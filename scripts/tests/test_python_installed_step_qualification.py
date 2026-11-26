# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Exercise installed qualification guards without importing or running PAACS."""
import base64
from copy import deepcopy
import hashlib
import errno
import importlib.util
import io
import json
import os
from pathlib import Path
from types import SimpleNamespace
import sys
import subprocess
import tempfile
import time
import unittest
from unittest.mock import Mock, patch


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("installed_step_qualification", ROOT / "benchmarks/qualify_installed_step.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class InstalledStepQualificationTests(unittest.TestCase):
    def test_development_versions_remain_distinct_from_release_and_rc_pins(self):
        self.assertEqual(runner.normalized_version("1.3.0-dev.0"), "1.3.0.dev0")
        self.assertEqual(runner.normalized_version("1.3.0.dev0"), "1.3.0.dev0")
        self.assertEqual(runner.normalized_version("1.2.0-rc.4"), "1.2.0rc4")
        self.assertEqual(runner.normalized_version("1.2.0"), "1.2.0")
        self.assertNotEqual(runner.normalized_version("1.3.0.dev0"), runner.normalized_version("1.3.0"))
        self.assertNotEqual(runner.normalized_version("1.3.0.dev0"), runner.normalized_version("1.3.0rc0"))
        for invalid in ("1.3.0-dev", "1.3.0.dev0+untracked", "1.3.0-dev.0ignored", "1.3.0.dev0rc1"):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                runner.normalized_version(invalid)

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cad = self.root / "owned.step"
        self.cad.write_bytes(b"controlled CAD identity; no native mesher is invoked")
        self.row = {"id": "owned", "path": self.cad.name, "sha256": runner.sha(self.cad), "size_mm": 1.}
        self.plan = {"schema": runner.SCHEMA, "package_version": "1.2.0", "files": [self.row]}
        self.plan_path = self.root / "plan.json"

    def jobs(self):
        runner.save(self.plan_path, self.plan)
        return runner.read_plan(self.plan_path, self.root)

    def report(self):
        return {"status": "succeeded", "reason": "none", "source_sha256": self.row["sha256"], "mesh_published": True,
            "audit": {"complete": True, "valid": True, **{name: 0 for name in ("missing_boundary", "unexpected_boundary",
                "nonmanifold_faces", "inconsistent_orientation", "invalid_connectivity", "nonfinite_nodes")}},
            "quality": {"complete": True, "valid": True, "inverted_cells": 0}, "seconds": {"meshing_total": 1.}}

    def mesh(self):
        report = self.report()
        return SimpleNamespace(report=report, body_tags=[1], face_tags=[2, 3], node_count=4, tetrahedron_count=1,
            write_archive=lambda path: Path(path).write_bytes(b"controlled archive"),
            write_vtk=lambda path: Path(path).write_text("controlled VTK export"))

    def summary(self, archive=None):
        return {"mesh_sha256": hashlib.sha256(b"controlled archive").hexdigest(),
            "orientation_counts": {"inverted": 0, "degenerate": 0, "positive": 1}, "unclassified_tet4_count": 0,
            "body_tet4_counts": {"1": 1}, "sizing": {"achieved": {"boundary_faces": {"2": {}, "3": {}}}},
            "quality": {"mean_ratio": {"minimum": .5}}, "metric_conventions": {"reference": "independent controlled test"},
            "bounds": [[0, 0, 0], [1, 1, 1]]}

    @staticmethod
    def sharing_failure():
        error = PermissionError(errno.EACCES, "controlled Windows sharing failure")
        error.winerror = 32
        return error

    def test_atomic_save_retries_transient_sharing_and_preserves_exact_contents(self):
        target = self.root / "receipt.json"
        prior = b'{"prior":true}\n'
        target.write_bytes(prior)
        value, calls, clock = {"complete": True, "values": [1, 2, 3]}, [], [100.]
        replace, before = Path.replace, runner.io_snapshot()
        def replacing(source, destination):
            calls.append(destination)
            if len(calls) == 1:
                self.assertEqual(target.read_bytes(), prior)
                raise self.sharing_failure()
            return replace(source, destination)
        with patch.object(Path, "replace", new=replacing), patch.object(runner.time, "perf_counter", side_effect=lambda: clock[0]), \
                patch.object(runner.time, "sleep", side_effect=lambda delay: clock.__setitem__(0, clock[0]+delay)):
            runner.save(target, value)
        self.assertEqual(target.read_bytes(), (json.dumps(value, indent=2)+"\n").encode())
        self.assertEqual(len(calls), 2)
        self.assertEqual(runner.IO_METRICS["sharing_retries"]-before["sharing_retries"], 1)
        self.assertAlmostEqual(runner.IO_METRICS["retry_wait_seconds"]-before["retry_wait_seconds"], .02)

    def test_atomic_save_exhaustion_is_bounded_and_nonsharing_errors_are_not_retried(self):
        target = self.root / "receipt.json"
        prior = b'{"prior":true}\n'
        target.write_bytes(prior)
        before = runner.io_snapshot()
        with patch.object(Path, "replace", side_effect=self.sharing_failure()) as replace, patch.object(runner.time, "sleep") as sleep:
            with self.assertRaises(PermissionError) as caught:
                runner.save(target, {"new": True})
        self.assertEqual(replace.call_count, 5)
        self.assertEqual([item.args[0] for item in sleep.call_args_list], list(runner.SAVE_RETRY_DELAYS))
        self.assertLess(sum(runner.SAVE_RETRY_DELAYS), 1)
        self.assertEqual(caught.exception.paacs_atomic_io["attempts"], 5)
        self.assertEqual(target.read_bytes(), prior)
        self.assertEqual(runner.IO_METRICS["sharing_exhaustions"]-before["sharing_exhaustions"], 1)
        with patch.object(Path, "replace", side_effect=OSError(errno.EIO, "controlled disk failure")) as replace, \
                patch.object(runner.time, "sleep") as sleep:
            with self.assertRaises(OSError):
                runner.save(target, {"new": True})
        self.assertEqual(replace.call_count, 1)
        sleep.assert_not_called()
        self.assertEqual(target.read_bytes(), prior)

    def test_receipt_retains_its_own_retry_wait_and_final_exhaustion_cost(self):
        target = self.root / "worker.json"
        prior = b'{"prior":true}\n'
        target.write_bytes(prior)
        replace, before = Path.replace, runner.io_snapshot()
        clock = [100.]
        def replacing(source, destination):
            if destination == target:
                raise self.sharing_failure()
            return replace(source, destination)
        value = {"schema": "paacs-installed-step-attempt-v1", "status": "succeeded", "pid": 1}
        with patch.object(Path, "replace", new=replacing), patch.object(runner.time, "perf_counter", side_effect=lambda: clock[0]), \
                patch.object(runner.time, "sleep", side_effect=lambda delay: clock.__setitem__(0, clock[0]+delay)):
            final = runner.retain_final(target, value)
        self.assertEqual(target.read_bytes(), prior)
        self.assertEqual(final.name, "worker-final.json")
        retained = runner.load(final)
        self.assertEqual(retained["status"], "infrastructure_failed")
        self.assertEqual(retained["io"]["sharing_exhaustions"]-before["sharing_exhaustions"], 1)
        self.assertEqual(retained["io"]["failed_io_calls"]-before["failed_io_calls"], 1)
        self.assertAlmostEqual(retained["io"]["retry_wait_seconds"]-before["retry_wait_seconds"], .3)

    def test_successful_final_receipt_includes_its_own_transient_retry(self):
        target, clock, before = self.root / "worker.json", [100.], runner.io_snapshot()
        replace, attempts = Path.replace, []
        def replacing(source, destination):
            attempts.append(destination)
            if len(attempts) == 1:
                raise self.sharing_failure()
            return replace(source, destination)
        value = {"schema": "paacs-installed-step-attempt-v1", "status": "succeeded"}
        with patch.object(Path, "replace", new=replacing), patch.object(runner.time, "perf_counter", side_effect=lambda: clock[0]), \
                patch.object(runner.time, "sleep", side_effect=lambda delay: clock.__setitem__(0, clock[0]+delay)):
            runner.retain_final(target, value)
        retained = runner.load(target)
        self.assertEqual(retained["status"], "succeeded")
        self.assertEqual(retained["io"]["sharing_retries"]-before["sharing_retries"], 1)
        self.assertAlmostEqual(retained["io"]["retry_wait_seconds"]-before["retry_wait_seconds"], .02)

    def test_explicit_plan_is_finite_pinned_and_serial(self):
        job, = self.jobs()
        self.assertEqual(job["package_version"], "1.2.0")
        self.assertEqual(job["options"]["workers"], 1)
        self.assertEqual(job["options"]["validation_workers"], 1)
        self.assertNotIn("max_work", job["options"])
        self.assertEqual(job["options"]["max_seconds"], job["allowances"]["native_seconds"])
        self.cad.write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "hash differs"):
            runner.read_plan(self.plan_path, self.root)

    def test_plan_rejects_escape_ambiguous_sizing_and_invalid_workers(self):
        for path in ("../owned.step", "C:/owned.step", "/owned.step", "owned.step:stream"):
            with self.subTest(path=path), self.assertRaises(ValueError):
                runner.inside(self.root, path)
        mutations = ({"size_factor": .15}, {"size_mm": True}, {"options": {"workers": True}},
            {"options": {"validation_workers": 2}}, {"options": {"max_cells": 0}},
            {"options": {"growth": float("nan")}}, {"expected_original_solids": True},
            {"allowances": {"native_seconds": 20, "hard_seconds": 10}})
        for change in mutations:
            self.plan["files"] = [{**self.row, **change}]
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.jobs()

    def test_symlinked_root_accepts_owned_file_but_rejects_outside_symlink(self):
        actual = self.root / "actual-cad-root"
        actual.mkdir()
        alias = self.root / "aliased-cad-root"
        inside = actual / "inside.step"
        inside.write_bytes(b"owned inside file")
        outside = self.root / "outside.step"
        outside.write_bytes(b"outside the admitted CAD root")
        escape = actual / "escaped.step"
        try:
            alias.symlink_to(actual, target_is_directory=True)
            escape.symlink_to(outside)
        except (OSError, NotImplementedError) as error:
            self.skipTest("Filesystem symlink privilege unavailable: " + str(error))
        self.assertEqual(runner.inside(alias, "inside.step"), inside.resolve())
        self.plan["files"] = [{**self.row, "path": "inside.step", "sha256": runner.sha(inside)}]
        runner.save(self.plan_path, self.plan)
        job, = runner.read_plan(self.plan_path, alias)
        self.assertEqual(Path(job["source"]), inside.resolve())
        with self.assertRaisesRegex(ValueError, "escapes its root"):
            runner.inside(alias, "escaped.step")
        self.plan["files"] = [{**self.row, "path": "escaped.step", "sha256": runner.sha(outside)}]
        runner.save(self.plan_path, self.plan)
        with self.assertRaisesRegex(ValueError, "escapes its root"):
            runner.read_plan(self.plan_path, alias)

    def test_relative_root_is_canonicalized_before_comparing_resolved_target(self):
        actual = self.root / "relative-cad-root"
        actual.mkdir()
        inside = actual / "inside.step"
        inside.write_bytes(b"owned relative-root file")
        self.plan["files"] = [{**self.row, "path": "inside.step", "sha256": runner.sha(inside)}]
        runner.save(self.plan_path, self.plan)
        previous = Path.cwd()
        try:
            os.chdir(self.root)
            self.assertEqual(runner.inside(Path("relative-cad-root"), "inside.step"), inside.resolve())
            job, = runner.read_plan(self.plan_path, Path("relative-cad-root"))
            self.assertEqual(Path(job["source"]), inside.resolve())
        finally:
            os.chdir(previous)

    def test_python_invocation_validates_resolved_target_but_preserves_selected_path(self):
        physical = self.root / "base-python"
        physical.write_bytes(b"controlled interpreter target")
        selected = self.root / "selected-venv-python"
        # Simulate a different resolution without requiring Windows link rights.
        with patch.object(Path, "resolve", return_value=physical) as resolution:
            invocation = runner.python_invocation(selected)
        self.assertEqual(invocation, selected.absolute())
        resolution.assert_called_once_with(strict=True)
        self.assertNotEqual(invocation, physical)

    def test_cli_plan_only_preserves_symlinked_selected_interpreter(self):
        physical = self.root / "base-python"
        physical.write_bytes(b"controlled interpreter target; never executed")
        selected = self.root / "venv/bin/python"
        selected.parent.mkdir(parents=True)
        try:
            selected.symlink_to(physical)
        except (OSError, NotImplementedError) as error:
            self.skipTest("Filesystem symlink privilege unavailable: " + str(error))
        runner.save(self.plan_path, self.plan)
        output = self.root / "plan-only"
        with patch.object(sys, "argv", ["qualify_installed_step.py", "--python", str(selected),
                "--plan", str(self.plan_path), "--root", str(self.root), "--output", str(output), "--plan-only"]), \
                patch.object(sys, "stdout", new=io.StringIO()), patch.object(runner.subprocess, "Popen") as launched:
            self.assertEqual(runner.main(), 0)
        launched.assert_not_called()
        self.assertEqual(runner.load(output / "plan.json")["python"], str(selected.absolute()))
        self.assertNotEqual(selected.absolute(), selected.resolve())

    def test_union_diagonal_sizing_paginates_and_includes_separated_bodies(self):
        pages = {0: {"offset": 0, "total": 2, "next_offset": 1, "entities": [{"tag": 1,
            "bounds": {"lower": [0, 0, 0], "upper": [3, 4, 0]}}]},
            1: {"offset": 1, "total": 2, "next_offset": None, "entities": [{"tag": 2,
                "bounds": {"lower": [6, 0, 0], "upper": [6, 8, 0]}}]}}
        model = SimpleNamespace(body_tags=[1, 2], inspect=lambda kind, offset, limit: pages[offset])
        size, basis = runner.resolve_size(model, {"size_mm": None, "size_factor": .15})
        self.assertAlmostEqual(size, 1.5)  # union diagonal=10; longest extent=8 would be wrong.
        self.assertEqual(basis["diagonal_mm"], 10)
        pages[1]["entities"][0]["tag"] = 1
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            runner.resolve_size(model, {"size_mm": None, "size_factor": .15})

    def test_loader_and_cad_overrides_are_removed(self):
        with patch.dict(os.environ, {"PYTHONPATH": "source", "LD_PRELOAD": "inject", "LD_AUDIT": "inject",
                "CSF_STEPDefaults": "sdk", "OpenCASCADE_DIR": "sdk", "DYLD_LIBRARY_PATH": "sdk"}):
            clean = runner.clean_environment(Path(sys.executable))
        for key in ("PYTHONPATH", "LD_PRELOAD", "LD_AUDIT", "CSF_STEPDefaults", "OpenCASCADE_DIR", "DYLD_LIBRARY_PATH"):
            self.assertNotIn(key, clean)
        self.assertEqual(clean["PYTHONNOUSERSITE"], "1")

    def test_native_and_archive_gates_preserve_all_validity_requirements(self):
        mesh = self.mesh()
        runner.validate_native(mesh.report, self.row["sha256"], [1], [2, 3], mesh)
        for section, name, value in [("audit", "complete", False), ("audit", "valid", False),
                ("quality", "complete", False), ("quality", "valid", False), ("quality", "inverted_cells", 1)] + [
                ("audit", name, 1) for name in ("missing_boundary", "unexpected_boundary", "nonmanifold_faces",
                    "inconsistent_orientation", "invalid_connectivity", "nonfinite_nodes")]:
            bad = deepcopy(mesh.report)
            bad[section][name] = value
            with self.subTest(section=section, name=name), self.assertRaises(ValueError):
                runner.validate_native(bad, self.row["sha256"], [1], [2, 3], mesh)
        archive = self.root / "mesh.paacsmesh"
        mesh.write_archive(archive)
        runner.validate_summary(self.summary(), archive, [1], [2, 3])
        for change in ({"mesh_sha256": "0"*64}, {"orientation_counts": {"inverted": 1, "degenerate": 0}},
                {"unclassified_tet4_count": 1}, {"body_tet4_counts": {"2": 1}}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                runner.validate_summary({**self.summary(), **change}, archive, [1], [2, 3])

    def installed_fixture(self):
        base = self.root / "site-packages"
        package = base / "paacs_mesh/__init__.py"
        native = base / "paacs_mesh/_native.pyd"
        package.parent.mkdir(parents=True)
        package.write_bytes(b"installed public wrapper")
        native.write_bytes(b"installed public extension")
        runtime = base / "paacs_mesh/cad.py"
        runtime.write_bytes(b"installed public CAD module")
        library = base / "paacs_mesh.libs/kernel.dll"
        library.parent.mkdir()
        library.write_bytes(b"installed CAD runtime")
        files = []
        for path in (package, native, runtime, library):
            value = base64.urlsafe_b64encode(hashlib.sha256(path.read_bytes()).digest()).decode().rstrip("=")
            class File:
                def __init__(self, name, digest):
                    self.name, self.hash = name, SimpleNamespace(mode="sha256", value=digest)
                def __str__(self):
                    return self.name
            files.append(File(path.relative_to(base).as_posix(), value))
        metadata = {"METADATA": "Name: paacs-mesh\nVersion: 1.2.0", "WHEEL": "Wheel-Version: 1.0",
            "RECORD": "controlled records", "INSTALLER": "pip", "direct_url.json": None}
        distribution = SimpleNamespace(version="1.2.0", files=files, locate_file=lambda name: base,
            read_text=lambda name: metadata.get(name))
        pm = SimpleNamespace(__file__=str(package), _native=SimpleNamespace(__file__=str(native)),
            __version__="1.2.0", native_version="1.2.0")
        return pm, distribution, metadata, native

    def test_distribution_provenance_refuses_stale_native_editable_and_record_mutation(self):
        pm, distribution, metadata, native = self.installed_fixture()
        with patch.object(runner.importlib.metadata, "distribution", return_value=distribution), patch.object(sys, "prefix", str(self.root)):
            identity = runner.package_identity(pm, "1.2.0")
            self.assertEqual(identity["native_sha256"], runner.sha(native))
            pm.native_version = "1.2.0-rc.4"
            with self.assertRaisesRegex(ValueError, "version"):
                runner.package_identity(pm, "1.2.0")
            pm.native_version = "1.2.0"
            metadata["direct_url.json"] = json.dumps({"url": "file:///source", "dir_info": {"editable": True}})
            with self.assertRaisesRegex(ValueError, "Editable"):
                runner.package_identity(pm, "1.2.0")
            metadata["direct_url.json"] = None
            native.write_bytes(b"mutated extension")
            with self.assertRaisesRegex(ValueError, "RECORD"):
                runner.package_identity(pm, "1.2.0")

    def test_all_hashed_python_and_cad_runtime_files_are_guarded(self):
        pm, distribution, metadata, native = self.installed_fixture()
        with patch.object(runner.importlib.metadata, "distribution", return_value=distribution), patch.object(sys, "prefix", str(self.root)):
            self.assertEqual(runner.package_identity(pm, "1.2.0")["verified_record_files"], 4)
            path = Path(pm.__file__).with_name("cad.py")
            path.write_bytes(b"modified CAD behavior")
            with self.assertRaisesRegex(ValueError, "runtime differs"):
                runner.package_identity(pm, "1.2.0")

    def test_missing_archive_or_worker_success_is_never_a_qualified_result(self):
        good = {"status": "succeeded", "native_status": "qualified"}
        self.assertEqual(runner.final_status(good, 0, None, self.root), "incomplete_result")
        for name in ("mesh.paacsmesh", "mesh.vtk", "native-report.json", "quality-summary.json", "inspection/inspection.json"):
            path = self.root / name
            path.parent.mkdir(exist_ok=True)
            path.write_bytes(b"controlled")
        self.assertEqual(runner.final_status(good, 0, None, self.root), "succeeded")
        self.assertEqual(runner.final_status(good, 0, "hard_timeout", self.root), "hard_timeout_native_outcome_not_recertified")
        self.assertEqual(runner.final_status({}, 0, None, self.root), "incomplete_result")
        self.assertEqual(runner.final_status({}, 3, None, self.root), "process_failed")
        self.assertEqual(runner.final_status({"status": "native_refused"}, 0, None, self.root), "native_refused")

    def test_request_errors_are_not_counted_as_geometry_failures(self):
        control = RuntimeError("unsupported parameter")
        control.reason = "unsupported-control"
        self.assertEqual(runner.refusal_kind(control, {}), ("request_rejected", "not_admitted"))
        self.assertEqual(runner.refusal_kind(TypeError("bad parameter"), {}), ("request_rejected", "not_admitted"))
        self.assertEqual(runner.refusal_kind(RuntimeError("boundary failed"), {"reason": "boundary_recovery_failed"}),
            ("native_refused", "refused"))

    class ControlledProcesses:
        def __init__(self, launcher):
            self.launcher, self.worker = launcher, None
            self.closed = self.stopped = False
        def contain_launcher(self):
            pass
        def establish(self, pid):
            if pid != self.launcher.pid:
                raise ValueError("controlled fixture rejects unknown PID")
            self.worker = self.launcher
            return {"worker_pid": pid, "launcher_pid": pid, "chain": [pid], "kind": "controlled direct fixture"}
        def alive(self):
            return self.launcher.poll() is None
        def exit_code(self):
            return self.launcher.poll()
        def cleanup(self, stop=False):
            self.stopped = stop
            if stop:
                self.launcher.kill()
            return []
        def close(self):
            self.closed = True

    def controlled_child(self, command, **kwargs):
        self.assertEqual(command[1], "-I")
        output = kwargs["cwd"]
        self.assertEqual(kwargs["env"]["TEMP"], str(output / "scratch"))
        job = runner.load(command[-1])
        runner.save(output / "worker-hello.json", {"schema": "paacs-installed-worker-hello-v1", "pid": 8123,
            "parent_pid": job["controller_pid"], "token": job["handshake_token"], "job_sha256": runner.sha(command[-1])})
        for name in ("mesh.paacsmesh", "mesh.vtk", "native-report.json", "quality-summary.json", "inspection/inspection.json"):
            target = output / name
            target.parent.mkdir(exist_ok=True)
            target.write_bytes(b"controlled fixture output")
        runner.save(output / "worker.json", {"pid": 8123, "status": "succeeded", "native_status": "qualified"})
        child = SimpleNamespace(pid=8123, forced=False, returncode=None)
        child.poll = lambda: 1 if child.forced else 0 if (output / "controller-ack.json").is_file() else None
        child.wait = lambda timeout=None: child.poll()
        child.kill = lambda: setattr(child, "forced", True)
        return child

    def test_controller_checks_final_disk_even_for_immediately_completed_worker(self):
        job, = self.jobs()
        job["allowances"]["max_disk_bytes"] = 1
        with patch.object(runner.subprocess, "Popen", side_effect=self.controlled_child), \
                patch.object(runner, "WindowsMemory", return_value=SimpleNamespace(sample=lambda child: None)), \
                patch.object(runner, "ProcessControl", self.ControlledProcesses):
            result = runner.run_attempt(Path(sys.executable), job, self.root / "attempt")
        self.assertEqual(result["status"], "output_disk_stop_native_outcome_not_recertified")
        self.assertGreater(result["durable_output_bytes"], 1)
        self.assertFalse(result["process_memory_limit_enforced"])

    def test_controller_retains_start_failure_and_rejects_missing_pid(self):
        job, = self.jobs()
        with patch.object(runner.subprocess, "Popen", side_effect=OSError("controlled loader failure")):
            failed = runner.run_attempt(Path(sys.executable), job, self.root / "loader")
        self.assertEqual(failed["status"], "process_start_failed")
        self.assertTrue((self.root / "loader/controller.json").is_file())
        def missing_pid(command, **kwargs):
            child = self.controlled_child(command, **kwargs)
            runner.save(kwargs["cwd"] / "worker.json", {"status": "succeeded", "native_status": "qualified"})
            return child
        with patch.object(runner.subprocess, "Popen", side_effect=missing_pid), \
                patch.object(runner, "WindowsMemory", return_value=SimpleNamespace(sample=lambda child: None)), \
                patch.object(runner, "ProcessControl", self.ControlledProcesses):
            result = runner.run_attempt(Path(sys.executable), job, self.root / "missing-pid")
        self.assertEqual(result["status"], "input_or_worker_identity_failed")

    def test_post_spawn_receipt_exhaustion_cleans_up_controlled_process(self):
        job, = self.jobs()
        replace, controls = Path.replace, []
        def replacing(source, destination):
            if destination.name == "controller.json" and json.loads(source.read_text()).get("pid") is not None:
                raise self.sharing_failure()
            return replace(source, destination)
        def supervise(child):
            control = self.ControlledProcesses(child)
            controls.append(control)
            return control
        with patch.object(Path, "replace", new=replacing), patch.object(runner.time, "sleep"), \
                patch.object(runner.subprocess, "Popen", side_effect=self.controlled_child), \
                patch.object(runner, "WindowsMemory", return_value=SimpleNamespace(sample=lambda child: None)), \
                patch.object(runner, "ProcessControl", side_effect=supervise):
            result = runner.run_attempt(Path(sys.executable), job, self.root / "locked-controller")
        self.assertTrue(controls[0].stopped)
        self.assertTrue(controls[0].closed)
        self.assertFalse(result["worker_acknowledged"])
        self.assertEqual(result["status"], "infrastructure_failed")
        self.assertTrue((self.root / "locked-controller/controller-final.json").is_file())

    def test_post_spawn_io_failure_stays_infrastructure_when_final_primary_write_recovers(self):
        job, = self.jobs()
        replace, failures = Path.replace, [0]
        def replacing(source, destination):
            if destination.name == "controller.json" and json.loads(source.read_text()).get("pid") is not None and failures[0] < 5:
                failures[0] += 1
                raise self.sharing_failure()
            return replace(source, destination)
        with patch.object(Path, "replace", new=replacing), patch.object(runner.time, "sleep"), \
                patch.object(runner.subprocess, "Popen", side_effect=self.controlled_child), \
                patch.object(runner, "WindowsMemory", return_value=SimpleNamespace(sample=lambda child: None)), \
                patch.object(runner, "ProcessControl", self.ControlledProcesses):
            result = runner.run_attempt(Path(sys.executable), job, self.root / "controller-recovers")
        self.assertEqual(result["status"], "startup_infrastructure_failed")
        self.assertFalse(result["worker_acknowledged"])
        self.assertEqual(runner.load(self.root / "controller-recovers/controller.json")["status"], "startup_infrastructure_failed")

    def test_worker_checkpoint_exhaustion_is_infrastructure_before_mesh_entry(self):
        for phase in ("controller_handshake", "meshing"):
            with self.subTest(phase=phase):
                job, = self.jobs()
                output = self.root / ("locked-" + phase)
                output.mkdir()
                job["output"] = str(output)
                path = output / "job.json"
                runner.save(path, job)
                model = SimpleNamespace(body_tags=[1], face_tags=[2, 3], info={}, source_sha256=self.row["sha256"], mesh=Mock())
                pm = SimpleNamespace(read_step=Mock(return_value=model))
                replace = Path.replace
                def replacing(source, destination):
                    value = json.loads(source.read_text())
                    if destination.name == "worker.json" and value.get("phase") == phase and value.get("status") == "running":
                        raise self.sharing_failure()
                    return replace(source, destination)
                with patch.object(Path, "replace", new=replacing), patch.object(runner.time, "sleep"), \
                        patch.dict(sys.modules, {"paacs_mesh": pm}), patch.object(runner, "await_controller"), \
                        patch.object(runner, "package_identity", return_value={}) as identity:
                    self.assertEqual(runner.worker(path), 1)
                retained = runner.load(output / "worker.json")
                self.assertEqual(retained["status"], "infrastructure_failed")
                self.assertEqual(retained["native_status"], "not_attempted")
                self.assertFalse(retained["mesh_call_started"])
                self.assertEqual(retained["io_failure"]["attempts"], 5)
                model.mesh.assert_not_called()
                if phase == "controller_handshake":
                    identity.assert_not_called()

    def test_disk_sampling_tolerates_atomic_temporary_bundle_rename(self):
        gone = SimpleNamespace(is_file=lambda: True, stat=lambda: (_ for _ in ()).throw(FileNotFoundError()))
        present = SimpleNamespace(is_file=lambda: True, stat=lambda: SimpleNamespace(st_size=17))
        root = SimpleNamespace(rglob=lambda pattern: [gone, present])
        self.assertEqual(runner.disk_bytes(root), 17)

    def test_creator_chain_rejects_foreign_cycles_and_reused_process_ids(self):
        chain = runner.ancestry_chain({30: 20, 20: 10}, 30, 10)
        self.assertEqual(chain, [30, 20, 10])
        runner.validate_creation_chain(chain, {10: (100, 500), 20: (200, 0), 30: (300, 0)})
        for parents in ({30: 99, 99: 0}, {30: 20, 20: 30}, {}):
            with self.subTest(parents=parents), self.assertRaises(ValueError):
                runner.ancestry_chain(parents, 30, 10)
        for times in ({10: (100, 500), 20: (200, 250), 30: (300, 0)},
                      {10: (400, 0), 20: (200, 0), 30: (300, 0)}):
            with self.subTest(times=times), self.assertRaises(ValueError):
                runner.validate_creation_chain(chain, times)

    def test_worker_handle_has_job_assignment_and_monitoring_rights(self):
        open_process = Mock(return_value=987)
        api = SimpleNamespace(kernel=SimpleNamespace(OpenProcess=open_process, CloseHandle=Mock()),
            checked=lambda value: value, lifetime=lambda handle: (100, 0))
        reference = runner.WindowsProcess(api, 8123)
        mask, inherit, pid = open_process.call_args.args
        self.assertEqual(pid, 8123)
        self.assertFalse(inherit)
        self.assertEqual(mask & 0x0101, 0x0101)  # SET_QUOTA + TERMINATE for job assignment.
        self.assertEqual(mask & 0x00100010, 0x00100010)  # SYNCHRONIZE + VM_READ.
        reference.close()
        api.kernel.CloseHandle.assert_called_once_with(987)

    def test_windows_creator_check_refreshes_exited_launcher_identity_before_assignment(self):
        control = runner.ProcessControl.__new__(runner.ProcessControl)
        control.api = SimpleNamespace(parents=lambda: {30: 10}, lifetime=lambda handle: (100, 150))
        control.root = SimpleNamespace(_handle=110, lifetime=(100, 0))  # Cached while originally alive.
        control.launcher = SimpleNamespace(pid=10)
        control.worker = None
        control.job = SimpleNamespace(assign=Mock())
        candidate = SimpleNamespace(pid=30, _handle=330, lifetime=(200, 0), alive=lambda: True, close=Mock())
        with patch.object(runner, "WindowsProcess", return_value=candidate):
            with self.assertRaisesRegex(ValueError, "reused PID"):
                control.establish(30)
        control.job.assign.assert_not_called()
        candidate.close.assert_called_once()

    def test_handshake_requires_matching_token_job_hash_and_pid(self):
        job = {"handshake_token": "controlled-token"}
        hello = {"schema": "paacs-installed-worker-hello-v1", "token": "controlled-token", "pid": 8123, "job_sha256": "a"*64}
        runner.validate_hello(hello, job, "a"*64)
        for change in ({"token": "foreign-token"}, {"job_sha256": "b"*64}, {"pid": True}, {"pid": 0}, {"schema": "foreign"}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                runner.validate_hello({**hello, **change}, job, "a"*64)

    def test_worker_acknowledgement_is_a_barrier_before_any_paacs_work(self):
        job, = self.jobs()
        output = self.root / "no-ack"
        output.mkdir()
        job.update(output=str(output), handshake_token="controlled", controller_pid=os.getpid())
        path = output / "job.json"
        runner.save(path, job)
        pm = SimpleNamespace(read_step=Mock(side_effect=AssertionError("Native must never run before acknowledgement")))
        with patch.dict(sys.modules, {"paacs_mesh": pm}), patch.object(runner, "package_identity") as identity, \
                patch.object(runner.time, "monotonic", side_effect=[0, 31]):
            self.assertEqual(runner.worker(path), 1)
        pm.read_step.assert_not_called()
        identity.assert_not_called()
        self.assertTrue((output / "worker-hello.json").is_file())
        self.assertEqual(runner.load(output / "worker.json")["native_status"], "not_attempted")

    def test_cleanup_uses_owned_job_and_retained_worker_if_launcher_exits(self):
        control = runner.ProcessControl.__new__(runner.ProcessControl)
        control.api = True
        control.job = SimpleNamespace(terminate=Mock(), close=Mock())
        control.worker = SimpleNamespace(terminate=Mock(), wait=Mock(), close=Mock())
        control.cleanup_members = []
        control.contain_pending_descendants = Mock(return_value=[])
        control.launcher = SimpleNamespace(pid=10, poll=lambda: 0, kill=Mock())
        with patch.object(runner.subprocess, "run") as taskkill:
            self.assertEqual(control.cleanup(stop=True), [])
        taskkill.assert_not_called()  # Never target a dead/reusable launcher PID.
        control.job.terminate.assert_called_once()
        control.worker.terminate.assert_called_once()
        control.worker.wait.assert_called_once()
        control.close()
        control.job.close.assert_called_once()
        control.worker.close.assert_called_once()

    @unittest.skipUnless(os.name == "nt", "Real Windows process/job startup-cleanup regression")
    def test_real_pre_handshake_child_is_contained_and_terminated_without_native_import(self):
        pid_file = self.root / "sleeping-child.json"
        sleeper_code = "import json,os,sys,time; from pathlib import Path; Path(sys.argv[1]).write_text(json.dumps({'pid':os.getpid()})); time.sleep(60)"
        parent_code = "import subprocess,sys; child=subprocess.Popen([sys.executable,'-I','-c',sys.argv[1],sys.argv[2]],creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0)); child.wait()"
        launcher = subprocess.Popen([sys.executable, "-I", "-c", parent_code, sleeper_code, str(pid_file)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, creationflags=subprocess.CREATE_NO_WINDOW)
        control = reference = None
        try:
            deadline = time.monotonic() + 10
            while not pid_file.is_file() and time.monotonic() < deadline:
                self.assertIsNone(launcher.poll(), "Controlled stand-in launcher exited unexpectedly")
                time.sleep(.02)
            self.assertTrue(pid_file.is_file(), "Controlled child did not start")
            pid = json.loads(pid_file.read_text())["pid"]
            # Deliberately assign the root only AFTER its child exists, and never
            # emit a worker hello. This reproduces the venv startup cleanup race.
            control = runner.ProcessControl(launcher)
            reference = runner.WindowsProcess(control.api, pid)
            control.contain_launcher()
            present = control.api.w.BOOL()
            control.api.checked(control.api.kernel.IsProcessInJob(reference._handle, control.job._handle, runner.ctypes.byref(present)))
            self.assertFalse(present.value, "Pre-existing stand-in child must begin outside the new job")
            self.assertIsNone(control.worker)
            warnings = control.cleanup(stop=True)
            self.assertTrue(all(message.startswith("Descendant capture refused:") for message in warnings), warnings)
            self.assertFalse(reference.alive(), "Pre-hello child survived controlled cleanup")
            self.assertIsNotNone(launcher.poll())
            self.assertTrue(control.cleanup_members, "No verified pre-existing descendants were captured")
            self.assertTrue(all(not child.alive() for child in control.cleanup_members), "A verified captured process survived")
        finally:
            if control is not None:
                try:
                    control.cleanup(stop=True)
                finally:
                    control.close()
            else:
                runner.terminate_launcher(launcher)
            if reference is not None:
                try:
                    reference.terminate()
                    reference.wait()
                finally:
                    reference.close()

    def test_controller_cleans_up_when_handshake_never_arrives(self):
        job, = self.jobs()
        controls = []
        def start(command, **kwargs):
            child = SimpleNamespace(pid=8123, forced=False)
            child.poll = lambda: 1 if child.forced else None
            child.kill = lambda: setattr(child, "forced", True)
            return child
        def supervise(child):
            control = self.ControlledProcesses(child)
            controls.append(control)
            return control
        ticks = iter(range(0, 1000, 31))
        with patch.object(runner.subprocess, "Popen", side_effect=start), \
                patch.object(runner, "WindowsMemory", return_value=SimpleNamespace(sample=lambda child: None)), \
                patch.object(runner, "ProcessControl", side_effect=supervise), \
                patch.object(runner.time, "perf_counter", side_effect=lambda: next(ticks)):
            result = runner.run_attempt(Path(sys.executable), job, self.root / "never-hello")
        self.assertEqual(result["stop"], "worker_handshake_timeout")
        self.assertFalse(result["worker_acknowledged"])
        self.assertTrue(controls[0].stopped)
        self.assertTrue(controls[0].closed)

    def test_completed_worker_after_hard_deadline_is_not_recertified(self):
        job, = self.jobs()
        job["allowances"]["hard_seconds"] = 1
        clock = [0.]
        class LateCleanup(self.ControlledProcesses):
            def cleanup(inner, stop=False):
                result = super().cleanup(stop=stop)
                clock[0] = 1.5
                return result
        with patch.object(runner.subprocess, "Popen", side_effect=self.controlled_child), \
                patch.object(runner, "WindowsMemory", return_value=SimpleNamespace(sample=lambda child: None)), \
                patch.object(runner, "ProcessControl", LateCleanup), \
                patch.object(runner.time, "perf_counter", side_effect=lambda: clock[0]):
            result = runner.run_attempt(Path(sys.executable), job, self.root / "late")
        self.assertEqual(result["status"], "hard_timeout_native_outcome_not_recertified")

    def test_startup_final_write_exhaustion_retains_bounded_log_infrastructure_signal(self):
        job, = self.jobs()
        output = self.root / "all-worker-receipts-locked"
        output.mkdir()
        job.update(output=str(output), handshake_token="controlled-startup-token", controller_pid=os.getpid())
        path = output / "job.json"
        runner.save(path, job)
        replace, before, stream = Path.replace, runner.io_snapshot(), io.StringIO()
        def replacing(source, destination):
            if destination.name in ("worker.json", "worker-final.json"):
                raise self.sharing_failure()
            return replace(source, destination)
        with patch.object(Path, "replace", new=replacing), patch.object(runner.time, "sleep"), \
                patch.object(sys, "stderr", new=stream), patch.object(runner, "package_identity") as identity:
            with self.assertRaises(PermissionError):
                runner.worker(path)
        identity.assert_not_called()
        (output / "worker.log").write_text(stream.getvalue(), encoding="utf-8")
        retained = runner.worker_receipt(output, None, job)
        self.assertEqual(retained["status"], "infrastructure_failed")
        self.assertFalse(retained["mesh_call_started"])
        self.assertEqual(retained["native_status"], "not_attempted")
        self.assertEqual(retained["io"]["sharing_exhaustions"]-before["sharing_exhaustions"], 3)

    def test_terminal_log_requires_exact_job_token_source_and_owned_pid(self):
        output = self.root / "terminal-auth"
        output.mkdir()
        job = {"handshake_token": "owned-token", "source_sha256": "a"*64}
        signal = {"schema": "paacs-installed-worker-terminal-io-failure-v1", "pid": 8123,
            "status": "infrastructure_failed", "mesh_call_started": True, "token": "owned-token", "source_sha256": "a"*64}
        fallback = {"pid": 8123, "status": "running"}
        runner.save(output / "worker.json", fallback)
        for changes in ({"token": "foreign"}, {"source_sha256": "b"*64}, {"pid": 9999}):
            (output / "worker.log").write_text(json.dumps({**signal, **changes})+"\n")
            self.assertEqual(runner.worker_receipt(output, 8123, job), fallback)
        (output / "worker.log").write_text(json.dumps(signal)+"\n")
        self.assertEqual(runner.worker_receipt(output, 8123, job), signal)

    def test_worker_executes_public_contract_and_retains_native_refusal(self):
        job, = self.jobs()
        output = self.root / "attempt"
        output.mkdir()
        job["output"] = str(output)
        path = output / "job.json"
        runner.save(path, job)
        mesh = self.mesh()
        model = SimpleNamespace(body_tags=[1], face_tags=[2, 3], info={}, source_sha256=self.row["sha256"],
            mesh=lambda *args, **kwargs: mesh)
        class Inspection:
            def __enter__(inner): return inner
            def __exit__(inner, *args): return None
            def summary(inner): return self.summary()
            def save(inner, path):
                path.mkdir()
                (path / "inspection.json").write_text("{}")
        pm = SimpleNamespace(read_step=lambda source: model, diagnose=lambda report: {}, inspect_mesh=lambda *a, **k: Inspection())
        with patch.dict(sys.modules, {"paacs_mesh": pm}), patch.object(runner, "await_controller"), patch.object(runner, "package_identity", return_value={"native_path": str(self.cad), "native_sha256": runner.sha(self.cad), "runtime_contents_sha256": "0"*64}):
            self.assertEqual(runner.worker(path), 0)
            self.assertEqual(runner.load(output / "worker.json")["status"], "succeeded")
            error = RuntimeError("controlled geometry refusal")
            error.reason, error.report = "boundary_recovery_failed", {"reason": "boundary_recovery_failed"}
            def refused(*args, **kwargs): raise error
            model.mesh = refused
            self.assertEqual(runner.worker(path), 0)
            receipt = runner.load(output / "worker.json")
            self.assertEqual(receipt["status"], "native_refused")
            self.assertEqual(receipt["reason"], "boundary_recovery_failed")
            self.assertTrue((output / "native-report.json").is_file())


if __name__ == "__main__":
    unittest.main()
