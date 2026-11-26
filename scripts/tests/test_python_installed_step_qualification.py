# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Exercise installed qualification guards without importing or running PAACS."""
import base64
from copy import deepcopy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("installed_step_qualification", ROOT / "benchmarks/qualify_installed_step.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class InstalledStepQualificationTests(unittest.TestCase):
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

    def controlled_child(self, command, **kwargs):
        self.assertEqual(command[1], "-I")
        output = kwargs["cwd"]
        self.assertEqual(kwargs["env"]["TEMP"], str(output / "scratch"))
        for name in ("mesh.paacsmesh", "mesh.vtk", "native-report.json", "quality-summary.json", "inspection/inspection.json"):
            target = output / name
            target.parent.mkdir(exist_ok=True)
            target.write_bytes(b"controlled fixture output")
        runner.save(output / "worker.json", {"pid": 8123, "status": "succeeded", "native_status": "qualified"})
        return SimpleNamespace(pid=8123, returncode=0, poll=lambda: 0, wait=lambda: 0)

    def test_controller_checks_final_disk_even_for_immediately_completed_worker(self):
        job, = self.jobs()
        job["allowances"]["max_disk_bytes"] = 1
        with patch.object(runner.subprocess, "Popen", side_effect=self.controlled_child), \
                patch.object(runner, "WindowsMemory", return_value=SimpleNamespace(sample=lambda child: None)):
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
                patch.object(runner, "WindowsMemory", return_value=SimpleNamespace(sample=lambda child: None)):
            result = runner.run_attempt(Path(sys.executable), job, self.root / "missing-pid")
        self.assertEqual(result["status"], "input_or_worker_identity_failed")

    def test_disk_sampling_tolerates_atomic_temporary_bundle_rename(self):
        gone = SimpleNamespace(is_file=lambda: True, stat=lambda: (_ for _ in ()).throw(FileNotFoundError()))
        present = SimpleNamespace(is_file=lambda: True, stat=lambda: SimpleNamespace(st_size=17))
        root = SimpleNamespace(rglob=lambda pattern: [gone, present])
        self.assertEqual(runner.disk_bytes(root), 17)

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
        with patch.dict(sys.modules, {"paacs_mesh": pm}), patch.object(runner, "package_identity", return_value={"native_path": str(self.cad), "native_sha256": runner.sha(self.cad), "runtime_contents_sha256": "0"*64}):
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
