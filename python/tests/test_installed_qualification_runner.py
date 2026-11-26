"""Real installed-wheel regression for subprocess identity and STEP qualification."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import paacs_mesh as pm


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "benchmarks/qualify_installed_step.py"


class InstalledQualificationRunnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        expected = os.environ.get("PAACS_EXPECT_INSTALLED") == "1"
        try:
            distribution = importlib.metadata.distribution("paacs-mesh")
            base = Path(distribution.locate_file("")).resolve()
            direct = distribution.read_text("direct_url.json")
            direct = json.loads(direct) if direct else {}
            installed = (Path(pm.__file__).resolve().is_relative_to(base)
                and bool(distribution.read_text("WHEEL")) and not direct.get("dir_info") and not direct.get("vcs_info"))
        except importlib.metadata.PackageNotFoundError:
            installed = False
        if not installed:
            if expected:
                raise AssertionError("CI expected an installed wheel; source/development imports are not an integration result")
            raise unittest.SkipTest("Real subprocess qualification requires an installed wheel, not a source Native development module")
        if not RUNNER.is_file():
            raise AssertionError("Installed-wheel qualification runner is missing from tested source")

    def test_owned_cylinder_passes_through_selected_installed_interpreter(self):
        spec = importlib.util.spec_from_file_location("installed_qualification_integration_helpers", RUNNER)
        runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        fixture = Path(os.environ.get("PAACS_TEST_STEP_FIXTURE", ROOT / "tests/fixtures/cad/cylinder.step"))
        self.assertTrue(fixture.is_file(), "Owned cylinder STEP fixture is required")
        with tempfile.TemporaryDirectory(prefix="paacs-installed-runner-") as temporary:
            root = Path(temporary)
            copied = root / "cylinder.step"
            shutil.copyfile(fixture, copied)
            plan = root / "plan.json"
            runner.save(plan, {"schema": runner.SCHEMA, "package_version": pm.__version__, "files": [{
                "id": "owned-cylinder", "path": copied.name, "sha256": runner.sha(copied), "size_mm": .6,
                "expected_original_solids": 1,
                "allowances": {"native_seconds": 120., "hard_seconds": 300., "inspection_seconds": 120.}}]})
            output = root / "results"
            # Exercise the actual CLI parsing/path-admission route and real worker.
            # Keeping the controller here avoids an uncontained outer redirector;
            # argv is the only mocked input, never PAACS or its Native execution.
            invocation = Path(sys.executable).absolute()
            with patch.object(sys, "argv", [str(RUNNER), "--python", str(invocation), "--plan", str(plan),
                    "--root", str(root), "--output", str(output)]):
                self.assertEqual(runner.main(), 0)
            admitted = runner.load(output / "plan.json")
            self.assertEqual(admitted["python"], str(invocation))
            controller = runner.load(output / "owned-cylinder/controller.json")
            self.assertEqual(controller["status"], "succeeded", json.dumps(controller, indent=2)[-16000:])
            self.assertTrue(controller["worker_acknowledged"])
            self.assertTrue(controller["worker_pid_verified"])
            self.assertEqual(controller["worker_pid"], controller["worker"]["pid"])
            self.assertEqual(controller["memory_measurement_pid"], controller["worker_pid"])
            self.assertEqual(controller["process_exit_code"], 0)
            self.assertEqual(controller["cleanup_errors"], [])
            self.assertEqual(controller["command"][0], str(invocation))
            self.assertEqual(Path(controller["worker"]["package"]["prefix"]).resolve(), Path(sys.prefix).resolve())
            self.assertTrue(Path(controller["worker"]["package"]["package_path"]).resolve().is_relative_to(Path(sys.prefix).resolve()))
            if os.name == "nt":
                self.assertTrue(controller["worker_identity"]["kill_on_close_job"])
                self.assertTrue(controller["process_memory_limit_enforced"])
                self.assertGreater(controller["peak_process_working_set_bytes"], 0)
                self.assertEqual(controller["worker_identity"]["chain"][0], controller["worker_pid"])
                self.assertEqual(controller["worker_identity"]["chain"][-1], controller["launcher_pid"])
            else:
                self.assertEqual(controller["worker_pid"], controller["launcher_pid"])
            report = runner.load(output / "owned-cylinder/native-report.json")
            self.assertTrue(report["audit"]["complete"] and report["audit"]["valid"])
            self.assertTrue(report["quality"]["complete"] and report["quality"]["valid"])
            self.assertEqual(report["quality"]["inverted_cells"], 0)
            self.assertEqual(report["source_sha256"], runner.sha(copied))
            self.assertEqual(controller["worker"]["native_status"], "qualified")
            for name in ("mesh.paacsmesh", "mesh.vtk", "quality-summary.json", "inspection/inspection.json"):
                self.assertTrue((output / "owned-cylinder" / name).is_file(), name)


if __name__ == "__main__":
    unittest.main()
