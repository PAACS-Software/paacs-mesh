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
            output.mkdir()
            job, = runner.read_plan(plan, root)
            # Run the controller in this test process, and the real Native worker
            # through sys.executable. An extra outer Windows redirector would
            # itself need containment if the test's outer timeout fired.
            runner.run_attempt(Path(sys.executable), job, output / "owned-cylinder")
            controller = runner.load(output / "owned-cylinder/controller.json")
            self.assertEqual(controller["status"], "succeeded", json.dumps(controller, indent=2)[-16000:])
            self.assertTrue(controller["worker_acknowledged"])
            self.assertTrue(controller["worker_pid_verified"])
            self.assertEqual(controller["worker_pid"], controller["worker"]["pid"])
            self.assertEqual(controller["memory_measurement_pid"], controller["worker_pid"])
            self.assertEqual(controller["process_exit_code"], 0)
            self.assertEqual(controller["cleanup_errors"], [])
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
