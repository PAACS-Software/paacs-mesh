"""Installed-command checks for discovery, structured output and saved jobs."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import paacs_mesh as pm

from test_native import fixture_path


class DiscoveryCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = fixture_path()
        cls.temporary = tempfile.TemporaryDirectory(prefix="paacs-cli-tests-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.directory = Path(cls.temporary.name)
        cls.model = pm.read_step(cls.fixture)
        cls.mesh = cls.model.mesh(.7)
        cls.report = cls.directory / "native-report.json"
        cls.report.write_text(json.dumps(cls.mesh.report, allow_nan=False), encoding="utf-8")
        cls.recipe = cls.directory / "recipe.json"
        pm.save_recipe(cls.model, pm.create_recipe(cls.model, .7), cls.recipe)

    def command(self, *arguments, expected=0, machine=True):
        child = subprocess.run(
            [sys.executable, "-I", "-m", "paacs_mesh", *map(str, arguments)],
            capture_output=True, text=True, timeout=90,
        )
        self.assertEqual(child.returncode, expected, child.stdout + child.stderr)
        if machine:
            # Parse the complete stream, including any leading OCCT diagnostics.
            # A caller must never need to search for the JSON's opening brace.
            return json.loads(child.stdout)
        return child

    def test_capabilities_are_discoverable_as_json_and_text(self):
        report = self.command("capabilities", "--json")
        self.assertEqual(report, pm.capabilities())
        self.assertEqual(report["schema"], "paacs-mesh-capabilities-v1")
        child = self.command("capabilities", machine=False)
        self.assertIn("PAACS Mesh capabilities", child.stdout)
        self.assertNotIn("Traceback", child.stderr)

    def test_inspection_is_paginated_and_source_bound(self):
        first = self.command("inspect", self.fixture, "--limit", 1, "--topology", "--json")
        expected = self.model.inspect(limit=1, topology=True)
        self.assertEqual(first, expected)
        second = self.command("inspect", self.fixture, "--offset", first["next_offset"],
                              "--limit", 1, "--json")
        self.assertNotEqual(first["entities"][0]["tag"], second["entities"][0]["tag"])
        self.assertEqual(second["source_sha256"], first["source_sha256"])
        text = self.command("inspect", self.fixture, "--limit", 1, machine=False).stdout
        self.assertIn("CAD faces (mm)", text)
        self.assertIn("Next page: --offset 1", text)

    def test_invalid_inspection_arguments_and_missing_source_have_stable_codes(self):
        invalid = self.command("inspect", self.fixture, "--limit", 1001, "--json", expected=2)
        self.assertEqual(invalid["reason"], "invalid_argument")
        missing = self.command("inspect", self.directory / "absent.step", "--json", expected=1)
        self.assertEqual(missing["schema"], "paacs-mesh-cli-error-v1")
        self.assertEqual(missing["stage"], "import")
        usage = self.command("inspect", "--kind", "nonsense", "--json", expected=2)
        self.assertEqual(usage["reason"], "invalid_argument")

    def test_diagnose_uses_same_authoritative_report_as_python(self):
        report = self.command("diagnose", self.report, "--json")
        self.assertEqual(report, pm.diagnose(self.mesh.report))
        self.assertIn(report["status"], ("valid", "valid_with_warnings"))
        text = self.command("diagnose", self.report, machine=False)
        self.assertEqual(text.stdout.strip(), pm.format_diagnostics(report).strip())

    def test_diagnose_refuses_invalid_json_and_does_not_reaudit(self):
        malformed = self.directory / "malformed.json"
        for payload in ('{"status":1,"status":2}', '{"quality":NaN}', '[]', '{}'):
            with self.subTest(payload=payload):
                malformed.write_text(payload, encoding="utf-8")
                report = self.command("diagnose", malformed, "--json", expected=2)
                self.assertEqual(report["status"], "failed")
        incomplete = self.mesh.report
        incomplete["audit"]["complete"] = False
        malformed.write_text(json.dumps(incomplete, allow_nan=False), encoding="utf-8")
        report = self.command("diagnose", malformed, "--json", expected=1)
        self.assertEqual(report["status"], "incomplete")

    def test_saved_recipe_run_creates_source_checked_atomic_bundle(self):
        output = self.directory / "job"
        report = self.command("run", self.recipe, "--source", self.fixture,
                              "--output", output, "--json")
        self.assertEqual(report["schema"], "paacs-mesh-run-v1")
        self.assertEqual(report["status"], "succeeded")
        self.assertEqual((output / "source.step").read_bytes(), self.fixture.read_bytes())
        for name in ("mesh.paacsmesh", "mesh.vtk", "recipe.json", "intent.json", "report.json"):
            self.assertTrue((output / name).is_file(), name)
        saved_report = json.loads((output / "report.json").read_text(encoding="utf-8"))
        self.assertEqual(report["diagnostics"], pm.diagnose(saved_report))
        before = (output / "report.json").read_bytes()
        refused = self.command("run", self.recipe, "--source", self.fixture,
                               "--output", output, "--json", expected=1)
        self.assertEqual(refused["reason"], "output_exists")
        self.assertEqual(before, (output / "report.json").read_bytes())

    def test_saved_recipe_mismatch_never_publishes(self):
        mismatched = self.directory / "mismatched-recipe.json"
        recipe = json.loads(self.recipe.read_text(encoding="utf-8"))
        recipe["source"]["sha256"] = "0" * 64
        mismatched.write_text(json.dumps(recipe), encoding="utf-8")
        output = self.directory / "refused-job"
        report = self.command("run", mismatched, "--source", self.fixture,
                              "--output", output, "--json", expected=1)
        self.assertEqual(report["reason"], "source_mismatch")
        self.assertFalse(output.exists())

    def test_legacy_direct_meshing_keeps_native_report_contract(self):
        output = self.directory / "legacy.paacsmesh"
        report = self.command(self.fixture, output, "--size", .7)
        self.assertTrue(report["mesh_published"])
        self.assertTrue(report["audit"]["valid"])
        self.assertEqual(report["source_sha256"], self.model.source_sha256)
        self.assertTrue(output.is_file())

    def test_native_and_python_stdout_are_routed_then_restored(self):
        code = """
import ctypes
import json
import os
from paacs_mesh.__main__ import _native_output_to_stderr
with _native_output_to_stderr():
    print('python progress')
    os.write(1, b'native progress\\n')
    runtime = ctypes.CDLL('ucrtbase' if os.name == 'nt' else None)
    runtime.puts(b'C buffered progress')
print(json.dumps({'restored': True}))
"""
        child = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True,
                               text=True, timeout=30)
        self.assertEqual(child.returncode, 0, child.stderr)
        self.assertEqual(json.loads(child.stdout), {"restored": True})
        for message in ("python progress", "native progress", "C buffered progress"):
            self.assertIn(message, child.stderr)


if __name__ == "__main__":
    unittest.main()
