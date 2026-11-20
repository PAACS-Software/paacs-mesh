#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Process-boundary regressions for the standalone CAD corpus evidence runner."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

RUNNER = Path(__file__).resolve().parents[1] / "scripts" / "run_cad_corpus.py"
spec = importlib.util.spec_from_file_location("cad_runner", RUNNER)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)

CHILD = r'''
import json, os, sys, time
from pathlib import Path
name = Path(sys.argv[1]).stem
output = Path(sys.argv[3])
payload = dict(schema="paacs-cad-report-v1", status="refused", stage="import",
               reason="unsupported-surface", mesh_published=False, archive_published=False,
               import_seconds=.01, preparation_seconds=0, generation_seconds=0)
if os.name == "nt":
    import ctypes
    from ctypes import wintypes
    get_mode = ctypes.WinDLL("kernel32").GetErrorMode
    get_mode.argtypes = []
    get_mode.restype = wintypes.UINT
    payload["inherited_error_mode"] = get_mode()
code = 3
if name == "success":
    # Force a visible working-set measurement and deterministic archive content.
    memory = bytearray(4 * 1024 * 1024)
    time.sleep(.03)
    output.write_bytes(b"PAACSM8\0" + b"fake child publication")
    payload.update(status="complete", stage="archive", reason="none",
                   mesh_published=True, archive_published=True,
                   preparation_seconds=.02, generation_seconds=.03)
    code = 0
elif name == "malformed":
    payload["reason"] = "parse-failed"
elif name == "timeout":
    time.sleep(5)
elif name == "crash":
    os._exit(17)
elif name == "partial":
    output.write_bytes(b"partial archive")
    payload.update(status="failed", stage="archive", reason="archive-write-failed",
                   mesh_published=True)
    code = 6
elif name == "bad-success":
    payload.update(status="complete", archive_published=True)
    code = 0
elif name == "missing-report":
    sys.exit(3)
print("OCCT transfer diagnostic before machine-readable report")
print(json.dumps(payload), flush=True)
sys.exit(code)
'''


class CorpusRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="paacs-cad-runner-test-")
        self.directory = Path(self.temporary.name)
        self.corpus = self.directory / "corpus"
        self.corpus.mkdir()
        self.executable = self.directory / "child.py"
        self.executable.write_text(CHILD, encoding="utf-8")
        self.output = self.directory / "evidence"

    def tearDown(self):
        self.temporary.cleanup()

    def invoke(self, *extra):
        return subprocess.run([sys.executable, str(RUNNER), str(self.executable), str(self.corpus),
                               "--output", str(self.output), "--limit", "30", "--smallest", "0",
                               "--repeats", "2", "--timeout", "2", *extra],
                              capture_output=True, text=True, timeout=30)

    def test_process_outcomes_memory_and_resume(self):
        for name in ("success", "unsupported", "malformed", "timeout", "crash",
                     "partial", "bad-success", "missing-report"):
            (self.corpus / (name + ".step")).write_text("ISO-10303-21;", encoding="utf-8")
        first = self.invoke()
        self.assertEqual(first.returncode, 1, first.stderr)
        rows_path = self.output / "results.jsonl"
        previous = rows_path.read_bytes()
        rows = [json.loads(line) for line in previous.decode().splitlines()]
        self.assertEqual(len(rows), 16)
        expected = {"success": "succeeded", "unsupported": "unsupported",
                    "malformed": "malformed-input", "timeout": "timeout", "crash": "crash",
                    "partial": "publication-violation", "bad-success": "invalid-result",
                    "missing-report": "missing-report"}
        for row in rows:
            self.assertEqual(row["outcome"], expected[Path(row["path"]).stem])
            self.assertFalse(row["process_memory_limit_enforced"])
            self.assertTrue((self.output / row["run_directory"] / "stdout.log").is_file())
        successes = [row for row in rows if row["outcome"] == "succeeded"]
        self.assertEqual(len({row["output_sha256"] for row in successes}), 1)
        if os.name == "nt":
            self.assertTrue(all(row["peak_process_working_set_bytes"] > 4 * 1024 * 1024
                                for row in successes))
            self.assertTrue(all(row["result"]["inherited_error_mode"] & 0x8003 == 0x8003
                                for row in successes))
        report = json.loads((self.output / "report.json").read_text())
        success_summary = next(row for row in report["cases"] if row["path"] == "success.step")
        self.assertTrue(success_summary["archive_hash_stable"])
        self.assertEqual(success_summary["timings_seconds"]["generation_seconds"]["median"], .03)
        self.assertEqual(self.invoke("--resume").returncode, 1)
        self.assertEqual(rows_path.read_bytes(), previous)
        # A truncated final append is recovered without replaying completed jobs.
        with rows_path.open("ab") as stream:
            stream.write(b'{"interrupted":')
        self.assertEqual(self.invoke("--resume").returncode, 1)
        self.assertEqual(rows_path.read_bytes(), previous)
        (self.corpus / "success.step").write_text("changed", encoding="utf-8")
        changed = self.invoke("--resume")
        self.assertEqual(changed.returncode, 2)
        self.assertIn("input changed", changed.stderr)
        self.assertEqual(rows_path.read_bytes(), previous)

    def test_deterministic_selection_and_input_safety(self):
        for i in range(7):
            (self.corpus / f"{i}.STEP").write_bytes(b"x" * (8 - i))
        count, first = runner.select_files(self.corpus, 2, 1, 3, [])
        _, second = runner.select_files(self.corpus, 2, 1, 3, [])
        self.assertEqual(count, 7)
        self.assertEqual(first, second)
        self.assertEqual([row["path"] for row in first], ["0.STEP", "1.STEP", "3.STEP", "6.STEP"])
        with self.assertRaises(ValueError):
            runner.select_files(self.corpus, 0, 0, 0, ["../child.py"])
        rejected = self.invoke("--output", str(self.corpus / "evidence"))
        self.assertEqual(rejected.returncode, 2)
        self.assertFalse((self.corpus / "evidence").exists())

    def test_report_parser_rejects_nonfinite_numbers(self):
        self.assertIsNone(runner.last_report('{"schema":"paacs-cad-report-v1","n":NaN}'))
        self.assertEqual(runner.last_report('noise\n{"schema":"paacs-cad-report-v1"}')
                         ["schema"], "paacs-cad-report-v1")

    def test_budget_and_requirement_taxonomy(self):
        expected = {"file-budget": "resource-limit", "meshing_budget": "resource-limit",
                    "allocation-failed": "allocation-failed",
                    "quality_requirement_not_met": "requirement-refusal",
                    "size_requirement_not_met": "requirement-refusal",
                    "cancelled": "cancelled", "kernel-failure": "failed"}
        for reason, outcome in expected.items():
            self.assertEqual(runner.classify(3, False, {"reason": reason}, False, False), outcome)
        self.assertEqual(runner.classify(-1073741819, False, None, False, False), "crash")
        self.assertEqual(runner.classify(0, False, {"status": "complete",
                                                  "archive_published": True,
                                                  "mesh_published": False}, True, True),
                         "invalid-result")


if __name__ == "__main__":
    unittest.main()
