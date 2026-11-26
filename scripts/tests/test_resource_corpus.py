"""Corpus resource policies forward explicitly without weakening qualification."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("resource_corpus_test_module", ROOT / "scripts/run_cad_corpus.py")
CORPUS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CORPUS)


class ResourceCorpusTests(unittest.TestCase):
    def request(self, extra=(), reported_resources=None):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            inputs, output = directory / "inputs", directory / "output"
            inputs.mkdir()
            (inputs / "fixture.step").write_bytes(b"Synthetic source for driver forwarding; never imported")
            executable = directory / "unused-cli.py"
            executable.write_text("# Synthetic executable identity; never executed\n")
            commands = []

            def refusal(command, *args):
                commands.append(command)
                return {"outcome": "resource-limit", "wall_seconds": .25,
                        "peak_process_working_set_bytes": None, "output_sha256": None,
                        "result": {"schema": "paacs-native-cad-report-v1", "status": "refused",
                                   "reason": "budget-exceeded", "mesh_published": False,
                                   "archive_published": False, "resources": reported_resources}}

            with patch.object(CORPUS, "run_case", side_effect=refusal), \
                 patch.object(CORPUS, "WindowsMemory", return_value=object()), \
                 patch.object(CORPUS, "suppress_native_error_dialogs"), \
                 contextlib.redirect_stdout(io.StringIO()):
                result = CORPUS.main([str(executable), str(inputs), "--output", str(output),
                                      "--repeats", "1", *extra])
            return result, commands, json.loads((output / "manifest.json").read_text()), \
                [json.loads(line) for line in (output / "results.jsonl").read_text().splitlines()]

    def test_non_native_defaults_keep_historical_work_flag_and_no_policy_fields(self):
        code, commands, manifest, rows = self.request()
        self.assertEqual(code, 0)
        command = commands[0]
        index = command.index("--max-work")
        self.assertEqual(command[index + 1], "2000000")
        self.assertNotIn("--resource-policy", command)
        self.assertNotIn("--recovery-mode", command)
        self.assertNotIn("--elements", command)
        self.assertNotIn("resource_policy", manifest["settings"])
        self.assertNotIn("recovery_mode", manifest["settings"])
        self.assertNotIn("elements", manifest["settings"])
        for name in CORPUS.OPTIONAL_NATIVE_FLAGS:
            self.assertNotIn(name, manifest["settings"])
        self.assertNotIn("resources_reported", rows[0])

    def test_inherited_normal_work_and_refusal_receipt_are_preserved(self):
        effective = {"policy": "normal", "requested_max_work": None,
                     "effective_max_work": 2**64 - 1, "logical_work_cap_disabled": True,
                     "max_working_bytes": 512 * 1024**2}
        code, commands, manifest, rows = self.request(
            ["--native-cad", "--resource-policy", "normal", "--max-work", "none"], effective)
        self.assertEqual(code, 0)
        self.assertNotIn("--max-work", commands[0])
        self.assertEqual(commands[0][commands[0].index("--resource-policy") + 1], "normal")
        self.assertIsNone(manifest["settings"]["max_work"])
        self.assertEqual(manifest["settings"]["resource_policy"], "normal")
        self.assertEqual(rows[0]["resource_policy_requested"], "normal")
        self.assertIsNone(rows[0]["max_work_requested"])
        self.assertEqual(rows[0]["resources_reported"], effective)
        self.assertEqual(rows[0]["outcome"], "resource-limit")
        self.assertFalse(rows[0]["result"]["mesh_published"])
        self.assertNotIn("--allow-centroids-within-chord-tolerance", commands[0])

    def test_explicit_work_cap_still_wins_and_no_effective_limits_are_invented(self):
        _, commands, manifest, rows = self.request(
            ["--native-cad", "--resource-policy", "normal", "--max-work", "300000000",
             "--max-working-bytes", "536870912"])
        command = commands[0]
        self.assertEqual(command[command.index("--max-work") + 1], "300000000")
        self.assertEqual(command[command.index("--max-working-bytes") + 1], "536870912")
        self.assertEqual(manifest["settings"]["max_work"], 300_000_000)
        self.assertIsNone(rows[0]["resources_reported"])

    def test_legacy_native_request_omits_new_flags_and_keeps_settings_identity(self):
        _, commands, manifest, rows = self.request(["--native-cad"])
        self.assertNotIn("--resource-policy", commands[0])
        self.assertNotIn("resource_policy", manifest["settings"])
        self.assertIsNone(rows[0]["resource_policy_requested"])
        self.assertIsNone(rows[0]["resources_reported"])

    def test_pinned_native_profile_flags_forward_only_when_selected(self):
        selected = ["--native-cad", "--recovery-mode", "advanced", "--elements", "tet4",
                    "--local-partition-recovery", "--trace-native-stages", "--adaptive-tet-recovery",
                    "--threads", "1", "--native-validation-threads", "1", "--jobs", "1"]
        _, commands, manifest, _ = self.request(selected)
        command = commands[0]
        for flag in ("--local-partition-recovery", "--trace-native-stages", "--adaptive-tet-recovery"):
            self.assertEqual(command.count(flag), 1)
        self.assertEqual(command[command.index("--recovery-mode") + 1], "advanced")
        self.assertEqual(command[command.index("--elements") + 1], "tet4")
        self.assertEqual(command[command.index("--threads") + 1], "1")
        self.assertEqual(command[command.index("--native-validation-threads") + 1], "1")
        self.assertNotIn("--jobs", command)
        self.assertEqual(manifest["settings"]["jobs"], 1)
        self.assertTrue(manifest["settings"]["adaptive_tet_recovery"])

    def test_non_native_policy_or_recovery_controls_are_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "must-not-exist"
            for flags in (["--resource-policy", "normal"], ["--recovery-mode", "advanced"],
                          ["--elements", "tet4"], ["--adaptive-tet-recovery"],
                          ["--local-partition-recovery"], ["--trace-native-stages"]):
                with self.subTest(flags=flags), contextlib.redirect_stderr(io.StringIO()), \
                     self.assertRaises(SystemExit) as caught:
                    CORPUS.main(["missing-cli", "missing-corpus", "--output", str(output), *flags])
                self.assertEqual(caught.exception.code, 2)
                self.assertFalse(output.exists())

    def test_resource_caps_are_typed_positive_or_inherited(self):
        self.assertIsNone(CORPUS.resource_cap("none"))
        self.assertEqual(CORPUS.resource_cap("300000000"), 300_000_000)
        for value in ("0", "-1", "1.5", "True", str(2**64)):
            with self.subTest(value=value), self.assertRaises(argparse.ArgumentTypeError):
                CORPUS.resource_cap(value)
        self.assertEqual(CORPUS.work_arguments(argparse.Namespace(
            max_work=2_000_000, resource_policy=None, native_cad=False)), ["--max-work", "2000000"])
        with self.assertRaises(ValueError):
            CORPUS.work_arguments(argparse.Namespace(max_work=None, resource_policy="normal", native_cad=False))

    def test_native_and_legacy_report_self_tests_still_pass_without_an_sdk(self):
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(CORPUS.main(["--self-test"]), 0)


if __name__ == "__main__":
    unittest.main()
