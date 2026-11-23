"""Resource experiment plans retain refusals and require comparable repetitions."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
import contextlib
import copy
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "benchmarks/construction_budget_trials.py"
SPEC = importlib.util.spec_from_file_location("resource_trial_test_module", SCRIPT)
TRIALS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TRIALS)


def settings(**changes):
    values = dict(work_values=None, working_byte_values=None, max_work=100_000_000,
                  max_working_bytes=1024**3, repetitions=1, wall_sizes=[.1, .08],
                  profiles=["legacy", "quadruple"], max_cases=256,
                  source=ROOT / "tests/fixtures/cad/cylinder.step", size=.5,
                  max_seconds=90, hard_seconds=150, development_native=None,
                  occt_sdk=None, plan_only=False)
    values.update(changes)
    return argparse.Namespace(**values)


def measured_repetition(number, **changes):
    values = dict(profile="quadruple", wall_size_mm=.1, global_size_mm=.5,
                  global_max_work=100_000_000, global_max_working_bytes=1024**3,
                  repetition=number, status="succeeded", mesh_sha256="a" * 64,
                  source_sha256="b" * 64, native_extension_sha256="c" * 64,
                  options=TRIALS.PROFILES["quadruple"], native_seconds=float(number),
                  process_seconds=float(number + 1))
    values.update(changes)
    return values


class ResourceTrialPlanTests(unittest.TestCase):
    def test_baseline_profiles_and_case_names_survive_a_new_public_default(self):
        for profile in ("legacy", "explicit-legacy", "double", "quadruple", "fair-legacy", "auto", "thorough"):
            self.assertEqual(TRIALS.PROFILES[profile]["resource_policy"], "bounded")
        self.assertEqual(TRIALS.PROFILES["normal"]["resource_policy"], "normal")
        self.assertEqual(TRIALS.PROFILES["normal"]["construction_effort"], "legacy")
        self.assertNotIn("max_surface_points", TRIALS.PROFILES["normal"])
        old = TRIALS.trial_plan(settings())
        self.assertEqual(len(old), 4)
        self.assertEqual(old[0]["name"], "legacy-wall-0p1")
        self.assertEqual(old[-1]["name"], "quadruple-wall-0p08")
        quadruple = TRIALS.PROFILES["quadruple"]
        normal = TRIALS.PROFILES["normal-quadruple"]
        self.assertEqual({key: value for key, value in normal.items() if key != "resource_policy"},
                         {key: value for key, value in quadruple.items() if key != "resource_policy"})

    def test_independent_grid_has_unique_bounded_repetition_paths(self):
        args = settings(work_values=[50_000_000, 100_000_000, 500_000_000],
                        working_byte_values=[512 * 1024**2, 1024**3], repetitions=3)
        rows = TRIALS.trial_plan(args)
        self.assertEqual(len(rows), 72)
        self.assertEqual(len({row["name"] for row in rows}), 72)
        self.assertEqual({(row["max_work"], row["max_working_bytes"]) for row in rows},
                         {(work, memory) for work in args.work_values for memory in args.working_byte_values})
        args.max_cases = 71
        with self.assertRaises(ValueError):
            TRIALS.trial_plan(args)

    def test_resource_caps_preserve_none_and_reject_invalid_native_values(self):
        self.assertIsNone(TRIALS.resource_cap("none"))
        self.assertEqual(TRIALS.resource_cap("50000000"), 50_000_000)
        for value in ("0", "-1", "1.5", "True", str(2**63)):
            with self.subTest(value=value), self.assertRaises(argparse.ArgumentTypeError):
                TRIALS.resource_cap(value)

    def test_cli_plan_only_does_not_load_native_and_records_inherited_limits(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "plan"
            child = subprocess.run([sys.executable, "-I", str(SCRIPT), "--output", str(output),
                "--profiles", "normal", "--wall-sizes", ".1", "--max-work", "none",
                "--max-working-bytes", "none", "--repetitions", "3", "--plan-only"],
                capture_output=True, text=True, timeout=30)
            self.assertEqual(child.returncode, 0, child.stderr)
            self.assertEqual(json.loads(child.stdout)["native_calls"], 0)
            plan = json.loads((output / "trial-plan.json").read_text())
            self.assertEqual(len(plan["cases"]), 3)
            self.assertTrue(all(row["max_work"] is None and row["max_working_bytes"] is None
                                for row in plan["cases"]))
            self.assertEqual(list(output.iterdir()), [output / "trial-plan.json"])

    def test_duplicate_cli_values_refuse_before_creating_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "plan"
            child = subprocess.run([sys.executable, "-I", str(SCRIPT), "--output", str(output),
                "--profiles", "quadruple", "--work-values", "50000000", "50000000", "--plan-only"],
                capture_output=True, text=True, timeout=30)
            self.assertEqual(child.returncode, 2)
            self.assertFalse(output.exists())


class ResourceTrialEvidenceTests(unittest.TestCase):
    def test_native_refusal_keeps_reported_work_and_elapsed_call_time(self):
        report = {"status": "refused", "reason": "budget_exceeded", "mesh_published": False,
                  "seconds": {"generation": .25, "meshing_total": .3},
                  "resources": {"policy": "bounded", "effective_max_work": 50_000_000},
                  "bodies": [{"solid_tag": 1, "native_surface_work": 17,
                              "initial_front": {"work": 123, "budget_reason": "work_limit"}}]}

        class Refusal(RuntimeError):
            stage, reason = "meshing", "budget_exceeded"

        error = Refusal("synthetic refusal for evidence-accounting test")
        error.report = report
        model = types.SimpleNamespace(select_faces=lambda **kwargs: types.SimpleNamespace(tags=[1]))
        model.mesh = lambda *args, **kwargs: (_ for _ in ()).throw(error)
        package = types.ModuleType("paacs_mesh")
        package.__version__ = "synthetic-test-only"
        package._native = types.SimpleNamespace(__file__=str(SCRIPT))
        package.read_step = lambda source: model
        package.FaceSize = lambda tag, size: (tag, size)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "case"
            args = settings(output=output, profile="quadruple", wall_size=.1, repetition=1)
            with patch.dict(sys.modules, {"paacs_mesh": package}), \
                 patch.object(TRIALS, "peak_process_bytes", return_value=64), \
                 patch.object(TRIALS.time, "perf_counter", side_effect=iter(range(20))), \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(TRIALS.worker(args), 1)
            saved = json.loads((output / "case.json").read_text())
            self.assertEqual(saved["call_seconds"], 1)
            self.assertEqual(saved["native_status"], "refused")
            self.assertEqual(saved["native_reported_seconds"], report["seconds"])
            self.assertEqual(saved["native_stage_counters"][0]["initial_front"]["work"], 123)
            self.assertIsNone(saved["native_reported_validation_work"])
            self.assertEqual(json.loads((output / "refusal-report.json").read_text()), report)

    def test_medians_require_complete_identical_valid_repetitions(self):
        rows = [measured_repetition(i) for i in (1, 2, 3)]
        median = TRIALS.repeated_summary(rows, 3)[0]
        self.assertEqual(median["native_mesh_median_seconds"], 2)
        self.assertEqual(median["attempt_process_median_seconds"], 3)
        self.assertIsNone(TRIALS.repeated_summary(rows[:2], 2)[0]["native_mesh_median_seconds"])
        for changed in ({"status": "failed"}, {"mesh_sha256": "d" * 64},
                        {"native_extension_sha256": "d" * 64}, {"source_sha256": None},
                        {"global_size_mm": .25}, {"repetition": 1}):
            altered = copy.deepcopy(rows)
            altered[1].update(changed)
            with self.subTest(changed=changed):
                self.assertIsNone(TRIALS.repeated_summary(altered, 3)[0]["native_mesh_median_seconds"])
        rows[1]["status"] = "failed"
        summary = TRIALS.repeated_summary(rows, 3)[0]
        self.assertEqual(summary["failed_attempts"], 1)
        self.assertEqual(summary["all_attempt_process_seconds"], [2, 3, 4])

    def test_controller_retains_partial_receipt_and_failed_process_after_result(self):
        calls = []

        def simulated_worker(command, **kwargs):
            self.assertEqual(command[command.index("--max-seconds") + 1], "180")
            self.assertEqual(command[command.index("--hard-seconds") + 1], "240")
            self.assertEqual(kwargs["timeout"], 240)
            case = Path(command[command.index("--output") + 1])
            case.mkdir()
            calls.append(case)
            if len(calls) == 1:
                (case / "case.json").write_text("{")
            else:
                (case / "case.json").write_text(json.dumps({"status": "succeeded", "native_status": "succeeded"}))
            return types.SimpleNamespace(returncode=1)

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "controller"
            args = settings(output=output, profiles=["quadruple"], wall_sizes=[.1], repetitions=3,
                            max_seconds=180, hard_seconds=240)
            with patch.object(TRIALS.subprocess, "run", side_effect=simulated_worker), \
                 patch.object(TRIALS.time, "perf_counter", side_effect=iter(range(20))), \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(TRIALS.controller(args), 1)
            rows = json.loads((output / "summary.json").read_text())
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[0]["reason"], "worker_receipt_unreadable")
            self.assertEqual(rows[1]["status"], "process_failed_after_result")
            self.assertEqual(rows[1]["worker_reported_status"], "succeeded")
            self.assertEqual((calls[0] / "case.json").read_text(), "{")
            summary = json.loads((output / "repetitions.json").read_text())[0]
            self.assertEqual(summary["failed_attempts"], 3)
            self.assertEqual(len(summary["all_attempt_process_seconds"]), 3)
            self.assertIsNone(summary["native_mesh_median_seconds"])


if __name__ == "__main__":
    unittest.main()
