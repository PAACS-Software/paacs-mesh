"""Owned scaling plans preserve identities and do not invent policy gains."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "benchmarks/resource_policy_scaling.py"
SPEC = importlib.util.spec_from_file_location("resource_scaling_test_module", SCRIPT)
SCALING = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SCALING)


def attempt(policy, repetition, **changes):
    row = {"case_id": "box-coarse", "profile": policy, "repetition": repetition,
           "wall_size_mm": None, "global_size_mm": .5, "global_max_work": 300_000_000,
           "global_max_working_bytes": 512 * 1024**2, "status": "succeeded",
           "mesh_sha256": "a" * 64, "source_sha256": "b" * 64,
           "native_extension_sha256": "c" * 64, "native_seconds": float(repetition),
           "process_seconds": float(repetition + 1), "nodes": 48, "tetrahedra": 100,
           "quality": {"mean_ratio": {"minimum": .2}}, "achieved_sizing": {},
           "initial_front_attempted": False}
    row.update(changes)
    return row


class ResourceScalingTests(unittest.TestCase):
    def test_default_plan_has_six_cases_two_policies_and_three_repetitions(self):
        rows = SCALING.plan(argparse.Namespace(cases=list(SCALING.CASES),
            policies=list(SCALING.POLICIES), repetitions=3))
        self.assertEqual(len(rows), 36)
        self.assertEqual(len({row["directory"] for row in rows}), 36)
        self.assertEqual({row["global_size_mm"] for row in rows if row["case_id"].startswith("box")}, {.5, .25})
        self.assertEqual({row["global_size_mm"] for row in rows if row["case_id"].startswith("sphere")}, {.5, .25})
        self.assertTrue(all(len(row["source_sha256"]) == 64 for row in rows))

    def test_tag_checks_compare_source_live_and_archive_identities(self):
        self.assertTrue(SCALING.verify_tags([2, 1], [7], [1, 2], [7], "test")["matches"])
        for faces, bodies in (([1], [7]), ([1, 2], [8]), ([1, 1, 2], [7])):
            with self.subTest(faces=faces, bodies=bodies), self.assertRaises(ValueError):
                SCALING.verify_tags([1, 2], [7], faces, bodies, "test")

    def test_reported_resources_must_honor_policy_and_identical_explicit_caps(self):
        resources = {"policy": "normal", "requested_max_work": 300_000_000,
                     "effective_max_work": 300_000_000, "max_working_bytes": 512 * 1024**2,
                     "logical_work_cap_disabled": False}
        self.assertEqual(SCALING.verify_resources({"resources": resources}, "normal", 300_000_000,
                         512 * 1024**2), resources)
        for change in ({"policy": "bounded"}, {"requested_max_work": None},
                       {"effective_max_work": 50_000_000}, {"max_working_bytes": 1024**3},
                       {"logical_work_cap_disabled": True}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                SCALING.verify_resources({"resources": dict(resources, **change)}, "normal",
                                         300_000_000, 512 * 1024**2)

    def test_inactive_identical_routes_are_preservation_evidence(self):
        rows = [attempt(policy, repetition) for policy in ("bounded", "normal") for repetition in (1, 2, 3)]
        summary = SCALING.summarize(rows, 3)
        pair = summary["paired_cases"][0]
        self.assertTrue(pair["complete_valid_comparison"])
        self.assertTrue(pair["initial_front_inactive_for_both"])
        self.assertTrue(pair["meshes_identical_across_policies"])
        self.assertFalse(pair["quality_or_simulation_improvement_inferred"])
        self.assertEqual(summary["repeated_groups"][0]["native_mesh_median_seconds"], 2)
        self.assertIsNone(SCALING.initial_front_activity({"bodies": []}))
        self.assertFalse(SCALING.initial_front_activity({"bodies": [{"initial_front": {"diagnostic": "not_attempted"}}]}))

    def test_failed_or_changed_identity_attempt_cannot_form_a_valid_pair(self):
        original = [attempt(policy, repetition) for policy in ("bounded", "normal") for repetition in (1, 2, 3)]
        for changed in ({"status": "failed"}, {"source_sha256": "d" * 64},
                        {"native_extension_sha256": "d" * 64}, {"mesh_sha256": "d" * 64}):
            rows = copy.deepcopy(original)
            rows[-1].update(changed)
            with self.subTest(changed=changed):
                self.assertFalse(SCALING.summarize(rows, 3)["paired_cases"][0]["complete_valid_comparison"])

    def test_cli_plan_only_requires_no_native_and_none_caps_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "plan"
            child = subprocess.run([sys.executable, "-I", str(SCRIPT), "--output", str(output), "--plan-only"],
                capture_output=True, text=True, timeout=30)
            self.assertEqual(child.returncode, 0, child.stderr)
            self.assertEqual(json.loads(child.stdout)["native_calls"], 0)
            self.assertEqual(len(json.loads((output / "plan.json").read_text())["cases"]), 36)
        with self.assertRaises(argparse.ArgumentTypeError):
            SCALING.explicit_cap("none")


if __name__ == "__main__":
    unittest.main()
