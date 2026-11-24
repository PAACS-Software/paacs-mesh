"""Pure planner and receipt checks; no native extension is imported."""

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
import tempfile
import types
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("initial_generation_plan_tests", ROOT / "benchmarks/qualify_initial_generation.py")
TRIALS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TRIALS)


def manifest():
    fixture = ROOT / "tests/fixtures/cad/thin_plate.step"
    sha = TRIALS.SHARED.sha256(fixture)
    common = {"kind": "owned", "fixture": "tests/fixtures/cad/thin_plate.step", "source_sha256": sha, "size": .1}
    return {"cases": [
        {**common, "case_id": "thin-plate-010", "baseline_action": "reuse_exact_artifact_and_settings"},
        {**common, "case_id": "thin-plate-020", "baseline_action": "pending_baseline_run"},
        {"case_id": "original51-00000003", "kind": "downloaded_reference_gap"},
        {"case_id": "holdout-example", "kind": "downloaded_holdout"}],
        "limits": {"max_nodes": 500000, "max_cells": 3000000, "max_work": None,
                   "max_working_bytes": 1024**3, "max_output_bytes": 256*1024**2, "max_audit_bytes": 512*1024**2},
        "growth": 1.3}


def settings(path, **changes):
    result = dict(case_manifest=path, cases=["thin-plate-010"], repetitions=1, native_seconds=600., hard_seconds=720.,
                  size_factor=None, max_chord_error=None, max_relative_volume_error=None,
                  curvature_points_per_turn=None, initial_topology_quality=None, advancing_interior_front=None,
                  include_thin_forced_pair=False, development_native=ROOT / "unused.pyd", occt_sdk=ROOT)
    result.update(changes)
    return argparse.Namespace(**result)


class InitialGenerationPlanTests(unittest.TestCase):
    def test_public_owned_manifest_is_fresh_reproducible_and_contains_no_local_corpus(self):
        path = ROOT / "benchmarks/manifests/initial-generation-owned.json"
        data = json.loads(path.read_text())
        self.assertEqual(len(data["cases"]), 14)
        self.assertNotIn("corpus_root", data)
        self.assertEqual(data["limits"]["max_seconds"], 600.)
        self.assertEqual(data["hard_process_seconds"], 720.)
        self.assertEqual((data["limits"]["workers"], data["limits"]["validation_workers"]), (1,1))
        self.assertEqual(data["limits"]["max_working_bytes"], 1024**3)
        for case in data["cases"]:
            self.assertEqual(case["kind"], "owned")
            self.assertEqual(case["baseline_action"], "fresh_run_no_frozen_baseline")
            self.assertNotIn("baseline_record", case)
            self.assertNotIn("baseline_mesh_sha256", case)
            self.assertFalse(Path(case["fixture"]).is_absolute())
        with patch.object(TRIALS.HELPERS, "development_identity", return_value={"test": "identity"}), \
             patch.object(TRIALS.SHARED, "load_development", side_effect=AssertionError("No native load")):
            plan = TRIALS.build_plan(settings(path, cases=None))
        self.assertEqual(len(plan["cases"]), 14)
        self.assertEqual(len({row["directory"] for row in plan["cases"]}), 14)

    def test_material_bounds_use_named_coordinates_and_union_all_bodies(self):
        first = {"bounds": {"lower": [-1., -2., 0.], "upper": [3., 1., 2.]}}
        self.assertEqual(TRIALS.material_bounds([first]), ([-1., -2., 0.], [3., 1., 2.], 29**.5))
        second = {"bounds": {"lower": [1., -4., -1.], "upper": [5., 2., 3.]}}
        self.assertEqual(TRIALS.material_bounds([first, second]), ([-1., -4., -1.], [5., 2., 3.], 88**.5))
        for boxes in ([], [{"bounds": {"lower": [0,0,0], "upper": [0,0,0]}}],
                      [{"bounds": {"lower": [0,0,0], "upper": [-1,1,1]}}],
                      [{"bounds": {"lower": [0,0,0], "upper": [float("inf"),1,1]}}]):
            with self.subTest(boxes=boxes), self.assertRaises(ValueError):
                TRIALS.material_bounds(boxes)

    def test_default_selection_keeps_gap_cases_optional_and_prioritizes_pairs(self):
        data = manifest()
        selected = TRIALS.select_cases(data)
        self.assertEqual([row["case_id"] for row in selected], ["thin-plate-010", "thin-plate-020", "holdout-example"])
        self.assertEqual(TRIALS.select_cases(data, ["original51-00000003"])[0]["case_id"], "original51-00000003")
        for invalid in ([], ["missing"], ["thin-plate-010", "thin-plate-010"]):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                TRIALS.select_cases(data, invalid)
        data["cases"].append(copy.deepcopy(data["cases"][0]))
        with self.assertRaises(ValueError):
            TRIALS.select_cases(data)

    def test_candidate_default_kwargs_keep_native_defaults_and_explicit_limits(self):
        data = manifest()
        args = settings(Path("not-read"))
        options = TRIALS.mesh_options(args, data)
        self.assertEqual((options["workers"], options["validation_workers"], options["max_seconds"]), (1, 1, 600.))
        self.assertEqual(options["max_working_bytes"], 1024**3)
        self.assertIsNone(options["max_work"])
        self.assertEqual(options["resource_policy"], "normal")
        self.assertEqual(options["construction_effort"], "legacy")
        self.assertTrue(all(options[key] is None for key in TRIALS.OPTION_KEYS))
        self.assertNotIn("minimum_quality", options)
        self.assertNotIn("warning_quality", options)
        args.curvature_points_per_turn = 32
        args.max_chord_error = .01
        self.assertEqual(TRIALS.mesh_options(args, data)["curvature_points_per_turn"], 32)
        self.assertEqual(TRIALS.mesh_options(args, data)["max_chord_error"], .01)

    def test_forced_pair_is_derived_without_mutating_the_pinned_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            data = manifest(); encoded = json.dumps(data)
            path.write_text(encoded)
            args = settings(path, include_thin_forced_pair=True, repetitions=2)
            with patch.object(TRIALS.HELPERS, "development_identity", return_value={"test": "identity"}), \
                 patch.object(TRIALS.SHARED, "load_development", side_effect=AssertionError("No native load")):
                plan = TRIALS.build_plan(args)
            self.assertEqual(path.read_text(), encoded)
            self.assertEqual(len(plan["cases"]), 6)
            self.assertEqual(len({row["directory"] for row in plan["cases"]}), 6)
            off, on = plan["cases"][2], plan["cases"][4]
            self.assertEqual(off["case_id"], "thin-plate-010-front-quality-off")
            self.assertEqual(on["case_id"], "thin-plate-010-front-quality-on")
            self.assertEqual(off["source_sha256"], plan["cases"][0]["source_sha256"])
            self.assertEqual(on["size"], .1)
            self.assertTrue(TRIALS.case_options(plan, off)["advancing_interior_front"])
            self.assertIs(TRIALS.case_options(plan, off)["initial_topology_quality"], False)
            self.assertIs(TRIALS.case_options(plan, on)["initial_topology_quality"], True)
            self.assertIsNone(plan["options"]["advancing_interior_front"])
            self.assertIsNone(plan["options"]["initial_topology_quality"])

    def test_changed_source_and_excessive_plan_refuse_before_identity_or_native(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            data = manifest(); data["cases"][0]["source_sha256"] = "0"*64
            path.write_text(json.dumps(data))
            with patch.object(TRIALS.HELPERS, "development_identity", side_effect=AssertionError("Too late")):
                with self.assertRaisesRegex(ValueError, "Pinned source changed"):
                    TRIALS.build_plan(settings(path))
            path.write_text(json.dumps(manifest()))
            with self.assertRaisesRegex(ValueError, "attempts"):
                TRIALS.build_plan(settings(path, repetitions=TRIALS.MAX_ATTEMPTS+1))
            with self.assertRaisesRegex(ValueError, "only to imported"):
                TRIALS.build_plan(settings(path, size_factor=.075))

    def test_source_paths_and_variant_overrides_cannot_escape_the_plan(self):
        with self.assertRaises(ValueError):
            TRIALS.relative_source(ROOT, "../outside.step")
        with self.assertRaises(ValueError):
            TRIALS.relative_source(ROOT, str((ROOT / "tests/fixtures/cad/box.step").resolve()))
        for overrides in ({"max_work": None}, {"initial_topology_quality": 1}):
            with self.subTest(overrides=overrides), self.assertRaises(ValueError):
                TRIALS.case_options({"options": {}}, {"option_overrides": overrides})

    def test_normal_receipt_requires_actual_requested_cap_and_mode(self):
        options = TRIALS.mesh_options(settings(Path("not-read")), manifest())
        good = {"resources": {"policy": "normal", "requested_max_work": None,
                              "effective_max_work": 2**62-1, "logical_work_cap_disabled": True,
                              "max_working_bytes": 1024**3}}
        TRIALS.check_resources(good, options)
        for key, value in (("policy", "bounded"), ("logical_work_cap_disabled", False),
                           ("max_working_bytes", 512*1024**2), ("effective_max_work", 0)):
            changed = copy.deepcopy(good); changed["resources"][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                TRIALS.check_resources(changed, options)
        explicit = {**options, "max_work": 300000000}
        good["resources"].update(requested_max_work=300000000, effective_max_work=300000000, logical_work_cap_disabled=False)
        TRIALS.check_resources(good, explicit)
        good["resources"]["effective_max_work"] -= 1
        with self.assertRaises(ValueError):
            TRIALS.check_resources(good, explicit)

    @unittest.skipUnless(importlib.util.find_spec("numpy"), "Optional NumPy is absent")
    def test_final_interior_nodes_are_not_reported_as_initial_front_sites(self):
        import numpy as np
        snapshot = {"points": np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1], [.2,.2,.2], [99,99,99.]]),
                    "cells": {"tri3": {"connectivity": np.array([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]),
                                       "entity_tags": np.array([1,2,3,4])},
                              "tet4": {"connectivity": np.array([[0,1,2,4], [0,1,3,4], [0,2,3,4], [1,2,3,4]])}}}
        class Mesh:
            def to_numpy(self, *, max_bytes):
                self.max_bytes = max_bytes
                return snapshot
        mesh = Mesh()
        report = {"bodies": [{"solid_tag": 1, "initial_front": {"nodes": 0, "diagnostic": "stalled"}}]}
        first = TRIALS.node_receipt(mesh, report)
        self.assertEqual(mesh.max_bytes, 256*1024**2)
        self.assertEqual((first["final_material_nodes"], first["final_boundary_nodes"], first["final_interior_nodes"]), (5,4,1))
        self.assertEqual(first["nodes_unused_by_tet4"], 1)
        self.assertEqual(first["initial_front_reported_added_sites"][0]["reported_added_sites"], 0)
        self.assertIsNone(first["initial_admitted_points_snapshot"])
        snapshot["points"][4,0] = .21
        second = TRIALS.node_receipt(mesh, report)
        self.assertEqual(first["final_boundary_ordered_points_sha256"], second["final_boundary_ordered_points_sha256"])
        self.assertNotEqual(first["final_interior_ordered_points_sha256"], second["final_interior_ordered_points_sha256"])

    def test_failed_attempt_costs_and_native_status_are_retained(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            attempt = output / "thin-repeat-001"; attempt.mkdir()
            (attempt / "case.json").write_text(json.dumps({"case_id": "thin", "status": "native_refused",
                "native_status": "refused", "call_seconds": 600.03, "native_report_file": "refusal-report.json"}))
            identity = {"native_extension_sha256": "a"*64, "source_snapshot_sha256": "b"*64}
            plan = {"cases": [{"case_id": "thin", "directory": attempt.name, "repetition": 1}],
                    "development_identity": identity}
            args = settings(Path(directory) / "unused-manifest.json", output=output)
            process = types.SimpleNamespace(pid=123, wait=lambda **kwargs: 1)
            with patch.object(TRIALS.HELPERS, "development_identity", return_value=identity), \
                 patch.object(TRIALS.subprocess, "Popen", return_value=process) as launch, contextlib.redirect_stdout(io.StringIO()):
                code = TRIALS.execute(args, plan)
            self.assertEqual(code, 1)
            rows = json.loads((output / "results.json").read_text())
            self.assertEqual(rows[0]["status"], "native_refused")
            self.assertEqual(rows[0]["call_seconds"], 600.03)
            self.assertEqual(rows[0]["native_report_file"], "refusal-report.json")
            self.assertFalse(rows[0]["process_timeout"])
            self.assertGreaterEqual(rows[0]["process_seconds"], 0)
            launch.assert_called_once()

    def test_invalid_fraction_or_curvature_is_rejected(self):
        self.assertEqual(TRIALS.curvature("32"), 32)
        self.assertEqual(TRIALS.fraction(".01"), .01)
        for value in ("-1", "65", "1.5", "True"):
            with self.subTest(value=value), self.assertRaises(argparse.ArgumentTypeError):
                TRIALS.curvature(value)
        for value in ("0", "1", "nan", "inf"):
            with self.subTest(value=value), self.assertRaises(argparse.ArgumentTypeError):
                TRIALS.fraction(value)


if __name__ == "__main__":
    unittest.main()
