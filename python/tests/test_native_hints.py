"""Check effective construction and bounded sizing diagnostics in native reports."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import json
import math
import os
from pathlib import Path
import unittest

import paacs_mesh as pm


def fixture_path():
    override = os.environ.get("PAACS_TEST_STEP_FIXTURE")
    return Path(override) if override else Path(__file__).resolve().parents[2] / "tests/fixtures/cad/cylinder.step"


class NativeHintsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = pm.read_step(fixture_path())
        cls.wall = cls.model.select_faces(surface="cylinder", radius=1, expect=1)
        options = dict(assembly="conforming", workers=1, validation_workers=1,
                       max_seconds=90, max_work=100_000_000)
        cls.uniform = cls.model.mesh(.5, **options)
        cls.controlled = cls.model.mesh(.5, face_sizes=[pm.FaceSize(cls.wall.tags[0], .3)], **options)
        cls.disabled = cls.model.mesh(.5, advancing_interior_front=False, **options)
        cls.wrapper_disabled = cls.model.mesh(.5, native_sizing=False, **options)

    def assert_native_routes(self, report):
        policy = report["effective_policy"]
        self.assertEqual(policy["body_routes"], [{"solid_tag": body["solid_tag"], "route": body["route"],
                                                 "reason": body["reason"]} for body in report["bodies"]])
        self.assertEqual(policy["wrapper_native_sizing"], report["sizing"]["native_sizing"])
        self.assertEqual(policy["wrapper_sizing_mode"], report["sizing"]["mode"])
        self.assertEqual(policy["requested_native_sizing"], report["sizing"]["requested_native_sizing"])
        self.assertEqual(policy["requested_advancing_interior_front"], report["requested_advancing_interior_front"])
        json.dumps(report, allow_nan=False)

    def test_auto_route_differs_from_wrapper_policy(self):
        report = self.uniform.report
        self.assert_native_routes(report)
        policy = report["effective_policy"]
        self.assertIsNone(policy["requested_native_sizing"])
        self.assertIsNone(policy["requested_advancing_interior_front"])
        self.assertEqual(policy["wrapper_sizing_mode"], "standard")
        self.assertFalse(policy["wrapper_native_sizing"])
        self.assertEqual(policy["body_routes"][0]["route"], "native_spatial_sizing")
        self.assertGreater(report["bodies"][0]["initial_front"]["nodes"], 0)
        # A uniform request uses initial construction without adaptive sizing.
        interior = report["bodies"][0]["interior_front"]
        self.assertEqual(interior["diagnostic"], "not_attempted")
        self.assertEqual(interior["nodes"], 0)
        self.assertEqual(interior["remaining_oversized"], 0)

    def test_selected_control_exposes_adaptive_completion_diagnostic(self):
        report = self.controlled.report
        self.assert_native_routes(report)
        self.assertTrue(report["audit"]["valid"])
        self.assertTrue(report["quality"]["valid"])
        self.assertEqual(report["effective_policy"]["wrapper_sizing_mode"], "spatial")
        self.assertEqual(report["sizing"]["face_sizes"][0]["tag"], self.wall.tags[0])
        interior = report["bodies"][0]["interior_front"]
        self.assertIn(interior["diagnostic"], {"budget_limited", "scaffold_refused", "targets_satisfied", "targets_unmet"})
        for key in ("nodes", "candidates", "work", "queries", "remaining_oversized",
                    "front_activations", "void_activations", "front_advances"):
            self.assertIs(type(interior[key]), int)
            self.assertGreaterEqual(interior[key], 0)
        self.assertIs(type(interior["seconds"]), float)
        self.assertTrue(math.isfinite(interior["seconds"]))
        self.assertGreater(interior["seconds"], 0)
        # Completion is the native stage diagnostic; no new guarantee about
        # final per-face achieved sizes is fabricated by these report fields.
        if interior["diagnostic"] == "targets_satisfied":
            self.assertEqual(interior["remaining_oversized"], 0)

    def test_explicit_front_disable_is_preserved(self):
        report = self.disabled.report
        self.assert_native_routes(report)
        self.assertFalse(report["effective_policy"]["requested_advancing_interior_front"])
        self.assertEqual(report["bodies"][0]["initial_front"]["diagnostic"], "not_attempted")
        self.assertEqual(report["bodies"][0]["interior_front"]["diagnostic"], "not_attempted")
        self.assertTrue(report["audit"]["valid"])
        # native_sizing=False also disables AUTO in the wrapper, but that
        # internal resolution must not rewrite the caller's unset front flag.
        wrapper = self.wrapper_disabled.report
        self.assert_native_routes(wrapper)
        self.assertFalse(wrapper["effective_policy"]["requested_native_sizing"])
        self.assertIsNone(wrapper["effective_policy"]["requested_advancing_interior_front"])
        self.assertEqual(wrapper["bodies"][0]["initial_front"]["diagnostic"], "not_attempted")

    def test_refusal_retains_only_observed_routes(self):
        with self.assertRaises(pm.MeshingError) as caught:
            self.model.mesh(.5, assembly="conforming", max_nodes=1, max_work=100_000_000, max_seconds=90)
        report = caught.exception.report
        self.assertFalse(report["mesh_published"])
        self.assert_native_routes(report)
        self.assertNotEqual(report["reason"], "none")
        self.assertIsNone(report["effective_policy"]["requested_advancing_interior_front"])


if __name__ == "__main__":
    unittest.main()
