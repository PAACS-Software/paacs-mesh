"""Explicit warning and strict quality thresholds keep native validity checks."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import json
import math
from pathlib import Path
import tempfile
import unittest

import paacs_mesh as pm
from test_native import fixture_path


class NativeQualityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = pm.read_step(fixture_path())
        cls.good = cls.model.mesh(.3, max_work=100_000_000)
        cls.controlled_options = dict(
            face_sizes=[pm.FaceSize(1, .16)],
            volume_controls=[pm.BoxSize((-.3, -.3, .7), (.3, .3, 1.3), .15)],
            min_size=.12, growth=1.55, native_sizing=True,
            assembly="conforming", max_work=100_000_000,
            workers=1, validation_workers=1,
        )
        cls.controlled = cls.model.mesh(.7, **cls.controlled_options)
        cls.warned = cls.model.mesh(.7, warning_quality=.25, **cls.controlled_options)

    def assert_valid(self, mesh):
        report = mesh.report
        self.assertTrue(report["mesh_published"])
        self.assertEqual(report["reason"], "none")
        self.assertTrue(report["audit"]["complete"])
        self.assertTrue(report["audit"]["valid"])
        self.assertTrue(report["quality"]["complete"])
        self.assertTrue(report["quality"]["valid"])
        self.assertEqual(report["quality"]["inverted_cells"], 0)
        self.assertGreater(report["quality"]["minimum"], 0)
        json.dumps(report, allow_nan=False)
        return report

    def archive(self, mesh, directory, name):
        path = Path(directory) / name
        mesh.write_archive(path)
        return path.read_bytes()

    def test_default_warning_policy_is_explicitly_disabled(self):
        for mesh in (self.good, self.controlled):
            with self.subTest(controlled=mesh is self.controlled):
                report = self.assert_valid(mesh)
                quality = report["quality"]
                self.assertEqual(quality["warning_threshold"], 0)
                self.assertIs(quality["warning_enabled"], False)
                self.assertEqual(quality["minimum_required"], 0)
                self.assertIsNone(quality["requested_warning_threshold"])
                self.assertIsNone(quality["requested_minimum_required"])
                self.assertIs(report["quality_warning"], False)
        self.assertGreaterEqual(self.good.report["quality"]["minimum"], .25)

    def test_warning_policy_does_not_change_the_controlled_mesh(self):
        report = self.assert_valid(self.warned)
        quality = report["quality"]
        self.assertEqual(quality["warning_threshold"], .25)
        self.assertIs(quality["warning_enabled"], True)
        self.assertEqual(quality["minimum_required"], 0)
        self.assertEqual(quality["requested_warning_threshold"], .25)
        self.assertIsNone(quality["requested_minimum_required"])
        # Future mesher improvements may remove this fixture's slivers. The
        # warning policy must follow the measured result, without locking a
        # particular low-quality mesh into the regression suite.
        self.assertIs(report["quality_warning"], quality["minimum"] < .25)
        self.assertEqual(self.warned.node_count, self.controlled.node_count)
        self.assertEqual(self.warned.cell_count, self.controlled.cell_count)
        self.assertEqual(quality["minimum"], self.controlled.report["quality"]["minimum"])
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(self.archive(self.controlled, directory, "default.paacsmesh"),
                             self.archive(self.warned, directory, "warning.paacsmesh"))

    def test_zero_thresholds_preserve_default_archive_and_requested_context(self):
        mesh = self.model.mesh(.3, warning_quality=0, minimum_quality=0, max_work=100_000_000)
        quality = self.assert_valid(mesh)["quality"]
        self.assertIs(quality["warning_enabled"], False)
        self.assertEqual(quality["requested_warning_threshold"], 0)
        self.assertEqual(quality["requested_minimum_required"], 0)
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(self.archive(self.good, directory, "default.paacsmesh"),
                             self.archive(mesh, directory, "explicit-zero.paacsmesh"))

    def test_a_high_warning_threshold_never_becomes_a_strict_requirement(self):
        mesh = self.model.mesh(.3, warning_quality=1, max_work=100_000_000)
        report = self.assert_valid(mesh)
        self.assertEqual(report["quality"]["warning_threshold"], 1)
        self.assertIs(report["quality_warning"], True)
        self.assertEqual(report["quality"]["minimum_required"], 0)
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(self.archive(self.good, directory, "default.paacsmesh"),
                             self.archive(mesh, directory, "warning.paacsmesh"))

    def test_good_cylinder_meets_an_explicit_requirement_without_warning(self):
        mesh = self.model.mesh(.3, warning_quality=.25, minimum_quality=.25, max_work=100_000_000)
        report = self.assert_valid(mesh)
        self.assertGreaterEqual(report["quality"]["minimum"], .25)
        self.assertIs(report["quality_warning"], False)
        self.assertIs(report["quality"]["warning_enabled"], True)
        self.assertEqual(report["quality"]["minimum_required"], .25)

    def test_strict_refusal_preserves_an_existing_valid_mesh(self):
        with tempfile.TemporaryDirectory() as directory:
            before = self.archive(self.controlled, directory, "held.paacsmesh")
            with self.assertRaises(pm.MeshingError) as caught:
                self.model.mesh(.3, minimum_quality=1, warning_quality=.25, max_work=100_000_000)
            self.assertEqual(caught.exception.reason, "quality_requirement_not_met")
            report = caught.exception.report
            self.assertFalse(report["mesh_published"])
            quality = report["quality"]
            self.assertEqual(quality["warning_threshold"], .25)
            self.assertEqual(quality["minimum_required"], 1)
            self.assertEqual(quality["requested_warning_threshold"], .25)
            self.assertEqual(quality["requested_minimum_required"], 1)
            json.dumps(report, allow_nan=False)
            self.assertEqual(self.archive(self.controlled, directory, "after-refusal.paacsmesh"), before)

    def test_thresholds_reject_nonfinite_out_of_range_and_boolean_scalars(self):
        for name in ("warning_quality", "minimum_quality"):
            for value in (-.001, 1.001, math.nan, math.inf, -math.inf, True, False, ".25", [], {}):
                with self.subTest(name=name, value=value), self.assertRaises((ValueError, TypeError, OverflowError)):
                    self.model.mesh(.3, **{name: value})

    def test_report_threshold_context_is_owned_and_json_safe(self):
        report = self.warned.report
        report["quality"]["warning_threshold"] = 0
        report["quality"]["requested_warning_threshold"] = 0
        report["quality"]["minimum_required"] = 1
        report["quality"]["requested_minimum_required"] = 1
        fresh = self.warned.report
        self.assertEqual(fresh["quality"]["warning_threshold"], .25)
        self.assertEqual(fresh["quality"]["requested_warning_threshold"], .25)
        self.assertEqual(fresh["quality"]["minimum_required"], 0)
        self.assertIsNone(fresh["quality"]["requested_minimum_required"])
        json.dumps(fresh, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
