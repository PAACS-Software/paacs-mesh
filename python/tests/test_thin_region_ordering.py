"""Opt-in phase-order controls and unchanged thick-part construction."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
from pathlib import Path
import tempfile
import unittest

import paacs_mesh as pm

ROOT = Path(__file__).resolve().parents[2]


class ThinRegionOrderingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = pm.read_step(ROOT / "tests/fixtures/cad/box.step")

    def test_capability_is_opt_in_and_bounded(self):
        capability = pm.capabilities()["controls"]["thin_region_ordering"]
        self.assertIs(capability["default"], False)
        self.assertIs(capability["experimental"], True)
        self.assertEqual(capability["source_samples"], 64)
        self.assertFalse(capability["forces_initial_front"])
        self.assertFalse(capability["changes_native_validity_gates"])

    def test_recipe_roundtrip_retains_independent_toggle(self):
        recipe = pm.create_recipe(self.model, .5, thin_region_ordering=True)
        self.assertTrue(recipe["options"]["thin_region_ordering"])
        self.assertFalse(recipe["options"]["boundary_first"])
        self.assertFalse(recipe["options"]["native_source_fallback"])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"recipe.json"
            pm.save_recipe(self.model, recipe, path)
            loaded = pm.load_recipe(self.model, path)
        self.assertEqual(loaded["options"], recipe["options"])

    def test_non_boolean_rejected_before_native_job(self):
        for value in (1, "auto", None, []):
            with self.subTest(value=value):
                with self.assertRaises((TypeError, ValueError)):
                    self.model.mesh(.5, thin_region_ordering=value)
                with self.assertRaises(pm.WorkflowError):
                    pm.create_recipe(self.model, .5, thin_region_ordering=value)

    def test_thick_front_keeps_mesh_and_native_validity(self):
        options = dict(assembly="conforming", advancing_interior_front=True,
            max_seconds=30, max_nodes=50000, max_cells=250000,
            max_working_bytes=1024**3, max_chord_error=.025)
        before = self.model.mesh(.5, **options)
        after = self.model.mesh(.5, thin_region_ordering=True, **options)
        self.assertFalse(before.report["requested_thin_region_ordering"])
        report = after.report
        self.assertTrue(report["requested_thin_region_ordering"])
        ordering = report["bodies"][0]["initial_front"]["ordering"]
        self.assertTrue(ordering["attempted"])
        self.assertTrue(ordering["complete"])
        self.assertFalse(ordering["selected_boundary_first"])
        for group in ("audit", "quality"):
            self.assertTrue(report[group]["complete"] and report[group]["valid"])
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory)/name for name in ("before.paacsmesh", "after.paacsmesh")]
            before.write_archive(paths[0])
            after.write_archive(paths[1])
            self.assertEqual(*(hashlib.sha256(path.read_bytes()).hexdigest() for path in paths))

    def test_face_sizing_abstains_without_changing_the_working_control(self):
        tag = self.model.face_tags[0]
        options = dict(assembly="conforming", advancing_interior_front=True,
            face_sizes=[pm.FaceSize(tag, .4)], workers=1, validation_workers=1,
            max_seconds=30, max_nodes=50000, max_cells=250000,
            max_working_bytes=1024**3, max_chord_error=.025)
        before = self.model.mesh(.5, **options)
        after = self.model.mesh(.5, thin_region_ordering=True, **options)
        for mesh in (before, after):
            report = mesh.report
            self.assertTrue(report["mesh_published"])
            self.assertEqual(report["source_sha256"], self.model.source_sha256)
            self.assertEqual(mesh.face_tags, list(self.model.face_tags))
            self.assertEqual(mesh.body_tags, list(self.model.body_tags))
            self.assertEqual(report["sizing"]["face_sizes"][0]["tag"], tag)
            self.assertEqual(report["sizing"]["face_sizes"][0]["size"], .4)
            for group in ("audit", "quality"):
                self.assertTrue(report[group]["complete"] and report[group]["valid"])
        self.assertFalse(before.report["requested_thin_region_ordering"])
        self.assertTrue(after.report["requested_thin_region_ordering"])
        ordering = after.report["bodies"][0]["initial_front"]["ordering"]
        self.assertFalse(ordering["attempted"])
        self.assertFalse(ordering["complete"])
        self.assertFalse(ordering["selected_boundary_first"])
        self.assertEqual(ordering["diagnostic"], "native_face_sizing_not_qualified")
        for counter in ("work", "samples", "distinct_samples", "blocked_samples",
                        "unknown_samples", "maximum_modeled_bytes"):
            self.assertEqual(ordering[counter], 0)
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory)/name for name in ("controlled-before.paacsmesh", "controlled-after.paacsmesh")]
            before.write_archive(paths[0])
            after.write_archive(paths[1])
            self.assertEqual(*(hashlib.sha256(path.read_bytes()).hexdigest() for path in paths))


if __name__ == "__main__":
    unittest.main()
