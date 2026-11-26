"""Explicit Native validation mode preserves mesh payloads and saved requests."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import paacs_mesh as pm
from paacs_mesh import __main__ as cli
from test_native import archive_data, fixture_path


class RayValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = pm.read_step(fixture_path().with_name("assembly.step"))

    def test_explicit_mode_keeps_payload_and_final_original_body_accounting(self):
        ordinary = self.model.mesh(.3, trimmed_ray_witnesses=False)
        fast = self.model.mesh(.3, trimmed_ray_witnesses=True)
        with tempfile.TemporaryDirectory() as directory:
            before, after = (Path(directory) / name for name in ("before.paacsmesh", "after.paacsmesh"))
            ordinary.write_archive(before)
            fast.write_archive(after)
            self.assertEqual(before.read_bytes(), after.read_bytes())
            _, _, cell_tags = archive_data(after)
        self.assertFalse(ordinary.report["requested_trimmed_ray_witnesses"])
        self.assertTrue(fast.report["requested_trimmed_ray_witnesses"])
        self.assertEqual(fast.report["effective_policy"]["auto_dispatch_scope"], "per_original_body")
        rays = 0
        for body in fast.report["bodies"]:
            population = sum(dim == 3 and tag == body["solid_tag"] for dim, tag in cell_tags)
            self.assertEqual(body["native_centroids_checked"] + body["native_final_ray_admissions"] +
                             body["native_centroid_exact_reuses"], population)
            self.assertTrue(body["native_centroid_source_guard"])
            rays += body["native_final_ray_admissions"]
        self.assertGreater(rays, 0)
        self.assertTrue(all(body["native_final_ray_admissions"] == 0 for body in ordinary.report["bodies"]))
        self.assertTrue(fast.report["audit"]["valid"])
        self.assertTrue(fast.report["quality"]["valid"])

    def test_recipe_roundtrip_and_historical_absence(self):
        recipe = pm.create_recipe(self.model, .5, trimmed_ray_witnesses=True, max_seconds=30)
        self.assertTrue(recipe["options"]["trimmed_ray_witnesses"])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "recipe.json"
            pm.save_recipe(self.model, recipe, path)
            loaded = pm.load_recipe(self.model, path)
            result = pm.run_recipe(self.model, loaded)
        self.assertTrue(result.mesh.report["requested_trimmed_ray_witnesses"])
        self.assertTrue(result.report["attempts"][0]["options"]["trimmed_ray_witnesses"])
        del recipe["options"]["trimmed_ray_witnesses"]
        historical = pm.run_recipe(self.model, recipe)
        self.assertFalse(historical.mesh.report["requested_trimmed_ray_witnesses"])

    def test_bool_admission_and_cli_mode_discovery(self):
        for value in (None, 0, 1, "true"):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    self.model.mesh(.5, trimmed_ray_witnesses=value)
                with self.assertRaises(TypeError):
                    self.model._native.mesh(.5, trimmed_ray_witnesses=value)
                with self.assertRaises(pm.WorkflowError):
                    pm.create_recipe(self.model, .5, trimmed_ray_witnesses=value)
        with tempfile.TemporaryDirectory() as directory:
            for flags, expected in (([], False), (["--trimmed-ray-witnesses"], True),
                                    (["--no-trimmed-ray-witnesses"], False)):
                output = io.StringIO()
                destination = Path(directory) / (str(expected) + str(len(flags)) + ".vtk")
                with contextlib.redirect_stdout(output):
                    code = cli.main([str(fixture_path().with_name("assembly.step")),
                                     str(destination), "--size", ".5", *flags])
                self.assertEqual(code, 0)
                self.assertTrue(destination.is_file())
                self.assertEqual(json.loads(output.getvalue())["requested_trimmed_ray_witnesses"], expected)


if __name__ == "__main__":
    unittest.main()
