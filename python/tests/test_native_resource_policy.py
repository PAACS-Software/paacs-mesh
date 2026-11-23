"""Strict native resource policy admission, precedence and old-recipe parity."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import json
from pathlib import Path
import tempfile
import unittest

import paacs_mesh as pm
from paacs_mesh import workflow as wf
from paacs_mesh.cad import _CONSTRUCTION_FIELDS

from test_native import fixture_path


class NativeResourcePolicyTests(unittest.TestCase):
    """An older binding must fail these checks instead of skipping new semantics."""

    @classmethod
    def setUpClass(cls):
        cls.model = pm.read_step(fixture_path())

    def assert_valid(self, mesh):
        report = mesh.report
        self.assertTrue(report["mesh_published"])
        for name in ("audit", "quality"):
            self.assertTrue(report[name]["complete"])
            self.assertTrue(report[name]["valid"])
        self.assertEqual(report["quality"]["inverted_cells"], 0)
        self.assertEqual(mesh.body_tags, list(self.model.body_tags))
        self.assertEqual(mesh.face_tags, list(self.model.face_tags))
        json.dumps(report, allow_nan=False)

    def test_omitted_work_uses_policy_while_memory_remains_bounded(self):
        for policy in ("bounded", "normal"):
            with self.subTest(policy=policy):
                mesh = self.model.mesh(.7, resource_policy=policy, workers=1,
                                       validation_workers=1, max_seconds=90)
                self.assert_valid(mesh)
                resources = mesh.report["resources"]
                self.assertEqual(resources["policy"], policy)
                self.assertIsNone(resources["requested_max_work"])
                self.assertEqual(resources["max_working_bytes"], 512 * 1024**2)
                self.assertEqual(resources["logical_work_cap_disabled"], policy == "normal")
                if policy == "bounded":
                    self.assertEqual(resources["effective_max_work"], 50_000_000)
                else:
                    self.assertGreaterEqual(resources["effective_max_work"], 2**60)

    def test_explicit_work_overrides_normal_and_tiny_allowance_refuses(self):
        mesh = self.model.mesh(.7, resource_policy="normal", max_work=100_000_000,
                               max_working_bytes=1024**3, workers=1, validation_workers=1,
                               max_seconds=90)
        self.assert_valid(mesh)
        resources = mesh.report["resources"]
        self.assertEqual(resources["requested_max_work"], 100_000_000)
        self.assertEqual(resources["effective_max_work"], 100_000_000)
        self.assertFalse(resources["logical_work_cap_disabled"])
        self.assertEqual(resources["max_working_bytes"], 1024**3)
        with self.assertRaises(pm.MeshingError) as caught:
            self.model.mesh(.7, resource_policy="normal", max_work=1,
                            workers=1, validation_workers=1, max_seconds=90)
        error = caught.exception
        self.assertEqual(error.reason, "budget_exceeded")
        self.assertFalse(error.report["mesh_published"])
        self.assertEqual(error.report["resources"]["effective_max_work"], 1)
        self.assertFalse(error.report["resources"]["logical_work_cap_disabled"])
        self.assert_valid(mesh)

    def test_old_recipe_absence_replays_the_bounded_legacy_mesh_exactly(self):
        recipe = wf.create_recipe(self.model, .7, resource_policy="normal",
                                  construction_effort="auto", max_seconds=90)
        del recipe["options"]["resource_policy"]
        for name in _CONSTRUCTION_FIELDS:
            del recipe["options"][name]
        expected = self.model.mesh(.7, resource_policy="bounded", construction_effort="legacy",
                                   workers=1, validation_workers=1, max_seconds=90)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "old-recipe.json"
            path.write_text(json.dumps(recipe), encoding="utf-8")
            loaded = wf.load_recipe(self.model, path)
            self.assertEqual(loaded["options"]["resource_policy"], "bounded")
            self.assertEqual(loaded["options"]["construction_effort"], "legacy")
            replayed = wf.run_recipe(self.model, loaded).mesh
            self.assert_valid(replayed)
            expected.write_archive(root / "expected.paacsmesh")
            replayed.write_archive(root / "replayed.paacsmesh")
            self.assertEqual((root / "expected.paacsmesh").read_bytes(),
                             (root / "replayed.paacsmesh").read_bytes())
            self.assertEqual(replayed.report["resources"], expected.report["resources"])

    def test_normal_resource_refusals_preserve_existing_output_and_source(self):
        fixture = fixture_path()
        source_bytes = fixture.read_bytes()
        retained = self.model.mesh(.7, resource_policy="bounded", workers=1,
                                   validation_workers=1, max_seconds=90)
        self.assert_valid(retained)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "retained.paacsmesh"
            retained.write_archive(output)
            expected = output.read_bytes()
            requests = (({"max_work": 1}, "budget_exceeded"),
                        ({"max_working_bytes": 1024}, "budget_exceeded"),
                        ({"max_nodes": 5}, "budget_exceeded"),
                        ({"max_cells": 8}, "budget_exceeded"),
                        ({"max_working_bytes": 1}, "budget_exceeded"))
            for limits, reason in requests:
                with self.subTest(limits=limits):
                    with self.assertRaises(pm.MeshingError) as caught:
                        candidate = self.model.mesh(.7, resource_policy="normal", workers=1,
                                                    validation_workers=1, max_seconds=90, **limits)
                        candidate.write_archive(output)
                    self.assertEqual(caught.exception.reason, reason)
                    self.assertFalse(caught.exception.report["mesh_published"])
                    self.assertEqual(caught.exception.report["resources"]["policy"], "normal")
                    if "max_work" in limits:
                        self.assertEqual(caught.exception.report["resources"]["effective_max_work"], 1)
                        self.assertFalse(caught.exception.report["resources"]["logical_work_cap_disabled"])
                    if "max_working_bytes" in limits:
                        self.assertEqual(caught.exception.report["resources"]["max_working_bytes"], limits["max_working_bytes"])
                    self.assertEqual(output.read_bytes(), expected)
                    self.assertEqual(fixture.read_bytes(), source_bytes)
            for call in (lambda: self.model.mesh(.7, resource_policy="unexpected"),
                         lambda: self.model._native.mesh(.7, resource_policy="unexpected")):
                with self.assertRaises(ValueError):
                    call()
                self.assertEqual(output.read_bytes(), expected)
                self.assertEqual(fixture.read_bytes(), source_bytes)
        self.assert_valid(retained)

    def test_normal_cooperative_cancellation_refuses_without_partial_publication(self):
        fixture = fixture_path()
        source_bytes = fixture.read_bytes()
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "retained.paacsmesh"
            sentinel = b"Existing output must survive a native cancellation"
            output.write_bytes(sentinel)
            with self.assertRaises(pm.MeshingError) as caught:
                candidate = self.model.mesh(.7, resource_policy="normal", workers=1,
                                            validation_workers=1, max_seconds=1e-9)
                candidate.write_archive(output)
            error = caught.exception
            self.assertEqual(error.reason, "cancelled")
            self.assertFalse(error.report["mesh_published"])
            self.assertEqual(error.report["resources"]["policy"], "normal")
            self.assertTrue(error.report["resources"]["logical_work_cap_disabled"])
            json.dumps(error.report, allow_nan=False)
            self.assertEqual(output.read_bytes(), sentinel)
            self.assertEqual(fixture.read_bytes(), source_bytes)


if __name__ == "__main__":
    unittest.main()
