"""Python discovery and guarded controls retain native requests and validity."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import copy
from dataclasses import replace
import inspect
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import paacs_mesh as pm
from paacs_mesh import cad

from test_native import fixture_path


class ApiUsabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = pm.read_step(fixture_path())
        cls.wall = cls.model.select_faces(surface="cylinder", radius=1, expect=1)
        cls.baseline = cls.model.mesh(.7)
        cls.temporary = tempfile.TemporaryDirectory(prefix="paacs-api-usability-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.directory = Path(cls.temporary.name)

    def archive(self, mesh, name):
        path = self.directory / name
        mesh.write_archive(path)
        return path.read_bytes()

    def test_mesh_has_a_real_python_signature_and_native_defaults(self):
        signature = inspect.signature(pm.CadModel.mesh)
        expected = dict(recovery="standard", assembly="preserve_parts", workers=1,
            validation_workers=1, max_nodes=None, max_cells=None, max_work=None,
            max_output_bytes=None, max_working_bytes=None, max_audit_bytes=None,
            max_seconds=None, face_sizes=None, selected_face_sizes=(),
            volume_controls=None, min_size=None, growth=1.3, native_sizing=None,
            warning_quality=None, minimum_quality=None, advancing_interior_front=None,
            construction_effort="legacy", max_surface_points=None,
            max_surface_chord_points=None, max_interior_candidates=None, max_interior_sites=None)
        self.assertEqual(set(signature.parameters), {"self", "size", *expected})
        for name, default in expected.items():
            self.assertEqual(signature.parameters[name].default, default)
            self.assertEqual(signature.parameters[name].kind, inspect.Parameter.KEYWORD_ONLY)
        self.assertIn("FaceSelection", inspect.getdoc(pm.CadModel.mesh))
        self.assertIn("in mm", inspect.getdoc(pm.CadModel.mesh))
        with self.assertRaises(TypeError):
            self.model.mesh(.7, unknown_option=True)

    def test_default_wrapper_keeps_the_native_constructor_output(self):
        native = self.model._native.mesh(.7)
        self.assertEqual(self.archive(self.baseline, "wrapper.paacsmesh"),
                         self.archive(native, "native.paacsmesh"))

    def test_direct_selection_explicit_face_and_recipe_produce_the_same_mesh(self):
        options = dict(assembly="conforming", growth=1.4, min_size=.18, max_work=100_000_000)
        selected = [(self.wall, .3, 1.2)]
        direct = self.model.mesh(.7, selected_face_sizes=selected, **options)
        explicit = self.model.mesh(.7, face_sizes=[pm.FaceSize(self.wall.tags[0], .3, growth=1.2)], **options)
        result = pm.run_recipe(self.model, pm.create_recipe(self.model, .7,
                              selected_face_sizes=selected, **options))
        expected = self.archive(explicit, "explicit.paacsmesh")
        self.assertEqual(self.archive(direct, "selected.paacsmesh"), expected)
        self.assertEqual(self.archive(result.mesh, "recipe.paacsmesh"), expected)
        for mesh in (direct, explicit, result.mesh):
            self.assertEqual(mesh.face_tags, list(self.model.face_tags))
            self.assertEqual(mesh.body_tags, list(self.model.body_tags))
            self.assertTrue(mesh.report["audit"]["complete"] and mesh.report["audit"]["valid"])
            self.assertTrue(mesh.report["quality"]["complete"] and mesh.report["quality"]["valid"])
            self.assertEqual(mesh.report["source_sha256"], self.model.source_sha256)
            self.assertEqual(mesh.report["sizing"]["face_sizes"][0]["growth"], 1.2)

    def test_none_means_no_optional_recipe_controls(self):
        expected = pm.create_recipe(self.model, .7)
        actual = pm.create_recipe(self.model, .7, face_sizes=None,
                                  selected_face_sizes=None, volume_controls=None)
        self.assertEqual(actual, expected)

    def test_direct_selection_keeps_the_requested_assembly_policy(self):
        with self.assertRaises(pm.MeshingError) as caught:
            self.model.mesh(.7, selected_face_sizes=[(self.wall, .3)])
        self.assertEqual(caught.exception.reason, "unsupported_control")
        self.assertEqual(caught.exception.report["requested_assembly"], "preserve_parts")
        self.assertFalse(caught.exception.report["mesh_published"])

    def test_foreign_and_changed_body_selections_are_rejected_consistently(self):
        foreign = pm.read_step(fixture_path()).select_faces(surface="cylinder", expect=1)
        for selected, reason in ((foreign, "foreign_selection"),
                                 (replace(self.wall, _body=999), "unknown_body")):
            with self.subTest(reason=reason):
                with self.assertRaises(pm.SelectionError) as direct:
                    self.model.mesh(.7, selected_face_sizes=[(selected, .3)])
                with self.assertRaises(pm.WorkflowError) as recipe:
                    pm.create_recipe(self.model, .7, selected_face_sizes=[(selected, .3)])
                self.assertEqual(direct.exception.reason, reason)
                self.assertEqual(recipe.exception.reason, reason)

    def test_property_guards_are_checked_for_direct_selection_sizing(self):
        guards = self.wall.guards
        guards[0]["area"] *= 2
        changed = replace(self.wall, _guards_json=json.dumps(guards))
        with self.assertRaises(pm.SelectionError) as caught:
            self.model.mesh(.7, selected_face_sizes=[(changed, .3)])
        self.assertEqual(caught.exception.reason, "guard_mismatch")

    def test_selected_control_errors_name_the_entry_and_keep_empty_selection_validation(self):
        empty = self.model.select_faces(surface="sphere", expect=(0, 0))
        for selection in (self.wall, empty):
            with self.subTest(empty=len(selection) == 0), self.assertRaises(pm.SelectionError) as caught:
                self.model.mesh(.7, selected_face_sizes=[(selection, 0)])
            self.assertEqual(caught.exception.stage, "sizing")
            self.assertEqual(caught.exception.reason, "invalid_argument")
            self.assertEqual(caught.exception.report["control_index"], 0)
            self.assertIn("selected_face_sizes[0]", str(caught.exception))
            diagnostic = pm.diagnose(caught.exception)
            issue = next(item for item in diagnostic["issues"] if item["code"] == "control_request_failed")
            self.assertEqual(issue["evidence"]["details"]["control_index"], 0)
            self.assertFalse(any(item["code"] == "mesh_refused" for item in diagnostic["issues"]))

    def test_selection_expansion_has_a_bound_before_native_meshing(self):
        with patch.object(cad, "_MAX_FACE_CONTROLS", 0):
            with self.assertRaises(pm.SelectionError) as caught:
                self.model.mesh(.7, selected_face_sizes=[(self.wall, .3)])
        self.assertEqual(caught.exception.reason, "control_budget")
        self.assertEqual(caught.exception.report["requested_controls"], 1)

    def test_existing_direct_face_control_refusals_still_come_from_native_validation(self):
        with self.assertRaises(pm.MeshingError) as caught:
            self.model.mesh(.7, assembly="conforming",
                            face_sizes=[pm.FaceSize(1, .3), pm.FaceSize(1, .4)])
        self.assertEqual(caught.exception.reason, "invalid_request")
        with self.assertRaises(pm.MeshingError) as combined:
            self.model.mesh(.7, assembly="conforming", face_sizes=[pm.FaceSize(1, .3)],
                            selected_face_sizes=[(self.wall, .4)])
        self.assertEqual(combined.exception.reason, "invalid_request")

    def test_recipe_errors_identify_the_duplicate_face_and_invalid_box(self):
        with self.assertRaises(pm.WorkflowError) as caught:
            pm.create_recipe(self.model, .7, face_sizes=[pm.FaceSize(1, .3)],
                             selected_face_sizes=[(self.wall, .4)])
        self.assertIn("face tag 1", str(caught.exception))
        recipe = pm.create_recipe(self.model, .7,
                                  volume_controls=[pm.BoxSize((0, 0, 0), (1, 1, 1), .3)])
        recipe["options"]["volume_controls"][0]["upper"] = [0, 1, 1]
        with patch.object(pm.CadModel, "mesh") as native:
            with self.assertRaises(pm.WorkflowError) as box:
                pm.run_recipe(self.model, recipe)
            native.assert_not_called()
        self.assertIn("volume_controls[0]", str(box.exception))
        diagnostic = pm.diagnose(box.exception)
        self.assertEqual(diagnostic["status"], "failed")
        self.assertIn("volume_controls[0]", pm.format_diagnostics(diagnostic))
        self.assertFalse(any(item["code"] == "mesh_refused" for item in diagnostic["issues"]))
        for collection, controls in (
            ("face_sizes", {"face_sizes": [pm.FaceSize(1, .3)]}),
            ("volume_controls", {"volume_controls": [pm.BoxSize((0, 0, 0), (1, 1, 1), .3)]}),
        ):
            malformed = pm.create_recipe(self.model, .7, **controls)
            del malformed["options"][collection][0]["growth"]
            with self.subTest(collection=collection), self.assertRaises(pm.WorkflowError) as missing:
                pm.run_recipe(self.model, malformed)
            self.assertIn(f"{collection}[0]", str(missing.exception))

    def test_selection_errors_and_cli_envelopes_keep_actual_candidates(self):
        with self.assertRaises(pm.SelectionError) as caught:
            self.model.select_faces(surface="plane", expect=1)
        error = caught.exception
        original = copy.deepcopy(error.report)
        diagnostic = pm.diagnose(error)
        self.assertEqual(diagnostic["status"], "failed")
        issue = next(item for item in diagnostic["issues"] if item["code"] == "selection_failed")
        self.assertEqual(issue["evidence"]["reason"], "cardinality")
        self.assertEqual(issue["evidence"]["details"]["matched"], 2)
        self.assertEqual(issue["evidence"]["details"]["candidate_tags"], original["candidate_tags"])
        self.assertEqual(issue["cad_face_tags"], [])
        self.assertEqual(issue["cell_ids"], [])
        self.assertIsNone(diagnostic["source_sha256"])
        self.assertIn("Request status: failed", pm.format_diagnostics(diagnostic))
        self.assertFalse(any(item["code"] == "mesh_refused" for item in diagnostic["issues"]))
        envelope = {"schema": "paacs-mesh-cli-error-v1", "status": "failed", "stage": error.stage,
                    "reason": error.reason, "message": str(error), "report": error.report}
        self.assertEqual(pm.diagnose(envelope), diagnostic)
        issue["evidence"]["details"]["candidate_tags"].append(999)
        self.assertEqual(error.report, original)

    def test_native_hint_is_preserved_without_changing_the_native_report(self):
        with self.assertRaises(pm.MeshingError) as caught:
            self.model.mesh(.7, face_sizes=[pm.FaceSize(1, .3)])
        error = caught.exception
        original = copy.deepcopy(error.report)
        diagnostic = pm.diagnose(error)
        self.assertEqual(diagnostic["summary"]["native_hint"], original["hint"])
        refused = next(item for item in diagnostic["issues"] if item["code"] == "mesh_refused")
        self.assertEqual(refused["evidence"]["hint"], original["hint"])
        self.assertIn("Native hint: " + original["hint"], pm.format_diagnostics(diagnostic))
        self.assertEqual(error.report, original)
        self.assertEqual(diagnostic["status"], "failed")

    def test_diagnostics_do_not_create_hints_or_accept_unstructured_exceptions(self):
        self.assertIsNone(pm.diagnose(self.baseline)["summary"]["native_hint"])
        native = self.baseline.report
        native["hint"] = {"invented": "not a native string"}
        with self.assertRaises(ValueError):
            pm.diagnose(native)
        with self.assertRaises(ValueError):
            pm.diagnose(ValueError("an unrelated exception"))


if __name__ == "__main__":
    unittest.main()
