"""Source-provider requests and witnesses are distinct from mesh certification."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import contextlib
import copy
import inspect
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import paacs_mesh as pm
from paacs_mesh import cad, diagnostics, workflow as wf, __main__ as cli
from test_construction_budgets import model_stub, synthetic_refusal
from test_construction_diagnostics import synthetic_report
from test_native_hints import fixture_path


def fallback_shell():
    return {"available": True, "accepted": True, "diagnostic": "accepted",
        "native_fallback_requested": True, "native_fallback_attempted": True,
        "native_fallback_diagnostic": "native_source_accepted",
        "attempts": [
            {"attempt": 1, "provider": "regular_chart", "surface_published": False,
             "global_sizing_applied": True, "curvature_points_per_turn": 32,
             "angular_deflection": None, "surface_diagnostic": "nonconforming_edge",
             "witness": {"stage": "curve_sampling", "detail": "parameter_schedule_refused",
                         "edge_tag": 9, "face_tag": None, "use_id": 0, "geometric_refusal": True}},
            {"attempt": 2, "provider": "native_tessellation", "surface_published": True,
             "global_sizing_applied": False, "curvature_points_per_turn": 32,
             "angular_deflection": .19634954084936207, "volume_criterion_met": True,
             "witness": {"stage": "not_reported", "detail": "not_reported", "edge_tag": None,
                         "face_tag": None, "use_id": None, "geometric_refusal": False}}]}


def material_receipt(*, complete=False):
    return {"diagnostic": "none" if complete else "cancelled", "complete": complete,
            "source_model_verified_unchanged": complete, "attempts": 2,
            "requested_points": 100, "attempted_points": 23,
            "native_queries": 10, "strict_reuses": 7, "ray_admissions": 6,
            "work": 300, "seconds": .12,
            "profiled_batches": 1, "ray_witness_attempts": 8,
            "classifier_setup_seconds": .01, "native_query_seconds": .04, "ray_witness_seconds": .02,
            "scope": "cumulative material classification attempts; last batch status; not final publication"}


class SourceFallbackTests(unittest.TestCase):
    """Pure stub/metadata tests; these do not claim a Native mesh was produced."""

    def setUp(self):
        self.model, self.native = model_stub()

    def test_boolean_default_and_explicit_forwarding_do_not_force_a_front(self):
        for function in (cad.CadModel.mesh, wf.create_recipe):
            parameter = inspect.signature(function).parameters["native_source_fallback"]
            self.assertIs(parameter.default, False)
            self.assertEqual(parameter.kind, inspect.Parameter.KEYWORD_ONLY)
        self.model.mesh(.5)
        self.assertIs(self.native.mesh.call_args.kwargs["native_source_fallback"], False)
        for enabled in (False, True):
            for front in (None, False, True):
                self.model.mesh(.5, native_source_fallback=enabled, advancing_interior_front=front)
                options = self.native.mesh.call_args.kwargs
                self.assertIs(options["native_source_fallback"], enabled)
                self.assertIs(options["advancing_interior_front"], front)
        for invalid in (None, 0, 1, "true", [], {}):
            with self.subTest(value=invalid), self.assertRaises((TypeError, ValueError)):
                self.model.mesh(.5, native_source_fallback=invalid)
            with self.subTest(recipe=invalid), self.assertRaises(wf.WorkflowError):
                wf.create_recipe(self.model, .5, native_source_fallback=invalid)

    def test_recipe_save_reopen_and_old_absence_keep_explicit_meaning(self):
        recipe = wf.create_recipe(self.model, .5, native_source_fallback=True)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "recipe.json"
            wf.save_recipe(self.model, recipe, path)
            self.assertIs(wf.load_recipe(self.model, path)["options"]["native_source_fallback"], True)
            old = copy.deepcopy(recipe)
            del old["options"]["native_source_fallback"]
            path.write_text(json.dumps(old), encoding="utf-8")
            self.assertIs(wf.load_recipe(self.model, path)["options"]["native_source_fallback"], False)

    def test_explicit_retries_preserve_provider_choice_and_controls(self):
        success = self.native.mesh.return_value
        self.native.mesh.side_effect = [synthetic_refusal("native-plc-cavity-material-refusal"),
                                       synthetic_refusal("native-plc-final-centroid-outside"), success]
        recipe = wf.create_recipe(self.model, .5, native_source_fallback=True,
            max_chord_error=.0025, max_relative_volume_error=.005,
            curvature_points_per_turn=32, allow_advanced=True, finer_factors=[.5], max_attempts=3)
        result = wf.run_recipe(self.model, recipe)
        self.assertEqual([call.args[0] for call in self.native.mesh.call_args_list], [.5, .5, .25])
        for call, attempt in zip(self.native.mesh.call_args_list, result.report["attempts"]):
            for options in (call.kwargs, attempt["options"]):
                self.assertIs(options["native_source_fallback"], True)
                self.assertEqual(options["max_chord_error"], .0025)
                self.assertEqual(options["max_relative_volume_error"], .005)
                self.assertEqual(options["curvature_points_per_turn"], 32)

    def test_cli_default_and_last_explicit_switch_are_forwarded(self):
        for flags, enabled in (([], False), (["--native-source-fallback"], True),
                (["--native-source-fallback", "--no-native-source-fallback"], False)):
            with self.subTest(flags=flags), patch.object(cli, "read_step", return_value=self.model), \
                 patch.object(cli, "_native_output_to_stderr", contextlib.nullcontext), \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cli.main(["source.step", "out.vtk", "--size", ".5", *flags]), 0)
            self.assertIs(self.native.mesh.call_args.kwargs["native_source_fallback"], enabled)

    def test_discovery_is_experimental_and_keeps_validity_and_control_requirements(self):
        control = pm.capabilities()["controls"]["native_source_fallback"]
        self.assertIs(control["default"], False)
        self.assertIs(control["experimental"], True)
        for key in ("forces_initial_front", "changes_native_validity_gates",
                    "explicit_surface_controls_bypassed", "native_source_spacing_certified"):
            self.assertIs(control[key], False)
        self.assertEqual(control["source_attempt_limits"],
                         {"normal": 3, "bounded_enabled": 2, "bounded_default": 1})

    def test_boundary_first_is_independent_strict_and_does_not_force_a_constructor(self):
        for function in (cad.CadModel.mesh, wf.create_recipe):
            parameter = inspect.signature(function).parameters["boundary_first"]
            self.assertIs(parameter.default, False)
            self.assertEqual(parameter.kind, inspect.Parameter.KEYWORD_ONLY)
        for provider in (False, True):
            for order in (False, True):
                for front in (None, False, True):
                    self.model.mesh(.5, native_source_fallback=provider, boundary_first=order,
                                    advancing_interior_front=front)
                    options = self.native.mesh.call_args.kwargs
                    self.assertIs(options["boundary_first"], order)
                    self.assertIs(options["native_source_fallback"], provider)
                    self.assertIs(options["advancing_interior_front"], front)
        for invalid in (None, 0, 1, "true", [], {}):
            with self.subTest(value=invalid), self.assertRaises((TypeError, ValueError)):
                self.model.mesh(.5, boundary_first=invalid)
            with self.subTest(recipe=invalid), self.assertRaises(wf.WorkflowError):
                wf.create_recipe(self.model, .5, boundary_first=invalid)

    def test_boundary_first_recipes_keep_order_through_reopen_and_every_retry(self):
        recipe = wf.create_recipe(self.model, .5, boundary_first=True,
            native_source_fallback=True, advancing_interior_front=True,
            allow_advanced=True, finer_factors=[.5], max_attempts=3)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "recipe.json"
            wf.save_recipe(self.model, recipe, path)
            loaded = wf.load_recipe(self.model, path)
            self.assertIs(loaded["options"]["boundary_first"], True)
            old = copy.deepcopy(recipe)
            del old["options"]["boundary_first"]
            path.write_text(json.dumps(old), encoding="utf-8")
            normalized = wf.load_recipe(self.model, path)
            self.assertIs(normalized["options"]["boundary_first"], False)
            self.assertIs(normalized["options"]["native_source_fallback"], True)
        success = self.native.mesh.return_value
        self.native.mesh.side_effect = [synthetic_refusal("native-plc-cavity-material-refusal"),
                                       synthetic_refusal("native-plc-final-centroid-outside"), success]
        result = wf.run_recipe(self.model, loaded)
        for call, attempt in zip(self.native.mesh.call_args_list, result.report["attempts"]):
            for options in (call.kwargs, attempt["options"]):
                self.assertIs(options["boundary_first"], True)
                self.assertIs(options["native_source_fallback"], True)
                self.assertIs(options["advancing_interior_front"], True)

    def test_boundary_first_cli_and_discovery_describe_an_explicit_order_override(self):
        for flags, expected in (([], False), (["--boundary-first"], True),
                               (["--boundary-first", "--no-boundary-first"], False)):
            with self.subTest(flags=flags), patch.object(cli, "read_step", return_value=self.model), \
                 patch.object(cli, "_native_output_to_stderr", contextlib.nullcontext), \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cli.main(["source.step", "out.vtk", "--size", ".5", *flags]), 0)
            self.assertIs(self.native.mesh.call_args.kwargs["boundary_first"], expected)
            self.assertIs(self.native.mesh.call_args.kwargs["native_source_fallback"], False)
        control = pm.capabilities()["controls"]["boundary_first"]
        self.assertIs(control["default"], False)
        self.assertIs(control["experimental"], True)
        self.assertIs(control["overrides_enabled_front_order"], True)
        self.assertIs(control["forces_initial_front"], False)
        self.assertIs(control["enables_source_fallback"], False)
        self.assertIs(control["changes_native_validity_gates"], False)

    def test_witnesses_are_bounded_and_info_does_not_recertify_a_mesh(self):
        report = synthetic_report()
        report["requested_native_source_fallback"] = True
        report["requested_boundary_first"] = True
        report["bodies"][0]["source_shell"] = fallback_shell()
        before = copy.deepcopy(report)
        diagnostics = pm.diagnose(report, max_issues=1)
        self.assertEqual(diagnostics["status"], "valid")
        self.assertIs(diagnostics["summary"]["requested_native_source_fallback"], True)
        self.assertIs(diagnostics["summary"]["requested_boundary_first"], True)
        issue = diagnostics["issues"][0]
        self.assertEqual(issue["code"], "source_fallback_spacing_unverified")
        self.assertEqual(issue["severity"], "info")
        self.assertIs(issue["evidence"]["sizing_completion_certified"], False)
        shell = diagnostics["bodies"][0]["source_shell"]
        self.assertEqual(shell["attempts_total"], 2)
        self.assertIs(shell["attempts_truncated"], True)
        self.assertEqual(shell["attempts"][0]["witness"]["edge_tag"], 9)
        self.assertEqual(shell["attempts"][0]["witness"]["use_id"], 0)
        self.assertEqual(report, before)
        report["audit"]["valid"] = False
        self.assertEqual(pm.diagnose(report)["status"], "failed")

    def test_no_spacing_advisory_for_ordinary_native_or_unaccepted_source(self):
        for change in ({"native_fallback_attempted": False}, {"accepted": False}):
            report = synthetic_report()
            report["bodies"][0]["source_shell"] = {**fallback_shell(), **change}
            result = pm.diagnose(report)
            self.assertNotIn("source_fallback_spacing_unverified", [row["code"] for row in result["issues"]])

    def test_malformed_witness_or_provider_receipt_is_rejected(self):
        for invalid in ({"native_fallback_requested": 1}, {"native_fallback_diagnostic": 4},
                        {"attempts": [{"provider": []}]}, {"attempts": [{"global_sizing_applied": 1}]},
                        {"attempts": [{"angular_deflection": 0}]},
                        {"attempts": [{"witness": {"edge_tag": -1}}]},
                        {"attempts": [{"witness": {"geometric_refusal": 1}}]}):
            report = synthetic_report()
            report["bodies"][0]["source_shell"] = {**fallback_shell(), **invalid}
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                pm.diagnose(report)

    def test_material_progress_survives_a_reported_stop_without_certifying_output(self):
        report = synthetic_report()
        report.update(status="refused", reason="cancelled", mesh_published=False)
        report["audit"]["complete"] = report["quality"]["complete"] = False
        body = report["bodies"][0]
        body.update(material_classification=material_receipt(), native_candidate_queries=35,
                    native_candidate_ray_admissions=6)
        before = copy.deepcopy(report)
        result = pm.diagnose(report, max_issues=1)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["bodies"][0]["material_classification"], material_receipt())
        self.assertEqual(result["bodies"][0]["native_candidate_queries"], 35)
        self.assertEqual(result["bodies"][0]["native_candidate_ray_admissions"], 6)
        self.assertEqual(report, before)

    def test_completed_material_batch_does_not_override_final_native_failure(self):
        report = synthetic_report()
        report["audit"]["valid"] = False
        report["bodies"][0]["material_classification"] = material_receipt(complete=True)
        result = pm.diagnose(report)
        self.assertEqual(result["status"], "failed")
        self.assertIs(result["bodies"][0]["material_classification"]["complete"], True)
        control = pm.capabilities()["controls"]["boundary_first"]
        self.assertIn("original_trimmed_CAD", control["material_admission"])
        self.assertIs(control["changes_native_validity_gates"], False)
        self.assertIs(pm.capabilities()["diagnostics"]["material_progress_certifies_final_mesh"], False)

    def test_material_and_candidate_progress_counter_types_and_coherence_are_checked(self):
        for invalid in ({"complete": 1}, {"source_model_verified_unchanged": "true"},
                        {"work": -1}, {"seconds": -1}, {"requested_points": 2},
                        {"native_queries": 11}, {"ray_admissions": True}):
            report = synthetic_report()
            report["bodies"][0]["material_classification"] = {**material_receipt(), **invalid}
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                pm.diagnose(report)
        for name in ("native_candidate_queries", "native_candidate_ray_admissions"):
            report = synthetic_report()
            report["bodies"][0][name] = True
            with self.subTest(name=name), self.assertRaises(ValueError):
                pm.diagnose(report)

    def test_material_timing_fields_reject_nonfinite_values_without_json_sanitizing(self):
        fields = ("seconds", "classifier_setup_seconds", "native_query_seconds", "ray_witness_seconds")
        for name in fields:
            for invalid in (float("nan"), float("inf"), -float("inf"), -1, True, "1", 10**1000):
                receipt = {**material_receipt(), name: invalid}
                with self.subTest(field=name, invalid=repr(invalid)), self.assertRaises(ValueError):
                    diagnostics._material_classification(receipt)
                report = synthetic_report()
                report["bodies"][0]["material_classification"] = receipt
                with self.subTest(public_field=name, invalid=repr(invalid)), self.assertRaises(ValueError):
                    pm.diagnose(report)

    def test_material_profile_counts_and_elapsed_partition_are_checked(self):
        for invalid in ({"profiled_batches": True}, {"profiled_batches": 3},
                        {"ray_witness_attempts": -1}, {"classifier_setup_seconds": None},
                        {"native_query_seconds": .2}):
            report = synthetic_report()
            report["bodies"][0]["material_classification"] = {**material_receipt(), **invalid}
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                pm.diagnose(report)
        unmeasured = {**material_receipt(), "profiled_batches": 0, "ray_witness_attempts": 0,
                      "classifier_setup_seconds": 0, "native_query_seconds": 0, "ray_witness_seconds": 0}
        report = synthetic_report()
        report["bodies"][0]["material_classification"] = unmeasured
        self.assertEqual(pm.diagnose(report)["bodies"][0]["material_classification"], unmeasured)
        self.assertIs(diagnostics._material_classification({"seconds": None})["seconds"], None)


class NativeSourceFallbackContractTests(unittest.TestCase):
    """Real binding/primitive checks; run only after the matching Native build."""

    def test_native_boolean_validation_and_explicit_off_constructor_are_preserved(self):
        source = fixture_path()
        original = source.read_bytes()
        model = pm.read_step(source)
        for invalid in (None, 0, 1, "true", []):
            with self.subTest(value=invalid), self.assertRaises(TypeError):
                model._native.mesh(.7, native_source_fallback=invalid)
            with self.subTest(order=invalid), self.assertRaises(TypeError):
                model._native.mesh(.7, boundary_first=invalid)
        options = {"assembly": "conforming", "advancing_interior_front": False,
                   "workers": 1, "validation_workers": 1, "max_work": 50_000_000}
        ordinary = model.mesh(.7, native_source_fallback=False, boundary_first=False, **options)
        enabled = model.mesh(.7, native_source_fallback=True, boundary_first=True, **options)
        for result, requested in ((ordinary, False), (enabled, True)):
            report = result.report
            self.assertIs(report["requested_native_source_fallback"], requested)
            self.assertIs(report["effective_policy"]["requested_native_source_fallback"], requested)
            self.assertIs(report["requested_boundary_first"], requested)
            self.assertIs(report["effective_policy"]["requested_boundary_first"], requested)
            self.assertIs(report["requested_advancing_interior_front"], False)
            self.assertTrue(report["audit"]["valid"] and report["quality"]["valid"])
            self.assertTrue(all(body["source_shell"]["native_fallback_attempted"] is False
                                and body["initial_front"]["diagnostic"] == "not_attempted"
                                for body in report["bodies"]))
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory) / name for name in ("ordinary.paacsmesh", "enabled.paacsmesh")]
            ordinary.write_archive(paths[0]);enabled.write_archive(paths[1])
            self.assertEqual(paths[0].read_bytes(), paths[1].read_bytes())
        self.assertEqual(ordinary.face_tags, enabled.face_tags)
        self.assertEqual(ordinary.body_tags, enabled.body_tags)
        self.assertEqual(source.read_bytes(), original)

    def test_boundary_first_live_route_reports_deferred_then_actual_interior_work(self):
        source = fixture_path()
        original = source.read_bytes()
        model = pm.read_step(source)
        mesh = model.mesh(.7, assembly="conforming", advancing_interior_front=True,
                          boundary_first=True, max_work=100_000_000, max_seconds=60,
                          workers=1, validation_workers=1)
        report = mesh.report
        self.assertIs(report["requested_boundary_first"], True)
        self.assertIs(report["requested_native_source_fallback"], False)
        self.assertTrue(report["audit"]["valid"] and report["quality"]["valid"])
        body = report["bodies"][0]
        self.assertEqual(body["initial_front"]["diagnostic"], "deferred_until_boundary_recovery")
        self.assertEqual(body["initial_front"]["nodes"], 0)
        self.assertGreater(body["constraint_recovery"]["seed_cells"], 0)
        self.assertNotEqual(body["interior_front"]["diagnostic"], "not_attempted")
        material = body["material_classification"]
        self.assertGreater(material["attempts"], 0)
        self.assertIs(material["complete"], True)
        self.assertIs(material["source_model_verified_unchanged"], True)
        self.assertLessEqual(material["attempted_points"], material["requested_points"])
        self.assertEqual(material["attempted_points"],
            material["native_queries"]+material["strict_reuses"]+material["ray_admissions"])
        self.assertGreater(material["profiled_batches"], 0)
        self.assertLessEqual(material["profiled_batches"], material["attempts"])
        self.assertGreaterEqual(material["ray_witness_attempts"], material["ray_admissions"])
        phases = ("classifier_setup_seconds", "native_query_seconds", "ray_witness_seconds")
        self.assertTrue(all(material[name] >= 0 for name in phases))
        self.assertLessEqual(sum(material[name] for name in phases), material["seconds"]+1e-9)
        for name in ("native_candidate_queries", "native_candidate_ray_admissions"):
            self.assertIs(type(body[name]), int)
            self.assertGreaterEqual(body[name], 0)
        self.assertNotEqual(pm.diagnose(report)["status"], "failed")
        self.assertEqual(mesh.face_tags, model.face_tags)
        self.assertEqual(mesh.body_tags, model.body_tags)
        self.assertEqual(source.read_bytes(), original)

    def test_boundary_first_caller_respects_explicit_interior_upper_caps(self):
        source = fixture_path()
        original = source.read_bytes()
        model = pm.read_step(source)
        try:
            mesh = model.mesh(.7, assembly="conforming", advancing_interior_front=True,
                boundary_first=True, max_interior_candidates=1, max_interior_sites=1,
                max_work=100_000_000, max_seconds=60, workers=1, validation_workers=1)
        except pm.MeshingError as error:
            report = error.report
            self.assertIs(report["mesh_published"], False)
            self.assertNotEqual(report["reason"], "none")
        else:
            report = mesh.report
            self.assertTrue(report["audit"]["valid"] and report["quality"]["valid"])
            self.assertEqual(mesh.face_tags, model.face_tags)
            self.assertEqual(mesh.body_tags, model.body_tags)
        self.assertEqual(report["construction"]["max_interior_candidates"], 1)
        self.assertEqual(report["construction"]["max_interior_sites"], 1)
        body = report["bodies"][0]
        self.assertEqual(body["initial_front"]["diagnostic"], "deferred_until_boundary_recovery")
        self.assertGreater(body["constraint_recovery"]["seed_cells"], 0)
        self.assertLessEqual(body["interior_front"]["candidates"], 1)
        self.assertLessEqual(body["interior_front"]["nodes"], 1)
        self.assertEqual(source.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
