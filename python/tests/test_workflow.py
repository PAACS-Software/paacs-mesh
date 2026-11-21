"""Installed-package checks for bounded recipes and local publication."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import copy
import json
from pathlib import Path
import shutil
import tempfile
import time
import unittest
from unittest.mock import patch

import paacs_mesh as pm
from paacs_mesh import workflow as wf

from test_native import fixture_path


class WorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = fixture_path()
        cls.temporary = tempfile.TemporaryDirectory(prefix="paacs-workflow-tests-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.directory = Path(cls.temporary.name)
        cls.model = pm.read_step(cls.fixture)
        cls.model.define_group("wall", cls.model.select_faces(surface="cylinder"), "wall-id")
        cls.baseline = cls.model.mesh(.7)
        cls.recipe = pm.create_recipe(cls.model, .7)
        cls.result = pm.run_recipe(cls.model, cls.recipe)

    def assert_valid(self, mesh):
        report = mesh.report
        self.assertTrue(report["mesh_published"])
        self.assertEqual(report["reason"], "none")
        for stage in ("audit", "quality"):
            self.assertTrue(report[stage]["complete"])
            self.assertTrue(report[stage]["valid"])
        self.assertEqual(report["quality"]["inverted_cells"], 0)

    def archive_bytes(self, mesh, name):
        path = self.directory / name
        mesh.write_archive(path)
        return path.read_bytes()

    def synthetic_error(self, reason, diagnostic):
        """Exercise route policy without requiring a slow failing CAD fixture."""
        report = copy.deepcopy(self.baseline.report)
        report.update(status="refused", mesh_published=False, reason=reason)
        report["bodies"][0].update(reason=reason, volume_reason=diagnostic)
        error = pm.MeshingError("Synthetic failure for workflow policy qualification")
        error.stage, error.reason, error.report = "mesh", reason, report
        return error

    def test_default_run_matches_native_mesh_and_keeps_part_scope(self):
        self.assertIs(pm.create_recipe, wf.create_recipe)
        self.assert_valid(self.result.mesh)
        self.assertEqual(self.archive_bytes(self.result.mesh, "workflow-default.paacsmesh"),
                         self.archive_bytes(self.baseline, "native-default.paacsmesh"))
        report = self.result.report
        self.assertEqual(len(report["attempts"]), 1)
        self.assertEqual(report["attempts"][0]["options"]["recovery"], "standard")
        self.assertEqual(report["attempts"][0]["options"]["assembly"], "preserve_parts")
        self.assertEqual(self.result.mesh.report["mesh_scope"], "independent_parts")
        self.assertFalse(self.result.mesh.report["interbody_nonoverlap_certified"])
        self.assertTrue(report["cooperative_deadline"])
        self.assertFalse(report["publication_included"])
        for key in ("warning_quality", "minimum_quality"):
            self.assertIsNone(self.recipe["options"][key])
            self.assertIsNone(report["attempts"][0]["options"][key])

    def test_selected_face_and_box_controls_replay_through_native_engine(self):
        model = pm.read_step(self.fixture)
        wall = model.select_faces(surface="cylinder", radius=1)
        recipe = pm.create_recipe(model, .7, selected_face_sizes=[(wall, .16)],
            volume_controls=[pm.BoxSize((-.3, -.3, .7), (.3, .3, 1.3), .15)],
            min_size=.12, growth=1.55, assembly="conforming", max_work=100_000_000)
        saved = self.directory / "controlled-recipe.json"
        pm.save_recipe(model, recipe, saved)
        reopened = pm.read_step(self.fixture)
        loaded = pm.load_recipe(reopened, saved)
        self.assertEqual(loaded, recipe)
        result = pm.run_recipe(reopened, loaded)
        self.assert_valid(result.mesh)
        actual = result.report["attempts"][0]["options"]
        self.assertEqual(actual["face_sizes"], [{"tag": wall.tags[0], "size": .16, "growth": None}])
        self.assertEqual(actual["volume_controls"][0]["size"], .15)
        self.assertEqual(actual["min_size"], .12)
        self.assertEqual(actual["growth"], 1.55)
        self.assertEqual(actual["max_work"], 100_000_000)
        sizing = result.mesh.report["sizing"]
        self.assertTrue(sizing["native_sizing"])
        self.assertEqual(sizing["face_sizes"][0]["growth"], 1.55)
        self.assertEqual(result.mesh.face_tags, [1, 2, 3])
        self.assertEqual(result.mesh.body_tags, [1])

    def test_native_budget_failure_is_not_retried_or_inflated(self):
        with self.assertRaises(pm.MeshingError) as direct:
            self.model.mesh(.7, max_nodes=1)
        recipe = pm.create_recipe(self.model, .7, max_nodes=1, max_attempts=3,
                                  allow_advanced=True, finer_factors=[.5])
        with self.assertRaises(wf.WorkflowError) as caught:
            pm.run_recipe(self.model, recipe)
        self.assertEqual(caught.exception.reason, direct.exception.reason)
        attempts = caught.exception.report["attempts"]
        self.assertEqual(len(attempts), 1)
        self.assertEqual(attempts[0]["reason"], direct.exception.reason)
        self.assertEqual(attempts[0]["native_report"]["reason"], direct.exception.report["reason"])
        self.assertEqual(attempts[0]["options"]["max_nodes"], 1)
        self.assertEqual(attempts[0]["options"]["recovery"], "standard")
        self.assertFalse(attempts[0]["native_report"]["mesh_published"])
        self.assert_valid(self.result.mesh)

    def test_recipe_inputs_and_result_metadata_are_independent_copies(self):
        faces = [pm.FaceSize(1, .2)]
        boxes = [pm.BoxSize((-.2, -.2, .5), (.2, .2, 1.5), .15)]
        recipe = pm.create_recipe(self.model, .7, face_sizes=faces, volume_controls=boxes,
                                  warning_quality=.6, minimum_quality=.2)
        faces.append(pm.FaceSize(2, .1))
        boxes.clear()
        self.assertEqual(len(recipe["options"]["face_sizes"]), 1)
        self.assertEqual(len(recipe["options"]["volume_controls"]), 1)
        self.assertEqual(recipe["options"]["warning_quality"], .6)
        self.assertEqual(recipe["options"]["minimum_quality"], .2)
        copied_recipe, copied_report = self.result.recipe, self.result.report
        copied_recipe["options"]["size"] = 999
        copied_recipe["options"]["warning_quality"] = 1
        copied_report["attempts"].clear()
        self.assertEqual(self.result.recipe["options"]["size"], .7)
        self.assertIsNone(self.result.recipe["options"]["warning_quality"])
        self.assertEqual(len(self.result.report["attempts"]), 1)
        with self.assertRaises(AttributeError):
            self.result.mesh = self.baseline
        foreign = pm.read_step(self.fixture).select_faces(surface="cylinder")
        with self.assertRaises(wf.WorkflowError) as caught:
            pm.create_recipe(self.model, .7, selected_face_sizes=[(foreign, .1)])
        self.assertEqual(caught.exception.reason, "foreign_selection")

    def test_opted_recovery_records_standard_advanced_and_finer_attempts(self):
        model = pm.read_step(self.fixture)
        recipe = pm.create_recipe(model, .7, allow_advanced=True, finer_factors=[.5], max_attempts=3,
                                  max_nodes=1_000_000, max_work=5_000_000,
                                  warning_quality=1, minimum_quality=0)
        failures = [self.synthetic_error("geometry_error", "native-plc-cavity-material-refusal"),
                    self.synthetic_error("geometry_error", "native-plc-final-centroid-outside")]
        original = pm.CadModel.mesh
        def attempt(actual_model, size, **options):
            if failures:
                raise failures.pop(0)
            return original(actual_model, size, **options)
        with patch.object(pm.CadModel, "mesh", autospec=True, side_effect=attempt):
            result = pm.run_recipe(model, recipe)
        self.assert_valid(result.mesh)
        attempts = result.report["attempts"]
        self.assertEqual([item["status"] for item in attempts], ["failed", "failed", "succeeded"])
        self.assertEqual([item["options"]["recovery"] for item in attempts], ["standard", "advanced", "advanced"])
        self.assertEqual([item["options"]["size"] for item in attempts], [.7, .7, .35])
        self.assertEqual(attempts[0]["native_report"]["bodies"][0]["volume_reason"], "native-plc-cavity-material-refusal")
        for item in attempts:
            self.assertEqual(item["options"]["max_nodes"], 1_000_000)
            self.assertEqual(item["options"]["max_work"], 5_000_000)
            self.assertEqual(item["options"]["assembly"], "preserve_parts")
            self.assertEqual(item["options"]["warning_quality"], 1)
            self.assertEqual(item["options"]["minimum_quality"], 0)
            self.assertGreaterEqual(item["seconds"], 0)
        self.assertGreaterEqual(result.report["seconds"], sum(item["seconds"] for item in attempts))

    def test_quality_thresholds_persist_forward_and_legacy_recipes_keep_defaults(self):
        model = pm.read_step(self.fixture)
        minimum = self.baseline.report["quality"]["minimum"] / 2
        recipe = pm.create_recipe(model, .7, warning_quality=1, minimum_quality=minimum,
                                  allow_advanced=True, finer_factors=[.5], max_attempts=3)
        saved = self.directory / "quality-recipe.json"
        pm.save_recipe(model, recipe, saved)
        loaded = pm.load_recipe(model, saved)
        self.assertEqual(loaded["options"]["warning_quality"], 1)
        self.assertEqual(loaded["options"]["minimum_quality"], minimum)
        original, observed = pm.CadModel.mesh, []
        def forwarded(actual_model, size, **options):
            observed.append({key: options[key] for key in ("warning_quality", "minimum_quality")})
            return original(actual_model, size, **options)
        with patch.object(pm.CadModel, "mesh", autospec=True, side_effect=forwarded):
            result = pm.run_recipe(model, loaded)
        self.assert_valid(result.mesh)
        self.assertTrue(result.mesh.report["quality_warning"])
        self.assertEqual(observed, [{"warning_quality": 1, "minimum_quality": minimum}])
        self.assertEqual(len(result.report["attempts"]), 1)
        actual = result.report["attempts"][0]["options"]
        self.assertEqual(actual["warning_quality"], 1)
        self.assertEqual(actual["minimum_quality"], minimum)
        # Recipe v1 originally had neither quality option. Loading adds only
        # their disabled/default requests and does not alter native output.
        legacy = copy.deepcopy(self.recipe)
        for key in ("warning_quality", "minimum_quality"):
            del legacy["options"][key]
        old_path = self.directory / "legacy-quality-recipe.json"
        old_path.write_text(json.dumps(legacy), encoding="utf-8")
        restored = pm.load_recipe(self.model, old_path)
        self.assertIsNone(restored["options"]["warning_quality"])
        self.assertIsNone(restored["options"]["minimum_quality"])
        replayed = pm.run_recipe(self.model, restored)
        self.assertEqual(self.archive_bytes(replayed.mesh, "legacy-quality-default.paacsmesh"),
                         self.archive_bytes(self.baseline, "quality-native-default.paacsmesh"))
        self.assertEqual(legacy["options"].keys(), self.recipe["options"].keys() - {"warning_quality", "minimum_quality"})
        for changed in (dict(legacy["options"], unknown_quality_option=.5),
                        {key: value for key, value in legacy["options"].items() if key != "growth"}):
            invalid = copy.deepcopy(legacy)
            invalid["options"] = changed
            old_path.write_text(json.dumps(invalid), encoding="utf-8")
            with self.assertRaises(wf.WorkflowError):
                pm.load_recipe(self.model, old_path)

    def test_quality_threshold_validation_precedes_native_meshing(self):
        with patch.object(pm.CadModel, "mesh") as calls:
            for key in ("warning_quality", "minimum_quality"):
                for value in (-.01, 1.01, float("nan"), float("inf"), True, "0.5"):
                    with self.subTest(option=key, value=value), self.assertRaises(wf.WorkflowError):
                        pm.create_recipe(self.model, .7, **{key: value})
                for value in (None, 0, .5, 1):
                    with self.subTest(option=key, value=value):
                        recipe = pm.create_recipe(self.model, .7, **{key: value})
                        self.assertEqual(recipe["options"][key], value)
            self.assertEqual(calls.call_count, 0)

    def test_total_wall_budget_counts_preparation_failures_and_overruns(self):
        model = pm.read_step(self.fixture)
        recipe = pm.create_recipe(model, .7, max_seconds=.1, allow_advanced=True, max_attempts=3)
        failure = self.synthetic_error("geometry_error", "native-plc-cavity-material-refusal")
        def slow_failure(*args, **kwargs):
            time.sleep(.15)
            raise failure
        with patch.object(pm.CadModel, "mesh", side_effect=slow_failure) as calls:
            with self.assertRaises(wf.WorkflowError) as caught:
                pm.run_recipe(model, recipe)
        self.assertEqual(caught.exception.reason, "time_budget")
        self.assertEqual(calls.call_count, 1)
        self.assertEqual(len(caught.exception.report["attempts"]), 1)
        self.assertGreaterEqual(caught.exception.report["seconds"], .15)
        original_options = wf._native_options
        def slow_preparation(*args, **kwargs):
            time.sleep(.15)
            return original_options(*args, **kwargs)
        with patch.object(wf, "_native_options", side_effect=slow_preparation), patch.object(pm.CadModel, "mesh") as calls:
            with self.assertRaises(wf.WorkflowError) as caught:
                pm.run_recipe(model, recipe)
        self.assertEqual(caught.exception.reason, "time_budget")
        self.assertEqual(calls.call_count, 0)
        self.assertEqual(caught.exception.report["attempts"], [])
        def slow_success(*args, **kwargs):
            time.sleep(.15)
            return self.baseline
        with patch.object(pm.CadModel, "mesh", side_effect=slow_success):
            overrun = pm.run_recipe(model, recipe)
        self.assert_valid(overrun.mesh)
        self.assertTrue(overrun.report["time_budget_overrun"])

    def test_published_source_intent_recipe_and_archive_reopen_together(self):
        destination = self.directory / "roundtrip"
        publication = pm.publish_result(self.result, destination)
        self.assertEqual(set(publication["files"]), {path.name for path in destination.iterdir()})
        self.assertEqual((destination / "source.step").read_bytes(), self.fixture.read_bytes())
        reopened = pm.read_step(destination / "source.step", intent=destination / "intent.json")
        expected_groups = {group["name"]: {key: value for key, value in group.items() if key != "name"}
                           for group in self.result.recipe["intent"]["groups"]}
        self.assertEqual(reopened.groups, expected_groups)
        recipe = pm.load_recipe(reopened, destination / "recipe.json")
        self.assertEqual(recipe, self.result.recipe)
        replayed = wf.run_saved_recipe(destination / "recipe.json", destination / "source.step")
        self.assert_valid(replayed.mesh)
        self.assertGreater(replayed.report["import_seconds"], 0)
        self.assertEqual(self.archive_bytes(replayed.mesh, "replayed.paacsmesh"),
                         (destination / "mesh.paacsmesh").read_bytes())
        self.assertEqual(json.loads((destination / "report.json").read_text())["source_sha256"], reopened.source_sha256)
        self.assertTrue((destination / "mesh.vtk").read_text().startswith("# vtk DataFile"))

    def test_publication_uses_captured_groups_and_refuses_existing_outputs(self):
        model = pm.read_step(self.fixture)
        model.define_group("wall", model.select_faces(surface="cylinder"), "wall-id")
        result = pm.run_recipe(model, pm.create_recipe(model, .7))
        model.define_group("cap", model.select_faces(tags=[2]), "cap-id")
        destination = self.directory / "captured-groups"
        pm.publish_result(result, destination)
        intent = json.loads((destination / "intent.json").read_text())
        self.assertEqual([group["name"] for group in intent["groups"]], ["wall"])
        original = {path.name: path.read_bytes() for path in destination.iterdir()}
        for output in (destination, self.directory / "existing-empty"):
            output.mkdir(exist_ok=True)
            with self.assertRaises(wf.WorkflowError) as caught:
                pm.publish_result(result, output)
            self.assertEqual(caught.exception.reason, "output_exists")
        self.assertEqual({path.name: path.read_bytes() for path in destination.iterdir()}, original)
        self.assertEqual(list((self.directory / "existing-empty").iterdir()), [])

    def test_failed_export_removes_staging_and_preserves_owned_mesh(self):
        class FailedExport:
            def __init__(self, mesh):
                self.mesh = mesh
            def write_archive(self, path):
                self.mesh.write_archive(path)
            def write_vtk(self, path):
                raise OSError("Deliberate publication export failure")
        result = wf.WorkflowResult(self.result._model, FailedExport(self.result.mesh),
                                    self.result._recipe_json, self.result._report_json)
        destination = self.directory / "failed-export"
        with self.assertRaises(wf.WorkflowError) as caught:
            pm.publish_result(result, destination)
        self.assertEqual(caught.exception.stage, "publication")
        self.assertFalse(destination.exists())
        self.assertEqual(list(self.directory.glob(".failed-export.paacs-*")), [])
        self.assert_valid(self.result.mesh)
        self.assertEqual(self.archive_bytes(self.result.mesh, "still-owned.paacsmesh"),
                         self.archive_bytes(self.baseline, "still-native.paacsmesh"))

    def test_source_changes_embedding_budget_and_retained_bytes(self):
        source = self.directory / "copied.step"
        shutil.copyfile(self.fixture, source)
        model = pm.read_step(source)
        result = pm.run_recipe(model, pm.create_recipe(model, .7))
        saved = self.directory / "copied-recipe.json"
        pm.save_recipe(model, result.recipe, saved)
        with self.assertRaises(wf.WorkflowError) as caught:
            pm.publish_result(result, self.directory / "source-budget", max_source_bytes=1)
        self.assertEqual(caught.exception.reason, "source_budget")
        with self.assertRaises(wf.WorkflowError):
            pm.save_recipe(model, result.recipe, source)
        self.assertEqual(source.read_bytes(), self.fixture.read_bytes())
        source.write_bytes(source.read_bytes() + b"\n")
        with self.assertRaises(wf.WorkflowError) as caught:
            pm.publish_result(result, self.directory / "changed-source")
        self.assertEqual(caught.exception.reason, "source_changed")
        with self.assertRaises(wf.WorkflowError) as caught:
            pm.load_recipe(pm.read_step(source), saved)
        self.assertEqual(caught.exception.reason, "source_mismatch")
        self.assertFalse((self.directory / "changed-source").exists())
        retained = pm.read_step(source, retain_source=True)
        kept = pm.run_recipe(retained, pm.create_recipe(retained, .7))
        payload = source.read_bytes()
        source.unlink()
        pm.publish_result(kept, self.directory / "retained-source")
        self.assertEqual((self.directory / "retained-source/source.step").read_bytes(), payload)
        self.assert_valid(result.mesh)

    def test_invalid_recipes_refuse_before_any_native_attempt(self):
        invalids = ({"max_attempts": 9}, {"max_seconds": float("nan")}, {"finer_factors": [.5, .7]},
                    {"allow_advanced": 1}, {"growth": .9})
        with patch.object(pm.CadModel, "mesh") as calls:
            for options in invalids:
                with self.subTest(options=options), self.assertRaises(wf.WorkflowError):
                    pm.create_recipe(self.model, .7, **options)
            recipe = copy.deepcopy(self.recipe)
            recipe["schema"] = "unsupported"
            with self.assertRaises(wf.WorkflowError) as caught:
                pm.run_recipe(self.model, recipe)
            self.assertEqual(caught.exception.reason, "unknown_schema")
            self.assertEqual(calls.call_count, 0)
        duplicate = self.directory / "duplicate-recipe.json"
        duplicate.write_text('{"schema":"one","schema":"two"}', encoding="utf-8")
        with self.assertRaises(wf.WorkflowError):
            pm.load_recipe(self.model, duplicate)
        self.assert_valid(self.result.mesh)

    def test_report_budget_preserves_attempts_and_marks_omitted_reports(self):
        original = self.archive_bytes(self.baseline, "before-report-budget.paacsmesh")
        for label, payload_bytes, expected_calls in (("single", 9 * 1024 * 1024, 1),
                                                    ("aggregate", 3 * 1024 * 1024, 3)):
            with self.subTest(report=label):
                failure = self.synthetic_error("geometry_error", "native-plc-cavity-material-refusal")
                failure.report["large_diagnostic"] = "x" * payload_bytes
                recipe = pm.create_recipe(self.model, .7, allow_advanced=True, max_attempts=8,
                                          finer_factors=[.5, .25, .125, .0625, .03125, .015625, .0078125])
                with patch.object(pm.CadModel, "mesh", side_effect=failure) as calls:
                    with self.assertRaises(wf.WorkflowError) as caught:
                        pm.run_recipe(self.model, recipe)
                error, report = caught.exception, caught.exception.report
                self.assertEqual(error.reason, "report_budget")
                self.assertEqual(calls.call_count, expected_calls)
                self.assertEqual(len(report["attempts"]), expected_calls)
                self.assertFalse(report["diagnostics_complete"])
                self.assertLess(len(json.dumps(report).encode("utf-8")), 8 * 1024 * 1024)
                omitted = report["attempts"][-1]
                self.assertTrue(omitted["options_complete"])
                self.assertFalse(omitted["native_report_complete"])
                self.assertIsNone(omitted["native_report"])
                self.assertEqual(omitted["native_report_summary"]["reason"], "geometry_error")
                self.assertGreaterEqual(omitted["seconds"], 0)
                self.assertGreaterEqual(report["seconds"], sum(item["seconds"] for item in report["attempts"]))
                for retained in report["attempts"][:-1]:
                    self.assertTrue(retained["native_report_complete"])
                    self.assertEqual(retained["native_report"], failure.report)
        class LargeValidReport:
            @property
            def report(inner_self):
                report = self.baseline.report
                report["large_diagnostic"] = "x" * (9 * 1024 * 1024)
                return report
        with patch.object(pm.CadModel, "mesh", return_value=LargeValidReport()) as calls:
            with self.assertRaises(wf.WorkflowError) as caught:
                pm.run_recipe(self.model, self.recipe)
        self.assertEqual(caught.exception.reason, "report_budget")
        self.assertEqual(calls.call_count, 1)
        attempt = caught.exception.report["attempts"][0]
        self.assertEqual(attempt["status"], "succeeded")
        self.assertTrue(attempt["native_report_summary"]["mesh_published"])
        self.assertTrue(attempt["native_report_summary"]["audit"]["valid"])
        fallback = wf.WorkflowError("original", "Report fallback", {
            "attempts": [{"number": 1, "seconds": 1.25, "native_report": {"detail": "x" * (9 * 1024 * 1024)}}]})
        self.assertEqual(fallback.reason, "report_budget")
        self.assertEqual(fallback.report["attempts"][0]["seconds"], 1.25)
        self.assertFalse(fallback.report["diagnostics_complete"])
        self.assertEqual((self.directory / "before-report-budget.paacsmesh").read_bytes(), original)
        self.assert_valid(self.result.mesh)

    def test_racing_empty_directory_is_preserved_without_partial_output(self):
        destination = self.directory / "racing-output"
        original_rename = wf._rename_new_directory
        def racing_writer(staging, target):
            target.mkdir()
            original_rename(staging, target)
        with patch.object(wf, "_rename_new_directory", side_effect=racing_writer):
            with self.assertRaises(wf.WorkflowError) as caught:
                pm.publish_result(self.result, destination)
        self.assertEqual(caught.exception.reason, "output_exists")
        self.assertEqual(list(destination.iterdir()), [])
        self.assertEqual(list(self.directory.glob(".racing-output.paacs-*")), [])


if __name__ == "__main__":
    unittest.main()
