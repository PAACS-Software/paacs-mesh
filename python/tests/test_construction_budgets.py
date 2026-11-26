"""Construction option validation, guarded replay and native report acceptance."""

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
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import paacs_mesh as pm
from paacs_mesh import __main__ as cli
from paacs_mesh import workflow as wf
from paacs_mesh.cad import CadModel, _CONSTRUCTION_FIELDS, _CONSTRUCTION_LIMITS


SETTINGS = {"construction_effort": "thorough", "max_surface_points": 400,
            "max_surface_chord_points": 128, "max_interior_candidates": 20000,
            "max_interior_sites": 256}


def good_report():
    return {"status": "succeeded", "reason": "none", "mesh_published": True,
            "audit": {"complete": True, "valid": True},
            "quality": {"complete": True, "valid": True}, "bodies": []}


def model_stub(sha256="a" * 64):
    native = SimpleNamespace(info={"source_path": "source.step", "file_bytes": 12},
        model_identity=1, source_sha256=sha256, kernel_version="8.0.1",
        face_tags=[], body_tags=[1], edge_tags=[], vertex_tags=[],
        mesh=Mock(return_value=SimpleNamespace(report=good_report(), write=Mock())))
    return CadModel(native), native


def synthetic_refusal(diagnostic):
    report = good_report()
    report.update(status="refused", reason="geometry_error", mesh_published=False,
                  bodies=[{"solid_tag": 1, "reason": "geometry_error", "surface_reason": "none",
                           "volume_reason": diagnostic}])
    error = pm.MeshingError("Synthetic refusal for retry-forwarding qualification")
    error.stage, error.reason, error.report = "mesh", "geometry_error", report
    return error


class ConstructionBudgetTests(unittest.TestCase):
    """Pure Python contract checks; no meshing capability is simulated as native evidence."""

    def setUp(self):
        self.model, self.native = model_stub()

    def test_explicit_signatures_and_unchanged_legacy_defaults(self):
        for function in (CadModel.mesh, wf.create_recipe):
            parameters = inspect.signature(function).parameters
            self.assertEqual(parameters["construction_effort"].default, "legacy")
            self.assertEqual(parameters["assembly"].default, "preserve_parts")
            for name in _CONSTRUCTION_LIMITS:
                self.assertIsNone(parameters[name].default)
                self.assertEqual(parameters[name].kind, inspect.Parameter.KEYWORD_ONLY)
            self.assertFalse(any(value.kind == inspect.Parameter.VAR_KEYWORD for value in parameters.values()))
        self.model.mesh(.5)
        actual = self.native.mesh.call_args.kwargs
        self.assertEqual(actual["construction_effort"], "legacy")
        self.assertTrue(all(actual[name] is None for name in _CONSTRUCTION_LIMITS))
        self.assertEqual(actual["assembly"], "preserve_parts")
        self.assertEqual(actual["workers"], 1)
        self.assertEqual(actual["validation_workers"], 1)

    def test_exact_native_forwarding_preserves_existing_controls(self):
        self.model.mesh(.5, **SETTINGS, max_work=50_000_000, workers=2,
                        validation_workers=1, growth=1.4, min_size=.05,
                        recovery="standard", advancing_interior_front=None)
        actual = self.native.mesh.call_args.kwargs
        self.assertEqual(self.native.mesh.call_args.args, (.5,))
        self.assertEqual({name: actual[name] for name in _CONSTRUCTION_FIELDS}, SETTINGS)
        self.assertEqual(actual["max_work"], 50_000_000)
        self.assertEqual(actual["workers"], 2)
        self.assertEqual(actual["growth"], 1.4)
        self.assertEqual(actual["min_size"], .05)
        self.assertEqual(actual["recovery"], "standard")
        self.assertIsNone(actual["advancing_interior_front"])

    def test_invalid_effort_and_caps_never_enter_native_meshing(self):
        for value in (None, True, 1, "", "AUTO", "unbounded"):
            with self.subTest(effort=value), self.assertRaises(ValueError):
                self.model.mesh(.5, construction_effort=value)
        for name in _CONSTRUCTION_LIMITS:
            for value in (True, False, 0, -1, 1., "2", 2**63):
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    self.model.mesh(.5, **{name: value})
        self.native.mesh.assert_not_called()
        for name in _CONSTRUCTION_LIMITS:
            with self.subTest(recipe_field=name), self.assertRaises(wf.WorkflowError) as caught:
                wf.create_recipe(self.model, .5, **{name: 0})
            self.assertEqual(caught.exception.reason, "invalid_recipe")

    def test_guarded_recipe_save_load_and_old_missing_fields(self):
        recipe = wf.create_recipe(self.model, .5, **SETTINGS, max_work=50_000_000)
        self.assertEqual({name: recipe["options"][name] for name in _CONSTRUCTION_FIELDS}, SETTINGS)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "recipe.json"
            wf.save_recipe(self.model, recipe, path)
            self.assertEqual(wf.load_recipe(self.model, path), recipe)
            foreign, _ = model_stub("b" * 64)
            with self.assertRaises(wf.WorkflowError) as caught:
                wf.load_recipe(foreign, path)
            self.assertEqual(caught.exception.reason, "source_mismatch")
            old = copy.deepcopy(recipe)
            for name in _CONSTRUCTION_FIELDS:
                del old["options"][name]
            path.write_text(json.dumps(old), encoding="utf-8")
            restored = wf.load_recipe(self.model, path)
            self.assertEqual(restored["options"]["construction_effort"], "legacy")
            self.assertTrue(all(restored["options"][name] is None for name in _CONSTRUCTION_LIMITS))
            self.assertEqual(restored["options"]["max_work"], 50_000_000)
            # A partially extended recipe also has explicit historical defaults.
            old["options"]["max_surface_points"] = 123
            path.write_text(json.dumps(old), encoding="utf-8")
            restored = wf.load_recipe(self.model, path)
            self.assertEqual(restored["options"]["construction_effort"], "legacy")
            self.assertEqual(restored["options"]["max_surface_points"], 123)
            old["options"]["max_surface_point"] = 999
            path.write_text(json.dumps(old), encoding="utf-8")
            with self.assertRaises(wf.WorkflowError) as caught:
                wf.load_recipe(self.model, path)
            self.assertEqual(caught.exception.reason, "invalid_recipe")

    def test_every_retry_preserves_effort_caps_and_global_budget(self):
        success = self.native.mesh.return_value
        self.native.mesh.side_effect = [synthetic_refusal("native-plc-cavity-material-refusal"),
                                        synthetic_refusal("native-plc-final-centroid-outside"), success]
        recipe = wf.create_recipe(self.model, .5, **SETTINGS, max_work=50_000_000,
                                  allow_advanced=True, finer_factors=[.5], max_attempts=3)
        result = wf.run_recipe(self.model, recipe)
        calls = self.native.mesh.call_args_list
        self.assertEqual([call.args[0] for call in calls], [.5, .5, .25])
        self.assertEqual([call.kwargs["recovery"] for call in calls], ["standard", "advanced", "advanced"])
        for call, attempt in zip(calls, result.report["attempts"]):
            self.assertEqual({name: call.kwargs[name] for name in _CONSTRUCTION_FIELDS}, SETTINGS)
            self.assertEqual({name: attempt["options"][name] for name in _CONSTRUCTION_FIELDS}, SETTINGS)
            self.assertEqual(call.kwargs["max_work"], 50_000_000)
            self.assertEqual(attempt["options"]["max_work"], 50_000_000)
        self.assertEqual(recipe["options"]["size"], .5)

    def test_direct_cli_flags_and_recipe_run_cannot_override_saved_caps(self):
        args = ["source.step", "out.vtk", "--size", ".5", "--construction-effort", "thorough",
                "--max-surface-points", "400", "--max-surface-chord-points", "128",
                "--max-interior-candidates", "20000", "--max-interior-sites", "256", "--max-work", "50000000"]
        with patch.object(cli, "read_step", return_value=self.model), \
             patch.object(cli, "_native_output_to_stderr", contextlib.nullcontext), \
             contextlib.redirect_stdout(io.StringIO()) as captured:
            self.assertEqual(cli.main(args), 0)
        json.loads(captured.getvalue())
        actual = self.native.mesh.call_args.kwargs
        self.assertEqual({name: actual[name] for name in _CONSTRUCTION_FIELDS}, SETTINGS)
        self.assertEqual(actual["max_work"], 50_000_000)
        with self.assertRaises(ValueError):
            cli._command_parser("run").parse_args(["recipe.json", "--source", "source.step", "--output", "new", "--max-interior-sites", "999"])
        for value in ("0", "-1", "1.5", str(2**63)):
            with self.subTest(value=value), patch.object(cli, "read_step") as read, \
                 contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as caught:
                cli.main(["source.step", "out.vtk", "--size", ".5", "--max-surface-points", value])
            self.assertEqual(caught.exception.code, 2)
            read.assert_not_called()


class NativeConstructionBudgetTests(unittest.TestCase):
    """Strict rebuilt-binding acceptance; intentionally not skipped on older wheels."""

    def test_real_native_caps_replay_and_report_without_relaxing_validity(self):
        from test_native import fixture_path
        model = pm.read_step(fixture_path())
        options = dict(SETTINGS, assembly="conforming", workers=1, validation_workers=1,
                       advancing_interior_front=True, max_work=100_000_000, max_seconds=90)
        direct = model.mesh(.5, **options)
        report = direct.report
        construction = report["construction"]
        self.assertEqual(construction["effort"], SETTINGS["construction_effort"])
        for name in _CONSTRUCTION_LIMITS:
            self.assertEqual(construction[name], SETTINGS[name])
        self.assertEqual(report["recovery"], "standard")
        self.assertTrue(report["mesh_published"])
        for name in ("audit", "quality"):
            self.assertTrue(report[name]["complete"])
            self.assertTrue(report[name]["valid"])
        self.assertEqual(report["quality"]["inverted_cells"], 0)
        self.assertEqual(direct.body_tags, [1])
        self.assertEqual(direct.face_tags, [1, 2, 3])
        front = report["bodies"][0]["initial_front"]
        self.assertEqual(report["route_resolution"]["reason"], "explicit_front_enabled")
        self.assertGreater(front["work"], 0)
        self.assertLessEqual(front["candidates"], front["candidate_limit"])
        self.assertLessEqual(front["candidate_limit"], SETTINGS["max_interior_candidates"])
        self.assertLessEqual(front["nodes"], front["site_limit"])
        self.assertLessEqual(front["site_limit"], SETTINGS["max_interior_sites"])
        self.assertTrue(report["bodies"][0]["surface_spacing"])
        for spacing in report["bodies"][0]["surface_spacing"]:
            self.assertLessEqual(spacing["added_points"], spacing["point_limit"])
            self.assertLessEqual(spacing["point_limit"], SETTINGS["max_surface_points"])
            self.assertLessEqual(spacing["chord_added_points"], spacing["chord_point_limit"])
            self.assertLessEqual(spacing["chord_point_limit"], SETTINGS["max_surface_chord_points"])
        recipe = wf.create_recipe(model, .5, **options)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            wf.save_recipe(model, recipe, root / "recipe.json")
            replayed = pm.run_saved_recipe(root / "recipe.json", fixture_path())
            direct.write_archive(root / "direct.paacsmesh")
            replayed.mesh.write_archive(root / "replayed.paacsmesh")
            self.assertEqual((root / "direct.paacsmesh").read_bytes(), (root / "replayed.paacsmesh").read_bytes())
        self.assertEqual(replayed.mesh.report["construction"], construction)

    def test_old_recipe_fields_preserve_real_native_legacy_policy(self):
        from test_native import fixture_path
        model = pm.read_step(fixture_path())
        recipe = wf.create_recipe(model, .5, construction_effort="auto", max_seconds=90)
        for name in _CONSTRUCTION_FIELDS:
            del recipe["options"][name]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "old-recipe.json"
            path.write_text(json.dumps(recipe), encoding="utf-8")
            loaded = wf.load_recipe(model, path)
            result = wf.run_recipe(model, loaded)
        self.assertEqual(result.mesh.report["construction"]["effort"], "legacy")
        for name in _CONSTRUCTION_LIMITS:
            self.assertIsNone(result.mesh.report["construction"][name])


if __name__ == "__main__":
    unittest.main()
