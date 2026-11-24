"""Pure forwarding and recipe contracts; no native mesh is constructed here."""

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
from paacs_mesh import cad, workflow as wf, __main__ as cli
from test_construction_budgets import model_stub, synthetic_refusal
from test_construction_diagnostics import synthetic_report


class FidelityControlTests(unittest.TestCase):
    def setUp(self):
        self.model, self.native = model_stub()

    def test_nullable_public_parameters_and_explicit_forwarding(self):
        for function in (cad.CadModel.mesh, wf.create_recipe):
            for name in (*cad._FIDELITY_FIELDS, *cad._CURVATURE_FIELDS, "initial_topology_quality"):
                parameter = inspect.signature(function).parameters[name]
                self.assertIsNone(parameter.default)
                self.assertEqual(parameter.kind, inspect.Parameter.KEYWORD_ONLY)
        self.model.mesh(.5)
        self.assertIsNone(self.native.mesh.call_args.kwargs["max_chord_error"])
        self.assertIsNone(self.native.mesh.call_args.kwargs["max_relative_volume_error"])
        self.model.mesh(.5, max_chord_error=.0025, max_relative_volume_error=.005)
        self.assertEqual(self.native.mesh.call_args.kwargs["max_chord_error"], .0025)
        self.assertEqual(self.native.mesh.call_args.kwargs["max_relative_volume_error"], .005)

    def test_invalid_fidelity_is_rejected_before_native_meshing(self):
        for name in cad._FIDELITY_FIELDS:
            values = (0, -1, True, "0.1", float("nan"), float("inf"))
            if name == "max_relative_volume_error":
                values += (1, 2)
            for value in values:
                with self.subTest(name=name, value=value), self.assertRaises((TypeError, ValueError)):
                    self.model.mesh(.5, **{name: value})
                with self.subTest(recipe=name, value=value), self.assertRaises(wf.WorkflowError):
                    wf.create_recipe(self.model, .5, **{name: value})
        self.native.mesh.assert_not_called()

    def test_recipes_roundtrip_and_old_absence_preserves_automatic_defaults(self):
        recipe = wf.create_recipe(self.model, .5, max_chord_error=.0025, max_relative_volume_error=.005,
                                  curvature_points_per_turn=32, initial_topology_quality=False)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "recipe.json"
            wf.save_recipe(self.model, recipe, path)
            self.assertEqual(wf.load_recipe(self.model, path), recipe)
            old = copy.deepcopy(recipe)
            for field in (*cad._FIDELITY_FIELDS, *cad._CURVATURE_FIELDS, "initial_topology_quality"):
                del old["options"][field]
            path.write_text(json.dumps(old), encoding="utf-8")
            loaded = wf.load_recipe(self.model, path)
            self.assertTrue(all(loaded["options"][field] is None for field in
                                (*cad._FIDELITY_FIELDS, *cad._CURVATURE_FIELDS, "initial_topology_quality")))

    def test_every_retry_retains_fidelity_when_global_size_changes(self):
        success = self.native.mesh.return_value
        self.native.mesh.side_effect = [synthetic_refusal("native-plc-cavity-material-refusal"),
                                       synthetic_refusal("native-plc-final-centroid-outside"), success]
        recipe = wf.create_recipe(self.model, .5, max_chord_error=.0025,
            max_relative_volume_error=.005, curvature_points_per_turn=32, initial_topology_quality=False,
            allow_advanced=True, finer_factors=[.5], max_attempts=3)
        result = wf.run_recipe(self.model, recipe)
        self.assertEqual([call.args[0] for call in self.native.mesh.call_args_list], [.5, .5, .25])
        for call, attempt in zip(self.native.mesh.call_args_list, result.report["attempts"]):
            for options in (call.kwargs, attempt["options"]):
                self.assertEqual(options["max_chord_error"], .0025)
                self.assertEqual(options["max_relative_volume_error"], .005)
                self.assertEqual(options["curvature_points_per_turn"], 32)
                self.assertIs(options["initial_topology_quality"], False)

    def test_cli_forwarding_and_invalid_values_do_not_import(self):
        with patch.object(cli, "read_step", return_value=self.model), \
             patch.object(cli, "_native_output_to_stderr", contextlib.nullcontext), \
             contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(cli.main(["source.step", "out.vtk", "--size", ".5",
                "--max-chord-error", ".0025", "--max-relative-volume-error", ".005"]), 0)
        self.assertEqual(self.native.mesh.call_args.kwargs["max_chord_error"], .0025)
        self.assertEqual(self.native.mesh.call_args.kwargs["max_relative_volume_error"], .005)
        for flag, value in (("--max-chord-error", "0"), ("--max-chord-error", "nan"),
                            ("--max-relative-volume-error", "1"), ("--max-relative-volume-error", "-1")):
            with self.subTest(flag=flag, value=value), patch.object(cli, "read_step") as read, \
                 contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                cli.main(["source.step", "out.vtk", "--size", ".5", flag, value])
            read.assert_not_called()

    def test_discovery_reports_defaults_without_claiming_achieved_fidelity(self):
        controls = pm.capabilities()["controls"]
        self.assertEqual(controls["max_chord_error"]["units"], "mm")
        self.assertIsNone(controls["max_relative_volume_error"]["default"])
        self.assertEqual(controls["max_relative_volume_error"]["native_default"], .01)
        self.assertEqual(controls["curvature_points_per_turn"]["native_default"], 0)

    def test_curvature_integer_bounds_default_and_cli_forwarding(self):
        for value in (None, 0, 32, 64):
            self.model.mesh(.5, curvature_points_per_turn=value)
            self.assertEqual(self.native.mesh.call_args.kwargs["curvature_points_per_turn"], value)
        for value in (-1, 65, True, 1.5, "32"):
            with self.subTest(value=value), self.assertRaises((TypeError, ValueError)):
                self.model.mesh(.5, curvature_points_per_turn=value)
            with self.subTest(recipe=value), self.assertRaises(wf.WorkflowError):
                wf.create_recipe(self.model, .5, curvature_points_per_turn=value)
        with patch.object(cli, "read_step", return_value=self.model), \
             patch.object(cli, "_native_output_to_stderr", contextlib.nullcontext), \
             contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(cli.main(["source.step", "out.vtk", "--size", ".5",
                                       "--curvature-points-per-turn", "32"]), 0)
        self.assertEqual(self.native.mesh.call_args.kwargs["curvature_points_per_turn"], 32)

    def test_topology_toggle_is_nullable_and_does_not_force_a_front(self):
        for value in (None, True, False):
            for policy in ("normal", "bounded"):
                self.model.mesh(.5, resource_policy=policy, initial_topology_quality=value)
                options = self.native.mesh.call_args.kwargs
                self.assertIs(options["initial_topology_quality"], value)
                self.assertIsNone(options["advancing_interior_front"])
        for value in (0, 1, "false", [], {}):
            with self.subTest(value=value), self.assertRaises((TypeError, ValueError)):
                self.model.mesh(.5, initial_topology_quality=value)
            with self.subTest(recipe=value), self.assertRaises(wf.WorkflowError):
                wf.create_recipe(self.model, .5, initial_topology_quality=value)
        for flag, expected in (("--initial-topology-quality", True), ("--no-initial-topology-quality", False)):
            with patch.object(cli, "read_step", return_value=self.model), \
                 patch.object(cli, "_native_output_to_stderr", contextlib.nullcontext), \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cli.main(["source.step", "out.vtk", "--size", ".5", flag]), 0)
            self.assertIs(self.native.mesh.call_args.kwargs["initial_topology_quality"], expected)
        capability = pm.capabilities()["controls"]["initial_topology_quality"]
        self.assertEqual(capability["policy_defaults"], {"normal": True, "bounded": False})
        self.assertFalse(capability["forces_initial_front"])

    def test_route_and_fidelity_receipts_do_not_recertify_native_validity(self):
        report = synthetic_report()
        route = {"reason": "verified_whole_body_primitive", "candidate": "analytic",
                 "preflight_attempted": True, "whole_body_verified": True, "preflight_work": 20,
                 "attempts": [{"method": "legacy_dispatch", "work": 30, "published": False,
                               "reason": "geometry_error", "seconds": .1},
                              {"method": "native_front", "work": 40, "published": True,
                               "reason": "none", "seconds": .2}]}
        report["route_resolution"] = report["bodies"][0]["route_resolution"] = route
        report["fidelity"] = {"requested_max_chord_error": .0025,
            "requested_max_relative_volume_error": .005, "effective_max_relative_volume_error": .005}
        before = copy.deepcopy(report)
        diagnostics = pm.diagnose(report)
        self.assertEqual(diagnostics["bodies"][0]["route_resolution"], route)
        self.assertEqual(diagnostics["summary"]["fidelity"], report["fidelity"])
        self.assertEqual(diagnostics["status"], "valid")
        self.assertEqual(report, before)
        report["audit"]["valid"] = False
        self.assertEqual(pm.diagnose(report)["status"], "failed")
        for invalid in ({"whole_body_verified": "true"}, {"preflight_work": -1},
                        {"attempts": [{}, {}, {}]}, {"attempts": [{"published": 1}]}):
            report["bodies"][0]["route_resolution"] = invalid
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                pm.diagnose(report)
