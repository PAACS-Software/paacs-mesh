"""Resource policy intent, explicit overrides and historical recipe replay."""

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
from paacs_mesh import __main__ as cli, cad, workflow as wf

from test_construction_budgets import model_stub, synthetic_refusal
from test_construction_diagnostics import synthetic_report


class ResourcePolicyTests(unittest.TestCase):
    """Stub tests verify forwarding only; native acceptance is tested separately."""

    def setUp(self):
        self.model, self.native = model_stub()

    def test_explicit_signatures_defaults_and_omitted_limits_are_preserved(self):
        for function in (cad.CadModel.mesh, wf.create_recipe):
            parameter = inspect.signature(function).parameters["resource_policy"]
            self.assertEqual(parameter.default, cad._DEFAULT_RESOURCE_POLICY)
            self.assertEqual(parameter.kind, inspect.Parameter.KEYWORD_ONLY)
        self.model.mesh(.5)
        options = self.native.mesh.call_args.kwargs
        self.assertEqual(options["resource_policy"], cad._DEFAULT_RESOURCE_POLICY)
        self.assertEqual(options["construction_effort"], "legacy")
        self.assertIsNone(options["max_work"])
        self.assertIsNone(options["max_working_bytes"])
        for policy in ("normal", "bounded"):
            self.model.mesh(.5, resource_policy=policy, max_work=123456,
                            max_working_bytes=987654, max_nodes=123, max_cells=456)
            options = self.native.mesh.call_args.kwargs
            self.assertEqual(options["resource_policy"], policy)
            for name, value in (("max_work", 123456), ("max_working_bytes", 987654),
                                ("max_nodes", 123), ("max_cells", 456)):
                self.assertEqual(options[name], value)

    def test_invalid_policies_do_not_enter_native_meshing(self):
        for value in (None, True, False, 0, 1, "", "Normal", "unlimited", []):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.model.mesh(.5, resource_policy=value)
            with self.subTest(recipe=value), self.assertRaises(wf.WorkflowError) as caught:
                wf.create_recipe(self.model, .5, resource_policy=value)
            self.assertEqual(caught.exception.reason, "invalid_recipe")
        self.native.mesh.assert_not_called()

    def test_recipe_roundtrip_and_missing_policy_keep_literal_bounded_meaning(self):
        recipe = wf.create_recipe(self.model, .5, resource_policy="normal")
        self.assertEqual(recipe["options"]["resource_policy"], "normal")
        self.assertIsNone(recipe["options"]["max_work"])
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "recipe.json"
            wf.save_recipe(self.model, recipe, path)
            self.assertEqual(wf.load_recipe(self.model, path), recipe)
            old = copy.deepcopy(recipe)
            del old["options"]["resource_policy"]
            path.write_text(json.dumps(old), encoding="utf-8")
            # Normalizing old files must not depend on a new-request default.
            with patch.object(cad, "_DEFAULT_RESOURCE_POLICY", "normal"):
                loaded = wf.load_recipe(self.model, path)
            self.assertEqual(loaded["options"]["resource_policy"], "bounded")
            self.assertIsNone(loaded["options"]["max_work"])
            old["options"]["resource_policie"] = "normal"
            path.write_text(json.dumps(old), encoding="utf-8")
            with self.assertRaises(wf.WorkflowError) as caught:
                wf.load_recipe(self.model, path)
            self.assertEqual(caught.exception.reason, "invalid_recipe")

    def test_all_retries_keep_policy_and_explicit_global_limits(self):
        success = self.native.mesh.return_value
        self.native.mesh.side_effect = [synthetic_refusal("native-plc-cavity-material-refusal"),
                                       synthetic_refusal("native-plc-final-centroid-outside"), success]
        recipe = wf.create_recipe(self.model, .5, resource_policy="normal", max_work=345678,
            max_working_bytes=987654, allow_advanced=True, finer_factors=[.5], max_attempts=3)
        result = wf.run_recipe(self.model, recipe)
        self.assertEqual([call.args[0] for call in self.native.mesh.call_args_list], [.5, .5, .25])
        for call, attempt in zip(self.native.mesh.call_args_list, result.report["attempts"]):
            for options in (call.kwargs, attempt["options"]):
                self.assertEqual(options["resource_policy"], "normal")
                self.assertEqual(options["max_work"], 345678)
                self.assertEqual(options["max_working_bytes"], 987654)
        self.assertEqual(recipe["options"]["size"], .5)

    def test_cli_forwards_resource_intent_and_run_has_no_override(self):
        args = ["source.step", "out.vtk", "--size", ".5", "--resource-policy", "normal",
                "--max-work", "345678", "--max-working-bytes", "987654",
                "--max-nodes", "123", "--max-cells", "456"]
        with patch.object(cli, "read_step", return_value=self.model), \
             patch.object(cli, "_native_output_to_stderr", contextlib.nullcontext), \
             contextlib.redirect_stdout(io.StringIO()) as captured:
            self.assertEqual(cli.main(args), 0)
        json.loads(captured.getvalue())
        options = self.native.mesh.call_args.kwargs
        self.assertEqual(options["resource_policy"], "normal")
        for name, value in (("max_work", 345678), ("max_working_bytes", 987654),
                            ("max_nodes", 123), ("max_cells", 456)):
            self.assertEqual(options[name], value)
        for flag, value in (("--resource-policy", "normal"), ("--max-working-bytes", "987654"),
                            ("--max-nodes", "123"), ("--max-cells", "456"), ("--max-work", "345678")):
            with self.subTest(flag=flag), self.assertRaises(ValueError):
                cli._command_parser("run").parse_args([
                    "recipe.json", "--source", "source.step", "--output", "new", flag, value])

    def test_cli_rejects_invalid_limits_before_importing(self):
        for flag in ("--max-work", "--max-working-bytes", "--max-nodes", "--max-cells"):
            for value in ("0", "-1", "1.5", str(2**63)):
                with self.subTest(flag=flag, value=value), patch.object(cli, "read_step") as read, \
                     contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as caught:
                    cli.main(["source.step", "out.vtk", "--size", ".5", flag, value])
                self.assertEqual(caught.exception.code, 2)
                read.assert_not_called()

    def test_discovery_and_diagnostics_preserve_reported_resource_evidence(self):
        capabilities = pm.capabilities()["resources"]
        self.assertEqual(capabilities["default"], cad._DEFAULT_RESOURCE_POLICY)
        self.assertEqual(set(capabilities["policies"]), {"normal", "bounded"})
        self.assertTrue(capabilities["explicit_limits_override_policy"])
        report = synthetic_report()
        resources = {"policy": "normal", "requested_max_work": None,
            "effective_max_work": 2**63 - 1, "logical_work_cap_disabled": True,
            "max_working_bytes": 512 * 1024**2, "scope": "reported native admission"}
        report["resources"] = resources
        diagnosed = pm.diagnose(report)
        self.assertEqual(diagnosed["summary"]["resources"], resources)
        self.assertEqual(diagnosed["validation"]["basis"], "reported_native_checks")
        diagnosed["summary"]["resources"]["policy"] = "bounded"
        self.assertEqual(report["resources"]["policy"], "normal")
        del report["resources"]
        self.assertIsNone(pm.diagnose(report)["summary"]["resources"])


if __name__ == "__main__":
    unittest.main()
