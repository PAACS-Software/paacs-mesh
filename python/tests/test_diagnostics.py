"""Reported validity is preserved; suggestions never invent local mesh evidence."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import copy
import json
import unittest

import paacs_mesh as pm
from test_native import fixture_path


class DiagnosticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = pm.read_step(fixture_path())
        cls.mesh = cls.model.mesh(.7)

    def test_native_mesh_interpretation_retains_identity_and_evidence(self):
        native = self.mesh.report
        result = pm.diagnose(self.mesh)
        self.assertEqual(result["schema"], "paacs-mesh-diagnostics-v1")
        self.assertIn(result["status"], ("valid", "valid_with_warnings"))
        self.assertEqual(result["source_sha256"], self.model.source_sha256)
        self.assertEqual(result["summary"]["quality"]["minimum"], native["quality"]["minimum"])
        self.assertTrue(result["validation"]["native_audit_complete"])
        self.assertFalse(result["validation"]["interbody_nonoverlap_certified"])
        self.assertEqual(result["bodies"][0]["maximum_boundary_sample_error"],
                         native["bodies"][0]["maximum_boundary_sample_error"])
        json.dumps(result, allow_nan=False)
        result["summary"]["nodes"] = -1
        self.assertNotEqual(pm.diagnose(self.mesh)["summary"]["nodes"], -1)

    def test_quality_warning_does_not_refuse_valid_mesh_or_invent_localization(self):
        native = self.mesh.report
        native["quality_warning"] = True
        native["quality"]["warning_threshold"] = 1
        result = pm.diagnose(native)
        self.assertEqual(result["status"], "valid_with_warnings")
        issue = next(item for item in result["issues"] if item["code"] == "preferred_quality_not_met")
        self.assertEqual(issue["cad_face_tags"], [])
        self.assertEqual(issue["cell_ids"], [])
        self.assertEqual(issue["body_tags"], [])
        self.assertIn("quality_to_CAD_face_association", result["unavailable"])
        self.assertIn("warning: preferred_quality_not_met", pm.format_diagnostics(result))
        self.assertEqual(issue["id"], pm.diagnose(copy.deepcopy(native))["issues"][0]["id"])

    def test_assembly_overlap_keeps_part_validity_and_reports_exact_body_pair(self):
        native = self.mesh.report
        native["interferences"] = [{"first_solid_tag": 3, "second_solid_tag": 2,
                                     "state": "material_overlap", "overlap_volume": .125}]
        result = pm.diagnose(native)
        self.assertEqual(result["status"], "valid_with_warnings")
        issue = next(item for item in result["issues"] if item["code"] == "assembly_interference")
        self.assertEqual(issue["body_tags"], [2, 3])
        self.assertEqual(issue["evidence"]["overlap_volume_mm3"], .125)
        self.assertFalse(result["validation"]["interbody_nonoverlap_certified"])

    def test_native_failures_and_missing_checks_cannot_become_success(self):
        for stage in ("audit", "quality"):
            native = self.mesh.report
            native[stage]["complete"] = False
            self.assertNotIn(pm.diagnose(native)["status"], ("valid", "valid_with_warnings"))
        native = self.mesh.report
        native["audit"]["invalid_connectivity"] = 1
        self.assertEqual(pm.diagnose(native)["status"], "failed")
        native = self.mesh.report
        native.update(mesh_published=False, reason="resource_limit")
        native["bodies"][0].update(reason="resource_limit", volume_reason="work_budget")
        result = pm.diagnose(native)
        self.assertEqual(result["status"], "failed")
        body = next(item for item in result["issues"] if item["scope"] == "body")
        self.assertEqual(body["body_tags"], [1])
        self.assertEqual(body["evidence"]["volume_reason"], "work_budget")

    def test_workflow_failure_and_omitted_final_attempt_are_never_hidden(self):
        workflow = {"schema": "paacs-mesh-workflow-v1", "status": "failed",
                    "source_sha256": self.model.source_sha256, "seconds": 3,
                    "attempts": [{"native_report": self.mesh.report}, {"native_report": None}]}
        result = pm.diagnose(workflow)
        self.assertEqual(result["status"], "failed")
        self.assertIsNone(result["summary"]["nodes"])
        self.assertEqual(result["timings"]["attempts"], 2)
        workflow.update(status="succeeded", diagnostics_complete=False)
        workflow["attempts"][-1]["native_report"] = self.mesh.report
        self.assertNotIn(pm.diagnose(workflow)["status"], ("valid", "valid_with_warnings"))

    def test_report_limits_count_hidden_issues_and_do_not_drop_failure_status(self):
        native = self.mesh.report
        native["quality_warning"] = True
        native["bodies"][0]["reason"] = "meshing_failed"
        result = pm.diagnose(native, max_issues=1)
        self.assertEqual(len(result["issues"]), 1)
        self.assertTrue(result["issues_truncated"])
        self.assertGreater(result["issue_count"], 1)
        self.assertEqual(result["status"], "failed")
        for maximum in (True, 0, 257):
            with self.assertRaises(ValueError):
                pm.diagnose(native, max_issues=maximum)

    def test_untrusted_report_shape_and_nonfinite_values_are_rejected(self):
        for value in ({}, {"status": "success"}, [], {"mesh_published": "true"}):
            with self.assertRaises(ValueError):
                pm.diagnose(value)
        native = self.mesh.report
        native["quality"]["minimum"] = float("nan")
        with self.assertRaises(ValueError):
            pm.diagnose(native)
        for count in (-1, True, "0"):
            native = self.mesh.report
            native["audit"]["nonfinite_nodes"] = count
            with self.assertRaises(ValueError):
                pm.diagnose(native)
        for key in ("reason", "units", "source_sha256", "status"):
            native = self.mesh.report
            del native[key]
            self.assertNotIn(pm.diagnose(native)["status"], ("valid", "valid_with_warnings"))
        workflow = {"schema": "paacs-mesh-workflow-v1", "status": "succeeded",
                    "source_sha256": "0" * 64, "attempts": [{"native_report": self.mesh.report}]}
        with self.assertRaisesRegex(ValueError, "different CAD sources"):
            pm.diagnose(workflow)

    def test_capability_discovery_exposes_limitations_without_numpy(self):
        result = pm.capabilities()
        self.assertEqual(result["schema"], "paacs-mesh-capabilities-v1")
        self.assertEqual(result["controls"]["spatial_sizing_requires"]["solids"], 1)
        self.assertTrue(result["diagnostics"]["local_cell_quality"])
        self.assertEqual(result["diagnostics"]["local_cell_quality_requires"], "opt_in_inspect_mesh")
        self.assertFalse(result["diagnostics"]["diagnose_includes_per_cell_quality"])
        self.assertEqual(result["step_volume_cell_families"], ["Tet4"])
        self.assertEqual(result["assemblies"]["default"], "preserve_parts")
        self.assertIn("diagnose", result["cli_commands"])
        json.dumps(result, allow_nan=False)

    def test_nested_cli_workflow_failure_keeps_final_native_identity_and_failure(self):
        wrapped = {"status": "failed", "stage": "publication", "reason": "io_error", "report": {
            "schema": "paacs-mesh-workflow-v1", "status": "succeeded",
            "source_sha256": self.model.source_sha256, "attempts": [{"native_report": self.mesh.report}]}}
        result = pm.diagnose(wrapped)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["source_sha256"], self.model.source_sha256)
        self.assertEqual(result["summary"]["nodes"], self.mesh.node_count)
        self.assertEqual(result["timings"]["attempts"], 1)

    def test_contradictory_native_claims_cannot_strengthen_validity(self):
        for change in ("assembly", "status", "strict_quality"):
            native = self.mesh.report
            if change == "assembly":
                native["interbody_nonoverlap_certified"] = True
            elif change == "status":
                native["status"] = "refused"
            else:
                native["quality"]["minimum_required"] = 1
            result = pm.diagnose(native)
            self.assertEqual(result["status"], "failed")
            self.assertTrue(result["issues"])
            self.assertFalse(result["validation"]["interbody_nonoverlap_certified"])


if __name__ == "__main__":
    unittest.main()
