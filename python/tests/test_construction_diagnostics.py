"""Budget evidence is distinct from validity, quality and solver adequacy."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import copy
import unittest

import paacs_mesh as pm


def synthetic_report():
    # Synthetic reported-check metadata tests interpretation, not certification.
    return {"status": "succeeded", "reason": "none", "mesh_published": True,
            "source_sha256": "a" * 64, "units": "mm", "nodes": 4, "cells": 1, "tetrahedra": 1,
            "mesh_scope": "conforming_domain", "interbody_nonoverlap_certified": True,
            "audit": {"complete": True, "valid": True, "missing_boundary": 0,
                      "unexpected_boundary": 0, "nonmanifold_faces": 0, "inconsistent_orientation": 0,
                      "invalid_connectivity": 0, "nonfinite_nodes": 0},
            "quality": {"complete": True, "valid": True, "inverted_cells": 0,
                        "minimum": .5, "warning_threshold": 0, "minimum_required": 0},
            "construction": {"effort": "auto", "max_surface_points": 512},
            "bodies": [{"solid_tag": 7, "reason": "none", "surface_spacing": [], "initial_front": {}}],
            "interferences": []}


class ConstructionDiagnosticTests(unittest.TestCase):
    def test_surface_limit_reports_measured_face_without_refusing_validity(self):
        report = synthetic_report()
        row = {"face_tag": 11, "point_limit": 512, "added_points": 512,
               "pending_entries": 17, "budget_limited": True, "queue_exhausted": False,
               "termination": "point_limit"}
        report["bodies"][0]["surface_spacing"] = [row]
        before = copy.deepcopy(report)
        result = pm.diagnose(report)
        self.assertEqual(result["status"], "valid_with_warnings")
        self.assertTrue(result["validation"]["native_audit_valid"])
        warning = next(issue for issue in result["issues"] if issue["code"] == "construction_budget_limited")
        self.assertEqual(warning["body_tags"], [7])
        self.assertEqual(warning["cad_face_tags"], [11])
        self.assertEqual(warning["cell_ids"], [])
        self.assertFalse(warning["evidence"]["sizing_completion_certified"])
        self.assertEqual(result["summary"]["construction"], report["construction"])
        self.assertFalse(warning["evidence"]["native_validity_gates_changed"])
        self.assertEqual(report, before)

    def test_initial_limit_and_incomplete_source_have_explicit_evidence(self):
        report = synthetic_report()
        report["bodies"][0]["initial_front"] = {"diagnostic": "budget_limited",
            "budget_reason": "candidate_limit", "candidate_limit": 16384, "candidates": 16384,
            "reserved_bulk_candidates": 4096, "source_phase_deferred": True}
        issue = pm.diagnose(report)["issues"][0]
        self.assertEqual(issue["code"], "construction_budget_limited")
        self.assertEqual(issue["evidence"]["initial_stage"]["candidate_limit"], 16384)
        self.assertEqual(issue["cad_face_tags"], [])
        self.assertNotIn("automatically_refine", issue["candidate_actions"])

    def test_blocked_candidates_and_old_reports_are_not_invented_budget_events(self):
        report = synthetic_report()
        report["bodies"][0]["surface_spacing"] = [{"termination": "blocked_candidates",
            "budget_limited": False, "queue_exhausted": False}]
        self.assertEqual(pm.diagnose(report)["status"], "valid")
        report["bodies"][0] = {"solid_tag": 7, "reason": "none"}
        self.assertEqual(pm.diagnose(report)["status"], "valid")

    def test_budget_warning_cannot_hide_a_native_refusal(self):
        report = synthetic_report()
        report.update(status="refused", reason="budget_exceeded", mesh_published=False)
        report["bodies"][0]["initial_front"] = {"diagnostic": "budget_limited", "budget_reason": "work_limit"}
        result = pm.diagnose(report)
        self.assertEqual(result["status"], "failed")
        self.assertIn("mesh_refused", [issue["code"] for issue in result["issues"]])

    def test_malformed_optional_budget_metadata_is_rejected(self):
        for changed in ({"surface_spacing": [1]}, {"surface_spacing": [{"budget_limited": "true"}]},
                        {"initial_front": []}, {"initial_front": {"budget_reason": 42}}):
            with self.subTest(changed=changed):
                report = synthetic_report()
                report["bodies"][0].update(changed)
                with self.assertRaises(ValueError):
                    pm.diagnose(report)

    def test_surface_evidence_and_body_rows_have_explicit_output_bounds(self):
        report = synthetic_report()
        report["bodies"][0]["surface_spacing"] = [
            {"face_tag": tag, "budget_limited": True, "queue_exhausted": False}
            for tag in range(1, 101)]
        result = pm.diagnose(report, max_issues=3)
        issue = result["issues"][0]
        self.assertEqual(issue["cad_face_tags"], [1, 2, 3])
        self.assertTrue(issue["evidence"]["cad_face_tags_truncated"])
        self.assertEqual(len(issue["evidence"]["surface_stages"]), 3)
        self.assertEqual(issue["evidence"]["surface_stages_total"], 100)
        self.assertTrue(issue["evidence"]["surface_stages_truncated"])
        body = result["bodies"][0]
        self.assertEqual(len(body["surface_spacing"]), 3)
        self.assertEqual(body["surface_spacing_total"], 100)
        self.assertTrue(body["surface_spacing_truncated"])
        self.assertEqual(len(report["bodies"][0]["surface_spacing"]), 100)


if __name__ == "__main__":
    unittest.main()
