"""Interpret reported Native summaries; these tests perform no meshing."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import copy
import json
import unittest

import paacs_mesh as pm

from test_construction_diagnostics import synthetic_report


def thin_plate_report():
    """Retained owned thin-plate .25 Native measurements, without mesh data."""
    report = synthetic_report()
    report.update(source_sha256="5cf8ac4351360f0fa25d6be84b25a01fed3d071c3c1753b89730450f795c6362",
                  nodes=120, cells=548, tetrahedra=312, quality_warning=False)
    report["quality"].update(metric="subtet_mean_ratio", minimum=7.93372169838619e-05,
        percentile_1=.0005074100815680177, percentile_5=.007416299214871934,
        median=.0558126048202823, mean=.05281869667154538, warning_threshold=0,
        minimum_required=0, warning_enabled=False)
    return report


class ShapeReviewTests(unittest.TestCase):
    def test_owned_thin_plate_requires_review_with_configured_warnings_disabled(self):
        report = thin_plate_report()
        before = copy.deepcopy(report)
        result = pm.diagnose(report)
        self.assertEqual(result["status"], "valid")
        self.assertTrue(result["validation"]["native_audit_valid"])
        self.assertTrue(result["validation"]["native_quality_valid"])
        self.assertEqual(result["shape_review"]["status"], "recommended")
        self.assertEqual(result["shape_review"]["solver_adequacy"], "unassessed")
        self.assertEqual({i["code"] for i in result["issues"]},
                         {"near_flat_shape_review", "low_shape_distribution_review"})
        for issue in result["issues"]:
            self.assertEqual(issue["severity"], "info")
            self.assertEqual(issue["evidence"]["metric"], "subtet_mean_ratio")
            self.assertEqual(issue["evidence"]["minimum"], before["quality"]["minimum"])
            self.assertFalse(issue["evidence"]["native_validity_changed"])
            self.assertEqual(issue["cell_ids"], [])
            self.assertFalse(any("refine" in action for action in issue["candidate_actions"]))
        self.assertNotIn("preferred_quality_not_met", [i["code"] for i in result["issues"]])
        self.assertEqual(report, before)
        self.assertIn("Shape review: recommended; solver adequacy unassessed", pm.format_diagnostics(result))

    def test_good_native_shape_has_no_automatic_trigger_or_solver_approval(self):
        report = synthetic_report()
        report["quality"].update(metric="subtet_mean_ratio", minimum=.55, percentile_5=.65, median=.83)
        result = pm.diagnose(report)
        self.assertEqual(result["issues"], [])
        self.assertEqual(result["shape_review"]["status"], "no_automatic_trigger")
        self.assertEqual(result["shape_review"]["solver_adequacy"], "unassessed")

    def test_low_native_tail_is_reviewed_without_a_near_flat_minimum(self):
        report = synthetic_report()
        report["quality"].update(metric="subtet_mean_ratio", minimum=.0212261737804,
                                 percentile_5=.0337287874208, median=.113512167539)
        result = pm.diagnose(report)
        self.assertEqual([i["code"] for i in result["issues"]], ["low_shape_distribution_review"])
        self.assertEqual(result["status"], "valid")

    def test_missing_or_unknown_metric_does_not_get_an_invented_review(self):
        for metric in (None, "unknown_convention"):
            report = thin_plate_report()
            report["quality"]["metric"] = metric
            result = pm.diagnose(report)
            self.assertEqual(result["shape_review"]["status"], "unavailable")
            self.assertEqual(result["issues"], [])

    def test_review_is_bounded_and_cannot_hide_a_native_failure(self):
        report = thin_plate_report()
        result = pm.diagnose(report, max_issues=1)
        self.assertEqual(len(result["issues"]), 1)
        self.assertEqual(result["issue_count"], 2)
        self.assertTrue(result["issues_truncated"])
        self.assertLess(len(json.dumps(result, allow_nan=False)), 16384)
        report.update(status="refused", reason="geometry_error", mesh_published=False)
        report["audit"]["valid"] = False
        result = pm.diagnose(report, max_issues=1)
        self.assertEqual(result["status"], "failed")
        self.assertFalse(result["validation"]["native_audit_valid"])
        self.assertEqual(result["shape_review"]["status"], "recommended")

    def test_review_policy_is_discoverable_without_altering_quality_controls(self):
        policy = pm.capabilities()["diagnostics"]["automatic_shape_review"]
        self.assertEqual(policy["metric"], "subtet_mean_ratio")
        self.assertEqual(policy["near_flat_minimum"], .001)
        self.assertFalse(pm.capabilities()["diagnostics"]["shape_review_changes_native_acceptance"])
        report = thin_plate_report()
        report["quality"]["warning_threshold"] = .25
        result = pm.diagnose(report)
        self.assertIn("preferred_quality_not_met", [i["code"] for i in result["issues"]])
        self.assertEqual(result["summary"]["quality"]["warning_threshold"], .25)
        self.assertEqual(result["shape_review"]["status"], "recommended")


if __name__ == "__main__":
    unittest.main()
