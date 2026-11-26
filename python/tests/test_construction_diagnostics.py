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

    def test_original_source_volume_failure_and_unsized_fallback_stay_distinct(self):
        report = synthetic_report()
        attempts = [{"attempt": i + 1, "chord_error_limit": .01 / 2**i,
                     "surface_published": True, "surface_diagnostic": "none",
                     "relative_volume_error": error, "source_model_verified_unchanged": True,
                     "volume_criterion_met": False, "work": 100 + i,
                     "construction_budget_limited": False, "construction_stalled": True}
                    for i, error in enumerate((.03, .02, .011))]
        body = report["bodies"][0]
        body["initial_front_fallback"] = {"attempted": True, "reason": "geometry_error",
            "surface_diagnostic": "none", "volume_diagnostic": "native-plc-source-volume-error",
            "work": 303, "seconds": .25, "source_shell": {"available": True, "accepted": False,
            "diagnostic": "native-plc-source-volume-error", "maximum_relative_volume_error": .01,
            "attempts": attempts}}
        body["source_shell"] = {"available": True, "accepted": True, "diagnostic": "accepted",
            "attempts": [{"attempt": 1, "surface_published": True, "surface_diagnostic": "none",
                          "source_model_verified_unchanged": True, "volume_criterion_met": True,
                          "relative_volume_error": .005}]}
        body["volume_sizing"] = {"applied": False, "diagnostic": "kernel_fan_has_no_volume_size_field",
            "requested_global_size": .5, "sizing_completion_certified": False,
            "basis": "reported_native_construction"}
        before = copy.deepcopy(report)
        result = pm.diagnose(report, max_issues=1)
        self.assertEqual(result["status"], "valid_with_warnings")
        self.assertTrue(result["validation"]["native_audit_valid"])
        self.assertEqual(result["issue_count"], 2)
        self.assertTrue(result["issues_truncated"])
        self.assertEqual(result["issues"][0]["code"], "volume_sizing_not_applied")
        self.assertFalse(result["issues"][0]["evidence"]["sizing_completion_certified"])
        original = result["bodies"][0]["initial_front_fallback"]
        self.assertEqual(original["volume_diagnostic"], "native-plc-source-volume-error")
        self.assertEqual(original["surface_diagnostic"], "none")
        self.assertEqual(original["source_shell"]["attempts_total"], 3)
        self.assertEqual(len(original["source_shell"]["attempts"]), 1)
        self.assertTrue(original["source_shell"]["attempts_truncated"])
        self.assertEqual(report, before)

    def test_source_geometry_receipt_is_not_a_sizing_or_solver_certificate(self):
        report = synthetic_report()
        body = report["bodies"][0]
        body["source_shell"] = {"available": True, "accepted": True, "diagnostic": "accepted",
            "attempts": [{"attempt": 1, "surface_published": True, "surface_diagnostic": "none",
                          "source_model_verified_unchanged": True, "volume_criterion_met": True,
                          "relative_volume_error": .005}]}
        body["volume_sizing"] = {"applied": True, "diagnostic": "sized_interior_construction",
                                 "sizing_completion_certified": False}
        result = pm.diagnose(report)
        self.assertEqual(result["status"], "valid")
        self.assertEqual(result["issues"], [])
        self.assertFalse(result["bodies"][0]["volume_sizing"]["sizing_completion_certified"])
        self.assertEqual(result["shape_review"]["solver_adequacy"], "unassessed")

    def test_malformed_source_and_sizing_receipts_are_rejected(self):
        for changed in ({"source_shell": {"available": "yes"}},
                        {"source_shell": {"accepted": 1}},
                        {"source_shell": {"attempts": [{}, {}, {}, {}]}},
                        {"source_shell": {"attempts": [{"surface_published": "true"}]}},
                        {"volume_sizing": {"applied": "true"}},
                        {"initial_front_fallback": {"source_shell": {"attempts": [{}, {}, {}, {}]}}},
                        {"initial_front_fallback": {"volume_diagnostic": 1}}):
            with self.subTest(changed=changed):
                report = synthetic_report()
                report["bodies"][0].update(changed)
                with self.assertRaises(ValueError):
                    pm.diagnose(report)

    def test_retained_initial_topology_is_copied_without_recertifying_geometry(self):
        report = synthetic_report()
        receipt = {"diagnostic": "complete", "work": 40,
                   "initial_topology_used": True, "initial_topology_diagnostic": "imported",
                   "initial_topology_cells": 12, "initial_topology_import_work": 24,
                   "initial_topology_peak_modeled_bytes": 4096}
        report["bodies"][0]["constraint_recovery"] = receipt
        before = copy.deepcopy(report)
        result = pm.diagnose(report)
        self.assertEqual(result["bodies"][0]["constraint_recovery"], receipt)
        self.assertEqual(result["status"], "valid")
        self.assertEqual(result["issues"], [])
        self.assertEqual(result["shape_review"]["solver_adequacy"], "unassessed")
        self.assertEqual(report, before)
        report["audit"]["valid"] = False
        report["status"], report["reason"], report["mesh_published"] = "refused", "geometry_error", False
        failed = pm.diagnose(report)
        self.assertEqual(failed["status"], "failed")
        self.assertTrue(failed["bodies"][0]["constraint_recovery"]["initial_topology_used"])

    def test_missing_or_refused_initial_topology_has_no_invented_success(self):
        report = synthetic_report()
        result = pm.diagnose(report)
        self.assertIsNone(result["bodies"][0]["constraint_recovery"])
        receipt = {"initial_topology_used": False, "initial_topology_diagnostic": "import_refused",
                   "initial_topology_cells": 12, "initial_topology_import_work": 24,
                   "initial_topology_peak_modeled_bytes": 4096}
        report["bodies"][0]["constraint_recovery"] = receipt
        self.assertEqual(pm.diagnose(report)["bodies"][0]["constraint_recovery"], receipt)

    def test_malformed_initial_topology_receipts_are_rejected(self):
        for receipt in ([1], {"initial_topology_used": "true"},
                        {"initial_topology_diagnostic": None}, {"initial_topology_cells": -1},
                        {"initial_topology_import_work": True}, {"initial_topology_peak_modeled_bytes": 1.5}):
            with self.subTest(receipt=receipt):
                report = synthetic_report()
                report["bodies"][0]["constraint_recovery"] = receipt
                with self.assertRaises(ValueError):
                    pm.diagnose(report)

    def test_cancelled_source_volume_pass_is_interrupted_validation_not_a_volume_defect(self):
        report = synthetic_report()
        report.update(status="refused", reason="cancelled", mesh_published=False)
        report["bodies"][0]["source_shell"] = {"available": True, "accepted": False,
            "diagnostic": "volume_criterion_met", "maximum_relative_volume_error": .01,
            "attempts": [{"attempt": 1, "surface_published": True, "volume_criterion_met": True,
                          "source_model_verified_unchanged": True, "relative_volume_error": 3.861482873600608e-5}]}
        result = pm.diagnose(report)
        self.assertEqual(result["status"], "failed")
        interruption = next(i for i in result["issues"] if i["code"] == "source_shell_validation_interrupted")
        self.assertEqual(interruption["severity"], "info")
        self.assertEqual(interruption["evidence"]["native_reason"], "cancelled")
        self.assertNotIn("source_shell_not_accepted", [i["code"] for i in result["issues"]])

    def test_partition_stage_evidence_is_preserved_and_output_bounded(self):
        report = synthetic_report()
        body = report["bodies"][0]
        body.update(native_partition_seconds=.3, native_assembly_seconds=.05,
                    partition_reason="none", partition_policy="native_planes",
                    partition_attempts=5, partition_pieces=3, partition_interfaces=2,
                    partition_work=120, partition_volume_change=1e-12,
                    partition_attempt_reports=[{"policy": "native_planes", "diagnostic": "piece_refused",
                        "failed_piece_tag": i, "pieces": 3, "interfaces": 2, "work": 20,
                        "seconds": .1, "surface_deflection": .005} for i in range(1, 6)])
        before = copy.deepcopy(report)
        result = pm.diagnose(report, max_issues=2)
        measured = result["bodies"][0]
        self.assertEqual(measured["native_partition_seconds"], .3)
        self.assertEqual(measured["native_assembly_seconds"], .05)
        self.assertEqual(measured["partition_attempts"], 5)
        self.assertEqual(measured["partition_attempt_reports_total"], 5)
        self.assertEqual(len(measured["partition_attempt_reports"]), 2)
        self.assertTrue(measured["partition_attempt_reports_truncated"])
        self.assertEqual(measured["partition_attempt_reports"][0]["failed_piece_tag"], 1)
        self.assertEqual(result["status"], "valid")
        self.assertEqual(report, before)
        body["partition_attempt_reports"] = [{"work": -1}]
        with self.assertRaises(ValueError):
            pm.diagnose(report)


if __name__ == "__main__":
    unittest.main()
