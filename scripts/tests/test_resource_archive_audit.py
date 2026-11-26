"""Synthetic/stub guards for corpus archive screening; no meshing evidence."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import copy
from contextlib import redirect_stderr, redirect_stdout
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("resource_archive_audit_tests", ROOT / "benchmarks/audit_resource_corpus.py")
AUDIT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = AUDIT
SPEC.loader.exec_module(AUDIT)


def native_report():
    return {"schema": AUDIT.NATIVE_SCHEMA, "status": "succeeded", "reason": "none",
            "mesh_published": True, "archive_published": True, "audit_valid": True,
            "nodes": 4, "cells": 1, "solids": 1, "validated_parts": 1,
            "boundary_missing": 0, "boundary_unexpected": 0,
            "nonmanifold_faces": 0, "inconsistent_face_orientation": 0,
            "mesh_scope": "independent_parts", "interbody_nonoverlap_certified": False,
            "bodies": [{"solid_tag": 1, "reason": "none", "route": "synthetic_stub"}],
            "interferences": []}


def quality_summary():
    metrics = {}
    for name in AUDIT.METRICS:
        value = 30. if name == "minimum_dihedral_deg" else .5
        metrics[name] = {"count": 1, "unavailable": 0,
                         **{key: value for key in ("minimum", "maximum", "mean", "median", "p01", "p05", "p95")}}
    return {"schema": "paacs-mesh-quality-summary-v1", "scope": "entire_snapshot",
            "node_count": 4, "cell_count": 1, "cell_counts": {"tet4": 1},
            "orientation_counts": {"positive": 1, "inverted": 0, "degenerate": 0, "unsupported": 0},
            "range_limited_cells": 0, "quality": metrics, "body_tet4_counts": {"1": 1},
            "unclassified_tet4_count": 0, "metric_conventions": {"fixture": "synthetic test metadata"}}


def python_report():
    report = native_report()
    for key in ("schema", "archive_published", "solids", *AUDIT.AUDIT_COUNTERS):
        del report[key]
    report.update(source_sha256="a" * 64, units="mm", kernel_version="8.0.1", body_tags=[1],
                  tetrahedra=1,
                  audit={"complete": True, "valid": True, "missing_boundary": 0, "unexpected_boundary": 0,
                         "nonmanifold_faces": 0, "inconsistent_orientation": 0, "invalid_connectivity": 0,
                         "nonfinite_nodes": 0},
                  quality={"complete": True, "valid": True, "inverted_cells": 0})
    return report


class StubInspection:
    def __init__(self, path, summary=None, **limits):
        self.mesh_sha256 = hashlib.sha256(Path(path).read_bytes()).hexdigest()
        self.max_seconds = limits["max_seconds"]
        self._summary = quality_summary() if summary is None else summary
        self.validation = {"basis": "unavailable", "archive_structure_checked": True,
                           "independent_native_reaudit_performed": False}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def summary(self, *, spatial_bins):
        if spatial_bins != 1:
            raise ValueError("Unexpected expensive spatial summary")
        return copy.deepcopy(self._summary)

    def query(self, **options):
        if options["associations"] is not False or options["limit"] > 1000:
            raise ValueError("Query exceeded intended audit scope")
        return {"total": 0, "cells": [], "next_offset": None}


class ResourceArchiveAuditTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="paacs-resource-audit-test-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.run = self.root / "run"
        self.run.mkdir()
        self.archive = self.run / "mesh.paacsmesh"
        self.archive.write_bytes(b"Synthetic archive bytes; inspected only by a unit-test stub")
        self.source = {"path": "00000001/fixture.step", "bytes": 12, "sha256": "a" * 64}
        self.row = {"path": self.source["path"], "input_sha256": self.source["sha256"],
                    "input_unchanged": True, "repeat": 0, "outcome": "succeeded",
                    "wall_seconds": 2.5, "output": "mesh.paacsmesh", "output_exists": True,
                    "output_sha256": hashlib.sha256(self.archive.read_bytes()).hexdigest(),
                    "result": native_report()}
        self.inspection = Mock(side_effect=StubInspection)
        self.limits = AUDIT.Limits(max_bytes=64 * 1024**2, max_cells=100, max_nodes=100,
                                   max_work=10000, max_seconds=10, max_total_seconds=20)

    def check_row(self, row=None, expected=None, inspection=None):
        return AUDIT.audit_record(self.run, self.source, row or self.row, expected,
            self.limits, inspection or self.inspection)

    def write_inputs(self, rows=None, repeats=1):
        manifest = {"schema": "paacs-cad-corpus-v1", "files": [self.source], "settings": {"repeats": repeats}}
        (self.run / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        (self.run / "results.jsonl").write_text("".join(json.dumps(row) + "\n" for row in (rows or [self.row])), encoding="utf-8")

    def write_trial(self, report=None):
        report = python_report() if report is None else report
        report_path = self.run / "native-report.json"
        report_path.write_text(json.dumps(report), encoding="utf-8")
        case = {"schema": "paacs-construction-budget-case-v1", "status": "succeeded", "native_status": "succeeded",
                "source": "fixture.step", "source_sha256": "a" * 64, "native_report_file": "native-report.json",
                "native_report_sha256": hashlib.sha256(report_path.read_bytes()).hexdigest(),
                "mesh_sha256": self.row["output_sha256"], "call_seconds": 2.5, "repetition": 1}
        (self.run / "case.json").write_text(json.dumps(case), encoding="utf-8")
        return case

    def test_changed_archive_hash_is_refused_before_inspection(self):
        self.archive.write_bytes(b"changed")
        result = self.check_row()
        self.assertEqual((result["status"], result["reason"]), ("failed", "archive_hash"))
        self.inspection.assert_not_called()

    def test_missing_native_report_cannot_qualify_an_archive(self):
        self.row["result"] = None
        result = self.check_row()
        self.assertEqual((result["status"], result["reason"]), ("failed", "missing_native_report"))
        self.inspection.assert_not_called()

    def test_snapshot_hash_is_rechecked_and_all_limits_forwarded(self):
        def changed_snapshot(path, **limits):
            inspected = StubInspection(path, **limits)
            inspected.mesh_sha256 = "b" * 64
            return inspected
        result = self.check_row(inspection=Mock(side_effect=changed_snapshot))
        self.assertEqual(result["reason"], "archive_hash")
        result = self.check_row()
        self.assertEqual(result["status"], "passed")
        kwargs = self.inspection.call_args.kwargs
        for name in ("max_bytes", "max_cells", "max_nodes", "max_work"):
            self.assertEqual(kwargs[name], getattr(self.limits, name))
        self.assertLessEqual(kwargs["max_seconds"], self.limits.max_seconds)
        self.assertTrue(result["archive_screen_performed"])
        self.assertFalse(result["independent_native_reaudit_performed"])
        self.assertEqual(result["saved_native_checks"]["basis"], AUDIT.SAVED_BASIS)
        self.assertEqual(result["mesh_wall_seconds"], 2.5)
        self.assertIn("independent_audit_seconds", result)

    def test_whole_snapshot_aggregate_failures_are_retained(self):
        variants = []
        inverted = quality_summary()
        inverted["orientation_counts"].update(positive=0, inverted=1)
        variants.append((inverted, "nonpositive_orientation"))
        degenerate = quality_summary()
        degenerate["orientation_counts"].update(positive=0, degenerate=1)
        variants.append((degenerate, "nonpositive_orientation"))
        unavailable = quality_summary()
        unavailable["quality"]["mean_ratio"].update(count=0, unavailable=1)
        variants.append((unavailable, "metric_unavailable"))
        zero = quality_summary()
        zero["quality"]["scaled_jacobian"]["minimum"] = 0
        variants.append((zero, "metric_nonpositive_or_out_of_range"))
        nonfinite = quality_summary()
        nonfinite["quality"]["minimum_dihedral_deg"]["mean"] = float("nan")
        variants.append((nonfinite, "metric_nonfinite"))
        limited = quality_summary()
        limited["range_limited_cells"] = 1
        variants.append((limited, "metric_range_loss"))
        unclassified = quality_summary()
        unclassified["unclassified_tet4_count"] = 1
        variants.append((unclassified, "body_classification"))
        wrong_body = quality_summary()
        wrong_body["body_tet4_counts"] = {"2": 1}
        variants.append((wrong_body, "body_classification"))
        for summary, reason in variants:
            with self.subTest(reason=reason):
                inspection = Mock(side_effect=lambda path, **kw: StubInspection(path, summary, **kw))
                result = self.check_row(inspection=inspection)
                self.assertEqual((result["status"], result["reason"]), ("failed", reason))

    def test_known_original_body_count_is_required_without_claiming_geometry(self):
        expected = {**self.source, "expected_original_solids": 2}
        self.assertEqual(self.check_row(expected=expected)["reason"], "original_body_count")
        expected["expected_original_solids"] = 1
        expected["expected_body_tags"] = [1]
        result = self.check_row(expected=expected)
        self.assertEqual(result["status"], "passed")
        ownership = result["body_ownership"]
        self.assertTrue(ownership["original_count_checked"])
        self.assertTrue(ownership["original_ids_checked"])
        self.assertFalse(ownership["original_body_geometry_reclassification_performed"])
        expected["expected_body_tags"] = [7]
        self.assertEqual(self.check_row(expected=expected)["reason"], "original_body_tags")

    def test_resource_refusal_is_not_silently_counted_as_a_pass(self):
        self.inspection.side_effect = AUDIT.AuditFailure("memory_budget", "Synthetic API refusal")
        result = self.check_row()
        self.assertEqual((result["status"], result["reason"]), ("resource_limit", "memory_budget"))

    def test_native_failure_retained_and_unexpected_output_rejected(self):
        self.row["outcome"] = "resource-limit"
        self.assertEqual(self.check_row()["reason"], "publication_violation")
        self.archive.unlink()
        self.row.update(output_exists=False, output_sha256=None, result=None)
        result = self.check_row()
        self.assertEqual(result["status"], "not_published")
        self.assertEqual(result["reason"], "resource-limit")
        self.assertFalse(result["archive_screen_performed"])
        self.inspection.assert_not_called()

    def test_path_escape_is_refused(self):
        for output in ("../elsewhere.paacsmesh", "C:/elsewhere.paacsmesh", "/elsewhere.paacsmesh"):
            with self.subTest(output=output):
                self.row["output"] = output
                self.assertEqual(self.check_row()["reason"], "path_invalid")
        self.inspection.assert_not_called()

    def test_atomic_complete_receipt_retains_missing_and_duplicate_attempts(self):
        self.write_inputs(rows=[self.row, self.row], repeats=2)
        output = self.root / "audit.json"
        result = AUDIT.run_audit(self.run, output, limits=self.limits, inspect_mesh=self.inspection)
        self.assertEqual(result["status"], "incomplete")
        self.assertEqual(result["counts"], {"passed": 1, "failed": 2})
        self.assertEqual({row.get("reason") for row in result["cases"]}, {None, "missing_result", "duplicate_result"})
        self.assertEqual(json.loads(output.read_text()), result)
        held = output.read_bytes()
        with self.assertRaises(AUDIT.AuditFailure) as caught:
            AUDIT.run_audit(self.run, output, limits=self.limits, inspect_mesh=self.inspection)
        self.assertEqual(caught.exception.reason, "output_exists")
        self.assertEqual(output.read_bytes(), held)

    def test_progress_flushes_every_attempt_without_hiding_refusals_or_duplicates(self):
        class ProgressStream(io.StringIO):
            flushes = 0

            def flush(self):
                self.flushes += 1
                super().flush()

        refused = {**self.row, "repeat": 1, "outcome": "resource-limit", "output": "absent.paacsmesh",
                   "output_exists": False, "output_sha256": None, "result": None}
        self.write_inputs(rows=[self.row, refused, self.row], repeats=2)
        progress = ProgressStream()
        output = self.root / "audit.json"
        result = AUDIT.run_audit(self.run, output, limits=self.limits,
                                 inspect_mesh=self.inspection, progress=progress)
        events = [json.loads(line) for line in progress.getvalue().splitlines()]
        self.assertEqual(len(events), 6)
        self.assertEqual(progress.flushes, 6)
        self.assertEqual([event["event"] for event in events], ["start", "end"] * 3)
        self.assertEqual([event["status"] for event in events[1::2]], ["passed", "not_published", "failed"])
        self.assertEqual(events[-1]["reason"], "duplicate_result")
        for record, pair in enumerate(zip(events[::2], events[1::2]), 1):
            start, end = pair
            for event in pair:
                self.assertEqual(event["schema"], AUDIT.PROGRESS_SCHEMA)
                self.assertEqual(event["record"], record)
                self.assertEqual(event["source_sha256"], self.source["sha256"])
                self.assertEqual(event["path"], self.source["path"])
                self.assertNotIn("worst_cells", event)
            self.assertGreaterEqual(end["independent_audit_seconds"], start["independent_audit_seconds"])
        self.assertEqual(result["counts"], {"passed": 1, "not_published": 1, "failed": 1})
        self.assertEqual(json.loads(output.read_text()), result)

    def test_progress_bounds_identifiers_and_details(self):
        progress = io.StringIO()
        source = {**self.source, "path": "\u2603" * 100_000}
        AUDIT._progress(progress, "end", source, 0, 1, AUDIT.time.perf_counter(),
                        {"status": "failed", "reason": "x" * 100_000})
        event = json.loads(progress.getvalue())
        self.assertTrue(event["path_truncated"])
        self.assertEqual(len(event["path"]), 1024)
        self.assertEqual(len(event["reason"]), 128)
        self.assertLess(len(progress.getvalue().encode("utf-8")), 8192)

    def test_cli_progress_stays_on_stderr_with_one_atomic_result_on_stdout(self):
        self.write_inputs()
        output = self.root / "audit.json"
        stdout, stderr = io.StringIO(), io.StringIO()
        real_audit = AUDIT.run_audit

        def inspect_stub(*args, **options):
            return real_audit(*args, **options, inspect_mesh=self.inspection)

        with patch.object(AUDIT, "run_audit", side_effect=inspect_stub), \
                redirect_stdout(stdout), redirect_stderr(stderr):
            status = AUDIT.main([str(self.run), "--output", str(output), "--progress"])
        self.assertEqual(status, 0)
        self.assertEqual(len(stdout.getvalue().splitlines()), 1)
        self.assertEqual(json.loads(stdout.getvalue())["status"], "complete")
        self.assertEqual([json.loads(line)["event"] for line in stderr.getvalue().splitlines()], ["start", "end"])
        self.assertEqual(json.loads(output.read_text())["counts"], {"passed": 1})

    def test_plan_source_hash_mismatch_refuses_before_inspection(self):
        self.write_inputs()
        plan = self.root / "plan.json"
        plan.write_text(json.dumps({"files": [{**self.source, "sha256": "b" * 64}]}), encoding="utf-8")
        output = self.root / "audit.json"
        with self.assertRaises(AUDIT.AuditFailure) as caught:
            AUDIT.run_audit(self.run, output, plan=plan, limits=self.limits, inspect_mesh=self.inspection)
        self.assertEqual(caught.exception.reason, "source_hash")
        self.assertFalse(output.exists())
        self.inspection.assert_not_called()

    def test_explicit_trial_keeps_python_checks_saved_and_archive_checks_independent(self):
        self.write_trial()
        result = AUDIT.audit_trial(self.run, self.root / "trial-audit.json", limits=self.limits,
                                  inspect_mesh=self.inspection, expected_original_solids=1,
                                  expected_body_tags=[1])
        self.assertEqual(result["mode"], "construction_trial")
        self.assertEqual(result["status"], "complete")
        record = result["cases"][0]
        self.assertEqual(record["saved_native_checks"]["saved_report_format"], "python")
        self.assertEqual(record["saved_native_checks"]["basis"], "saved_native_python_report_not_recertified")
        self.assertFalse(record["independent_native_reaudit_performed"])
        self.assertTrue(record["archive_structure"]["archive_structure_checked"])
        self.assertTrue(record["body_ownership"]["original_count_checked"])

    def test_trial_report_substitution_is_refused_before_inspection(self):
        self.write_trial()
        (self.run / "native-report.json").write_text(json.dumps({"status": "substituted"}), encoding="utf-8")
        output = self.root / "trial-audit.json"
        with self.assertRaises(AUDIT.AuditFailure) as caught:
            AUDIT.audit_trial(self.run, output, limits=self.limits, inspect_mesh=self.inspection)
        self.assertEqual(caught.exception.reason, "report_hash")
        self.assertFalse(output.exists())
        self.inspection.assert_not_called()

    def test_trial_archive_substitution_cannot_pass(self):
        self.write_trial()
        self.archive.write_bytes(b"substituted after the saved trial")
        result = AUDIT.audit_trial(self.run, self.root / "trial-audit.json", limits=self.limits,
                                  inspect_mesh=self.inspection)
        self.assertEqual(result["status"], "incomplete")
        self.assertEqual(result["cases"][0]["reason"], "archive_hash")
        self.inspection.assert_not_called()

    def test_trial_source_fingerprint_and_native_checks_cannot_be_substituted(self):
        report = python_report()
        report["source_sha256"] = "b" * 64
        self.write_trial(report)
        with self.assertRaises(AUDIT.AuditFailure) as caught:
            AUDIT.audit_trial(self.run, self.root / "trial-audit.json", limits=self.limits,
                              inspect_mesh=self.inspection)
        self.assertEqual(caught.exception.reason, "source_hash")
        report = python_report()
        report["audit"]["complete"] = False
        self.write_trial(report)
        with self.assertRaises(AUDIT.AuditFailure) as caught:
            AUDIT.audit_trial(self.run, self.root / "trial-audit.json", limits=self.limits,
                              inspect_mesh=self.inspection)
        self.assertEqual(caught.exception.reason, "native_audit_not_valid")
        self.inspection.assert_not_called()

    def test_trial_native_cancellation_is_retained_without_an_archive_certificate(self):
        report = python_report()
        report.update(status="refused", reason="cancelled", mesh_published=False)
        report["audit"].update(complete=False, valid=False)
        report["quality"].update(complete=False, valid=False)
        case = self.write_trial(report)
        case.update(status="failed", native_status="refused", mesh_sha256=None)
        (self.run / "case.json").write_text(json.dumps(case), encoding="utf-8")
        self.archive.unlink()
        result = AUDIT.audit_trial(self.run, self.root / "trial-audit.json", limits=self.limits,
                                  inspect_mesh=self.inspection)
        self.assertEqual(result["status"], "complete")
        self.assertFalse(result["all_planned_meshes_passed"])
        self.assertEqual(result["independent_archive_passes"], 0)
        record = result["cases"][0]
        self.assertEqual((record["status"], record["reason"]), ("not_published", "cancelled"))
        self.assertFalse(record["archive_screen_performed"])
        self.inspection.assert_not_called()


if __name__ == "__main__":
    unittest.main()
