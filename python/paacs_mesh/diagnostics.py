"""Bounded interpretations of native reports, without rerunning or relaxing audits."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import json

from ._cad_common import integer, json_copy


DIAGNOSTIC_SCHEMA = "paacs-mesh-diagnostics-v1"
CAPABILITY_SCHEMA = "paacs-mesh-capabilities-v1"
_AUDIT_ERRORS = ("missing_boundary", "unexpected_boundary", "nonmanifold_faces",
                 "inconsistent_orientation", "invalid_connectivity", "nonfinite_nodes")
_QUALITY = ("metric", "complete", "valid", "minimum", "percentile_1", "percentile_5",
            "median", "mean", "inverted_cells", "warning_threshold", "minimum_required")
_BODY = ("solid_tag", "route", "reason", "surface_reason", "volume_reason", "cad_volume",
         "mesh_volume", "relative_volume_error", "maximum_boundary_node_error",
         "maximum_boundary_sample_error", "surface_chord_error_limit",
         "boundary_embedding_certified", "centroid_outside_distance_limit", "native_centroids_checked")


def capabilities():
    """Describe implemented operations and restrictions; this is not a CAD guarantee."""
    from . import __version__, native_version

    return {
        "schema": CAPABILITY_SCHEMA, "package_version": __version__, "engine_version": native_version,
        "units": {"length": "mm", "volume": "mm3"},
        "operations": ["inspect", "select_faces", "preview", "mesh", "diagnose", "run_recipe", "publish_result"],
        "cli_commands": ["capabilities", "inspect", "run", "diagnose"],
        "input_formats": ["STEP"], "output_formats": ["VTK", "paacsmesh"],
        "step_volume_cell_families": ["Tet4"],
        "selection": {"predicates": ["tags", "name", "group", "surface", "body", "area", "radius",
                                      "centroid_box", "normal", "adjacent_to"],
                      "ambiguity": "explicit_expectation", "identity": "source_sha256_and_guarded_local_tags",
                      "inspection_page_limit": 1000},
        "controls": {"global_size": True, "face_size": True, "box_size": True,
                     "minimum_size": True, "growth": True,
                     "spatial_sizing_requires": {"solids": 1, "assembly": "conforming", "cell_family": "Tet4"}},
        "recovery": {"default": "standard", "experimental": "advanced",
                     "recipe_retries": "explicit_bounded_diagnostic_selection"},
        "assemblies": {"default": "preserve_parts", "interference": "reported_without_fusing",
                       "independent_parts_certify_global_nonoverlap": False},
        "diagnostics": {"schema": DIAGNOSTIC_SCHEMA, "native_report_summary": True,
                        "body_cad_error": True, "local_cell_quality": False,
                        "face_achieved_size": False, "solver_error_estimation": False,
                        "maximum_issues": 256, "report_budget_bytes": 8 * 1024 * 1024},
        "limitations": ["arbitrary_STEP_geometry_is_not_guaranteed",
                        "shape_quality_is_not_simulation_convergence",
                        "diagnostics_interpret_reported_checks_and_do_not_recertify_meshes",
                        "CAD_tags_are_local_to_the_guarded_source_and_kernel",
                        "native_calls_are_serialized_with_internal_workers",
                        "NumPy_is_optional_for_array_extraction"],
    }


def _mapping(value, name):
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _records(value, name):
    if not isinstance(value, list) or any(not isinstance(item, dict) for item in value):
        raise ValueError(f"{name} must be a list of JSON objects")
    return value


def _source(value):
    if value is None:
        return None
    if (not isinstance(value, str) or len(value) != 64 or
            any(char not in "0123456789abcdef" for char in value)):
        raise ValueError("source_sha256 must be a lowercase SHA-256 or null")
    return value


def _positive(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and value > 0


def _native_report(report):
    """Select the final attempt, never a successful earlier attempt from a failed run."""
    workflow = None
    envelopes = []
    while report.get("status") == "failed" and "report" in report and "mesh_published" not in report:
        if len(envelopes) == 4:
            raise ValueError("Too many nested failure envelopes")
        envelopes.append(report)
        report = _mapping(report["report"], "report")
    if report.get("schema") == "paacs-mesh-workflow-v1":
        workflow = report
        attempts = _records(report.get("attempts", []), "attempts")
        native = report.get("native_report")
        if attempts:
            native = attempts[-1].get("native_report")
        report = {} if native is None else _mapping(native, "native_report")
    if envelopes:
        workflow = dict(workflow or {})
        workflow.update(status="failed", outer_failure={key: envelopes[0].get(key) for key in ("stage", "reason")})
    if "mesh_published" not in report and workflow is None:
        raise ValueError("Expected a native mesh report or paacs-mesh-workflow-v1 report")
    if "mesh_published" in report and not isinstance(report["mesh_published"], bool):
        raise ValueError("mesh_published must be bool")
    return report, workflow


def _validate_measurements(native, audit, quality):
    for record, fields in ((audit, ("complete", "valid")), (quality, ("complete", "valid"))):
        for key in fields:
            if key in record and not isinstance(record[key], bool):
                raise ValueError(f"{key} must be bool")
    for record, fields in ((audit, _AUDIT_ERRORS), (quality, ("inverted_cells",)),
                           (native, ("nodes", "cells", "tetrahedra"))):
        for key in fields:
            if key in record:
                integer(record[key], key, minimum=0)
    for key in ("minimum", "percentile_1", "percentile_5", "median", "mean", "warning_threshold", "minimum_required"):
        value = quality.get(key)
        if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 1):
            raise ValueError(f"quality.{key} must be a number in [0,1] or null")


def diagnose(report_or_mesh, *, max_issues=32):
    """Return owned JSON evidence, stable issue IDs and suggestions, with no retries.

    Accept a native ``Mesh``, its report, a workflow result/report, or a CLI
    failure envelope. Saved reports are interpreted as supplied, not authenticated
    against a mesh file. Missing local evidence is reported explicitly.
    """
    maximum = integer(max_issues, "max_issues")
    if maximum > 256:
        raise ValueError("max_issues must not exceed 256")
    raw = report_or_mesh if isinstance(report_or_mesh, dict) else getattr(report_or_mesh, "report", None)
    raw = json_copy(_mapping(raw, "report"))
    native, workflow = _native_report(raw)
    source = _source(native.get("source_sha256", (workflow or {}).get("source_sha256")))
    if native.get("units", "mm") != "mm":
        raise ValueError("Native mesh reports must use mm")
    audit = _mapping(native.get("audit", {}), "audit")
    quality = _mapping(native.get("quality", {}), "quality")
    _validate_measurements(native, audit, quality)
    workflow_source = _source((workflow or {}).get("source_sha256"))
    if workflow_source is not None and source is not None and workflow_source != source:
        raise ValueError("Workflow and final native attempt identify different CAD sources")
    bodies = _records(native.get("bodies", []), "bodies")
    interferences = _records(native.get("interferences", []), "interferences")
    issues, total, severities = [], 0, set()

    def issue(code, severity, evidence, *, scope="mesh", tags=(), actions=()):
        nonlocal total
        total += 1
        severities.add(severity)
        if len(issues) >= maximum:
            return
        body_tags = list(tags)
        identity = json.dumps([source, code, scope, body_tags], separators=(",", ":"))
        issues.append({"id": hashlib.sha256(identity.encode("utf-8")).hexdigest()[:24],
                       "code": code, "severity": severity, "scope": scope,
                       "body_tags": body_tags, "cad_face_tags": [], "cell_ids": [],
                       "evidence": evidence, "candidate_actions": list(actions)})

    published = native.get("mesh_published") is True
    missing = [key for key in ("status", "reason", "units", "source_sha256") if native.get(key) is None]
    missing += [f"audit.{key}" for key in ("complete", "valid", *_AUDIT_ERRORS) if key not in audit]
    missing += [f"quality.{key}" for key in ("complete", "valid", "inverted_cells") if key not in quality]
    complete = not missing and audit.get("complete") is True and quality.get("complete") is True
    valid = (audit.get("valid") is True and quality.get("valid") is True and
             not any(_positive(audit.get(key)) for key in _AUDIT_ERRORS) and
             not _positive(quality.get("inverted_cells")))
    reason = native.get("reason")
    if published and native.get("status") not in (None, "succeeded"):
        issue("native_report_inconsistent", "error", {"status": native.get("status"), "mesh_published": True})
    if not published or reason not in (None, "none"):
        issue("mesh_refused", "error", {"reason": reason or (workflow or {}).get("reason", "unavailable"),
                                       "mesh_published": published},
              actions=("review_native_failure",))
    if published and not complete:
        issue("validation_incomplete", "error", {"audit_complete": audit.get("complete"),
                                                 "quality_complete": quality.get("complete"),
                                                 "missing_fields": missing},
              actions=("complete_native_validation",))
    if published and not valid:
        issue("native_validity_failed", "error", {"audit_valid": audit.get("valid"),
                                                  "quality_valid": quality.get("valid"),
                                                  "audit_counts": {key: audit.get(key) for key in _AUDIT_ERRORS},
                                                  "inverted_cells": quality.get("inverted_cells")},
              actions=("review_native_failure",))
    minimum, threshold = quality.get("minimum"), quality.get("warning_threshold")
    below = (isinstance(minimum, (int, float)) and not isinstance(minimum, bool) and
             _positive(threshold) and minimum < threshold)
    if native.get("quality_warning") is True or below:
        issue("preferred_quality_not_met", "warning",
              {key: quality.get(key) for key in ("metric", "minimum", "percentile_1", "percentile_5", "median", "warning_threshold")},
              actions=("review_quality_distribution", "review_sizing_without_assuming_quality_improves"))
    required = quality.get("minimum_required")
    if minimum is not None and _positive(required) and minimum < required:
        issue("strict_quality_not_met", "error", {"minimum": minimum, "minimum_required": required},
              actions=("review_explicit_quality_requirement",))
    for body in bodies:
        tag = body.get("solid_tag")
        if tag is not None:
            integer(tag, "solid_tag")
        body_reason = body.get("reason")
        if body_reason and body_reason != "none":
            issue("body_meshing_failed", "error",
                  {key: body.get(key) for key in ("reason", "surface_reason", "volume_reason", "route")},
                  scope="body", tags=() if tag is None else (tag,),
                  actions=("inspect_body_and_native_stage", "review_explicit_recipe_retry_policy"))
    for pair in interferences:
        state = pair.get("state")
        if state not in ("material_overlap", "no_material_overlap", "unresolved"):
            raise ValueError("Unknown assembly interference state")
        tags = [integer(pair.get(key), key) for key in ("first_solid_tag", "second_solid_tag")]
        if state != "no_material_overlap":
            issue("assembly_interference" if state == "material_overlap" else "assembly_interference_unresolved",
                  "warning", {"state": state, "overlap_volume_mm3": pair.get("overlap_volume")},
                  scope="assembly", tags=sorted(tags), actions=("review_preserved_part_interference",))
    nonoverlap = native.get("interbody_nonoverlap_certified") is True
    if nonoverlap and (native.get("mesh_scope") != "conforming_domain" or
                      any(pair.get("state") != "no_material_overlap" for pair in interferences)):
        issue("native_report_inconsistent", "error",
              {"mesh_scope": native.get("mesh_scope"), "interbody_nonoverlap_certified": True,
               "interference_states": sorted({pair.get("state") for pair in interferences})}, scope="assembly")
        nonoverlap = False
    workflow_failed = workflow is not None and workflow.get("status") != "succeeded"
    workflow_incomplete = workflow is not None and workflow.get("diagnostics_complete") is False
    if workflow_failed or workflow_incomplete:
        issue("workflow_incomplete" if workflow_incomplete else "workflow_failed", "error",
              {"status": workflow.get("status"), "reason": workflow.get("reason"),
               "diagnostics_complete": workflow.get("diagnostics_complete"),
               "outer_failure": workflow.get("outer_failure")}, scope="workflow",
              actions=("review_recorded_attempts_and_budgets",))
    if (workflow_failed or not published or not valid or reason not in (None, "none") or
            native.get("status") not in (None, "succeeded") or "error" in severities and complete):
        status = "failed"
    elif not complete or workflow_incomplete:
        status = "incomplete"
    else:
        status = "valid_with_warnings" if "warning" in severities else "valid"
    result = {
        "schema": DIAGNOSTIC_SCHEMA, "status": status, "source_sha256": source, "units": "mm",
        "validation": {"basis": "reported_native_checks", "native_audit_complete": audit.get("complete") is True,
                       "native_audit_valid": audit.get("valid") is True,
                       "native_quality_complete": quality.get("complete") is True,
                       "native_quality_valid": quality.get("valid") is True,
                       "mesh_scope": native.get("mesh_scope"),
                       "interbody_nonoverlap_certified": nonoverlap},
        "summary": {"nodes": native.get("nodes"), "cells": native.get("cells"),
                    "tetrahedra": native.get("tetrahedra"),
                    "quality": {key: quality.get(key) for key in _QUALITY}},
        "bodies": [{key: body.get(key) for key in _BODY} for body in bodies[:maximum]],
        "bodies_truncated": len(bodies) > maximum,
        "timings": {"native_seconds": native.get("seconds"),
                    "workflow_seconds": None if workflow is None else workflow.get("seconds"),
                    "attempts": None if workflow is None else len(workflow.get("attempts", []))},
        "issues": issues, "issue_count": total, "issues_truncated": total > maximum,
        "unavailable": ["local_cell_quality", "quality_to_CAD_face_association", "achieved_face_size",
                        "solver_error_indicator"],
    }
    return json_copy(result)


def format_diagnostics(report):
    """Render the same structured evidence for people, without adding conclusions."""
    report = json_copy(_mapping(report, "diagnostics"))
    if report.get("schema") != DIAGNOSTIC_SCHEMA:
        raise ValueError("Expected paacs-mesh-diagnostics-v1")
    summary = report["summary"]
    lines = [f"Mesh status: {report['status']}",
             f"Cells: {summary['cells']}; nodes: {summary['nodes']}; units: mm",
             f"Validation: reported native checks; scope: {report['validation']['mesh_scope']}"]
    for item in report["issues"]:
        body = f" (bodies {', '.join(map(str, item['body_tags']))})" if item["body_tags"] else ""
        lines.append(f"{item['severity']}: {item['code']}{body}")
        lines.append("  Evidence: " + json.dumps(item["evidence"], sort_keys=True, allow_nan=False))
        if item["candidate_actions"]:
            lines.append("  Consider: " + ", ".join(item["candidate_actions"]))
    if report["issues_truncated"]:
        lines.append(f"Showing {len(report['issues'])} of {report['issue_count']} issues.")
    lines.append("Local cell/face quality and solver accuracy are not assessed by this summary.")
    return "\n".join(lines)
