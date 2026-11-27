"""Bounded interpretations of native reports, without rerunning or relaxing audits."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import json

from ._cad_common import integer, json_copy, number


DIAGNOSTIC_SCHEMA = "paacs-mesh-diagnostics-v1"
CAPABILITY_SCHEMA = "paacs-mesh-capabilities-v1"
_AUDIT_ERRORS = ("missing_boundary", "unexpected_boundary", "nonmanifold_faces",
                 "inconsistent_orientation", "invalid_connectivity", "nonfinite_nodes")
_QUALITY = ("metric", "complete", "valid", "minimum", "percentile_1", "percentile_5",
            "median", "mean", "inverted_cells", "warning_threshold", "minimum_required")
_BODY = ("solid_tag", "route", "reason", "surface_reason", "volume_reason", "cad_volume",
         "mesh_volume", "relative_volume_error", "maximum_boundary_node_error",
         "maximum_boundary_sample_error", "surface_chord_error_limit",
         "boundary_embedding_certified", "centroid_outside_distance_limit", "native_centroids_checked",
         "initial_front", "surface_spacing", "source_shell", "initial_front_fallback", "volume_sizing",
         "constraint_recovery", "fidelity", "route_resolution", "initial_topology_quality",
         "material_classification", "native_candidate_queries", "native_candidate_ray_admissions")
_PARTITION = ("native_partition_seconds", "native_assembly_seconds", "partition_reason", "partition_policy",
              "partition_attempts", "partition_pieces", "partition_interfaces", "partition_work",
              "partition_volume_change", "partition_volume_admission", "sampled_partition_trim_witness",
              "partition_trim_samples", "maximum_partition_trim_drift", "maximum_partition_trim_tolerance_ratio",
              "surface_diagonal_repairs", "exact_surface_pairs", "kernel_trials", "partition_attempt_reports")
_SHAPE_REVIEW = {"metric": "subtet_mean_ratio", "near_flat_minimum": .001,
                 "low_median": .1, "low_tail_median": .2, "low_tail_percentile_5": .05}


def capabilities():
    """Describe implemented operations and restrictions; this is not a CAD guarantee."""
    from . import __version__, native_version
    from .cad import _DEFAULT_RESOURCE_POLICY, _RESOURCE_POLICIES

    return {
        "schema": CAPABILITY_SCHEMA, "package_version": __version__, "engine_version": native_version,
        "units": {"length": "mm", "volume": "mm3"},
        "operations": ["inspect", "select_faces", "preview", "body_relations", "mesh", "diagnose", "run_recipe", "publish_result",
                       "inspect_mesh", "query_facial_proximity", "visualize_mesh", "compare_meshes"],
        "cli_commands": ["capabilities", "inspect", "run", "diagnose", "inspect-mesh", "visualize", "compare"],
        "input_formats": ["STEP"], "output_formats": ["VTK", "paacsmesh"],
        "step_volume_cell_families": ["Tet4"],
        "selection": {"predicates": ["tags", "name", "group", "surface", "body", "area", "radius",
                                      "centroid_box", "normal", "adjacent_to"],
                      "ambiguity": "explicit_expectation", "identity": "source_sha256_and_guarded_local_tags",
                      "inspection_page_limit": 1000},
        "controls": {"global_size": True, "face_size": True, "box_size": True,
                     "minimum_size": True, "growth": True,
                     "max_chord_error": {"default": None, "units": "mm", "explicit": "positive"},
                     "max_relative_volume_error": {"default": None, "native_default": .01, "explicit": "0 < value < 1"},
                     "curvature_points_per_turn": {"default": None, "native_default": 0, "integer_range": [0, 64]},
                     "initial_topology_quality": {"default": None, "policy_defaults": {"normal": True, "bounded": False},
                                                  "forces_initial_front": False, "changes_native_validity_gates": False},
                     "native_source_fallback": {"default": False, "experimental": True,
                         "forces_initial_front": False, "changes_native_validity_gates": False,
                         "trigger": "remaining_named_geometric_chart_or_wire_refusal",
                         "explicit_surface_controls_bypassed": False,
                         "native_source_spacing_certified": False,
                         "source_attempt_limits": {"normal": 3, "bounded_enabled": 2, "bounded_default": 1}},
                     "boundary_first": {"default": False, "experimental": True,
                         "applies_to": "live_spatial_Tet4_constructor", "forces_initial_front": False,
                         "overrides_enabled_front_order": True, "enables_source_fallback": False,
                         "material_admission": "positive_original_trimmed_CAD_ray_else_classifier_in_serial_strict_mode",
                         "changes_native_validity_gates": False},
                     "construction_efforts": ["legacy", "auto", "thorough"],
                     "thin_region_ordering": {"default": False, "experimental": True,
                         "applies_to": "normal_native_spatial_Tet4_constructor", "forces_initial_front": False,
                         "source_samples": 64, "maximum_work": 5000000, "maximum_working_bytes": 33554432,
                         "changes_native_validity_gates": False, "changes_cad": False,
                         "qualified_sizing": ["global", "analytic_volume_controls"],
                         "face_controls": "abstain_native_face_sizing_not_qualified",
                         "scope": "advisory_phase_order; ambiguous_or_limited_queries_keep_original_order"},
                     "construction_default": "legacy",
                     "construction_stage_caps": ["max_surface_points", "max_surface_chord_points",
                                                 "max_interior_candidates", "max_interior_sites"],
                     "construction_caps_override_overall_resources": False,
                     "spatial_sizing_requires": {"solids": 1, "assembly": "conforming", "cell_family": "Tet4"}},
        "recovery": {"default": "standard", "experimental": "advanced",
                     "recipe_retries": "explicit_bounded_diagnostic_selection"},
        "validation": {"trimmed_ray_witnesses": {"default": False, "opt_in": True,
                           "eligible": "strict_serial_native_validation",
                           "ambiguous_or_unsupported": "ordinary_original_CAD_classifier",
                           "changes_geometry_tolerance": False, "certifies_solver_accuracy": False}},
        "resources": {"policies": list(_RESOURCE_POLICIES), "default": _DEFAULT_RESOURCE_POLICY,
                      "explicit_limits_override_policy": True,
                      "logical_work_is_cpu_instructions": False,
                      "modeled_memory_is_OS_hard_limit": False},
        "assemblies": {"default": "preserve_parts", "interference": "reported_without_fusing",
                       "independent_parts_certify_global_nonoverlap": False,
                       "normal_auto_dispatch": "per_original_body",
                       "bounded_auto_dispatch": "historical_whole_request",
                       "body_relations": {"explicit_pairs_required": True, "maximum_pairs": 64,
                           "relations": ["disconnected", "face_contact", "edge_contact", "point_contact",
                                         "material_overlap", "unresolved"],
                           "certifies_mesh_interfaces": False, "changes_assembly_policy": False}},
        "diagnostics": {"schema": DIAGNOSTIC_SCHEMA, "native_report_summary": True,
                        "body_cad_error": True, "local_cell_quality": True,
                        "local_cell_quality_requires": "opt_in_inspect_mesh",
                        "diagnose_includes_per_cell_quality": False,
                        "automatic_shape_review": dict(_SHAPE_REVIEW),
                        "shape_review_changes_native_acceptance": False,
                        "source_shell_attempt_limit": 3,
                        "material_classification_progress": True,
                        "material_progress_certifies_final_mesh": False,
                        "face_achieved_size": False, "solver_error_estimation": False,
                        "maximum_issues": 256, "report_budget_bytes": 8 * 1024 * 1024},
        "mesh_inspection": {"api": "inspect_mesh", "cell_id_scheme": "zero_based_native_archive_cell_index_v1",
                            "mesh_identity": "exact_archive_sha256", "page_limit": 1000,
                            "all_cells_available_by_pagination": True, "numpy_required": False,
                            "metrics": ["mean_ratio", "minimum_dihedral_deg", "scaled_jacobian"],
                            "metric_scope": "linear_Tet4_only; other_families_return_null",
                            "archive_only_provenance": "units_source_and_native_validation_unavailable",
                            "saved_reports_are_reaudited": False, "quality_field_export": "VTK"},
        "cad_proximity": {"method": "occt_trimmed_face_distance", "sample": "cell_node_centroid",
                          "requires_matching_source_kernel_and_units": True,
                          "ownership_or_quality_causation_inferred": False,
                          "maximum_face_evaluations": 100000,
                          "limit_scope": "bounded_inputs_outputs; OCCT_workspace_and_single_solve_not_OS_capped"},
        "mesh_visualization": {"format": "SVG", "deterministic_camera": True,
                               "interior_visibility": "explicit_xray_highlights",
                               "cut_rule": "whole_cells_by_centroid", "cad_context": "source_matched_display_tessellation"},
        "mesh_comparison": {"api": "compare_meshes", "requires_reported_native_validity": True,
                            "requires_same_source_and_units": True, "cell_id_correspondence_assumed": False,
                            "spatial_comparison": "common_world_grid", "solver_accuracy_inferred": False},
        "pipeline": {"module": "paacs_mesh.pipeline", "introduced": "1.2.1", "opt_in": True,
                     "compatibility_module": "paacs_mesh.experimental",
                     "functions": ["prepare_boundary", "load_boundary", "load_initial_points", "initial_points"],
                     "boundary_methods": ["surface", "edges", "cad_topology", "generate_initial_points",
                                          "mesh", "save", "to_numpy"],
                     "point_cloud_methods": ["points", "to_numpy", "save"],
                     "scope": {"single_original_solid": True, "uniform_sizing": True,
                               "recovery": "standard", "constructor": "normal_advancing_interior",
                               "volume_cell_family": "Tet4", "local_sizing_controls": False},
                     "artifacts": {"surface": "qualified_oriented_tagged_Tri3",
                                   "cad_edges": "retained_parameter_ordered_CAD_Line2_chains",
                                   "cad_incidence": "explicit_owner_lists_senses_and_wire_uses",
                                   "formats": ["VTK", "paacsmesh", "hash_bound_directory_bundle"]},
                     "initial_points": {"before_constraint_recovery_and_volume_fill": True,
                                        "includes_auxiliary_topology_finalization": True,
                                        "volume_certified": False, "immutable_owned_bytes": True,
                                        "none_and_empty_supplied_are_distinct": True,
                                        "invalid_supplied_points": "refuse_without_automatic_replacement"},
                     "identity": ["source_sha256", "boundary_sha256", "native_extension_sha256",
                                  "frozen_normalized_options", "millimetre_units"],
                     "cross_build_bundle_replay": False, "numpy_required": False,
                     "numpy_array_copies": "optional_owned_initially_readonly",
                     "changes_default_meshing": False, "changes_native_validity_gates": False,
                     "includes_AI_or_training": False},
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


def _source_shell(value):
    source = _mapping(value, "source_shell")
    if "available" in source and not isinstance(source["available"], bool):
        raise ValueError("source_shell.available must be bool")
    if source.get("accepted") is not None and not isinstance(source["accepted"], bool):
        raise ValueError("source_shell.accepted must be bool or null")
    for name in ("native_fallback_requested", "native_fallback_attempted"):
        if name in source and not isinstance(source[name], bool):
            raise ValueError(f"source_shell.{name} must be bool")
    if "native_fallback_diagnostic" in source and not isinstance(source["native_fallback_diagnostic"], str):
        raise ValueError("source_shell.native_fallback_diagnostic must be a string")
    attempts = _records(source.get("attempts", []), "source_shell.attempts")
    if len(attempts) > 3:
        raise ValueError("source_shell has more than three reported attempts")
    for attempt in attempts:
        for name in ("surface_published", "source_model_verified_unchanged",
                     "construction_budget_limited", "construction_stalled", "global_sizing_applied"):
            if name in attempt and not isinstance(attempt[name], bool):
                raise ValueError(f"source_shell.attempts.{name} must be bool")
        if attempt.get("volume_criterion_met") is not None and not isinstance(attempt["volume_criterion_met"], bool):
            raise ValueError("source_shell.attempts.volume_criterion_met must be bool or null")
        for name in ("attempt", "work", "nodes", "triangles", "curvature_points_per_turn"):
            if name in attempt:
                integer(attempt[name], "source_shell.attempts." + name, minimum=0)
        if "provider" in attempt and not isinstance(attempt["provider"], str):
            raise ValueError("source_shell.attempts.provider must be a string")
        angle = attempt.get("angular_deflection")
        if angle is not None and not _positive(angle):
            raise ValueError("source_shell.attempts.angular_deflection must be positive or null")
        if "witness" in attempt:
            witness = _mapping(attempt["witness"], "source_shell.attempts.witness")
            for name in ("stage", "detail"):
                if name in witness and not isinstance(witness[name], str):
                    raise ValueError(f"source_shell.attempts.witness.{name} must be a string")
            for name in ("edge_tag", "face_tag", "use_id"):
                if witness.get(name) is not None:
                    integer(witness[name], "source_shell.attempts.witness." + name, minimum=0)
            if "geometric_refusal" in witness and not isinstance(witness["geometric_refusal"], bool):
                raise ValueError("source_shell.attempts.witness.geometric_refusal must be bool")
    return source, attempts


def _route_resolution(value):
    receipt = _mapping(value, "route_resolution")
    for name in ("reason", "candidate"):
        if name in receipt and not isinstance(receipt[name], str):
            raise ValueError(f"route_resolution.{name} must be a string")
    for name in ("preflight_attempted", "whole_body_verified"):
        if name in receipt and not isinstance(receipt[name], bool):
            raise ValueError(f"route_resolution.{name} must be bool")
    if "preflight_work" in receipt:
        integer(receipt["preflight_work"], "route_resolution.preflight_work", minimum=0)
    attempts = _records(receipt.get("attempts", []), "route_resolution.attempts")
    if len(attempts) > 2:
        raise ValueError("route_resolution has more than two reported construction attempts")
    for attempt in attempts:
        if "published" in attempt and not isinstance(attempt["published"], bool):
            raise ValueError("route_resolution.attempts.published must be bool")
        if "work" in attempt:
            integer(attempt["work"], "route_resolution.attempts.work", minimum=0)
    return receipt


def _material_classification(value):
    receipt = _mapping(value, "material_classification")
    if "diagnostic" in receipt and not isinstance(receipt["diagnostic"], str):
        raise ValueError("material_classification.diagnostic must be a string")
    for name in ("complete", "source_model_verified_unchanged"):
        if name in receipt and not isinstance(receipt[name], bool):
            raise ValueError(f"material_classification.{name} must be bool")
    counts = ("attempts", "requested_points", "attempted_points", "native_queries",
              "strict_reuses", "ray_admissions", "work", "profiled_batches", "ray_witness_attempts")
    for name in counts:
        if name in receipt:
            integer(receipt[name], "material_classification." + name, minimum=0)
    if ("requested_points" in receipt and "attempted_points" in receipt
            and receipt["attempted_points"] > receipt["requested_points"]):
        raise ValueError("material_classification attempted points exceed requested points")
    partitions = ("native_queries", "strict_reuses", "ray_admissions")
    if all(name in receipt for name in (*partitions, "attempted_points")):
        if sum(receipt[name] for name in partitions) != receipt["attempted_points"]:
            raise ValueError("material_classification point counters disagree")
    if ("profiled_batches" in receipt and "attempts" in receipt
            and receipt["profiled_batches"] > receipt["attempts"]):
        raise ValueError("material_classification profiled batches exceed attempts")
    phases = ("classifier_setup_seconds", "native_query_seconds", "ray_witness_seconds")
    for name in ("seconds", *phases):
        value = receipt.get(name)
        if value is None and name in phases and name in receipt:
            raise ValueError(f"material_classification.{name} must be a finite nonnegative number")
        if value is not None:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"material_classification.{name} must be a finite nonnegative number or null")
            number(value, "material_classification." + name, minimum=0)
    if receipt.get("seconds") is not None and all(receipt.get(name) is not None for name in phases):
        measured = sum(receipt[name] for name in phases)
        total = receipt["seconds"]
        if measured > total + max(1e-9, abs(total)*1e-9):
            raise ValueError("material_classification profiled phases exceed elapsed seconds")
    return receipt


def _report_input(value):
    """Copy structured exceptions into the same envelope used by the CLI."""
    if isinstance(value, Exception):
        stage, reason = getattr(value, "stage", None), getattr(value, "reason", None)
        if not isinstance(stage, str) or not isinstance(reason, str):
            raise ValueError("Expected a PAACS exception with stage, reason and report")
        return {"schema": "paacs-mesh-cli-error-v1", "status": "failed", "stage": stage,
                "reason": reason, "message": str(value), "report": getattr(value, "report", None)}
    return value if isinstance(value, dict) else getattr(value, "report", None)


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
        operation = report if workflow is None and "mesh_published" not in report else None
        workflow = dict(workflow or {})
        workflow.update(status="failed", outer_failure={key: envelopes[0].get(key)
                        for key in ("stage", "reason", "message")})
        if operation is not None:
            workflow["operation_failure"] = {**workflow["outer_failure"], "details": operation}
            if operation.get("source_sha256") is not None:
                workflow["source_sha256"] = operation["source_sha256"]
            report = {}
    if "mesh_published" not in report and workflow is None:
        raise ValueError("Expected a native mesh report or paacs-mesh-workflow-v1 report")
    if "mesh_published" in report and not isinstance(report["mesh_published"], bool):
        raise ValueError("mesh_published must be bool")
    return report, workflow


def _validate_measurements(native, audit, quality):
    for name in ("requested_native_source_fallback", "requested_boundary_first", "requested_thin_region_ordering"):
        if name in native and not isinstance(native[name], bool):
            raise ValueError(f"{name} must be bool")
    if native.get("hint") is not None and not isinstance(native["hint"], str):
        raise ValueError("hint must be a string or null")
    if quality.get("metric") is not None and not isinstance(quality["metric"], str):
        raise ValueError("quality.metric must be a string or null")
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

    Accept a native ``Mesh``, its report, a workflow result/report, a PAACS
    exception, or a CLI failure envelope. Selection and request errors retain
    their stage and actual evidence; they do not imply a volume mesh was tried.
    Saved reports are interpreted as supplied, not authenticated against a mesh
    file. Missing local evidence is reported explicitly.
    """
    maximum = integer(max_issues, "max_issues")
    if maximum > 256:
        raise ValueError("max_issues must not exceed 256")
    raw = _report_input(report_or_mesh)
    raw = json_copy(_mapping(raw, "report"))
    native, workflow = _native_report(raw)
    operation_failure = (workflow or {}).get("operation_failure")
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
    if native.get("route_resolution") is not None:
        _route_resolution(native["route_resolution"])
    interferences = _records(native.get("interferences", []), "interferences")
    issues, total, severities = [], 0, set()

    def issue(code, severity, evidence, *, scope="mesh", tags=(), face_tags=(), actions=()):
        nonlocal total
        total += 1
        severities.add(severity)
        if len(issues) >= maximum:
            return
        body_tags = list(tags)
        identity = json.dumps([source, code, scope, body_tags], separators=(",", ":"))
        issues.append({"id": hashlib.sha256(identity.encode("utf-8")).hexdigest()[:24],
                       "code": code, "severity": severity, "scope": scope,
                       "body_tags": body_tags, "cad_face_tags": list(face_tags), "cell_ids": [],
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
    if operation_failure is not None:
        stage = operation_failure.get("stage")
        code = "selection_failed" if stage in ("selection", "inspection") else (
            "control_request_failed" if stage == "sizing" else "operation_failed")
        issue(code, "error", operation_failure, scope="selection" if code == "selection_failed" else "request",
              actions=("review_selection_candidates" if code == "selection_failed" else "review_request_and_report",))
    elif not published or reason not in (None, "none"):
        issue("mesh_refused", "error", {"reason": reason or (workflow or {}).get("reason", "unavailable"),
                                       "mesh_published": published, "hint": native.get("hint")},
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
    population = max((value for value in (audit.get("volume_cells"), native.get("tetrahedra"))
                      if isinstance(value, int) and not isinstance(value, bool) and value > 0), default=0)
    minimum, threshold = quality.get("minimum"), quality.get("warning_threshold")
    below = (population and isinstance(minimum, (int, float)) and not isinstance(minimum, bool) and
             _positive(threshold) and minimum < threshold)
    if population and (native.get("quality_warning") is True or below):
        issue("preferred_quality_not_met", "warning",
              {key: quality.get(key) for key in ("metric", "minimum", "percentile_1", "percentile_5", "median", "warning_threshold")},
              actions=("review_quality_distribution", "review_sizing_without_assuming_quality_improves"))
    required = quality.get("minimum_required")
    if population and minimum is not None and _positive(required) and minimum < required:
        issue("strict_quality_not_met", "error", {"minimum": minimum, "minimum_required": required},
              actions=("review_explicit_quality_requirement",))
    # These conservative review heuristics do not change preferred/strict
    # thresholds or certify solver adequacy, including intentional anisotropy.
    shape_review = {"status": "unavailable", "basis": "reported_native_summary",
                    "metric": quality.get("metric"), "solver_adequacy": "unassessed",
                    "native_validity_changed": False, "policy": dict(_SHAPE_REVIEW), "triggers": []}
    shape_review["reported_volume_population"] = population or None
    if not population:
        shape_review["unavailable_reason"] = "no_reported_volume_population"
    if population and quality.get("complete") is True and quality.get("metric") == _SHAPE_REVIEW["metric"]:
        median, p05 = quality.get("median"), quality.get("percentile_5")
        if minimum is not None or median is not None:
            shape_review["status"] = "no_automatic_trigger"
        evidence = {key: quality.get(key) for key in
                    ("metric", "minimum", "percentile_1", "percentile_5", "median", "mean",
                     "warning_threshold", "minimum_required")}
        evidence.update(basis="reported_native_summary", solver_adequacy="unassessed",
                        intentional_anisotropy_possible=True, native_validity_changed=False,
                        local_cells_identified=False)
        if minimum is not None and minimum <= _SHAPE_REVIEW["near_flat_minimum"]:
            shape_review["triggers"].append("near_flat_minimum")
            issue("near_flat_shape_review", "info", {**evidence,
                  "review_minimum_threshold": _SHAPE_REVIEW["near_flat_minimum"]},
                  actions=("inspect_low_quality_cells_and_geometry", "assess_solver_and_element_formulation"))
        if median is not None and (median <= _SHAPE_REVIEW["low_median"] or
                p05 is not None and median <= _SHAPE_REVIEW["low_tail_median"] and
                p05 <= _SHAPE_REVIEW["low_tail_percentile_5"]):
            shape_review["triggers"].append("low_shape_distribution")
            issue("low_shape_distribution_review", "info", {**evidence,
                  "review_median_threshold": _SHAPE_REVIEW["low_median"],
                  "review_tail_median_threshold": _SHAPE_REVIEW["low_tail_median"],
                  "review_percentile_5_threshold": _SHAPE_REVIEW["low_tail_percentile_5"]},
                  actions=("inspect_shape_distribution_and_regions", "assess_intended_anisotropy_and_solver_needs"))
        if shape_review["triggers"]:
            shape_review["status"] = "recommended"
    for body in bodies:
        for name in ("native_candidate_queries", "native_candidate_ray_admissions"):
            if name in body:
                integer(body[name], name, minimum=0)
        if body.get("material_classification") is not None:
            _material_classification(body["material_classification"])
        if body.get("partition_attempt_reports") is not None:
            partition_attempts = _records(body["partition_attempt_reports"], "partition_attempt_reports")
            for attempt in partition_attempts:
                for name in ("policy", "diagnostic"):
                    if name in attempt and not isinstance(attempt[name], str):
                        raise ValueError(f"partition_attempt_reports.{name} must be a string")
                for name in ("work", "pieces", "interfaces"):
                    if name in attempt:
                        integer(attempt[name], "partition_attempt_reports." + name, minimum=0)
                if "failed_piece_tag" in attempt:
                    integer(attempt["failed_piece_tag"], "partition_attempt_reports.failed_piece_tag", minimum=-1)
        tag = body.get("solid_tag")
        if tag is not None:
            integer(tag, "solid_tag")
        body_reason = body.get("reason")
        if body_reason and body_reason != "none":
            issue("body_meshing_failed", "error",
                  {key: body.get(key) for key in ("reason", "surface_reason", "volume_reason", "route")},
                  scope="body", tags=() if tag is None else (tag,),
                  actions=("inspect_body_and_native_stage", "review_explicit_recipe_retry_policy"))
        if body.get("route_resolution") is not None:
            _route_resolution(body["route_resolution"])
        source_shell = body.get("source_shell")
        if body.get("constraint_recovery") is not None:
            recovery = _mapping(body["constraint_recovery"], "constraint_recovery")
            if ("initial_topology_used" in recovery
                    and not isinstance(recovery["initial_topology_used"], bool)):
                raise ValueError("constraint_recovery.initial_topology_used must be bool")
            if ("initial_topology_diagnostic" in recovery
                    and not isinstance(recovery["initial_topology_diagnostic"], str)):
                raise ValueError("constraint_recovery.initial_topology_diagnostic must be a string")
            for name in ("initial_topology_cells", "initial_topology_import_work",
                         "initial_topology_peak_modeled_bytes"):
                if name in recovery:
                    integer(recovery[name], "constraint_recovery." + name, minimum=0)
        if source_shell is not None:
            source_shell, source_attempts = _source_shell(source_shell)
            native_sources = [row for row in source_attempts if row.get("surface_published") is True
                              and row.get("global_sizing_applied") is False]
            if (source_shell.get("native_fallback_attempted") is True
                    and source_shell.get("accepted") is True and native_sources):
                issue("source_fallback_spacing_unverified", "info", {
                    "diagnostic": source_shell.get("native_fallback_diagnostic"),
                    "attempts": native_sources[:maximum], "attempts_total": len(native_sources),
                    "attempts_truncated": len(native_sources) > maximum,
                    "sizing_completion_certified": False,
                    "native_validity_gates_changed": False}, scope="body", tags=() if tag is None else (tag,),
                    actions=("inspect_actual_route_and_achieved_sizes", "review_source_provider_and_controls"))
            if source_shell.get("available") is True and source_shell.get("accepted") is False:
                interrupted = (native.get("reason") in ("cancelled", "budget_exceeded", "allocation_failed")
                               and source_attempts and source_attempts[-1].get("volume_criterion_met") is True)
                issue("source_shell_validation_interrupted" if interrupted else "source_shell_not_accepted",
                    "info" if interrupted else "warning", {
                    "diagnostic": source_shell.get("diagnostic"),
                    "native_reason": native.get("reason"),
                    "maximum_relative_volume_error": source_shell.get("maximum_relative_volume_error"),
                    "attempts": source_attempts[:maximum], "attempts_total": len(source_attempts),
                    "attempts_truncated": len(source_attempts) > maximum,
                    "sizing_completion_certified": False}, scope="body", tags=() if tag is None else (tag,),
                    actions=("inspect_interrupted_source_validation", "review_cancellation_or_resource_reason") if interrupted
                    else ("inspect_source_geometry_attempts", "review_geometry_and_resource_limits"))
        volume_sizing = body.get("volume_sizing")
        if volume_sizing is not None:
            volume_sizing = _mapping(volume_sizing, "volume_sizing")
            if volume_sizing.get("applied") is not None and not isinstance(volume_sizing["applied"], bool):
                raise ValueError("volume_sizing.applied must be bool or null")
            if volume_sizing.get("applied") is False:
                issue("volume_sizing_not_applied", "warning", dict(volume_sizing), scope="body",
                      tags=() if tag is None else (tag,),
                      actions=("inspect_actual_route_and_achieved_sizes", "review_sized_volume_construction"))
        fallback = body.get("initial_front_fallback")
        if fallback is not None:
            fallback = _mapping(fallback, "initial_front_fallback")
            if fallback.get("volume_diagnostic") is not None and not isinstance(fallback["volume_diagnostic"], str):
                raise ValueError("initial_front_fallback.volume_diagnostic must be a string or null")
            if fallback.get("source_shell") is not None:
                _source_shell(fallback["source_shell"])
            if fallback.get("attempted") is True:
                issue("initial_construction_fallback", "info", {
                    key: fallback.get(key) for key in ("reason", "surface_diagnostic", "volume_diagnostic", "work", "seconds")},
                    scope="body", tags=() if tag is None else (tag,),
                    actions=("inspect_original_construction_failure", "review_source_shell_attempts"))
        surface = _records(body.get("surface_spacing", []), "surface_spacing")
        for row in surface:
            for name in ("budget_limited", "queue_exhausted"):
                if name in row and not isinstance(row[name], bool):
                    raise ValueError(f"surface_spacing.{name} must be bool")
        limited, limited_count, limited_faces, faces_truncated = [], 0, set(), False
        for row in surface:
            if row.get("budget_limited") is not True:
                continue
            limited_count += 1
            if len(limited) < maximum:
                limited.append(row)
            face = row.get("face_tag")
            if isinstance(face, int) and not isinstance(face, bool) and face > 0:
                if face in limited_faces or len(limited_faces) < maximum:
                    limited_faces.add(face)
                else:
                    faces_truncated = True
        initial = _mapping(body.get("initial_front", {}), "initial_front")
        if initial.get("budget_reason") is not None and not isinstance(initial["budget_reason"], str):
            raise ValueError("initial_front.budget_reason must be a string or null")
        initial_limited = (initial.get("diagnostic") == "budget_limited"
                           and initial.get("budget_reason") not in (None, "none"))
        if limited_count or initial_limited:
            face_tags = sorted(limited_faces)
            issue("construction_budget_limited", "warning", {
                "surface_stages": limited,
                "surface_stages_total": limited_count,
                "surface_stages_truncated": limited_count > len(limited),
                "cad_face_tags_truncated": faces_truncated,
                "initial_stage": {key: initial.get(key) for key in ("budget_reason", "candidate_limit", "site_limit",
                    "query_limit", "work_limit", "candidates", "nodes", "work", "queries",
                    "reserved_bulk_candidates", "reserved_bulk_sites", "source_phase_deferred")}
                    if initial_limited else None,
                "sizing_completion_certified": False,
                "native_validity_gates_changed": False}, scope="body", tags=() if tag is None else (tag,),
                face_tags=face_tags, actions=("review_construction_limits_and_requested_size",
                                            "measure_achieved_sizing_before_raising_limits"))
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
    if (workflow_failed or workflow_incomplete) and operation_failure is None:
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
    def body_summary(body):
        summary = {key: body.get(key) for key in _BODY}
        summary.update({key: body[key] for key in _PARTITION if key in body})
        if body.get("partition_attempt_reports") is not None:
            attempts = _records(body["partition_attempt_reports"], "partition_attempt_reports")
            summary["partition_attempt_reports"] = attempts[:maximum]
            summary["partition_attempt_reports_total"] = len(attempts)
            summary["partition_attempt_reports_truncated"] = len(attempts) > maximum
        def bounded_shell(shell):
            if not isinstance(shell, dict):
                return shell
            rows = _records(shell.get("attempts", []), "source_shell.attempts")
            return {**shell, "attempts": rows[:maximum], "attempts_total": len(rows),
                    "attempts_truncated": len(rows) > maximum}
        if body.get("source_shell") is not None:
            summary["source_shell"] = bounded_shell(body["source_shell"])
        if isinstance(body.get("initial_front_fallback"), dict):
            summary["initial_front_fallback"] = dict(body["initial_front_fallback"])
            if "source_shell" in summary["initial_front_fallback"]:
                summary["initial_front_fallback"]["source_shell"] = bounded_shell(summary["initial_front_fallback"]["source_shell"])
        if "surface_spacing" in body:
            rows = body["surface_spacing"]
            summary.update(surface_spacing=rows[:maximum], surface_spacing_total=len(rows),
                           surface_spacing_truncated=len(rows) > maximum)
        return summary

    result = {
        "schema": DIAGNOSTIC_SCHEMA, "status": status, "source_sha256": source, "units": "mm",
        "shape_review": shape_review,
        "validation": {"basis": "reported_native_checks", "native_audit_complete": audit.get("complete") is True,
                       "native_audit_valid": audit.get("valid") is True,
                       "native_quality_complete": quality.get("complete") is True,
                       "native_quality_valid": quality.get("valid") is True,
                       "mesh_scope": native.get("mesh_scope"),
                       "interbody_nonoverlap_certified": nonoverlap},
        "summary": {"nodes": native.get("nodes"), "cells": native.get("cells"),
                    "tetrahedra": native.get("tetrahedra"),
                    "native_hint": native.get("hint"),
                    "construction": native.get("construction"),
                    "resources": native.get("resources"), "fidelity": native.get("fidelity"),
                    "route_resolution": native.get("route_resolution"),
                    "requested_native_source_fallback": native.get("requested_native_source_fallback"),
                    "requested_boundary_first": native.get("requested_boundary_first"),
                    "requested_thin_region_ordering": native.get("requested_thin_region_ordering"),
                    "quality": {key: quality.get(key) for key in _QUALITY}},
        "bodies": [body_summary(body) for body in bodies[:maximum]],
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
    request_only = summary["cells"] is None and any(
        item["scope"] in ("selection", "request") for item in report["issues"])
    label = "Request" if request_only else "Mesh"
    lines = [f"{label} status: {report['status']}",
             f"Cells: {summary['cells']}; nodes: {summary['nodes']}; units: mm",
             f"Validation: reported native checks; scope: {report['validation']['mesh_scope']}"]
    if "shape_review" in report:
        lines.append("Shape review: " + report["shape_review"]["status"] + "; solver adequacy unassessed")
    if summary.get("native_hint"):
        lines.append("Native hint: " + summary["native_hint"])
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
