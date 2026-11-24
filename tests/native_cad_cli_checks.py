#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Verify native CAD request boundaries and existing-archive preservation."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile


def invalid_number(token):
    raise ValueError(f"nonfinite report number: {token}")


def report_from(stdout):
    reports = []
    for line in stdout.splitlines():
        try:
            value = json.loads(line, parse_constant=invalid_number)
        except ValueError:
            continue
        if isinstance(value, dict) and value.get("schema") == "paacs-native-cad-report-v1":
            reports.append(value)
    if len(reports) != 1:
        raise RuntimeError(f"Expected one native report, got {len(reports)}: {stdout[-4000:]}")
    return reports[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--fixtures", type=Path,
                        default=Path(__file__).resolve().parent / "fixtures" / "cad")
    parser.add_argument("--resource-checks-only", action="store_true",
                        help="run only resource-policy admission and preservation regressions")
    args = parser.parse_args()
    executable = args.executable.resolve()
    rows = []
    with tempfile.TemporaryDirectory(prefix="paacs-native-cad-cli-") as directory:
        output = Path(directory) / "held.paacsmesh"
        output.write_bytes(b"existing archive sentinel")

        def run(name, fixture, options=(), success=False, reason=None, size=".3"):
            previous = output.read_bytes()
            source = args.fixtures / fixture
            source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
            command = [str(executable), str(source), size, str(output),
                       "--native-cad", "--report-json", *options]
            child = subprocess.run(command, capture_output=True, text=True, timeout=60)
            receipt = report_from(child.stdout)
            if hashlib.sha256(source.read_bytes()).hexdigest() != source_hash:
                raise RuntimeError(f"{name}: modified the imported source file")
            if (child.returncode == 0) != success:
                raise RuntimeError(f"{name}: code={child.returncode}, report={receipt}, stderr={child.stderr}")
            if success:
                if receipt.get("status") != "succeeded" or receipt.get("reason") != "none":
                    raise RuntimeError(f"{name}: wrong success report: {receipt}")
                if receipt.get("mesh_published") is not True or receipt.get("archive_published") is not True:
                    raise RuntimeError(f"{name}: missing publication: {receipt}")
                if not output.read_bytes().startswith(b"PAACSM8\0"):
                    raise RuntimeError(f"{name}: archive magic missing")
                if receipt.get("nodes", 0) <= 0 or receipt.get("cells", 0) <= 0:
                    raise RuntimeError(f"{name}: empty successful mesh")
                if receipt.get("millimetres_per_model_unit") != 1:
                    raise RuntimeError(f"{name}: STEP transfer units missing")
            else:
                if receipt.get("mesh_published") is not False or receipt.get("archive_published") is not False:
                    raise RuntimeError(f"{name}: refusal published mesh/archive: {receipt}")
                if output.read_bytes() != previous:
                    raise RuntimeError(f"{name}: replaced existing archive")
                normalized = str(receipt.get("reason", "")).replace("_", "-")
                if normalized != reason:
                    raise RuntimeError(f"{name}: expected {reason}, got {normalized}: {receipt}")
            rows.append({"case": name, "passed": True, "returncode": child.returncode,
                         "archive_preserved": not success, "result": receipt})
            return receipt

        normal = run("default_normal_resource_admission", "cylinder.step", success=True, size=".7")
        resource = normal.get("resources", {})
        if (resource.get("policy") != "normal" or resource.get("requested_max_work") is not None
                or resource.get("logical_work_cap_disabled") is not True
                or resource.get("effective_max_work", 0) < 2**60
                or resource.get("max_working_bytes") != 512 * 1024**2
                or normal.get("audit_valid") is not True
                or normal.get("construction", {}).get("effort") != "legacy"):
            raise RuntimeError(f"Normal admission did not retain explicit resource evidence: {normal}")
        default_normal_bytes = output.read_bytes()
        explicit_normal = run("explicit_normal_resource_admission", "cylinder.step",
                              ("--resource-policy", "normal"), success=True, size=".7")
        if (explicit_normal.get("resources") != resource or output.read_bytes() != default_normal_bytes):
            raise RuntimeError("Default and explicit Normal requests differ in resource intent or archive")
        bounded = run("explicit_bounded_resource_admission", "cylinder.step",
                      ("--resource-policy", "bounded"), success=True, size=".7")
        bounded_resource = bounded.get("resources", {})
        if (bounded_resource.get("policy") != "bounded"
                or bounded_resource.get("requested_max_work") is not None
                or bounded_resource.get("effective_max_work") != 50_000_000
                or bounded_resource.get("logical_work_cap_disabled") is not False
                or bounded_resource.get("max_working_bytes") != 512 * 1024**2
                or bounded.get("construction", {}).get("effort") != "legacy"
                or bounded.get("audit_valid") is not True):
            raise RuntimeError("Explicit Bounded admission lost its historical limits or native checks")
        for name, options in (
                ("normal_explicit_work_before_policy", ("--max-work", "100000000", "--resource-policy", "normal")),
                ("normal_explicit_work_after_policy", ("--resource-policy", "normal", "--max-work", "100000000"))):
            explicit = run(name, "cylinder.step", (*options, "--max-working-bytes", str(1024**3)),
                           success=True, size=".7")
            resource = explicit.get("resources", {})
            if (resource.get("policy") != "normal" or resource.get("requested_max_work") != 100_000_000
                    or resource.get("effective_max_work") != 100_000_000
                    or resource.get("logical_work_cap_disabled") is not False
                    or resource.get("max_working_bytes") != 1024**3):
                raise RuntimeError(f"Explicit work did not override Normal independent of flag order: {explicit}")
        for name, options in (
                ("normal_tiny_work", ("--max-work", "1")),
                ("normal_tiny_working_memory", ("--max-working-bytes", "1024")),
                ("normal_one_byte_working_memory", ("--max-working-bytes", "1")),
                ("normal_tiny_nodes", ("--max-nodes", "5")),
                ("normal_tiny_cells", ("--max-cells", "8"))):
            refused = run(name, "cylinder.step", options,
                          reason="budget-exceeded", size=".7")
            if refused.get("resources", {}).get("policy") != "normal":
                raise RuntimeError("Resource refusal lost the requested Normal policy")
            if name == "normal_tiny_work" and (
                    refused["resources"].get("effective_max_work") != 1
                    or refused["resources"].get("logical_work_cap_disabled") is not False):
                raise RuntimeError("Normal bypassed an explicit tiny work ceiling")
        run("invalid_resource_policy", "cylinder.step", ("--resource-policy", "unexpected"), reason="invalid-request")
        run("missing_resource_policy", "cylinder.step", ("--resource-policy",), reason="invalid-request")
        # End every otherwise valid prefix with an unknown flag so these checks
        # exercise parsing and receipts without importing CAD or constructing a mesh.
        parser_cases = (
            ("parse_default", ("--unknown-resource-test-option",), "normal", None),
            ("parse_normal", ("--resource-policy", "normal", "--unknown-resource-test-option"), "normal", None),
            ("parse_normal_then_bounded", ("--resource-policy", "normal", "--resource-policy", "bounded", "--unknown-resource-test-option"), "bounded", None),
            ("parse_bounded_then_normal", ("--resource-policy", "bounded", "--resource-policy", "normal", "--unknown-resource-test-option"), "normal", None),
            ("parse_work_before_bounded", ("--max-work", "123456", "--resource-policy", "bounded", "--unknown-resource-test-option"), "bounded", 123456),
            ("parse_work_after_bounded", ("--resource-policy", "bounded", "--max-work", "123456", "--unknown-resource-test-option"), "bounded", 123456),
            ("parse_work_before_normal", ("--max-work", "123456", "--resource-policy", "normal", "--unknown-resource-test-option"), "normal", 123456),
            ("parse_work_after_normal", ("--resource-policy", "normal", "--max-work", "123456", "--unknown-resource-test-option"), "normal", 123456),
            ("parse_work_survives_policy_switch", ("--resource-policy", "normal", "--max-work", "123456", "--resource-policy", "bounded", "--unknown-resource-test-option"), "bounded", 123456),
            ("parse_last_valid_work", ("--max-work", "123456", "--resource-policy", "normal", "--max-work", "654321", "--unknown-resource-test-option"), "normal", 654321),
            ("parse_invalid_policy_preserves_normal", ("--resource-policy", "normal", "--resource-policy", "unexpected"), "normal", None),
            ("parse_missing_policy_preserves_normal", ("--resource-policy", "normal", "--resource-policy"), "normal", None),
            ("parse_zero_work_preserves_default", ("--max-work", "0"), "normal", None),
            ("parse_zero_work_preserves_normal", ("--resource-policy", "normal", "--max-work", "0"), "normal", None),
            ("parse_zero_work_preserves_explicit", ("--max-work", "123456", "--resource-policy", "normal", "--max-work", "0"), "normal", 123456),
            ("parse_invalid_work_preserves_explicit", ("--max-work", "123456", "--resource-policy", "normal", "--max-work", "invalid"), "normal", 123456),
            ("parse_negative_work_preserves_explicit", ("--max-work", "123456", "--resource-policy", "normal", "--max-work", "-1"), "normal", 123456),
            ("parse_missing_work_preserves_explicit", ("--max-work", "123456", "--resource-policy", "normal", "--max-work"), "normal", 123456),
        )
        for name, options, policy, requested in parser_cases:
            parsed = run(name, "cylinder.step", options, reason="invalid-request")
            resources = parsed.get("resources", {})
            expected_work = requested if requested is not None else (
                (2**64 - 1) // 4 if policy == "normal" else 50_000_000)
            if (resources.get("policy") != policy
                    or resources.get("requested_max_work") != requested
                    or resources.get("effective_max_work") != expected_work
                    or resources.get("logical_work_cap_disabled") is not (policy == "normal" and requested is None)
                    or resources.get("max_working_bytes") != 512 * 1024**2
                    or parsed.get("import_seconds") != 0 or parsed.get("solids") != 0):
                raise RuntimeError(f"{name}: parser error changed recognized resource intent or imported CAD: {parsed}")
        for alias in ("--max-volume-error", "--max-relative-volume-error"):
            parsed = run("parse_fidelity_" + alias.lstrip("-"), "cylinder.step",
                         ("--max-chord-error", ".0025", alias, ".005",
                          "--curvature-points-per-turn", "32", "--unknown-resource-test-option"),
                         reason="invalid-request")
            fidelity = parsed.get("fidelity", {})
            if (fidelity.get("requested_max_chord_error") != .0025
                    or fidelity.get("requested_max_relative_volume_error") != .005
                    or fidelity.get("effective_max_relative_volume_error") != .005
                    or fidelity.get("requested_curvature_points_per_turn") != 32
                    or fidelity.get("effective_curvature_points_per_turn") != 32
                    or parsed.get("import_seconds") != 0 or parsed.get("solids") != 0):
                raise RuntimeError(f"Fidelity alias lost intent or imported CAD on parser refusal: {parsed}")
        if args.resource_checks_only:
            print(json.dumps({"schema": "paacs-native-cad-cli-checks-v1", "passed": len(rows),
                              "cases": rows}, indent=2, allow_nan=False))
            return

        cylinder = run("legacy_native_cylinder", "cylinder.step",
                       ("--no-advancing-interior-front",), success=True)
        if cylinder.get("curved_faces") != 1 or cylinder["bodies"][0].get("route") != "cylinder":
            raise RuntimeError(f"cylinder did not consume native curved CAD: {cylinder}")
        capable = run("plain_default_cylinder_capability", "cylinder.step", success=True)
        capability = capable.get("route_resolution", {})
        if (capable.get("audit_valid") is not True or capability.get("whole_body_verified") is not True
                or capability.get("reason") != "verified_whole_body_primitive"
                or capable["bodies"][0].get("route") != "cylinder"
                or capability.get("preflight_work", 0) <= 0
                or capability.get("preflight_work", 0) > capable.get("validation_work", 0)):
            raise RuntimeError(f"automatic constant-field cylinder lost its verified capability: {capable}")
        plain = run("explicit_initial_cylinder_phases", "cylinder.step",
                    ("--advancing-interior-front",), success=True)
        plain_body = plain["bodies"][0]
        plain_front = plain_body.get("initial_front", {})
        if (plain.get("advancing_interior_front_policy") != "enabled"
                or plain.get("advancing_interior_front") is not True
                or plain.get("audit_valid") is not True
                or plain_front.get("wall_nodes", 0) == 0 or plain_front.get("void_nodes", 0) == 0
                or plain_front.get("source_sweep_complete") is not True
                or plain_front.get("delaunay_commits") != plain_front.get("nodes")
                or plain_front.get("full_certifications", 0) < plain_front.get("delaunay_commits", 0)
                or plain_body.get("boundary_embedding_certified") is not True
                or plain_body.get("native_centroid_source_guard") is not True
                or plain_body.get("native_centroids_checked", 0) == 0
                or plain_body.get("centroid_outside_distance_limit") != 0
                or plain_body.get("centroids_within_chord_tolerance") != 0
                or plain_body.get("interior_front", {}).get("diagnostic") != "not_attempted"):
            raise RuntimeError(f"plain defaults did not complete the source sweep and populate the interior voids: {plain}")
        # Historical Bounded requests retain their original single source-build
        # refusal and byte-identical analytic fallback. Normal may retry the
        # source chord under the same volume tolerance before placing points.
        default_sphere = run("bounded_automatic_surface_sphere", "sphere.step",
                             ("--resource-policy", "bounded"), success=True)
        default_sphere_body = default_sphere["bodies"][0]
        if (default_sphere.get("audit_valid") is not True or default_sphere.get("solids") != 1
                or default_sphere.get("advancing_interior_front_policy") != "auto"
                or default_sphere.get("advancing_interior_front_selected") is not True
                or default_sphere.get("advancing_interior_front") is not False
                or default_sphere_body.get("solid_tag") != 1
                or default_sphere_body.get("boundary_embedding_certified") is not True
                or default_sphere_body.get("native_volume_reason") != "none"
                or default_sphere_body.get("centroid_outside_distance_limit") != 0
                or default_sphere_body.get("centroids_within_chord_tolerance") != 0
                or default_sphere_body.get("maximum_outside_centroid_distance") != 0
                or default_sphere_body.get("relative_volume_error", 1) > .01):
            raise RuntimeError(f"automatic defaults regressed the established sphere's strict native publication: {default_sphere}")
        fallback = default_sphere.get("initial_front_fallback", {})
        body_fallback = default_sphere_body.get("initial_front_fallback", {})
        rejected_source = fallback.get("source_shell", {})
        rejected_attempts = rejected_source.get("attempts", [])
        if (fallback.get("attempted") is not True or fallback.get("reason") != "geometry_error"
                or fallback.get("work", 0) <= 0 or fallback.get("seconds", -1) < 0
                or fallback.get("surface_diagnostic") in {None, "not_attempted"}
                or rejected_source.get("available") is not True or rejected_source.get("accepted") is not False
                or len(rejected_attempts) != 1
                or default_sphere_body.get("volume_sizing", {}).get("applied") is not False
                or default_sphere_body.get("volume_sizing", {}).get("sizing_completion_certified") is not False
                or default_sphere.get("validation_work", 0) <= fallback.get("work", 0)
                or default_sphere.get("total_seconds", 0) < fallback.get("seconds", 0)
                or body_fallback != fallback
                or default_sphere_body.get("initial_front", {}).get("diagnostic") != "not_attempted"):
            raise RuntimeError(f"automatic sphere fallback hid its failed pre-point surface work: {default_sphere}")
        rejected = rejected_attempts[0]
        if rejected.get("surface_published") is True:
            if (rejected.get("volume_criterion_met") is not False
                    or rejected.get("source_model_verified_unchanged") is not True
                    or fallback.get("volume_diagnostic") != "native-plc-source-volume-error"):
                raise RuntimeError(f"published sphere source lost its volume-refusal certificate: {fallback}")
        elif (rejected.get("surface_published") is not False
                or rejected.get("volume_criterion_met") is not None
                or rejected.get("surface_diagnostic") != "geometry_error"
                or fallback.get("volume_diagnostic") != rejected.get("surface_diagnostic")):
            raise RuntimeError(f"unpublished sphere source invented a volume check or lost its original refusal: {fallback}")
        default_sphere_bytes = output.read_bytes()
        forced_sphere = run("forced_initial_sphere_surface_refusal", "sphere.step",
                            ("--resource-policy", "bounded", "--advancing-interior-front"), reason="geometry-error")
        if forced_sphere.get("initial_front_fallback", {}).get("attempted") is not False:
            raise RuntimeError("explicit initial-front ON silently fell back after source refusal")
        sphere_reference = run("automatic_sphere_matches_legacy_off", "sphere.step",
                               ("--resource-policy", "bounded", "--no-advancing-interior-front"), success=True)
        if output.read_bytes() != default_sphere_bytes:
            raise RuntimeError("automatic pre-point fallback changed the established sphere archive")
        reference_body = sphere_reference["bodies"][0]
        for key in ("route", "solid_tag", "boundary_embedding_certified", "native_surface_reason", "native_volume_reason",
                    "native_centroid_threads", "native_centroids_checked", "native_centroid_exact_reuses",
                    "native_centroid_reference_queries", "native_centroid_geometry_pairs", "native_centroid_source_guard",
                    "centroid_outside_distance_limit", "centroids_within_chord_tolerance", "native_centroid_distance_queries",
                    "maximum_outside_centroid_distance", "native_chord_samples_checked", "chord_samples_within_tolerance",
                    "maximum_outside_chord_sample_distance", "cad_volume", "mesh_volume", "relative_volume_error",
                    "boundary_node_error", "boundary_sample_error", "surface_chord_error_limit"):
            if key not in default_sphere_body or default_sphere_body[key] != reference_body.get(key):
                raise RuntimeError(f"automatic sphere fallback changed the established certificate field {key}")
        sphere_sizes = []
        for size in (".5", ".25"):
            current = run(f"normal_sphere_source_retry_{size}", "sphere.step",
                          ("--max-working-bytes", str(1024**3)), success=True, size=size)
            body = current["bodies"][0]
            shell = body.get("source_shell", {})
            attempts = shell.get("attempts", [])
            front = body.get("initial_front", {})
            sizing = body.get("volume_sizing", {})
            # This single-solid route publishes Tri3 boundary + Tet4 volume.
            # Raw CLI reports total cells; the Python wrapper adds tetrahedra.
            tetrahedra = current.get("cells", 0) - body.get("surface_triangles", 0)
            if (current.get("resources", {}).get("policy") != "normal"
                    or current.get("initial_front_fallback", {}).get("attempted") is not False
                    or current.get("advancing_interior_front") is not True
                    or current.get("audit_valid") is not True
                    or current.get("solids") != 1 or current.get("auxiliary_faces") != 0
                    or tetrahedra <= 0
                    or body.get("route") != "native_spatial_sizing"
                    or body.get("boundary_embedding_certified") is not True
                    or body.get("native_volume_reason") != "none"
                    or body.get("native_centroids_checked", 0) != tetrahedra
                    or body.get("centroid_outside_distance_limit") != 0
                    or body.get("maximum_outside_centroid_distance") != 0
                    or body.get("relative_volume_error", 1) > .01
                    or shell.get("accepted") is not True or shell.get("diagnostic") != "accepted"
                    or shell.get("requested_global_size") != float(size)
                    or shell.get("requested_chord_error") is not None
                    or shell.get("maximum_relative_volume_error") != .01
                    or not 2 <= len(attempts) <= 3
                    or attempts[0].get("volume_criterion_met") is not False
                    or attempts[0].get("relative_volume_error", 0) <= .01
                    or attempts[-1].get("volume_criterion_met") is not True
                    or attempts[-1].get("relative_volume_error", 1) > .01
                    or front.get("work", 0) <= 0 or front.get("nodes", 0) <= 0
                    or sizing.get("applied") is not True
                    or sizing.get("requested_global_size") != float(size)
                    or sizing.get("sizing_completion_certified") is not False):
                raise RuntimeError(f"Normal source retry lost geometry, work or sizing evidence: {current}")
            if (any(attempt.get("work", 0) <= 0 or attempt.get("seconds", -1) < 0
                    or attempt.get("source_model_verified_unchanged") is not True for attempt in attempts)
                    or any(after["chord_error_limit"] >= before["chord_error_limit"]
                           for before, after in zip(attempts, attempts[1:]))
                    or sum(attempt["work"] for attempt in attempts) > body.get("native_surface_work", 0)
                    or body.get("native_surface_work", 0) > current.get("validation_work", 0)):
                raise RuntimeError(f"Source retries were not progressively tighter and fully charged: {current}")
            sphere_sizes.append((tetrahedra, front["nodes"], output.read_bytes()))
        # This fixed owned sphere used to publish the same fan at both sizes.
        # Count and initial-site responses demonstrate use of the volume field;
        # they do not certify every achieved edge length or solver accuracy.
        if (sphere_sizes[1][0] <= sphere_sizes[0][0] or sphere_sizes[1][1] <= sphere_sizes[0][1]
                or sphere_sizes[1][2] == sphere_sizes[0][2]):
            raise RuntimeError("Normal sphere construction still ignores the requested global size")
        sized = run("curved_face_and_box", "cylinder.step",
                    ("--no-advancing-interior-front", "--face-size", "1", ".38", "--size-box", "-.3", "-.3", ".7", ".3", ".3", "1.3", ".31",
                     "--min-size", ".34", "--growth", "1.55", "--max-work", "100000000"),
                    success=True, size=".7")
        sized_body = sized["bodies"][0]
        if (sized.get("native_spatial_sizing") is not True or sized.get("audit_valid") is not True
                or sized.get("face_size_controls") != 1 or sized.get("size_box_controls") != 1
                or sized_body.get("route") != "native_spatial_sizing"
                or sized_body.get("local_delaunay_refinement_nodes", 0) == 0):
            raise RuntimeError(f"curved controls did not refine and audit the native volume: {sized}")
        advancing = run("advancing_interior_front", "cylinder.step",
                        ("--face-size", "1", ".16", "--size-box", "-.3", "-.3", ".7", ".3", ".3", "1.3", ".15",
                         "--min-size", ".12", "--growth", "1.55", "--advancing-interior-front",
                         "--improve-tet-quality", "--max-work", "100000000"),
                        success=True, size=".7")
        front = advancing["bodies"][0].get("interior_front", {})
        placement = advancing["bodies"][0].get("initial_front", {})
        if (advancing.get("advancing_interior_front") is not True or advancing.get("audit_valid") is not True
                or placement.get("nodes", 0) == 0 or placement.get("source_activations", 0) == 0
                or placement.get("wall_nodes", 0) == 0
                or front.get("nodes", 0) == 0
                or front.get("front_advances", 0) + front.get("void_activations", 0) == 0
                or front.get("diagnostic") not in {"targets_satisfied", "targets_unmet", "budget_limited"}):
            raise RuntimeError(f"advancing interior front did not report audited construction: {advancing}")
        uniform_options = ("--resource-policy", "bounded", "--max-work", "100000000",
                           "--max-chord-error", ".008660254037844387",
                           "--min-quality", "0", "--warning-quality", ".25")
        initial = run("bounded_uniform_initial_interior_front", "cylinder.step",
                      uniform_options, success=True)
        initial_body = initial["bodies"][0]
        seed = initial_body.get("initial_front", {})
        if (initial.get("advancing_interior_front") is not True
                or initial.get("advancing_interior_front_policy") != "auto"
                or initial.get("advancing_interior_front_selected") is not True
                or initial.get("face_size_controls") != 0 or initial.get("size_box_controls") != 0
                or initial.get("improve_tet_quality") is not False
                or initial.get("audit_valid") is not True
                or seed.get("nodes", 0) == 0 or seed.get("propagated_activations", 0) == 0
                or seed.get("source_sweep_complete") is not True
                or seed.get("wall_nodes", 0) == 0 or seed.get("void_nodes", 0) == 0
                or seed.get("wall_nodes", 0) + seed.get("void_nodes", 0) > seed.get("nodes", 0)
                or seed.get("void_candidates", 0) < seed.get("void_nodes", 0)
                or not isinstance(seed.get("void_queue_exhausted"), bool)
                or (seed.get("maximum_unresolved_void_ratio", -1) is None
                    and seed.get("void_queue_exhausted") is not False)
                or (seed.get("maximum_unresolved_void_ratio", -1) is not None
                    and seed.get("maximum_unresolved_void_ratio", -1) < 0)
                or seed.get("delaunay_commits") != seed.get("nodes")
                or seed.get("full_certifications", 0) < seed.get("delaunay_commits", 0)
                or seed.get("maximum_planned_cavity", 0) <= 64
                or seed.get("compactions", 0) == 0 or seed.get("retired_cells_reclaimed", 0) == 0
                or seed.get("storage_rejections") != 0
                or not 0 < seed.get("compaction_peak_modeled_bytes", 0) <= seed.get("auxiliary_working_bytes", 0)
                or initial.get("quality_percentile_1", 0) < .25
                or initial_body.get("interior_front", {}).get("diagnostic") != "not_attempted"):
            raise RuntimeError(f"uniform initial construction did not run independently: {initial}")
        if (initial_body.get("boundary_embedding_certified") is not True
                or initial_body.get("native_centroid_source_guard") is not True
                or initial_body.get("native_centroids_checked", 0) <= 0
                or initial_body.get("centroid_outside_distance_limit") != 0
                or initial_body.get("centroids_within_chord_tolerance") != 0
                or initial_body.get("maximum_outside_centroid_distance") != 0
                or initial_body.get("recovery_steiner_nodes") != 0
                or initial_body.get("local_delaunay_refinement_nodes") != 0
                or initial_body.get("local_delaunay_quality_edits") != 0):
            raise RuntimeError(f"bounded initial construction lost its strict native publication guards: {initial}")
        automatic_bytes = output.read_bytes()
        forced = run("initial_front_last_flag_on", "cylinder.step",
                     (*uniform_options, "--no-advancing-interior-front", "--advancing-interior-front"), success=True)
        if (forced.get("advancing_interior_front") is not True
                or forced.get("advancing_interior_front_policy") != "enabled"
                or output.read_bytes() != automatic_bytes):
            raise RuntimeError("automatic and last-explicit-on initial construction differ")
        minimum = run("automatic_uniform_minimum_size", "cylinder.step",
                      (*uniform_options, "--min-size", ".1"), success=True)
        if minimum.get("advancing_interior_front") is not True or output.read_bytes() != automatic_bytes:
            raise RuntimeError("a uniform minimum-size bound changed or disabled automatic initial construction")
        single_part = run("automatic_single_independent_part", "cylinder.step",
                          (*uniform_options, "--independent-parts"), success=True)
        single_body = single_part["bodies"][0]
        if (single_part.get("solids") != 1 or single_part.get("mesh_scope") != "independent_parts"
                or single_part.get("advancing_interior_front") is not True
                or single_part.get("advancing_interior_front_policy") != "auto"
                or single_part.get("advancing_interior_front_selected") is not True
                or single_body.get("solid_tag") != 1
                or single_body.get("initial_front", {}).get("nodes", 0) == 0
                or single_body.get("initial_front", {}).get("source_sweep_complete") is not True
                or single_body.get("boundary_embedding_certified") is not True
                or single_body.get("native_centroid_source_guard") is not True
                or single_body.get("native_centroids_checked", 0) == 0
                or single_body.get("centroid_outside_distance_limit") != 0
                or single_body.get("centroids_within_chord_tolerance") != 0):
            raise RuntimeError("single-part automatic construction lost the independent wrapper's native tags/audits")
        run("forced_single_independent_front_refused", "cylinder.step",
            ("--advancing-interior-front", "--independent-parts"), reason="unsupported-control")
        legacy = run("initial_front_last_flag_off", "cylinder.step",
                     (*uniform_options, "--advancing-interior-front", "--no-advancing-interior-front"), success=True)
        if (legacy.get("advancing_interior_front") is not False
                or legacy.get("advancing_interior_front_policy") != "disabled"
                or legacy.get("advancing_interior_front_selected") is not False
                or legacy["bodies"][0].get("route") != "cylinder"
                or legacy["bodies"][0].get("initial_front", {}).get("diagnostic") != "not_attempted"):
            raise RuntimeError("last-explicit-off did not restore the legacy cylinder route")
        no_fallback = run("automatic_without_native_fallback", "cylinder.step",
                          ("--no-native-fallback",), success=True)
        if (no_fallback.get("advancing_interior_front") is not False
                or no_fallback["bodies"][0].get("route") != "cylinder"):
            raise RuntimeError("automatic policy displaced a recognized route without native fallback")
        run("forced_front_without_native_fallback_refused", "cylinder.step",
            ("--no-native-fallback", "--advancing-interior-front"), reason="unsupported-control")
        native = run("automatic_native_linear_elements", "box.step", ("--elements", "native"), success=True)
        native_bytes = output.read_bytes()
        native_off = run("legacy_native_linear_elements", "box.step",
                         ("--elements", "native", "--no-advancing-interior-front"), success=True)
        if (native.get("advancing_interior_front") is not False
                or native.get("advancing_interior_front_selected") is not False
                or native.get("advancing_interior_front_policy") != "auto"
                or native_off.get("advancing_interior_front") is not False
                or native_off.get("advancing_interior_front_policy") != "disabled"
                or any(body.get("initial_front", {}).get("diagnostic") != "not_attempted" for body in native["bodies"])
                or output.read_bytes() != native_bytes):
            raise RuntimeError("automatic policy changed the existing non-Tet4 route")
        default_box = run("default_uniform_box", "box.step", success=True)
        box_resolution = default_box.get("route_resolution", {})
        if (default_box.get("advancing_interior_front") is not False
                or default_box.get("advancing_interior_front_policy") != "auto"
                or default_box.get("audit_valid") is not True
                or box_resolution.get("reason") != "verified_whole_body_mapped_or_sweep"
                or box_resolution.get("candidate") != "mapped_hex8"
                or box_resolution.get("whole_body_verified") is not True
                or box_resolution.get("preflight_work", 0) <= 0
                or default_box["bodies"][0].get("initial_front", {}).get("diagnostic") != "not_attempted"
                or default_box["bodies"][0].get("solid_tag") != 1
                or default_box["bodies"][0].get("relative_volume_error", 1) > .01):
            raise RuntimeError(f"automatic box lost verified whole-body selection or Native validity: {default_box}")
        default_box_bytes = output.read_bytes()
        run("advancing_front_independent_parts_refused", "assembly.step",
            ("--advancing-interior-front", "--independent-parts"), reason="unsupported-control")
        standard = run("legacy_planar_box", "box.step", ("--no-advancing-interior-front",), success=True, size="3e-1")
        standard_bytes = output.read_bytes()
        if standard_bytes != default_box_bytes:
            raise RuntimeError("verified automatic box changed the established mapped construction")
        advanced = run("advanced_recognized_box", "box.step", ("--no-advancing-interior-front", "--recovery-mode", "advanced"), success=True)
        quality_disabled = run("advanced_box_with_regular_quality_disabled", "box.step",
                               ("--no-advancing-interior-front", "--recovery-mode", "advanced", "--no-improve-regular-quality"), success=True)
        quality_enabled = run("standard_box_with_regular_quality_enabled", "box.step",
                              ("--no-advancing-interior-front", "--improve-regular-quality"), success=True)
        if standard.get("improve_regular_quality") is not False or \
                advanced.get("improve_regular_quality") is not True or \
                quality_disabled.get("improve_regular_quality") is not False or \
                quality_enabled.get("improve_regular_quality") is not True or \
                output.read_bytes() != standard_bytes:
            raise RuntimeError("regular quality policies did not preserve the recognized box route")
        reset = run("explicit_standard_box", "box.step",
                    ("--no-advancing-interior-front", "--recovery-mode", "advanced", "--recovery-mode", "standard"), success=True)
        if standard["recovery_mode"] != "standard" or advanced["recovery_mode"] != "advanced" or \
                reset["recovery_mode"] != "standard" or output.read_bytes() != standard_bytes or \
                advanced["bodies"][0]["constraint_recovery"]["diagnostic"] != "not_attempted":
            raise RuntimeError("recovery mode did not preserve the standard or recognized route")
        local_delaunay = run("advanced_box_with_local_delaunay", "box.step",
                             ("--no-advancing-interior-front", "--recovery-mode", "advanced", "--live-material-recovery",
                              "--local-delaunay-refinement"), success=True)
        body = local_delaunay["bodies"][0]
        if (standard.get("local_delaunay_refinement") is not False
                or advanced.get("local_delaunay_refinement") is not False
                or local_delaunay.get("local_delaunay_refinement") is not True
                or body.get("local_delaunay_support_nodes") != 0
                or body.get("local_delaunay_refinement_nodes") != 0
                or body.get("local_delaunay_quality_edits") != 0
                or body.get("local_delaunay_used_global_fallback") is not False
                or output.read_bytes() != standard_bytes):
            raise RuntimeError("local Delaunay option did not preserve the recognized box route")
        assembly = run("separated_assembly", "assembly.step", success=True)
        if assembly.get("solids") != 2 or len(assembly.get("bodies", [])) != 2:
            raise RuntimeError(f"assembly bodies were lost: {assembly}")
        if (assembly.get("advancing_interior_front") is not False
                or assembly.get("advancing_interior_front_selected") is not False
                or assembly.get("advancing_interior_front_policy") != "auto"
                or any(body.get("initial_front", {}).get("diagnostic") != "not_attempted" for body in assembly["bodies"])):
            raise RuntimeError("automatic policy leaked into assembly children")
        assembly_bytes = output.read_bytes()
        run("legacy_separated_assembly", "assembly.step", ("--no-advancing-interior-front",), success=True)
        if output.read_bytes() != assembly_bytes:
            raise RuntimeError("automatic policy changed an existing conforming assembly archive")
        independent = run("automatic_independent_parts", "assembly.step", ("--independent-parts",), success=True)
        if (independent.get("advancing_interior_front") is not False or independent.get("solids") != 2
                or any(body.get("initial_front", {}).get("diagnostic") != "not_attempted" for body in independent["bodies"])):
            raise RuntimeError("automatic policy leaked into independent assembly children")
        sphere = run("legacy_native_surface_sphere", "sphere.step",
                     ("--no-advancing-interior-front", "--warning-quality", ".25"), success=True)
        if sphere["bodies"][0].get("route") != "native_surface_star" or sphere["bodies"][0].get("native_volume_reason") != "none":
            raise RuntimeError(f"sphere fallback certificate missing: {sphere}")
        refusals = [
            ("disabled_sphere_fallback", "sphere.step", ("--resource-policy", "bounded", "--no-native-fallback"), "unsupported-body"),
            ("malformed_brep", "malformed.brep", ("--brep",), "parse-failed"),
            ("unknown_face", "cylinder.step", ("--face-size", "999999", ".1"), "invalid-request"),
            ("unknown_spatial_face", "cylinder.step",
             ("--native-sizing", "--face-size", "999999", ".1"), "invalid-request"),
            ("inverted_size_box", "box.step",
             ("--size-box", ".4", ".2", ".2", ".3", ".8", ".8", ".1"), "invalid-request"),
            ("nonfinite_size_box", "box.step",
             ("--size-box", ".2", ".2", ".2", ".8", ".8", ".8", "nan"), "invalid-request"),
            ("spatial_sizing_assembly", "assembly.step",
             ("--size-box", ".2", ".2", ".2", ".8", ".8", ".8", ".1"), "unsupported-control"),
            ("strict_quality", "cylinder.step", ("--min-quality", "1"), "quality-requirement-not-met"),
            ("zero_work", "cylinder.step", ("--max-work", "0"), "invalid-request"),
            ("one_work", "cylinder.step", ("--max-work", "1"), "budget-exceeded"),
            ("working_storage", "cylinder.step", ("--max-working-bytes", "1"), "budget-exceeded"),
            ("audit_storage", "cylinder.step", ("--max-audit-bytes", "1"), "budget-exceeded"),
            ("node_budget", "cylinder.step", ("--max-nodes", "1"), "budget-exceeded"),
            ("cell_budget", "cylinder.step", ("--max-cells", "1"), "budget-exceeded"),
            ("output_budget", "cylinder.step", ("--max-bytes", "1"), "budget-exceeded"),
            ("fallback_cell_budget", "sphere.step", ("--no-advancing-interior-front", "--max-cells", "1"), "budget-exceeded"),
            ("fallback_node_budget", "sphere.step", ("--no-advancing-interior-front", "--max-nodes", "1"), "budget-exceeded"),
            ("fallback_working_storage", "sphere.step", ("--no-advancing-interior-front", "--max-working-bytes", "1"), "budget-exceeded"),
            ("fallback_strict_quality", "sphere.step", ("--no-advancing-interior-front", "--min-quality", "1"), "quality-requirement-not-met"),
            ("file_budget", "box.step", ("--max-file-bytes", "1"), "file-budget"),
            ("invalid_budget_number", "box.step", ("--max-nodes", "-1"), "invalid-request"),
            ("invalid_element_flag", "box.step", ("--elements", "tet10"), "invalid-request"),
            ("missing_element_flag", "box.step", ("--elements",), "invalid-request"),
            ("invalid_recovery_mode", "box.step", ("--recovery-mode", "fast"), "invalid-request"),
            ("missing_recovery_mode", "box.step", ("--recovery-mode",), "invalid-request"),
            ("nan_growth", "cylinder.step", ("--growth", "nan"), "invalid-request"),
            ("underflow_box_coordinate", "box.step",
             ("--size-box", "1e-999", ".2", ".2", ".8", ".8", ".8", ".1"), "invalid-request"),
            ("legacy_minimum_without_control", "cylinder.step", ("--no-advancing-interior-front", "--min-size", ".1"), "unsupported-control"),
            ("independent_minimum_without_control", "cylinder.step", ("--independent-parts", "--min-size", ".1"), "unsupported-control"),
        ]
        for name, fixture, options, reason in refusals:
            run(name, fixture, options, reason=reason)
        for name, size in (("nan_size", "nan"), ("trailing_size", ".3junk"),
                           ("leading_space", " .3"), ("trailing_space", ".3 "),
                           ("leading_plus", "+.3"), ("hex_size", "0x1p-2"),
                           ("decimal_comma", "0,3"), ("overflow_size", "1e999"),
                           ("zero_size", "0")):
            run(name, "box.step", reason="invalid-request", size=size)
        # Use a private copy: the input itself must survive output alias errors.
        source = Path(directory) / 'source.step'
        source_bytes = (args.fixtures / 'cylinder.step').read_bytes()
        source.write_bytes(source_bytes)
        alias = Path(directory) / 'source-alias.step'
        os.link(source, alias)
        for name, destination in [('same_input_output', source), ('hardlink_input_output', alias)]:
            child = subprocess.run([str(executable), str(source), '.3', str(destination),
                                    '--native-cad', '--report-json'],
                                   capture_output=True, text=True, timeout=60)
            receipt = report_from(child.stdout)
            if child.returncode != 2 or receipt.get('reason') != 'invalid_request' or \
                    receipt.get('mesh_published') or receipt.get('archive_published') or \
                    source.read_bytes() != source_bytes or destination.read_bytes() != source_bytes:
                raise RuntimeError(f'{name}: source alias was not preserved: {receipt}')
            rows.append({'case': name, 'passed': True, 'returncode': child.returncode,
                         'source_preserved': True, 'result': receipt})
    print(json.dumps({"schema": "paacs-native-cad-cli-checks-v1", "passed": len(rows),
                      "cases": rows}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
