"""Run a pinned thin/curved/holdout STEP plan serially against a development build.

Planning imports no CAD or native extension. Each attempt gets its own process,
all refusals remain in the results, and inspection time is separate from meshing.
The input manifest contains paths and hashes, not downloaded geometry contents.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("initial_generation_helpers", Path(__file__).with_name("qualify_thin_curved_step.py"))
HELPERS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPERS)
SHARED = HELPERS.SHARED
LIMIT_KEYS = ("max_nodes", "max_cells", "max_work", "max_working_bytes", "max_output_bytes", "max_audit_bytes")
OPTION_KEYS = ("max_chord_error", "max_relative_volume_error", "curvature_points_per_turn",
               "initial_topology_quality", "advancing_interior_front")
MAX_ATTEMPTS = 128


def positive(value):
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise argparse.ArgumentTypeError("Expected a finite positive number")
    return result


def curvature(value):
    try:
        result = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("Expected an integer from zero to 64") from error
    if str(result) != value or not 0 <= result <= 64:
        raise argparse.ArgumentTypeError("Expected an integer from zero to 64")
    return result


def fraction(value):
    result = positive(value)
    if result >= 1:
        raise argparse.ArgumentTypeError("Expected a positive fraction less than one")
    return result


def relative_source(root, relative):
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError("Manifest source paths must stay relative to their selected root")
    result = (root / path).resolve(strict=True)
    if not result.is_relative_to(root.resolve()) or not result.is_file():
        raise ValueError("Manifest source escaped its selected root or is not a file")
    return result


def source_path(spec, manifest):
    if "fixture" in spec:
        return relative_source(ROOT, spec["fixture"])
    return relative_source(Path(manifest["corpus_root"]), spec["path"])


def material_bounds(entities):
    """Union the named lower/upper bounds from complete body inspection."""
    if not entities:
        raise ValueError("No material body bounds are available for size-factor resolution")
    boxes = [entity["bounds"] for entity in entities]
    for box in boxes:
        if (len(box["lower"]) != 3 or len(box["upper"]) != 3
            or any(not math.isfinite(value) for key in ("lower", "upper") for value in box[key])
            or any(low > high for low, high in zip(box["lower"], box["upper"]))):
            raise ValueError("Imported body bounds must be finite ordered3D boxes")
    lower = [min(box["lower"][axis] for box in boxes) for axis in range(3)]
    upper = [max(box["upper"][axis] for box in boxes) for axis in range(3)]
    diagonal = math.dist(lower, upper)
    if not math.isfinite(diagonal) or diagonal <= 0:
        raise ValueError("Imported material bounds have no finite positive diagonal")
    return lower, upper, diagonal


def select_cases(manifest, selected=None):
    cases = manifest["cases"]
    names = [case["case_id"] for case in cases]
    if len(names) != len(set(names)):
        raise ValueError("Manifest case IDs must be unique")
    if selected is None:
        # Exact baseline pairs first, then owned additions, then the fixed holdouts.
        return sorted((case for case in cases if case["kind"] != "downloaded_reference_gap"),
                      key=lambda case: (0 if case.get("baseline_action") == "reuse_exact_artifact_and_settings"
                                        else 1 if case["kind"] == "owned" else 2, names.index(case["case_id"])))
    if not selected or len(selected) != len(set(selected)) or not set(selected) <= set(names):
        raise ValueError("Selected case IDs must be unique entries in the pinned manifest")
    lookup = {case["case_id"]: case for case in cases}
    return [lookup[name] for name in selected]


def mesh_options(args, manifest):
    limits = {key: manifest["limits"][key] for key in LIMIT_KEYS}
    options = {**limits, "max_seconds": args.native_seconds, "workers": 1, "validation_workers": 1,
            "resource_policy": "normal", "construction_effort": "legacy", "recovery": "standard",
            "growth": manifest["growth"], **{key: getattr(args, key) for key in OPTION_KEYS}}
    # Keep the established default call unchanged. A frozen older package may
    # not expose these experimental options; their defaults are still disabled.
    for name in ("native_source_fallback", "boundary_first"):
        enabled = getattr(args, name, False)
        if enabled is not None and type(enabled) is not bool:
            raise ValueError(f"{name} must be bool or omitted")
        if enabled:
            options[name] = True
    return options


def build_plan(args):
    manifest = json.loads(args.case_manifest.read_text(encoding="utf-8-sig"))
    cases = select_cases(manifest, args.cases)
    if args.include_thin_forced_pair:
        base = next(case for case in manifest["cases"] if case["case_id"] == "thin-plate-010")
        variants = [{**base, "case_id": "thin-plate-010-front-quality-" + label,
                     "base_case_id": base["case_id"], "baseline_action": "explicit_forced_front_variant",
                     "variant_scope": "Controlled forced-front pair; differs from AUTO baseline route selection",
                     "option_overrides": {"advancing_interior_front": True, "initial_topology_quality": enabled}}
                    for label, enabled in (("off", False), ("on", True))]
        insertion = max((index+1 for index, case in enumerate(cases)
                         if case.get("baseline_action") == "reuse_exact_artifact_and_settings"), default=0)
        cases = cases[:insertion] + variants + cases[insertion:]
    if args.size_factor is not None and any("fixture" in case for case in cases):
        raise ValueError("A size-factor override applies only to imported corpus cases")
    if not cases or not 1 <= args.repetitions <= MAX_ATTEMPTS or len(cases) * args.repetitions > MAX_ATTEMPTS:
        raise ValueError(f"Execution plan must contain between one and {MAX_ATTEMPTS} attempts")
    rows = []
    for spec in cases:
        path = source_path(spec, manifest)
        expected = spec.get("source_sha256", spec.get("sha256"))
        if SHARED.sha256(path) != expected:
            raise ValueError(f"Pinned source changed: {spec['case_id']}")
        for repetition in range(1, args.repetitions + 1):
            rows.append({**spec, "source_file": str(path), "source_sha256": expected,
                         "repetition": repetition, "directory": f"{spec['case_id']}-repeat-{repetition:03d}",
                         **({"size_factor": args.size_factor} if args.size_factor is not None else {})})
    identity = HELPERS.development_identity(args)
    return {"schema": "paacs-initial-generation-execution-plan-v1", "input_manifest": str(args.case_manifest),
            "input_manifest_sha256": SHARED.sha256(args.case_manifest), "development_identity": identity,
            "harness_sha256": SHARED.sha256(Path(__file__)), "cases": rows,
            "options": mesh_options(args, manifest), "hard_process_seconds": args.hard_seconds,
            "requested_native_source_fallback": getattr(args, "native_source_fallback", False),
            "requested_boundary_first": getattr(args, "boundary_first", False),
            "inspection_limits": {"max_cells": manifest["limits"]["max_cells"],
                                  "max_nodes": manifest["limits"]["max_nodes"],
                                  "max_bytes": 1024**3, "max_work": 100_000_000, "max_seconds": 60.},
            "automatic_external_retries": False, "repetitions": args.repetitions,
            "timing_scope": "quiet serial single-worker development qualification; every internal retry is included in API/process time",
            "initial_points_scope": "Reported initial-front added sites and final interior nodes are separate. No initial-point coordinates are inferred from final nodes.",
            "publication_scope": "Downloaded CAD sources/figures stay local; this is not installed-wheel qualification."}


def case_options(plan, spec):
    overrides = spec.get("option_overrides", {})
    if set(overrides) - {"advancing_interior_front", "initial_topology_quality"}:
        raise ValueError("Unsupported per-case option override")
    if any(type(value) is not bool for value in overrides.values()):
        raise ValueError("Per-case toggle overrides must be booleans")
    return {**plan["options"], **overrides}


def check_boolean_intent(report, options, name):
    expected = options.get(name, False)
    key = "requested_" + name
    available = key in report
    if (expected or available) and report.get(key) is not expected:
        raise ValueError(f"Native {name.replace('_', '-')} intent differs from the pinned request")
    return available


def check_source_fallback_intent(report, options):
    return check_boolean_intent(report, options, "native_source_fallback")


def check_resources(report, options):
    resource = report.get("resources", {})
    if (resource.get("policy") != "normal" or resource.get("requested_max_work") != options["max_work"]
        or resource.get("max_working_bytes") != options["max_working_bytes"]
        or resource.get("logical_work_cap_disabled") is not (options["max_work"] is None)):
        raise ValueError("Native resource receipt differs from the pinned request")
    if options["max_work"] is not None and resource.get("effective_max_work") != options["max_work"]:
        raise ValueError("Native explicit work ceiling differs from the pinned request")
    if type(resource.get("effective_max_work")) is not int or resource["effective_max_work"] <= 0:
        raise ValueError("Native effective work counter must be a positive integer")
    check_source_fallback_intent(report, options)
    check_boolean_intent(report, options, "boundary_first")


def node_receipt(mesh, report):
    """Bounded final-node diagnostics; these do not reconstruct initial sites."""
    import numpy as np
    snapshot = mesh.to_numpy(max_bytes=256 * 1024**2)
    points = snapshot["points"]
    triangles = snapshot["cells"]["tri3"]
    boundary = np.unique(triangles["connectivity"])
    material = np.unique(snapshot["cells"]["tet4"]["connectivity"])
    interior = np.setdiff1d(material, boundary, assume_unique=True)
    def signature(array, dtype):
        return hashlib.sha256(np.asarray(array, dtype=dtype, order="C").tobytes()).hexdigest()
    return {"all_nodes": len(points), "final_material_nodes": len(material), "final_boundary_nodes": len(boundary),
            "final_interior_nodes": len(interior), "nodes_unused_by_tet4": len(points)-len(material),
            "final_boundary_ordered_points_sha256": signature(points[boundary], "<f8"),
            "final_tri3_connectivity_sha256": signature(triangles["connectivity"], "<i8"),
            "final_tri3_tags_sha256": signature(triangles["entity_tags"], "<i4"),
            "final_interior_ordered_points_sha256": signature(points[interior], "<f8"),
            "initial_front_reported_added_sites": [{"body_tag": body["solid_tag"],
                "diagnostic": body.get("initial_front", {}).get("diagnostic"),
                "reported_added_sites": body.get("initial_front", {}).get("nodes")}
                for body in report.get("bodies", [])],
            "initial_admitted_points_snapshot": None,
            "scope": "Final material-node sets and order-sensitive signatures. Initial-front counts are reported separately; final interior nodes can include later closure/recovery points."}


def worker(args):
    plan = json.loads(args.execution_plan.read_text())
    spec = next(row for row in plan["cases"] if row["directory"] == args.attempt)
    options = case_options(plan, spec)
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter(); call_started = inspection_started = None
    record = {"schema": "paacs-initial-generation-attempt-v1", "case_id": spec["case_id"],
              "repetition": spec["repetition"], "specification": spec,
              "status": "incomplete", "native_status": "not_attempted", "phase": "package_import",
              "source_sha256": spec["source_sha256"], "options": options,
              "development_identity": plan["development_identity"], "python": sys.version,
              "platform": platform.platform(), "external_attempts": 1,
              "quality_threshold_arguments": "omitted; retain Native defaults"}
    def checkpoint():
        record["worker_seconds"] = time.perf_counter()-started
        SHARED.save(output / "case.json", record)
    checkpoint()
    try:
        identity = HELPERS.development_identity(args)
        expected = plan["development_identity"]
        if any(identity[key] != expected[key] for key in ("native_extension_sha256", "source_snapshot_sha256")):
            raise ValueError("Selected source/Native identity changed after plan freeze")
        if SHARED.sha256(Path(spec["source_file"])) != spec["source_sha256"]:
            raise ValueError("Pinned STEP source changed before import")
        SHARED.load_development(args)
        import paacs_mesh as pm
        if (Path(pm.__file__).resolve() != ROOT / "python/paacs_mesh/__init__.py"
            or Path(pm._native.__file__).resolve() != args.development_native):
            raise ValueError("Loaded source package/Native extension differs from the selected candidate")
        record.update(package_path=str(Path(pm.__file__).resolve()), native_path=str(Path(pm._native.__file__).resolve()),
                      package_version=pm.__version__)
        record["phase"] = "step_import"; checkpoint(); imported = time.perf_counter()
        model = pm.read_step(spec["source_file"])
        record["import_seconds"] = time.perf_counter()-imported
        faces, bodies = list(model.face_tags), list(model.body_tags)
        if spec.get("original_solids") is not None and len(bodies) != spec["original_solids"]:
            raise ValueError("Imported body count differs from the pinned original")
        assembly = "conforming" if len(bodies) == 1 else "preserve_parts"
        size = spec.get("size")
        if size is None:
            entities = model.inspect(kind="bodies", limit=1000)["entities"]
            if len(entities) != len(bodies):
                raise ValueError("Complete material-body inspection is required for size-factor resolution")
            lower, upper, diagonal = material_bounds(entities)
            size = diagonal * spec["size_factor"]
            record.update(material_bounds={"lower": lower, "upper": upper}, imported_diagonal_mm=diagonal)
        selected = None if "face_size" not in spec else model.select_faces(surface="cylinder", radius=spec["radius"], expect=1)
        controls = [] if selected is None else [pm.FaceSize(tag, spec["face_size"]) for tag in selected.tags]
        record.update(requested_size_mm=size, expected_face_tags=faces, expected_body_tags=bodies,
                      assembly=assembly, face_sizes=[] if selected is None else
                      [{"tag": tag, "size_mm": spec["face_size"]} for tag in selected.tags])
        record["phase"] = "meshing"; record["native_status"] = "running"; checkpoint()
        call_started = time.perf_counter()
        mesh = model.mesh(size, face_sizes=controls, assembly=assembly, **options)
        record.update(call_seconds=time.perf_counter()-call_started, native_status="succeeded")
        report = mesh.report
        SHARED.retain_native_report(output, record, report, "native-report.json")
        record["phase"] = "native_checks"; checkpoint()
        record["diagnose"] = pm.diagnose(report)
        if (report.get("source_sha256") != spec["source_sha256"] or not report.get("mesh_published")
            or not report["audit"]["complete"] or not report["audit"]["valid"]
            or not report["quality"]["complete"] or not report["quality"]["valid"]):
            raise ValueError("Native identity/publication/audit/quality gate failed")
        check_resources(report, options)
        HELPERS.HELPERS.verify_tags(faces, bodies, mesh.face_tags, mesh.body_tags, "Live mesh")
        record.update(nodes=report["nodes"], tetrahedra=report["tetrahedra"],
                      mesh_scope=report.get("mesh_scope"), interferences=report.get("interferences"),
                      interbody_nonoverlap_certified=report.get("interbody_nonoverlap_certified"))
        record["phase"] = "export"; checkpoint(); exported = time.perf_counter()
        mesh.write_archive(output / "mesh.paacsmesh")
        record.update(export_seconds=time.perf_counter()-exported, mesh_sha256=SHARED.sha256(output / "mesh.paacsmesh"))
        record["phase"] = "inspection"; checkpoint(); inspection_started = time.perf_counter()
        remaining = min(plan["inspection_limits"]["max_seconds"], args.hard_seconds-(time.perf_counter()-started)-4.)
        if remaining <= 0:
            raise ValueError("No process deadline remains for independent archive inspection")
        with pm.inspect_mesh(output / "mesh.paacsmesh", **{**plan["inspection_limits"], "max_seconds": remaining}) as inspected:
            summary = inspected.summary()
        SHARED.save(output / "quality-summary.json", summary)
        if (summary["orientation_counts"]["inverted"] or summary["orientation_counts"]["degenerate"]
            or summary["unclassified_tet4_count"] or summary["mesh_sha256"] != record["mesh_sha256"]):
            raise ValueError("Independent archive orientation/ownership/identity screen failed")
        record["archive_tags"] = HELPERS.HELPERS.verify_tags(faces, bodies,
            [int(tag) for tag in summary["sizing"]["achieved"]["boundary_faces"]],
            [int(tag) for tag in summary["body_tet4_counts"]], "Decoded archive")
        record.update(quality=summary["quality"], achieved_sizing=summary["sizing"], bounds=summary["bounds"],
                      orientation_counts=summary["orientation_counts"], metric_conventions=summary["metric_conventions"])
        record["inspection_seconds"] = time.perf_counter()-inspection_started
        record["phase"] = "optional_node_diagnostics"; checkpoint(); diagnostic_started = time.perf_counter()
        try:
            record["node_receipt"] = node_receipt(mesh, report)
            if "thickness" in spec:
                record["thickness_capture"] = HELPERS.thickness_measure(mesh, spec["thickness"])
        except Exception as diagnostic_error:
            record["optional_node_diagnostics_error"] = {"type": type(diagnostic_error).__name__, "message": str(diagnostic_error)}
        record["optional_node_diagnostics_seconds"] = time.perf_counter()-diagnostic_started
        record.update(status="succeeded", phase="complete")
    except Exception as error:
        if record["phase"] == "meshing" and call_started is not None:
            record["call_seconds"] = time.perf_counter()-call_started
        if record["phase"] == "inspection" and inspection_started is not None:
            record["inspection_seconds"] = time.perf_counter()-inspection_started
        record.update(status="inspection_or_harness_failed", failed_phase=record["phase"],
                      reason=getattr(error, "reason", None), exception_type=type(error).__name__, exception=str(error))
        if isinstance(getattr(error, "report", None), dict):
            SHARED.retain_native_report(output, record, error.report, "refusal-report.json")
            record.update(status="native_refused", native_status="refused")
            try:
                record["source_fallback_intent_verified"] = check_source_fallback_intent(error.report, options)
                record["boundary_first_intent_verified"] = check_boolean_intent(error.report, options, "boundary_first")
            except ValueError as receipt_error:
                # Keep the original refusal, elapsed cost and full report even
                # when the selected build did not receive the intended option.
                record.update(status="native_receipt_mismatch", native_receipt_error=str(receipt_error))
            try:
                record["diagnose"] = pm.diagnose(error.report)
            except Exception as diagnostic_error:
                record["diagnose_error"] = str(diagnostic_error)
    record["source_unchanged"] = SHARED.sha256(Path(spec["source_file"])) == spec["source_sha256"]
    record["native_unchanged"] = SHARED.sha256(args.development_native) == plan["development_identity"]["native_extension_sha256"]
    current_identity = HELPERS.development_identity(args)
    record["source_snapshot_unchanged"] = current_identity["source_snapshot_sha256"] == plan["development_identity"]["source_snapshot_sha256"]
    record["peak_process_working_set_bytes"] = SHARED.peak_process_bytes()
    if not record["source_unchanged"] or not record["native_unchanged"] or not record["source_snapshot_unchanged"]:
        record["status"] = "identity_changed"
    checkpoint()
    print(json.dumps({key: record.get(key) for key in ("case_id", "status", "native_status", "tetrahedra", "call_seconds", "reason")}), flush=True)
    return 0 if record["status"] == "succeeded" else 1


def execute(args, plan):
    records = []
    for spec in plan["cases"]:
        identity = HELPERS.development_identity(args)
        if any(identity[key] != plan["development_identity"][key] for key in ("native_extension_sha256", "source_snapshot_sha256")):
            raise ValueError("Candidate identity changed between attempts; remaining cases were not run")
        directory = args.output / spec["directory"]
        command = [sys.executable, "-I", str(Path(__file__).resolve()), "--worker", "--attempt", spec["directory"],
                   "--output", str(directory), "--execution-plan", str(args.output / "plan.json"),
                   "--case-manifest", str(args.case_manifest), "--development-native", str(args.development_native),
                   "--occt-sdk", str(args.occt_sdk), "--native-seconds", str(args.native_seconds), "--hard-seconds", str(args.hard_seconds)]
        began = time.perf_counter(); timed_out = False
        with (args.output / (spec["directory"]+".stdout.log")).open("w") as stdout, (args.output / (spec["directory"]+".stderr.log")).open("w") as stderr:
            process = subprocess.Popen(command, cwd=ROOT, stdout=stdout, stderr=stderr, start_new_session=os.name != "nt")
            print(json.dumps({"case_id": spec["case_id"], "repetition": spec["repetition"], "phase": "started", "pid": process.pid}), flush=True)
            try:
                code = process.wait(timeout=args.hard_seconds)
            except subprocess.TimeoutExpired:
                timed_out = True; code = None
                if os.name == "nt":
                    subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"], stdout=stderr, stderr=stderr, timeout=10, check=False)
                else:
                    os.killpg(process.pid, signal.SIGKILL)
                process.wait(timeout=10)
        path = directory / "case.json"
        try:
            result = json.loads(path.read_text()) if path.is_file() else {"case_id": spec["case_id"], "status": "missing_worker_record"}
        except (OSError, ValueError) as error:
            result = {"case_id": spec["case_id"], "status": "unreadable_worker_record", "exception": str(error)}
        result.update(process_seconds=time.perf_counter()-began, process_exit_code=code, process_timeout=timed_out,
                      command=command, archive_exists=(directory / "mesh.paacsmesh").is_file())
        if timed_out:
            result["status"] = "process_timeout"
        elif code != 0 and result.get("status") == "succeeded":
            result["status"] = "process_failed_after_worker_record"
        records.append(result); SHARED.save(args.output / "results.json", records)
        print(json.dumps({key: result.get(key) for key in ("case_id", "status", "native_status", "tetrahedra", "call_seconds", "process_seconds", "reason")}), flush=True)
    return 0 if all(record["status"] == "succeeded" for record in records) else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--development-native", type=Path, required=True)
    parser.add_argument("--occt-sdk", type=Path, required=True)
    parser.add_argument("--cases", nargs="+")
    parser.add_argument("--repetitions", type=int, default=1)
    parser.add_argument("--native-seconds", type=positive, default=600.)
    parser.add_argument("--hard-seconds", type=positive, default=720.)
    parser.add_argument("--size-factor", type=positive)
    parser.add_argument("--max-chord-error", type=positive)
    parser.add_argument("--max-relative-volume-error", type=fraction)
    parser.add_argument("--curvature-points-per-turn", type=curvature)
    parser.add_argument("--initial-topology-quality", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--advancing-interior-front", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--native-source-fallback", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--boundary-first", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--include-thin-forced-pair", action="store_true")
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--execution-plan", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--attempt", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.repetitions < 1 or args.hard_seconds <= args.native_seconds:
        parser.error("Positive repetitions and a process deadline above the Native deadline are required")
    for name in ("case_manifest", "development_native", "occt_sdk"):
        setattr(args, name, getattr(args, name).resolve(strict=True))
    if not args.development_native.is_file() or not args.occt_sdk.is_dir():
        parser.error("Development extension and OCCT SDK must exist")
    args.output = args.output.resolve()
    if args.worker:
        if args.execution_plan is None or args.attempt is None:
            parser.error("Worker needs its frozen execution plan and attempt ID")
        return worker(args)
    plan = build_plan(args)
    args.output.mkdir(parents=True, exist_ok=False)
    SHARED.save(args.output / "plan.json", plan)
    if args.plan_only:
        print(json.dumps({"planned": len(plan["cases"]), "native_calls": 0, "output": str(args.output)}))
        return 0
    return execute(args, plan)


if __name__ == "__main__":
    raise SystemExit(main())
