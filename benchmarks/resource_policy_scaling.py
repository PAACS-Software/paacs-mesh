"""Compare resource policies on six owned STEP smoke/scaling cases.

The default plan contains 36 isolated attempts: two box sizes, two sphere
sizes, a through-hole and an oblique L, each under Bounded and Normal with
three repetitions. Both policies receive the same explicit work and memory
ceilings. --plan-only imports no CAD and starts no mesher.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
import importlib.util
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("resource_scaling_helpers", Path(__file__).with_name("construction_budget_trials.py"))
SHARED = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SHARED)
CASES = {
    "box-coarse": ("tests/fixtures/cad/box.step", .5),
    "box-fine": ("tests/fixtures/cad/box.step", .25),
    "sphere-coarse": ("tests/fixtures/cad/sphere.step", .5),
    "sphere-fine": ("tests/fixtures/cad/sphere.step", .25),
    "through-hole-smoke": ("tests/fixtures/cad/through_hole.step", .5),
    "oblique-smoke": ("tests/fixtures/cad/oblique_l.step", .5),
}
POLICIES = {"bounded": "legacy", "normal": "normal"}
FIXED_LIMITS = {"max_nodes": 500_000, "max_cells": 3_000_000,
                "max_output_bytes": 256 * 1024**2, "max_audit_bytes": 512 * 1024**2}
INSPECTION_LIMITS = {"max_nodes": 500_000, "max_cells": 3_000_000,
                     "max_bytes": 512 * 1024**2, "max_work": 20_000_000, "max_seconds": 60}


def plan(args):
    rows = []
    for name in args.cases:
        relative, size = CASES[name]
        source = ROOT / relative
        if not source.is_file():
            raise ValueError(f"Owned fixture is missing: {relative}")
        for policy in args.policies:
            for repetition in range(1, args.repetitions + 1):
                rows.append({"case_id": name, "fixture": relative, "source_sha256": SHARED.sha256(source),
                    "global_size_mm": size, "policy": policy, "repetition": repetition,
                    "directory": f"{name}-{policy}-repeat-{repetition:03d}"})
    return rows


def initial_front_activity(report):
    bodies = report.get("bodies", [])
    diagnostics = [(body.get("initial_front") or {}).get("diagnostic") for body in bodies]
    return None if not diagnostics or any(value is None for value in diagnostics) else any(
        value != "not_attempted" for value in diagnostics)


def verify_tags(expected_faces, expected_bodies, actual_faces, actual_bodies, label):
    receipt = {"expected_face_tags": sorted(expected_faces), "expected_body_tags": sorted(expected_bodies),
               "face_tags": sorted(actual_faces), "body_tags": sorted(actual_bodies)}
    receipt["matches"] = (receipt["expected_face_tags"] == receipt["face_tags"]
                          and receipt["expected_body_tags"] == receipt["body_tags"])
    if not receipt["matches"]:
        raise ValueError(f"{label} face/body identities differ from the imported owned fixture")
    return receipt


def verify_resources(report, policy, work, memory):
    resources = report.get("resources")
    if (not isinstance(resources, dict) or resources.get("policy") != policy
            or resources.get("requested_max_work") != work or resources.get("effective_max_work") != work
            or resources.get("max_working_bytes") != memory
            or resources.get("logical_work_cap_disabled") is not False):
        raise ValueError("Native resource receipt differs from the requested fixed comparison limits")
    return resources


def worker(args):
    SHARED.load_development(args)
    import paacs_mesh as pm

    relative, size = CASES[args.case]
    source, output = ROOT / relative, args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    options = dict(SHARED.PROFILES[POLICIES[args.policy]])
    started, meshing_started = time.perf_counter(), None
    record = {"schema": "paacs-construction-budget-case-v1", "study": "owned-resource-policy-scaling",
        "case_id": args.case, "profile": args.policy, "wall_size_mm": None, "global_size_mm": size,
        "source": relative, "source_sha256": SHARED.sha256(source), "python": sys.version,
        "host": platform.platform(), "package_version": pm.__version__,
        "native_extension_sha256": SHARED.sha256(Path(pm._native.__file__)), "options": options,
        "global_max_work": args.max_work, "global_max_working_bytes": args.max_working_bytes,
        "other_global_limits": dict(FIXED_LIMITS, max_seconds=args.max_seconds),
        "requested_controls": {"global_size_mm": size, "face_sizes": [], "volume_controls": [], "growth": 1.3},
        "workers": 1, "validation_workers": 1, "retry_count": 0, "repetition": args.repetition,
        "status": "incomplete", "native_status": "not_attempted", "phase": "import"}
    try:
        imported = time.perf_counter()
        model = pm.read_step(source)
        record["import_seconds"] = time.perf_counter() - imported
        expected_faces, expected_bodies = list(model.face_tags), list(model.body_tags)
        record["expected_tags"] = {"face_tags": expected_faces, "body_tags": expected_bodies}
        record["phase"] = "meshing"
        meshing_started = time.perf_counter()
        mesh = model.mesh(size, assembly="conforming", workers=1, validation_workers=1,
            growth=1.3, max_work=args.max_work, max_working_bytes=args.max_working_bytes,
            max_seconds=args.max_seconds, **FIXED_LIMITS, **options)
        record["call_seconds"] = time.perf_counter() - meshing_started
        record["peak_process_bytes_after_mesh"] = SHARED.peak_process_bytes()
        record["phase"] = "native_checks"
        report = mesh.report
        SHARED.retain_native_report(output, record, report, "native-report.json")
        record["native_status"] = "checks_failed"
        if (report.get("source_sha256") != record["source_sha256"] or not report["mesh_published"]
                or not report["audit"]["complete"] or not report["audit"]["valid"]
                or not report["quality"]["complete"] or not report["quality"]["valid"]):
            raise ValueError("Native source identity or unchanged publication checks failed")
        verify_resources(report, args.policy, args.max_work, args.max_working_bytes)
        record["native_status"] = "succeeded"
        record["native_seconds"] = report["seconds"]["meshing_total"]
        record["nodes"], record["tetrahedra"] = report["nodes"], report["tetrahedra"]
        record["initial_front_attempted"] = initial_front_activity(report)
        record["body_routes"] = [{"solid_tag": body.get("solid_tag"), "route": body.get("route"),
                                 "initial_front": (body.get("initial_front") or {}).get("diagnostic")}
                                for body in report["bodies"]]
        record["live_tags"] = verify_tags(expected_faces, expected_bodies, mesh.face_tags, mesh.body_tags, "Live mesh")
        record["phase"] = "export"
        exporting = time.perf_counter()
        archive = output / "mesh.paacsmesh"
        mesh.write_archive(archive)
        record["export_seconds"] = time.perf_counter() - exporting
        record["mesh_sha256"] = SHARED.sha256(archive)
        record["phase"] = "analysis"
        analyzing = time.perf_counter()
        with pm.inspect_mesh(mesh, **INSPECTION_LIMITS) as inspected:
            measured = inspected.summary()
            record["below_mean_ratio_0p1"] = inspected.query(quality_below=.1, limit=1, associations=False)["total"]
        SHARED.save(output / "quality-summary.json", measured)
        record["quality"] = measured["quality"]
        record["achieved_sizing"] = measured["sizing"]
        if measured["orientation_counts"]["inverted"] or measured["orientation_counts"]["degenerate"]:
            raise ValueError("Independent inspection found an inverted or degenerate Tet4")
        with pm.open_mesh(archive, **INSPECTION_LIMITS) as reopened:
            archived = reopened.summary()
        record["archive_tags"] = verify_tags(expected_faces, expected_bodies,
            [int(tag) for tag in archived["sizing"]["achieved"]["boundary_faces"]],
            [int(tag) for tag in archived["body_tet4_counts"]], "Reopened archive")
        if archived["unclassified_tet4_count"] or archived["mesh_sha256"] != record["mesh_sha256"]:
            raise ValueError("Reopened archive lost body ownership or archive identity")
        record["archive_tag_basis"] = "Decoded archive classifications; no repeated native CAD audit"
        record["measurement_seconds"] = time.perf_counter() - analyzing
        record["phase"], record["status"] = "complete", "succeeded"
    except Exception as error:
        if record["phase"] == "meshing" and meshing_started is not None:
            record["call_seconds"] = time.perf_counter() - meshing_started
            record["peak_process_bytes_after_mesh"] = SHARED.peak_process_bytes()
            record["native_status"] = "failed_without_native_report"
        record.update(status="failed", reason=getattr(error, "reason", None),
            stage=getattr(error, "stage", record["phase"]), exception=str(error), exception_type=type(error).__name__)
        if isinstance(getattr(error, "report", None), dict):
            SHARED.retain_native_report(output, record, error.report, "refusal-report.json")
            if record["phase"] == "meshing":
                record["native_status"] = "refused" if error.report.get("mesh_published") is False else "failed_with_report"
    record["worker_seconds"] = time.perf_counter() - started
    record["peak_process_bytes_after_analysis"] = SHARED.peak_process_bytes()
    SHARED.save(output / "case.json", record)
    print(json.dumps(record, allow_nan=False))
    return 0 if record["status"] == "succeeded" else 1


def summarize(rows, repetitions):
    groups = {}
    for row in rows:
        groups.setdefault((row["case_id"], row["profile"]), []).append(row)
    repeated, paired = [], []
    for (name, policy), attempts in groups.items():
        summary = SHARED.repeated_summary(attempts, repetitions)[0]
        repeated.append({"case_id": name, "policy": policy, **summary})
    for name in sorted({row["case_id"] for row in rows}):
        sides = {policy: groups.get((name, policy), []) for policy in POLICIES}
        all_passed = all(len(side) == repetitions and all(row["status"] == "succeeded" for row in side)
                         for side in sides.values())
        all_attempts = [attempt for attempts in sides.values() for attempt in attempts]
        identities = all(bool(all_attempts) and None not in {attempt.get(key) for attempt in all_attempts}
                         and len({attempt.get(key) for attempt in all_attempts}) == 1
                         for key in ("source_sha256", "native_extension_sha256"))
        stable = all(len({attempt.get("mesh_sha256") for attempt in attempts}) == 1
                     for attempts in sides.values())
        complete = all_passed and identities and stable
        row = {"case_id": name, "complete_valid_comparison": complete,
               "all_native_attempts_passed": all_passed,
               "quality_or_simulation_improvement_inferred": False}
        if complete:
            hashes = {side: sorted({attempt["mesh_sha256"] for attempt in attempts}) for side, attempts in sides.items()}
            row.update(meshes_identical_across_policies=hashes["bounded"] == hashes["normal"] and len(hashes["bounded"]) == 1,
                initial_front_inactive_for_both=all(attempt.get("initial_front_attempted") is False
                                                   for attempts in sides.values() for attempt in attempts),
                counts={side: {"nodes": attempts[0]["nodes"], "tetrahedra": attempts[0]["tetrahedra"]}
                        for side, attempts in sides.items()},
                quality={side: attempts[0]["quality"] for side, attempts in sides.items()},
                achieved_sizing={side: attempts[0]["achieved_sizing"] for side, attempts in sides.items()},
                interpretation="Same mesh is preservation evidence; inactive initial-front routes are not evidence of allocation improvement.")
        paired.append(row)
    return {"schema": "paacs-owned-resource-scaling-summary-v1", "attempts": len(rows),
        "successful_attempts": sum(row["status"] == "succeeded" for row in rows),
        "repeated_groups": repeated, "paired_cases": paired,
        "measurement": "One isolated serial process per attempt; native time excludes import/export/inspection; process time includes all work. No retries or validity relaxations. Repetition medians require matching source, binary and mesh identity."}


def controller(args):
    planned = plan(args)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    SHARED.save(output / "plan.json", {"schema": "paacs-owned-resource-scaling-plan-v1", "cases": planned,
        "work_cap": args.max_work, "working_bytes": args.max_working_bytes, "fixed_limits": FIXED_LIMITS,
        "native_seconds": args.max_seconds, "process_seconds": args.hard_seconds,
        "repetitions": args.repetitions, "maximum_attempts": 36,
        "case_design": "Two box/sphere scaling pairs and two broader smoke cases, not an unseen set."})
    if args.plan_only:
        print(json.dumps({"status": "planned_only", "attempts": len(planned), "native_calls": 0}))
        return 0
    rows = []
    for entry in planned:
        case = output / entry["directory"]
        command = [sys.executable, "-I", str(Path(__file__).resolve()), "--worker", "--case", entry["case_id"],
            "--policy", entry["policy"], "--repetition", str(entry["repetition"]), "--output", str(case),
            "--max-work", str(args.max_work), "--max-working-bytes", str(args.max_working_bytes),
            "--max-seconds", str(args.max_seconds), "--hard-seconds", str(args.hard_seconds)]
        if args.development_native is not None:
            command += ["--development-native", str(args.development_native.resolve()), "--occt-sdk", str(args.occt_sdk.resolve())]
        started = time.perf_counter()
        with (output / (entry["directory"] + ".stdout.log")).open("wb") as stdout, \
             (output / (entry["directory"] + ".stderr.log")).open("wb") as stderr:
            try:
                child = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=args.hard_seconds, check=False)
                code, timed_out = child.returncode, False
            except subprocess.TimeoutExpired:
                code, timed_out = None, True
        row = {"status": "timeout" if timed_out else "process_failed"}
        if (case / "case.json").is_file():
            try:
                row = json.loads((case / "case.json").read_text())
                if not isinstance(row, dict):
                    raise ValueError("Worker receipt must be an object")
            except (ValueError, UnicodeError, OSError) as error:
                row = {"status": "process_failed", "reason": "worker_receipt_unreadable", "exception": str(error)}
        if (timed_out or code != 0) and row.get("status") == "succeeded":
            row.update(worker_reported_status=row["status"], status="timeout" if timed_out else "process_failed_after_result")
        if row.get("status") == "succeeded" and row.get("source_sha256") != entry["source_sha256"]:
            row.update(worker_reported_status=row["status"], status="source_identity_changed", reason="fixed_fixture_changed")
        row.update(case_id=entry["case_id"], profile=entry["policy"], wall_size_mm=None,
            global_size_mm=entry["global_size_mm"], global_max_work=args.max_work,
            global_max_working_bytes=args.max_working_bytes, repetition=entry["repetition"],
            process_seconds=time.perf_counter() - started, process_exit_code=code, process_timeout=timed_out)
        rows.append(row)
        SHARED.save(output / "results.json", rows)
        SHARED.save(output / "summary.json", summarize(rows, args.repetitions))
        print(json.dumps({key: row.get(key) for key in ("case_id", "profile", "repetition", "status", "tetrahedra",
            "native_seconds", "below_mean_ratio_0p1", "initial_front_attempted", "reason")}), flush=True)
    return 0 if all(row["status"] == "succeeded" for row in rows) else 1


def explicit_cap(text):
    value = SHARED.resource_cap(text)
    if value is None:
        raise argparse.ArgumentTypeError("Scaling comparisons require identical explicit positive limits")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cases", choices=tuple(CASES), nargs="+", default=list(CASES))
    parser.add_argument("--policies", choices=tuple(POLICIES), nargs="+", default=list(POLICIES))
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--max-work", type=explicit_cap, default=300_000_000)
    parser.add_argument("--max-working-bytes", type=explicit_cap, default=512 * 1024**2)
    parser.add_argument("--max-seconds", type=float, default=90)
    parser.add_argument("--hard-seconds", type=float, default=150)
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--development-native", type=Path)
    parser.add_argument("--occt-sdk", type=Path)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--case", choices=tuple(CASES), default="box-coarse", help=argparse.SUPPRESS)
    parser.add_argument("--policy", choices=tuple(POLICIES), default="bounded", help=argparse.SUPPRESS)
    parser.add_argument("--repetition", type=int, default=1, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not (1 <= args.repetitions <= 3 and 1 <= args.repetition <= 3):
        parser.error("Repetition counts must be 1..3; the complete plan is bounded to 36 attempts")
    if len(set(args.cases)) != len(args.cases) or len(set(args.policies)) != len(args.policies):
        parser.error("Duplicate case or policy would overwrite evidence")
    if not (math.isfinite(args.max_seconds) and math.isfinite(args.hard_seconds)
            and 0 < args.max_seconds < args.hard_seconds):
        parser.error("Deadlines must be finite and the hard deadline must exceed the native deadline")
    if args.development_native is not None and args.occt_sdk is None:
        parser.error("Development-native loading requires --occt-sdk")
    return worker(args) if args.worker else controller(args)


if __name__ == "__main__":
    raise SystemExit(main())
