"""Serial baseline on eight owned thin/curved STEP requests.

This is a measured qualification probe, not a simulation-adequacy certificate.
Every Native refusal, process deadline and incomplete inspection is retained.
No geometry-specific plate helper or explicit primitive mesher is called.
The default uses the frozen installed wheel. Explicit development loading has
separate binary/source receipts and never presents itself as a wheel test.
Thin-part thickness measurements require NumPy in the worker interpreter.
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
SPEC = importlib.util.spec_from_file_location("thin_curved_helpers", Path(__file__).with_name("resource_policy_scaling.py"))
HELPERS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPERS)
SHARED = HELPERS.SHARED
CASES = {
    "thin-plate-025": {"fixture": "tests/fixtures/cad/thin_plate.step", "size": .25, "thickness": .002},
    "thin-plate-010": {"fixture": "tests/fixtures/cad/thin_plate.step", "size": .10, "thickness": .002},
    "plate-hole-4-bore-1": {"fixture": "examples/elasticity_plate_hole/plate-hole.step", "size": 4., "radius": 5., "face_size": 1., "thickness": 2.},
    "plate-hole-2-bore-05": {"fixture": "examples/elasticity_plate_hole/plate-hole.step", "size": 2., "radius": 5., "face_size": .5, "thickness": 2.},
    "sphere-05": {"fixture": "tests/fixtures/cad/sphere.step", "size": .5},
    "sphere-025": {"fixture": "tests/fixtures/cad/sphere.step", "size": .25},
    "cylinder-wall-015": {"fixture": "tests/fixtures/cad/cylinder.step", "size": .5, "radius": 1., "face_size": .15},
    "torus-rotated-05": {"fixture": "tests/fixtures/curved-reference/torus_rotated.step", "size": .5},
}
DEFAULT_CASES = tuple(CASES)
CASES.update({
    "thin-plate-020": {"fixture": "tests/fixtures/cad/thin_plate.step", "size": .20, "thickness": .002},
    "thin-plate-015": {"fixture": "tests/fixtures/cad/thin_plate.step", "size": .15, "thickness": .002},
})
LIMITS = {"max_nodes": 500_000, "max_cells": 3_000_000, "max_output_bytes": 256 * 1024**2,
          "max_audit_bytes": 512 * 1024**2, "max_working_bytes": 1024**3, "max_work": None,
          "max_seconds": 90., "workers": 1, "validation_workers": 1}
PROVENANCE = {
    "thin_plate": "tests/occt_corpus.cpp: OCCT box2x1x.002; planar/sweep corpus previously checked analytic volume, classifications and12-layer grading.",
    "plate_hole": "examples/elasticity_plate_hole/make_step.cpp:100x50x2 box minusR5 bore. Prior FE convergence uses its specialized planar-Delaunay sweep, not this main STEP route.",
    "curved_reference": "tests/primitive_cad_reference.cpp + tests/primitive_fixtures.hpp: generated OCCT torusR2/r1,origin(7,-3,4),axis(1,2,3),xDirection(3,-1,.2). Prior checks mesh explicitly parameterized primitives, not automatic STEP meshing.",
    "sphere_gap": "Prior resource-policy scaling returned the same2834-Tet native_surface_star mesh for sphere.5/.25; those were valid Native outputs, not sizing qualification.",
}


def wheel_identity(args):
    receipt = json.loads(args.wheel_receipt.read_text(encoding="utf-8-sig"))
    wheel = args.wheel_receipt.parent / receipt["wheel"]
    if SHARED.sha256(wheel) != receipt["sha256"]:
        raise ValueError("Frozen wheel differs from its validation receipt")
    direct_paths = list(args.installed_env.glob("Lib/site-packages/paacs_mesh-*.dist-info/direct_url.json"))
    if not direct_paths:
        direct_paths = list(args.installed_env.glob("lib/python*/site-packages/paacs_mesh-*.dist-info/direct_url.json"))
    if len(direct_paths) != 1:
        raise ValueError("Expected one installed PAACS direct-url identity")
    direct = json.loads(direct_paths[0].read_text())
    if direct["archive_info"]["hashes"]["sha256"] != receipt["sha256"]:
        raise ValueError("Installed package came from a different wheel")
    return {"wheel_sha256": receipt["sha256"], "wheel_file": wheel.name,
            "wheel_receipt_sha256": SHARED.sha256(args.wheel_receipt),
            "installed_direct_url_sha256": SHARED.sha256(direct_paths[0]),
            "source_receipt": receipt.get("source_revision", {})}


def development_identity(args):
    """Hash the selected binary and the current development source separately."""
    paths = set()
    for directory in (ROOT / "src", ROOT / "include", ROOT / "python/paacs_mesh"):
        paths.update(path for path in directory.rglob("*") if path.is_file()
                     and path.suffix in (".cpp", ".hpp", ".h", ".py"))
    paths.update(path for path in (ROOT / "python").iterdir()
                 if path.is_file() and path.suffix in (".cpp", ".hpp", ".h"))
    files = {path.relative_to(ROOT).as_posix(): SHARED.sha256(path) for path in sorted(paths)}
    snapshot = hashlib.sha256(json.dumps(files, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                            text=True, check=True, timeout=10).stdout.strip()
    status = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True,
                            text=True, check=True, timeout=10).stdout.splitlines()
    return {"native_extension_sha256": SHARED.sha256(args.development_native),
        "native_file": args.development_native.relative_to(ROOT).as_posix()
            if args.development_native.is_relative_to(ROOT) else str(args.development_native),
        "source_commit": commit, "source_worktree_clean": not status,
        "source_status": status, "source_snapshot_sha256": snapshot, "source_files": files,
        "occt_sdk": str(args.occt_sdk),
        "note": "Development binary selected explicitly. Source snapshot and binary identity are recorded separately; this is not an installed-wheel or publication receipt."}


def artifact_identity(args):
    if args.development_native is None:
        return {"artifact_mode": "installed_wheel", "wheel_identity": wheel_identity(args)}
    identity = development_identity(args)
    for key, expected in (("native_extension_sha256", args.expected_native_sha256),
                          ("source_snapshot_sha256", args.expected_source_sha256)):
        if expected is not None and identity[key] != expected:
            raise ValueError(f"Development {key} changed after the serial plan was frozen")
    return {"artifact_mode": "development", "development_identity": identity}


def case_plan(args):
    return [{"case_id": name, **CASES[name], "source_sha256": SHARED.sha256(ROOT / CASES[name]["fixture"])}
            for name in args.cases]


def thickness_measure(mesh, thickness):
    import numpy as np
    snapshot = mesh.to_numpy(max_bytes=256 * 1024**2)
    points = snapshot["points"]
    z = points[:, 2]
    low, high = float(z.min()), float(z.max())
    tolerance = max(1e-12, thickness * 1e-8)
    connectivity = snapshot["cells"]["tet4"]["connectivity"]
    spans = np.ptp(points[connectivity, 2], axis=1) / thickness
    return {"axis": "z", "known_CAD_thickness_mm": thickness, "mesh_extent_mm": high-low,
        "mesh_lower_z": low, "mesh_upper_z": high,
        "interior_z_nodes": int(np.count_nonzero((z > low+tolerance) & (z < high-tolerance))),
        "distinct_z_coordinates_rounded_12dp": int(len(np.unique(np.round(z, 12)))),
        "tet_z_span_fraction": {"median": float(np.median(spans)), "p95": float(np.quantile(spans, .95)),
                                "maximum": float(spans.max())},
        "tets_spanning_at_least_98pct_thickness": int(np.count_nonzero(spans >= .98)),
        "note": "Distinct coordinates are not structured layers. This measures geometry capture, not thin-part solver adequacy."}


def worker(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    spec = CASES[args.case]
    limits = {**LIMITS, "max_seconds": args.native_seconds}
    started, call_started, inspection_started = time.perf_counter(), None, None
    record = {"schema": "paacs-thin-curved-baseline-case-v1", "case_id": args.case, "specification": spec,
        "source_sha256": SHARED.sha256(ROOT / spec["fixture"]), "status": "incomplete", "native_status": "not_attempted",
        "phase": "import_package", "limits": limits, "quality_threshold_arguments": "omitted; installed defaults",
        "timing_scope": args.timing_scope,
        "python": sys.version, "platform": platform.platform(), "retry_count": 0,
        **artifact_identity(args)}

    def save():
        record["worker_seconds"] = time.perf_counter()-started
        SHARED.save(output / "case.json", record)

    save()
    try:
        SHARED.load_development(args)
        import paacs_mesh as pm
        if args.development_native is None:
            if not Path(pm.__file__).resolve().is_relative_to(args.installed_env):
                raise ValueError("Imported a source/development package instead of the selected installed wheel")
        elif (Path(pm.__file__).resolve() != ROOT / "python/paacs_mesh/__init__.py"
              or Path(pm._native.__file__).resolve() != args.development_native):
            raise ValueError("Loaded package or Native extension differs from the explicit development selection")
        record.update(package_version=pm.__version__, package_path=str(Path(pm.__file__).resolve()),
                      native_extension_sha256=SHARED.sha256(Path(pm._native.__file__)))
        record["phase"] = "import"
        imported = time.perf_counter()
        model = pm.read_step(ROOT / spec["fixture"])
        record["import_seconds"] = time.perf_counter()-imported
        faces, bodies = list(model.face_tags), list(model.body_tags)
        selected = None if "face_size" not in spec else model.select_faces(surface="cylinder", radius=spec["radius"], expect=1)
        controls = [] if selected is None else [pm.FaceSize(tag, spec["face_size"]) for tag in selected.tags]
        record["controls"] = {"global_size_mm": spec["size"], "face_sizes": [] if selected is None else
            [{"tag": tag, "size_mm": spec["face_size"]} for tag in selected.tags], "growth": 1.3,
            "resource_policy": "normal", "construction_effort": "legacy", "recovery": "standard", "assembly": "conforming"}
        record["phase"] = "meshing"
        save()
        call_started = time.perf_counter()
        mesh = model.mesh(spec["size"], face_sizes=controls, growth=1.3, recovery="standard", assembly="conforming",
                          resource_policy="normal", construction_effort="legacy", **limits)
        record["call_seconds"] = time.perf_counter()-call_started
        record["peak_process_bytes_after_mesh"] = SHARED.peak_process_bytes()
        record["phase"] = "native_checks"
        report = mesh.report
        SHARED.retain_native_report(output, record, report, "native-report.json")
        if args.development_native is not None and SHARED.sha256(args.development_native) != record["native_extension_sha256"]:
            raise ValueError("Development Native extension changed during meshing")
        record["diagnose"] = pm.diagnose(report)
        record["raw_native_flags"] = {key: report.get(key) for key in
            ("status", "reason", "mesh_published", "audit_valid", "quality_warning", "audit", "quality", "resources", "effective_policy")}
        if (report.get("source_sha256") != record["source_sha256"] or not report.get("mesh_published")
            or not report["audit"]["complete"] or not report["audit"]["valid"]
            or not report["quality"]["complete"] or not report["quality"]["valid"]):
            raise ValueError("Native identity/publication/audit/quality gate failed")
        resource = report["resources"]
        if (resource["policy"] != "normal" or resource["requested_max_work"] is not None
            or resource["logical_work_cap_disabled"] is not True or resource["max_working_bytes"] != limits["max_working_bytes"]):
            raise ValueError("Native resource receipt differs from baseline request")
        HELPERS.verify_tags(faces, bodies, mesh.face_tags, mesh.body_tags, "Live mesh")
        record.update(native_status="succeeded", native_seconds=report["seconds"]["meshing_total"],
                      nodes=report["nodes"], tetrahedra=report["tetrahedra"])
        record["phase"] = "export"
        mesh.write_archive(output / "mesh.paacsmesh")
        record["mesh_sha256"] = SHARED.sha256(output / "mesh.paacsmesh")
        save()
        record["phase"] = "inspection"
        inspection_started = time.perf_counter()
        remaining = min(25., args.hard_seconds-(time.perf_counter()-started)-4.)
        if remaining <= 0:
            raise ValueError("No process deadline remains for independent inspection")
        with pm.inspect_mesh(output / "mesh.paacsmesh", max_cells=limits["max_cells"], max_nodes=limits["max_nodes"],
                             max_bytes=1024**3, max_work=100_000_000, max_seconds=remaining) as inspected:
            summary = inspected.summary()
        SHARED.save(output / "quality-summary.json", summary)
        if (summary["orientation_counts"]["inverted"] or summary["orientation_counts"]["degenerate"]
            or summary["unclassified_tet4_count"] or summary["mesh_sha256"] != record["mesh_sha256"]):
            raise ValueError("Independent archive orientation/ownership/identity screen failed")
        record["archive_tags"] = HELPERS.verify_tags(faces, bodies,
            [int(tag) for tag in summary["sizing"]["achieved"]["boundary_faces"]],
            [int(tag) for tag in summary["body_tet4_counts"]], "Decoded archive")
        record.update(quality=summary["quality"], achieved_sizing=summary["sizing"], bounds=summary["bounds"],
                      orientation_counts=summary["orientation_counts"], metric_conventions=summary["metric_conventions"])
        if "thickness" in spec:
            record["thickness_capture"] = thickness_measure(mesh, spec["thickness"])
        record["inspection_seconds"] = time.perf_counter()-inspection_started
        record["status"], record["phase"] = "succeeded", "complete"
    except Exception as error:
        if record["phase"] == "meshing" and call_started is not None:
            record["call_seconds"] = time.perf_counter()-call_started
        if record["phase"] == "inspection" and inspection_started is not None:
            record["inspection_seconds"] = time.perf_counter()-inspection_started
        record.update(status="meshing_exception_without_report" if record["phase"] == "meshing" else "inspection_or_harness_failed",
                      reason=getattr(error, "reason", None), exception_type=type(error).__name__, exception=str(error))
        if isinstance(getattr(error, "report", None), dict):
            SHARED.retain_native_report(output, record, error.report, "refusal-report.json")
            record["native_status"] = "refused" if error.report.get("mesh_published") is False else "failed_with_report"
            if record["native_status"] == "refused":
                record["status"] = "native_refused"
            try:
                record["diagnose"] = pm.diagnose(error.report)
            except Exception as diagnostic_error:
                record["diagnose_error"] = str(diagnostic_error)
        record["failed_phase"] = record["phase"]
    save()
    print(json.dumps({key: record.get(key) for key in ("case_id", "status", "native_status", "tetrahedra", "call_seconds", "reason")}, allow_nan=False), flush=True)
    return 0 if record["status"] == "succeeded" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--installed-env", type=Path, default=ROOT / "build/quality-feedback/wheel-build/installed-test")
    parser.add_argument("--wheel-receipt", type=Path, default=ROOT / "build/resource-policy/wheel-final/wheel-validation.json")
    parser.add_argument("--development-native", type=Path)
    parser.add_argument("--occt-sdk", type=Path)
    parser.add_argument("--cases", choices=tuple(CASES), nargs="+", default=list(DEFAULT_CASES))
    parser.add_argument("--native-seconds", type=float, default=90.)
    parser.add_argument("--hard-seconds", type=float, default=120.)
    parser.add_argument("--timing-scope", choices=("quiet_single_attempt", "contention_diagnostic"),
                        default="quiet_single_attempt")
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--case", choices=tuple(CASES), help=argparse.SUPPRESS)
    parser.add_argument("--expected-native-sha256", help=argparse.SUPPRESS)
    parser.add_argument("--expected-source-sha256", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.development_native is None:
        if args.occt_sdk is not None:
            parser.error("--occt-sdk requires an explicit --development-native")
        args.installed_env = args.installed_env.resolve(strict=True)
        args.wheel_receipt = args.wheel_receipt.resolve(strict=True)
    else:
        if args.occt_sdk is None:
            parser.error("--development-native requires --occt-sdk")
        args.development_native = args.development_native.resolve(strict=True)
        args.occt_sdk = args.occt_sdk.resolve(strict=True)
        if not args.development_native.is_file() or not args.occt_sdk.is_dir():
            parser.error("Development Native extension and OCCT SDK must exist")
    if (not math.isfinite(args.native_seconds) or not math.isfinite(args.hard_seconds)
        or args.native_seconds <= 0 or args.hard_seconds <= args.native_seconds
        or len(set(args.cases)) != len(args.cases)):
        raise ValueError("Cases must be unique and finite positive process deadline must exceed native deadline")
    if args.worker:
        if args.case is None:
            raise ValueError("Worker case missing")
        return worker(args)
    plan = {"schema": "paacs-thin-curved-baseline-plan-v1", "cases": case_plan(args),
        "limits": {**LIMITS, "max_seconds": args.native_seconds}, "harness_sha256": SHARED.sha256(Path(__file__)),
        "hard_process_seconds": args.hard_seconds, **artifact_identity(args), "provenance": PROVENANCE,
        "timing_scope": args.timing_scope,
        "repetitions": 1, "timing_note": "No medians or repeated timing claim; every process and Native refusal retained."}
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    SHARED.save(args.output / "plan.json", plan)
    if args.plan_only:
        print(json.dumps({"planned": len(plan["cases"]), "native_calls": 0, "output": str(args.output)}, allow_nan=False))
        return 0
    executable = (Path(sys.executable) if args.development_native is not None else
                  args.installed_env / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python"))
    records = []
    for case in plan["cases"]:
        directory = args.output / case["case_id"]
        command = [str(executable), "-I", str(Path(__file__).resolve()), "--worker", "--case", case["case_id"],
            "--output", str(directory), "--installed-env", str(args.installed_env), "--wheel-receipt", str(args.wheel_receipt),
            "--native-seconds", str(args.native_seconds), "--hard-seconds", str(args.hard_seconds),
            "--timing-scope", args.timing_scope]
        if args.development_native is not None:
            identity = plan["development_identity"]
            command += ["--development-native", str(args.development_native), "--occt-sdk", str(args.occt_sdk),
                "--expected-native-sha256", identity["native_extension_sha256"],
                "--expected-source-sha256", identity["source_snapshot_sha256"]]
        began, timed_out = time.perf_counter(), False
        with (args.output / (case["case_id"]+".stdout.log")).open("w") as stdout, (args.output / (case["case_id"]+".stderr.log")).open("w") as stderr:
            try:
                process = subprocess.Popen(command, cwd=ROOT, stdout=stdout, stderr=stderr,
                    start_new_session=sys.platform != "win32")
                code = process.wait(timeout=args.hard_seconds)
            except subprocess.TimeoutExpired:
                timed_out, code = True, None
                # A Windows venv interpreter can have a redirector child. Stop
                # the whole owned process tree before the next serial attempt.
                if sys.platform == "win32":
                    subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                        stdout=stderr, stderr=stderr, check=False, timeout=10)
                else:
                    os.killpg(process.pid, signal.SIGKILL)
                process.wait(timeout=10)
        path = directory / "case.json"
        try:
            result = json.loads(path.read_text()) if path.is_file() else {"case_id": case["case_id"], "status": "missing_worker_record"}
        except (OSError, ValueError) as error:
            result = {"case_id": case["case_id"], "status": "unreadable_worker_record", "exception": str(error)}
        result.update(process_seconds=time.perf_counter()-began, process_exit_code=code, process_timeout=timed_out,
                      command=command, archive_exists=(directory / "mesh.paacsmesh").is_file())
        if timed_out:
            result["status"] = "process_timeout"
        elif code != 0 and result.get("status") == "succeeded":
            result["status"] = "process_failed_after_worker_record"
        records.append(result)
        SHARED.save(args.output / "results.json", records)
        print(json.dumps({key: result.get(key) for key in ("case_id", "status", "native_status", "tetrahedra", "call_seconds", "process_seconds", "reason")}, allow_nan=False), flush=True)
    return 0 if all(row["status"] == "succeeded" for row in records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
