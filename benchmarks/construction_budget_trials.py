"""Compare bounded construction policies on the owned cylinder STEP.

Every case uses one native meshing worker and one validation worker, makes one
attempt, and retains refusals. The optional development-native path loads an
explicit locally compiled extension with its SDK; ordinary use measures the
installed wheel. No validity or geometric tolerance is changed.
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
import statistics
import subprocess
import sys
import time
import types


PROFILES = {
    "legacy": {"construction_effort": "legacy"},
    "explicit-legacy": {"construction_effort": "legacy", "max_surface_points": 256,
                        "max_surface_chord_points": 96, "max_interior_candidates": 16384, "max_interior_sites": 512},
    "double": {"construction_effort": "legacy", "max_surface_points": 512,
               "max_surface_chord_points": 192, "max_interior_candidates": 32768, "max_interior_sites": 1024},
    "quadruple": {"construction_effort": "legacy", "max_surface_points": 1024,
                  "max_surface_chord_points": 384, "max_interior_candidates": 65536, "max_interior_sites": 2048},
    "fair-legacy": {"construction_effort": "auto", "max_surface_points": 256,
                    "max_surface_chord_points": 96, "max_interior_candidates": 16384, "max_interior_sites": 512},
    "auto": {"construction_effort": "auto"},
    "thorough": {"construction_effort": "thorough"},
}
# Explicitly freeze baseline policy when the public default later changes.
for _profile_options in PROFILES.values():
    _profile_options["resource_policy"] = "bounded"
PROFILES["normal"] = {"construction_effort": "legacy", "resource_policy": "normal"}
PROFILES["normal-quadruple"] = dict(PROFILES["quadruple"], resource_policy="normal")


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8192), b""):
            digest.update(block)
    return digest.hexdigest()


def peak_process_bytes():
    """Process high-water observation, distinct from modeled native admission."""
    if os.name != "nt":
        import resource
        value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return int(value if sys.platform == "darwin" else value * 1024)
    import ctypes
    from ctypes import wintypes
    class Counters(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("faults", wintypes.DWORD)] + [
            (name, ctypes.c_size_t) for name in ("peak_working_set", "working_set", "peak_paged", "paged",
                                                "peak_nonpaged", "nonpaged", "pagefile", "peak_pagefile")]
    counters = Counters()
    counters.cb = ctypes.sizeof(counters)
    if ctypes.windll.psapi.GetProcessMemoryInfo(ctypes.c_void_p(-1), ctypes.byref(counters), counters.cb):
        return counters.peak_working_set
    return None


def load_development(args):
    if args.development_native is None:
        return
    if args.occt_sdk is None or not args.development_native.is_file():
        raise ValueError("Development loading requires an existing extension and --occt-sdk")
    package = Path(__file__).resolve().parents[1] / "python/paacs_mesh"
    if "paacs_mesh" in sys.modules:
        raise RuntimeError("Development extension must be selected before importing paacs_mesh")
    handles = []
    if os.name == "nt":
        for directory in (args.occt_sdk / "bin", args.development_native.parent.parent.parent / "Release"):
            if directory.is_dir():
                handles.append(os.add_dll_directory(str(directory.resolve())))
    resources = args.occt_sdk / "share/opencascade/resources"
    for variable, suffix in (("CSF_XSMessage", "XSMessage"), ("CSF_SHMessage", "SHMessage"),
                             ("CSF_STEPDefaults", "XSTEPResource"), ("CSF_IGESDefaults", "XSTEPResource")):
        if (resources / suffix).is_dir():
            os.environ[variable] = str((resources / suffix).resolve())
    module = types.ModuleType("paacs_mesh")
    module.__file__ = str(package / "__init__.py")
    module.__path__ = [str(package)]
    sys.modules["paacs_mesh"] = module
    spec = importlib.util.spec_from_file_location("paacs_mesh._native", args.development_native)
    native = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = native
    spec.loader.exec_module(native)
    module._native = native
    exec(compile((package / "__init__.py").read_text(encoding="utf-8"), module.__file__, "exec"), module.__dict__)
    module._development_dll_handles = handles


def retain_native_report(output, record, report, filename):
    """Keep the complete receipt on success and refusal; never infer missing work."""
    save(output / filename, report)
    record["native_report_file"] = filename
    record["native_report_sha256"] = sha256(output / filename)
    record["native_reported_status"] = report.get("status")
    record["native_reported_reason"] = report.get("reason")
    record["native_reported_seconds"] = report.get("seconds")
    record["native_reported_validation_work"] = report.get("validation_work")
    record["native_resources"] = report.get("resources")
    record["native_audit"] = report.get("audit")
    record["native_quality"] = report.get("quality")
    record["native_stage_counters"] = [{key: body.get(key) for key in (
        "solid_tag", "route", "reason", "native_surface_work", "native_volume_work",
        "native_surface_seconds", "native_volume_seconds", "initial_front", "surface_spacing")}
        for body in report.get("bodies", [])]


def worker(args):
    load_development(args)
    import paacs_mesh as pm

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    meshing_started = None
    record = {"schema": "paacs-construction-budget-case-v1", "profile": args.profile,
              "wall_size_mm": args.wall_size, "global_size_mm": args.size,
              "source_sha256": sha256(args.source), "source": str(args.source.resolve()),
              "python": sys.version, "host": platform.platform(), "package_version": pm.__version__,
              "native_extension_sha256": sha256(Path(pm._native.__file__)),
              "options": PROFILES[args.profile], "workers": 1, "validation_workers": 1,
              "retry_count": 0, "status": "incomplete", "native_status": "not_attempted",
              "global_max_work": args.max_work, "global_max_working_bytes": args.max_working_bytes,
              "repetition": args.repetition,
              "requested_controls": {"growth": 1.3, "wall_size_mm": args.wall_size,
                  "global_size_mm": args.size},
              "other_global_limits": {"max_nodes": 500000, "max_cells": 3000000,
                  "max_output_bytes": 256 * 1024**2, "max_audit_bytes": 512 * 1024**2,
                  "max_seconds": args.max_seconds}}
    try:
        record["phase"] = "import"
        imported = time.perf_counter()
        model = pm.read_step(args.source)
        record["import_seconds"] = time.perf_counter() - imported
        record["phase"] = "selection"
        wall = None if args.wall_size is None else model.select_faces(surface="cylinder", radius=1, expect=1)
        controls = [] if wall is None else [pm.FaceSize(tag, args.wall_size) for tag in wall.tags]
        record["phase"] = "meshing"
        meshing_started = time.perf_counter()
        mesh = model.mesh(args.size, face_sizes=controls, assembly="conforming",
            workers=1, validation_workers=1, growth=1.3,
            max_nodes=500000, max_cells=3000000, max_work=args.max_work,
            max_output_bytes=256 * 1024**2, max_working_bytes=args.max_working_bytes,
            max_audit_bytes=512 * 1024**2, max_seconds=args.max_seconds, **PROFILES[args.profile])
        record["call_seconds"] = time.perf_counter() - meshing_started
        record["peak_process_bytes_after_mesh"] = peak_process_bytes()
        record["phase"] = "native_checks"
        record["native_status"] = "returned_unchecked"
        report = mesh.report
        retain_native_report(output, record, report, "native-report.json")
        record["native_status"] = "checks_failed"
        if report.get("source_sha256") != record["source_sha256"]:
            raise RuntimeError("Imported source identity changed during the trial")
        if not (report["mesh_published"] and report["audit"]["complete"] and report["audit"]["valid"]
                and report["quality"]["complete"] and report["quality"]["valid"]):
            raise RuntimeError("Native mesh did not complete the unchanged publication checks")
        record["native_status"] = "succeeded"
        record["phase"] = "export"
        exported = time.perf_counter()
        mesh.write_archive(output / "mesh.paacsmesh")
        record["export_seconds"] = time.perf_counter() - exported
        record["mesh_sha256"] = sha256(output / "mesh.paacsmesh")
        record["nodes"], record["tetrahedra"] = report["nodes"], report["tetrahedra"]
        record["native_seconds"] = report["seconds"]["meshing_total"]
        record["phase"] = "analysis"
        measured = time.perf_counter()
        with pm.inspect_mesh(mesh, max_cells=3000000, max_bytes=512 * 1024**2) as inspected:
            summary = inspected.summary()
            save(output / "quality-summary.json", summary)
            record["quality"] = summary["quality"]
            record["below_mean_ratio_0p1"] = inspected.query(quality_below=.1, limit=1, associations=False)["total"]
        reference_path = Path(__file__).with_name("refinement_diagnostics.py")
        spec = importlib.util.spec_from_file_location("construction_trial_reference", reference_path)
        reference = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reference)
        measurements = reference.measure(output / "mesh.paacsmesh")
        save(output / "independent-measurements.json", measurements)
        record["measurement_seconds"] = time.perf_counter() - measured
        record["wall_sizes"] = None if wall is None else measurements["face_unique_edge_lengths_mm"][str(wall.tags[0])]
        record["status"] = "succeeded"
        record["phase"] = "complete"
    except Exception as error:
        if record["phase"] == "meshing" and meshing_started is not None:
            record["call_seconds"] = time.perf_counter() - meshing_started
            record["peak_process_bytes_after_mesh"] = peak_process_bytes()
            record["native_status"] = "failed_without_native_report"
        record.update(status="failed", reason=getattr(error, "reason", None),
                      stage=getattr(error, "stage", None), exception=str(error), exception_type=type(error).__name__)
        if isinstance(getattr(error, "report", None), dict):
            retain_native_report(output, record, error.report, "refusal-report.json")
            if record["phase"] == "meshing":
                record["native_status"] = "refused" if error.report.get("mesh_published") is False else "failed_with_report"
    record["worker_seconds"] = time.perf_counter() - started
    record["peak_process_bytes_after_analysis"] = peak_process_bytes()
    save(output / "case.json", record)
    print(json.dumps(record, allow_nan=False))
    return 0 if record["status"] == "succeeded" else 1


def trial_plan(args):
    """Preserve old single-run names; qualify independent grids explicitly."""
    work_values = args.work_values or [args.max_work]
    memory_values = args.working_byte_values or [args.max_working_bytes]
    expanded = bool(args.work_values or args.working_byte_values or args.repetitions > 1)
    count = len(args.wall_sizes) * len(args.profiles) * len(work_values) * len(memory_values) * args.repetitions
    if count > args.max_cases:
        raise ValueError(f"Trial plan contains {count} cases, exceeding --max-cases {args.max_cases}")
    rows = []
    for wall in args.wall_sizes:
        for profile in args.profiles:
            for work in work_values:
                for memory in memory_values:
                    for repetition in range(1, args.repetitions + 1):
                        name = profile + "-wall-" + ("none" if wall is None else str(wall).replace(".", "p"))
                        if expanded:
                            name += (f"-work-{'default' if work is None else work}"
                                     f"-memory-{'default' if memory is None else memory}-repeat-{repetition:03d}")
                        rows.append(dict(name=name, profile=profile, wall_size_mm=wall,
                                         max_work=work, max_working_bytes=memory, repetition=repetition))
    return rows


def repeated_summary(rows, requested_repetitions):
    """Report all attempts; mesh-time medians require three identical valid meshes."""
    groups = {}
    for row in rows:
        key = (row["profile"], row["wall_size_mm"], row["global_max_work"], row["global_max_working_bytes"])
        groups.setdefault(key, []).append(row)
    result = []
    for key, attempts in groups.items():
        successes = [r for r in attempts if r.get("status") == "succeeded"]
        complete = (len(attempts) == requested_repetitions
                    and sorted(r.get("repetition", 0) for r in attempts) == list(range(1, requested_repetitions + 1)))
        hashes = sorted({r["mesh_sha256"] for r in successes if r.get("mesh_sha256")})
        identities_match = all(len({r.get(field) for r in attempts}) == 1
                               and attempts[0].get(field) is not None
                               for field in ("source_sha256", "native_extension_sha256"))
        same_requests = all(len({json.dumps(r.get(field), sort_keys=True, allow_nan=False) for r in attempts}) == 1
                            for field in ("global_size_mm", "options", "requested_controls", "other_global_limits"))
        comparable = (complete and len(attempts) >= 3 and len(successes) == len(attempts)
                      and len(hashes) == 1 and all(r.get("mesh_sha256") for r in attempts)
                      and identities_match and same_requests)
        result.append({"profile": key[0], "wall_size_mm": key[1], "global_max_work": key[2],
            "global_max_working_bytes": key[3], "attempts": len(attempts), "successful_attempts": len(successes),
            "failed_attempts": len(attempts) - len(successes), "requested_repetitions_complete": complete,
            "all_attempt_process_seconds": [r["process_seconds"] for r in attempts],
            "all_successful_mesh_sha256": hashes,
            "attempt_process_median_seconds": statistics.median(r["process_seconds"] for r in attempts)
                if complete and len(attempts) >= 3 else None,
            "native_mesh_median_seconds": statistics.median(r["native_seconds"] for r in attempts)
                if comparable else None,
            "native_median_basis": "three_or_more_complete_identical_valid_mesh_repetitions" if comparable else "unavailable",
            "timing_limitations": "Serial observations; no hardware isolation, randomized order or simulation accuracy claim."})
    return result


def controller(args):
    plan = trial_plan(args)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    rows = []
    save(output / "trial-plan.json", {"schema": "paacs-resource-trial-plan-v1",
        "source_sha256": sha256(args.source), "source": str(args.source.resolve()),
        "profiles": {name: PROFILES[name] for name in args.profiles}, "repetitions": args.repetitions,
        "order": "serial wall, profile, work, memory, repetition", "cases": plan})
    if args.plan_only:
        print(json.dumps({"status": "planned_only", "cases": len(plan), "native_calls": 0,
                          "plan": str(output / "trial-plan.json")}))
        return 0
    for planned in plan:
        name, wall, profile = planned["name"], planned["wall_size_mm"], planned["profile"]
        case = output / name
        command = [sys.executable, "-I", str(Path(__file__).resolve()), "--worker", "--profile", profile,
                   "--source", str(args.source.resolve()), "--size", str(args.size),
                   "--output", str(case), "--wall-size", "none" if wall is None else str(wall),
                   "--max-seconds", str(args.max_seconds), "--hard-seconds", str(args.hard_seconds),
                   "--max-work", str(planned["max_work"]),
                   "--max-working-bytes", str(planned["max_working_bytes"]),
                   "--repetition", str(planned["repetition"])]
        if args.development_native is not None:
            command += ["--development-native", str(args.development_native.resolve()),
                        "--occt-sdk", str(args.occt_sdk.resolve())]
        started = time.perf_counter()
        with (output / (name + ".stdout.log")).open("wb") as stdout, (output / (name + ".stderr.log")).open("wb") as stderr:
            try:
                process = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=args.hard_seconds, check=False)
                code, timed_out = process.returncode, False
            except subprocess.TimeoutExpired:
                code, timed_out = None, True
        elapsed = time.perf_counter() - started
        row = {"status": "timeout" if timed_out else "process_failed"}
        if (case / "case.json").is_file():
            try:
                loaded = json.loads((case / "case.json").read_text())
                if not isinstance(loaded, dict):
                    raise ValueError("Worker receipt must be a JSON object")
                row = loaded
            except (ValueError, UnicodeError, OSError) as error:
                row.update(reason="worker_receipt_unreadable", exception=str(error))
        if (timed_out or code != 0) and row.get("status") == "succeeded":
            row["worker_reported_status"] = row["status"]
            row["status"] = "timeout" if timed_out else "process_failed_after_result"
        row.update(profile=profile, wall_size_mm=wall, global_max_work=planned["max_work"],
                   global_max_working_bytes=planned["max_working_bytes"], repetition=planned["repetition"],
                   process_seconds=elapsed, process_exit_code=code, process_timeout=timed_out)
        rows.append(row)
        save(output / "summary.json", rows)
        save(output / "repetitions.json", repeated_summary(rows, args.repetitions))
        print(json.dumps({key: row.get(key) for key in
              ("profile", "wall_size_mm", "global_max_work", "global_max_working_bytes", "repetition", "status",
               "tetrahedra", "native_seconds", "below_mean_ratio_0p1", "reason")}), flush=True)
    return 0 if all(row["status"] == "succeeded" for row in rows) else 1


def wall_size(text):
    if text.lower() == "none":
        return None
    value = float(text)
    if not math.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError("wall size must be positive and finite, or none")
    return value


def resource_cap(text):
    if text.lower() == "none":
        return None
    try:
        value = int(text)
    except ValueError as error:
        raise argparse.ArgumentTypeError("resource cap must be a positive integer or none") from error
    if not 1 <= value <= 2**63 - 1:
        raise argparse.ArgumentTypeError("resource cap must fit a positive signed 64-bit integer")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1] / "tests/fixtures/cad/cylinder.step")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profiles", choices=tuple(PROFILES), nargs="+", default=["legacy", "double", "quadruple", "auto"])
    parser.add_argument("--wall-sizes", type=wall_size, nargs="+", default=[.2, .1, .08])
    parser.add_argument("--size", type=float, default=.5)
    parser.add_argument("--max-seconds", type=float, default=90)
    parser.add_argument("--hard-seconds", type=float, default=150)
    parser.add_argument("--max-work", type=resource_cap, default=100000000,
                        help="Overall work cap; none inherits native policy, not necessarily unlimited")
    parser.add_argument("--max-working-bytes", type=resource_cap, default=1024**3,
                        help="Modeled working-memory cap; none inherits native policy")
    parser.add_argument("--work-values", type=resource_cap, nargs="+", help="Independent grid of overall logical-work caps")
    parser.add_argument("--working-byte-values", type=resource_cap, nargs="+", help="Independent grid of modeled working-memory caps")
    parser.add_argument("--repetitions", type=int, default=1)
    parser.add_argument("--max-cases", type=int, default=256)
    parser.add_argument("--plan-only", action="store_true", help="Write the bounded run plan without importing or meshing CAD")
    parser.add_argument("--development-native", type=Path)
    parser.add_argument("--occt-sdk", type=Path)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--profile", choices=tuple(PROFILES), default="legacy", help=argparse.SUPPRESS)
    parser.add_argument("--wall-size", type=wall_size, default=.1, help=argparse.SUPPRESS)
    parser.add_argument("--repetition", type=int, default=1, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if (not args.source.is_file() or not all(math.isfinite(v) for v in (args.size, args.max_seconds, args.hard_seconds))
            or args.size <= 0 or args.max_seconds <= 0 or args.hard_seconds <= args.max_seconds
            or any(value is not None and value <= 0 for value in (args.max_work, args.max_working_bytes))):
        parser.error("source must exist, size/deadlines must be finite and positive, and hard deadline must exceed native deadline")
    if args.development_native is not None and args.occt_sdk is None:
        parser.error("--development-native requires --occt-sdk")
    if not (1 <= args.repetitions <= 100 and 1 <= args.repetition <= 100 and 1 <= args.max_cases <= 4096):
        parser.error("Repetition counts must be 1..100 and max-cases must be 1..4096")
    for name, values in (("profiles", args.profiles), ("wall sizes", args.wall_sizes),
                         ("work values", args.work_values or []), ("working-byte values", args.working_byte_values or [])):
        if len(set(values)) != len(values):
            parser.error(f"Duplicate {name} would overwrite case evidence")
    if not args.worker:
        try:
            trial_plan(args)
        except ValueError as error:
            parser.error(str(error))
    return worker(args) if args.worker else controller(args)


if __name__ == "__main__":
    raise SystemExit(main())
