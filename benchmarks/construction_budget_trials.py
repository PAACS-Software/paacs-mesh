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


def worker(args):
    load_development(args)
    import paacs_mesh as pm

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    record = {"schema": "paacs-construction-budget-case-v1", "profile": args.profile,
              "wall_size_mm": args.wall_size, "global_size_mm": args.size,
              "source_sha256": sha256(args.source), "source": str(args.source.resolve()),
              "python": sys.version, "host": platform.platform(), "package_version": pm.__version__,
              "native_extension_sha256": sha256(Path(pm._native.__file__)),
              "options": PROFILES[args.profile], "workers": 1, "validation_workers": 1,
              "retry_count": 0, "status": "incomplete", "native_status": "not_attempted",
              "global_max_work": args.max_work, "global_max_working_bytes": args.max_working_bytes}
    try:
        imported = time.perf_counter()
        model = pm.read_step(args.source)
        record["import_seconds"] = time.perf_counter() - imported
        wall = None if args.wall_size is None else model.select_faces(surface="cylinder", radius=1, expect=1)
        controls = [] if wall is None else [pm.FaceSize(tag, args.wall_size) for tag in wall.tags]
        meshing = time.perf_counter()
        mesh = model.mesh(args.size, face_sizes=controls, assembly="conforming",
            workers=1, validation_workers=1, growth=1.3,
            max_nodes=500000, max_cells=3000000, max_work=args.max_work,
            max_output_bytes=256 * 1024**2, max_working_bytes=args.max_working_bytes,
            max_audit_bytes=512 * 1024**2, max_seconds=args.max_seconds, **PROFILES[args.profile])
        record["call_seconds"] = time.perf_counter() - meshing
        report = mesh.report
        assert report["mesh_published"] and report["audit"]["complete"] and report["audit"]["valid"]
        assert report["quality"]["complete"] and report["quality"]["valid"]
        record["native_status"] = "succeeded"
        record["peak_process_bytes_after_mesh"] = peak_process_bytes()
        save(output / "native-report.json", report)
        exported = time.perf_counter()
        mesh.write_archive(output / "mesh.paacsmesh")
        record["export_seconds"] = time.perf_counter() - exported
        record["mesh_sha256"] = sha256(output / "mesh.paacsmesh")
        record["nodes"], record["tetrahedra"] = report["nodes"], report["tetrahedra"]
        record["native_seconds"] = report["seconds"]["meshing_total"]
        record["native_quality"] = report["quality"]
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
    except Exception as error:
        record.update(status="failed", reason=getattr(error, "reason", None),
                      stage=getattr(error, "stage", None), exception=str(error), exception_type=type(error).__name__)
        if isinstance(getattr(error, "report", None), dict):
            save(output / "refusal-report.json", error.report)
    record["worker_seconds"] = time.perf_counter() - started
    record["peak_process_bytes_after_analysis"] = peak_process_bytes()
    save(output / "case.json", record)
    print(json.dumps(record, allow_nan=False))
    return 0 if record["status"] == "succeeded" else 1


def controller(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    rows = []
    for wall in args.wall_sizes:
        for profile in args.profiles:
            name = profile + "-wall-" + ("none" if wall is None else str(wall).replace(".", "p"))
            case = output / name
            command = [sys.executable, "-I", str(Path(__file__).resolve()), "--worker", "--profile", profile,
                       "--source", str(args.source.resolve()), "--size", str(args.size),
                       "--output", str(case), "--wall-size", "none" if wall is None else str(wall),
                       "--max-seconds", str(args.max_seconds), "--max-work", str(args.max_work),
                       "--max-working-bytes", str(args.max_working_bytes)]
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
            row = json.loads((case / "case.json").read_text()) if (case / "case.json").is_file() else {
                "profile": profile, "wall_size_mm": wall, "status": "timeout" if timed_out else "process_failed"}
            row.update(process_seconds=elapsed, process_exit_code=code, process_timeout=timed_out)
            rows.append(row)
            save(output / "summary.json", rows)
            print(json.dumps({key: row.get(key) for key in
                  ("profile", "wall_size_mm", "status", "tetrahedra", "native_seconds", "below_mean_ratio_0p1", "reason")}))
    return 0 if all(row["status"] == "succeeded" for row in rows) else 1


def wall_size(text):
    if text.lower() == "none":
        return None
    value = float(text)
    if not math.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError("wall size must be positive and finite, or none")
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
    parser.add_argument("--max-work", type=int, default=100000000)
    parser.add_argument("--max-working-bytes", type=int, default=1024**3)
    parser.add_argument("--development-native", type=Path)
    parser.add_argument("--occt-sdk", type=Path)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--profile", choices=tuple(PROFILES), default="legacy", help=argparse.SUPPRESS)
    parser.add_argument("--wall-size", type=wall_size, default=.1, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if (not args.source.is_file() or not all(math.isfinite(v) for v in (args.size, args.max_seconds, args.hard_seconds))
            or args.size <= 0 or args.max_seconds <= 0 or args.hard_seconds <= args.max_seconds
            or args.max_work <= 0 or args.max_working_bytes <= 0):
        parser.error("source must exist, size/deadlines must be finite and positive, and hard deadline must exceed native deadline")
    if args.development_native is not None and args.occt_sdk is None:
        parser.error("--development-native requires --occt-sdk")
    return worker(args) if args.worker else controller(args)


if __name__ == "__main__":
    raise SystemExit(main())
