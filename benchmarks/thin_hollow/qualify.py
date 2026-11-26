"""Serial public-Python meshing trials on the owned hollow-extrusion fixtures.

Every attempt, including a refusal or timeout, has a separate receipt. The
controller does not import the native extension. Lengths and sizes are in mm.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

PROFILES = {
    "auto-order": {"thin_region_ordering": True},
    "default": {},
    "boundary-first": {"boundary_first": True},
    "front": {"advancing_interior_front": True},
    "front-boundary-first": {"advancing_interior_front": True, "boundary_first": True},
    "front-auto-order": {"advancing_interior_front": True, "thin_region_ordering": True},
    "front-source-fallback": {"advancing_interior_front": True, "native_source_fallback": True},
    "legacy": {"advancing_interior_front": False},
    "advanced": {"recovery": "advanced", "boundary_first": True, "native_source_fallback": True},
}

def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()

def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))

def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(".json.tmp")
    for attempt in range(8):
        try:
            temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n", encoding="utf-8")
            os.replace(temporary, path)
            return
        except OSError as error:
            if getattr(error, "winerror", None) not in (5, 32, 33) or attempt == 7:
                raise
            time.sleep(.01*(attempt+1))

def now():
    return datetime.now(timezone.utc).isoformat()

def guards(job):
    for path, digest in job["input_sha256"].items():
        if sha(path) != digest:
            raise ValueError("Pinned input changed: "+path)
    for name, digest in job["stage_files"].items():
        if sha(Path(job["stage"])/name) != digest:
            raise ValueError("Staged package changed: "+name)

def peak_bytes(process):
    if os.name != "nt":
        return None
    import ctypes
    class Counters(ctypes.Structure):
        _fields_ = [("cb", ctypes.c_ulong), ("faults", ctypes.c_ulong)]+[(n, ctypes.c_size_t)for n in
            ("peak", "working", "peak_paged", "paged", "peak_nonpaged", "nonpaged", "pagefile", "peak_pagefile", "private")]
    counters = Counters()
    counters.cb = ctypes.sizeof(counters)
    function = ctypes.WinDLL("psapi").GetProcessMemoryInfo
    function.argtypes = [ctypes.c_void_p, ctypes.POINTER(Counters), ctypes.c_ulong]
    function.restype = ctypes.c_int
    return int(counters.peak) if function(int(process._handle), ctypes.byref(counters), counters.cb) else None

def worker(path):
    job = load(path)
    folder = Path(job["output"])
    receipt = {"pid": os.getpid(), "status": "running", "phase": "guards", "started_utc": now(),
        "source_sha256": job["source_sha256"], "native_sha256": job["native_sha256"], "native_status": "not_started"}
    began = time.perf_counter()
    dll_scope = None
    save(job["worker_receipt"], receipt)
    try:
        guards(job)
        if job["sdk"]:
            sdk = Path(job["sdk"])
            if os.name == "nt":
                dll_scope = os.add_dll_directory(str(sdk/"bin"))
            for variable, directory in (("CSF_XSMessage", "XSMessage"), ("CSF_SHMessage", "SHMessage"),
                    ("CSF_STEPDefaults", "XSTEPResource"), ("CSF_IGESDefaults", "XSTEPResource")):
                os.environ[variable] = str(sdk/"share/opencascade/resources"/directory)
        sys.path.insert(0, job["stage"])
        at = time.perf_counter()
        import paacs_mesh as pm
        if Path(pm.__file__).resolve().parent != Path(job["stage"])/"paacs_mesh":
            raise ValueError("Unexpected Python package")
        receipt.update(package_import_seconds=time.perf_counter()-at, phase="step_import")
        save(job["worker_receipt"], receipt)
        at = time.perf_counter()
        model = pm.read_step(job["source"])
        receipt.update(step_import_seconds=time.perf_counter()-at, original_body_tags=model.body_tags,
            original_face_tags=model.face_tags, model_identity=model.model_identity)
        if model.source_sha256 != job["source_sha256"] or len(model.body_tags) != 1:
            raise ValueError("Expected the unchanged single-solid fixture")
        save(folder/"model-info.json", model.info)
        save(folder/"body-inspection.json", model.inspect("bodies", limit=100))
        receipt["phase"] = "meshing"
        save(job["worker_receipt"], receipt)
        at = time.perf_counter()
        report = None
        try:
            mesh = model.mesh(job["global_h_mm"], **job["kwargs"])
            report = mesh.report
            receipt.update(native_status="succeeded", status="native_accepted")
        except Exception as error:
            report = getattr(error, "report", None)
            receipt.update(native_status="refused" if isinstance(report, dict) else "unknown",
                status="native_refused" if isinstance(report, dict) else "infrastructure_failure",
                exception_type=type(error).__name__, exception=str(error), reason=getattr(error, "reason", None))
        finally:
            receipt["mesh_call_seconds"] = time.perf_counter()-at
        if isinstance(report, dict):
            save(folder/"native-report.json", report)
            receipt.update(native_report_file="native-report.json", native_report_sha256=sha(folder/"native-report.json"))
        if receipt["native_status"] == "succeeded":
            if not (report["mesh_published"] and report["source_sha256"] == job["source_sha256"]
                    and all(report[s][k] is True for s in ("audit", "quality") for k in ("complete", "valid"))):
                raise ValueError("Native publication or validity gate failed")
            if sorted(mesh.body_tags) != sorted(model.body_tags) or sorted(mesh.face_tags) != sorted(model.face_tags):
                raise ValueError("Published mesh lost CAD tags")
            receipt["phase"] = "export"
            save(job["worker_receipt"], receipt)
            at = time.perf_counter()
            mesh.write_archive(folder/"mesh.paacsmesh")
            mesh.write_vtk(folder/"mesh.vtk")
            receipt.update(export_seconds=time.perf_counter()-at, archive_sha256=sha(folder/"mesh.paacsmesh"))
        guards(job)
        model._assert_identity()
        receipt.update(input_guards_unchanged=True, phase="complete")
    except BaseException as error:
        receipt.update(status="infrastructure_failure", error_type=type(error).__name__, error=str(error))
    finally:
        receipt.update(worker_seconds=time.perf_counter()-began, finished_utc=now())
        save(job["worker_receipt"], receipt)
        if dll_scope is not None:
            dll_scope.close()
    return 0 if receipt["status"] in ("native_accepted", "native_refused") else 1

def run(job):
    folder = Path(job["output"])
    folder.mkdir(parents=True, exist_ok=False)
    save(folder/"job.json", job)
    command = [sys.executable, "-I", str(Path(__file__).resolve()), "--worker", str(folder/"job.json")]
    receipt = {"schema": "paacs-hollow-extrusion-controller-v1", "job": job, "command": command,
        "runner_sha256": sha(__file__), "started_utc": now(), "state": "running"}
    save(folder/"controller-receipt.json", receipt)
    began = time.perf_counter()
    process, peak, hard = None, None, False
    try:
        with (folder/"worker.log").open("wb") as stream:
            process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            receipt["pid"] = process.pid
            save(folder/"controller-receipt.json", receipt)
            while process.poll() is None:
                measured = peak_bytes(process)
                if measured is not None:
                    peak = max(peak or 0, measured)
                if time.perf_counter()-began > job["hard_seconds"]:
                    hard = True
                    process.kill()
                    break
                time.sleep(.2)
            process.wait()
    finally:
        if process is not None and process.poll() is None:
            process.kill()
            process.wait()
    observed = load(job["worker_receipt"]) if Path(job["worker_receipt"]).exists() else {}
    receipt.update(worker=observed, state="complete", process_seconds=time.perf_counter()-began,
        returncode=process.returncode, hard_timeout=hard, worker_pid_verified=observed.get("pid") == process.pid,
        peak_process_working_set_bytes=peak, finished_utc=now())
    guards(job)
    receipt["input_guards_unchanged"] = True
    receipt["status"] = "hard_timeout_native_outcome_unknown" if hard else observed.get("status", "missing_receipt")
    save(folder/"controller-receipt.json", receipt)
    return receipt

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--stage", type=Path)
    parser.add_argument("--sdk", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cases", nargs="+")
    parser.add_argument("--profiles", nargs="+", choices=PROFILES, default=["default"])
    parser.add_argument("--assembly", choices=["preserve_parts", "conforming"], default="preserve_parts")
    parser.add_argument("--sizes", nargs="+", type=float, default=[4.])
    parser.add_argument("--chord", type=float)
    parser.add_argument("--memory-gib", type=float, default=2.)
    parser.add_argument("--seconds", type=float, default=180.)
    parser.add_argument("--plan-only", action="store_true")
    args = parser.parse_args()
    if args.worker:
        return worker(args.worker)
    if not all((args.manifest, args.stage, args.output)) or not args.cases:
        parser.error("Supply --manifest, --stage, --output and --cases")
    if any(not 0 < v < float("inf") for v in (*args.sizes, args.memory_gib, args.seconds)):
        parser.error("Sizes and resource limits must be finite and positive")
    if args.chord is not None and not 0 < args.chord < float("inf"):
        parser.error("Chord tolerance must be finite and positive")
    manifest_path, stage = args.manifest.resolve(strict=True), args.stage.resolve(strict=True)
    manifest = load(manifest_path)
    records = manifest["fixtures"]
    by_id = {c["id"]: c for c in records}
    native_files = [p for p in (stage/"paacs_mesh").iterdir() if p.name.startswith("_native.") and p.suffix in (".pyd", ".so")]
    if len(native_files) != 1:
        raise ValueError("Expected one staged native extension")
    stage_files = {str(p.relative_to(stage)): sha(p)for p in sorted((stage/"paacs_mesh").rglob("*"))
        if p.is_file() and p.suffix in (".py", ".pyd", ".so", ".dll")}
    runtime_files = {str(p.resolve()): sha(p)for p in args.sdk.joinpath("bin").glob("*.dll")} if args.sdk else {}
    output = args.output.resolve()
    if output.exists():
        raise ValueError("Choose a fresh output directory")
    jobs = []
    for name in args.cases:
        if Path(name).name != name or any(c in name for c in ("/", "\\", ":")):
            raise ValueError("Fixture IDs must be single path components")
        case = by_id[name]
        source = (manifest_path.parent/case["step"]).resolve(strict=True)
        if not source.is_relative_to(manifest_path.parent):
            raise ValueError("Fixture STEP must be beside its owned manifest")
        for size in args.sizes:
            for profile in args.profiles:
                directory = output/f"{name}-h{size:g}-{profile}"
                kwargs = {"assembly": args.assembly, "workers": 1, "validation_workers": 1,
                    "resource_policy": "normal", "construction_effort": "legacy", "max_seconds": args.seconds,
                    "max_nodes": 500000, "max_cells": 3000000, "max_output_bytes": 256*1024**2,
                    "max_audit_bytes": 512*1024**2, "max_working_bytes": int(args.memory_gib*1024**3),
                    **PROFILES[profile]}
                if args.chord is not None:
                    kwargs.update(max_chord_error=args.chord, max_relative_volume_error=.01)
                job = {"schema": "paacs-hollow-extrusion-job-v1", "id": name, "profile": profile, "specification": case,
                    "source": str(source), "source_sha256": sha(source), "global_h_mm": size, "kwargs": kwargs,
                    "stage": str(stage), "stage_files": stage_files, "native_sha256": sha(native_files[0]),
                    "sdk": str(args.sdk.resolve(strict=True)) if args.sdk else None,
                    "input_sha256": {str(manifest_path): sha(manifest_path), str(source): sha(source), str(Path(__file__).resolve()): sha(__file__), **runtime_files},
                    "hard_seconds": args.seconds+60, "output": str(directory), "worker_receipt": str(directory/"worker-receipt.json")}
                jobs.append(job)
    if args.plan_only:
        print(json.dumps({"native_calls": 0, "jobs": jobs}, indent=2))
        return 0
    output.mkdir(parents=True, exist_ok=False)
    save(output/"plan.json", {"schema": "paacs-hollow-extrusion-plan-v1", "jobs": jobs})
    results = []
    for job in jobs:
        receipt = run(job)
        results.append({"id": job["id"], "profile": job["profile"], "size": job["global_h_mm"],
            "status": receipt["status"], "process_seconds": receipt["process_seconds"], "receipt": str(Path(job["output"])/"controller-receipt.json")})
        save(output/"results.json", {"complete": len(results) == len(jobs), "attempts": results,
            "total_process_seconds": sum(r["process_seconds"]for r in results)})
        print(json.dumps(results[-1]), flush=True)
        if receipt["status"] in ("infrastructure_failure", "missing_receipt"):
            raise RuntimeError("Trial infrastructure failed; remaining requests were not executed")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
