#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Qualify explicit STEP inputs through an installed PAACS Python distribution.

This controller imports no mesher. Each attempt runs serially under the selected
Python with -I, without source paths, CAD SDK paths or loader overrides. The plan
pins paths, input hashes, version and all attempt options. Native checks remain
authoritative for CAD embedding; independent inspection screens encoded Tet4
shape, ownership and export identity, not a new CAD audit.
"""
import argparse
import base64
import ctypes
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import platform
import re
import subprocess
import sys
import time


SCHEMA = "paacs-installed-step-plan-v1"
MAX_JSON = 32 * 1024**2
DEFAULT_OPTIONS = {
    "workers": 1, "validation_workers": 1, "assembly": "preserve_parts",
    "resource_policy": "normal", "construction_effort": "legacy", "recovery": "standard",
    "max_nodes": 2_000_000, "max_cells": 5_000_000,
    "max_working_bytes": 2 * 1024**3, "max_output_bytes": 512 * 1024**2,
    "max_audit_bytes": 512 * 1024**2,
}
DEFAULT_ALLOWANCES = {"native_seconds": 300., "hard_seconds": 600., "inspection_seconds": 180.,
    "inspection_bytes": 2 * 1024**3, "inspection_work": 500_000_000,
    "max_process_bytes": 8 * 1024**3, "max_disk_bytes": 4 * 1024**3}


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, value):
    data = (json.dumps(value, indent=2, allow_nan=False) + "\n").encode("utf-8")
    if len(data) > MAX_JSON:
        raise ValueError("JSON output exceeds 32 MiB allowance")
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(data)
    temporary.replace(path)


def load(path):
    path = Path(path)
    if path.stat().st_size > MAX_JSON:
        raise ValueError("JSON input exceeds 32 MiB allowance")
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("Duplicate JSON key")
            result[key] = value
        return result
    def constant(value):
        raise ValueError("Nonfinite JSON constant: " + value)
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=pairs, parse_constant=constant)


def finite_tree(value):
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("Nonfinite plan option")
    if isinstance(value, dict):
        for key, child in value.items():
            if not isinstance(key, str):
                raise ValueError("Plan object keys must be strings")
            finite_tree(child)
    elif isinstance(value, list):
        for child in value:
            finite_tree(child)
    elif value is not None and not isinstance(value, (str, bool, int, float)):
        raise ValueError("Plan must contain JSON data only")


def normalized_version(value):
    match = re.fullmatch(r"(\d+)\.(\d+)\.(\d+)(?:-?rc\.?([0-9]+))?", value)
    if not match:
        raise ValueError("Expected a stable or RC PAACS version")
    return ".".join(match.group(i) for i in (1, 2, 3)) + ("rc" + match.group(4) if match.group(4) else "")


def inside(root, relative):
    if not isinstance(relative, str) or not relative or "\0" in relative:
        raise ValueError("Expected a relative CAD path")
    windows, posix = PureWindowsPath(relative), PurePosixPath(relative)
    if windows.drive or windows.root or posix.is_absolute() or ".." in windows.parts or ".." in posix.parts:
        raise ValueError("CAD path escapes its root")
    if any(":" in item for item in windows.parts):
        raise ValueError("CAD path contains a drive or alternate stream")
    target = (root / relative).resolve(strict=True)
    if not target.is_relative_to(root) or not target.is_file():
        raise ValueError("CAD path escapes its root or is not a regular file")
    return target


def read_plan(path, root):
    plan = load(path)
    finite_tree(plan)
    if plan.get("schema") != SCHEMA:
        raise ValueError("Unknown installed qualification plan schema")
    version = normalized_version(plan["package_version"])
    files = plan.get("files")
    if not isinstance(files, list) or not 1 <= len(files) <= 10_000:
        raise ValueError("Plan requires 1..10000 explicit files")
    ids, jobs = set(), []
    for row in files:
        identifier = row["id"]
        if not isinstance(identifier, str) or not re.fullmatch(r"[A-Za-z0-9_-]{1,96}", identifier) or identifier in ids:
            raise ValueError("Case IDs must be unique safe directory names")
        ids.add(identifier)
        source = inside(root, row["path"])
        digest = row["sha256"]
        if not isinstance(digest, str) or not re.fullmatch("[0-9a-f]{64}", digest) or sha(source) != digest:
            raise ValueError("Pinned CAD hash differs: " + identifier)
        if ("size_mm" in row) == ("size_factor" in row):
            raise ValueError("Specify exactly one of size_mm or size_factor")
        size = row.get("size_mm", row.get("size_factor"))
        if isinstance(size, bool) or not isinstance(size, (int, float)) or not math.isfinite(size) or size <= 0:
            raise ValueError("Sizing must be finite and positive")
        options = {**DEFAULT_OPTIONS, **row.get("options", {})}
        if any(type(options[name]) is not int or options[name] != 1 for name in ("workers", "validation_workers")):
            raise ValueError("Qualification requires one meshing and one validation worker")
        if any(name in options for name in ("callback", "on_progress", "cancel", "source")):
            raise ValueError("Plan must not inject executable options")
        allowances = {**DEFAULT_ALLOWANCES, **row.get("allowances", {})}
        if set(allowances) != set(DEFAULT_ALLOWANCES):
            raise ValueError("Unknown allowance")
        for name, value in allowances.items():
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError("Allowances must be finite and positive")
            if name.endswith(("_bytes", "_work")) and not isinstance(value, int):
                raise ValueError("Count and byte allowances must be integers")
        if allowances["hard_seconds"] <= allowances["native_seconds"]:
            raise ValueError("Hard deadline must exceed native deadline")
        if "max_seconds" in options and options["max_seconds"] != allowances["native_seconds"]:
            raise ValueError("Use native_seconds to specify the native deadline")
        options["max_seconds"] = allowances["native_seconds"]
        for name in ("max_nodes", "max_cells", "max_working_bytes", "max_output_bytes", "max_audit_bytes"):
            value = options[name]
            if isinstance(value, bool) or not isinstance(value, int) or not 0 < value < 2**64:
                raise ValueError(name + " must be a positive uint64")
        expected_solids = row.get("expected_original_solids")
        if expected_solids is not None and (type(expected_solids) is not int or expected_solids < 1):
            raise ValueError("expected_original_solids must be a positive integer")
        profile = row.get("profile", "unspecified")
        if not isinstance(profile, str) or len(profile) > 128:
            raise ValueError("profile must be a bounded string")
        jobs.append({"id": identifier, "profile": profile, "relative_path": row["path"], "source": str(source), "source_sha256": digest,
            "size_mm": row.get("size_mm"), "size_factor": row.get("size_factor"),
            "options": options, "allowances": allowances, "package_version": version,
            "expected_original_solids": expected_solids})
    return jobs


def clean_environment(python):
    env = dict(os.environ)
    for name in list(env):
        upper = name.upper()
        if upper.startswith(("PYTHON", "CSF_", "CASROOT", "OCCT", "OPENCASCADE", "PAACS_", "DYLD_", "LD_")):
            del env[name]
    paths = [str(Path(python).parent)]
    if os.name == "nt":
        system = Path(env.get("SystemRoot", r"C:\Windows"))
        paths += [str(system / "System32"), str(system)]
    else:
        paths += ["/usr/bin", "/bin"]
    env["PATH"] = os.pathsep.join(paths)
    env["PYTHONNOUSERSITE"] = "1"
    return env


def package_identity(pm, expected):
    distribution = importlib.metadata.distribution("paacs-mesh")
    version = normalized_version(distribution.version)
    if version != expected or normalized_version(pm.__version__) != expected or normalized_version(pm.native_version) != expected:
        raise ValueError("Installed package version differs from plan")
    base = Path(distribution.locate_file("")).resolve()
    package = Path(pm.__file__).resolve()
    native = Path(pm._native.__file__).resolve()
    if not package.is_relative_to(base) or not native.is_relative_to(package.parent):
        raise ValueError("Package/native paths differ from the installed distribution")
    if not base.is_relative_to(Path(sys.prefix).resolve()):
        raise ValueError("Distribution is outside the selected isolated Python prefix")
    wheel_metadata = distribution.read_text("WHEEL")
    if not wheel_metadata or "Wheel-Version:" not in wheel_metadata:
        raise ValueError("Installed distribution has no wheel metadata")
    files = distribution.files or []
    if len(files) > 10_000:
        raise ValueError("Installed wheel RECORD exceeds 10000 entries")
    recorded = {str(item).replace("\\", "/"): item for item in files}
    for path in (package, native):
        entry = recorded.get(path.relative_to(base).as_posix())
        if entry is None or entry.hash is None or entry.hash.mode != "sha256":
            raise ValueError("Installed package/native is absent from hashed wheel RECORD")
        actual = base64.urlsafe_b64encode(bytes.fromhex(sha(path))).decode().rstrip("=")
        if entry.hash.value != actual:
            raise ValueError("Installed package/native differs from wheel RECORD")
    checked, checked_bytes, contents = 0, 0, hashlib.sha256()
    for name, entry in sorted(recorded.items()):
        if entry.hash is None:
            continue  # RECORD itself and generated caches are normally unhashed.
        target = (base / name).resolve(strict=True)
        if not target.is_relative_to(Path(sys.prefix).resolve()) or not target.is_file() or entry.hash.mode != "sha256":
            raise ValueError("Wheel RECORD contains an unowned file or unsupported hash")
        checked_bytes += target.stat().st_size
        if checked_bytes > 512 * 1024**2:
            raise ValueError("Installed wheel hash checks exceed 512 MiB")
        digest = sha(target)
        if entry.hash.value != base64.urlsafe_b64encode(bytes.fromhex(digest)).decode().rstrip("="):
            raise ValueError("Installed runtime differs from wheel RECORD: " + name)
        contents.update((name + "\0" + digest + "\n").encode())
        checked += 1
    direct = distribution.read_text("direct_url.json")
    direct = json.loads(direct) if direct else None
    if direct and (direct.get("dir_info") or direct.get("vcs_info")):
        raise ValueError("Editable/source/VCS installations are not wheel qualification")
    if direct and not str(direct.get("url", "")).split("?", 1)[0].endswith(".whl"):
        raise ValueError("Direct installation must refer to a wheel")
    return {"distribution_version": distribution.version, "package_version": pm.__version__,
        "engine_version": pm.native_version, "package_path": str(package), "native_path": str(native),
        "native_sha256": sha(native), "distribution_metadata_sha256": hashlib.sha256(distribution.read_text("METADATA").encode()).hexdigest(),
        "direct_url": direct, "wheel_record_sha256": hashlib.sha256((distribution.read_text("RECORD") or "").encode()).hexdigest(),
        "verified_record_files": checked, "verified_record_bytes": checked_bytes, "runtime_contents_sha256": contents.hexdigest(),
        "installer": distribution.read_text("INSTALLER"), "python": sys.version, "prefix": sys.prefix,
        "platform": platform.platform()}


def resolve_size(model, job):
    if job["size_mm"] is not None:
        return job["size_mm"], None
    lower, upper, tags, offset = [math.inf]*3, [-math.inf]*3, set(), 0
    while True:
        page = model.inspect("bodies", offset=offset, limit=100)
        if page.get("offset") != offset or page.get("total") != len(model.body_tags):
            raise ValueError("Incomplete body-bound pagination")
        entities = page["entities"]
        if not entities or len(tags) + len(entities) > len(model.body_tags):
            raise ValueError("Incomplete or duplicate body-bound page")
        for entity in entities:
            tag = entity["tag"]
            if tag in tags:
                raise ValueError("Duplicate body-bound record")
            tags.add(tag)
            box = entity["bounds"]
            lo, hi = box["lower"], box["upper"]
            if len(lo) != 3 or len(hi) != 3 or any(not math.isfinite(v) for v in lo + hi) or any(lo[i] > hi[i] for i in range(3)):
                raise ValueError("Invalid body bounds")
            lower = [min(lower[i], lo[i]) for i in range(3)]
            upper = [max(upper[i], hi[i]) for i in range(3)]
        end = page.get("next_offset")
        if end is None:
            break
        if type(end) is not int or end != offset + len(entities):
            raise ValueError("Non-progressing body-bound pagination")
        offset = end
    if tags != set(model.body_tags):
        raise ValueError("Body bounds did not cover original CAD")
    diagonal = math.dist(lower, upper)
    value = diagonal * job["size_factor"]
    if not math.isfinite(value) or value <= 0:
        raise ValueError("Invalid union-bounds size")
    return value, {"rule": "factor_times_Euclidean_diagonal_of_union_body_bounds", "factor": job["size_factor"],
        "lower": lower, "upper": upper, "diagonal_mm": diagonal, "body_count": len(tags)}


def validate_native(report, digest, body_tags, face_tags, mesh):
    if (not isinstance(report, dict) or report.get("source_sha256") != digest or report.get("mesh_published") is not True
            or report.get("status") != "succeeded" or report.get("reason") != "none"):
        raise ValueError("Native publication/source identity gate failed")
    for section in ("audit", "quality"):
        if not all(report.get(section, {}).get(name) is True for name in ("complete", "valid")):
            raise ValueError("Native " + section + " gate failed")
    if report["quality"].get("inverted_cells") != 0:
        raise ValueError("Native inversion gate failed")
    for name in ("missing_boundary", "unexpected_boundary", "nonmanifold_faces", "inconsistent_orientation", "invalid_connectivity", "nonfinite_nodes"):
        if report["audit"].get(name) != 0:
            raise ValueError("Native audit defect: " + name)
    if set(mesh.body_tags) != set(body_tags) or set(mesh.face_tags) != set(face_tags):
        raise ValueError("Published mesh lost CAD body/face tags")


def validate_summary(summary, archive, body_tags, face_tags):
    if summary.get("mesh_sha256") != sha(archive):
        raise ValueError("Exported archive differs from independently inspected snapshot")
    if any(summary.get("orientation_counts", {}).get(name) != 0 for name in ("inverted", "degenerate")):
        raise ValueError("Independent archive orientation failed")
    if summary.get("unclassified_tet4_count") != 0:
        raise ValueError("Independent archive has unclassified tetrahedra")
    if set(map(int, summary["body_tet4_counts"])) != set(body_tags):
        raise ValueError("Encoded body tags differ")
    if set(map(int, summary["sizing"]["achieved"]["boundary_faces"])) != set(face_tags):
        raise ValueError("Encoded face tags differ")


def refusal_kind(error, report):
    reason = str(getattr(error, "reason", report.get("reason", ""))).replace("-", "_")
    if isinstance(error, (TypeError, ValueError, OverflowError)) or reason.startswith(("invalid_request", "unsupported_control", "control_request")):
        return "request_rejected", "not_admitted"
    return "native_refused", "refused"


def worker(job_path):
    job, began = load(job_path), time.perf_counter()
    output = Path(job["output"])
    record = {"schema": "paacs-installed-step-attempt-v1", "id": job["id"], "pid": os.getpid(),
        "started_utc": now(), "phase": "package_import", "status": "running", "native_status": "not_attempted",
        "source_sha256": job["source_sha256"], "relative_path": job["relative_path"], "profile": job["profile"],
        "size_mm": job["size_mm"], "size_factor": job["size_factor"], "options": job["options"],
        "allowances": job["allowances"], "external_retry_count": 0,
        "cad_validation_basis": "native CAD checks retained; inspection does not re-audit CAD"}
    def checkpoint(phase=None):
        if phase:
            record["phase"] = phase
        record["worker_seconds"] = time.perf_counter() - began
        save(output / "worker.json", record)
    checkpoint()
    try:
        import paacs_mesh as pm
        record["package"] = package_identity(pm, job["package_version"])
        if sha(job["source"]) != job["source_sha256"]:
            raise ValueError("Source changed before import")
        checkpoint("step_import")
        at = time.perf_counter()
        model = pm.read_step(job["source"])
        record["import_seconds"] = time.perf_counter() - at
        bodies, faces = list(model.body_tags), list(model.face_tags)
        record.update(original_body_tags=bodies, original_face_tags=faces, model_info=model.info)
        if model.source_sha256 != job["source_sha256"]:
            raise ValueError("Imported CAD source differs")
        if job["expected_original_solids"] is not None and len(bodies) != job["expected_original_solids"]:
            raise ValueError("Original body count differs from pinned corpus")
        checkpoint("size_resolution")
        at = time.perf_counter()
        size, sizing_basis = resolve_size(model, job)
        record.update(effective_size_mm=size, sizing_basis=sizing_basis, size_resolution_seconds=time.perf_counter()-at)
        checkpoint("meshing")
        at = time.perf_counter()
        try:
            mesh = model.mesh(size, **job["options"])
        except Exception as error:
            record["mesh_call_seconds"] = time.perf_counter() - at
            report = getattr(error, "report", None)
            if not isinstance(report, dict) and not isinstance(error, (TypeError, ValueError, OverflowError)):
                raise
            report = report or {}
            if report:
                save(output / "native-report.json", report)
            status, native_status = refusal_kind(error, report)
            record.update(status=status, native_status=native_status, reason=getattr(error, "reason", report.get("reason")),
                error_type=type(error).__name__, error=str(error)[:4096], diagnostic=pm.diagnose(error) if report else None)
            return 0
        record["mesh_call_seconds"] = time.perf_counter() - at
        report = mesh.report
        save(output / "native-report.json", report)
        record["native_report_sha256"] = sha(output / "native-report.json")
        record["native_status"] = "published_pending_qualification"
        validate_native(report, job["source_sha256"], bodies, faces, mesh)
        record.update(nodes=mesh.node_count, tetrahedra=mesh.tetrahedron_count, native_seconds=report.get("seconds"),
            native_quality=report.get("quality"), body_routes=report.get("effective_policy"), diagnostic=pm.diagnose(report))
        checkpoint("exports")
        at = time.perf_counter()
        archive = output / "mesh.paacsmesh"
        mesh.write_archive(archive)
        mesh.write_vtk(output / "mesh.vtk")
        record["export_seconds"] = time.perf_counter() - at
        checkpoint("independent_inspection")
        at = time.perf_counter()
        limits = job["allowances"]
        with pm.inspect_mesh(mesh, max_nodes=job["options"]["max_nodes"], max_cells=job["options"]["max_cells"],
                max_bytes=limits["inspection_bytes"], max_work=limits["inspection_work"], max_seconds=limits["inspection_seconds"]) as inspected:
            summary = inspected.summary()
            validate_summary(summary, archive, bodies, faces)
            inspected.save(output / "inspection")
        save(output / "quality-summary.json", summary)
        record.update(inspection_seconds=time.perf_counter()-at, archive_sha256=sha(archive), vtk_sha256=sha(output / "mesh.vtk"),
            quality=summary["quality"], achieved_sizing=summary["sizing"], orientation_counts=summary["orientation_counts"],
            metric_conventions=summary["metric_conventions"], bounds=summary["bounds"],
            body_tet4_counts=summary["body_tet4_counts"], native_status="qualified")
        if (sha(job["source"]) != job["source_sha256"]
                or package_identity(pm, job["package_version"])["runtime_contents_sha256"] != record["package"]["runtime_contents_sha256"]):
            raise ValueError("Pinned source/native binary changed during attempt")
        record["status"] = "succeeded"
        record["phase"] = "complete"
        return 0
    except Exception as error:
        report = getattr(error, "report", None)
        if record["phase"] == "step_import" and isinstance(report, dict):
            save(output / "import-report.json", report)
            record.update(status="import_refused", reason=getattr(error, "reason", report.get("reason")),
                error_type=type(error).__name__, error=str(error)[:4096])
            return 0
        record.update(status="qualification_failed", error_type=type(error).__name__, error=str(error)[:4096])
        return 1
    finally:
        record["finished_utc"] = now()
        checkpoint()


class WindowsMemory:
    def __init__(self):
        self.query = None
        if os.name == "nt":
            from ctypes import wintypes
            class Counters(ctypes.Structure):
                _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + [
                    (name, ctypes.c_size_t) for name in ("PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage",
                    "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage", "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage")]
            self.counters = Counters
            self.query = ctypes.WinDLL("psapi", use_last_error=True).GetProcessMemoryInfo
            self.query.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
            self.query.restype = wintypes.BOOL

    def sample(self, process):
        if self.query is None:
            return None
        values = self.counters()
        values.cb = ctypes.sizeof(values)
        return int(values.PeakWorkingSetSize) if self.query(int(process._handle), ctypes.byref(values), values.cb) else None


def final_status(receipt, returncode, stop, output):
    if stop:
        return stop + "_native_outcome_not_recertified"
    if returncode != 0:
        return "qualification_failed" if receipt.get("status") == "qualification_failed" else "process_failed"
    if receipt.get("status") in ("native_refused", "request_rejected", "import_refused"):
        return receipt["status"]
    required = ("mesh.paacsmesh", "mesh.vtk", "native-report.json", "quality-summary.json", "inspection/inspection.json")
    if receipt.get("status") != "succeeded" or receipt.get("native_status") != "qualified" or not all((output / name).is_file() for name in required):
        return "incomplete_result"
    return "succeeded"


def run_attempt(python, job, output):
    output.mkdir()
    scratch = output / "scratch"
    scratch.mkdir()
    job = {**job, "output": str(output)}
    save(output / "job.json", job)
    command = [str(python), "-I", str(Path(__file__).resolve()), "--worker", str(output / "job.json")]
    began, memory, peak, stop = time.perf_counter(), WindowsMemory(), None, None
    receipt = {"schema": "paacs-installed-step-controller-v1", "id": job["id"], "job": job,
        "command": command, "runner_sha256": sha(__file__), "started_utc": now(), "status": "running",
        "process_memory_limit_enforced": os.name == "nt",
        "memory_scope": "modeled native memory is enforced by mesher; sampled peak process working-set watchdog is Windows-only; unavailable process memory is null",
        "vtk_validation_scope": "native export and byte hash; VTK tags are not independently parsed"}
    save(output / "controller.json", receipt)
    environment = clean_environment(python)
    environment.update({key: str(scratch) for key in ("TEMP", "TMP", "TMPDIR")})
    with (output / "worker.log").open("wb") as stream:
        try:
            process = subprocess.Popen(command, cwd=output, env=environment, stdout=stream, stderr=subprocess.STDOUT,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        except OSError as error:
            receipt.update(status="process_start_failed", error=str(error)[:4096], process_seconds=time.perf_counter()-began, finished_utc=now())
            save(output / "controller.json", receipt)
            return receipt
        receipt["pid"] = process.pid
        save(output / "controller.json", receipt)
        try:
            last_disk = began
            while process.poll() is None:
                sampled = memory.sample(process)
                peak = max(peak or 0, sampled) if sampled is not None else peak
                elapsed = time.perf_counter() - began
                if elapsed > job["allowances"]["hard_seconds"]:
                    stop = "hard_timeout"
                elif peak is not None and peak > job["allowances"]["max_process_bytes"]:
                    stop = "process_memory_stop"
                if time.perf_counter() - last_disk >= 2:
                    if disk_bytes(output) > job["allowances"]["max_disk_bytes"]:
                        stop = "output_disk_stop"
                    last_disk = time.perf_counter()
                if stop:
                    process.kill()
                    break
                time.sleep(.1)
            process.wait()
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
    sampled = memory.sample(process)
    peak = max(peak or 0, sampled) if sampled is not None else peak
    final_disk = disk_bytes(output)
    if final_disk > job["allowances"]["max_disk_bytes"]:
        stop = stop or "output_disk_stop"
    if peak is not None and peak > job["allowances"]["max_process_bytes"]:
        stop = stop or "process_memory_stop"
    try:
        observed = load(output / "worker.json") if (output / "worker.json").is_file() else {}
    except (OSError, ValueError) as error:
        observed = {"status": "malformed_worker_receipt", "error": str(error)[:4096]}
    receipt.update(worker=observed, worker_pid_verified=observed.get("pid") == process.pid,
        process_exit_code=process.returncode, process_seconds=time.perf_counter()-began,
        peak_process_working_set_bytes=peak, durable_output_bytes=final_disk, stop=stop, finished_utc=now())
    receipt["process_memory_limit_enforced"] = os.name == "nt" and peak is not None
    receipt["status"] = final_status(observed, process.returncode, stop, output)
    if (observed and not receipt["worker_pid_verified"]) or sha(job["source"]) != job["source_sha256"]:
        receipt["status"] = "input_or_worker_identity_failed"
    save(output / "controller.json", receipt)
    return receipt


def disk_bytes(root):
    total = 0
    for path in root.rglob("*"):
        try:
            if path.is_file():
                total += path.stat().st_size
        except FileNotFoundError:
            # Public inspection bundle publication renames a temporary directory.
            continue
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", type=Path)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        return worker(args.worker)
    if not all((args.python, args.plan, args.root, args.output)):
        parser.error("Supply --python, --plan, --root and --output")
    python, root = args.python.resolve(strict=True), args.root.resolve(strict=True)
    if not python.is_file() or not root.is_dir() or args.output.exists():
        parser.error("Use an existing interpreter/CAD directory and a fresh output directory")
    jobs = read_plan(args.plan, root)
    output = args.output.resolve()
    output.mkdir(parents=True)
    save(output / "plan.json", {"schema": SCHEMA, "plan_sha256": sha(args.plan), "python": str(python), "jobs": jobs,
        "execution": "serial isolated processes; all retries inside native calls are retained in native reports"})
    if args.plan_only:
        print(json.dumps({"planned": len(jobs), "native_calls": 0}))
        return 0
    rows = []
    for job in jobs:
        result = run_attempt(python, job, output / job["id"])
        rows.append(result)
        save(output / "results.json", {"schema": "paacs-installed-step-results-v1", "attempts": rows,
            "attempt_count": len(rows), "complete": len(rows) == len(jobs)})
        print(json.dumps({"id": job["id"], "status": result["status"], "process_seconds": result["process_seconds"]}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
