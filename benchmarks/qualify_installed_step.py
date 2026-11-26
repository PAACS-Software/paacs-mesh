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
import signal
import subprocess
import sys
import time
import uuid


SCHEMA = "paacs-installed-step-plan-v1"
MAX_JSON = 32 * 1024**2
# SYNCHRONIZE, QUERY(_LIMITED)_INFORMATION, VM_READ, SET_QUOTA and TERMINATE.
# SET_QUOTA + TERMINATE are required to assign an existing redirector child to our job.
WINDOWS_WORKER_RIGHTS = 0x00101511
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
    # Callers may supply platform aliases such as macOS /var -> /private/var.
    # Compare canonical root and target paths while retaining the escape guard.
    root = Path(root).resolve(strict=True)
    if not root.is_dir():
        raise ValueError("CAD root must be a directory")
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
    root = Path(root).resolve(strict=True)
    if not root.is_dir():
        raise ValueError("CAD root must be a directory")
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
        "started_utc": now(), "phase": "controller_handshake", "status": "running", "native_status": "not_attempted",
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
        await_controller(job_path, job)
        checkpoint("package_import")
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


def await_controller(job_path, job):
    """Do not import PAACS until our live process identity is acknowledged."""
    output = Path(job["output"])
    identity = {"schema": "paacs-installed-worker-hello-v1", "pid": os.getpid(), "parent_pid": os.getppid(),
        "token": job["handshake_token"], "job_sha256": sha(job_path)}
    save(output / "worker-hello.json", identity)
    deadline = time.monotonic() + min(30., job["allowances"]["hard_seconds"])
    while time.monotonic() < deadline:
        path = output / "controller-ack.json"
        if path.is_file():
            ack = load(path)
            if (ack.get("schema") != "paacs-installed-controller-ack-v1" or ack.get("worker_pid") != os.getpid()
                    or ack.get("controller_pid") != job["controller_pid"] or ack.get("token") != identity["token"]
                    or ack.get("job_sha256") != identity["job_sha256"]):
                raise ValueError("Controller acknowledgement does not match this worker/job")
            return
        time.sleep(.02)
    raise TimeoutError("Controller did not verify/acknowledge this worker before startup deadline")


def ancestry_chain(parents, worker_pid, launcher_pid):
    """Resolve only a finite creator chain ending at the controlled launcher."""
    chain = [worker_pid]
    while chain[-1] != launcher_pid:
        parent = parents.get(chain[-1])
        if type(parent) is not int or parent <= 0 or parent in chain or len(chain) >= 32:
            raise ValueError("Worker PID is not a verified descendant of the controlled launcher")
        chain.append(parent)
    return chain


def validate_creation_chain(chain, lifetimes):
    for child, parent in zip(chain, chain[1:]):
        created = lifetimes[child][0]
        start, end = lifetimes[parent]
        if created < start or (end and created > end):
            raise ValueError("Process creator chain contains a reused PID or inconsistent lifetime")


class WindowsAPI:
    """Narrow process/job API; process handles remain owned until cleanup."""
    def __init__(self):
        from ctypes import wintypes as w
        self.w = w
        self.kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        declarations = {
            "OpenProcess": ([w.DWORD, w.BOOL, w.DWORD], w.HANDLE),
            "CloseHandle": ([w.HANDLE], w.BOOL),
            "GetProcessTimes": ([w.HANDLE, ctypes.POINTER(w.FILETIME), ctypes.POINTER(w.FILETIME), ctypes.POINTER(w.FILETIME), ctypes.POINTER(w.FILETIME)], w.BOOL),
            "WaitForSingleObject": ([w.HANDLE, w.DWORD], w.DWORD),
            "TerminateProcess": ([w.HANDLE, w.UINT], w.BOOL),
            "GetExitCodeProcess": ([w.HANDLE, ctypes.POINTER(w.DWORD)], w.BOOL),
            "CreateToolhelp32Snapshot": ([w.DWORD, w.DWORD], w.HANDLE),
            "CreateJobObjectW": ([ctypes.c_void_p, w.LPCWSTR], w.HANDLE),
            "SetInformationJobObject": ([w.HANDLE, ctypes.c_int, ctypes.c_void_p, w.DWORD], w.BOOL),
            "AssignProcessToJobObject": ([w.HANDLE, w.HANDLE], w.BOOL),
            "IsProcessInJob": ([w.HANDLE, w.HANDLE, ctypes.POINTER(w.BOOL)], w.BOOL),
            "TerminateJobObject": ([w.HANDLE, w.UINT], w.BOOL),
        }
        for name, (args, result) in declarations.items():
            function = getattr(self.kernel, name)
            function.argtypes, function.restype = args, result

    def checked(self, value):
        if not value:
            raise ctypes.WinError(ctypes.get_last_error())
        return value

    def lifetime(self, handle):
        values = [self.w.FILETIME() for _ in range(4)]
        self.checked(self.kernel.GetProcessTimes(handle, *(ctypes.byref(value) for value in values)))
        return tuple((value.dwHighDateTime << 32) | value.dwLowDateTime for value in values[:2])

    def parents(self):
        w = self.w
        class Entry(ctypes.Structure):
            _fields_ = [("dwSize", w.DWORD), ("cntUsage", w.DWORD), ("pid", w.DWORD),
                ("heap", ctypes.c_size_t), ("module", w.DWORD), ("threads", w.DWORD),
                ("parent", w.DWORD), ("priority", w.LONG), ("flags", w.DWORD), ("name", w.WCHAR * 260)]
        first, next_entry = self.kernel.Process32FirstW, self.kernel.Process32NextW
        for function in (first, next_entry):
            function.argtypes, function.restype = [w.HANDLE, ctypes.POINTER(Entry)], w.BOOL
        handle = self.kernel.CreateToolhelp32Snapshot(2, 0)
        if handle == ctypes.c_void_p(-1).value:
            raise ctypes.WinError(ctypes.get_last_error())
        try:
            entry, result = Entry(), {}
            entry.dwSize = ctypes.sizeof(entry)
            present = first(handle, ctypes.byref(entry))
            while present:
                result[int(entry.pid)] = int(entry.parent)
                present = next_entry(handle, ctypes.byref(entry))
            return result
        finally:
            self.kernel.CloseHandle(handle)


class WindowsProcess:
    def __init__(self, api, pid, handle=None):
        self.api, self.pid, self.owned = api, pid, handle is None
        self._handle = api.checked(api.kernel.OpenProcess(WINDOWS_WORKER_RIGHTS, False, pid)) if handle is None else int(handle)
        self.lifetime = api.lifetime(self._handle)

    def alive(self):
        result = self.api.kernel.WaitForSingleObject(self._handle, 0)
        if result not in (0, 258):
            raise ctypes.WinError(ctypes.get_last_error())
        return result == 258

    def terminate(self):
        if self.alive():
            if not self.api.kernel.TerminateProcess(self._handle, 1):
                error = ctypes.get_last_error()
                # Taskkill/job termination can already be completing on another
                # thread. Windows then reports ACCESS_DENIED before the process
                # handle becomes signalled; accept only observed termination.
                if self.api.kernel.WaitForSingleObject(self._handle, 1000) != 0:
                    raise ctypes.WinError(error)

    def wait(self, milliseconds=5000):
        if self.api.kernel.WaitForSingleObject(self._handle, milliseconds) != 0:
            raise TimeoutError("Verified worker did not terminate within cleanup deadline")

    def close(self):
        if self.owned and self._handle:
            self.api.kernel.CloseHandle(self._handle)
            self._handle = None

    def exit_code(self):
        if self.alive():
            return None
        code = self.api.w.DWORD()
        self.api.checked(self.api.kernel.GetExitCodeProcess(self._handle, ctypes.byref(code)))
        return int(code.value)


class WindowsJob:
    def __init__(self, api):
        self.api, self._handle = api, api.checked(api.kernel.CreateJobObjectW(None, None))
        w = api.w
        class Basic(ctypes.Structure):
            _fields_ = [("process_time", ctypes.c_int64), ("job_time", ctypes.c_int64), ("flags", w.DWORD),
                ("minimum", ctypes.c_size_t), ("maximum", ctypes.c_size_t), ("active", w.DWORD),
                ("affinity", ctypes.c_size_t), ("priority", w.DWORD), ("scheduling", w.DWORD)]
        class IO(ctypes.Structure):
            _fields_ = [(name, ctypes.c_uint64) for name in ("read_ops", "write_ops", "other_ops", "read", "write", "other")]
        class Extended(ctypes.Structure):
            _fields_ = [("basic", Basic), ("io", IO), ("process_memory", ctypes.c_size_t), ("job_memory", ctypes.c_size_t),
                ("peak_process", ctypes.c_size_t), ("peak_job", ctypes.c_size_t)]
        info = Extended()
        info.basic.flags = 0x2000  # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE, non-inherited controller handle.
        try:
            api.checked(api.kernel.SetInformationJobObject(self._handle, 9, ctypes.byref(info), ctypes.sizeof(info)))
        except BaseException:
            self.close()
            raise

    def assign(self, process):
        present = self.api.w.BOOL()
        self.api.checked(self.api.kernel.IsProcessInJob(process._handle, self._handle, ctypes.byref(present)))
        if not present.value:
            self.api.checked(self.api.kernel.AssignProcessToJobObject(self._handle, process._handle))

    def terminate(self):
        self.api.checked(self.api.kernel.TerminateJobObject(self._handle, 1))

    def close(self):
        if self._handle:
            self.api.kernel.CloseHandle(self._handle)
            self._handle = None


class ProcessControl:
    def __init__(self, launcher):
        self.launcher, self.worker, self.chain, self.job = launcher, None, None, None
        self.cleanup_members = []
        self.api = WindowsAPI() if os.name == "nt" else None
        self.root = WindowsProcess(self.api, launcher.pid, launcher._handle) if self.api else None
        if self.api:
            self.job = WindowsJob(self.api)

    def contain_launcher(self):
        if self.job:
            self.job.assign(self.root)

    def establish(self, pid):
        if type(pid) is not int or pid <= 0:
            raise ValueError("Worker handshake has no valid PID")
        if not self.api:
            if pid != self.launcher.pid:
                raise ValueError("Non-Windows worker PID differs from controlled launcher")
            self.worker, self.chain = self.launcher, [pid]
            return {"worker_pid": pid, "launcher_pid": pid, "chain": [pid], "kind": "direct_pid"}
        candidate, intermediates = WindowsProcess(self.api, pid), []
        try:
            if not candidate.alive():
                raise ValueError("Worker terminated before identity verification")
            chain = ancestry_chain(self.api.parents(), pid, self.launcher.pid)
            # Refresh the retained root handle: a cached live (exit=0) stamp
            # must not authenticate a child created after a reused launcher PID.
            self.root.lifetime = self.api.lifetime(self.root._handle)
            lifetimes = {pid: candidate.lifetime, self.launcher.pid: self.root.lifetime}
            for parent in chain[1:-1]:
                reference = WindowsProcess(self.api, parent)
                intermediates.append(reference)
                lifetimes[parent] = reference.lifetime
            validate_creation_chain(chain, lifetimes)
            if ancestry_chain(self.api.parents(), pid, self.launcher.pid) != chain:
                raise ValueError("Worker ancestry changed during verification")
            self.worker, self.chain = candidate, chain
            self.job.assign(candidate)  # Required before acknowledgement; no uncontained Native work.
            return {"worker_pid": pid, "launcher_pid": self.launcher.pid, "chain": chain,
                "creation_and_exit_filetimes": {str(key): list(value) for key, value in lifetimes.items()},
                "kind": "windows_verified_process_handles", "kill_on_close_job": True}
        except BaseException:
            if self.worker is not candidate:
                candidate.close()
            raise
        finally:
            for reference in intermediates:
                reference.close()

    def alive(self):
        return self.worker.alive() if self.api and self.worker else self.launcher.poll() is None

    def exit_code(self):
        return self.worker.exit_code() if self.api and self.worker else self.launcher.poll()

    def contain_pending_descendants(self):
        """Capture pre-job redirector children even when no hello ever arrived.

        Parents are contained before children, so later children inherit our
        job. Extra bounded snapshots cover children born during that capture.
        Every accepted handle passes the same creator/lifetime checks as the
        acknowledged worker; unrelated or reused PIDs are never assigned/killed.
        """
        failures, deadline = [], time.monotonic() + 5
        for _ in range(3):
            parents = self.api.parents()
            scoped = []
            for pid in parents:
                if pid == self.launcher.pid:
                    continue
                try:
                    chain = ancestry_chain(parents, pid, self.launcher.pid)
                except ValueError:
                    continue
                scoped.append((len(chain), pid))
            for _, pid in sorted(scoped):
                if len(self.cleanup_members) >= 384 or time.monotonic() >= deadline:
                    return failures + ["Bounded descendant capture exhausted; controlled tree/job fallback remains active"]
                previous, previous_chain = self.worker, self.chain
                try:
                    self.establish(pid)
                except (OSError, ValueError) as error:
                    failures.append("Descendant capture refused: " + str(error)[:1024])
                finally:
                    if self.worker is not None and self.worker is not previous:
                        self.cleanup_members.append(self.worker)
                    self.worker, self.chain = previous, previous_chain
        return failures

    def cleanup(self, stop=False):
        """Terminate only the owned job/group and opened, verified worker handle."""
        errors = []
        if stop:
            if self.api:
                try:
                    errors.extend(self.contain_pending_descendants())
                except OSError as error:
                    errors.append(str(error))
            if self.api and self.launcher.poll() is None:
                # Controlled Popen root only. The job and retained handle avoid orphaning a redirector child.
                try:
                    subprocess.run([str(Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32/taskkill.exe"),
                        "/PID", str(self.launcher.pid), "/T", "/F"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=10, check=False)
                except (OSError, subprocess.TimeoutExpired) as error:
                    errors.append(str(error))
            elif not self.api and self.launcher.poll() is None:
                try:
                    os.killpg(self.launcher.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            if self.job:
                try:
                    self.job.terminate()
                except OSError as error:
                    errors.append(str(error))
            for reference in self.cleanup_members:
                reference.terminate()
                reference.wait()
            if self.api and self.worker:
                self.worker.terminate()
                self.worker.wait()
            if self.launcher.poll() is None:
                self.launcher.kill()
        if self.launcher.poll() is None:
            self.launcher.wait(timeout=5)
        return errors

    def close(self):
        if self.job:
            self.job.close()  # Also protects against interrupted controller cleanup.
        if self.api and self.worker:
            self.worker.close()
        for reference in self.cleanup_members:
            reference.close()


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


def validate_hello(hello, job, job_digest):
    if (hello.get("schema") != "paacs-installed-worker-hello-v1" or hello.get("token") != job["handshake_token"]
            or hello.get("job_sha256") != job_digest or type(hello.get("pid")) is not int or hello["pid"] <= 0):
        raise ValueError("Worker handshake token/job/PID differs from the controlled request")


def terminate_launcher(launcher):
    """Startup fallback targets only the Popen-owned launcher and its tree."""
    if launcher.poll() is None:
        if os.name == "nt":
            try:
                subprocess.run([str(Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32/taskkill.exe"),
                    "/PID", str(launcher.pid), "/T", "/F"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=10, check=False)
            finally:
                if launcher.poll() is None:
                    launcher.kill()
        else:
            try:
                os.killpg(launcher.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    launcher.wait(timeout=5)


def run_attempt(python, job, output):
    output.mkdir()
    scratch = output / "scratch"
    scratch.mkdir()
    job = {**job, "output": str(output), "handshake_token": uuid.uuid4().hex, "controller_pid": os.getpid()}
    job_path = output / "job.json"
    save(job_path, job)
    job_digest = sha(job_path)
    command = [str(python), "-I", str(Path(__file__).resolve()), "--worker", str(job_path)]
    began, memory, peak, stop = time.perf_counter(), WindowsMemory(), None, None
    control, process, acknowledged, pending_exception = None, None, False, None
    worker_exit, cleanup_errors = None, []
    receipt = {"schema": "paacs-installed-step-controller-v1", "id": job["id"], "job": job,
        "command": command, "runner_sha256": sha(__file__), "started_utc": now(), "status": "running",
        "process_memory_limit_enforced": False,
        "memory_scope": "modeled native memory is enforced by mesher; sampled actual-worker peak working-set watchdog is Windows-only; unavailable process memory is null",
        "vtk_validation_scope": "native export and byte hash; VTK tags are not independently parsed"}
    save(output / "controller.json", receipt)
    environment = clean_environment(python)
    environment.update({key: str(scratch) for key in ("TEMP", "TMP", "TMPDIR")})
    with (output / "worker.log").open("wb") as stream:
        try:
            process = subprocess.Popen(command, cwd=output, env=environment, stdout=stream, stderr=subprocess.STDOUT,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0), start_new_session=os.name != "nt")
        except OSError as error:
            receipt.update(status="process_start_failed", error=str(error)[:4096], process_seconds=time.perf_counter()-began, finished_utc=now())
            save(output / "controller.json", receipt)
            return receipt
        receipt.update(pid=process.pid, launcher_pid=process.pid)
        save(output / "controller.json", receipt)
        try:
            control = ProcessControl(process)
            control.contain_launcher()
            last_disk = began
            while True:
                hello_path = output / "worker-hello.json"
                if not acknowledged and hello_path.is_file():
                    hello = load(hello_path)
                    validate_hello(hello, job, job_digest)
                    identity = control.establish(hello["pid"])
                    peak = None  # The benchmark RSS belongs to the verified worker, not a venv redirector.
                    expected_parent = identity["chain"][1] if len(identity["chain"]) > 1 else job["controller_pid"]
                    if hello.get("parent_pid") != expected_parent:
                        raise ValueError("Worker handshake creator differs from verified process chain")
                    receipt.update(worker_identity=identity, worker_pid=hello["pid"], memory_measurement_pid=hello["pid"])
                    save(output / "controller.json", receipt)
                    if process.poll() is not None:
                        stop = "launcher_exited_before_acknowledgement"
                        break
                    save(output / "controller-ack.json", {"schema": "paacs-installed-controller-ack-v1",
                        "controller_pid": job["controller_pid"], "worker_pid": hello["pid"],
                        "token": job["handshake_token"], "job_sha256": job_digest})
                    acknowledged = True
                reference = control.worker if control.worker is not None else process
                sampled = memory.sample(reference)
                peak = max(peak or 0, sampled) if sampled is not None else peak
                elapsed = time.perf_counter()-began
                if acknowledged and not control.alive():
                    worker_exit = control.exit_code()
                    process.wait(timeout=5)
                    break
                if process.poll() is not None:
                    stop = "launcher_exited_before_worker" if acknowledged else "launcher_exited_before_handshake"
                elif not acknowledged and elapsed > min(30., job["allowances"]["hard_seconds"]):
                    stop = "worker_handshake_timeout"
                elif elapsed > job["allowances"]["hard_seconds"]:
                    stop = "hard_timeout"
                elif peak is not None and peak > job["allowances"]["max_process_bytes"]:
                    stop = "process_memory_stop"
                if time.perf_counter()-last_disk >= 2:
                    if disk_bytes(output) > job["allowances"]["max_disk_bytes"]:
                        stop = "output_disk_stop"
                    last_disk = time.perf_counter()
                if stop:
                    break
                time.sleep(.02 if not acknowledged else .1)
        except BaseException as error:
            stop = "controller_exception"
            receipt.update(controller_error_type=type(error).__name__, controller_error=str(error)[:4096])
            if isinstance(error, (KeyboardInterrupt, SystemExit)):
                pending_exception = error
        finally:
            try:
                if control is not None:
                    cleanup_errors.extend(control.cleanup(stop=stop is not None))
                    if control.worker is not None:
                        sampled = memory.sample(control.worker)
                        peak = max(peak or 0, sampled) if sampled is not None else peak
                        worker_exit = control.exit_code()
                else:
                    terminate_launcher(process)
            except BaseException as error:
                cleanup_errors.append(type(error).__name__ + ": " + str(error)[:4096])
                stop = stop or "process_cleanup_failed"
            finally:
                if control is not None:
                    control.close()  # Windows KILL_ON_JOB_CLOSE contains abrupt controller cancellation too.
    final_disk = disk_bytes(output)
    final_elapsed = time.perf_counter()-began
    if final_elapsed > job["allowances"]["hard_seconds"]:
        stop = stop or "hard_timeout"
    if final_disk > job["allowances"]["max_disk_bytes"]:
        stop = stop or "output_disk_stop"
    if peak is not None and peak > job["allowances"]["max_process_bytes"]:
        stop = stop or "process_memory_stop"
    try:
        observed = load(output / "worker.json") if (output / "worker.json").is_file() else {}
    except (OSError, ValueError) as error:
        observed = {"status": "malformed_worker_receipt", "error": str(error)[:4096]}
    receipt.update(worker=observed, worker_pid_verified=acknowledged and observed.get("pid") == receipt.get("worker_pid"),
        worker_acknowledged=acknowledged, launcher_exit_code=process.poll(), process_exit_code=worker_exit,
        process_seconds=final_elapsed, peak_process_working_set_bytes=peak,
        durable_output_bytes=final_disk, stop=stop, cleanup_errors=cleanup_errors, finished_utc=now())
    receipt["process_memory_limit_enforced"] = os.name == "nt" and peak is not None
    receipt["status"] = final_status(observed, worker_exit, stop, output)
    try:
        source_unchanged = sha(job["source"]) == job["source_sha256"]
    except OSError:
        source_unchanged = False
    if (observed and not receipt["worker_pid_verified"]) or not source_unchanged:
        receipt["status"] = "input_or_worker_identity_failed"
    save(output / "controller.json", receipt)
    if pending_exception is not None:
        raise pending_exception
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


def python_invocation(value):
    """Validate the executable target without changing the venv invocation path.

    On POSIX, venv/bin/python commonly symlinks to a base interpreter. Executing
    the resolved base path loses pyvenv.cfg discovery and selects another prefix.
    """
    invocation = Path(value).absolute()
    if not invocation.resolve(strict=True).is_file():
        raise ValueError("Selected Python target must be a regular executable file")
    return invocation


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
    python, root = python_invocation(args.python), args.root.resolve(strict=True)
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
