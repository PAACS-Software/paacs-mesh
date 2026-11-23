#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Run a bounded, reproducible CAD sample in isolated processes.

Example from the repository root with the corpus unpacked beside the repository:
  python scripts/run_cad_corpus.py build/windows-occt/Release/paacs-cad-mesh.exe \
    ../abc_0000_step_v00 --output evidence/local-cad \
    --limit 30 --smallest 6 --spread 12 --repeats 3 --timeout 60

File/work/output budgets are mesher limits. The timeout kills the direct child;
peak process memory is measured, not limited. Refusal is reported separately
from a successful mesh and does not qualify unsupported geometry.
Use --native-cad for the native BRep adapter; recovery extensions are explicit.
The historical runner work cap remains 2,000,000. Use --resource-policy normal
with --max-work none to inherit the selected CLI's Normal work policy.
Run --self-test alone to check both report schemas without a CAD SDK.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import ctypes
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time
import uuid

SCHEMA = "paacs-cad-corpus-v1"
MAGIC = b"PAACSM8\0"
REPORT_SCHEMAS = {"paacs-cad-report-v1", "paacs-native-cad-report-v1"}
MODE_FLAGS = ("native_cad", "planar_recovery", "general_recovery", "separate_contact_bodies", "separate_matched_box_faces", "reuse_native_queries", "defer_discarded_native_facet_checks", "separate_matched_native_faces", "retain_auxiliary_geometry", "native_kernel_admission", "longitudinal_recovery", "independent_parts", "allow_centroids_within_chord_tolerance")
OPTIONAL_NATIVE_FLAGS = ("local_partition_recovery", "adaptive_tet_recovery", "trace_native_stages")
STAGES = ("import_seconds", "preparation_seconds", "generation_seconds",
          "validation_seconds", "audit_seconds", "quality_seconds",
          "archive_seconds", "total_seconds")


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def suppress_native_error_dialogs():
    if os.name != "nt":
        return
    from ctypes import wintypes
    # Children inherit this mode. Missing DLL/entry-point and native fault
    # dialogs must not turn a loader failure into an apparent meshing timeout;
    # the child still returns its native failure status for classification.
    set_mode = ctypes.WinDLL("kernel32", use_last_error=True).SetErrorMode
    set_mode.argtypes = [wintypes.UINT]
    set_mode.restype = wintypes.UINT
    set_mode(0x0001 | 0x0002 | 0x8000)


def select_files(root, limit, smallest, spread, includes):
    files = sorted((p for p in root.rglob("*") if p.is_file() and
                    p.suffix.lower() in {".step", ".stp", ".brep"}),
                   key=lambda p: p.relative_to(root).as_posix())
    chosen = set(files[:limit])
    chosen.update(sorted(files, key=lambda p: (p.stat().st_size,
                                             p.relative_to(root).as_posix()))[:smallest])
    if files and spread:
        for i in range(spread):
            chosen.add(files[i * (len(files) - 1) // max(1, spread - 1)])
    for value in includes:
        path = (root / value).resolve()
        if not path.is_relative_to(root) or path not in files:
            raise ValueError(f"included CAD file is not in the corpus: {value}")
        chosen.add(path)
    rows = []
    for path in sorted(chosen, key=lambda p: p.relative_to(root).as_posix()):
        rows.append({"path": path.relative_to(root).as_posix(),
                     "bytes": path.stat().st_size, "sha256": digest(path)})
    return len(files), rows


class WindowsMemory:
    """Read the direct child's peak working set; no process-memory cap."""
    def __init__(self):
        self.query = None
        if os.name != "nt":
            return
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + [
                (name, ctypes.c_size_t) for name in
                ("PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage",
                 "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage",
                 "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage")]
        self.counters = Counters
        self.query = ctypes.WinDLL("psapi", use_last_error=True).GetProcessMemoryInfo
        self.query.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        self.query.restype = wintypes.BOOL

    def peak(self, process):
        if self.query is None:
            return None
        counters = self.counters()
        counters.cb = ctypes.sizeof(counters)
        if self.query(int(process._handle), ctypes.byref(counters), counters.cb):
            return int(counters.PeakWorkingSetSize)
        return None


def log_tail(path, maximum=1024 * 1024):
    with path.open("rb") as stream:
        stream.seek(0, 2)
        size = stream.tell()
        stream.seek(max(0, size - maximum))
        return stream.read().decode("utf-8", errors="replace"), size > maximum


def last_report(text):
    def invalid_constant(value):
        raise ValueError(f"nonfinite JSON number: {value}")
    for line in reversed(text.splitlines()):
        try:
            data = json.loads(line, parse_constant=invalid_constant)
            if isinstance(data, dict) and data.get("schema") in REPORT_SCHEMAS:
                return data
        except (ValueError, TypeError):
            pass
    return None


def classify(code, timed_out, payload, output_exists, valid_archive):
    if timed_out:
        return "timeout"
    if code not in range(7):
        return "crash"
    if code == 0:
        success = "succeeded" if payload and payload.get("schema") == \
            "paacs-native-cad-report-v1" else "complete"
        if payload and payload.get("status") == success and valid_archive and \
                payload.get("archive_published") is True and payload.get("mesh_published") is True:
            return "succeeded"
        return "invalid-result"
    if output_exists:
        return "publication-violation"
    if not payload:
        return "missing-report"
    reason = str(payload.get("reason", "")).replace("_", "-")
    if reason.startswith("unsupported-") or reason in {
            "unsupported-chart", "no-matching-strategy", "unresolved-body-interface"}:
        return "unsupported"
    if "budget" in reason:
        return "resource-limit"
    if reason == "allocation-failed":
        return "allocation-failed"
    if reason in {"size-requirement-not-met", "quality-requirement-not-met"}:
        return "requirement-refusal"
    if reason == "parse-failed":
        return "malformed-input"
    if reason in {"invalid-shape", "polyhedral-validation-failed", "volume-mismatch"}:
        return "invalid-geometry"
    if reason == "cancelled":
        return "cancelled"
    if reason in {"invalid-request", "file-unavailable"}:
        return "rejected"
    return "failed"


def mode_arguments(args):
    return ["--" + key.replace("_", "-") for key in MODE_FLAGS if getattr(args, key, False)]


def native_profile_arguments(args):
    flags = []
    for key in ("recovery_mode", "elements"):
        value = getattr(args, key, None)
        if value is not None:
            flags += ["--" + key.replace("_", "-"), value]
    flags += ["--" + key.replace("_", "-") for key in OPTIONAL_NATIVE_FLAGS if getattr(args, key, False)]
    if flags and not args.native_cad:
        raise ValueError("native recovery profile requires --native-cad")
    return flags


def resource_cap(text):
    """None inherits CLI policy; an explicit cap must fit native uint64."""
    if text.lower() == "none":
        return None
    try:
        value = int(text)
    except ValueError as error:
        raise argparse.ArgumentTypeError("resource cap must be a positive integer or none") from error
    if not 1 <= value <= 2**64 - 1:
        raise argparse.ArgumentTypeError("resource cap must fit a positive native 64-bit integer")
    return value


def work_arguments(args):
    """Keep the historical default flag; omit only an explicitly inherited cap."""
    flags = []
    if args.max_work is not None:
        flags += ["--max-work", str(args.max_work)]
    policy = getattr(args, "resource_policy", None)
    if policy is not None:
        if not args.native_cad:
            raise ValueError("resource policy requires --native-cad")
        flags += ["--resource-policy", policy]
    return flags


def self_test():
    """Independent report fixtures retain legacy classification guarantees."""
    legacy = {"schema": "paacs-cad-report-v1", "status": "complete",
              "mesh_published": True, "archive_published": True}
    native = {**legacy, "schema": "paacs-native-cad-report-v1", "status": "succeeded"}
    for payload in (legacy, native):
        assert last_report("OCCT diagnostic\n" + json.dumps(payload)) == payload
        assert classify(0, False, payload, True, True) == "succeeded"
        assert classify(0, False, payload, False, False) == "invalid-result"
        assert classify(0, False, {**payload, "mesh_published": False}, True, True) == "invalid-result"
    assert classify(0, False, {**native, "status": "complete"}, True, True) == "invalid-result"
    assert classify(0, False, {**legacy, "status": "succeeded"}, True, True) == "invalid-result"
    for reason, expected in (("unsupported_body", "unsupported"),
                             ("unresolved_body_interface", "unsupported"),
                             ("budget_exceeded", "resource-limit"),
                             ("allocation_failed", "allocation-failed"),
                             ("cancelled", "cancelled"),
                             ("invalid_shape", "invalid-geometry")):
        refusal = {**native, "status": "refused", "reason": reason,
                   "mesh_published": False, "archive_published": False}
        assert classify(3, False, refusal, False, False) == expected
        assert classify(3, False, refusal, True, False) == "publication-violation"
    assert last_report('{"schema":"paacs-native-cad-report-v1","n":NaN}') is None
    assert mode_arguments(argparse.Namespace(native_cad=False, planar_recovery=False,
                                            general_recovery=False)) == []
    assert mode_arguments(argparse.Namespace(native_cad=True, planar_recovery=True,
                                            general_recovery=True)) == [
                                                "--native-cad", "--planar-recovery", "--general-recovery"]
    assert work_arguments(argparse.Namespace(max_work=2_000_000, resource_policy=None,
                                              native_cad=False)) == ["--max-work", "2000000"]
    assert work_arguments(argparse.Namespace(max_work=None, resource_policy="normal",
                                              native_cad=True)) == ["--resource-policy", "normal"]
    assert native_profile_arguments(argparse.Namespace(native_cad=False)) == []
    manifest = {"files": [{"path": "fixture.step"}], "corpus_files": 1,
                "settings": {"repeats": 1}}
    rows = [{"path": "fixture.step", "outcome": "succeeded", "wall_seconds": .7,
             "result": {**native, "validation_seconds": .02},
             "peak_process_working_set_bytes": None, "output_sha256": "fixture"}]
    assert summary(manifest, rows)["cases"][0]["timings_seconds"]["validation_seconds"]["median"] == .02
    print("CAD corpus runner legacy/native schema, refusal, budget, publication and validation-timing fixtures passed.")
    return 0


def run_case(command, timeout, stdout_path, stderr_path, archive_path, memory):
    started = time.perf_counter()
    peak = None
    timed_out = False
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr,
                                   creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        try:
            while True:
                measured = memory.peak(process)
                if measured is not None:
                    peak = max(peak or 0, measured)
                if process.poll() is not None:
                    break
                remaining = timeout - (time.perf_counter() - started)
                if remaining <= 0:
                    timed_out = True
                    process.kill()
                    break
                time.sleep(min(.02, remaining))
        except BaseException:
            process.kill()
            process.wait()
            raise
        process.wait()
        measured = memory.peak(process)
        if measured is not None:
            peak = max(peak or 0, measured)
    stdout_text, stdout_truncated = log_tail(stdout_path)
    stderr_text, stderr_truncated = log_tail(stderr_path)
    payload = last_report(stdout_text)
    exists = archive_path.is_file()
    valid = False
    if exists:
        with archive_path.open("rb") as stream:
            valid = stream.read(len(MAGIC)) == MAGIC
    outcome = classify(process.returncode, timed_out, payload, exists, valid)
    return {"outcome": outcome,
            "returncode": process.returncode, "timed_out": timed_out,
            "wall_seconds": time.perf_counter() - started, "timeout_seconds": timeout,
            "peak_process_working_set_bytes": peak,
            "process_memory_limit_enforced": False,
            "output_exists": exists, "archive_header_valid": valid,
            "output_publication_verified": outcome == "succeeded",
            "unexpected_output": exists and outcome != "succeeded",
            "output_bytes": archive_path.stat().st_size if exists else 0,
            "output_sha256": digest(archive_path) if exists else None,
            "mesh_published": payload.get("mesh_published") if payload else None,
            "archive_published": payload.get("archive_published") if payload else None,
            "result": payload, "stderr_tail": stderr_text,
            "stdout_tail": stdout_text if not payload else None,
            "stdout_tail_truncated": stdout_truncated,
            "stderr_tail_truncated": stderr_truncated}


def load_records(path):
    if not path.exists():
        return []
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    rows = []
    for index, line in enumerate(lines):
        try:
            rows.append(json.loads(line))
        except ValueError:
            if index + 1 == len(lines) and not line.endswith("\n"):
                # A killed driver can leave only its final JSONL write incomplete.
                path.write_text("".join(lines[:index]), encoding="utf-8")
                break
            raise ValueError(f"invalid evidence record at line {index + 1}") from None
    return rows


def summary(manifest, rows):
    cases = []
    for source in manifest["files"]:
        selected = [row for row in rows if row["path"] == source["path"]]
        good = [row for row in selected if row["outcome"] == "succeeded"]
        timings = {}
        for key in ("wall_seconds", *STAGES):
            timed_rows = selected if key == "wall_seconds" else good
            values = [row.get(key) if key == "wall_seconds" else
                      (row.get("result") or {}).get(key) for row in timed_rows]
            values = [value for value in values if isinstance(value, (int, float)) and
                      math.isfinite(value)]
            if values:
                timings[key] = {"median": statistics.median(values),
                                "minimum": min(values), "maximum": max(values)}
        peaks = [row["peak_process_working_set_bytes"] for row in selected
                 if row["peak_process_working_set_bytes"] is not None]
        cases.append({"path": source["path"], "runs": len(selected),
                      "outcomes": dict(Counter(row["outcome"] for row in selected)),
                      "archive_hash_stable": len({row["output_sha256"] for row in good}) == 1
                      if len(good) > 1 else None, "timings_seconds": timings,
                      "peak_process_working_set_bytes": max(peaks) if peaks else None})
    return {"schema": SCHEMA, "corpus_files": manifest["corpus_files"],
            "selected_files": len(manifest["files"]), "completed_runs": len(rows),
            "planned_runs": len(manifest["files"]) * manifest["settings"]["repeats"],
            "outcomes": dict(Counter(row["outcome"] for row in rows)), "cases": cases,
            "measurement": "Fresh sequential child per repeat; no warmup. Windows peak working set includes the complete direct child. Non-Windows peak is unavailable. Wall summaries include all attempts; stage summaries include successful runs only. Stage timers come from the CLI; generation combines surface and volume work. The timeout kills the direct child. Memory is measured, not limited. Archive magic checks publication presence, not independent mesh validity."}


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if argv == ["--self-test"]:
        return self_test()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("executable", type=Path)
    parser.add_argument("corpus", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=30, help="first N sorted paths")
    parser.add_argument("--smallest", type=int, default=6, help="additional K smallest files")
    parser.add_argument("--spread", type=int, default=0, help="additional evenly spaced sorted paths")
    parser.add_argument("--include", action="append", default=[], help="additional path relative to corpus")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("--size-factor", type=float, default=.15, help="target / imported bbox diagonal")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--max-file-bytes", type=int, default=64 * 1024 * 1024)
    parser.add_argument("--max-work", type=resource_cap, default=2_000_000,
                        help="positive logical-work cap; none omits the flag and inherits the CLI policy")
    parser.add_argument("--resource-policy", choices=("bounded", "normal"), default=None,
                        help="optional native resource policy; omitted preserves the selected CLI's default")
    parser.add_argument("--max-nodes", type=int, default=200_000)
    parser.add_argument("--max-cells", type=int, default=1_000_000)
    parser.add_argument("--max-bytes", type=int, default=128 * 1024 * 1024)
    parser.add_argument("--max-audit-bytes", type=int, default=128 * 1024 * 1024)
    parser.add_argument("--max-native-pieces", type=int, default=None, help="explicit native partition-piece admission; preserves the CLI default when omitted")
    parser.add_argument("--max-working-bytes", type=resource_cap, default=None, help="optional positive native working-payload cap; none inherits policy")
    parser.add_argument("--recovery-patch-states", type=int, default=None, help="optional bounded local-region search admission")
    parser.add_argument("--native-validation-threads", type=int, default=None, help="1..8 independent native CAD validation workers; meshing thread count is separate")
    parser.add_argument("--jobs", type=int, default=1, help="1..4 independent CAD processes; timings record this concurrency")
    parser.add_argument("--reuse-native-queries", action="store_true", help="reuse one per-call surface projector, with every point and trim check retained")
    parser.add_argument("--defer-discarded-native-facet-checks", action="store_true", help="defer facet sampling for discarded direct-star candidates; accepted results still complete it")
    parser.add_argument("--separate-matched-native-faces", action="store_true", help="explicit separate, exactly matched native assembly contact surfaces")
    parser.add_argument("--retain-auxiliary-geometry", action="store_true", help="retain valid standalone document geometry, while meshing every solid body")
    parser.add_argument("--native-kernel-admission", action="store_true", help="additional bounded concave fallback with original-CAD checks of each candidate fan")
    parser.add_argument("--longitudinal-recovery", action="store_true", help="optional native lengthwise partitions and bounded local Tet4 recovery")
    parser.add_argument("--independent-parts", action="store_true", help="preserve original parts and report interference; no global nonoverlap claim")
    parser.add_argument("--allow-centroids-within-chord-tolerance", action="store_true", help="explicit CAD chord policy for PLC centroids and interior samples; interior nodes remain strictly IN/ON")
    parser.add_argument("--local-partition-recovery", action="store_true", help="optional native local partition recovery")
    parser.add_argument("--adaptive-tet-recovery", action="store_true", help="optional native adaptive tetrahedron recovery")
    parser.add_argument("--trace-native-stages", action="store_true", help="record native stage traces without changing mesh validity settings")
    parser.add_argument("--recovery-mode", choices=("standard", "advanced"), default=None,
                        help="optional native recovery mode; omitted preserves the CLI default")
    parser.add_argument("--elements", choices=("tet4",), default=None,
                        help="explicit Tet4 native corpus profile; omitted preserves the CLI default")
    parser.add_argument("--native-cad", action="store_true", help="use the optional native OCCT adapter")
    parser.add_argument("--planar-recovery", action="store_true", help="enable bounded planar surface recovery")
    parser.add_argument("--general-recovery", action="store_true", help="enable bounded general tet edge recovery")
    parser.add_argument("--separate-contact-bodies", action="store_true", help="retain separate point/edge contact boundaries; refuse shared faces and overlap")
    parser.add_argument("--separate-matched-box-faces", action="store_true", help="native coarse boxes only: retain separate body nodes at exactly matched complete contact faces")
    parser.add_argument("--warning-quality", type=float, default=0, help="native warning threshold; does not reject valid low-quality cells")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    if min(args.limit, args.smallest, args.spread) < 0 or args.repeats < 1 or args.threads < 1:
        parser.error("selection counts must be nonnegative; repeats and threads must be positive")
    if not (args.timeout > 0 and math.isfinite(args.timeout) and
            args.size_factor > 0 and math.isfinite(args.size_factor)):
        parser.error("timeout and size factor must be finite and positive")
    if min(args.max_file_bytes, args.max_nodes, args.max_cells,
           args.max_bytes, args.max_audit_bytes) < 0:
        parser.error("resource budgets must be nonnegative")
    if (args.resource_policy is not None or args.recovery_mode is not None or args.elements is not None
            or any(getattr(args, key) for key in OPTIONAL_NATIVE_FLAGS)) and not args.native_cad:
        parser.error("resource policy and native recovery controls require --native-cad")
    if not math.isfinite(args.warning_quality) or not 0 <= args.warning_quality <= 1 or (args.warning_quality > 0 and not args.native_cad):
        parser.error("warning quality must be in [0,1] and requires --native-cad")
    if args.max_native_pieces is not None and (args.max_native_pieces < 1 or not args.native_cad):
        parser.error("native piece admission must be positive and requires --native-cad")
    if args.native_validation_threads is not None and (not 1 <= args.native_validation_threads <= 8 or not args.native_cad):
        parser.error("native validation workers must be in [1,8] and require --native-cad")
    if args.max_working_bytes is not None and not args.native_cad:
        parser.error("native working admission requires --native-cad")
    if args.recovery_patch_states is not None and (not 0<=args.recovery_patch_states<=65536 or not args.native_cad):
        parser.error("region search admission must be in [0,65536] and requires --native-cad")
    if not 1 <= args.jobs <= 4:
        parser.error("jobs must be in [1,4]")
    if args.separate_matched_box_faces and not args.native_cad:
        parser.error("matched box-face contacts require --native-cad")
    if args.reuse_native_queries and not args.native_cad:
        parser.error("native query reuse requires --native-cad")
    if (args.defer_discarded_native_facet_checks or args.separate_matched_native_faces or args.retain_auxiliary_geometry or args.native_kernel_admission or args.longitudinal_recovery or args.independent_parts or args.allow_centroids_within_chord_tolerance) and not args.native_cad:
        parser.error("native fallback controls require --native-cad")
    executable, root, output = args.executable.resolve(), args.corpus.resolve(), args.output.resolve()
    if not executable.is_file() or not root.is_dir():
        parser.error("executable and corpus directory must exist")
    if output == root or output.is_relative_to(root):
        parser.error("evidence output must be outside the input corpus")
    settings = {key: value for key, value in vars(args).items()
                if key not in {"executable", "corpus", "output", "resume"}
                and (key not in (*MODE_FLAGS, *OPTIONAL_NATIVE_FLAGS) or value) and (key != "warning_quality" or value)
                and (key not in {"max_working_bytes","recovery_patch_states"} or value is not None)
                and (key not in {"resource_policy", "recovery_mode", "elements"} or value is not None)
                and (key != "max_native_pieces" or value is not None)}
    output.mkdir(parents=True, exist_ok=True)
    manifest_path, records_path = output / "manifest.json", output / "results.jsonl"
    identity = {"corpus_root": str(root), "executable": str(executable),
                "executable_sha256": digest(executable), "settings": settings}
    if args.resume:
        if not manifest_path.exists():
            parser.error("--resume requires an existing manifest.json")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("schema") != SCHEMA or any(manifest.get(k) != v for k, v in identity.items()):
            parser.error("resume requires identical executable, corpus and settings")
        for source in manifest["files"]:
            path = (root / source["path"]).resolve()
            if not path.is_relative_to(root) or digest(path) != source["sha256"]:
                parser.error(f"input changed since manifest creation: {source['path']}")
        rows = load_records(records_path)
    else:
        if manifest_path.exists() or records_path.exists():
            parser.error("evidence already exists; use --resume or a new output directory")
        count, files = select_files(root, args.limit, args.smallest, args.spread, args.include)
        if not files:
            parser.error("selection contains no CAD files")
        manifest = {"schema": SCHEMA, **identity, "corpus_files": count, "files": files,
                    "created_utc": datetime.now(timezone.utc).isoformat(),
                    "platform": platform.platform(), "python": platform.python_version()}
        save_json(manifest_path, manifest)
        rows = []
    completed = {(row["path"], row["repeat"]) for row in rows}
    suppress_native_error_dialogs()
    memory = WindowsMemory()
    prefix = [sys.executable, str(executable)] if executable.suffix.lower() == ".py" else [str(executable)]
    tasks = []
    for index, source in enumerate(manifest["files"]):
        if digest(root / source["path"]) != source["sha256"]:
            raise ValueError(f"input changed since manifest creation: {source['path']}")
        for repeat in range(args.repeats):
            if (source["path"], repeat) in completed:
                continue
            tasks.append((index, source, repeat))
    def execute(task):
            index, source, repeat = task
            run_id = f"{index:04d}-{repeat:02d}-{uuid.uuid4().hex[:8]}"
            directory = output / "runs" / run_id
            directory.mkdir(parents=True)
            archive = directory / "mesh.paacsmesh"
            command = prefix + [str(root / source["path"]), "1", str(archive), "--report-json",
                                "--size-factor", str(args.size_factor), "--threads", str(args.threads)]
            command += ["--max-file-bytes", str(args.max_file_bytes)]
            command += work_arguments(args)
            for key in ("max_nodes", "max_cells", "max_bytes", "max_audit_bytes"):
                command += ["--" + key.replace("_", "-"), str(getattr(args, key))]
            command += mode_arguments(args)
            if args.warning_quality > 0:
                command += ["--warning-quality", str(args.warning_quality)]
            if args.max_native_pieces is not None:
                command += ["--max-native-pieces", str(args.max_native_pieces)]
            if args.native_validation_threads is not None:
                command += ["--native-validation-threads", str(args.native_validation_threads)]
            command += native_profile_arguments(args)
            for key in ("max_working_bytes","recovery_patch_states"):
                if getattr(args,key) is not None:
                    command += ["--"+key.replace("_","-"),str(getattr(args,key))]
            row = run_case(command, args.timeout, directory / "stdout.log", directory / "stderr.log",
                           archive, memory)
            if args.native_cad:
                reported_resources = (row.get("result") or {}).get("resources")
                row.update(resource_policy_requested=args.resource_policy,
                           max_work_requested=args.max_work,
                           max_working_bytes_requested=args.max_working_bytes,
                           resources_reported=reported_resources)
            row.update(path=source["path"], input_sha256=source["sha256"], repeat=repeat,
                       command=command, run_directory=str(directory.relative_to(output)),
                       output=str(archive.relative_to(output)))
            row["input_unchanged"] = digest(root / source["path"]) == source["sha256"]
            if not row["input_unchanged"]:
                row["outcome"] = "publication-violation"
            return row
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        pending = [executor.submit(execute, task) for task in tasks]
        for future in as_completed(pending):
            row = future.result()
            rows.append(row)
            with records_path.open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(row, allow_nan=False) + "\n")
                stream.flush()
            save_json(output / "report.json", summary(manifest, rows))
            print(f"[{len(rows)}/{len(manifest['files']) * args.repeats}] {row['path']} "
                  f"repeat={row['repeat'] + 1}: {row['outcome']} ({row['wall_seconds']:.3f}s)", flush=True)
    report = summary(manifest, rows)
    save_json(output / "report.json", report)
    print(json.dumps(report["outcomes"], sort_keys=True))
    return 1 if any(row["outcome"] in {"crash", "timeout", "publication-violation",
                                      "invalid-result", "missing-report", "failed"} for row in rows) else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as error:
        sys.exit(f"CAD corpus runner: {error}")
