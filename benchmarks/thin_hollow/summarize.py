"""Pair retained hollow-fixture trial receipts without loading meshes or Native code.

Inputs may be qualifier results.json files or individual controller receipts.
Every distinct controller is retained, including refusals and timed-out work.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
import os
from pathlib import Path
import time

PHASES = ("preparation", "generation", "validation", "audit", "quality")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value):
    return value if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) and value >= 0 else None


class Reader:
    def __init__(self, args):
        self.args, self.started, self.bytes = args, time.perf_counter(), 0
        self.cache = {}

    def load(self, path):
        path = Path(path).resolve(strict=True)
        if path in self.cache:
            return self.cache[path]
        require(len(self.cache) < self.args.max_files, "Input file allowance exhausted")
        size = path.stat().st_size
        require(size <= self.args.max_file_bytes and self.bytes+size <= self.args.max_input_bytes, "JSON byte allowance exhausted")
        payload = path.read_bytes()
        require(len(payload) == size, "Input changed during its read")
        value = json.loads(payload.decode("utf-8-sig"), parse_constant=lambda v: require(False, "Nonfinite JSON: "+v))
        self.bytes += size
        result = (value, hashlib.sha256(payload).hexdigest())
        self.cache[path] = result
        self.deadline()
        return result

    def deadline(self):
        require(time.perf_counter()-self.started <= self.args.max_seconds, "Summary deadline exhausted")

    def unchanged(self):
        for path, (_, digest) in self.cache.items():
            self.deadline()
            with path.open("rb") as stream:
                require(hashlib.file_digest(stream, "sha256").hexdigest() == digest, "Input JSON changed: "+str(path))


def expand(paths, reader):
    controllers, indexes = [], []
    for path in paths:
        path = Path(path).resolve(strict=True)
        value, digest = reader.load(path)
        if "attempts" in value:
            require(len(value["attempts"]) <= reader.args.max_attempts, "Result index exceeds attempt allowance")
            indexes.append({"path": str(path), "sha256": digest, "complete": value.get("complete"),
                            "recorded_total_process_seconds": value.get("total_process_seconds")})
            controllers.extend(Path(row["receipt"]) for row in value["attempts"])
        else:
            require("job" in value, "Unsupported trial receipt")
            controllers.append(path)
    require(len(controllers) <= reader.args.max_attempts, "Attempt allowance exhausted")
    return controllers, indexes


def attempt(path, reader):
    path = Path(path).resolve(strict=True)
    controller, digest = reader.load(path)
    job, worker = controller["job"], controller.get("worker", {})
    report, report_ref, integrity = None, None, []
    report_path = Path(job["output"])/worker.get("native_report_file", "native-report.json")
    if report_path.exists():
        report, report_sha = reader.load(report_path)
        report_ref = {"path": str(report_path), "sha256": report_sha}
        if report_sha != worker.get("native_report_sha256"):
            integrity.append("native_report_hash_mismatch")
    if worker.get("source_sha256") != job.get("source_sha256") or (report and report.get("source_sha256") != job.get("source_sha256")):
        integrity.append("source_hash_mismatch")
    if worker.get("native_sha256") != job.get("native_sha256"):
        integrity.append("native_identity_mismatch")
    guarded = all(value is True for value in (controller.get("input_guards_unchanged"), worker.get("input_guards_unchanged")))
    if not guarded:
        integrity.append("completed_input_guards_unavailable")
    if controller.get("worker_pid_verified") is not True:
        integrity.append("worker_pid_unverified")
    report = report or {}
    if report.get("mesh_published") is True and any(sorted(worker.get("original_"+k) or []) != sorted(report.get(k) or []) for k in ("body_tags", "face_tags")):
        integrity.append("published_cad_tags_differ")
    native_gates = report.get("status") == "succeeded" and report.get("mesh_published") is True and all(
        report.get(section, {}).get(field) is True for section in ("audit", "quality") for field in ("complete", "valid"))
    accepted = controller.get("status") == "native_accepted" and controller.get("state") == "complete" and not controller.get("hard_timeout") and guarded and native_gates and not integrity
    quality = report.get("quality", {})
    quality_available = quality.get("complete") is True and quality.get("valid") is True
    quality_fields = ("minimum", "percentile_1", "percentile_5", "median", "mean", "minimum_volume", "total_volume", "inverted_cells")
    body_fields = ("solid_tag", "route", "reason", "surface_reason", "volume_reason", "surface_triangles", "surface_attempts",
                   "relative_volume_error", "maximum_boundary_node_error", "maximum_boundary_sample_error", "surface_chord_error_limit",
                   "boundary_embedding_certified", "native_centroid_source_guard", "volume_sizing", "initial_front", "interior_front",
                   "route_resolution", "initial_front_fallback", "source_shell", "partition_attempt_reports", "constraint_recovery",
                   "native_surface_seconds", "native_volume_seconds", "native_partition_seconds", "native_assembly_seconds",
                   "initial_topology_quality", "material_classification", "region_recovery",
                   "native_cap_sizing", "local_delaunay_quality_edits")
    bodies = [{key: body.get(key) for key in body_fields} for body in report.get("bodies", [])]
    seconds = report.get("seconds", {})
    native_total = finite(seconds.get("meshing_total"))
    phase_sum = math.fsum(seconds[k] for k in PHASES) if all(finite(seconds.get(k)) is not None for k in PHASES) else None
    kwargs = job.get("kwargs", {})
    size = finite(job.get("global_h_mm"))
    require(size is not None and size > 0, "Positive requested size unavailable")
    key = [job["id"], size, kwargs.get("assembly"), kwargs.get("max_chord_error")]
    return {"attempt_id": digest, "controller": str(path), "key": key, "profile": job.get("profile"),
            "status": controller.get("status", "unknown"), "native_status": worker.get("native_status", "unknown"),
            "reason": report.get("reason", worker.get("reason")), "worker_phase": worker.get("phase"),
            "qualified_native_acceptance": bool(accepted), "native_gates_passed": bool(native_gates), "integrity_flags": integrity,
            "hard_timeout": controller.get("hard_timeout"), "worker_pid_verified": controller.get("worker_pid_verified"),
            "source_sha256": job.get("source_sha256"), "native_sha256": job.get("native_sha256"),
            "archive_sha256_recorded": worker.get("archive_sha256"), "native_report": report_ref,
            "original_body_tags": worker.get("original_body_tags"), "original_face_tags": worker.get("original_face_tags"),
            "reported_body_tags": report.get("body_tags"), "reported_face_tags": report.get("face_tags"),
            "fixture_specification": job.get("specification"), "requested_controls": kwargs,
            "requested_size_mm": size, "effective_sizing_report": report.get("sizing"),
            "achieved_size_report": report.get("achieved_size"),
            "achieved_size_available": isinstance(report.get("achieved_size"), dict),
            "published": report.get("mesh_published"), "counts": {k: report.get(k) for k in ("nodes", "cells", "tetrahedra")},
            "quality": {"available": quality_available, "metric": quality.get("metric"),
                        "reported_histogram": quality.get("histogram"),
                        "measurements": {k: finite(quality.get(k)) if quality_available else None for k in quality_fields}},
            "route_resolution": report.get("route_resolution"), "initial_front_fallback": report.get("initial_front_fallback"),
            "clocks": {"controller_process_seconds": finite(controller.get("process_seconds")),
                       "worker_seconds": finite(worker.get("worker_seconds")), "package_import_seconds": finite(worker.get("package_import_seconds")),
                       "step_import_seconds": finite(worker.get("step_import_seconds")), "public_mesh_call_seconds": finite(worker.get("mesh_call_seconds")),
                       "export_seconds": finite(worker.get("export_seconds")), "native_meshing_total_seconds": native_total,
                       "native_reported_seconds": {k: finite(v) for k, v in seconds.items()}, "native_top_phase_sum_seconds": phase_sum,
                       "native_unassigned_top_seconds": native_total-phase_sum if native_total is not None and phase_sum is not None else None},
            "peak_worker_bytes": controller.get("peak_process_working_set_bytes") if controller.get("worker_pid_verified") is True else None,
            "body_reports": bodies}


def comparison(first, second):
    equal = lambda field: None if first.get(field) is None or second.get(field) is None else first.get(field) == second.get(field)
    controls = sorted(set(first["requested_controls"]) | set(second["requested_controls"]))
    changes = {k: {"baseline": first["requested_controls"].get(k), "candidate": second["requested_controls"].get(k)}
               for k in controls if first["requested_controls"].get(k) != second["requested_controls"].get(k)}
    qc = ("minimum", "percentile_1", "percentile_5", "median", "mean")
    quality_delta = {k: second["quality"]["measurements"][k]-first["quality"]["measurements"][k]
                     if first["quality"]["measurements"][k] is not None and second["quality"]["measurements"][k] is not None else None for k in qc}
    clock_deltas = {k: second["clocks"][k]-first["clocks"][k] if first["clocks"][k] is not None and second["clocks"][k] is not None else None
                    for k in ("controller_process_seconds", "public_mesh_call_seconds", "native_meshing_total_seconds")}
    cad_fields = ("solid_tag", "relative_volume_error", "maximum_boundary_node_error", "maximum_boundary_sample_error", "surface_chord_error_limit")
    cad = lambda row: [{k: b.get(k) for k in cad_fields} for b in row["body_reports"]]
    a, b = first["qualified_native_acceptance"], second["qualified_native_acceptance"]
    outcome = "preserved_acceptance" if a and b else "gained_acceptance" if b else "regressed_acceptance" if a else "neither_accepted"
    return {"baseline_attempt": first["attempt_id"], "candidate_attempt": second["attempt_id"], "outcome": outcome,
            "both_qualified": a and b, "same_attempt_reference": first["attempt_id"] == second["attempt_id"],
            "source_sha256_equal": equal("source_sha256"), "fixture_specification_equal": equal("fixture_specification"),
            "original_body_tags_equal": equal("original_body_tags"), "original_face_tags_equal": equal("original_face_tags"),
            "recorded_archive_sha256_equal": equal("archive_sha256_recorded"),
            "counts_equal": first["counts"] == second["counts"] if all(v is not None for r in (first, second) for v in r["counts"].values()) else None,
            "requested_control_changes": changes, "cad_measurements_equal": cad(first) == cad(second) if a and b and cad(first) and cad(second) else None,
            "effective_sizing_report_equal": equal("effective_sizing_report"),
            "cad_measurements": {"baseline": cad(first), "candidate": cad(second)},
            "achieved_size_report_equal": equal("achieved_size_report"), "achieved_size_comparison_available": first["achieved_size_available"] and second["achieved_size_available"],
            "quality_candidate_minus_baseline": quality_delta, "clock_candidate_minus_baseline_seconds": clock_deltas,
            "comparison_scope": "Receipt-only measurements. Equal requested h does not establish equal achieved size/CAD error. Distinct meshes need independent archive size/quality inspection before a comparable-mesh runtime claim."}


def ledger(rows):
    clocks = ("controller_process_seconds", "public_mesh_call_seconds", "native_meshing_total_seconds")
    return {"attempts": len(rows), "accepted": sum(r["qualified_native_acceptance"] for r in rows), "statuses": dict(Counter(r["status"] for r in rows)),
            "clocks": {k: {"known_subtotal_seconds": math.fsum(r["clocks"][k] for r in rows if r["clocks"][k] is not None),
                            "known_attempts": sum(r["clocks"][k] is not None for r in rows), "unknown_attempts": sum(r["clocks"][k] is None for r in rows)} for k in clocks}}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline", nargs="+", type=Path, required=True); p.add_argument("--candidate", nargs="+", type=Path, required=True)
    p.add_argument("--development", nargs="+", type=Path, default=[], help="Retained earlier trials: cost ledger only, not baseline/candidate pairing")
    p.add_argument("--output", type=Path, required=True); p.add_argument("--max-seconds", type=float, default=30)
    for name, default in (("max-files",2048),("max-attempts",1000),("max-file-bytes",16*1024**2),("max-input-bytes",128*1024**2),("max-output-bytes",16*1024**2)):
        p.add_argument("--"+name, type=int, default=default)
    args = p.parse_args()
    require(all(v > 0 and math.isfinite(v) for k, v in vars(args).items() if k.startswith("max_")), "Positive finite bounds required")
    output = args.output.resolve(); require(not output.exists(), "Choose a fresh summary JSON path")
    reader = Reader(args); groups, indexes, all_rows, aliases = {}, {}, {}, []
    for group, paths in (("baseline",args.baseline),("candidate",args.candidate),("development",args.development)):
        controllers, indexes[group] = expand(paths, reader); rows = {}
        for path in controllers:
            row = attempt(path, reader); digest = row["attempt_id"]
            if digest in rows:
                aliases.append({"group": group, "controller": str(path), "attempt_id": digest, "counted_again": False})
            rows[digest] = row; all_rows[digest] = row
        groups[group] = list(rows.values())
    by_key = defaultdict(lambda: {"baseline":[],"candidate":[]})
    for group, rows in groups.items():
        if group == "development":
            continue
        for row in rows:
            by_key[json.dumps(row["key"],separators=(",",":"),allow_nan=False)][group].append(row)
    pairs = []
    for key, rows in sorted(by_key.items()):
        paired = len(rows["baseline"]) == len(rows["candidate"]) == 1
        pairs.append({"key":json.loads(key), "baseline_attempt_ids":[r["attempt_id"] for r in rows["baseline"]],
                      "candidate_attempt_ids":[r["attempt_id"] for r in rows["candidate"]],
                      "pairing":"one_to_one" if paired else "missing_or_multiple_attempts_no_best_attempt_selected",
                      "comparison":comparison(rows["baseline"][0],rows["candidate"][0]) if paired else None})
    reader.unchanged()
    report = {"schema":"paacs-hollow-trial-comparison-v1", "native_execution":False, "archives_or_source_geometry_read":False,
              "pair_key":["fixture_id","requested_h_mm","assembly","requested_chord_mm"], "input_indexes":indexes,
              "groups":groups,"pairs":pairs,"duplicate_index_references":aliases,
              "shared_attempt_ids_between_groups":sorted({r["attempt_id"] for r in groups["baseline"]} & {r["attempt_id"] for r in groups["candidate"]}),
              "group_ledgers":{g:ledger(r) for g,r in groups.items()},"all_unique_attempts_ledger":ledger(list(all_rows.values())),
              "clock_scope":"All distinct attempts, including refused/time-out costs, are included once per ledger. Public call and Native total include recorded internal retries. Nested source/route/partition timers are retained descriptively and must not be added to totals.",
              "unknowns":"Missing clocks/quality/achieved sizing remain null, never zero. No best retry is chosen and no equal-mesh, solver-accuracy or global-CAD claim is inferred.",
              "source_files":{str(path):digest for path,(_,digest) in reader.cache.items()},"elapsed_seconds":time.perf_counter()-reader.started}
    encoded = (json.dumps(report,indent=2,allow_nan=False)+"\n").encode()
    require(len(encoded) <= args.max_output_bytes, "Summary output allowance exhausted")
    output.parent.mkdir(parents=True,exist_ok=True); temporary = output.with_name(output.name+".tmp")
    temporary.write_bytes(encoded); os.replace(temporary,output)
    print(json.dumps({"output":str(output),"baseline_attempts":len(groups["baseline"]),"candidate_attempts":len(groups["candidate"]),"pair_keys":len(pairs),"native_execution":False}))


if __name__ == "__main__":
    try:
        main()
    except (OSError,ValueError,KeyError) as error:
        raise SystemExit("Receipt summary refused: "+str(error))
