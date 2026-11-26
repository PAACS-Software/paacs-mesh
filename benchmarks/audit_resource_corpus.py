"""Independently screen every published Tet4 archive in a retained corpus run.

No CAD is imported and no meshing or native CAD re-audit is performed. Archive
coordinates/connectivity and dimensionless Tet4 metrics are independently
checked by the existing inspection API. CAD checks remain saved native report
evidence. Keep outputs local when source CAD redistribution is not cleared:
optional worst-cell details contain actual mesh coordinates.

The optional plan has a files list keyed by path and sha256, with optional
expected_original_solids and expected_body_tags. Counts/IDs are compared only
when recorded; matching classifications do not reclassify points against CAD.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
from collections import Counter
from dataclasses import dataclass
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import stat
import sys
import tempfile
import time


SCHEMA = "paacs-resource-corpus-archive-audit-v1"
PROGRESS_SCHEMA = "paacs-resource-corpus-archive-audit-progress-v1"
NATIVE_SCHEMA = "paacs-native-cad-report-v1"
SAVED_BASIS = "saved_native_cli_report_not_recertified"
METRICS = {"mean_ratio": 1., "minimum_dihedral_deg": 180., "scaled_jacobian": 1.}
AUDIT_COUNTERS = ("boundary_missing", "boundary_unexpected", "nonmanifold_faces",
                  "inconsistent_face_orientation")
VOLUME_FAMILIES = ("tet10", "prism6", "hex8", "pyramid5")


class AuditFailure(ValueError):
    def __init__(self, reason, message):
        self.reason = reason
        super().__init__(message)


@dataclass(frozen=True)
class Limits:
    max_bytes: int = 1024**3
    max_cells: int = 5_000_000
    max_nodes: int = 2_000_000
    max_work: int = 100_000_000
    max_seconds: float = 120
    max_total_seconds: float = 3600
    max_metadata_bytes: int = 32 * 1024**2
    max_records: int = 4096
    max_output_bytes: int = 64 * 1024**2
    detail_limit: int = 10

    def validate(self):
        for name, value in vars(self).items():
            minimum = 0 if name == "detail_limit" else 1
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < minimum:
                raise ValueError(f"{name} must be finite and at least {minimum}")
            if "seconds" not in name and not isinstance(value, int):
                raise ValueError(f"{name} must be an integer")
        if self.detail_limit > 1000:
            raise ValueError("detail_limit must not exceed the inspection page limit of 1000")


def _require(condition, reason, message):
    if not condition:
        raise AuditFailure(reason, message)


def _int(value, name, minimum=0):
    _require(isinstance(value, int) and not isinstance(value, bool) and value >= minimum,
             "metadata_invalid", f"{name} must be an integer >= {minimum}")
    return value


def _hash(value, name):
    _require(isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value),
             "metadata_invalid", f"{name} must be a lowercase SHA-256")
    return value


def _regular(path):
    _require(not path.is_symlink() and not getattr(path, "is_junction", lambda: False)(),
             "input_not_regular", "Inputs must not be symbolic links or junctions")
    _require(stat.S_ISREG(path.stat().st_mode), "input_not_regular", "Input must be a regular file")


def _inside(root, relative):
    _require(isinstance(relative, str) and relative and "\0" not in relative,
             "path_invalid", "A relative output path is required")
    windows, posix = PureWindowsPath(relative), PurePosixPath(relative)
    _require(not windows.drive and not windows.root and not posix.is_absolute()
             and ".." not in windows.parts and ".." not in posix.parts
             and not any(":" in part for part in windows.parts), "path_invalid", "Output path escapes the run")
    target = root / relative
    _require(target.resolve().is_relative_to(root.resolve()), "path_invalid", "Output path escapes the run")
    cursor = target
    while cursor != root:
        _require(not cursor.is_symlink() and not getattr(cursor, "is_junction", lambda: False)(),
                 "path_invalid", "Output path contains a link or junction")
        cursor = cursor.parent
    return target


def _deadline(started, seconds):
    _require(time.perf_counter() - started <= seconds, "time_budget", "Cooperative audit deadline exceeded")


def _digest(path, maximum, started=None, seconds=None):
    _regular(path)
    size = path.stat().st_size
    _require(size <= maximum, "input_budget", "Input exceeds the admitted byte allowance")
    result = hashlib.sha256()
    with path.open("rb") as stream:
        remaining = size
        while remaining:
            if started is not None:
                _deadline(started, seconds)
            block = stream.read(min(1024**2, remaining))
            _require(bool(block), "input_changed", "Input shortened while hashing")
            result.update(block)
            remaining -= len(block)
        _require(not stream.read(1), "input_changed", "Input grew while hashing")
    return result.hexdigest()


def _pairs(values):
    result = {}
    for key, value in values:
        _require(key not in result, "metadata_invalid", "Duplicate JSON field")
        result[key] = value
    return result


def _float(value):
    result = float(value)
    _require(math.isfinite(result), "metadata_invalid", "Nonfinite JSON number")
    return result


def _reject(value):
    raise AuditFailure("metadata_invalid", f"Nonfinite JSON constant {value}")


def _decode(payload):
    try:
        result = json.loads(payload.decode("utf-8"), object_pairs_hook=_pairs,
                            parse_float=_float, parse_constant=_reject)
    except (UnicodeError, json.JSONDecodeError, RecursionError) as error:
        raise AuditFailure("metadata_invalid", "Metadata must be finite UTF-8 JSON within parser depth") from error
    _require(isinstance(result, dict), "metadata_invalid", "Metadata must be a JSON object")
    return result


def _json_file(path, maximum):
    _regular(path)
    _require(path.stat().st_size <= maximum, "metadata_budget", "Metadata exceeds its parsing allowance")
    with path.open("rb") as stream:
        payload = stream.read(maximum + 1)
    _require(len(payload) <= maximum, "metadata_budget", "Metadata grew beyond its parsing allowance")
    return _decode(payload), hashlib.sha256(payload).hexdigest()


def check_native(report, report_kind="cli"):
    _require(isinstance(report, dict), "missing_native_report", "Published archive has no saved native report")
    _require(report_kind in ("cli", "python"), "native_report_schema", "Unknown saved report format")
    if report_kind == "cli":
        _require(report.get("schema") == NATIVE_SCHEMA, "native_report_schema", "Expected a NativeCLI report")
    _require(report.get("status") == "succeeded" and report.get("reason") == "none"
             and report.get("mesh_published") is True,
             "native_report_refusal", "Saved report does not declare successful publication")
    _require(report.get("audit_valid") is True, "native_audit_not_valid", "Saved native audit did not pass")
    _int(report.get("nodes"), "reported nodes", 1)
    _int(report.get("cells"), "reported cells", 1)
    if report_kind == "cli":
        _require(report.get("archive_published") is True, "native_report_refusal", "Saved CLI archive was not published")
        for name in AUDIT_COUNTERS:
            _require(_int(report.get(name), name) == 0, "native_audit_not_valid", f"Saved {name} is nonzero")
        solids = _int(report.get("solids"), "solids", 1)
    else:
        _hash(report.get("source_sha256"), "native source_sha256")
        _require(report.get("units") == "mm" and isinstance(report.get("kernel_version"), str)
                 and bool(report["kernel_version"]), "native_report_schema", "Python report lacks source units or kernel identity")
        audit, quality = report.get("audit"), report.get("quality")
        _require(isinstance(audit, dict) and audit.get("complete") is True and audit.get("valid") is True,
                 "native_audit_not_valid", "Saved Python native audit is incomplete or invalid")
        for name in ("missing_boundary", "unexpected_boundary", "nonmanifold_faces", "inconsistent_orientation",
                     "invalid_connectivity", "nonfinite_nodes"):
            _require(_int(audit.get(name), name) == 0, "native_audit_not_valid", f"Saved {name} is nonzero")
        _require(isinstance(quality, dict) and quality.get("complete") is True and quality.get("valid") is True
                 and _int(quality.get("inverted_cells"), "inverted_cells") == 0,
                 "native_quality_not_valid", "Saved Python native quality checks did not pass")
        declared_tags = report.get("body_tags")
        _require(isinstance(declared_tags, list) and bool(declared_tags), "body_report_mismatch", "Missing Python body identities")
        solids = len(declared_tags)
    bodies = report.get("bodies")
    _require(isinstance(bodies, list) and len(bodies) == solids and all(isinstance(b, dict) for b in bodies),
             "body_report_mismatch", "Saved bodies do not match the reported original solid count")
    tags = [_int(body.get("solid_tag"), "solid_tag", 1) for body in bodies]
    _require(len(set(tags)) == len(tags) and all(body.get("reason") == "none" for body in bodies),
             "body_report_mismatch", "Saved body identities are duplicated or refused")
    if report_kind == "python":
        _require(sorted(_int(tag, "body tag", 1) for tag in declared_tags) == sorted(tags),
                 "body_report_mismatch", "Python body identities differ from body reports")
    return sorted(tags)


def check_summary(summary, report, expected=None, report_kind="cli"):
    """Check whole-snapshot aggregates; never infer original CAD membership."""
    _require(summary.get("schema") == "paacs-mesh-quality-summary-v1" and summary.get("scope") == "entire_snapshot",
             "inspection_summary", "Expected a whole-snapshot inspection summary")
    counts = summary.get("cell_counts", {})
    _require(isinstance(counts, dict), "inspection_summary", "Missing cell counts")
    tet4 = _int(counts.get("tet4"), "tet4", 1)
    _require(not any(_int(counts.get(name, 0), name) for name in VOLUME_FAMILIES),
             "unsupported_volume_cells", "Published output contains unsupported volume cell families")
    _require(summary.get("node_count") == report.get("nodes") and summary.get("cell_count") == report.get("cells"),
             "report_counts", "Saved report counts differ from the archive")
    if report_kind == "python":
        _require(report.get("tetrahedra") == tet4, "report_counts", "Saved Python Tet4 count differs from the archive")
    orientation = summary.get("orientation_counts", {})
    _require(orientation.get("positive") == tet4 and orientation.get("inverted") == 0
             and orientation.get("degenerate") == 0, "nonpositive_orientation", "Some Tet4 orientations are not positive")
    _require(orientation.get("unsupported") == summary["cell_count"] - tet4,
             "inspection_summary", "Orientation and supported-cell counts differ")
    _require(_int(summary.get("range_limited_cells"), "range_limited_cells") == 0,
             "metric_range_loss", "Some nonzero dimensional metrics are unavailable or unrepresentable")
    quality = summary.get("quality", {})
    for name, upper in METRICS.items():
        population = quality.get(name, {})
        _require(population.get("count") == tet4 and population.get("unavailable") == 0,
                 "metric_unavailable", f"{name} is not available for every Tet4")
        for key in ("minimum", "maximum", "mean", "median", "p01", "p05", "p95"):
            value = population.get(key)
            _require(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value),
                     "metric_nonfinite", f"{name}.{key} is nonfinite or missing")
        _require(0 < population["minimum"] <= population["maximum"] <= upper + 1e-12,
                 "metric_nonpositive_or_out_of_range", f"{name} contains zero, negative or out-of-range values")
    tags = check_native(report, report_kind)
    body_counts = summary.get("body_tet4_counts", {})
    _require(isinstance(body_counts, dict) and summary.get("unclassified_tet4_count") == 0,
             "body_classification", "Tet4 body classifications are unavailable")
    measured = {int(tag): _int(count, "body tetrahedra", 1) for tag, count in body_counts.items()}
    _require(sorted(measured) == tags and sum(measured.values()) == tet4,
             "body_classification", "Archive classifications do not cover exactly the recorded body tags")
    expected = expected or {}
    original_count = expected.get("expected_original_solids")
    if original_count is not None:
        _require(len(tags) == _int(original_count, "expected_original_solids", 1),
                 "original_body_count", "Original body count was not preserved")
    original_tags = expected.get("expected_body_tags")
    if original_tags is not None:
        _require(isinstance(original_tags, list) and sorted(_int(t, "expected_body_tag", 1) for t in original_tags) == tags,
                 "original_body_tags", "Known original body IDs were not preserved")
    return {"recorded_body_tags": tags, "body_tet4_counts": body_counts,
            "expected_original_solids": original_count, "known_original_body_tags": original_tags,
            "original_count_checked": original_count is not None, "original_ids_checked": original_tags is not None,
            "classification_tags_match_saved_body_tags": True,
            "original_body_geometry_reclassification_performed": False}


def _saved_report(report, report_kind="cli"):
    result = {key: report.get(key) for key in ("schema", "status", "reason", "mesh_published", "archive_published",
        "solids", "validated_parts", "mesh_scope", "interbody_nonoverlap_certified", "audit_valid",
        *AUDIT_COUNTERS, "minimum_quality", "quality_mean", "native_validation_threads", "total_seconds",
        "preparation_seconds", "generation_seconds", "validation_seconds", "audit_seconds", "quality_seconds",
        "millimetres_per_model_unit", "interferences")}
    result["bodies"] = [{key: body.get(key) for key in ("solid_tag", "route", "reason", "relative_volume_error",
        "boundary_embedding_certified", "boundary_node_error", "boundary_sample_error", "centroid_outside_distance_limit",
        "maximum_outside_centroid_distance", "native_centroids_checked", "native_centroid_source_guard")}
        for body in report.get("bodies", [])]
    if report_kind == "python":
        result.update(audit=report["audit"], quality=report["quality"], units=report["units"],
                      source_sha256=report["source_sha256"], kernel_version=report["kernel_version"],
                      seconds=report.get("seconds"))
    result.update(basis=SAVED_BASIS if report_kind == "cli" else "saved_native_python_report_not_recertified",
                  saved_report_format=report_kind, independent_native_reaudit_performed=False)
    return result


def audit_record(run, source, row, expected, limits, inspect_mesh, *, report_kind="cli"):
    started = time.perf_counter()
    result = {"path": source["path"], "source_sha256": source["sha256"], "repeat": row.get("repeat"),
              "native_outcome": row.get("outcome"), "mesh_wall_seconds": row.get("wall_seconds"),
              "status": "failed", "independent_native_reaudit_performed": False,
              "global_intersection_audit_performed": False}
    try:
        _require(row.get("input_sha256") == source["sha256"], "source_hash", "Runner source hash differs from manifest")
        if expected is not None:
            _require(expected.get("sha256") == source["sha256"], "source_hash", "Source hash differs from the fixed plan")
        archive = _inside(run, row.get("output"))
        report = row.get("result")
        if row.get("outcome") != "succeeded":
            _require(not archive.exists() and not row.get("output_exists") and not row.get("output_sha256"),
                     "publication_violation", "A failed native attempt left an output archive")
            result.update(status="not_published", reason=row.get("outcome"),
                          saved_native_reason=report.get("reason") if isinstance(report, dict) else None,
                          archive_screen_performed=False)
            return result
        check_native(report, report_kind)
        if report_kind == "cli":
            _require(row.get("input_unchanged") is True, "source_changed", "Source was not retained unchanged during meshing")
        target_hash = _hash(row.get("output_sha256"), "output_sha256")
        _require(_digest(archive, limits.max_bytes // 2, started, limits.max_seconds) == target_hash,
                 "archive_hash", "Published archive hash changed")
        report_hash = hashlib.sha256(json.dumps(report, sort_keys=True, separators=(",", ":"),
                                               allow_nan=False).encode("utf-8")).hexdigest()
        result.update(mesh_sha256=target_hash, native_report_sha256=report_hash,
                      native_report_hash_encoding="sorted_compact_UTF8_JSON", saved_native_checks=_saved_report(report, report_kind))
        remaining = limits.max_seconds - (time.perf_counter() - started)
        _require(remaining > 0, "time_budget", "Audit deadline exceeded before snapshot")
        with inspect_mesh(archive, max_cells=limits.max_cells, max_nodes=limits.max_nodes,
                          max_bytes=limits.max_bytes, max_work=limits.max_work, max_seconds=remaining) as inspected:
            _require(inspected.mesh_sha256 == target_hash, "archive_hash", "Immutable inspection snapshot hash differs")
            _require(inspected.validation.get("archive_structure_checked") is True,
                     "archive_structure_unchecked", "Inspection did not complete archive structure checks")
            inspected.max_seconds = limits.max_seconds - (time.perf_counter() - started)
            _require(inspected.max_seconds > 0, "time_budget", "Audit deadline exceeded before summary")
            summary = inspected.summary(spatial_bins=1)
            ownership = check_summary(summary, report, expected, report_kind)
            result.update(node_count=summary["node_count"], cell_count=summary["cell_count"],
                cell_counts=summary["cell_counts"], quality=summary["quality"],
                orientation_counts=summary["orientation_counts"], range_limited_cells=summary["range_limited_cells"],
                metric_conventions=summary.get("metric_conventions"), body_ownership=ownership,
                archive_structure=inspected.validation,
                signed_volume_check={"basis": "positive exact orientation and no dimensional metric range loss",
                                     "direct_volume_pagination_performed": False},
                achieved_tet_sizes=summary.get("sizing", {}).get("achieved", {}).get("tet4"),
                bounds=summary.get("bounds"), archive_screen_performed=True)
            if limits.detail_limit:
                inspected.max_seconds = limits.max_seconds - (time.perf_counter() - started)
                _require(inspected.max_seconds > 0, "time_budget", "Audit deadline exceeded before quality query")
                page = inspected.query(quality_below=.1, metric="mean_ratio", sort_by="mean_ratio",
                                       limit=limits.detail_limit, associations=False)
                result.update(below_mean_ratio_0p1=page["total"], worst_cells=page["cells"],
                              worst_cells_truncated=page["next_offset"] is not None)
            _deadline(started, limits.max_seconds)
        _require(_digest(archive, limits.max_bytes // 2, started, limits.max_seconds) == target_hash,
                 "archive_hash", "Original archive changed during the independent screen")
        result.update(status="passed", validation_scope="archive structure, positive Tet4 metrics and recorded body tags")
    except (AuditFailure, ValueError, TypeError, KeyError, OverflowError, OSError) as error:
        reason = getattr(error, "reason", "inspection_failed")
        result.update(status="resource_limit" if "budget" in reason else "failed", reason=reason,
                      detail=str(error)[:2048])
    finally:
        result["independent_audit_seconds"] = time.perf_counter() - started
    return result


def atomic_json(output, data, maximum):
    _require(not os.path.lexists(output), "output_exists", "Audit output must be a new file")
    payload = (json.dumps(data, indent=2, allow_nan=False) + "\n").encode("utf-8")
    _require(len(payload) <= maximum, "output_budget", "Audit output exceeds its byte allowance")
    _require(output.parent.is_dir(), "output_parent", "Audit output parent must exist")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=output.parent, prefix=".paacs-corpus-audit-", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        # Same-directory hard link publishes the complete file exclusively.
        # Unlike replace(), a racing pre-existing destination stays untouched.
        if os.name == "nt":
            # Windows rename is atomic on one filesystem and refuses an
            # existing destination. Hard links can be unavailable in sandboxes.
            os.rename(temporary, output)
        else:
            os.link(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _progress(stream, event, source, repeat, record, started, result=None):
    if stream is None:
        return
    path = source["path"]
    message = {"schema": PROGRESS_SCHEMA, "event": event, "record": record,
               "path": path[:1024], "path_truncated": len(path) > 1024,
               "source_sha256": source["sha256"], "repeat": repeat,
               "independent_audit_seconds": time.perf_counter() - started}
    if result is not None:
        message["status"] = result["status"]
        if result.get("reason") is not None:
            message["reason"] = str(result["reason"])[:128]
    # At most two short lines per retained input record; never include cells.
    print(json.dumps(message, allow_nan=False), file=stream, flush=True)


def run_audit(run, output, *, plan=None, limits=Limits(), inspect_mesh=None, progress=None):
    limits.validate()
    started = time.perf_counter()
    run, output = Path(run).resolve(strict=True), Path(output).absolute()
    _require(not os.path.lexists(output), "output_exists", "Audit output must be a new file")
    metadata_cap = min(limits.max_metadata_bytes, limits.max_bytes // 32)
    manifest, manifest_hash = _json_file(run / "manifest.json", metadata_cap)
    _require(manifest.get("schema") == "paacs-cad-corpus-v1", "manifest_schema", "Expected a corpus manifest")
    files = manifest.get("files")
    _require(isinstance(files, list) and 0 < len(files) <= limits.max_records,
             "manifest_files", "Manifest files are missing or exceed the record allowance")
    repeats = _int(manifest.get("settings", {}).get("repeats"), "repeats", 1)
    _require(len(files) * repeats <= limits.max_records, "record_budget", "Planned attempts exceed max_records")
    sources = {}
    for source in files:
        _require(isinstance(source, dict) and isinstance(source.get("path"), str), "manifest_files", "Invalid source row")
        _hash(source.get("sha256"), "source sha256")
        _require(source["path"] not in sources, "manifest_files", "Duplicate source path")
        sources[source["path"]] = source
    expected, plan_hash = {}, None
    if plan is not None:
        declared, plan_hash = _json_file(Path(plan), metadata_cap)
        _require(isinstance(declared.get("files"), list), "plan_files", "Fixed plan files are missing")
        for source in declared["files"]:
            _require(isinstance(source, dict) and source.get("path") not in expected, "plan_files", "Invalid or duplicate plan source")
            expected[source["path"]] = source
        _require(set(expected) == set(sources), "plan_files", "Fixed plan and run source sets differ")
        for key in sources:
            _require(expected[key].get("sha256") == sources[key]["sha256"], "source_hash", "Plan source hashes differ")
    if inspect_mesh is None:
        import paacs_mesh as pm
        inspect_mesh = pm.inspect_mesh
    results_path = run / "results.jsonl"
    _regular(results_path)
    _require(results_path.stat().st_size <= limits.max_metadata_bytes, "metadata_budget", "Results exceed input byte allowance")
    rows, duplicate_rows, consumed, records = {}, [], 0, 0
    with results_path.open("rb") as stream:
        while True:
            payload = stream.readline(metadata_cap + 1)
            if not payload:
                break
            consumed += len(payload)
            records += 1
            _require(len(payload) <= metadata_cap and consumed <= limits.max_metadata_bytes,
                     "metadata_budget", "Results JSONL exceeds its parsing allowance")
            _require(records <= limits.max_records, "record_budget", "Results exceed max_records")
            record = _decode(payload)
            key = (record.get("path"), record.get("repeat"))
            _require(key[0] in sources and isinstance(key[1], int) and not isinstance(key[1], bool)
                     and 0 <= key[1] < repeats, "result_identity", "Unexpected source or repeat in results")
            attempt_started = time.perf_counter()
            _progress(progress, "start", sources[key[0]], key[1], records, attempt_started)
            if key in rows:
                duplicate = {"path": key[0], "repeat": key[1], "status": "failed", "reason": "duplicate_result"}
                duplicate_rows.append(duplicate)
                _progress(progress, "end", sources[key[0]], key[1], records, attempt_started, duplicate)
                continue
            if time.perf_counter() - started >= limits.max_total_seconds:
                rows[key] = {"path": key[0], "repeat": key[1], "status": "resource_limit", "reason": "total_time_budget"}
            else:
                rows[key] = audit_record(run, sources[key[0]], record, expected.get(key[0]), limits, inspect_mesh)
            _progress(progress, "end", sources[key[0]], key[1], records, attempt_started, rows[key])
    ordered = []
    for source in files:
        for repeat in range(repeats):
            ordered.append(rows.get((source["path"], repeat), {"path": source["path"], "repeat": repeat,
                "status": "failed", "reason": "missing_result", "archive_screen_performed": False}))
    ordered.extend(duplicate_rows)
    counts = dict(Counter(row["status"] for row in ordered))
    receipt = {"schema": SCHEMA, "status": "complete" if not any(row["status"] in ("failed", "resource_limit") for row in ordered) else "incomplete",
        "manifest_sha256": manifest_hash, "plan_sha256": plan_hash, "planned_attempts": len(files) * repeats,
        "independent_archive_passes": counts.get("passed", 0),
        "all_planned_meshes_passed": len(ordered) == len(files) * repeats and counts.get("passed", 0) == len(ordered),
        "independent_native_reaudit_performed": False, "global_intersection_audit_performed": False,
        "saved_CAD_validation_basis": SAVED_BASIS, "source_identity_basis": "runner input_sha256 and input_unchanged marker matched to manifest and optional fixed plan",
        "limits": vars(limits), "limit_scope": "max_bytes/max_work bound each inspection snapshot/operation; metadata, record and output limits are separate; no process RSS cap; deadlines are cooperative",
        "timing_scope": "independent audit only; excluded from saved meshing deadline and timings",
        "counts": counts, "cases": ordered, "independent_total_seconds": time.perf_counter() - started,
        "contains_local_mesh_coordinates": bool(limits.detail_limit)}
    atomic_json(output, receipt, limits.max_output_bytes)
    return receipt


def audit_trial(case_directory, output, *, expected_original_solids=None, expected_body_tags=None,
                limits=Limits(), inspect_mesh=None):
    """Audit an explicitly selected construction trial and its saved Python report.

    The report is checked in its original format, without translating it into a
    NativeCLI report or claiming a new native CAD audit. Both recorded artifact
    hashes and the source fingerprint must agree before inspection starts.
    """
    limits.validate()
    started = time.perf_counter()
    root, output = Path(case_directory).resolve(strict=True), Path(output).absolute()
    _require(not os.path.lexists(output), "output_exists", "Audit output must be a new file")
    maximum = min(limits.max_metadata_bytes, limits.max_bytes // 32)
    case, case_hash = _json_file(root / "case.json", maximum)
    _require(case.get("schema") == "paacs-construction-budget-case-v1", "trial_schema", "Expected a construction trial record")
    _require(case.get("native_status") in ("succeeded", "refused"), "trial_native_status", "Trial lacks a completed native outcome")
    source_hash = _hash(case.get("source_sha256"), "trial source_sha256")
    report_path = _inside(root, case.get("native_report_file"))
    report_hash = _hash(case.get("native_report_sha256"), "native_report_sha256")
    try:
        report, actual_report_hash = _json_file(report_path, maximum)
    except FileNotFoundError as error:
        raise AuditFailure("missing_native_report", "Recorded Python report file is missing") from error
    _require(actual_report_hash == report_hash, "report_hash", "Recorded Python report file changed")
    _require(report.get("source_sha256") == source_hash, "source_hash", "Trial and Python report source fingerprints differ")
    source = {"path": PureWindowsPath(str(case.get("source", "trial-source"))).name, "sha256": source_hash}
    expected = {"sha256": source_hash, "expected_original_solids": expected_original_solids,
                "expected_body_tags": expected_body_tags}
    repeat = _int(case.get("repetition", 1), "repetition", 1) - 1
    if case["native_status"] == "succeeded":
        check_native(report, "python")
        if inspect_mesh is None:
            import paacs_mesh as pm
            inspect_mesh = pm.inspect_mesh
        row = {"outcome": "succeeded", "input_sha256": source_hash, "repeat": repeat,
               "output": "mesh.paacsmesh", "output_exists": True,
               "output_sha256": case.get("mesh_sha256"), "result": report,
               "wall_seconds": case.get("call_seconds")}
        result = audit_record(root, source, row, expected, limits, inspect_mesh, report_kind="python")
    else:
        _require(report.get("mesh_published") is False and report.get("status") != "succeeded"
                 and not os.path.lexists(root / "mesh.paacsmesh") and not case.get("mesh_sha256"),
                 "publication_violation", "Refused trial retained a published output")
        result = {"path": source["path"], "source_sha256": source_hash, "repeat": repeat,
                  "status": "not_published", "reason": report.get("reason"), "native_outcome": "refused",
                  "mesh_wall_seconds": case.get("call_seconds"), "saved_native_checks": _saved_report(report, "python"),
                  "archive_screen_performed": False, "independent_native_reaudit_performed": False,
                  "global_intersection_audit_performed": False,
                  "independent_audit_seconds": time.perf_counter() - started}
    result.update(original_worker_status=case.get("status"), trial_case_sha256=case_hash,
                  native_report_file_sha256=report_hash,
                  source_identity_basis="matching trial/Python-report fingerprints; CAD source was not imported")
    _require(_digest(report_path, maximum, started, limits.max_total_seconds) == report_hash,
             "report_hash", "Original Python report changed during the independent screen")
    receipt = {"schema": SCHEMA, "mode": "construction_trial", "status": "complete" if result["status"] in ("passed", "not_published") else "incomplete",
        "planned_attempts": 1, "independent_archive_passes": int(result["status"] == "passed"),
        "all_planned_meshes_passed": result["status"] == "passed", "trial_case_sha256": case_hash,
        "native_report_file_sha256": report_hash, "independent_native_reaudit_performed": False,
        "global_intersection_audit_performed": False,
        "saved_CAD_validation_basis": "saved_native_python_report_not_recertified",
        "limits": vars(limits), "limit_scope": "inspection limits per snapshot/operation; no process RSS cap; deadlines are cooperative",
        "timing_scope": "independent archive audit only; excluded from saved mesh timings",
        "counts": {result["status"]: 1}, "cases": [result],
        "independent_total_seconds": time.perf_counter() - started,
        "contains_local_mesh_coordinates": bool(limits.detail_limit)}
    atomic_json(output, receipt, limits.max_output_bytes)
    return receipt


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path, nargs="?")
    parser.add_argument("--trial-case", type=Path, help="explicit construction trial case directory instead of a corpus run")
    parser.add_argument("--expected-original-solids", type=int, help="known original body count for an explicit trial case")
    parser.add_argument("--expected-body-tags", type=int, nargs="+", help="known original body IDs for an explicit trial case")
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--progress", action="store_true",
                        help="flush bounded corpus attempt start/end JSONL to stderr")
    parser.add_argument("--development-native", type=Path)
    parser.add_argument("--occt-sdk", type=Path)
    defaults = Limits()
    for name, default in vars(defaults).items():
        parser.add_argument("--" + name.replace("_", "-"), type=float if "seconds" in name else int, default=default)
    args = parser.parse_args(argv)
    try:
        _require((args.run_directory is None) != (args.trial_case is None), "audit_mode",
                 "Choose exactly one corpus run or --trial-case")
        if args.development_native is not None:
            _require(args.occt_sdk is not None, "development_loader", "Development loading requires --occt-sdk")
            spec = importlib.util.spec_from_file_location("resource_audit_loader", Path(__file__).with_name("construction_budget_trials.py"))
            loader = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(loader)
            loader.load_development(args)
        limits = Limits(**{name: getattr(args, name) for name in vars(defaults)})
        if args.trial_case is not None:
            _require(args.plan is None, "audit_mode", "--plan applies to a corpus run, not an explicit trial")
            _require(not args.progress, "audit_mode", "--progress applies to corpus attempts, not an explicit trial")
            receipt = audit_trial(args.trial_case, args.output, expected_original_solids=args.expected_original_solids,
                                  expected_body_tags=args.expected_body_tags, limits=limits)
        else:
            _require(args.expected_original_solids is None and args.expected_body_tags is None,
                     "audit_mode", "Trial body expectations require --trial-case")
            receipt = run_audit(args.run_directory, args.output, plan=args.plan, limits=limits,
                                progress=sys.stderr if args.progress else None)
        print(json.dumps({"schema": SCHEMA, "status": receipt["status"], "counts": receipt["counts"],
                          "output": str(args.output.absolute()), "independent_total_seconds": receipt["independent_total_seconds"]}, allow_nan=False))
        return 0 if receipt["status"] == "complete" else 1
    except (AuditFailure, ValueError, OSError, KeyError, TypeError, ImportError) as error:
        print(json.dumps({"schema": SCHEMA, "status": "refused", "reason": getattr(error, "reason", "invalid_request"),
                          "detail": str(error)[:2048]}, allow_nan=False))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
