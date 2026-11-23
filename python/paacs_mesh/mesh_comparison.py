"""Bounded comparisons of source-bound inspected meshes; no cell correspondence."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from array import array
import json
import math
import time

from ._cad_common import MAX_JSON_BYTES, integer, json_text, number
from .mesh_archive import CELL_NAMES, InspectionError
from .mesh_inspection import METRICS, MeshInspection, _centroid, distribution
from .quality_metrics import METRIC_CONVENTIONS


SIZE_METRICS = ("edge_min", "edge_median", "edge_max", "edge_rms", "volume_equivalent_edge")
CAD_MEASUREMENTS = ("cad_volume", "mesh_volume", "relative_volume_error",
                    "maximum_boundary_node_error", "maximum_boundary_sample_error",
                    "surface_chord_error_limit", "centroid_outside_distance_limit")
_VALID = {"valid", "valid_with_warnings"}
_LOW_QUALITY = .2


class _Budget:
    def __init__(self, first, second, work, seconds, bins):
        self.max_work = integer(work, "max_work")
        self.max_seconds = number(seconds, "max_seconds", minimum=0)
        if self.max_seconds <= 0:
            raise ValueError("max_seconds must be positive")
        self.started = time.perf_counter()
        self.work = 0
        # Additional comparison-owned storage; existing snapshots are owned by
        # their callers. Respect the smaller admitted input memory allowance.
        self.max_bytes = min(first.max_bytes, second.max_bytes)
        self.retained = 65536 + 4 * MAX_JSON_BYTES + 4096 * bins ** 3
        self.retained += 2 * (first._archive.metadata_bytes + second._archive.metadata_bytes)
        self.peak = self.retained
        self._check_memory(self.retained)
        self.check()

    def _check_memory(self, total):
        if total > self.max_bytes:
            raise InspectionError("memory_budget", "Comparison scratch exceeds the input memory allowance",
                                  details={"scratch_reservation_bytes": total,
                                           "max_scratch_bytes": self.max_bytes})

    def reserve(self, size):
        self._check_memory(self.retained + size)
        self.retained += size
        self.peak = max(self.peak, self.retained)

    def temporary(self, size):
        self._check_memory(self.retained + size)
        self.peak = max(self.peak, self.retained + size)

    def take(self, amount=1):
        self.work += amount
        if self.work > self.max_work:
            raise InspectionError("work_budget", "Comparison exceeds max_work",
                                  details={"work": self.work, "max_work": self.max_work})
        self.check()

    def check(self):
        if time.perf_counter() - self.started > self.max_seconds:
            raise InspectionError("time_budget", "Comparison exceeds cooperative max_seconds")


def _difference(first, second):
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) for value in (first, second)):
        return None
    try:
        if not all(math.isfinite(value) for value in (first, second)):
            return None
        value = second - first
        return value if math.isfinite(value) else None
    except OverflowError:
        return None


def _distribution_delta(first, second):
    return {key: _difference(first.get(key), second.get(key))
            for key in ("count", "minimum", "p01", "p05", "median", "p95", "maximum", "mean")}


def _common_bounds(first, second, budget):
    result = None
    for inspected in (first, second):
        for node_id in range(inspected.node_count):
            budget.take(3)
            point = inspected._archive.node(node_id)[0]
            if result is None:
                result = [point[:], point[:]]
            else:
                for axis in range(3):
                    result[0][axis] = min(result[0][axis], point[axis])
                    result[1][axis] = max(result[1][axis], point[axis])
    return result


def _bin_index(point, extent, bins):
    index = []
    for axis in range(3):
        low, high = extent[0][axis], extent[1][axis]
        scale = max(abs(low), abs(high)) or 1.
        fraction = ((point[axis] / scale - low / scale) / (high / scale - low / scale)) if high > low else 0.
        index.append(max(0, min(bins - 1, int(fraction * bins))))
    return tuple(index)


def _interpolate(low, high, fraction):
    if fraction == 0:
        return low
    if fraction == 1:
        return high
    scale = max(abs(low), abs(high)) or 1.
    first, second = low / scale, high / scale
    value = first * (1 - fraction) + second * fraction
    return max(first, min(second, value)) * scale


def _bin_bounds(index, extent, bins):
    return [[_interpolate(extent[0][axis], extent[1][axis], (index[axis] + end) / bins)
             for axis in range(3)] for end in (0, 1)]


def _sizing(report, achieved, faces):
    sizing = report.get("sizing") or {}
    if not isinstance(sizing, dict):
        raise InspectionError("comparison_report", "Native sizing metadata must be an object")
    face_controls, boxes = sizing.get("face_sizes", []), sizing.get("volume_controls", [])
    if (not isinstance(face_controls, list) or not isinstance(boxes, list)
            or any(not isinstance(item, dict) for item in face_controls + boxes)):
        raise InspectionError("comparison_report", "Native sizing controls must be object lists")
    return {
        "requested": {"global_size": sizing.get("global_size", report.get("requested_size")),
            "min_size": sizing.get("requested_min_size"), "growth": sizing.get("growth"),
            "native_sizing": sizing.get("requested_native_sizing"),
            "face_sizes": [{"tag": item.get("tag"), "size": item.get("size"),
                            "growth": item.get("requested_growth")} for item in face_controls],
            "volume_controls": [{key: item.get(key) for key in
                ("kind", "lower", "upper", "size", "requested_growth")} for item in boxes]},
        "effective": {"minimum_target": sizing.get("min_size"), "global_target": sizing.get("global_size"),
            "face_targets": [{"tag": item.get("tag"), "size": item.get("effective_size"),
                              "growth": item.get("growth")} for item in face_controls],
            "volume_targets": [{key: item.get(key) for key in
                ("kind", "lower", "upper", "effective_size", "growth")} for item in boxes],
            "resolved_routing": report.get("effective_policy"), "sizing_completion_certificate": "unavailable"},
        "achieved": {"tet4": achieved, "boundary_faces": faces},
        "achieved_growth": {"available": False, "reason": "no_neighbor_growth_audit_performed"},
        "boundary_population": "triangle-edge incidences; shared edges appear in each incident triangle",
        "tet4_population": "per-cell edge and regular-volume-equivalent-edge distributions",
        "volume_equivalent_edge_caveat": "Volume proxy depends on shape; it does not certify maximum edge or accuracy",
        "target_is_not_a_guarantee_of_every_edge_length": True,
    }


def _profile(inspected, report, extent, bins, budget):
    cells = inspected.cell_count
    tetrahedra = inspected._archive.type_counts[5]
    face_lengths, spatial, body_counts = {}, {}, {}
    unknown_bodies = 0
    for cell_id in range(cells):
        budget.take(8)
        kind, nodes, classification = inspected._archive.cell(cell_id)
        if kind == 2 and classification["dimension"] == 2 and classification["tag"] > 0:
            budget.take(12)
            tag = str(classification["tag"])
            if tag not in face_lengths:
                budget.reserve(512)
                face_lengths[tag] = array("d")
            budget.reserve(48)
            points = [inspected._archive.node(node)[0] for node in nodes]
            face_lengths[tag].extend(value if math.isfinite(value) else math.nan for value in (
                math.dist(points[a], points[b]) for a, b in ((0, 1), (1, 2), (2, 0))))
        if kind != 5:
            continue
        budget.take(24)
        center = _centroid([inspected._archive.node(node)[0] for node in nodes])
        index = _bin_index(center, extent, bins)
        row = spatial.get(index)
        if row is None:
            row = {"count": 0, "minimum_mean_ratio": None, "below_0_2": 0, "unavailable_mean_ratio": 0,
                   "body_tags": set(), "unclassified_body_count": 0}
            spatial[index] = row
        value = inspected._metric(cell_id, "mean_ratio")
        row["count"] += 1
        if value is not None:
            row["minimum_mean_ratio"] = value if row["minimum_mean_ratio"] is None else min(row["minimum_mean_ratio"], value)
            row["below_0_2"] += value < _LOW_QUALITY
        else:
            row["unavailable_mean_ratio"] += 1
        if classification["dimension"] == 3 and classification["tag"] > 0:
            tag = classification["tag"]
            if tag not in row["body_tags"]:
                budget.reserve(256)
                row["body_tags"].add(tag)
            if tag not in body_counts:
                budget.reserve(256)
            body_counts[tag] = body_counts.get(tag, 0) + 1
        else:
            row["unclassified_body_count"] += 1
            unknown_bodies += 1
    quality, achieved = {}, {}
    budget.temporary(128 * max(tetrahedra, 1))
    sort_work = tetrahedra * max(1, (tetrahedra - 1).bit_length())
    for name in METRICS + SIZE_METRICS:
        budget.take(cells + sort_work)
        measured = distribution(inspected._metric(index, name) for index in range(cells)
                                if inspected._orientation[index] != 2)
        budget.check()
        (quality if name in METRICS else achieved)[name] = measured
    boundary_sizes = {}
    for tag, values in face_lengths.items():
        budget.temporary(128 * max(len(values), 1))
        budget.take(len(values) * (1 + max(1, (len(values) - 1).bit_length())))
        boundary_sizes[tag] = distribution(values)
        budget.check()
    return {"node_count": inspected.node_count, "cell_count": cells, "tetrahedron_count": tetrahedra,
            "cell_counts": dict(zip(CELL_NAMES, inspected._archive.type_counts)),
            "quality": quality, "sizing": _sizing(report, achieved, boundary_sizes),
            "unavailable_tet4_measurements": {name: tetrahedra - value["count"]
                                              for name, value in {**quality, **achieved}.items()},
            "assembly": {key: report.get(key) for key in
                         ("mesh_scope", "interbody_nonoverlap_certified", "interferences")},
            "body_tet4_counts": {str(tag): count for tag, count in sorted(body_counts.items())},
            "unclassified_tet4_count": unknown_bodies}, spatial


def _region_record(value):
    value = value or {}
    return {"count": value.get("count", 0), "minimum_mean_ratio": value.get("minimum_mean_ratio"),
            "below_0_2": value.get("below_0_2", 0), "body_tags": sorted(value.get("body_tags", ())),
            "unavailable_mean_ratio": value.get("unavailable_mean_ratio", 0),
            "unclassified_body_count": value.get("unclassified_body_count", 0)}


def _cad_errors(first, second, budget):
    maps = []
    for report in (first, second):
        found = {}
        for body in report.get("bodies", []):
            budget.take(8)
            tag = body.get("solid_tag")
            if tag is None:
                continue
            tag = integer(tag, "solid_tag")
            if tag in found:
                raise InspectionError("comparison_report", "Native report duplicates a CAD body tag")
            budget.reserve(2048)
            found[tag] = {key: body.get(key) for key in CAD_MEASUREMENTS}
        maps.append(found)
    kernels = [report.get("kernel_version") for report in (first, second)]
    paired = isinstance(kernels[0], str) and bool(kernels[0]) and kernels[0] == kernels[1]
    common = {"basis": "reported_native_body_measurements_not_independent_CAD_reaudit",
              "body_comparison_available": paired, "kernel_versions": kernels,
              "length_units": "mm", "volume_units": "mm3"}
    if not paired:
        return {**common, "body_identity_basis": "unavailable_without_matching_known_source_and_kernel",
                "reason": "CAD body tags cannot be paired across missing or different kernel identities",
                "first_body_measurements": [{"body_tag": tag, **value} for tag, value in sorted(maps[0].items())],
                "second_body_measurements": [{"body_tag": tag, **value} for tag, value in sorted(maps[1].items())],
                "bodies": []}
    bodies = []
    for tag in sorted(maps[0].keys() | maps[1].keys()):
        budget.take(8)
        before, after = maps[0].get(tag), maps[1].get(tag)
        bodies.append({"body_tag": tag, "first": before, "second": after,
                       "delta": {key: _difference((before or {}).get(key), (after or {}).get(key))
                                 for key in CAD_MEASUREMENTS}})
    return {**common, "body_identity_basis": "reported_positive_CAD_body_tags_with_same_source_hash_and_kernel",
            "bodies": bodies}


def _timing(report):
    seconds = report.get("seconds")
    if seconds is not None and not isinstance(seconds, dict):
        raise InspectionError("comparison_report", "Native seconds metadata must be an object")
    return {"reported_seconds": seconds, "meshing_scope": "final_native_attempt_only",
            "import_scope": "reported_model_import", "workflow_total_seconds": None,
            "failed_retry_seconds": None, "retry_count": None,
            "retry_timing_available": False, "reason": "no_source_and_final_report_paired_workflow_was_supplied"}


def compare_meshes(first, second, *, spatial_bins=4, max_work=20_000_000, max_seconds=60):
    """Compare reported-valid snapshots of the same CAD source and units.

    Spatial bins share one union bounding box. Cell IDs are never paired across
    snapshots. Native timings cover the final attempt only; this interface does
    not accept or guess workflow reports and therefore cannot total failed
    retries. Comparisons do not establish solver accuracy or CAD validity.
    """
    if not isinstance(first, MeshInspection) or not isinstance(second, MeshInspection):
        raise TypeError("compare_meshes requires two MeshInspection objects")
    bins = integer(spatial_bins, "spatial_bins")
    if bins > 16:
        raise ValueError("spatial_bins must not exceed 16 per axis")
    budget = _Budget(first, second, max_work, max_seconds, bins)
    identities = [item.identity() for item in (first, second)]
    source = identities[0]["source_sha256"]
    if source is None or source != identities[1]["source_sha256"]:
        raise InspectionError("comparison_source", "Comparison requires the same known CAD source SHA-256")
    units = identities[0]["units"]
    if units is None or units != identities[1]["units"]:
        raise InspectionError("comparison_units", "Comparison requires the same known length units")
    budget.take(1 + (first._archive.metadata_bytes + second._archive.metadata_bytes) // 1024)
    validation = [item.validation for item in (first, second)]
    if any(item["reported_native_status"] not in _VALID for item in validation):
        raise InspectionError("comparison_validation", "Both snapshots require reported valid native checks",
                              details={"first": validation[0], "second": validation[1]})
    reports = [first.native_report, second.native_report]
    budget.check()
    extent = _common_bounds(first, second, budget)
    profiles, spatial = [], []
    for item, report in zip((first, second), reports):
        profile, regions = _profile(item, report, extent, bins, budget)
        profiles.append(profile)
        spatial.append(regions)
    regions, changed = [], []
    for index in sorted(spatial[0].keys() | spatial[1].keys()):
        budget.take(16)
        before, after = (_region_record(item.get(index)) for item in spatial)
        delta = {key: _difference(before[key], after[key])
                 for key in ("count", "minimum_mean_ratio", "below_0_2", "unavailable_mean_ratio", "unclassified_body_count")}
        row = {"index": list(index), "bounds": _bin_bounds(index, extent, bins),
               "first": before, "second": after, "delta": delta}
        regions.append(row)
        if before != after:
            changed.append(list(index))
    timing = [_timing(report) for report in reports]
    timing_keys = set(timing[0]["reported_seconds"] or {}) | set(timing[1]["reported_seconds"] or {})
    result = {"schema": "paacs-mesh-comparison-v1", "source_sha256": source, "units": units,
        "first": {"identity": identities[0], "validation": validation[0], **profiles[0]},
        "second": {"identity": identities[1], "validation": validation[1], **profiles[1]},
        "delta_convention": "second_minus_first; null means at least one measurement was unavailable or out of range",
        "counts_delta": {key: profiles[1][key] - profiles[0][key]
                         for key in ("node_count", "cell_count", "tetrahedron_count")},
        "cell_counts_delta": {name: profiles[1]["cell_counts"][name] - profiles[0]["cell_counts"][name]
                              for name in CELL_NAMES},
        "quality_delta": {name: _distribution_delta(profiles[0]["quality"][name], profiles[1]["quality"][name])
                          for name in METRICS},
        "achieved_tet4_size_delta": {name: _distribution_delta(profiles[0]["sizing"]["achieved"]["tet4"][name],
                                         profiles[1]["sizing"]["achieved"]["tet4"][name]) for name in SIZE_METRICS},
        "cad_error": _cad_errors(*reports, budget),
        "spatial_comparison": {"kind": "common_uniform_axis_aligned_centroid_bins", "bounds": extent,
            "bins_per_axis": bins, "low_mean_ratio_threshold": _LOW_QUALITY,
            "assignment": "centroid; half-open bins with the global upper boundary in the final bin",
            "body_identity_basis": "positive_dimension_3_archive_classification_with_known_CAD_source",
            "regions": regions, "affected_region_indices": changed},
        "cell_correspondence": "not_assumed; cell IDs belong only to their own mesh hash",
        "timings": {"first": timing[0], "second": timing[1], "reported_seconds_delta": {
            key: _difference((timing[0]["reported_seconds"] or {}).get(key),
                             (timing[1]["reported_seconds"] or {}).get(key)) for key in sorted(timing_keys)},
            "end_to_end_runtime_comparison_available": False},
        "limits": {"max_work": budget.max_work, "work_used": budget.work,
            "max_seconds": budget.max_seconds, "cooperative_deadline": True,
            "max_scratch_bytes": budget.max_bytes, "scratch_reservation_bytes": budget.peak,
            "scratch_scope": "additional_comparison_storage_excluding_caller_owned_snapshots",
            "max_output_bytes": MAX_JSON_BYTES},
        "metric_conventions": METRIC_CONVENTIONS,
        "measurement_population": "finite Tet4 measurements; unavailable values are counted separately, not replaced by zero",
        "limitations": ["reported_native_validity_is_not_an_independent_reaudit",
                        "finer_elements_do_not_guarantee_better_shape_or_solver_accuracy",
                        "native_timing_does_not_total_failed_workflow_retries",
                        "Tet4_quality_only; other_cell_families_are_counted_but_not_measured",
                        "body_tags_are_reported_source_local_identities_not_independently_rebound"]}
    try:
        encoded = json_text(result)
    except ValueError as error:
        raise InspectionError("output_budget", "Comparison cannot fit finite JSON within the 8 MiB output budget") from error
    budget.take(3 * max(1, len(encoded.encode("utf-8")) // 64))
    result["limits"]["work_used"] = budget.work
    # The serialization charge covers encoding, accounting refresh and copying.
    # The final deadline check also applies to this small final population.
    result = json.loads(json_text(result))
    budget.check()
    return result
