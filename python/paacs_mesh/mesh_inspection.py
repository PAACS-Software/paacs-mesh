"""Opt-in Tet4 measurements and paginated access to immutable mesh snapshots."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from array import array
from collections import Counter
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
import time

from ._cad_common import bounds, integer, json_copy, json_text, number, sequence
from .mesh_archive import ArchiveView, CELL_NAMES, InspectionError, atomic_stream, file_hash
from .quality_metrics import METRIC_CONVENTIONS, tet4_metrics


METRICS = ("mean_ratio", "minimum_dihedral_deg", "scaled_jacobian")
MEASUREMENTS = METRICS + ("signed_volume", "edge_min", "edge_median", "edge_max", "edge_rms", "volume_equivalent_edge")
ORIENTATIONS = {"inverted": -1, "degenerate": 0, "positive": 1, "unsupported": 2}
TET_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3))


def distribution(values, *, edges=None):
    """Exact sorted linear-rank quantiles; caller admits sorting storage first."""
    finite, unavailable = [], 0
    for value in values:
        if value is None or not math.isfinite(value):
            unavailable += 1
        else:
            finite.append(value)
    values = finite
    values.sort()
    result = {"count": len(values), "unavailable": unavailable, "minimum": None, "p01": None, "p05": None,
              "median": None, "p95": None, "maximum": None, "mean": None}
    if values:
        def quantile(p):
            position = (len(values) - 1) * p
            low = int(position)
            fraction = position - low
            first, second = values[low], values[min(low + 1, len(values) - 1)]
            scale = max(abs(first), abs(second)) or 1.
            local = (first / scale) * (1 - fraction) + (second / scale) * fraction
            return max(first / scale, min(second / scale, local)) * scale
        try:
            mean = math.fsum(values) / len(values)
        except OverflowError:
            scale = max(abs(value) for value in values) or 1.
            local = math.fsum(value / scale for value in values) / len(values)
            mean = max(-1., min(1., local)) * scale
        result.update(minimum=values[0], p01=quantile(.01), p05=quantile(.05), median=quantile(.5),
                      p95=quantile(.95), maximum=values[-1], mean=mean)
    if edges is not None:
        counts = [0] * (len(edges) + 1)
        for value in values:
            counts[sum(value >= edge for edge in edges)] += 1
        result["histogram"] = {"edges": list(edges), "counts": counts,
                               "rule": "left closed; values equal to an edge enter the higher bin"}
    return result


def _centroid(points):
    result = []
    for axis in range(3):
        values = [point[axis] for point in points]
        try:
            value = math.fsum(values) / len(values)
        except OverflowError:
            scale = max(abs(v) for v in values) or 1.
            value = (math.fsum(v / scale for v in values) / len(values)) * scale
        result.append(min(max(value, min(values)), max(values)))
    return result


def _box(points):
    return [[min(p[axis] for p in points) for axis in range(3)],
            [max(p[axis] for p in points) for axis in range(3)]]


class MeshInspection:
    """Owned immutable snapshot, bounded compact metrics, stable native cell IDs.

    Inspection measures element shape. It never retries meshing, changes geometry
    or replaces the native audit. Bare archive v1 has no stored units, CAD hash or
    validation report; these stay unavailable unless paired in a saved bundle.
    """

    def __init__(self, source, *, max_bytes=256 * 1024 * 1024, max_cells=1_000_000,
                 max_nodes=1_000_000, max_work=20_000_000, max_seconds=60):
        self._construction_started = time.perf_counter()
        self.max_seconds = number(max_seconds, "max_seconds", minimum=0)
        if self.max_seconds <= 0:
            raise ValueError("max_seconds must be positive")
        self._archive = ArchiveView(source, max_bytes=max_bytes, max_cells=max_cells,
                                    max_nodes=max_nodes, max_work=max_work, max_seconds=self.max_seconds)
        self.max_bytes, self.max_work = self._archive.max_bytes, self._archive.max_work
        self._columns = {}
        self._orientation = array("b")
        self._range_limited = array("b")
        self._boundary = None
        self.storage_bytes = self._archive.storage_bytes
        try:
            self._analyze()
        except BaseException:
            self.close()
            raise

    @property
    def mesh_sha256(self):
        return self._archive.mesh_sha256

    @property
    def source_sha256(self):
        return self._archive.source_sha256

    @property
    def units(self):
        return self._archive.units

    @property
    def cell_count(self):
        return self._archive.cell_count

    @property
    def node_count(self):
        return self._archive.node_count

    @property
    def native_report(self):
        return json_copy(self._archive.report) if self._archive.report is not None else None

    @property
    def validation(self):
        report = self._archive.report
        status = "unavailable"
        if report is not None:
            from .diagnostics import diagnose
            try:
                status = diagnose(report)["status"]
            except (ValueError, TypeError):
                status = "invalid_report"
        return {"basis": self._archive.report_basis, "reported_native_status": status,
                "independent_native_reaudit_performed": False,
                "archive_structure_checked": True}

    def identity(self):
        return {"mesh_sha256": self.mesh_sha256, "source_sha256": self.source_sha256,
                "kernel_version": (self._archive.report or {}).get("kernel_version"),
                "units": self.units, "cell_id_scheme": "zero_based_native_archive_cell_index_v1",
                "node_id_scheme": "zero_based_native_node_index_v1"}

    def _reserve(self, scratch, operation):
        if scratch + self.storage_bytes > self.max_bytes:
            raise InspectionError("memory_budget", f"{operation} exceeds max_bytes",
                                  details={"owned_storage_bytes": self.storage_bytes,
                                           "scratch_reservation_bytes": scratch, "max_bytes": self.max_bytes})

    def _deadline(self, started, work):
        if work > self.max_work:
            raise InspectionError("work_budget", "Operation exceeds max_work")
        if time.perf_counter() - started > self.max_seconds:
            raise InspectionError("time_budget", "Operation exceeds its cooperative max_seconds")

    def _analyze(self):
        count = self.cell_count
        payload = count * (8 * len(MEASUREMENTS) + 2)
        self._reserve(payload + 4096, "Per-element measurements")
        self._columns = {name: array("d", [math.nan]) * count for name in MEASUREMENTS}
        self._orientation = array("b", [2]) * count
        self._range_limited = array("b", [0]) * count
        started = self._construction_started
        analysis_started = time.perf_counter()
        for cell_id in range(count):
            if cell_id % 128 == 0:
                self._deadline(started, cell_id + 1)
            kind, nodes, _ = self._archive.cell(cell_id)
            if kind != 5:
                continue
            measured = tet4_metrics([self._archive.node(n)[0] for n in nodes])
            for name in MEASUREMENTS:
                value = measured[name]
                self._columns[name][cell_id] = value if value is not None else math.nan
            self._orientation[cell_id] = ORIENTATIONS[measured["orientation"]]
            self._range_limited[cell_id] = measured["range_limited"]
        self._deadline(started, count)
        self.storage_bytes += payload + 4096
        self.analysis_seconds = time.perf_counter() - analysis_started

    def _metric(self, cell_id, name):
        value = self._columns[name][cell_id]
        return value if math.isfinite(value) else None

    def _boundary_index(self):
        if self._boundary is not None:
            return self._boundary
        # Tuple keys, small lists and Python integers are conservatively admitted
        # before construction; no full volume-face graph is needed for this map.
        self._reserve(512 * self._archive.type_counts[2], "CAD boundary association index")
        result = {}
        started = time.perf_counter()
        for cell_id in range(self.cell_count):
            if cell_id % 256 == 0:
                self._deadline(started, cell_id + 1)
            kind, nodes, classification = self._archive.cell(cell_id)
            if kind == 2 and classification["dimension"] == 2:
                result.setdefault(tuple(sorted(nodes)), []).append((cell_id, classification["dimension"], classification["tag"]))
        self._boundary = result
        self.storage_bytes += 512 * self._archive.type_counts[2]
        return result

    def _record(self, cell_id, *, associations=True):
        kind, nodes, classification = self._archive.cell(cell_id)
        points_and_classes = [self._archive.node(n) for n in nodes]
        points = [pair[0] for pair in points_and_classes]
        source, kernel = self.source_sha256, (self._archive.report or {}).get("kernel_version")
        guarded_cad = (isinstance(source, str) and len(source) == 64
                       and all(char in "0123456789abcdef" for char in source)
                       and isinstance(kernel, str) and bool(kernel))
        result = {"cell_id": cell_id, "cell_type": CELL_NAMES[kind], "node_ids": list(nodes),
                  "points": points, "node_classifications": [pair[1] for pair in points_and_classes],
                  "classification": classification,
                  "body_tag": classification["tag"] if classification["dimension"] == 3 else None,
                  "cad_body_identity_known": guarded_cad and classification["dimension"] == 3 and classification["tag"] > 0,
                  "centroid": _centroid(points), "bounds": _box(points),
                  "metrics": {name: self._metric(cell_id, name) for name in MEASUREMENTS},
                  "orientation": next(name for name, value in ORIENTATIONS.items() if value == self._orientation[cell_id]),
                  "range_limited": bool(self._range_limited[cell_id]),
                  "metric_support": "linear_tet4" if kind == 5 else "unsupported_cell_family"}
        faces = []
        if associations and kind == 5:
            boundary = self._boundary_index()
            for face in TET_FACES:
                for boundary_id, dimension, tag in boundary.get(tuple(sorted(nodes[i] for i in face)), ()):
                    faces.append({"boundary_cell_id": boundary_id, "dimension": dimension, "tag": tag,
                                  "evidence": "identical_boundary_triangle_node_ids"})
                    if len(faces) > 32:
                        raise InspectionError("association_budget", "Cell has more than 32 boundary associations")
        result["cad_boundary_associations"] = faces
        result["nearby_cad_surfaces"] = {"available": False,
                                          "reason": "no_exact_trimmed_CAD_projection_was_performed"}
        return result

    def query(self, *, cell_id=None, offset=0, limit=100, metric="mean_ratio", quality_below=None,
              region=None, spatial_mode="centroid", body_tags=None, cell_type="tet4",
              orientation=None, sort_by="cell_id", descending=False, associations=True,
              cad=None, cad_distance=None, cad_face_tags=None, max_cad_work=100000,
              max_cad_seconds=10, cad_include_faces=False):
        """Query all matching cells, with exact total and explicit resource refusal.

        Region modes are centroid inclusion or AABB overlap, not a claim that a
        cell intersects the exact CAD region. Sorting ties use the native cell ID.
        Offset pages refer to the filtered/sorted result, never renumbered cells.
        """
        offset, limit = integer(offset, "offset", minimum=0), integer(limit, "limit")
        if limit > 1000:
            raise ValueError("limit must not exceed 1000")
        if metric not in METRICS or sort_by not in ("cell_id",) + METRICS:
            raise ValueError("Unknown quality metric or sort key")
        if not isinstance(descending, bool) or not isinstance(associations, bool):
            raise ValueError("descending and associations must be bool")
        if spatial_mode not in ("centroid", "aabb_overlap"):
            raise ValueError("spatial_mode must be centroid or aabb_overlap")
        if cell_type is not None and cell_type not in CELL_NAMES:
            raise ValueError("Unknown cell_type")
        if orientation is not None and orientation not in ORIENTATIONS:
            raise ValueError("Unknown orientation")
        threshold = number(quality_below, "quality_below") if quality_below is not None else None
        cad_threshold = number(cad_distance, "cad_distance", minimum=0) if cad_distance is not None else None
        if cad is None and (cad_threshold is not None or cad_face_tags is not None or cad_include_faces):
            raise ValueError("CAD proximity options require the matching CadModel")
        box = bounds(region, "region") if region is not None else None
        bodies = None if body_tags is None else {integer(tag, "body_tag", minimum=0)
                                                for tag in sequence(body_tags, "body_tags", maximum=1000)}
        if cell_id is None:
            candidates = range(self.cell_count)
        else:
            cell_id = integer(cell_id, "cell_id", minimum=0)
            if cell_id >= self.cell_count:
                raise ValueError("cell_id is outside this mesh")
            candidates = (cell_id,)
            cell_type = None
        started = time.perf_counter()
        boundary_work = self.cell_count if associations and self._boundary is None else 0
        base_work = len(candidates) + boundary_work
        self._deadline(started, base_work)
        if associations:
            self._boundary_index()
        query_reservation = 224 * len(candidates) + 96 * 1024 * min(limit, len(candidates)) + 8192
        self._reserve(query_reservation, "Filtered ordering and output page")
        found = []
        for work, index in enumerate(candidates, 1):
            if work % 128 == 0:
                self._deadline(started, work)
            kind, nodes, classification = self._archive.cell(index)
            if cell_type is not None and CELL_NAMES[kind] != cell_type:
                continue
            if bodies is not None and (classification["dimension"] != 3 or classification["tag"] not in bodies):
                continue
            if orientation is not None and self._orientation[index] != ORIENTATIONS[orientation]:
                continue
            value = self._metric(index, metric)
            if threshold is not None and (value is None or value >= threshold):
                continue
            if box is not None:
                points = [self._archive.node(n)[0] for n in nodes]
                if spatial_mode == "centroid":
                    center = _centroid(points)
                    if not all(box[0][axis] <= center[axis] <= box[1][axis] for axis in range(3)):
                        continue
                else:
                    extent = _box(points)
                    if any(extent[1][axis] < box[0][axis] or extent[0][axis] > box[1][axis] for axis in range(3)):
                        continue
            found.append(index)
        if sort_by == "cell_id":
            if descending:
                found.reverse()
        else:
            # Missing values are always last; cell IDs resolve ties consistently.
            found.sort(key=lambda index: (self._metric(index, sort_by) is None,
                       (-1 if descending else 1) * (self._metric(index, sort_by) or 0), index))
        proximity = None
        proximity_rows = {}
        evaluations = 0
        if cad is not None:
            from .cad_proximity import query_facial_proximity
            cad_work = integer(max_cad_work, "max_cad_work")
            if cad_work > 100000:
                raise ValueError("max_cad_work must not exceed 100000 CAD distance evaluations")
            remaining_work = self.max_work - base_work
            if found and remaining_work <= 0:
                raise InspectionError("work_budget", "Query exhausted its work budget before CAD proximity")
            remaining_seconds = self.max_seconds - (time.perf_counter() - started)
            if remaining_seconds <= 0:
                raise InspectionError("time_budget", "Query exhausted its deadline before CAD proximity")
            proximity = query_facial_proximity(self, found, cad, face_tags=cad_face_tags,
                max_work=max(1, min(cad_work, remaining_work)), max_seconds=min(number(max_cad_seconds, "max_cad_seconds", minimum=0), remaining_seconds),
                include_faces=cad_include_faces, _reserved_bytes=query_reservation)
            proximity_rows = {row["cell_id"]: row for row in proximity["cells"]}
            evaluations = proximity["evaluations"]
            self._deadline(started, base_work + evaluations)
            if cad_threshold is not None:
                found = [index for index in found if proximity_rows[index]["minimum_distance"] <= cad_threshold]
        end = min(offset + limit, len(found))
        records = [self._record(index, associations=associations) for index in found[offset:end]]
        if proximity is not None:
            for record in records:
                record["nearby_cad_surfaces"] = {"available": True, "method": proximity["method"],
                    "sample": proximity["sample"], "sample_scope": proximity["sample_scope"],
                    "does_not_establish": proximity["does_not_establish"],
                    "limit_scope": proximity["limit_scope"], **proximity_rows[record["cell_id"]]}
        self._deadline(started, base_work + evaluations)
        result = {"schema": "paacs-mesh-cell-page-v1", **self.identity(), "total": len(found),
                  "offset": offset, "limit": limit, "next_offset": end if end < len(found) else None,
                  "metric_conventions": METRIC_CONVENTIONS,
                  "filter": {"metric": metric, "quality_below": threshold, "region": box,
                             "spatial_mode": spatial_mode, "body_tags": sorted(bodies) if bodies is not None else None,
                             "cell_type": cell_type, "orientation": orientation,
                             "cad_distance_inclusive": cad_threshold,
                             "cad_face_tags": proximity["selected_face_tags"] if proximity is not None else None},
                  "cad_proximity": {key: value for key, value in proximity.items() if key != "cells"} if proximity is not None else None,
                  "sort": {"by": sort_by, "descending": descending}, "cells": records}
        copied = json_copy(result)
        self._deadline(started, base_work + evaluations)
        return copied

    def summary(self, *, spatial_bins=4):
        bins = integer(spatial_bins, "spatial_bins")
        if bins > 16:
            raise ValueError("spatial_bins must not exceed 16 per axis")
        report = self._archive.report or {}
        sizing = report.get("sizing")
        sizing = {} if sizing is None else sizing
        if not isinstance(sizing, dict):
            raise InspectionError("inspection_report", "Native sizing metadata must be an object")
        for name in ("face_sizes", "volume_controls"):
            controls = sizing.get(name, [])
            if not isinstance(controls, list) or any(not isinstance(item, dict) for item in controls):
                raise InspectionError("inspection_report", f"Native sizing.{name} must be an object list")
        # A boundary triangle can introduce a distinct face group. Admit each
        # group's distribution, dictionary keys and JSON-copy overhead before
        # building the map, rather than budgeting only its three edge values.
        self._reserve(384 * self.cell_count + 4096 * self._archive.type_counts[2] + 8192 * bins**3,
                      "Quality distributions")
        started = time.perf_counter()
        self._deadline(started, 12 * self.cell_count + self.node_count)
        measured = [index for index in range(self.cell_count) if self._orientation[index] != 2]
        edges = {"mean_ratio": (.05, .1, .2, .25, .5, .8),
                 "minimum_dihedral_deg": (1, 5, 10, 20, 40, 60),
                 "scaled_jacobian": (0, .05, .1, .2, .5, .8)}
        quality = {name: distribution((self._metric(index, name) for index in measured), edges=edges[name]) for name in METRICS}
        size = {name: distribution(self._metric(index, name) for index in measured)
                for name in ("edge_min", "edge_median", "edge_max", "edge_rms", "volume_equivalent_edge")}
        extent = None
        for node in range(self.node_count):
            if node % 128 == 0:
                self._deadline(started, 12 * self.cell_count + node + 1)
            point = self._archive.node(node)[0]
            if extent is None:
                extent = [point[:], point[:]]
            else:
                for axis in range(3):
                    extent[0][axis] = min(extent[0][axis], point[axis])
                    extent[1][axis] = max(extent[1][axis], point[axis])
        spatial, body_counts, face_lengths = {}, Counter(), {}
        unclassified_tetrahedra = 0
        for index in range(self.cell_count):
            if index % 128 == 0:
                self._deadline(started, index + 1)
            kind, nodes, classification = self._archive.cell(index)
            if kind == 2 and classification["dimension"] == 2:
                points = [self._archive.node(n)[0] for n in nodes]
                target = face_lengths.setdefault(str(classification["tag"]), [])
                target.extend(value if math.isfinite(value) else None
                              for value in (math.dist(points[a], points[b]) for a, b in ((0, 1), (1, 2), (2, 0))))
            if kind != 5:
                continue
            center = _centroid([self._archive.node(n)[0] for n in nodes])
            key = []
            for axis in range(3):
                low, high = extent[0][axis], extent[1][axis]
                scale = max(abs(low), abs(high)) or 1.
                fraction = ((center[axis] / scale - low / scale) / (high / scale - low / scale)) if high > low else 0.
                key.append(max(0, min(bins - 1, int(fraction * bins))))
            key = tuple(key)
            record = spatial.setdefault(key, {"count": 0, "minimum_mean_ratio": 1., "below_0_2": 0})
            value = self._metric(index, "mean_ratio")
            record["count"] += 1
            record["minimum_mean_ratio"] = min(record["minimum_mean_ratio"], value)
            record["below_0_2"] += value < .2
            local = record.setdefault("quality", {})
            for name in METRICS:
                measured_value = self._metric(index, name)
                population = local.setdefault(name, {"count": 0, "unavailable": 0, "minimum": None,
                    "maximum": None, "_sum": 0., "_compensation": 0., "histogram": {"edges": list(edges[name]),
                    "counts": [0] * (len(edges[name]) + 1)}})
                if measured_value is None:
                    population["unavailable"] += 1
                    continue
                population["count"] += 1
                population["minimum"] = measured_value if population["minimum"] is None else min(population["minimum"], measured_value)
                population["maximum"] = measured_value if population["maximum"] is None else max(population["maximum"], measured_value)
                adjusted = measured_value - population["_compensation"]
                total = population["_sum"] + adjusted
                population["_compensation"] = (total - population["_sum"]) - adjusted
                population["_sum"] = total
                population["histogram"]["counts"][sum(measured_value >= edge for edge in edges[name])] += 1
            if classification["dimension"] == 3:
                body_counts[str(classification["tag"])] += 1
            else:
                unclassified_tetrahedra += 1
        for record in spatial.values():
            for population in record["quality"].values():
                total = population.pop("_sum")
                population["mean"] = total / population["count"] if population["count"] else None
                population.pop("_compensation")
        requested = {"global_size": sizing.get("global_size"), "min_size": sizing.get("requested_min_size"),
                     "growth": sizing.get("growth"), "native_sizing": sizing.get("requested_native_sizing"),
                     "face_sizes": [{"tag": item.get("tag"), "size": item.get("size"),
                                     "growth": item.get("requested_growth")} for item in sizing.get("face_sizes", [])],
                     "volume_controls": [{key: item.get(key) for key in ("kind", "lower", "upper", "size", "requested_growth")}
                                         for item in sizing.get("volume_controls", [])]}
        effective = {"minimum_target": sizing.get("min_size"), "global_target": sizing.get("global_size"),
                     "face_targets": [{"tag": item.get("tag"), "size": item.get("effective_size"),
                                       "growth": item.get("growth")} for item in sizing.get("face_sizes", [])],
                     "volume_targets": [{key: item.get(key) for key in ("kind", "lower", "upper", "effective_size", "growth")}
                                        for item in sizing.get("volume_controls", [])],
                     "resolved_routing": report.get("effective_policy"),
                     "sizing_completion_certificate": "unavailable",
                     "note": "Effective field targets and observed routes do not certify every achieved edge length"}
        result = {"schema": "paacs-mesh-quality-summary-v1", "scope": "entire_snapshot", **self.identity(), "validation": self.validation,
                  "metric_conventions": METRIC_CONVENTIONS, "node_count": self.node_count,
                  "cell_count": self.cell_count, "cell_counts": dict(zip(CELL_NAMES, self._archive.type_counts)),
                  "quality": quality, "orientation_counts": {name: self._orientation.count(value) for name, value in ORIENTATIONS.items()},
                  "range_limited_cells": sum(self._range_limited),
                  "sizing": {"requested": requested, "effective": effective,
                             "achieved": {"tet4": size, "boundary_faces": {tag: distribution(values) for tag, values in face_lengths.items()}},
                             "boundary_population": "triangle-edge incidences; shared edges appear in each incident triangle",
                             "tet4_population": "per-cell edge and regular-volume-equivalent-edge distributions",
                             "volume_equivalent_edge_caveat": "Volume proxy is affected by element shape; it is not a maximum-edge or accuracy certificate",
                             "target_is_not_a_guarantee_of_every_edge_length": True},
                  "body_tet4_counts": dict(body_counts), "unclassified_tet4_count": unclassified_tetrahedra, "bounds": extent,
                  "spatial_distribution": {"kind": "uniform_axis_aligned_centroid_bins", "bins_per_axis": bins,
                                           "bounds": extent, "bins": [{"index": list(key), **value} for key, value in sorted(spatial.items())]},
                  "native_timings": report.get("seconds"),
                  "inspection": {"analysis_seconds": self.analysis_seconds, "owned_storage_bytes": self.storage_bytes,
                                 "max_bytes": self.max_bytes, "work_budget_scope": "per operation",
                                 "not_a_solver_error_or_convergence_estimate": True}}
        copied = json_copy(result)
        self._deadline(started, 12 * self.cell_count + self.node_count)
        return copied

    def save(self, path):
        """Save a new hash-paired inspection bundle; never replace an existing path."""
        destination = Path(path).absolute()
        if os.path.lexists(destination):
            raise FileExistsError("Inspection bundles never replace existing paths")
        if not destination.parent.is_dir():
            raise ValueError("Bundle parent must already exist")
        report = self._archive.report
        if report is None:
            raise InspectionError("report_unavailable", "A bare archive has no native report to bind")
        with tempfile.TemporaryDirectory(dir=destination.parent, prefix=".paacs-inspection-") as temporary:
            staging = Path(temporary) / "bundle"
            staging.mkdir()
            shutil.copyfile(self._archive.path, staging / "mesh.paacsmesh")
            (staging / "report.json").write_text(json_text(report), encoding="utf-8")
            manifest = {"schema": "paacs-mesh-inspection-bundle-v1", **self.identity(),
                        "report_sha256": file_hash(staging / "report.json"), "validation_basis": self._archive.report_basis}
            (staging / "inspection.json").write_text(json_text(manifest), encoding="utf-8")
            os.rename(staging, destination)
        return destination

    def write_vtk(self, path, *, max_output_bytes=128 * 1024 * 1024):
        """Export original IDs/classifications and measured Tet4 fields as legacy VTK."""
        maximum = integer(max_output_bytes, "max_output_bytes")
        self._reserve(8192, "VTK streaming buffer")
        protected = self._archive.protected_paths
        started, output_bytes, work = time.perf_counter(), 0, 0
        self._deadline(started, 16 * self.cell_count + self.node_count)
        with atomic_stream(path, protected=protected) as stream:
            for line in self._vtk_lines():
                work += 1
                if work % 128 == 0:
                    self._deadline(started, work)
                data = line.encode("ascii")
                output_bytes += len(data)
                if output_bytes > maximum:
                    raise InspectionError("output_budget", "VTK output exceeds max_output_bytes")
                stream.write(data)
            self._deadline(started, work)
        return Path(path).absolute()

    def _vtk_lines(self):
        yield f"# vtk DataFile Version 3.0\nPAACS mesh_sha256={self.mesh_sha256}; source_sha256={self.source_sha256 or 'unknown'}\nASCII\nDATASET UNSTRUCTURED_GRID\n"
        yield f"POINTS {self.node_count} double\n"
        for node in range(self.node_count):
            yield " ".join(format(v, ".17g") for v in self._archive.node(node)[0]) + "\n"
        total = sum((1 + width) * count for width, count in zip((2, 3, 3, 6, 4, 4, 10, 6, 8, 5), self._archive.type_counts))
        yield f"CELLS {self.cell_count} {total}\n"
        for index in range(self.cell_count):
            _, nodes, _ = self._archive.cell(index)
            yield str(len(nodes)) + " " + " ".join(map(str, nodes)) + "\n"
        yield f"CELL_TYPES {self.cell_count}\n"
        types = (3, 21, 5, 22, 9, 10, 24, 13, 12, 14)
        for index in range(self.cell_count):
            yield str(types[self._archive.cell(index)[0]]) + "\n"
        yield f"CELL_DATA {self.cell_count}\n"
        for field, dtype, getter in (("native_cell_id", "long", lambda index: index),
                ("entity_dimension", "int", lambda index: self._archive.cell(index)[2]["dimension"]),
                ("entity_tag", "int", lambda index: self._archive.cell(index)[2]["tag"]),
                ("quality_supported", "int", lambda index: int(self._orientation[index] != 2)),
                ("orientation_code", "int", lambda index: self._orientation[index]),
                ("range_limited", "int", lambda index: self._range_limited[index])):
            yield f"SCALARS {field} {dtype} 1\nLOOKUP_TABLE default\n"
            for index in range(self.cell_count):
                yield str(getter(index)) + "\n"
        for field in MEASUREMENTS:
            yield f"SCALARS {field} double 1\nLOOKUP_TABLE default\n"
            for value in self._columns[field]:
                yield format(value, ".17g") + "\n"
        yield f"POINT_DATA {self.node_count}\n"
        for name, getter in (("native_node_id", lambda index: index),
                             ("entity_dimension", lambda index: self._archive.node(index)[1]["dimension"]),
                             ("entity_tag", lambda index: self._archive.node(index)[1]["tag"])):
            yield f"SCALARS {name} long 1\nLOOKUP_TABLE default\n"
            for index in range(self.node_count):
                yield str(getter(index)) + "\n"

    def visualize(self, path, **options):
        from .mesh_visualization import visualize_mesh
        return visualize_mesh(self, path, **options)

    def close(self):
        self._archive.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def inspect_mesh(source, **limits):
    """Measure a native Mesh, archive or saved bundle; NumPy is not required."""
    return MeshInspection(source, **limits)


def open_mesh(path, **limits):
    """Open an existing archive/bundle for bounded structural/shape inspection."""
    return MeshInspection(path, **limits)
