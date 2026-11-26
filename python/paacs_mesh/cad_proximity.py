"""Explicit cell-sample proximity to selected, truly trimmed CAD faces.

OCCT computes point-to-shape distance including face trims, edges and vertices:
https://occt3d.com/dev/doc/refman/html/class_b_rep_extrema___dist_shape_shape.html
These are sampled distances, not full-cell distances, boundary ownership,
quality causation, CAD error bounds, or a new mesh validity certificate.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from bisect import bisect_left
import time

from ._cad_common import integer, number, point, sequence
from .mesh_archive import InspectionError


MAX_CAD_EVALUATIONS = 100_000
_NATIVE_RESULT_BYTES = 8 * 1024 * 1024
_MAX_FACES = 1000


def proximity_scratch_bytes(candidate_count, face_count, *, include_faces=False):
    """Conservative owned result plus transient native chunk reservation."""
    candidates = integer(candidate_count, "candidate_count", minimum=0)
    faces = integer(face_count, "face_count", minimum=0)
    if not isinstance(include_faces, bool):
        raise ValueError("include_faces must be bool")
    if candidates == 0:
        return 65536
    return (4096 * candidates + (2048 if include_faces else 256) * candidates * faces
            + _NATIVE_RESULT_BYTES + 65536)


def _count(value, name, maximum):
    if isinstance(value, (str, bytes, bytearray, dict)):
        raise ValueError(f"{name} must be a bounded sequence")
    try:
        count = len(value)
    except TypeError as error:
        raise ValueError(f"{name} must be a bounded sequence") from error
    if count > maximum:
        raise InspectionError("cad_work_budget", f"{name} exceeds its diagnostic cap; restrict cells or face tags")
    return count


def _check_row(row, sample, tags, include_faces):
    """Check completeness before a row may participate in exact-total filtering."""
    try:
        nearest = sequence(row["nearest_face_tags"], "nearest_face_tags", maximum=len(tags))
        if (row["point"] != sample or row["measured_face_count"] != len(tags)
                or row["closest_face_tag"] not in tags or row["closest_face_tag"] not in nearest
                or len(nearest) != len(set(nearest)) or any(tag not in tags for tag in nearest)):
            raise ValueError("Incomplete selected-face receipt")
        number(row["minimum_distance"], "minimum_distance", minimum=0)
        number(row["face_tolerance"], "face_tolerance", minimum=0)
        point(row["closest_point"], "closest_point")
        if include_faces:
            faces = sequence(row["faces"], "faces", length=len(tags), maximum=len(tags))
            if sorted(face["face_tag"] for face in faces) != tags:
                raise ValueError("Selected-face measurements differ")
            for face in faces:
                number(face["distance"], "distance", minimum=0)
                number(face["tolerance"], "tolerance", minimum=0)
                point(face["closest_point"], "closest_point")
            if min(face["distance"] for face in faces) != row["minimum_distance"]:
                raise ValueError("Selected-face minimum differs")
    except (ValueError, TypeError, KeyError) as error:
        raise InspectionError("cad_query_receipt", "Native CAD proximity sample receipt is incomplete") from error


def query_facial_proximity(inspection, candidate_cell_ids, model, *, face_tags=None,
                          max_work=MAX_CAD_EVALUATIONS, max_seconds=10,
                          include_faces=False, _reserved_bytes=0):
    """Measure every selected face at each candidate's arithmetic node centroid.

    For Tet4 this sample is its geometric centroid. For nonlinear or other cell
    families it is explicitly a node-mean sample, not an integrated geometric
    centroid. IDs remain native archive indices. Source SHA-256, kernel version
    and millimetre units must match; bare archives have no such provenance.

    All selected-face evaluations must finish or the operation raises. An OCCT
    solve and its internal allocations are not preemptible or RSS bounded.
    PAACS limits inputs/output and checks its cooperative deadline between
    solves, including native operation-gate waiting. Restrict cells/regions or
    face tags when candidate_count * face_count exceeds the admitted budget.
    ``_reserved_bytes`` accounts for an enclosing query's simultaneous scratch.
    """
    from .cad import CadModel
    from .mesh_inspection import MeshInspection, _centroid

    if not isinstance(inspection, MeshInspection) or not isinstance(model, CadModel):
        raise TypeError("CAD proximity requires a MeshInspection and CadModel")
    maximum = integer(max_work, "max_work")
    if maximum > MAX_CAD_EVALUATIONS:
        raise ValueError("max_work must not exceed 100000 CAD distance evaluations")
    seconds = number(max_seconds, "max_seconds", minimum=0)
    if seconds <= 0:
        raise ValueError("max_seconds must be positive")
    if not isinstance(include_faces, bool):
        raise ValueError("include_faces must be bool")
    reserved = integer(_reserved_bytes, "_reserved_bytes", minimum=0)
    started = time.perf_counter()

    def remaining():
        value = seconds - (time.perf_counter() - started)
        if value <= 0:
            raise InspectionError("cad_time_budget", "CAD proximity exceeded cooperative max_seconds")
        return value

    model._assert_identity()
    report = inspection._archive.report or {}
    if (inspection.source_sha256 is None or report.get("kernel_version") is None
            or inspection.units != "mm"):
        raise InspectionError("cad_identity_unavailable", "CAD proximity requires guarded source, kernel and millimetre units")
    if inspection.source_sha256 != model.source_sha256 or report["kernel_version"] != model.kernel_version:
        raise InspectionError("cad_identity_mismatch", "Inspection and CAD model must have the same source and kernel")
    if inspection._archive._mapped is None:
        raise InspectionError("closed", "Inspection snapshot is closed")
    known_faces = model._tag_sets["faces"]
    selected = known_faces if face_tags is None else face_tags
    candidate_count = _count(candidate_cell_ids, "candidate_cell_ids", maximum)
    face_count = _count(selected, "face_tags", _MAX_FACES)
    if face_count == 0:
        raise ValueError("face_tags must not be empty")
    evaluations = candidate_count * face_count
    effective_maximum = min(maximum, inspection.max_work)
    if evaluations > effective_maximum:
        raise InspectionError("cad_work_budget", "CAD proximity exceeds max_work; restrict cells, regions or face tags",
                              details={"candidate_cells": candidate_count, "selected_faces": face_count,
                                       "evaluations": evaluations, "maximum": effective_maximum})
    inspection._reserve(reserved + proximity_scratch_bytes(candidate_count, face_count,
                                                          include_faces=include_faces), "CAD proximity")
    ids = [integer(value, "cell_id", minimum=0)
           for value in sequence(candidate_cell_ids, "candidate_cell_ids", maximum=maximum)]
    tags = [integer(value, "face_tag") for value in sequence(selected, "face_tags", maximum=_MAX_FACES)]
    if len(set(ids)) != len(ids) or len(set(tags)) != len(tags):
        raise ValueError("Cell IDs and face tags must be unique")
    if any(value >= inspection.cell_count for value in ids):
        raise ValueError("cell_id is outside this mesh")
    for tag in tags:
        position = bisect_left(known_faces, tag)
        if position == len(known_faces) or known_faces[position] != tag:
            raise ValueError(f"Unknown CAD face tag {tag}")
    tags.sort()
    # Native admission includes both DTOs and Python rows; shrink chunks before
    # entering native code instead of relying on a failed allocation to retry.
    bytes_per_point = 4096 + face_count * (1024 if include_faces else 160)
    chunk_size = min(64, (_NATIVE_RESULT_BYTES - 65536) // bytes_per_point)
    if chunk_size < 1:
        raise InspectionError("cad_memory_budget", "One selected-face sample exceeds the native output cap")
    records = []
    for offset in range(0, len(ids), chunk_size):
        remaining()
        chunk_ids = ids[offset:offset + chunk_size]
        samples = []
        for cell_id in chunk_ids:
            _, nodes, _ = inspection._archive.cell(cell_id)
            samples.append(_centroid([inspection._archive.node(node)[0] for node in nodes]))
        try:
            native = model._native._face_distances(samples, tags, max_work=len(samples) * face_count,
                                                    max_seconds=remaining(), include_faces=include_faces,
                                                    max_result_bytes=_NATIVE_RESULT_BYTES)
        except MemoryError as error:
            raise InspectionError("cad_kernel_allocation", "OCCT could not allocate its CAD proximity workspace") from error
        except (ValueError, RuntimeError, OverflowError) as error:
            if time.perf_counter() - started >= seconds:
                raise InspectionError("cad_time_budget", "CAD proximity exceeded cooperative max_seconds; no partial result returned",
                                      details={"native_message": str(error)[:1024]}) from error
            raise InspectionError("cad_query_failed", "Trimmed CAD proximity did not complete; no partial result returned",
                                  details={"native_message": str(error)[:1024]}) from error
        model._assert_identity()
        if (native.get("schema") != "paacs-cad-face-distances-v1"
                or native.get("source_sha256") != inspection.source_sha256
                or native.get("kernel_version") != model.kernel_version
                or native.get("evaluations") != len(samples) * face_count
                or len(native.get("rows", ())) != len(samples)):
            raise InspectionError("cad_query_identity", "Native CAD proximity receipt differs from its request")
        for cell_id, sample, row in zip(chunk_ids, samples, native["rows"]):
            _check_row(row, sample, tags, include_faces)
            result = dict(row)
            result.pop("point")
            result.update(cell_id=cell_id, sample_point=sample)
            records.append(result)
        remaining()
    remaining()
    return {"schema": "paacs-mesh-cad-proximity-v1", **inspection.identity(),
            "kernel_version": model.kernel_version, "cad_model_identity": model.model_identity,
            "method": "occt_trimmed_face_distance", "sample": "arithmetic_cell_node_centroid",
            "sample_scope": "Tet4 geometric centroid; other families use a node mean, not integrated curved geometry",
            "selected_face_tags": tags, "evaluations": evaluations, "complete": True,
            "nearest_tag_rule": "distance difference <= max(native face tolerances,32epsilon*max(1,minimum_distance))",
            "does_not_establish": ["boundary_ownership", "quality_causation", "full_cell_distance", "CAD_error_bound", "mesh_validity"],
            "limit_scope": "bounded PAACS inputs/output; one OCCT solve is not preemptible or RSS bounded",
            "cells": records}
