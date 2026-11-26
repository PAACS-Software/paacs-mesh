"""Opt-in frozen-boundary Tet4 research; saved claims never replace CAD checks.

Bundles contain exact source STEP geometry and are private CAD artifacts. Loading
reimports that STEP and revalidates the untrusted tagged boundary natively.
``points=None`` selects classical placement; ``points=[]`` selects supplied mode
with no interior seeds. Other supplied clouds are all-or-nothing, in mm.
"""
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
from __future__ import annotations

import hashlib
from functools import lru_cache
import json
import math
import os
from pathlib import Path
import stat
import struct
import tempfile
import time

from . import _native
from .cad import CadModel, read_step
from ._cad_common import integer, number

_SCHEMA = "paacs-experimental-boundary-bundle-v1"
_SOURCE_LIMIT = 64 * 1024**2
_METADATA_LIMIT = 1024**2
_ARCHIVE_LIMIT = 512 * 1024**2
_POINT_LIMIT = 1_000_000
_FILES = ("source.step", "boundary.paacsmesh", "boundary-report.json")
_DEFAULTS = {"max_nodes": 200000, "max_cells": 1000000,
             "max_work": (2**64-1)//4, "max_working_bytes": 512*1024**2,
             "max_output_bytes": 128*1024**2, "max_audit_bytes": 128*1024**2,
             "max_seconds": 300.}


def _hash(data):
    return hashlib.sha256(data).hexdigest()


def _json_bytes(value):
    result = (json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")
    if len(result) > _METADATA_LIMIT:
        raise ValueError("Boundary metadata exceeds 1 MiB")
    return result


def _json(data):
    if len(data) > _METADATA_LIMIT:
        raise ValueError("Boundary metadata exceeds 1 MiB")
    def pairs(items):
        value = {}
        for key, item in items:
            if key in value:
                raise ValueError("Duplicate boundary metadata key")
            value[key] = item
        return value
    def constant(value):
        raise ValueError("Nonfinite boundary metadata")
    try:
        result = json.loads(data, object_pairs_hook=pairs, parse_constant=constant)
        _json_bytes(result)
    except (RecursionError, UnicodeError) as error:
        raise ValueError("Malformed boundary metadata") from error
    if not isinstance(result, dict):
        raise ValueError("Boundary metadata must be an object")
    return result


def _snapshot(path, maximum):
    """Snapshot bounded regular files, including descriptor checks after open."""
    path = Path(path)
    if path.is_symlink():
        raise ValueError("Boundary inputs must not be symbolic links")
    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0) | getattr(os, "O_NONBLOCK", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    with os.fdopen(descriptor, "rb") as stream:
        info = os.fstat(stream.fileno())
        if not stat.S_ISREG(info.st_mode) or not 0 < info.st_size <= maximum:
            raise ValueError("Boundary input is not a bounded nonempty regular file")
        payload = stream.read(info.st_size + 1)
        if len(payload) != info.st_size:
            raise ValueError("Boundary input changed during snapshot")
    return payload


def _options(size, max_chord_error=None, max_relative_volume_error=.01, **limits):
    if set(limits) - set(_DEFAULTS):
        raise ValueError("Unsupported experimental control; only uniform single-solid Tet4 is admitted")
    result = {**_DEFAULTS, **limits}
    for name in _DEFAULTS:
        result[name] = number(result[name], name) if name == "max_seconds" else integer(result[name], name)
        if result[name] <= 0:
            raise ValueError(name + " must be positive")
    if result["max_seconds"] > 86400:
        raise ValueError("Cooperative deadline exceeds 86400 seconds")
    size = number(size, "size")
    chord = None if max_chord_error is None else number(max_chord_error, "max_chord_error")
    relative = number(max_relative_volume_error, "max_relative_volume_error")
    if size <= 0 or (chord is not None and chord <= 0) or not 0 < relative < 1:
        raise ValueError("Positive uniform size/chord and relative volume error below one are required")
    return {"size": size, "max_chord_error": chord, "max_relative_volume_error": relative, **result}


@lru_cache(maxsize=4)
def _extension_hash(path, revision, dirty):
    """Capture loaded-build provenance once, without reading an unbounded file."""
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_BINARY", 0) | getattr(os, "O_NONBLOCK", 0))
    with os.fdopen(descriptor, "rb") as stream:
        info = os.fstat(stream.fileno())
        if not stat.S_ISREG(info.st_mode) or not 0 < info.st_size <= 64*1024**2:
            raise ValueError("Native extension exceeds the 64 MiB provenance hash allowance")
        digest, remaining = hashlib.sha256(), info.st_size
        while remaining:
            block = stream.read(min(remaining, 65536))
            if not block:
                raise ValueError("Native extension changed during provenance hashing")
            digest.update(block)
            remaining -= len(block)
        after = os.fstat(stream.fileno())
        if stream.read(1) or (info.st_size, info.st_mtime_ns, info.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
            raise ValueError("Native extension changed during provenance hashing")
        return digest.hexdigest()


def _build():
    return {"engine_version": _native.engine_version,
            "source_revision": _native.source_revision, "source_dirty": _native.source_dirty,
            "native_extension_sha256": _extension_hash(str(Path(_native.__file__).resolve(strict=True)),
                                                       _native.source_revision, _native.source_dirty)}


def _topology(model):
    """Integer incidence/orientation fingerprint, without process-local tokens."""
    result, used = {}, 0
    for kind in ("bodies", "faces", "edges", "vertices"):
        records, offset = [], 0
        while True:
            page = model.inspect(kind, offset=offset, limit=100, topology=kind == "faces")
            total, entities = page["total"], page["entities"]
            if (type(total) is not int or total > 10000 or page["offset"] != offset
                    or not isinstance(entities, list) or len(entities) > 100
                    or (offset < total and not entities)):
                raise ValueError("CAD topology exceeds fingerprint allowance or is incomplete")
            for entity in entities:
                record = {"tag": entity["tag"]}
                if kind == "bodies":
                    record["face_uses"] = entity["face_uses"]
                elif kind == "faces":
                    record.update(surface=entity["surface"], owners=entity["owners"],
                        wires=[{"outer": wire["outer"], "edges": [{key: edge[key] for key in
                            ("edge_tag", "start_vertex_tag", "end_vertex_tag", "orientation", "sense", "seam", "degenerated")}
                            for edge in wire["edges"]]} for wire in entity["wires"]])
                    # Normal vectors are floating measurements rather than incidence.
                    record["owners"] = [{key: owner[key] for key in ("body_tag", "sense")} for owner in record["owners"]]
                elif kind == "edges":
                    record.update({key: entity[key] for key in ("curve_type", "start_vertex_tag", "end_vertex_tag", "face_tags", "periodic", "degenerated")})
                used += len(_json_bytes(record))
                if used > _METADATA_LIMIT - 65536:
                    raise ValueError("CAD topology fingerprint exceeds its metadata allowance")
                records.append(record)
            offset += len(entities)
            if page["next_offset"] is None:
                if offset != total:
                    raise ValueError("CAD topology pagination is incomplete")
                break
            if page["next_offset"] != offset:
                raise ValueError("CAD topology pagination is inconsistent")
        result[kind] = records
        _json_bytes(result)
    return _hash(_json_bytes(result))


def _point_bytes(points, options):
    if points is None:
        return None
    if isinstance(points, (str, bytes, bytearray)) or not hasattr(points, "__len__") or not hasattr(points, "__getitem__"):
        raise ValueError("points must be a bounded sequence of finite triples, or None")
    count = len(points)
    if count > min(_POINT_LIMIT, options["max_nodes"]) or count*24 > options["max_working_bytes"]//4:
        raise ValueError("Supplied point count/storage exceeds allowance")
    payload = bytearray(count*24)
    for i in range(count):
        point = points[i]
        if isinstance(point, (str, bytes)) or not hasattr(point, "__len__") or not hasattr(point, "__getitem__") or len(point) != 3:
            raise ValueError("Each supplied point must contain three finite numbers")
        struct.pack_into("<3d", payload, i*24, *(number(point[j], "point coordinate") for j in range(3)))
    return bytes(payload)


class PreparedBoundary:
    """Owned exact-certified Tri3 boundary and frozen uniform sizing contract.

    The report certifies boundary preparation only, never volume filling or
    solver adequacy. Saving includes source CAD. Resource limits are modeled
    Native/encoded-buffer allowances, not a process RSS/hard timeout guarantee.
    """
    __slots__ = ("_model", "_boundary", "_source", "_archive", "_options", "_metadata", "_report")

    def __init__(self, model, boundary, source, archive, options, metadata, report):
        self._model, self._boundary = model, boundary
        self._source, self._archive = bytes(source), bytes(archive)
        self._options, self._metadata, self._report = _json_bytes(options), _json_bytes(metadata), _json_bytes(report)

    @property
    def report(self):
        return _json(self._report)

    @property
    def source_sha256(self):
        return _hash(self._source)

    @property
    def boundary_sha256(self):
        return _hash(self._archive)

    @property
    def source_bytes(self):
        """Exact original STEP bytes; this is private CAD data."""
        return self._source

    def to_numpy(self, *, max_bytes=128*1024**2):
        """Copy owned read-only boundary points/Tri3/CAD-class arrays; NumPy optional."""
        return self._boundary.to_numpy(max_bytes=integer(max_bytes, "max_bytes"))

    def mesh(self, points=None):
        """Classical placement or strict supplied seeds through the same kernel.

        Native refusal retains its report. Invalid or unretained supplied points
        cannot fall back to ordinary placement. The frozen boundary is immutable.
        ``Mesh.initial_points(max_bytes=...)`` copies the actual seed cloud.
        """
        began = time.perf_counter()
        packing, called, supplied_hash, extension_hash = None, False, None, None
        try:
            self._model._assert_identity()
            extension_hash = _build()["native_extension_sha256"]
            options = _json(self._options)
            payload = _point_bytes(points, options)
            supplied_hash = None if payload is None else _hash(payload)
            packing = time.perf_counter() - began
            called = True
            mesh = _native._experimental_mesh_boundary(self._model._native, self._boundary, options, payload)
            mesh._experimental_set_python_timing(time.perf_counter() - began, packing, extension_hash)
            return mesh
        except BaseException as error:
            report = getattr(error, "report", None)
            report = dict(report) if isinstance(report, dict) else {}
            prior = report.get("experimental_boundary_volume")
            details = dict(prior) if isinstance(prior, dict) else {}
            details.update(python_total_seconds=time.perf_counter()-began, point_packing_seconds=packing,
                           supplied_input_sha256=supplied_hash,
                           native_extension_sha256=extension_hash,
                           native_binding_called=called,
                           python_timing_scope="Failed Python replay snapshot includes request admission/packing and any invoked native binding work, including post-Native hashes; precedes exception timing attachment/propagation. Native phases are incomplete separate spans.")
            report["experimental_boundary_volume"] = details
            error.report = report
            raise

    def save(self, directory):
        """Exclusively create a bundle; publish its manifest after complete files."""
        directory = Path(directory).absolute()
        files = {"source.step": self._source, "boundary.paacsmesh": self._archive,
                 "boundary-report.json": self._report}
        manifest = {"schema": _SCHEMA, "units": "mm", **_json(self._metadata),
                    "options": _json(self._options),
                    "files": {name: {"sha256": _hash(data), "bytes": len(data)} for name, data in files.items()}}
        encoded = _json_bytes(manifest)
        directory.mkdir(mode=0o700, exist_ok=False)
        owned = []
        try:
            for name, data in [*files.items(), ("manifest.json", encoded)]:
                path = directory / name
                with path.open("xb") as stream:
                    owned.append(path)
                    stream.write(data)
                    stream.flush()
                    os.fsync(stream.fileno())
        except BaseException:
            for path in reversed(owned):
                path.unlink(missing_ok=True)
            directory.rmdir()
            raise
        return directory


def initial_points(value, *, max_bytes=16*1024**2):
    """Copy a replay result/refusal's actual pre-constraint cloud, never a mesh.

    Refusals expose a separate at-most-64-MiB immutable IEEE64 payload. An
    unavailable capture refuses explicitly; a captured empty cloud returns ().
    ``max_bytes`` conservatively bounds decoded tuple storage as well as bytes.
    """
    maximum = integer(max_bytes, "max_bytes")
    if not isinstance(value, BaseException):
        return value.initial_points(max_bytes=maximum)
    payload = getattr(value, "initial_points_bytes", None)
    if payload is None:
        raise ValueError("Initial point snapshot unavailable: " + str(getattr(value, "initial_points_unavailable_reason", "not_captured")))
    if not isinstance(payload, bytes) or len(payload) > 64*1024**2 or len(payload) % 24:
        raise ValueError("Malformed initial point snapshot")
    report = getattr(value, "report", {})
    details = report.get("experimental_boundary_volume", {}) if isinstance(report, dict) else {}
    if (not isinstance(details, dict) or details.get("initial_points_available") is not True
            or details.get("initial_points") != len(payload)//24
            or details.get("initial_points_sha256") != _hash(payload)):
        raise ValueError("Initial point snapshot differs from its scalar receipt")
    if maximum < 65536 or len(payload)//24 > (maximum-65536)//256:
        raise ValueError("Initial point snapshot exceeds decoded max_bytes")
    if any(not math.isfinite(v) for point in struct.iter_unpack("<3d", payload) for v in point):
        raise ValueError("Initial point snapshot is nonfinite")
    return tuple(struct.iter_unpack("<3d", payload))


def prepare_boundary(model, *, size, max_chord_error=None, max_relative_volume_error=.01, **limits):
    """Prepare and retain a qualified boundary BEFORE any volume/seed attempt.

    First scope: one original solid, uniform sizing, conforming Tet4, fixed
    Standard recovery and the Normal classical front. Unsupported controls
    refuse. Prefer ``read_step(..., retain_source=True)``; otherwise the original
    regular file must still have exactly its imported source hash.
    """
    began = time.perf_counter()
    if not isinstance(model, CadModel):
        raise TypeError("model must be an imported CadModel")
    model._assert_identity()
    build = _build()
    options = _options(size, max_chord_error, max_relative_volume_error, **limits)
    if len(model.body_tags) != 1:
        raise ValueError("Experimental boundary replay requires exactly one original solid")
    source_size = model.info.get("file_bytes")
    if type(source_size) is not int or not 0 < source_size <= _SOURCE_LIMIT:
        raise ValueError("Imported source exceeds the experimental source byte allowance")
    source = model._native._source_bytes
    if source is None:
        source = _snapshot(model.info["source_path"], _SOURCE_LIMIT)
    if not isinstance(source, bytes) or len(source) != source_size or _hash(source) != model.source_sha256:
        raise ValueError("Exact imported source is unavailable/changed; reimport with retain_source=True")
    topology = _topology(model)
    boundary = _native._experimental_prepare_boundary(model._native, options)
    maximum = min(_ARCHIVE_LIMIT, options["max_output_bytes"], options["max_working_bytes"]//8)
    archive = boundary.archive(maximum)
    report = boundary.report
    metadata = {"source_sha256": _hash(source), "boundary_sha256": _hash(archive),
                "topology_sha256": topology, "kernel_version": model.kernel_version, "build": build}
    report["provenance"] = metadata
    report["python_prepare_seconds"] = time.perf_counter() - began
    report["timing_scope"] = "Preparation snapshot includes Python source/fingerprint/archive/hash work and Native preparation; precedes final metadata encoding/return; no volume call"
    return PreparedBoundary(model, boundary, source, archive, options, metadata, report)


def load_boundary(directory, *, max_source_bytes=_SOURCE_LIMIT, max_archive_bytes=64*1024**2):
    """Snapshot, reimport and natively rebind untrusted saved boundary data.

    Digests establish byte identity; the fresh CAD/node/face/facet/exact-shell
    checks establish boundary admission. A mismatched build, kernel, topology,
    sizing or malformed input refuses rather than remeshing or trusting flags.
    This first bundle version requires the exact same extension binary/build;
    portable archive encoding does not make a new build's CAD checks equivalent.
    """
    began = time.perf_counter()
    directory = Path(directory).absolute()
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError("Boundary bundle must be a regular local directory")
    source_cap, archive_cap = integer(max_source_bytes, "max_source_bytes"), integer(max_archive_bytes, "max_archive_bytes")
    if source_cap > _SOURCE_LIMIT or archive_cap > _ARCHIVE_LIMIT:
        raise ValueError("Bundle read allowances exceed experimental ceilings")
    manifest = _json(_snapshot(directory / "manifest.json", _METADATA_LIMIT))
    required = {"schema", "units", "source_sha256", "boundary_sha256", "topology_sha256", "kernel_version", "build", "options", "files"}
    if set(manifest) != required or manifest["schema"] != _SCHEMA or manifest["units"] != "mm" or _json_bytes(manifest["build"]) != _json_bytes(_build()):
        raise ValueError("Boundary bundle schema/build/units differ")
    for name in ("source_sha256", "boundary_sha256", "topology_sha256"):
        digest = manifest[name]
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("Boundary identities must be lowercase SHA-256 strings")
    if not isinstance(manifest["files"], dict) or set(manifest["files"]) != set(_FILES):
        raise ValueError("Boundary bundle must use the fixed local file names")
    options = _options(**manifest["options"])
    if options != manifest["options"]:
        raise ValueError("Boundary options are not the exact normalized contract")
    payloads = {}
    caps = (source_cap, min(archive_cap, options["max_output_bytes"], options["max_working_bytes"]//8), _METADATA_LIMIT)
    for name, cap in zip(_FILES, caps, strict=True):
        entry = manifest["files"][name]
        if not isinstance(entry, dict) or set(entry) != {"sha256", "bytes"} or type(entry["bytes"]) is not int or not 0 < entry["bytes"] <= cap:
            raise ValueError("Boundary member metadata exceeds its allowance")
        data = _snapshot(directory / name, cap)
        if len(data) != entry["bytes"] or _hash(data) != entry["sha256"]:
            raise ValueError("Boundary member hash/length differs")
        payloads[name] = data
    source, archive = payloads["source.step"], payloads["boundary.paacsmesh"]
    if manifest["source_sha256"] != _hash(source) or manifest["boundary_sha256"] != _hash(archive):
        raise ValueError("Boundary/source identities differ")
    original_report = _json(payloads["boundary-report.json"])
    if (original_report.get("schema") != "paacs-experimental-boundary-report-v1"
            or original_report.get("prepared") is not True or original_report.get("size") != options["size"]
            or original_report.get("maximum_relative_volume_error") != options["max_relative_volume_error"]):
        raise ValueError("Saved boundary receipt differs from its frozen contract")
    with tempfile.TemporaryDirectory(prefix="paacs-boundary-import-") as temporary:
        source_path = Path(temporary) / "source.step"
        with source_path.open("xb") as stream:
            stream.write(source)
        model = read_step(source_path, max_file_bytes=source_cap, retain_source=True)
    if model.source_sha256 != manifest["source_sha256"] or model.kernel_version != manifest["kernel_version"] or _topology(model) != manifest["topology_sha256"]:
        raise ValueError("Reimported CAD source/kernel/topology identity differs")
    effective = number(original_report.get("effective_chord_error"), "saved effective chord error")
    if options["max_chord_error"] is None:
        body = model.inspect("bodies", limit=2)
        if body["total"] != 1 or len(body["entities"]) != 1:
            raise ValueError("Boundary reload requires one original solid")
        box = body["entities"][0]["bounds"]
        lower, upper = box["lower"], box["upper"]
        if len(lower) != 3 or len(upper) != 3:
            raise ValueError("Missing original body bounds for default chord admission")
        diagonal = math.dist([number(x,"body bound") for x in lower], [number(x,"body bound") for x in upper])
        ceiling = min(options["size"]*.05, diagonal*.0025)
    else:
        ceiling = options["max_chord_error"]
    if not 0 < effective <= ceiling:
        raise ValueError("Saved effective chord error exceeds its requested/default source bound")
    rebind_options = {**options, "max_chord_error": effective}
    boundary = _native._experimental_rebind_boundary(model._native, archive, rebind_options)
    rebound_archive = boundary.archive(min(_ARCHIVE_LIMIT, options["max_output_bytes"], options["max_working_bytes"]//8))
    if rebound_archive != archive:
        raise ValueError("Native rebind changed the frozen boundary bytes")
    report = boundary.report
    if any(report.get(key) != original_report.get(key) for key in ("nodes", "triangles", "solid_tag", "size", "effective_chord_error", "maximum_relative_volume_error")):
        raise ValueError("Fresh boundary admission differs from its saved receipt")
    report["saved_preparation"] = original_report
    metadata = {key: manifest[key] for key in ("source_sha256", "boundary_sha256", "topology_sha256", "kernel_version", "build")}
    report["provenance"] = metadata
    report["python_load_seconds"] = time.perf_counter() - began
    report["timing_scope"] = "Load snapshot includes file admission/hashing, fresh STEP import, fingerprint and Native boundary revalidation; precedes final metadata encoding/return; no volume call"
    return PreparedBoundary(model, boundary, source, archive, options, metadata, report)
