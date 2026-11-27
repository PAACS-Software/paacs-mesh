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
import shutil
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
_POINT_SCHEMA = "paacs-experimental-initial-points-bundle-v1"
_POINT_REPORT_SCHEMA = "paacs-experimental-initial-points-report-v1"
_POINT_FILES = ("points.bin", "point-report.json")
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


def _save_bundle(directory, files, encoded_manifest, prefix):
    """Stage owned files beside a new destination and publish by directory rename."""
    destination = Path(directory).absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError("Bundle destination already exists")
    if not destination.parent.is_dir():
        raise ValueError("Bundle destination requires an existing parent directory")
    staging = Path(tempfile.mkdtemp(dir=destination.parent, prefix=prefix))
    try:
        for name, data in [*files.items(), ("manifest.json", encoded_manifest)]:
            with (staging/name).open("xb") as stream:
                stream.write(data);stream.flush();os.fsync(stream.fileno())
        if destination.exists() or destination.is_symlink():
            raise FileExistsError("Bundle destination appeared during staging")
        os.rename(staging, destination)
        return destination
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def _snapshot(path, maximum, *, allow_empty=False):
    """Snapshot bounded regular files, including descriptor checks after open."""
    path = Path(path)
    if path.is_symlink():
        raise ValueError("Boundary inputs must not be symbolic links")
    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0) | getattr(os, "O_NONBLOCK", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    with os.fdopen(descriptor, "rb") as stream:
        info = os.fstat(stream.fileno())
        if not stat.S_ISREG(info.st_mode) or not (0 if allow_empty else 1) <= info.st_size <= maximum:
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


def _topology_records(model):
    """Bounded explicit CAD incidence/orientation, without inferred ownership."""
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
    return result


def _topology(model):
    """Integer incidence/orientation fingerprint, without process-local tokens."""
    return _hash(_json_bytes(_topology_records(model)))


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


def _point_payload(payload, options):
    """Validate encoded ownership and finite coordinates before decoded copies."""
    if not isinstance(payload, bytes) or len(payload) % 24:
        raise ValueError("Initial points must be immutable IEEE64 triples")
    count = len(payload)//24
    if count > min(_POINT_LIMIT, options["max_nodes"]) or len(payload) > options["max_working_bytes"]//4:
        raise ValueError("Initial point count/storage exceeds allowance")
    if any(not math.isfinite(value) for point in struct.iter_unpack("<3d", payload) for value in point):
        raise ValueError("Initial point cloud is nonfinite")
    return count


def _copy_point_tuples(payload, maximum):
    maximum = integer(maximum, "max_bytes")
    if maximum < 65536 or len(payload)//24 > (maximum-65536)//256:
        raise ValueError("Initial point cloud exceeds decoded max_bytes")
    return tuple(struct.iter_unpack("<3d", payload))


def _generated_report(payload, report):
    if (not isinstance(report, dict) or report.get("schema") != _POINT_REPORT_SCHEMA
            or report.get("generated") is not True or report.get("reason") != "none"
            or report.get("initial_points_available") is not True
            or type(report.get("initial_points")) is not int or report["initial_points"] != len(payload)//24
            or report.get("initial_points_sha256") != _hash(payload)
            or report.get("boundary_unchanged") is not True
            or report.get("recovery_performed") is not False
            or report.get("volume_construction_performed") is not False
            or report.get("mesh_published", False) is not False
            or report.get("units", "mm") != "mm"):
        raise ValueError("Initial point cloud differs from its generation-stage receipt")
    _json_bytes(report)


class InitialPointCloud:
    """Immutable staged classical sites; this is not a completed volume mesh.

    Coordinates are little-endian IEEE64 millimetres. Tuple/NumPy extraction
    copies explicitly under a decoded-storage allowance. Saved clouds bind to
    the exact source, frozen boundary, build and resource/sizing contract.
    """
    __slots__ = ("_payload", "_options", "_metadata", "_report", "_protected_paths", "_sealed")

    def __init__(self, payload, options, metadata, report, *, protected_paths=()):
        options = _options(**options)
        _point_payload(payload, options)
        _generated_report(payload, report)
        if not isinstance(metadata, dict) or set(metadata) != {"source_sha256", "boundary_sha256", "build"}:
            raise ValueError("Initial cloud requires exact source/boundary/build provenance")
        for name in ("source_sha256", "boundary_sha256"):
            digest = metadata[name]
            if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
                raise ValueError("Initial cloud identities must be lowercase SHA-256 strings")
        if (not isinstance(metadata["build"], dict)
                or _json_bytes(report.get("provenance")) != _json_bytes(metadata)):
            raise ValueError("Initial cloud receipt provenance differs")
        self._payload, self._options = payload, _json_bytes(options)
        self._metadata, self._report = _json_bytes(metadata), _json_bytes(report)
        if not isinstance(protected_paths, (tuple, list)) or len(protected_paths) > 15:
            raise ValueError("Initial cloud accepts at most 15 protected inputs")
        self._protected_paths = tuple(dict.fromkeys(Path(path).absolute() for path in protected_paths))
        self._sealed = True

    def __setattr__(self, key, value):
        if getattr(self, "_sealed", False):
            raise AttributeError("InitialPointCloud is immutable")
        object.__setattr__(self, key, value)

    def __delattr__(self, key):
        raise AttributeError("InitialPointCloud is immutable")

    @property
    def report(self):
        return _json(self._report)

    @property
    def point_count(self):
        return len(self._payload)//24

    @property
    def sha256(self):
        return _hash(self._payload)

    @property
    def source_sha256(self):
        return _json(self._metadata)["source_sha256"]

    @property
    def boundary_sha256(self):
        return _json(self._metadata)["boundary_sha256"]

    @property
    def provenance(self):
        return _json(self._metadata)

    @property
    def options(self):
        return _json(self._options)

    @property
    def units(self):
        return "mm"

    def points(self, *, max_bytes=16*1024**2):
        """Copy immutable coordinate tuples; an available empty cloud returns ()."""
        return _copy_point_tuples(self._payload, max_bytes)

    def to_numpy(self, *, max_bytes=128*1024**2):
        """Return an independent, initially read-only float64 N×3 array."""
        maximum = integer(max_bytes, "max_bytes")
        if maximum < 65536 or len(self._payload) > maximum-65536:
            raise ValueError("Initial point array exceeds max_bytes")
        try:
            import numpy as np
        except ImportError as error:
            raise ImportError("Initial point arrays require NumPy; install 'paacs-mesh[numpy]'.") from error
        array = np.frombuffer(self._payload, dtype=np.dtype("<f8")).reshape((-1, 3)).copy()
        array.setflags(write=False)
        return array

    def _for_boundary(self, boundary, build, options):
        metadata = _json(self._metadata)
        if (_json_bytes(metadata) != _json_bytes({"source_sha256": boundary.source_sha256,
                         "boundary_sha256": boundary.boundary_sha256, "build": build})
                or self._options != _json_bytes(options)):
            raise ValueError("Initial point cloud source/boundary/build/options differ")
        return self._payload

    def save(self, directory):
        """Stage a bounded private cloud bundle, then rename to a new destination."""
        files = {"points.bin": self._payload, "point-report.json": self._report}
        manifest = {"schema": _POINT_SCHEMA, "units": "mm", **_json(self._metadata),
                    "options": _json(self._options), "point_count": self.point_count,
                    "files": {name: {"sha256": _hash(data), "bytes": len(data)} for name, data in files.items()}}
        return _save_bundle(directory, files, _json_bytes(manifest), ".paacs-points-")


class PreparedBoundary:
    """Owned exact-certified Tri3 boundary and frozen uniform sizing contract.

    The report certifies boundary preparation only, never volume filling or
    solver adequacy. Saving includes source CAD. Resource limits are modeled
    Native/encoded-buffer allowances, not a process RSS/hard timeout guarantee.
    """
    __slots__ = ("_model", "_boundary", "_source", "_archive", "_options", "_metadata", "_report", "_protected_paths")

    def __init__(self, model, boundary, source, archive, options, metadata, report, *, protected_paths=()):
        self._model, self._boundary = model, boundary
        self._source, self._archive = bytes(source), bytes(archive)
        self._options, self._metadata, self._report = _json_bytes(options), _json_bytes(metadata), _json_bytes(report)
        if not isinstance(protected_paths, (tuple, list)) or len(protected_paths) > 15:
            raise ValueError("Prepared boundary accepts at most 15 additional protected inputs")
        original = model.info.get("source_path")
        paths = ((original,) if original else ()) + tuple(protected_paths)
        self._protected_paths = tuple(dict.fromkeys(Path(path).absolute() for path in paths))

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

    def cad_topology(self, *, max_bytes=_METADATA_LIMIT):
        """Copy explicit CAD incidences and orientations under the frozen hash.

        Face owner lists and senses come from CAD inspection. Triangle edges or
        nearby coordinates are never substituted for CAD edges or ownership.
        """
        maximum = integer(max_bytes, "max_bytes")
        self._model._assert_identity()
        records = _topology_records(self._model)
        encoded = _json_bytes(records)
        if maximum <= 0 or len(encoded) > maximum:
            raise ValueError("CAD topology copy exceeds max_bytes")
        if _hash(encoded) != _json(self._metadata)["topology_sha256"]:
            raise ValueError("CAD topology differs from the prepared boundary")
        return _json(encoded)

    def surface(self):
        """Return the qualified tagged Tri3 surface as an exportable artifact."""
        from .pipeline_export import GeometryArtifact
        return GeometryArtifact(self._boundary, self._artifact_metadata("surface"),
                                protected_paths=self._protected_paths)

    def edges(self):
        """Return retained ordered CAD-edge chains, never all triangle edges."""
        from .pipeline_export import GeometryArtifact
        metadata = self._artifact_metadata("cad_edges")
        handle = _native._experimental_cad_edges(self._model._native, self._boundary, _json(self._options))
        return GeometryArtifact(handle, metadata, protected_paths=self._protected_paths)

    def _artifact_metadata(self, kind):
        return {"kind": kind, "units": "mm", "source_sha256": self.source_sha256,
                "boundary_sha256": self.boundary_sha256, "cad_topology": self.cad_topology(),
                "provenance": _json(self._metadata)}

    def generate_initial_points(self):
        """Generate classical sites before constraint recovery or volume fill.

        This performs the same initial placement and auxiliary topology work
        used by classical replay, then stops. It is a staged input artifact,
        not a completed/validated volume mesh or a learned point generator.
        """
        began = time.perf_counter()
        called = False
        try:
            self._model._assert_identity()
            build = _build();options = _json(self._options)
            called = True
            result = _native._experimental_generate_initial_points(self._model._native, self._boundary, options)
            if not isinstance(result, dict) or set(result) != {"points", "report"}:
                raise ValueError("Malformed initial-point generation response")
            payload, report = result["points"], result["report"]
            _point_payload(payload, options);_generated_report(payload, report)
            metadata = {"source_sha256": self.source_sha256, "boundary_sha256": self.boundary_sha256,
                        "build": build}
            report = {**report, "provenance": metadata, "units": "mm",
                      "python_generation_seconds": time.perf_counter()-began,
                      "timing_scope": "Generation includes original placement and auxiliary topology finalization; excludes constraint recovery, material fill and volume validation. Python snapshot includes provenance and encoded receipt checks; precedes immutable receipt encoding and return. Native spans are not an additive partition."}
            return InitialPointCloud(payload, options, metadata, report)
        except BaseException as error:
            report = getattr(error, "report", None)
            report = dict(report) if isinstance(report, dict) else {}
            report.update(python_generation_seconds=time.perf_counter()-began,
                          native_binding_called=called,
                          python_timing_scope="Failed generation snapshot includes Python admission and any invoked Native placement; no volume replay is requested by this operation.")
            error.report = report
            raise

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
            build = _build()
            extension_hash = build["native_extension_sha256"]
            options = _json(self._options)
            protected_paths = self._protected_paths
            if isinstance(points, InitialPointCloud):
                payload = points._for_boundary(self, build, options)
                protected_paths = tuple(dict.fromkeys((*protected_paths, *points._protected_paths)))
            else:
                payload = _point_bytes(points, options)
            if len(protected_paths) > 16:
                raise ValueError("Boundary and cloud exceed the protected-input path allowance")
            supplied_hash = None if payload is None else _hash(payload)
            packing = time.perf_counter() - began
            called = True
            mesh = _native._experimental_mesh_boundary(self._model._native, self._boundary, options, payload, protected_paths)
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
        """Stage complete files/manifest, then rename to a new bundle directory."""
        files = {"source.step": self._source, "boundary.paacsmesh": self._archive,
                 "boundary-report.json": self._report}
        manifest = {"schema": _SCHEMA, "units": "mm", **_json(self._metadata),
                    "options": _json(self._options),
                    "files": {name: {"sha256": _hash(data), "bytes": len(data)} for name, data in files.items()}}
        return _save_bundle(directory, files, _json_bytes(manifest), ".paacs-boundary-")


def initial_points(value, *, max_bytes=16*1024**2):
    """Copy a replay result/refusal's actual pre-constraint cloud, never a mesh.

    Refusals expose a separate at-most-64-MiB immutable IEEE64 payload. An
    unavailable capture refuses explicitly; a captured empty cloud returns ().
    ``max_bytes`` conservatively bounds decoded tuple storage as well as bytes.
    """
    maximum = integer(max_bytes, "max_bytes")
    if isinstance(value, InitialPointCloud):
        return value.points(max_bytes=maximum)
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


def load_initial_points(directory, *, max_bytes=32*1024**2):
    """Load a hash-bound staged cloud without claiming volume validity.

    No CAD or volume operation occurs here. Passing the returned cloud to a
    prepared boundary rechecks source/boundary/build/options identity and all
    ordinary Native supplied-point and publication requirements.
    """
    directory = Path(directory).absolute()
    maximum = integer(max_bytes, "max_bytes")
    if maximum <= 0 or maximum > 64*1024**2:
        raise ValueError("Initial cloud read allowance must be positive and at most 64 MiB")
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError("Initial cloud bundle must be a regular local directory")
    manifest = _json(_snapshot(directory/"manifest.json", _METADATA_LIMIT))
    keys = {"schema", "units", "source_sha256", "boundary_sha256", "build", "options", "point_count", "files"}
    if (set(manifest) != keys or manifest["schema"] != _POINT_SCHEMA or manifest["units"] != "mm"
            or _json_bytes(manifest["build"]) != _json_bytes(_build())):
        raise ValueError("Initial cloud bundle schema/build/units differ")
    for name in ("source_sha256", "boundary_sha256"):
        digest = manifest[name]
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("Initial cloud identities must be lowercase SHA-256 strings")
    options = _options(**manifest["options"])
    if options != manifest["options"] or type(manifest["point_count"]) is not int or not 0 <= manifest["point_count"] <= min(_POINT_LIMIT, options["max_nodes"]):
        raise ValueError("Initial cloud options/count differ from the bounded contract")
    if not isinstance(manifest["files"], dict) or set(manifest["files"]) != set(_POINT_FILES):
        raise ValueError("Initial cloud bundle must use the fixed local file names")
    payloads = {}
    for name, cap in (("points.bin", min(maximum, _POINT_LIMIT*24, options["max_working_bytes"]//4)),
                      ("point-report.json", _METADATA_LIMIT)):
        entry = manifest["files"][name]
        if (not isinstance(entry, dict) or set(entry) != {"sha256", "bytes"}
                or type(entry["bytes"]) is not int or not (0 if name == "points.bin" else 1) <= entry["bytes"] <= cap):
            raise ValueError("Initial cloud member metadata exceeds its allowance")
        payload = _snapshot(directory/name, cap, allow_empty=name == "points.bin")
        if len(payload) != entry["bytes"] or _hash(payload) != entry["sha256"]:
            raise ValueError("Initial cloud member hash/length differs")
        payloads[name] = payload
    if len(payloads["points.bin"]) != manifest["point_count"]*24:
        raise ValueError("Initial cloud count differs from encoded triples")
    metadata = {key: manifest[key] for key in ("source_sha256", "boundary_sha256", "build")}
    report = _json(payloads["point-report.json"])
    if _json_bytes(report.get("provenance")) != _json_bytes(metadata):
        raise ValueError("Initial cloud receipt provenance differs")
    return InitialPointCloud(payloads["points.bin"], options, metadata, report,
                             protected_paths=(directory/"manifest.json", *(directory/name for name in _POINT_FILES)))


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
    return PreparedBoundary(model, boundary, source, archive, options, metadata, report,
                            protected_paths=(directory/"manifest.json", *(directory/name for name in _FILES)))
