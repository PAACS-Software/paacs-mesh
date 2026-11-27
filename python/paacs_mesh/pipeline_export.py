"""Tagged surface/CAD-edge artifacts; these are not qualified volume meshes.

Single-file writes preserve native mesh data. ``save`` additionally retains
units, exact CAD incidence, oriented uses, and source/boundary provenance in a
hash-linked directory bundle. No STEP data is copied into an export bundle.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import math
import os
from pathlib import Path
import shutil
import struct
import tempfile
import time

from ._cad_common import integer, json_copy, json_text, number, protected_destination


_HEADER = struct.Struct("<8sIIQQQ")
_NODE = struct.Struct("<3dii")
_CELL = struct.Struct("<B3sii")
_METADATA_CAP = 1024 * 1024
_DEFAULT_BYTES = 128 * 1024 * 1024


def _json_bytes(value):
    payload = (json_text(value) + "\n").encode("utf-8")
    if len(payload) > _METADATA_CAP:
        raise ValueError("Pipeline artifact metadata exceeds its 1 MiB allowance")
    return payload


def _digest(value):
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError("Pipeline identities must be lowercase SHA-256 strings")
    return value


class _Budget:
    def __init__(self, max_bytes, max_work, max_seconds):
        self.maximum = integer(max_bytes, "max_output_bytes")
        self.remaining = integer(max_work, "max_work")
        self.seconds = number(max_seconds, "max_seconds", minimum=0)
        if self.seconds <= 0:
            raise ValueError("max_seconds must be positive")
        self.started = time.perf_counter()

    def tick(self, work=1):
        if work > self.remaining:
            raise ValueError("Pipeline artifact work allowance exceeded")
        self.remaining -= work
        if time.perf_counter() - self.started > self.seconds:
            raise ValueError("Pipeline artifact time allowance exceeded")


class GeometryArtifact:
    """Owned Tri3 surface or Line2 CAD-edge view of a prepared boundary.

    Node and cell IDs are zero-based native row indices. CAD-edge artifacts
    deliberately retain *all* original boundary nodes, including unused nodes;
    Line2 cells follow each source curve's increasing native parameter. They
    are not every triangle edge. Face wire orientation and shared body/face
    incidence remain explicit in ``report`` and the saved metadata.

    Work/byte/time allowances cover this Python export and native writer's
    admission, not the pre-existing mesh or an OS-enforced process limit.
    Time checks are cooperative between native writer calls and bounded loops.
    """

    def __init__(self, native_handle, metadata, *, protected_paths=()):
        if not isinstance(metadata, dict) or metadata.get("kind") not in ("surface", "cad_edges"):
            raise ValueError("Artifact kind must be surface or cad_edges")
        if metadata.get("units") != "mm":
            raise ValueError("Prepared CAD artifacts use millimetres")
        for key in ("source_sha256", "boundary_sha256"):
            _digest(metadata.get(key))
        _json_bytes(metadata)
        self._metadata = json_copy(metadata)
        report = native_handle.report
        _json_bytes(report)
        if (metadata["kind"] == "surface" and report.get("prepared") is not True) or (
                metadata["kind"] == "cad_edges" and report.get("complete") is not True):
            raise ValueError("Pipeline artifact requires a completed native stage")
        self._native, self._report = native_handle, json_copy(report)
        self._protected = tuple(Path(path) for path in protected_paths)

    @property
    def kind(self):
        return self._metadata["kind"]

    @property
    def units(self):
        return self._metadata["units"]

    @property
    def source_sha256(self):
        return self._metadata["source_sha256"]

    @property
    def boundary_sha256(self):
        return self._metadata["boundary_sha256"]

    @property
    def metadata(self):
        return json_copy(self._metadata)

    @property
    def report(self):
        return json_copy(self._report)

    def _snapshot(self, budget):
        budget.tick()
        payload = self._native.archive(budget.maximum)
        if not isinstance(payload, bytes) or len(payload) > budget.maximum or len(payload) < _HEADER.size:
            raise ValueError("Pipeline archive exceeds its encoded byte allowance")
        magic, version, reserved, nodes, cells, connectivity = _HEADER.unpack_from(payload)
        kind, width, dimension = (2, 3, 2) if self.kind == "surface" else (0, 2, 1)
        cell_bytes = _CELL.size + width * 8
        if (magic != b"PAACSM8\0" or version != 1 or reserved != 0 or not nodes or not cells
                or connectivity != cells * width
                or _HEADER.size + nodes * _NODE.size + cells * cell_bytes != len(payload)):
            raise ValueError("Pipeline archive header/type/count contract differs")
        budget.tick(nodes + cells * (width + 1))
        for index in range(nodes):
            x, y, z, node_dimension, tag = _NODE.unpack_from(payload, _HEADER.size + index * _NODE.size)
            if not all(math.isfinite(p) for p in (x, y, z)) or node_dimension not in (0, 1, 2) or tag <= 0:
                raise ValueError("Pipeline archive contains invalid CAD node data")
            if index % 1024 == 0:
                budget.tick(0)
        offset = _HEADER.size + nodes * _NODE.size
        indices = struct.Struct("<" + "Q" * width)
        for index in range(cells):
            cell_kind, padding, cell_dimension, tag = _CELL.unpack_from(payload, offset)
            ids = indices.unpack_from(payload, offset + _CELL.size)
            if (cell_kind != kind or padding != b"\0\0\0" or cell_dimension != dimension
                    or tag <= 0 or max(ids) >= nodes or len(set(ids)) != width):
                raise ValueError("Pipeline archive contains unexpected cells or connectivity")
            offset += cell_bytes
            if index % 1024 == 0:
                budget.tick(0)
        if self.kind == "surface" and hashlib.sha256(payload).hexdigest() != self.boundary_sha256:
            raise ValueError("Surface archive differs from the frozen boundary identity")
        budget.tick()
        return payload

    def archive(self, *, max_bytes=_DEFAULT_BYTES, max_work=20_000_000, max_seconds=60):
        """Return immutable native archive bytes, retaining original node IDs.

        The raw archive has CAD classes but no units/provenance fields. Use
        ``save`` when those accompanying identities must travel with the mesh.
        """
        return self._snapshot(_Budget(max_bytes, max_work, max_seconds))

    def to_numpy(self, *, max_bytes=_DEFAULT_BYTES):
        """Return owned native arrays plus copied stage/provenance metadata.

        NumPy is optional. Arrays keep every node and original native cell
        index; no node filtering, renumbering or CAD association is inferred.
        """
        result = self._native.to_numpy(max_bytes=integer(max_bytes, "max_bytes"))
        result["artifact_metadata"] = self.metadata
        result["stage_report"] = self.report
        result["volume_qualified"] = False
        return result

    def _destination(self, path):
        path = Path(path).absolute()
        if path.suffix.lower() not in (".vtk", ".paacsmesh"):
            raise ValueError("Pipeline exports support .vtk and .paacsmesh")
        for source in self._protected:
            protected_destination(path, source)
        return path

    def write(self, path, *, max_output_bytes=_DEFAULT_BYTES, max_work=20_000_000, max_seconds=60):
        """Atomically write one native tagged mesh file; return its path.

        VTK/archive output keeps full coordinate precision and CAD node/cell
        classes. Native zero-based IDs are their row order. A single file does
        not carry the complete metadata; ``save`` publishes a provenance bundle.
        Existing output is replaced only after a complete native write.
        """
        budget = _Budget(max_output_bytes, max_work, max_seconds)
        destination = self._destination(path)
        self._snapshot(budget)
        self._native.write(destination, budget.maximum)
        # Native writer owns its transactional byte admission. No late refusal
        # after it has committed a single-file output.
        return destination

    def save(self, directory, *, max_output_bytes=2 * _DEFAULT_BYTES, max_work=20_000_000, max_seconds=60):
        """Publish a new hash-linked directory bundle with archive, VTK and CAD metadata.

        Refusals remove only this call's private staging directory. Existing
        destinations are never replaced. The source STEP is not copied.
        ``volume_qualified=False`` is explicit, even if a later volume run
        succeeds. This operation is an export, not a volume validation.
        """
        budget = _Budget(max_output_bytes, max_work, max_seconds)
        destination = Path(directory).absolute()
        if destination.exists() or destination.is_symlink() or not destination.parent.is_dir():
            raise ValueError("Artifact bundle requires a new directory in an existing parent")
        payload = self._snapshot(budget)
        metadata = {**self.metadata, "schema": "paacs-pipeline-artifact-v1",
                    "volume_qualified": False, "node_id_space": "original_boundary_zero_based",
                    "cell_id_space": "native_zero_based_row_order",
                    "single_file_metadata_scope": "native coordinates/connectivity/CAD classes only"}
        members = {"metadata.json": _json_bytes(metadata), "report.json": _json_bytes(self.report),
                   "mesh.paacsmesh": payload}
        used = sum(len(value) for value in members.values())
        if used >= budget.maximum:
            raise ValueError("Artifact bundle exceeds its aggregate encoded byte allowance")
        staging = Path(tempfile.mkdtemp(dir=destination.parent, prefix=".paacs-artifact-"))
        try:
            for name, value in members.items():
                budget.tick()
                with (staging / name).open("xb") as stream:
                    stream.write(value)
                    stream.flush()
                    os.fsync(stream.fileno())
            self._native.write(staging / "mesh.vtk", budget.maximum - used)
            files = {}
            for name in (*members, "mesh.vtk"):
                file = staging / name
                size = file.stat().st_size
                if name == "mesh.vtk":
                    used += size
                if used > budget.maximum:
                    raise ValueError("Artifact bundle exceeds its aggregate encoded byte allowance")
                digest = hashlib.sha256()
                with file.open("rb") as stream:
                    while block := stream.read(65536):
                        budget.tick()
                        digest.update(block)
                files[name] = {"bytes": size, "sha256": digest.hexdigest()}
            manifest = _json_bytes({"schema": "paacs-pipeline-artifact-manifest-v1", "kind": self.kind,
                "units": self.units, "source_sha256": self.source_sha256,
                "boundary_sha256": self.boundary_sha256, "volume_qualified": False, "files": files})
            if used + len(manifest) > budget.maximum:
                raise ValueError("Artifact bundle manifest exceeds its aggregate encoded byte allowance")
            with (staging / "manifest.json").open("xb") as stream:
                stream.write(manifest)
                stream.flush()
                os.fsync(stream.fileno())
            budget.tick()
            if destination.exists() or destination.is_symlink():
                raise ValueError("Artifact destination appeared during staging")
            os.rename(staging, destination)
            return destination
        finally:
            if staging.exists():
                shutil.rmtree(staging)
