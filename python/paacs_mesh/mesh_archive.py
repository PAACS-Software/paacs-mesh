"""Bounded, read-only access to archive v1; no native validity certification."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from array import array
import hashlib
import json
import math
import mmap
import os
from pathlib import Path
import shutil
import struct
import tempfile
import time
from contextlib import contextmanager

from ._cad_common import MAX_JSON_BYTES, integer, json_copy, json_text, number


CELL_NAMES = ("line2", "line3", "tri3", "tri6", "quad4", "tet4", "tet10", "prism6", "hex8", "pyramid5")
CELL_WIDTHS = (2, 3, 3, 6, 4, 4, 10, 6, 8, 5)
HEADER = struct.Struct("<8sIIQQQ")
NODE = struct.Struct("<3dii")
CELL = struct.Struct("<B3sii")


class InspectionError(ValueError):
    """A bounded inspection was refused; no mesh validity gate is changed."""

    def __init__(self, reason, message, *, details=None):
        self.stage, self.reason = "mesh_inspection", reason
        self.report = {"schema": "paacs-mesh-inspection-error-v1", "reason": reason,
                       "details": details or {}}
        super().__init__(message)


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8192), b""):
            digest.update(block)
    return digest.hexdigest()


def _metadata(path, *, max_length, available):
    """Admit JSON object/copy overhead before a size-bounded immutable read."""
    size = Path(path).stat().st_size
    reservation = 64 * size + 65536
    if size > max_length or reservation > available:
        raise InspectionError("memory_budget", "Metadata parsing exceeds its byte allowance")
    with Path(path).open("rb") as stream:
        payload = stream.read(size + 1)
    if len(payload) != size:
        raise InspectionError("metadata_changed", "Metadata changed during its snapshot read")
    def unique_pairs(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON field {key!r}")
            result[key] = value
        return result
    def reject(value):
        raise ValueError(f"Nonfinite JSON constant {value}")
    try:
        result = json.loads(payload.decode("utf-8"), object_pairs_hook=unique_pairs, parse_constant=reject)
        if not isinstance(result, dict):
            raise ValueError("Metadata must be a JSON object")
        json_text(result)
    except (ValueError, RecursionError) as error:
        raise InspectionError("bundle_metadata", str(error)) from error
    return result, payload, reservation


def atomic_bytes(path, payload, *, protected=()):
    """Replace a complete export atomically, protecting all input files."""
    destination = Path(path).absolute()
    if not destination.name or not destination.parent.is_dir() or destination.is_dir():
        raise ValueError("Output must name a file in an existing directory")
    for original in protected:
        original = Path(original)
        if destination.resolve() == original.resolve() or (
                destination.exists() and original.exists() and os.path.samefile(destination, original)):
            raise ValueError("Inspection exports cannot replace their input files")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=".paacs-inspection-", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return destination


@contextmanager
def atomic_stream(path, *, protected=()):
    """Bounded writers can stream into a sibling file before atomic replacement."""
    destination = Path(path).absolute()
    if not destination.parent.is_dir() or destination.is_dir():
        raise ValueError("Output must name a file in an existing directory")
    for original in protected:
        original = Path(original)
        if destination.resolve() == original.resolve() or (destination.exists() and original.exists()
                and os.path.samefile(destination, original)):
            raise ValueError("Inspection exports cannot replace input files")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=".paacs-inspection-", delete=False) as stream:
            temporary = Path(stream.name)
            yield stream
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


class ArchiveView:
    """Immutable file snapshot with bounded offsets, coordinates and connectivity.

    The input is copied before mapping so changes to the caller's file cannot
    invalidate offsets or silently change IDs. Limits cover encoded bytes and
    owned storage, not the Python interpreter or an already existing native mesh.
    """

    def __init__(self, source, *, max_bytes=256 * 1024 * 1024, max_cells=1_000_000,
                 max_nodes=1_000_000, max_work=20_000_000, max_seconds=60):
        self._started = time.perf_counter()
        self.max_seconds = number(max_seconds, "max_seconds", minimum=0)
        if self.max_seconds <= 0:
            raise ValueError("max_seconds must be positive")
        self.max_bytes = integer(max_bytes, "max_bytes")
        if self.max_bytes < 65536:
            raise ValueError("max_bytes must be at least 65536 for bounded parser/hash overhead")
        self.max_cells = integer(max_cells, "max_cells")
        self.max_nodes = integer(max_nodes, "max_nodes")
        self.max_work = integer(max_work, "max_work")
        self._temporary = tempfile.TemporaryDirectory(prefix="paacs-inspection-")
        self._stream = self._mapped = None
        self.original_path = None
        self.protected_paths = []
        self.report = None
        self.report_basis = "unavailable"
        self._expected_mesh_sha256 = None
        self.source_sha256 = self.units = None
        self.metadata_bytes = 0
        try:
            self._load_source(source)
            self._parse()
            self.mesh_sha256 = file_hash(self.path)
            if self._expected_mesh_sha256 is not None and self.mesh_sha256 != self._expected_mesh_sha256:
                raise InspectionError("bundle_hash", "Immutable snapshot differs from the bundle mesh hash")
            self._deadline()
        except BaseException:
            self.close()
            raise

    def _load_source(self, source):
        from ._native import Mesh

        self.path = Path(self._temporary.name) / "mesh.paacsmesh"
        if isinstance(source, Mesh):
            if source.node_count > self.max_nodes or source.cell_count > self.max_cells:
                raise InspectionError("count_budget", "Native mesh counts exceed inspection limits")
            if source.node_count + 11 * source.cell_count > self.max_work:
                raise InspectionError("work_budget", "Native snapshot exceeds its conservative traversal allowance")
            # Every supported cell has at most ten connectivity indices. Bound
            # the native serialization before it allocates its encoded buffer.
            estimate = 40 + 32 * source.node_count + 92 * source.cell_count
            if estimate * 2 > self.max_bytes:
                raise InspectionError("memory_budget", "Native snapshot exceeds max_bytes",
                                      details={"snapshot_upper_bound": estimate})
            if 96 * source.cell_count > self.max_bytes:
                raise InspectionError("memory_budget", "Native body-tag discovery exceeds its conservative allowance")
            body_count = len(source.body_tags)
            report_reservation = 16 * 1024 * 1024 + 128 * 1024 * body_count
            if report_reservation + 2 * estimate > self.max_bytes:
                raise InspectionError("memory_budget", "Native report construction exceeds its conservative allowance")
            self.report = json_copy(source.report)
            self.metadata_bytes = 64 * len(json.dumps(self.report, ensure_ascii=True).encode("utf-8")) + 65536
            self.report_basis = "live_native_mesh_report"
            self.source_sha256 = self.report.get("source_sha256")
            self.units = self.report.get("units")
            if self.report.get("source_path"):
                self.protected_paths.append(Path(self.report["source_path"]))
            source.write_archive(self.path)
            self._deadline()
            return
        original = Path(source).absolute()
        if original.is_dir():
            manifest_path = original / "inspection.json"
            manifest, _, manifest_reservation = _metadata(manifest_path, max_length=65536,
                                                         available=self.max_bytes // 4)
            if manifest.get("schema") != "paacs-mesh-inspection-bundle-v1":
                raise InspectionError("bundle_schema", "Unsupported inspection bundle")
            mesh_path, report_path = original / "mesh.paacsmesh", original / "report.json"
            if mesh_path.is_symlink() or report_path.is_symlink():
                raise InspectionError("bundle_source", "Bundle inputs must be regular local files")
            if mesh_path.stat().st_size > self.max_bytes // 2 or report_path.stat().st_size > MAX_JSON_BYTES:
                raise InspectionError("memory_budget", "Bundle inputs exceed their byte budgets")
            self._expected_mesh_sha256 = manifest.get("mesh_sha256")
            if (not isinstance(self._expected_mesh_sha256, str) or len(self._expected_mesh_sha256) != 64
                    or any(char not in "0123456789abcdef" for char in self._expected_mesh_sha256)):
                raise InspectionError("bundle_hash", "Bundle requires a lowercase SHA-256 mesh identity")
            report, report_bytes, self.metadata_bytes = _metadata(report_path, max_length=MAX_JSON_BYTES,
                available=self.max_bytes - manifest_reservation - 2 * mesh_path.stat().st_size)
            if hashlib.sha256(report_bytes).hexdigest() != manifest.get("report_sha256"):
                raise InspectionError("bundle_hash", "Bundle report hash differs")
            self.report = report
            self.source_sha256 = self.report.get("source_sha256")
            self.units = self.report.get("units")
            if manifest.get("source_sha256") != self.source_sha256 or manifest.get("units") != self.units:
                raise InspectionError("bundle_metadata", "Bundle provenance differs from its report")
            self.report_basis = "hash_paired_saved_native_report_not_recertified"
            self.protected_paths.extend((manifest_path, mesh_path, report_path))
            if self.report.get("source_path"):
                self.protected_paths.append(Path(self.report["source_path"]))
            original = mesh_path
        self.original_path = original
        self.protected_paths.append(original)
        if not original.is_file():
            raise InspectionError("archive_input", "Input must be an archive or inspection bundle")
        size = original.stat().st_size
        if size > self.max_bytes // 2:
            raise InspectionError("memory_budget", "Encoded input exceeds the snapshot budget")
        # Read no more than the admitted size even if a file grows during copy.
        with original.open("rb") as incoming, self.path.open("xb") as outgoing:
            remaining = size
            while remaining:
                self._deadline()
                block = incoming.read(min(1024 * 1024, remaining))
                if not block:
                    raise InspectionError("archive_changed", "Archive changed while being copied")
                outgoing.write(block)
                remaining -= len(block)
            if incoming.read(1):
                raise InspectionError("archive_changed", "Archive grew while being copied")

    def _parse(self):
        self._stream = self.path.open("rb")
        size = self.path.stat().st_size
        if size < HEADER.size:
            raise InspectionError("archive_header", "Truncated archive header")
        self._mapped = mmap.mmap(self._stream.fileno(), 0, access=mmap.ACCESS_READ)
        self.protected_paths.append(self.path)
        magic, version, reserved, self.node_count, self.cell_count, connectivity = HEADER.unpack_from(self._mapped)
        if (magic, version, reserved) != (b"PAACSM8\0", 1, 0):
            raise InspectionError("archive_header", "Unsupported archive header")
        if self.node_count > self.max_nodes or self.cell_count > self.max_cells:
            raise InspectionError("count_budget", "Archive counts exceed inspection limits")
        expected = 40 + 32 * self.node_count + 12 * self.cell_count + 8 * connectivity
        if expected != size or max(self.node_count, self.cell_count, connectivity) > 2**63 - 1:
            raise InspectionError("archive_length", "Archive counts do not match its exact length")
        self.storage_bytes = size + 8 * self.cell_count + self.metadata_bytes + 16384
        if self.storage_bytes > self.max_bytes:
            raise InspectionError("memory_budget", "Archive plus cell offsets exceeds max_bytes")
        if self.node_count + self.cell_count + connectivity > self.max_work:
            raise InspectionError("work_budget", "Archive validation exceeds max_work")
        self.offsets = array("Q", [0]) * self.cell_count
        self.type_counts = [0] * 10
        for node in range(self.node_count):
            if node % 256 == 0:
                self._deadline()
            x, y, z, dim, tag = NODE.unpack_from(self._mapped, 40 + node * 32)
            if not all(math.isfinite(v) for v in (x, y, z)) or not 0 <= dim <= 3 or tag < 0:
                raise InspectionError("archive_node", f"Invalid coordinate/classification at node {node}")
        offset, actual_connectivity = 40 + self.node_count * 32, 0
        for cell_id in range(self.cell_count):
            if cell_id % 256 == 0:
                self._deadline()
            if offset + 12 > size:
                raise InspectionError("archive_cell", "Truncated cell record")
            kind, reserved, dim, tag = CELL.unpack_from(self._mapped, offset)
            if kind >= 10 or reserved != b"\0\0\0" or not 0 <= dim <= 3 or tag < 0:
                raise InspectionError("archive_cell", f"Invalid cell header at {cell_id}")
            width = CELL_WIDTHS[kind]
            if offset + 12 + width * 8 > size:
                raise InspectionError("archive_cell", "Truncated connectivity")
            nodes = struct.unpack_from("<" + "q" * width, self._mapped, offset + 12)
            if any(n < 0 or n >= self.node_count for n in nodes):
                raise InspectionError("archive_connectivity", f"Invalid node index at cell {cell_id}")
            self.offsets[cell_id] = offset
            self.type_counts[kind] += 1
            actual_connectivity += width
            offset += 12 + width * 8
        if actual_connectivity != connectivity or offset != size:
            raise InspectionError("archive_length", "Connectivity totals differ")
        if self.report is not None and (self.report.get("nodes") != self.node_count or self.report.get("cells") != self.cell_count):
            raise InspectionError("report_counts", "Report counts do not match the snapshot")

    def _deadline(self):
        if time.perf_counter() - self._started > self.max_seconds:
            raise InspectionError("time_budget", "Snapshot exceeds cooperative max_seconds")

    def node(self, node_id):
        if self._mapped is None:
            raise InspectionError("closed", "Inspection snapshot is closed")
        node_id = integer(node_id, "node_id", minimum=0)
        if node_id >= self.node_count:
            raise ValueError("node_id is outside this mesh")
        x, y, z, dim, tag = NODE.unpack_from(self._mapped, 40 + 32 * node_id)
        return [x, y, z], {"dimension": dim, "tag": tag}

    def cell(self, cell_id):
        if self._mapped is None:
            raise InspectionError("closed", "Inspection snapshot is closed")
        cell_id = integer(cell_id, "cell_id", minimum=0)
        if cell_id >= self.cell_count:
            raise ValueError("cell_id is outside this mesh")
        kind, _, dim, tag = CELL.unpack_from(self._mapped, self.offsets[cell_id])
        nodes = struct.unpack_from("<" + "q" * CELL_WIDTHS[kind], self._mapped, self.offsets[cell_id] + 12)
        return kind, nodes, {"dimension": dim, "tag": tag}

    def close(self):
        if self._mapped is not None:
            self._mapped.close()
            self._mapped = None
        if self._stream is not None:
            self._stream.close()
            self._stream = None
        self._temporary.cleanup()

    def __del__(self):
        try:
            self.close()
        except (AttributeError, OSError):
            pass
