"""Independent archive measurements for installed-wheel sizing checks."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from collections import Counter
from dataclasses import dataclass
import math
from pathlib import Path
import statistics
import struct


TET_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
TET_FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))


@dataclass(frozen=True)
class Cell:
    kind: int
    dimension: int
    tag: int
    nodes: tuple


@dataclass(frozen=True)
class Archive:
    points: tuple
    node_tags: tuple
    cells: tuple


def read_archive(path):
    """Decode the documented archive format, including cell connectivity."""
    data = Path(path).read_bytes()
    magic, version, flags, nodes, cells, connectivity = struct.unpack_from("<8sIIQQQ", data)
    if magic != b"PAACSM8\0" or version != 1 or flags != 0:
        raise ValueError("Unexpected mesh archive header")
    points, node_tags, records = [], [], []
    at = 40
    for _ in range(nodes):
        x, y, z, dimension, tag = struct.unpack_from("<dddii", data, at)
        if not all(math.isfinite(value) for value in (x, y, z)):
            raise ValueError("Nonfinite archived coordinate")
        points.append((x, y, z))
        node_tags.append((dimension, tag))
        at += 32
    widths = (2, 3, 3, 6, 4, 4, 10, 6, 8, 5)
    observed_connectivity = 0
    for _ in range(cells):
        kind, dimension, tag = struct.unpack_from("<B3xii", data, at)
        if kind >= len(widths):
            raise ValueError("Unexpected archived cell type")
        at += 12
        width = widths[kind]
        indices = struct.unpack_from(f"<{width}q", data, at)
        if not all(0 <= index < nodes for index in indices):
            raise ValueError("Invalid archived connectivity")
        records.append(Cell(kind, dimension, tag, indices))
        observed_connectivity += width
        at += 8 * width
    if at != len(data) or observed_connectivity != connectivity:
        raise ValueError("Archive count or byte-length mismatch")
    return Archive(tuple(points), tuple(node_tags), tuple(records))


def quantile(values, fraction):
    ordered = sorted(values)
    if not ordered:
        return None
    position = (len(ordered) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def stats(values):
    return {"count": len(values), "minimum": min(values) if values else None,
            "median": quantile(values, .5), "p95": quantile(values, .95),
            "maximum": max(values) if values else None}


def measure(archive, face_tag, lower, upper, global_size, box_target):
    """Measure achieved lengths; these are not strict sizing certificates.

    A tetrahedron belongs to the measurement box when its centroid is inside.
    Its characteristic size is the median of its six edge lengths. Adjacent
    ratios compare that size across a shared face. The minimum-size clamp and
    distance-growth setting do not impose hard bounds on these measurements.
    """
    points = archive.points
    selected_edges = set()
    tetrahedra, sizes, centers = [], [], []
    for cell in archive.cells:
        if cell.kind == 2 and cell.dimension == 2 and cell.tag == face_tag:
            selected_edges.update(tuple(sorted((cell.nodes[i], cell.nodes[(i + 1) % 3])))
                                  for i in range(3))
        if cell.kind == 5:
            xyz = [points[index] for index in cell.nodes]
            lengths = [math.dist(xyz[a], xyz[b]) for a, b in TET_EDGES]
            if min(lengths) <= 0:
                raise ValueError("Zero-length tetrahedron edge")
            tetrahedra.append(cell.nodes)
            sizes.append(statistics.median(lengths))
            centers.append(tuple(sum(point[axis] for point in xyz) / 4 for axis in range(3)))
    selected_lengths = [math.dist(points[a], points[b]) for a, b in selected_edges]
    inside, outside, near = [], [], []
    for center, size in zip(centers, sizes):
        in_box = all(lo <= value <= hi for value, lo, hi in zip(center, lower, upper))
        (inside if in_box else outside).append(size)
        distance = math.sqrt(sum(max(0, lo - value, value - hi) ** 2
                                 for value, lo, hi in zip(center, lower, upper)))
        near.append(distance <= global_size * .5)
    uses, adjacent, near_adjacent = {}, [], []
    for index, nodes in enumerate(tetrahedra):
        for face in TET_FACES:
            key = tuple(sorted(nodes[k] for k in face))
            if key not in uses:
                uses[key] = index
                continue
            previous = uses[key]
            ratio = max(sizes[index], sizes[previous]) / min(sizes[index], sizes[previous])
            adjacent.append(ratio)
            if near[index] or near[previous]:
                near_adjacent.append(ratio)
    tags = Counter(f"{cell.dimension}:{cell.tag}" for cell in archive.cells)
    return {
        "nodes": len(points), "tet4": len(tetrahedra),
        "selected_face_edges_mm": stats(selected_lengths),
        "box_tet_median_edges_mm": stats(inside),
        "outside_box_tet_median_edges_mm": stats(outside),
        "volume_tet_median_edges_mm": stats(sizes),
        "adjacent_tet_size_ratio": stats(adjacent),
        "near_box_adjacent_tet_size_ratio": stats(near_adjacent),
        "box_reference_target_mm": box_target,
        "box_tet_size_to_reference_ratio": stats([size / box_target for size in inside]),
        "box_tets_above_reference_target": sum(size > box_target for size in inside),
        "box_tets_above_twice_reference_target": sum(size > 2 * box_target for size in inside),
        "tags": dict(sorted(tags.items())),
    }
