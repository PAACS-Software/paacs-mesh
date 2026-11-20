# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Independently check the archived interface, boundary groups and cell maps."""

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from mesh_visual_data import FACES, read_archive, sha256


WALLS = ('x-', 'x+', 'y-', 'y+')
EXTERIOR = WALLS + ('inlet', 'outlet')
TOL = 1e-10


def require(condition, message):
    if not condition:
        raise ValueError(message)


def plane(points, lo, hi):
    matches = [2*axis + side for axis in range(len(lo)) for side, position in enumerate((lo[axis], hi[axis]))
               if np.all(np.abs(points[:, axis] - position) <= TOL)]
    require(len(matches) == 1, 'Expected one containing face plane')
    return matches[0]


def orientation(ids):
    start = ids.index(min(ids))
    return ids[start:] + ids[:start]


def opposite(a, b):
    return orientation(a) == orientation(b[::-1])


def surface_summary(faces, points, side, lo, hi):
    """Check one oriented rectangle, including all perimeter edge intervals."""
    axis, direction = divmod(side, 2)
    in_plane = [i for i in range(3) if i != axis]
    low, high = lo[in_plane], hi[in_plane]
    edges = defaultdict(list)
    total_area = 0.
    for face in faces:
        p = points[list(face)]
        require(plane(p, lo, hi) == side, 'Face appears on the wrong rectangle')
        spans = p[1:] - p[0]
        area = .5*np.cross(spans[:-1], spans[1:]).sum(axis=0)[axis]*(2*direction - 1)
        require(area > 0, 'Reversed or degenerate rectangle facet')
        total_area += area
        for i, p0 in enumerate(face):
            p1 = face[(i+1) % len(face)]
            edges[tuple(sorted((p0, p1)))].append((p0, p1))

    perimeter = [[] for _ in range(4)]
    for uses in edges.values():
        if len(uses) == 2:
            require(uses[0] == uses[1][::-1], 'Inconsistent rectangle edge winding')
            continue
        require(len(uses) == 1, 'Nonmanifold rectangle edge')
        ends = points[list(uses[0])][:, in_plane]
        edge_side = plane(ends, low, high)
        perimeter[edge_side].append(tuple(np.sort(ends[:, 1-edge_side//2])))

    for edge_side, segments in enumerate(perimeter):
        along = 1 - edge_side//2
        at = low[along]
        for left, right in sorted(segments):
            require(abs(left-at) <= TOL and right > left, 'Gap or overlap in rectangle perimeter')
            at = right
        require(abs(at-high[along]) <= TOL, 'Incomplete rectangle perimeter')
    np.testing.assert_allclose(total_area, np.prod(high-low), rtol=1e-11)
    return float(total_area)


def wedge_measures(wedges):
    """Integrate the exact determinant polynomial and minimise it on the reference wedge."""
    bottom, top = wedges[:, :3], wedges[:, 3:]
    base = bottom[:, 1:] - bottom[:, :1]
    translation = top - bottom
    delta = translation[:, 1:] - translation[:, :1]
    cross0 = np.cross(base[:, 0], base[:, 1])
    cross1 = np.cross(base[:, 0], delta[:, 1]) + np.cross(delta[:, 0], base[:, 1])
    cross2 = np.cross(delta[:, 0], delta[:, 1])
    # The determinant is affine in each triangular section, so its minimum
    # at fixed thickness lies at a triangle vertex. Each of the three
    # resulting thickness polynomials is quadratic on the unit interval.
    c0, c1, c2 = [np.einsum('ij,ikj->ik', vector, translation) for vector in (cross0, cross1, cross2)]
    minima = np.minimum(c0, c0+c1+c2)
    t = np.divide(-c1, 2*c2, out=np.zeros_like(c1), where=c2 > 0)
    middle = (t > 0) & (t < 1)
    minima[middle] = np.minimum(minima[middle], (c0+t*(c1+t*c2))[middle])
    return minima.min(axis=1), (c0 + .5*c1 + c2/3).sum(axis=1)/6


def wall_sectors(wedges):
    """Identify each prism's wall and slab, allowing matching mitered corner faces."""
    candidates = []
    for axis, sign in ((0, -1), (0, 1), (1, -1), (1, 1)):
        normal = sign*wedges[:, :, axis]
        axis_caps = normal.reshape(-1, 2, 3)
        # Both triangular caps have a constant wall-normal coordinate.
        candidates.append(np.all(normal >= np.abs(wedges[:, :, 1-axis])-TOL, axis=1)
                          & np.all(np.ptp(axis_caps, axis=2) <= TOL, axis=1)
                          & (np.ptp(normal, axis=1) > TOL))
    candidates = np.column_stack(candidates)
    require(np.all(candidates.sum(axis=1) == 1), 'A prism crosses a wall sector or has nonparallel caps')
    return np.argmax(candidates, axis=1)


def audit(archive, first_height, growth, wall_layers, exact_input):
    mesh = read_archive(archive)
    points = mesh.points
    lo, hi = points.min(axis=0), points.max(axis=0)
    np.testing.assert_allclose(lo, [-5., -5., 0.], atol=TOL)
    np.testing.assert_allclose(hi, [5., 5., 40.], atol=TOL)
    require(len(np.unique(points, axis=0)) == len(points), 'Coincident nodes have separate identities')
    volumes = [cell for cell in mesh.cells if cell.dim == 3]
    require({cell.kind for cell in volumes} == {5, 7} and len({cell.tag for cell in volumes}) == 1,
            'Expected a single fluid region with both Prism6 and Tet4')
    declared = {}
    for cell in mesh.cells:
        if cell.dim != 2:
            continue
        key = tuple(sorted(cell.ids))
        require(cell.kind in (2, 4) and key not in declared, 'Invalid or duplicate exterior face')
        declared[key] = cell

    incidence = defaultdict(list)
    prism_nodes, tetra_nodes = [], []
    for cell in volumes:
        ids = cell.ids
        require(len(set(ids)) == len(ids), 'Repeated node in volume element')
        own = prism_nodes if cell.kind == 7 else tetra_nodes
        index = len(own)
        own.append(ids)
        for face in FACES[cell.kind]:
            oriented = tuple(ids[i] for i in face)
            incidence[tuple(sorted(oriented))].append((cell.kind, index, oriented))
    require(len({node for row in prism_nodes+tetra_nodes for node in row}) == len(points), 'Unused node')

    prism_ids = np.asarray(prism_nodes)
    prism_xyz = points[prism_ids]
    tet_xyz = points[np.asarray(tetra_nodes)]
    core_lo, core_hi = tet_xyz.min(axis=(0, 1)), tet_xyz.max(axis=(0, 1))
    thickness = first_height*growth**np.arange(wall_layers)
    core_half = 5 - thickness.sum()
    np.testing.assert_allclose(core_lo, [-core_half, -core_half, 0], atol=TOL)
    np.testing.assert_allclose(core_hi, [core_half, core_half, 40], atol=TOL)
    wedge_min, prism_volumes = wedge_measures(prism_xyz)
    tet_det = np.linalg.det(tet_xyz[:, 1:] - tet_xyz[:, :1])
    require(np.all(wedge_min > 0) and np.all(tet_det > 0), 'Nonpositive whole-cell Jacobian')
    core_volume = tet_det.sum()/6
    shell_volume = prism_volumes.sum()
    expected_core = (2*core_half)**2*40
    np.testing.assert_allclose([core_volume, shell_volume], [expected_core, 4000-expected_core], rtol=1e-11)
    sectors = wall_sectors(prism_xyz)
    prism_counts, first_heights = [], []
    expected = 5 - np.concatenate(([0.], np.cumsum(thickness)))
    for wall in range(4):
        own = prism_xyz[sectors == wall]
        own_ids = prism_ids[sectors == wall]
        axis, direction = divmod(wall, 2)
        coordinate = (2*direction - 1)*own[:, :, axis]
        require(len(coordinate) > 0, 'Missing wall prism sector')
        # Each constant-normal slice is a scaled copy of the same wall
        # triangulation. This verifies the homothety and its axial invariance
        # even for the bilinear diagonal faces between paired prisms.
        tangent = own[:, :, 1-axis]/coordinate
        np.testing.assert_allclose(tangent[:, :3], tangent[:, 3:], atol=TOL, rtol=0)
        np.testing.assert_allclose(own[:, :3, 2], own[:, 3:, 2], atol=TOL, rtol=0)
        outer, inner = coordinate.max(axis=1), coordinate.min(axis=1)
        layer = np.argmin(abs(outer[:, None]-expected[:-1]), axis=1)
        np.testing.assert_allclose(outer, expected[layer], atol=TOL, rtol=0)
        np.testing.assert_allclose(inner, expected[layer+1], atol=TOL, rtol=0)
        measured_volume = np.bincount(layer, weights=prism_volumes[sectors == wall], minlength=wall_layers)
        np.testing.assert_allclose(measured_volume, (expected[:-1]**2-expected[1:]**2)*40, rtol=1e-11)
        for depth, radius in enumerate(expected[:-1]):
            on_layer = layer == depth
            bottom = abs(coordinate[on_layer, 0]-radius) <= TOL
            caps = [tuple(int(ids[i]) for i in FACES[7][0 if first else 1])
                    for ids, first in zip(own_ids[on_layer], bottom)]
            surface_summary(caps, points, wall, np.array([-radius, -radius, 0.]), np.array([radius, radius, 40.]))
        prism_counts.append(np.bincount(layer, minlength=wall_layers).tolist())
        first_heights.append(float(np.min((outer-inner)[layer == 0])))

    exposed = [[] for _ in EXTERIOR]
    interface = [[] for _ in WALLS]
    core_boundary = []
    exterior_tags = [set() for _ in EXTERIOR]
    corners = Counter()
    for key, uses in incidence.items():
        if len(uses) == 1:
            kind, _, oriented = uses[0]
            require(key in declared and orientation(oriented) == orientation(declared[key].ids),
                    'Missing or misoriented exterior face')
            surface = plane(points[list(oriented)], lo, hi)
            exposed[surface].append(oriented)
            exterior_tags[surface].add(declared.pop(key).tag)
            if kind == 5:
                core_boundary.append(oriented)
            continue
        require(len(uses) == 2 and opposite(uses[0][2], uses[1][2]), 'Invalid interior face incidence')
        (kind_a, index_a, a), (kind_b, index_b, b) = uses
        if kind_a != kind_b:
            require(len(key) == 3, 'Expected a shared triangle at the Prism6–Tet4 transition')
            face = a if kind_a == 5 else b
            side = plane(points[list(face)], core_lo, core_hi)
            require(side < 4, 'Transition triangle lies outside the four core walls')
            interface[side].append(face)
            core_boundary.append(face)
        elif kind_a == 7 and sectors[index_a] != sectors[index_b]:
            require(len(key) == 4 and sectors[index_a]//2 != sectors[index_b]//2,
                    'Invalid corner seam contact')
            corners[tuple(sorted((int(sectors[index_a]), int(sectors[index_b]))))] += 1

    require(not declared, 'Exterior declaration contains faces unused by volume cells')
    require(len(corners) == 4 and len(set(corners.values())) == 1, 'Incomplete miter seams')
    require(exterior_tags == [{1}, {2}, {3}, {4}, {5}, {6}], 'Unexpected STEP boundary groups')
    outer_areas = [surface_summary(faces, points, side, lo, hi) for side, faces in enumerate(exposed)]
    transition_areas = [surface_summary(faces, points, side, core_lo, core_hi)
                        for side, faces in enumerate(interface)]
    if exact_input is not None:
        with exact_input.open('w', encoding='ascii') as stream:
            used = np.unique(tetra_nodes)
            stream.write(f'{len(used)} {len(core_boundary)} {len(tetra_nodes)}\n')
            np.savetxt(stream, points[used], fmt='%.17g')
            np.savetxt(stream, np.searchsorted(used, core_boundary), fmt='%d')
            np.savetxt(stream, np.searchsorted(used, tetra_nodes), fmt='%d')
    return dict(archive_sha256=sha256(archive), nodes=len(points), prism6=len(prism_nodes), tet4=len(tetra_nodes),
                prism_minimum_jacobian_mm3=float(wedge_min.min()), tet_minimum_jacobian_mm3=float(tet_det.min()),
                prism_volume_mm3=float(shell_volume), tet_volume_mm3=float(core_volume),
                first_layer_height_mm=first_height, layer_growth=growth, wall_normal_widths_mm=thickness.tolist(),
                first_layer_height_by_wall_mm=dict(zip(WALLS, first_heights)),
                prism_cells_by_wall_and_layer=dict(zip(WALLS, prism_counts)),
                interface_triangles=sum(map(len, interface)), interface_triangles_by_wall=dict(zip(WALLS, map(len, interface))),
                interface_areas_mm2=dict(zip(WALLS, transition_areas)),
                corner_shared_quad_faces=sum(corners.values()), corner_shared_quad_faces_each=next(iter(corners.values())),
                exterior_faces_by_side=dict(zip(EXTERIOR, map(len, exposed))),
                exterior_area_by_side_mm2=dict(zip(EXTERIOR, outer_areas)),
                exterior_tags=dict(zip(EXTERIOR, map(lambda tags: next(iter(tags)), exterior_tags))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('receipt', type=Path)
    parser.add_argument('--first-height', type=float, required=True)
    parser.add_argument('--growth', type=float, required=True)
    parser.add_argument('--wall-layers', type=int, required=True)
    parser.add_argument('--tet-core-exact-input', type=Path)
    args = parser.parse_args()
    if not (0 < args.first_height < 5 and 1 <= args.growth <= 2 and 1 <= args.wall_layers <= 16):
        parser.error('Height, growth or layer count is outside this duct study')
    result = audit(args.archive, args.first_height, args.growth, args.wall_layers, args.tet_core_exact_input)
    args.receipt.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(f"{result['prism6']:,} Prism6, {result['tet4']:,} Tet4 and {result['interface_triangles']:,} shared triangles")


if __name__ == '__main__':
    main()
