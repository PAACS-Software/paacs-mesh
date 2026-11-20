# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Audit the two-wall normal Prism6 sweep and solve its fully developed flow."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from mesh_visual_data import read_archive, sha256


def projected_cells(mesh, first_height, wall_rows):
    """Recover the rectangular cross section of the wall-normal wedge grid."""
    wedges = np.array([cell.ids for cell in mesh.cells if cell.kind == 7], dtype=np.int64)
    if not len(wedges) or any(cell.kind != 7 for cell in mesh.cells if cell.dim == 3):
        raise ValueError("Expected an all-Prism6 duct")
    points = mesh.points[wedges]
    lower, upper = points[:, :3], points[:, 3:]
    translation = upper - lower
    if (not np.all(translation[:, :, (0, 2)] == 0) or np.any(translation[:, :, 1] <= 0) or
            not np.all(translation[:, :, 1] == translation[:, :1, 1])):
        raise ValueError("A prism does not sweep normally to the opposing y walls")
    jacobian = np.einsum("ij,ij->i", np.cross(lower[:, 1]-lower[:, 0], lower[:, 2]-lower[:, 0]),
                         translation[:, 0])
    if np.any(jacobian <= 0):
        raise ValueError("Nonpositive Prism6 Jacobian")
    stations = [np.unique(mesh.points[:, axis]) for axis in range(3)]
    widths = [np.diff(values) for values in stations]
    left, right = points.min(axis=1), points.max(axis=1)
    starts = np.column_stack([np.searchsorted(values, left[:, axis]) for axis, values in enumerate(stations)])
    ends = np.column_stack([np.searchsorted(values, right[:, axis]) for axis, values in enumerate(stations)])
    if not np.all(ends == starts + 1):
        raise ValueError("A swept prism spans nonadjacent station planes")
    voxels, multiplicity = np.unique(starts, axis=0, return_counts=True)
    if len(voxels) != np.prod([len(values) for values in widths]) or np.any(multiplicity != 2):
        raise ValueError("Expected two prisms in every Cartesian block")
    group_faces = np.bincount([cell.tag for cell in mesh.cells if cell.dim == 2], minlength=7)
    if len(group_faces) != 7 or group_faces[0] or not np.all(group_faces[1:] > 0):
        raise ValueError("The six boundary groups were not preserved")
    exterior_areas = np.zeros(6)
    for cell in (cell for cell in mesh.cells if cell.dim == 2):
        coordinates = mesh.points[list(cell.ids)]
        axis = (cell.tag - 1) // 2
        side = (cell.tag - 1) % 2
        if not np.all(coordinates[:, axis] == stations[axis][0 if side == 0 else -1]):
            raise ValueError("A boundary group appears on the wrong duct plane")
        vectors = coordinates[1:] - coordinates[0]
        exterior_areas[cell.tag - 1] += np.linalg.norm(np.cross(vectors[:-1], vectors[1:]), axis=1).sum()/2
    lx, ly, lz = [np.ptp(values) for values in stations]
    np.testing.assert_allclose(exterior_areas, [ly*lz, ly*lz, lx*lz, lx*lz, lx*ly, lx*ly], rtol=1e-12)
    volume = jacobian.sum()/2
    np.testing.assert_allclose(volume, lx*ly*lz, rtol=1e-12)
    thickness = widths[1]
    np.testing.assert_allclose(thickness, thickness[::-1], atol=1e-12)
    np.testing.assert_allclose(thickness[:wall_rows], first_height*1.22**np.arange(wall_rows), atol=1e-12)
    maximum_growth = np.maximum(thickness[1:]/thickness[:-1], thickness[:-1]/thickness[1:]).max()
    if maximum_growth > 1.22 + 1e-12:
        raise ValueError("The achieved wall-normal growth exceeds 1.22")
    return stations, dict(prism6=len(wedges), volume_mm3=float(volume),
                          minimum_full_cell_jacobian_mm3=float(jacobian.min()),
                          boundary_faces=dict(zip(("inlet", "outlet", "y-", "y+", "z-", "z+"),
                                                  group_faces[1:].tolist())),
                          first_height_mm=float(thickness[0]), maximum_growth=float(maximum_growth),
                          wall_normal_widths_mm=thickness.tolist())


def solve_section(stations):
    # Axially invariant Prism6 shape functions reduce to bilinear functions on
    # these rectangles: each is swept in y and has a linear tangential trace in z.
    axial, y, z = stations
    side = np.ptp(y)
    np.testing.assert_allclose(np.ptp(z), side)
    y_normal, z_normal = (y - y[0])/side, (z - z[0])/side
    yy, zz = np.meshgrid(y_normal, z_normal, indexing="ij")
    nodes = np.column_stack((yy.ravel(), zz.ravel()))
    grid = np.arange(len(nodes)).reshape(yy.shape)
    quads = np.column_stack([grid[:-1, :-1].ravel(), grid[1:, :-1].ravel(),
                             grid[1:, 1:].ravel(), grid[:-1, 1:].ravel()])
    dimensions = nodes[quads[:, 2]] - nodes[quads[:, 0]]
    areas = dimensions.prod(axis=1)
    corners = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
    stiffness = np.zeros((len(quads), 4, 4))
    for a in (-1/np.sqrt(3), 1/np.sqrt(3)):
        for b in (-1/np.sqrt(3), 1/np.sqrt(3)):
            basis_gradient = corners * (1 + corners[:, ::-1] * [b, a])/4
            physical_gradient = basis_gradient/dimensions[:, None, :]*2
            stiffness += (physical_gradient @ physical_gradient.transpose(0, 2, 1)) * (areas[:, None, None]/4)
    row = np.broadcast_to(quads[:, :, None], stiffness.shape).ravel()
    col = np.broadcast_to(quads[:, None, :], stiffness.shape).ravel()
    matrix = coo_matrix((stiffness.ravel(), (row, col)), shape=(len(nodes), len(nodes))).tocsr()
    load = np.bincount(quads.ravel(), weights=np.repeat(areas/4, 4), minlength=len(nodes))
    interior = ((nodes > 0) & (nodes < 1)).all(axis=1)
    field = np.zeros(len(nodes))
    field[interior] = spsolve(matrix[interior][:, interior], load[interior])
    conductance = float(load @ field)
    np.testing.assert_allclose(field @ (matrix @ field), conductance, rtol=1e-11)
    np.testing.assert_allclose((matrix @ field - load)[interior], 0, atol=1e-12)
    odd = np.arange(1., 8192., 2.)
    exact = (1 - 192/np.pi**5 * np.sum(np.tanh(odd*np.pi/2)/odd**5))/12
    viscosity = 1e-3
    length_m = np.ptp(axial)*1e-3
    side_m = side*1e-3
    mean_speed = .01
    scale = viscosity*length_m*mean_speed/side_m**2
    numerical_drop = scale/conductance
    exact_drop = scale/exact
    return dict(section_quads=len(quads), nondimensional_conductance=conductance,
                pressure_drop_pa=float(numerical_drop), analytical_pressure_drop_pa=float(exact_drop),
                pressure_error_percent=float(100*(numerical_drop/exact_drop - 1)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--first-height", type=float, required=True, help="Requested first layer height in mm")
    parser.add_argument("--wall-rows", type=int, required=True, help="Requested geometric rows at each wall")
    args = parser.parse_args()
    if not np.isfinite(args.first_height) or args.first_height <= 0 or not 1 <= args.wall_rows <= 16:
        parser.error("Provide a positive first height and 1 to 16 wall rows")
    stations, geometry = projected_cells(read_archive(args.archive), args.first_height, args.wall_rows)
    result = dict(archive_sha256=sha256(args.archive), **geometry, **solve_section(stations))
    args.receipt.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
