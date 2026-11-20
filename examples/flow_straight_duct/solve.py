# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Solve fully developed laminar flow on the section of a PAACS Prism6 duct."""

import argparse
import json
from pathlib import Path
import re
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.mesh_visual_data import read_archive, sha256


def row_keys(cells):
    return np.sort(cells, axis=1)


def ordered_rows(rows):
    return rows[np.lexsort(rows.T[::-1])]


def signed_double_area(xyz):
    a, b = xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]
    return a[:, 0]*b[:, 1] - a[:, 1]*b[:, 0]


def section_edges(xy, triangles):
    edges = triangles[:, ((0, 1), (1, 2), (2, 0))].reshape(-1, 2)
    _, first, inverse, count = np.unique(row_keys(edges), axis=0, return_index=True,
                                         return_inverse=True, return_counts=True)
    if np.any(count > 2):
        raise ValueError("The extracted section has nonmanifold edges")
    internal = count[inverse] == 2
    pair_order = np.argsort(inverse[internal], kind="stable")
    pairs = edges[internal][pair_order].reshape(-1, 2, 2)
    if not np.array_equal(pairs[:, 0], pairs[:, 1, ::-1]):
        raise ValueError("Interior section edges have inconsistent orientation")
    boundary = edges[first[count == 1]]
    wall_xy = xy[boundary]
    side_tests = [(abs(wall_xy[:, :, axis]-value) <= 1e-10).all(axis=1)
                  for axis in (0, 1) for value in (0, 1)]
    if not np.all(np.sum(side_tests, axis=0) == 1):
        raise ValueError("The section boundary must lie on the four square walls")
    lengths = np.linalg.norm(wall_xy[:, 1]-wall_xy[:, 0], axis=1)
    np.testing.assert_allclose(np.array(side_tests) @ lengths, 1., atol=1e-10, rtol=0)
    return boundary, first[count == 1]//3


def verify_exterior(points, caps, quads, columns, triangles, edges, bounds):
    cap_keys = np.concatenate([columns[i, triangles] for i in (0, -1)])
    if not np.array_equal(ordered_rows(row_keys(caps.ids)), ordered_rows(row_keys(cap_keys))):
        raise ValueError("The observed inlet/outlet differ from the prism section")
    a, b = edges.T
    wall_keys = np.stack((columns[:-1, a], columns[:-1, b],
                          columns[1:, b], columns[1:, a]), axis=-1).reshape(-1, 4)
    if not np.array_equal(ordered_rows(row_keys(quads.ids)), ordered_rows(row_keys(wall_keys))):
        raise ValueError("The observed walls differ from the straight prism sweep")
    lo, hi = bounds
    tol = 1e-10*(hi[0]-lo[0])
    tags = {}
    for axis, dimension in enumerate("xyz"):
        face_group = caps if axis == 2 else quads
        xyz = points[face_group.ids]
        normals = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
        for side, coordinate, direction in (("-", lo[axis], -1), ("+", hi[axis], 1)):
            on_plane = (abs(xyz[:, :, axis]-coordinate) <= tol).all(axis=1)
            found = np.unique(face_group.tags[on_plane])
            if len(found) != 1 or np.any(direction*normals[on_plane, axis] <= 0):
                raise ValueError(f"Expected consistently oriented CAD face {dimension}{side}")
            tags[dimension+side] = int(found[0])
    if len(set(tags.values())) != 6:
        raise ValueError("The CAD face tags do not distinguish all six duct sides")
    return tags


class BoundaryCells:
    def __init__(self, cells):
        self.ids = np.array([cell.ids for cell in cells], dtype=np.int64)
        self.tags = np.array([cell.tag for cell in cells], dtype=np.int64)


def extract_section(archive, metres_per_unit):
    points = archive.points*metres_per_unit
    volumes = [cell for cell in archive.cells if cell.dim == 3]
    surface = [cell for cell in archive.cells if cell.dim == 2]
    if (not volumes or not surface or any(cell.kind != 7 for cell in volumes)
            or any(cell.kind not in (2, 4) for cell in surface)):
        raise ValueError("Expected Prism6 volume cells with Tri3 caps and Quad4 walls")
    body = {cell.tag for cell in volumes}
    if len(body) != 1:
        raise ValueError("Expected a single fluid region")
    prisms = np.array([cell.ids for cell in volumes])
    caps = BoundaryCells([cell for cell in surface if cell.kind == 2])
    quads = BoundaryCells([cell for cell in surface if cell.kind == 4])
    if not len(caps.ids) or not len(quads.ids) or len(np.unique(prisms)) != len(points):
        raise ValueError("The volume, its exterior or its active node coverage is incomplete")

    bounds = np.stack((points.min(axis=0), points.max(axis=0)))
    lo, hi = bounds
    side, height, length = hi-lo
    if min(side, height, length) <= 0 or not np.isclose(side, height, atol=0, rtol=1e-10):
        raise ValueError("Expected a square section in xy and straight axial layers in z")
    z, axial, counts = np.unique(points[:, 2], return_inverse=True, return_counts=True)
    if len(z) < 2 or np.any(counts != counts[0]):
        raise ValueError("Each axial layer must contain the same section nodes")
    columns = np.lexsort((points[:, 1], points[:, 0], points[:, 2])).reshape(len(z), -1)
    section_points = points[columns[0], :2]
    projected = points[columns, :2]
    if (not np.allclose(projected, section_points, atol=side*1e-12, rtol=0)
            or len(np.unique(section_points, axis=0)) != len(section_points)):
        raise ValueError("The volume does not repeat individual section coordinates")
    index = np.empty(len(points), dtype=int)
    index[columns.ravel()] = np.tile(np.arange(len(section_points)), len(z))
    lower, upper = prisms[:, :3], prisms[:, 3:]
    prism_levels = axial[lower[:, 0]]
    if (not np.all(axial[lower] == prism_levels[:, None])
            or not np.all(axial[upper] == prism_levels[:, None]+1)
            or not np.array_equal(index[lower], index[upper])):
        raise ValueError("Found a prism that does not repeat the same ordered triangular section")
    triangles = index[lower[prism_levels == 0]]
    triangle_keys = ordered_rows(row_keys(triangles))
    if not len(triangles) or len(np.unique(triangle_keys, axis=0)) != len(triangles):
        raise ValueError("Missing or duplicate section triangles")
    for k in range(1, len(z)-1):
        if not np.array_equal(ordered_rows(row_keys(index[lower[prism_levels == k]])),
                              triangle_keys):
            raise ValueError("The triangular connectivity changed along the extrusion")
    xy = (section_points-lo[:2])/side
    if np.any(signed_double_area(xy[index[lower]]) <= 0):
        raise ValueError("The archive contains a nonpositive prism")
    edges, owners = section_edges(xy, triangles)
    tags = verify_exterior(points, caps, quads, columns, triangles, edges, bounds)
    return points, prisms, section_points, index, triangles, edges, owners, z, tags, next(iter(body)), bounds


def reference_conductance():
    n = np.arange(1, 1024, 2, dtype=float)
    value = (1 - 192/np.pi**5*np.sum(np.tanh(n*np.pi/2)/n**5))/12
    # Independently integrate the double-sine expansion of the unit forcing.
    squared = n*n
    product = squared[:, None]*squared[None, :]
    result = 64/np.pi**6*np.sum(1/(product*(squared[:, None]+squared[None, :])))
    np.testing.assert_allclose(value, result, rtol=1e-8)
    return float(value), float(result)


def analytical_velocity(xy):
    """Dimensionless velocity on [0,1]^2 for unit pressure drive."""
    n = np.arange(1, 384, 2, dtype=float)
    values = np.empty(len(xy))
    for start in range(0, len(xy), 4096):
        end = min(start+4096, len(xy))
        x, y = np.clip(xy[start:end], 0., 1.).T[:, :, None]
        k = n*np.pi
        ratio = (np.exp(-x*k) + np.exp(-(1-x)*k))/(1+np.exp(-k))
        values[start:end] = (4/np.pi**3)*np.sum(np.sin(y*k)*(1-ratio)/n**3, axis=1)
    return values


def solve_unit_section(xy, triangles, edges):
    corners = xy[triangles]
    jacobian = (corners[:, 1:]-corners[:, :1]).transpose(0, 2, 1)
    areas = signed_double_area(corners)/2
    np.testing.assert_allclose(areas.sum(), 1., rtol=0, atol=1e-10)
    inv = np.linalg.inv(jacobian)
    gradients = np.concatenate((-inv.sum(axis=1, keepdims=True), inv), axis=1)
    local = (gradients @ gradients.transpose(0, 2, 1))*areas[:, None, None]
    rows, cols = np.broadcast_arrays(triangles[:, :, None], triangles[:, None, :])
    matrix = coo_matrix((local.ravel(), (rows.ravel(), cols.ravel())),
                        shape=(len(xy), len(xy))).tocsc()
    load = np.bincount(triangles.ravel(), weights=np.repeat(areas/3, 3), minlength=len(xy))
    free = np.ones(len(xy), dtype=bool)
    free[edges] = False
    if not free.any():
        raise ValueError("The section has no interior degrees of freedom")
    velocity = np.zeros(len(xy))
    velocity[free] = spsolve(matrix[free][:, free], load[free])
    if not np.isfinite(velocity).all():
        raise ValueError("The section Poisson solve did not produce a finite solution")
    residual = matrix @ velocity - load
    error = np.linalg.norm(residual[free])/np.linalg.norm(load[free])
    conductance = float(load @ velocity)
    energy = float(velocity @ (matrix @ velocity))
    wall_reaction = float(-residual[~free].sum())
    if conductance<=0 or error>1e-9:
        raise ValueError("The section solve failed the flux or residual check")
    np.testing.assert_allclose(energy, conductance, rtol=1e-9)
    np.testing.assert_allclose(wall_reaction, 1., atol=1e-9, rtol=0)
    return velocity, conductance, areas, gradients, int(free.sum()), float(error), energy, wall_reaction


def velocity_error(xy, triangles, areas, velocity, reference_scale):
    q, w = np.polynomial.legendre.leggauss(6)
    q, w = (q+1)/2, w/2
    corners, values = xy[triangles], velocity[triangles]
    numerator = denominator = 0.
    for u, wu in zip(q, w):
        for v, wv in zip(q, w):
            shape = np.array([(1-u)*(1-v), u, (1-u)*v])
            reference = reference_scale*analytical_velocity(shape @ corners)
            finite_element = values @ shape
            weight = 2*(1-u)*wu*wv * areas
            numerator += np.dot(weight, (finite_element-reference)**2)
            denominator += np.dot(weight, reference**2)
    return float(np.sqrt(numerator/denominator))


def velocity_profile(xy, triangles, velocity, reference_scale):
    abscissas = np.linspace(0., 1., 101)
    positions = np.column_stack((abscissas, np.full(len(abscissas), .5)))
    corners = xy[triangles]
    lower, upper = corners.min(axis=1), corners.max(axis=1)
    profile = np.empty(len(positions))
    for i, point in enumerate(positions):
        ids = np.flatnonzero(np.all((lower-1e-12 <= point) & (point <= upper+1e-12), axis=1))
        xyz = corners[ids]
        jacobian = (xyz[:, 1:]-xyz[:, :1]).transpose(0, 2, 1)
        solved = np.linalg.solve(jacobian, (point-xyz[:, 0])[:, :, None]).squeeze(-1)
        barycentric = np.column_stack((1-solved.sum(axis=1), solved))
        found = np.min(barycentric, axis=1) >= -1e-11
        if not found.any():
            raise ValueError("The computed section does not cover a profile sample")
        values = np.sum(barycentric[found]*velocity[triangles[ids[found]]], axis=1)
        np.testing.assert_allclose(values, values[0], atol=1e-11, rtol=0)
        profile[i] = values[0]
    return abscissas, profile, reference_scale*analytical_velocity(positions)


def solve(label, path, output, metres_per_unit, mu, density, bulk_speed):
    start = time.perf_counter()
    archive_hash = sha256(path)
    archive = read_archive(path, max_bytes=256 << 20, max_nodes=300_000, max_cells=800_000)
    (points, prisms, section, section_index, triangles, edges, owners, z, tags, body,
     bounds) = extract_section(archive, metres_per_unit)
    lo, hi = bounds
    side, _, length = hi-lo
    xy = (section-lo[:2])/side
    prepared = time.perf_counter()
    velocity, conductance, areas, gradients, free, error, energy, reaction = solve_unit_section(
        xy, triangles, edges)
    solved = time.perf_counter()
    ref, independent_ref = reference_conductance()
    target_flow = bulk_speed*side*side
    slope = mu*bulk_speed/(conductance*side*side)
    exact_slope = mu*bulk_speed/(ref*side*side)
    velocity *= bulk_speed/conductance
    exact_scale = bulk_speed/ref
    exact_velocity = exact_scale*analytical_velocity(xy)
    achieved_flow = float(np.dot(areas, velocity[triangles].mean(axis=1))*side**2)
    np.testing.assert_allclose(achieved_flow, target_flow, rtol=1e-12)
    profile_x, profile, profile_exact = velocity_profile(xy, triangles, velocity, exact_scale)
    relative_error = velocity_error(xy, triangles, areas, velocity, exact_scale)

    delta = xy[edges[:, 1]]-xy[edges[:, 0]]
    edge_lengths = np.linalg.norm(delta, axis=1)
    outward = np.column_stack((delta[:, 1], -delta[:, 0]))/edge_lengths[:, None]
    cell_gradient = np.einsum("ei,eid->ed", velocity[triangles], gradients)
    wall_shear = -mu/side*np.sum(outward*cell_gradient[owners], axis=1)
    mean_wall_shear = float(np.average(wall_shear, weights=edge_lengths))
    reynolds = density*bulk_speed*side/mu
    reference = dict(conductance=ref, conductance_double_sine=independent_ref,
                     pressure_gradient_pa_m=float(exact_slope), pressure_drop_pa=float(exact_slope*length),
                     centerline_velocity_m_s=float(profile_exact[len(profile_exact)//2]),
                     darcy_friction_factor=float(2*exact_slope*side/(density*bulk_speed**2)),
                     mean_wall_shear_pa=float(exact_slope*side/4))
    # Measure the first normal node interval at each wall from the exported section.
    wall_spacing = []
    for axis in (0, 1):
        interior = xy[(xy[:, axis] > 1e-10) & (xy[:, axis] < 1-1e-10), axis]
        wall_spacing.extend([float(side*interior.min()), float(side*(1-interior.max()))])
    signature = output/f"{label}.npz"
    if sha256(path) != archive_hash:
        raise ValueError("The source Prism6 archive changed during the solve")
    np.savez_compressed(signature, points=points, prisms=prisms, velocity=velocity[section_index],
                       exact_velocity=exact_velocity[section_index], pressure=slope*(hi[2]-points[:, 2]),
                       section_points=section, section_triangles=triangles, section_velocity=velocity,
                       section_exact_velocity=exact_velocity, profile_coordinate_m=lo[0]+profile_x*side,
                       profile_velocity=profile, profile_exact_velocity=profile_exact,
                       wall_edges=edges, wall_shear_pa=wall_shear)
    result = dict(label=label, archive=str(path), archive_sha256=archive_hash,
                  solution=signature.name, solution_sha256=sha256(signature),
                  prism6=len(prisms), nodes=len(points), section_nodes=len(xy), section_tri3=len(triangles),
                  section_free_dofs=free, axial_layers=len(z)-1, wall_face_tags=tags, body_tag=body,
                  first_wall_spacing_m=min(wall_spacing), wall_spacings_m=wall_spacing,
                  axial_pitch_min_m=float(np.diff(z).min()), axial_pitch_max_m=float(np.diff(z).max()),
                  conductance=conductance, pressure_gradient_pa_m=float(slope),
                  pressure_drop_pa=float(slope*length),
                  pressure_drop_error_percent=100*float(slope/exact_slope-1),
                  relative_velocity_l2=relative_error, flow_m3_s=achieved_flow,
                  centerline_velocity_m_s=float(profile[len(profile)//2]),
                  darcy_friction_factor=float(2*slope*side/(density*bulk_speed**2)),
                  direct_wall_shear_pa=mean_wall_shear,
                  mean_wall_shear_at_fem_pressure_gradient_pa=float(slope*side/4),
                  force_balanced_wall_shear_pa=float(reaction*slope*side/4),
                  free_residual_relative=error, dimensionless_energy=energy, dimensionless_wall_reaction=reaction,
                  preparation_seconds=prepared-start, section_solve_seconds=solved-prepared)
    geometry = dict(side_m=float(side), length_m=float(length), bounds_m=bounds.tolist())
    fluid = dict(dynamic_viscosity_pa_s=mu, density_kg_m3=density, bulk_speed_m_s=bulk_speed,
                 reynolds=float(reynolds), target_flow_m3_s=float(target_flow))
    return geometry, fluid, reference, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", nargs=2, metavar=("LABEL", "ARCHIVE"), action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metres-per-unit", type=float, default=.001)
    parser.add_argument("--viscosity", type=float, default=.001)
    parser.add_argument("--density", type=float, default=1000.)
    parser.add_argument("--bulk-speed", type=float, default=.01)
    args = parser.parse_args()
    physical = (args.metres_per_unit, args.viscosity, args.density, args.bulk_speed)
    if not all(np.isfinite(physical)) or min(physical) <= 0:
        parser.error("Use finite positive units, viscosity, density and speed")
    labels = [label for label, _ in args.mesh]
    if len(set(labels)) != len(labels) or not all(re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_-]*", label)
                                                 for label in labels):
        parser.error("Choose distinct path-safe mesh labels")
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for label, archive in args.mesh:
        geometry, fluid, reference, row = solve(label, Path(archive), args.output, *physical)
        if rows:
            if (geometry != common_geometry or fluid != common_fluid or reference != common_reference):
                raise ValueError("Convergence levels must use identical duct bounds and physical inputs")
        else:
            common_geometry, common_fluid, common_reference = geometry, fluid, reference
        rows.append(row)
        print(f"{label}  {row['prism6']:,} Prism6  pressure drop {row['pressure_drop_pa']:.7f} Pa"
              f"  velocity L2 {100*row['relative_velocity_l2']:.3f}%")
    result = dict(schema="paacs-square-duct-flow-v1", geometry=common_geometry, fluid=common_fluid,
                  reference=common_reference, metres_per_archive_unit=args.metres_per_unit,
                  model="Fully developed steady axial flow, exact scalar reduction of incompressible Newtonian equations",
                  meshes=rows)
    (args.output/"summary.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n",
                                             encoding="utf-8")


if __name__ == "__main__":
    main()
