# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Solve periodic scalar viscous flow on an exported PAACS Prism6/Tet4 duct."""

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
from examples.flow_straight_duct.solve import analytical_velocity, reference_conductance
from scripts.mesh_visual_data import FACES, TET_EDGES, read_archive, sha256


def check_boundary(points, tets, prisms, surface):
    """Verify mixed face incidence and original CAD planes before applying BCs."""
    remaining = {}
    paired = set()
    for kind, elements in ((5, tets), (7, prisms)):
        for element in elements.tolist():
            for face in FACES[kind]:
                row = tuple(element[i] for i in face)
                key = tuple(sorted(row))
                if key in paired:
                    raise ValueError("A mixed face has more than two owners")
                previous = remaining.pop(key, None)
                if previous is None:
                    remaining[key] = row
                    continue
                index = previous.index(row[0])
                if row != tuple(previous[(index-i) % len(row)] for i in range(len(row))):
                    raise ValueError("An internal hybrid face has inconsistent winding")
                paired.add(key)
    faces = {name: [] for name in ("x-", "x+", "y-", "y+", "z-", "z+")}
    upper = points.max(axis=0)
    for cell in surface:
        row = cell.ids
        expected = remaining.pop(tuple(sorted(row)), None)
        if expected is None:
            raise ValueError("The archive boundary does not match the mixed cells")
        index = expected.index(row[0])
        if row != expected[index:] + expected[:index]:
            raise ValueError("A boundary facet has inconsistent winding")
        xyz = points[list(row)]
        sides = [axis+sign for d, axis in enumerate("xyz") for k, sign in enumerate("-+")
                 if np.allclose(xyz[:, d], upper[d] if k else 0., atol=1e-10, rtol=0)]
        if len(sides) != 1:
            raise ValueError("A boundary facet does not belong to a unique duct plane")
        faces[sides[0]].append(cell)
    if remaining:
        raise ValueError("Missing exported faces for the hybrid volume")
    tags = {side: {cell.tag for cell in cells} for side, cells in faces.items()}
    if any(len(value) != 1 for value in tags.values()) or len(set.union(*tags.values())) != 6:
        raise ValueError("Expected six distinct CAD side tags")
    return faces, {side: next(iter(value)) for side, value in tags.items()}


def periodic_dofs(points, faces):
    masters = np.arange(len(points))
    ends = []
    for side in ("z-", "z+"):
        nodes = np.unique([i for cell in faces[side] for i in cell.ids])
        x, y = points[nodes, :2].T
        ends.append(nodes[np.lexsort((y, x))])
    inlet, outlet = ends
    if not np.array_equal(points[inlet, :2], points[outlet, :2]):
        raise ValueError("Periodic caps must have identical transverse coordinates")
    if len(np.unique(points[inlet, :2], axis=0)) != len(inlet):
        raise ValueError("Periodic cap coordinates must identify unique nodes")
    masters[outlet] = inlet
    cap_keys = [sorted(tuple(sorted(masters[list(cell.ids)])) for cell in faces[side])
                for side in ("z-", "z+")]
    if cap_keys[0] != cap_keys[1]:
        raise ValueError("Periodic caps must have matching Tri3 and Quad4 connectivity")
    _, dofs = np.unique(masters, return_inverse=True)
    wall = np.unique([i for side in ("x-", "x+", "y-", "y+")
                      for cell in faces[side] for i in cell.ids])
    fixed = np.unique(dofs[wall])
    return dofs, fixed, wall, len(inlet)


def wedge_shapes(r, s, t):
    barycentric = np.array([1-r-s, r, s])
    triangle_gradient = np.array([[-1., -1.], [1., 0.], [0., 1.]])
    values = np.concatenate(((1-t)*barycentric, t*barycentric))
    gradient = np.empty((6, 3))
    gradient[:3, :2] = (1-t)*triangle_gradient
    gradient[3:, :2] = t*triangle_gradient
    gradient[:3, 2] = -barycentric
    gradient[3:, 2] = barycentric
    return values, gradient


def tet_elements(corners):
    jacobians = (corners[:, 1:] - corners[:, :1]).transpose(0, 2, 1)
    det = np.linalg.det(jacobians)
    if np.any(det <= 0):
        raise ValueError("The volume has a nonpositive Tet4 Jacobian")
    inv = np.linalg.inv(jacobians)
    gradient = np.concatenate((-inv.sum(axis=1, keepdims=True), inv), axis=1)
    volumes = det/6
    stiffness = (gradient @ gradient.transpose(0, 2, 1))*volumes[:, None, None]
    return stiffness, np.repeat(volumes[:, None]/4, 4, axis=1)


def wedge_quadrature(corners, order):
    q, w = np.polynomial.legendre.leggauss(order)
    q, w = (q+1)/2, w/2
    position = (corners-corners[:, :1]).transpose(0, 2, 1)
    for a, wa in zip(q, w):
        for b, wb in zip(q, w):
            for t, wt in zip(q, w):
                shape, parent_gradient = wedge_shapes(a, (1-a)*b, t)
                jacobian = position @ parent_gradient
                det = np.linalg.det(jacobian)
                if np.any(det <= 0):
                    raise ValueError("The volume has a nonpositive Prism6 Jacobian at quadrature")
                gradient = parent_gradient @ np.linalg.inv(jacobian)
                volume = det*(1-a)*wa*wb*wt
                yield shape, gradient, volume


def wedge_elements(corners, order):
    stiffness = np.zeros((len(corners), 6, 6))
    load = np.zeros((len(corners), 6))
    for shape, gradient, volume in wedge_quadrature(corners, order):
        stiffness += (gradient @ gradient.transpose(0, 2, 1))*volume[:, None, None]
        load += volume[:, None]*shape
    return stiffness, load


def assemble(points, tets, prisms, dofs, order):
    ndofs = int(dofs.max())+1
    load = np.zeros(ndofs)
    rows, cols, values = [], [], []
    for cells, local in ((tets, tet_elements(points[tets])),
                         (prisms, wedge_elements(points[prisms], order))):
        stiffness, force = local
        ids = dofs[cells]
        r, c = np.broadcast_arrays(ids[:, :, None], ids[:, None, :])
        rows.append(r.ravel())
        cols.append(c.ravel())
        values.append(stiffness.ravel())
        load += np.bincount(ids.ravel(), weights=force.ravel(), minlength=ndofs)
    matrix = coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))),
                        shape=(ndofs, ndofs)).tocsc()
    return matrix, load


def volume_diagnostics(points, tets, prisms, velocity, exact_scale, order):
    numerator = denominator = 0.

    def measure(corners, values, weights, shape):
        nonlocal numerator, denominator
        exact = exact_scale*analytical_velocity(shape @ corners[..., :2])
        numerical = values @ shape
        numerator += np.dot(weights, (exact-numerical)**2)
        denominator += np.dot(weights, exact**2)

    tet_corners, wedge_corners = points[tets], points[prisms]
    tet_values, wedge_values = velocity[tets], velocity[prisms]
    jacobian = (tet_corners[:, 1:] - tet_corners[:, :1]).transpose(0, 2, 1)
    det = np.linalg.det(jacobian)
    gradient_z = np.sum((tet_values[:, 1:]-tet_values[:, :1])*np.linalg.inv(jacobian)[:, :, 2], axis=1)
    axial_energy = np.dot(det/6, gradient_z**2)
    q, w = np.polynomial.legendre.leggauss(order)
    q, w = (q+1)/2, w/2
    for a, wa in zip(q, w):
        for b, wb in zip(q, w):
            for c, wc in zip(q, w):
                r, s, t = a, (1-a)*b, (1-a)*(1-b)*c
                measure(tet_corners, tet_values, det*(1-a)**2*(1-b)*wa*wb*wc,
                        np.array([1-r-s-t, r, s, t]))
    for shape, gradient, volume in wedge_quadrature(wedge_corners, order):
        gradient_z = np.sum((wedge_values-wedge_values[:, :1])*gradient[:, :, 2], axis=1)
        axial_energy += np.dot(volume, gradient_z**2)
        measure(wedge_corners, wedge_values, volume, shape)
    return float(np.sqrt(numerator/denominator)), float(axial_energy)


def intersect_edges(vertices, values, edges, z):
    intersections = []
    for a, b in edges:
        za, zb = vertices[a, 2], vertices[b, 2]
        if min(za, zb) < z < max(za, zb):
            fraction = (z-za)/(zb-za)
            intersections.append((vertices[a, :2]*(1-fraction) + vertices[b, :2]*fraction,
                                  values[a]*(1-fraction) + values[b]*fraction))
    return intersections


def cross_2d(a, b):
    return a[0]*b[1] - a[1]*b[0]


def integrate_cut(points, tets, prisms, velocity, z):
    area, flow = 0., 0.
    q, w = np.polynomial.legendre.leggauss(2)
    q, w = (q+1)/2, w/2
    for cells, is_wedge in ((tets, False), (prisms, True)):
        levels = points[cells, 2]
        intersects = (levels.min(axis=1) < z) & (z < levels.max(axis=1))
        for cell in cells[intersects]:
            vertices, values = points[cell], velocity[cell]
            if is_wedge:
                lower = intersect_edges(vertices[:3], values[:3], ((0, 1), (1, 2), (2, 0)), z)
                upper = intersect_edges(vertices[3:], values[3:], ((0, 1), (1, 2), (2, 0)), z)
                if len(lower) != 2 or len(upper) != 2:
                    raise ValueError("Expected two intersections on each paired prism triangle")
                a, b, c, d = [p for p, _ in lower+upper]
                va, vb, vc, vd = [u for _, u in lower+upper]
                # At fixed z the triangular basis restricts to a line and
                # its product with the layer-normal basis is bilinear.
                for s, ws in zip(q, w):
                    for t, wt in zip(q, w):
                        ds = (1-t)*(b-a) + t*(d-c)
                        dt = (1-s)*(c-a) + s*(d-b)
                        weight = abs(cross_2d(ds, dt))*ws*wt
                        u = (1-t)*((1-s)*va+s*vb) + t*((1-s)*vc+s*vd)
                        area += weight
                        flow += u*weight
            else:
                pairs = intersect_edges(vertices, values, TET_EDGES, z)
                if len(pairs) not in (3, 4):
                    raise ValueError("Expected a triangular or quadrilateral Tet4 section")
                center = np.mean([p for p, _ in pairs], axis=0)
                pairs.sort(key=lambda pair: np.arctan2(*(pair[0]-center)[::-1]))
                a, va = pairs[0]
                for (b, vb), (c, vc) in zip(pairs[1:], pairs[2:]):
                    triangle_area = abs(cross_2d(b-a, c-a))/2
                    area += triangle_area
                    flow += triangle_area*(va+vb+vc)/3
    return area, flow


def probe_profile(points, tets, prisms, velocity, z):
    abscissa = np.linspace(0., 1., 101)
    profile = []
    corner_groups = [(cells, points[cells]) for cells in (tets, prisms)]
    groups = [(cells, corners, corners.min(axis=1), corners.max(axis=1))
              for cells, corners in corner_groups]
    for x in abscissa:
        probe = np.array([x, .5, z])
        samples = []
        for k, (cells, all_corners, lower, upper) in enumerate(groups):
            possible = ((lower-1e-10 <= probe) & (probe <= upper+1e-10)).all(axis=1)
            for cell, corners in zip(cells[possible], all_corners[possible]):
                if k == 0:
                    coordinates = np.linalg.solve((corners[1:]-corners[0]).T, probe-corners[0])
                    shape = np.r_[1-coordinates.sum(), coordinates]
                else:
                    coordinates = np.array([1/3, 1/3, .5])
                    for _ in range(12):
                        shape, gradient = wedge_shapes(*coordinates)
                        difference = probe - shape @ corners
                        if np.linalg.norm(difference) <= 1e-12:
                            break
                        coordinates += np.linalg.solve(corners.T @ gradient, difference)
                    else:
                        continue
                if shape.min() >= -1e-9:
                    samples.append(shape @ velocity[cell])
        if not samples:
            raise ValueError(f"The hybrid mesh does not cover centerline sample x={x:g}")
        np.testing.assert_allclose(samples, samples[0], atol=1e-10, rtol=0)
        profile.append(samples[0])
    return abscissa, np.array(profile)


def solve(label, path, output, metres_per_unit, mu, density, bulk, order):
    began = time.perf_counter()
    source = sha256(path)
    mesh = read_archive(path, max_bytes=192 << 20, max_nodes=250_000, max_cells=450_000)
    physical = mesh.points*metres_per_unit
    lo, hi = physical.min(axis=0), physical.max(axis=0)
    side, height, length = hi-lo
    if min(side, height, length) <= 0 or not np.isclose(side, height, atol=0, rtol=1e-10):
        raise ValueError("Expected an xy square duct with axial coordinate z")
    points = (physical-lo)/side
    volumes = [cell for cell in mesh.cells if cell.dim == 3]
    surface = [cell for cell in mesh.cells if cell.dim == 2]
    if (any(cell.kind not in (5, 7) for cell in volumes) or not surface
            or any(cell.kind not in (2, 4) for cell in surface) or len({c.tag for c in volumes}) != 1):
        raise ValueError("Expected one region of classified Prism6 and Tet4 volume")
    tets = np.array([cell.ids for cell in volumes if cell.kind == 5])
    prisms = np.array([cell.ids for cell in volumes if cell.kind == 7])
    if not len(tets) or not len(prisms) or len(np.unique(np.r_[tets.ravel(), prisms.ravel()])) != len(points):
        raise ValueError("Expected a populated mixed volume with no unconnected nodes")
    if not np.array_equal(points[prisms[:, :3], 2], points[prisms[:, 3:], 2]):
        raise ValueError("Prism6 layer-normal edges must preserve axial coordinate for section integration")
    faces, tags = check_boundary(points, tets, prisms, surface)
    dofs, fixed, wall, cap_nodes = periodic_dofs(points, faces)
    prepared = time.perf_counter()
    matrix, load = assemble(points, tets, prisms, dofs, order)
    length_normalized = length/side
    np.testing.assert_allclose(load.sum(), length_normalized, atol=1e-10, rtol=0)
    free = np.ones(len(load), dtype=bool)
    free[fixed] = False
    if not free.any():
        raise ValueError("The hybrid system has no interior degrees of freedom")
    unit = np.zeros(len(load))
    unit[free] = spsolve(matrix[free][:, free], load[free])
    solved = time.perf_counter()
    if not np.isfinite(unit).all():
        raise ValueError("The periodic linear solve did not produce finite velocities")
    work = float(load @ unit)
    residual = matrix @ unit - load
    energy = float(unit @ (matrix @ unit))
    residual_relative = np.linalg.norm(residual[free])/np.linalg.norm(load[free])
    if work<=0 or energy<=0 or residual_relative>1e-9:
        raise ValueError("The periodic solve failed the flow or residual check")
    np.testing.assert_allclose(work, energy, atol=0, rtol=1e-9)
    np.testing.assert_allclose(-residual[fixed].sum(), length_normalized, rtol=1e-9)

    conductance = work/length_normalized
    ref, independent = reference_conductance()
    pressure_gradient = mu*bulk/(conductance*side**2)
    reference_gradient = mu*bulk/(ref*side**2)
    velocity = unit[dofs]*(bulk/conductance)
    exact_scale = bulk/ref
    analytical = analytical_velocity(points[:, :2])*exact_scale
    volume_error, axial_energy = volume_diagnostics(points, tets, prisms, velocity, exact_scale, order)
    target_flow = bulk*side**2
    net_flow = work*(bulk/conductance)*side**3/length
    pressure_drop = pressure_gradient*length
    dissipation = mu*side*(bulk/conductance)**2*energy
    np.testing.assert_allclose(dissipation, pressure_drop*net_flow, rtol=1e-9)
    if np.any(velocity[wall] != 0.):
        raise ValueError("Nonzero wall velocity")
    sections = []
    for fraction in (.137, .319, .571, .793, .917):
        z = fraction*length_normalized
        slab_area, flux = integrate_cut(points, tets, prisms, velocity, z)
        np.testing.assert_allclose(slab_area, 1., rtol=1e-8)
        flux *= side**2
        sections.append(dict(z_m=float(lo[2]+z*side), area_m2=float(slab_area*side**2),
                             flow_m3_s=float(flux), flux_error_percent=float(100*(flux/target_flow-1))))
    x, profile = probe_profile(points, tets, prisms, velocity, length_normalized/2)
    profile_exact = exact_scale*analytical_velocity(np.column_stack((x, np.full(len(x), .5))))

    xyz = points[prisms]
    first_spacing = []
    for axis in (0, 1):
        for coordinate in (0, 1):
            for offset, opposite in ((0, 3), (3, 0)):
                on_wall = np.isclose(xyz[:, offset:offset+3, axis], coordinate, atol=1e-10, rtol=0).all(axis=1)
                if on_wall.any():
                    first_spacing.extend(abs(xyz[on_wall, opposite:opposite+3, axis]-coordinate).ravel())
    if not first_spacing or min(first_spacing) <= 0:
        raise ValueError("Cannot measure the first prism interval")
    if sha256(path) != source:
        raise ValueError("The hybrid input changed during the solve")
    signature = output/f"{label}.npz"
    np.savez_compressed(signature, points=physical, prisms=prisms, tetrahedra=tets, velocity=velocity,
                       exact_velocity=analytical, pressure=pressure_gradient*(hi[2]-physical[:, 2]),
                       profile_coordinate_m=lo[0]+x*side, profile_velocity=profile,
                       profile_exact_velocity=profile_exact,
                       boundary_triangles=np.array([cell.ids for cell in surface if cell.kind == 2]),
                       boundary_quads=np.array([cell.ids for cell in surface if cell.kind == 4]))
    result = dict(label=label, archive=str(path), archive_sha256=source,
                  solution=signature.name, solution_sha256=sha256(signature),
                  prism6=len(prisms), tet4=len(tets), nodes=len(points), free_dofs=int(free.sum()),
                  periodic_node_pairs=cap_nodes, face_tags=tags,
                  first_wall_spacing_m=float(side*min(first_spacing)), quadrature_order=order,
                  conductance=conductance, pressure_gradient_pa_m=float(pressure_gradient),
                  pressure_drop_pa=float(pressure_drop),
                  pressure_drop_error_percent=float(100*(pressure_gradient/reference_gradient-1)),
                  centerline_velocity_m_s=float(profile[len(profile)//2]), flow_m3_s=float(net_flow),
                  section_fluxes=sections,
                  max_section_flux_error_percent=max(abs(row["flux_error_percent"]) for row in sections),
                  relative_axial_gradient_energy=float(np.sqrt(axial_energy/energy)*conductance/bulk),
                  relative_volume_velocity_l2=volume_error,
                  relative_profile_velocity_l2=float(np.linalg.norm(profile-profile_exact)/np.linalg.norm(profile_exact)),
                  free_residual_relative=float(residual_relative),
                  dissipation_w=float(dissipation), pressure_work_w=float(pressure_drop*net_flow),
                  wall_reaction=float(-residual[fixed].sum()), wall_velocity_max_m_s=float(np.max(abs(velocity[wall]))),
                  preparation_seconds=prepared-began, solve_seconds=solved-prepared)
    geometry = dict(side_m=float(side), length_m=float(length), bounds_m=[lo.tolist(), hi.tolist()])
    fluid = dict(dynamic_viscosity_pa_s=mu, density_kg_m3=density, bulk_speed_m_s=bulk,
                 reynolds=float(density*bulk*side/mu), target_flow_m3_s=float(target_flow))
    reference = dict(conductance=ref, conductance_double_sine=independent,
                     pressure_gradient_pa_m=float(reference_gradient), pressure_drop_pa=float(reference_gradient*length),
                     centerline_velocity_m_s=float(profile_exact[len(profile)//2]))
    return geometry, fluid, reference, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", nargs=2, metavar=("LABEL", "ARCHIVE"), action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metres-per-unit", type=float, default=.001)
    parser.add_argument("--viscosity", type=float, default=.001)
    parser.add_argument("--density", type=float, default=1000.)
    parser.add_argument("--bulk-speed", type=float, default=.01)
    parser.add_argument("--quadrature-order", type=int, choices=range(3, 7), default=4)
    args = parser.parse_args()
    physical = (args.metres_per_unit, args.viscosity, args.density, args.bulk_speed)
    if not all(np.isfinite(physical)) or min(physical) <= 0:
        parser.error("Use finite positive units, viscosity, density and bulk speed")
    labels = [label for label, _ in args.mesh]
    if len(set(labels)) != len(labels) or not all(re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_-]*", label)
                                                 for label in labels):
        parser.error("Choose distinct path-safe labels")
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for label, archive in args.mesh:
        geometry, fluid, reference, row = solve(label, Path(archive), args.output, *physical, args.quadrature_order)
        if rows and (geometry, fluid, reference) != common:
            raise ValueError("Convergence levels must use identical geometry, loading and fluid parameters")
        common = geometry, fluid, reference
        rows.append(row)
        print(f"{label}  {row['prism6']:,} Prism6 / {row['tet4']:,} Tet4"
              f"  pressure drop {row['pressure_drop_pa']:.7f} Pa"
              f"  worst section flux error {row['max_section_flux_error_percent']:.3f}%")
    result = dict(schema="paacs-hybrid-duct-flow-v1", geometry=geometry, fluid=fluid, reference=reference,
                  metres_per_archive_unit=args.metres_per_unit,
                  model="3D scalar axial-velocity Poisson equation with periodic end velocity and no-slip walls",
                  meshes=rows)
    (args.output/"summary.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n",
                                             encoding="utf-8")


if __name__ == "__main__":
    main()
