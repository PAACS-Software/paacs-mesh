# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Solve 3D linear elasticity on PAACS Tet4 plates with a circular through-hole."""

import argparse
import json
from pathlib import Path
import re
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from examples.elasticity_cantilever.solve import apply_constraints, assemble
from examples.elasticity_plate_hole.solver_input import read_solver_input
from scripts.mesh_visual_data import TET_EDGES, TET_FACES, read_archive, sha256

STRESS_WEIGHTS = np.array([1, 1, 1, 2, 2, 2])


def constitutive_matrix(young, poisson):
    shear = young / (2 * (1 + poisson))
    lame = young * poisson / ((1 + poisson) * (1 - 2 * poisson))
    elasticity = np.diag([2 * shear] * 3 + [shear] * 3)
    elasticity[:3, :3] += lame
    return elasticity


def kirsch_stress(xy, radius, remote):
    """In-plane Kirsch stress for remote x-tension, stored as xx yy zz xy yz xz."""
    x, y = np.moveaxis(xy, -1, 0)
    r2 = x*x + y*y
    c2, s2 = (x*x - y*y)/r2, 2*x*y/r2
    a2 = radius*radius/r2
    rr = .5*remote * ((1 - a2) + (1 - 4*a2 + 3*a2*a2)*c2)
    tt = .5*remote * ((1 + a2) - (1 + 3*a2*a2)*c2)
    rt = -.5*remote * (1 + 2*a2 - 3*a2*a2)*s2
    out = np.zeros((*r2.shape, 6))
    out[..., 0] = .5*((rr + tt) + (rr - tt)*c2) - rt*s2
    out[..., 1] = .5*((rr + tt) - (rr - tt)*c2) + rt*s2
    out[..., 3] = .5*(rr - tt)*s2 + rt*c2
    return out


def make_boundary(points, tets, triangles, tags, bounds, radius):
    """Identify exterior CAD faces and the cell supporting each boundary triangle."""
    faces = tets[:, TET_FACES].reshape(-1, 3)
    keys = np.sort(faces, axis=1)
    _, first, inverse, counts = np.unique(keys, axis=0, return_index=True,
                                         return_inverse=True, return_counts=True)
    if np.any(counts > 2):
        raise ValueError("The volume has a nonmanifold face")
    # Canonical key order allows the archive's surface records in any order.
    boundary_faces = faces[first[counts == 1]]
    triangle_keys = np.sort(triangles, axis=1)
    order = np.lexsort(triangle_keys.T[::-1])
    if not np.array_equal(triangle_keys[order], np.sort(boundary_faces, axis=1)):
        raise ValueError("The classified surface does not match the Tet4 boundary")
    owners = np.empty(len(triangles), dtype=int)
    owners[order] = first[counts == 1] // 4
    faces_group = faces[counts[inverse] == 2]
    if len(faces_group):
        interior_index = inverse[counts[inverse] == 2]
        paired = faces_group[np.argsort(interior_index, kind="stable")].reshape(-1, 2, 3)
        reverse = paired[:, 1, ::-1]
        if not np.any(np.stack([np.all(paired[:, 0] == np.roll(reverse, k, axis=1), axis=1)
                                for k in (0, 1, 2)]), axis=0).all():
            raise ValueError("The volume has inconsistent interior face winding")
    same_winding = [np.all(triangles[order] == np.roll(boundary_faces, k, axis=1), axis=1)
                    for k in (0, 1, 2)]
    if not np.any(np.stack(same_winding), axis=0).all():
        raise ValueError("The archive surface winding differs from the exterior Tet4 faces")

    lo, hi = bounds
    span = hi - lo
    tol = 1e-8*np.max(span)
    vertices = points[triangles]
    area = np.linalg.norm(np.cross(vertices[:, 1]-vertices[:, 0],
                                  vertices[:, 2]-vertices[:, 0]), axis=1)/2
    if np.any(area <= 0):
        raise ValueError("A classified boundary triangle has zero area")
    assigned = np.zeros(len(triangles), dtype=bool)
    sides = {}
    side_tags = {}
    for d, name in enumerate("xyz"):
        for k, sign in enumerate(("-", "+")):
            selected = (abs(vertices[:, :, d] - bounds[k, d]) <= tol).all(axis=1)
            native_tags = np.unique(tags[selected])
            if len(native_tags) != 1 or assigned[selected].any():
                raise ValueError(f"Expected one CAD face for outer plane {name}{sign}")
            if d < 2:
                planar_area = np.prod(np.delete(span, d))
                np.testing.assert_allclose(area[selected].sum(), planar_area, rtol=1e-8)
            assigned |= selected
            sides[name+sign] = selected
            side_tags[name+sign] = int(native_tags[0])
    hole = ~assigned
    hole_tags = np.unique(tags[hole])
    radial = np.linalg.norm(vertices[hole, :, :2], axis=2)
    if (len(hole_tags) != 1 or len(set(side_tags.values()) | set(hole_tags)) != 7
            or not np.allclose(radial, radius, rtol=0, atol=tol)):
        raise ValueError("Expected a single CAD-tagged cylindrical hole with vertices on its rim")
    return sides, side_tags, hole, int(hole_tags[0]), owners, area, tol


def exterior_traction(points, triangles, area, sides, radius, remote):
    """Integrate sigma_Kirsch n on x/y truncations. All other faces are free."""
    # The Kirsch traction is rational. A high-order tensor rule also
    # balances the force and moment on coarse, nonmatching outer faces.
    abscissa, weights = np.polynomial.legendre.leggauss(24)
    abscissa = (abscissa + 1)/2
    weights /= 2
    load = np.zeros_like(points)
    for direction, component in enumerate(("x", "y")):
        for sign, face in ((-1, component+"-"), (1, component+"+")):
            selected = sides[face]
            ids = triangles[selected]
            xyz = points[ids]
            forces = np.zeros_like(xyz)
            for u, wu in zip(abscissa, weights):
                for v, wv in zip(abscissa, weights):
                    shape = np.array([(1-u)*(1-v), u, (1-u)*v])
                    stress = kirsch_stress(shape @ xyz[..., :2], radius, remote)
                    traction = stress[:, [0, 3] if direction == 0 else [3, 1]]*sign
                    d_weight = area[selected]*(2*wu*wv*(1-u))
                    forces[:, :, :2] += traction[:, None, :] * (d_weight[:, None]*shape)[:, :, None]
            np.add.at(load, ids.ravel(), forces.reshape(-1, 3))
    return load


def corner_constraints(points, bounds, tol):
    # Three non-collinear corners supply the six rigid-motion constraints
    # while leaving the traction boundary unrestricted in deformation.
    lo, hi = bounds
    corners = np.array([lo, [hi[0], lo[1], lo[2]], [lo[0], hi[1], lo[2]]])
    distance = np.linalg.norm(points[None, :, :] - corners[:, None, :], axis=2)
    nodes = np.argmin(distance, axis=1)
    if np.any(distance[np.arange(3), nodes] > tol):
        raise ValueError("The plate lacks three non-collinear corners for rigid-mode constraints")
    return np.array([3*nodes[0], 3*nodes[0]+1, 3*nodes[0]+2,
                     3*nodes[1]+1, 3*nodes[1]+2, 3*nodes[2]+2])


def stress_sampler(points, tets, stress, tol):
    corners = points[tets]
    lower, upper = corners.min(axis=1), corners.max(axis=1)
    jacobian = (corners[:, 1:] - corners[:, :1]).transpose(0, 2, 1)
    reciprocal = np.linalg.inv(jacobian)
    centroid_index = cKDTree(corners.mean(axis=1))

    def sample(point):
        candidates = np.all((lower-tol <= point) & (point <= upper+tol), axis=1)
        ids = np.flatnonzero(candidates)
        solved = (reciprocal[ids] @ (point-corners[ids, 0])[:, :, None]).squeeze(-1)
        barycentric = np.column_stack((1-solved.sum(axis=1), solved))
        owners = ids[np.min(barycentric, axis=1) >= -1e-8]
        if len(owners) == 0:
            raise ValueError(f"No material tetrahedron at stress probe {point.tolist()}")
        return stress[owners].mean(axis=0), len(owners)
    def sample_many(probes):
        traces = np.empty((len(probes), 6))
        for start in range(0, len(probes), 1024):
            end = min(start+1024, len(probes))
            chunk = probes[start:end]
            _, nearby = centroid_index.query(chunk, k=min(32, len(tets)))
            if nearby.ndim == 1:
                nearby = nearby[:, None]
            v = chunk[:, None] - corners[nearby, 0]
            solved = (reciprocal[nearby] @ v[..., None]).squeeze(-1)
            inside = (solved.min(axis=2) > 1e-8) & (solved.sum(axis=2) < 1-1e-8)
            found = inside.any(axis=1)
            local = np.arange(len(chunk))
            owner = nearby[local, inside.argmax(axis=1)]
            chunk_traces = traces[start:end]
            chunk_traces[found] = stress[owner[found]]
            # A nearest-centroid candidate is only an acceleration. The
            # full bounding-box lookup handles elongated tets and edge traces.
            for i in np.flatnonzero(~found):
                traces[start+i], _ = sample(chunk[i])
        return traces

    return sample, sample_many


def radial_probes(sample, radius, remote, tol):
    probes = []
    for offset in (1., 1.05, 1.10, 1.20, 1.40):
        samples = []
        for sign in (-1, 1):
            nominal = np.array([0., sign*radius*offset, 0.])
            # Move into material by roundoff scale to avoid a lookup just
            # inside the CAD void when a rim node matches the ray.
            point = nominal.copy()
            if offset == 1:
                point[1] += sign*tol*.1
            trace, owner_count = sample(point)
            samples.append(dict(position_mm=nominal.tolist(), cells=owner_count,
                                fem_sigma_tt_mpa=float(trace[0])))
        value = np.mean([sample["fem_sigma_tt_mpa"] for sample in samples])
        reference = kirsch_stress(nominal[:2], radius, remote)[0]
        probes.append(dict(radius_over_hole_radius=offset, samples=samples,
                           fem_sigma_tt_mpa=float(value), kirsch_sigma_tt_mpa=float(reference),
                           fem_kt=float(value/remote), kirsch_kt=float(reference/remote)))
    return probes


def probe_grid(radius, half_height, minimum, maximum, angles, wa, radial_order, depth_order):
    q, wr = np.polynomial.legendre.leggauss(radial_order)
    z, wz = np.polynomial.legendre.leggauss(depth_order)
    radii = radius*((maximum-minimum)*q + maximum+minimum)/2
    rr, theta, zz = np.meshgrid(radii, angles, half_height*z, indexing="ij")
    xyz = np.column_stack(((rr*np.cos(theta)).ravel(), (rr*np.sin(theta)).ravel(), zz.ravel()))
    polar_weights = (radii*wr)[:, None, None] * wa[None, :, None] * wz[None, None, :]
    return xyz, theta.ravel(), polar_weights.ravel()


def strip_probe(sample_many, radius, remote, half_height):
    # Use identical physical samples on every mesh. This is finite
    # quadrature for a piecewise-constant stress field.
    q, w = np.polynomial.legendre.leggauss(24)
    half_angle = np.pi/12
    minimum, maximum = 1.05, 1.40
    angles = q*half_angle + np.array([-np.pi/2, np.pi/2])[:, None]
    xyz, theta, weights = probe_grid(radius, half_height, minimum, maximum, angles.ravel(),
                                     np.tile(w, 2), 16, 10)
    fem = sample_many(xyz)
    x, y = np.cos(theta), np.sin(theta)
    hoop = y*y*fem[:, 0] + x*x*fem[:, 1] - 2*x*y*fem[:, 3]
    # Integrate the analytic hoop stress with the polar area weight r dr.
    radial_measure = .5*(maximum**2-minimum**2)
    mean_a2_over_r2 = np.log(maximum/minimum)/radial_measure
    mean_a4_over_r4 = .5*(minimum**-2-maximum**-2)/radial_measure
    mean_cos2theta = -np.sin(2*half_angle)/(2*half_angle)
    analytical = .5*remote*(1 + mean_a2_over_r2 - (1 + 3*mean_a4_over_r4)*mean_cos2theta)
    return dict(radius_range_over_hole_radius=[minimum, maximum], half_angle_degrees=15,
                quadrature_samples=len(xyz), fem_sigma_tt_mpa=float(np.average(hoop, weights=weights)),
                kirsch_sigma_tt_mpa=float(analytical))


def fixed_annulus(sample_many, radius, remote, half_height):
    angles = 2*np.pi/192*(np.arange(192)+.5)
    minimum, maximum = 1.01, 2.
    xyz, _, weights = probe_grid(radius, half_height, minimum, maximum, angles,
                                 np.ones_like(angles), 16, 8)
    fem = sample_many(xyz)
    exact = kirsch_stress(xyz[:, :2], radius, remote)
    error2 = np.sum((fem-exact)**2 * STRESS_WEIGHTS, axis=1)
    exact2 = np.sum(exact**2 * STRESS_WEIGHTS, axis=1)
    relative = np.sqrt(np.dot(weights, error2)/np.dot(weights, exact2))
    return dict(radius_range_over_hole_radius=[minimum, maximum],
                quadrature_samples=len(xyz), relative_kirsch_stress_l2=float(relative))


def stress_diagnostics(stress, centroids, volumes, radius, remote):
    radial = np.linalg.norm(centroids[:, :2], axis=1)
    if (radial == 0).any():
        raise ValueError("The archive contains a tetrahedron centered on the hole axis")
    in_reference = radial >= radius
    exact = kirsch_stress(centroids[:, :2], radius, remote)
    # Circumferential unit direction in the xy plane.
    x, y = (centroids[:, :2]/radial[:, None]).T
    hoop = y*y*stress[:, 0] + x*x*stress[:, 1] - 2*x*y*stress[:, 3]
    vm2 = .5*((stress[:, 0]-stress[:, 1])**2 + (stress[:, 1]-stress[:, 2])**2 +
              (stress[:, 2]-stress[:, 0])**2) + 3*np.sum(stress[:, 3:]**2, axis=1)
    error2 = np.sum((stress-exact)**2 * STRESS_WEIGHTS, axis=1)
    exact2 = np.sum(exact**2 * STRESS_WEIGHTS, axis=1)
    relative = np.sqrt(np.dot(volumes[in_reference], error2[in_reference]) /
                       np.dot(volumes[in_reference], exact2[in_reference]))
    near = in_reference & (radial <= 2*radius)
    near_relative = np.sqrt(np.dot(volumes[near], error2[near]) /
                            np.dot(volumes[near], exact2[near])) if near.any() else None
    return hoop, np.sqrt(vm2), exact, dict(
        relative_kirsch_stress_l2=float(relative),
        relative_near_hole_kirsch_stress_l2=float(near_relative) if near_relative is not None else None,
        centroid_volume_inside_ideal_hole_mm3=float(volumes[~in_reference].sum()))


def unique_lengths(points, cells, edges):
    pairs = np.unique(np.sort(cells[:, edges].reshape(-1, 2), axis=1), axis=0)
    return np.linalg.norm(points[pairs[:, 1]] - points[pairs[:, 0]], axis=1)


def solve(label, path, output, radius, young, poisson, remote):
    start = time.perf_counter()
    original_hash = sha256(path)
    mesh = (read_solver_input(path) if path.suffix.lower() == ".inp"
            else read_archive(path, max_bytes=128 << 20, max_nodes=200_000, max_cells=240_000))
    bounds = np.array(mesh.bounds)
    lo, hi = bounds
    if (not np.allclose(lo, -hi, rtol=0, atol=max(hi-lo)*1e-8) or hi[2] <= 0
            or min(hi[:2]) <= 2*radius):
        raise ValueError("Expected a box centered at the origin enclosing the circular hole and its probes")
    volume_cells = [cell for cell in mesh.cells if cell.dim == 3]
    surface_cells = [cell for cell in mesh.cells if cell.dim == 2]
    if (not volume_cells or not surface_cells or any(cell.kind != 5 for cell in volume_cells)
            or any(cell.kind != 2 for cell in surface_cells)):
        raise ValueError("Expected classified Tri3 boundary and Tet4 volume")
    tets = np.array([cell.ids for cell in volume_cells])
    triangles = np.array([cell.ids for cell in surface_cells])
    tags = np.array([cell.tag for cell in surface_cells])
    body_tags = {cell.tag for cell in volume_cells}
    if len(body_tags) != 1 or len(np.unique(tets)) != len(mesh.points):
        raise ValueError("Expected a single Tet4 region with all nodes incident to volume cells")
    points = mesh.points
    sides, side_tags, hole, hole_tag, owners, area, tol = make_boundary(
        points, tets, triangles, tags, bounds, radius)
    elasticity = constitutive_matrix(young, poisson)
    matrix, strain, dofs, volume = assemble(points, tets, elasticity)
    load = exterior_traction(points, triangles, area, sides, radius, remote)
    fixed = corner_constraints(points, bounds, tol)
    assembled = time.perf_counter()
    displacement, free = apply_constraints(matrix, load.ravel(), fixed, 0.)
    solved = time.perf_counter()
    residual = matrix @ displacement - load.ravel()
    reaction = np.zeros_like(residual)
    reaction[fixed] = residual[fixed]
    reaction = reaction.reshape(-1, 3)
    load_scale = np.linalg.norm(load, axis=1).sum()
    moment_scale = load_scale * np.linalg.norm(hi-lo)
    force_balance = (load+reaction).sum(axis=0)
    moment_balance = np.cross(points, load+reaction).sum(axis=0)
    free_residual = np.linalg.norm(residual[free])/np.linalg.norm(load)
    if (free_residual > 1e-8 or np.linalg.norm(force_balance) > load_scale*1e-8
            or np.linalg.norm(moment_balance) > moment_scale*1e-8):
        raise ValueError("Sparse elasticity equilibrium check failed")
    support_reaction = np.linalg.norm(reaction)/load_scale
    if support_reaction > 1e-8:
        raise ValueError(f"Rigid-body constraints carried relative force {support_reaction:.3g}")
    work = float(load.ravel() @ displacement)
    energy = float(.5*(displacement @ (matrix @ displacement)))
    np.testing.assert_allclose(work, 2*energy, rtol=1e-8)

    stress = (strain @ displacement[dofs, None]).squeeze(-1) @ elasticity
    centroids = points[tets].mean(axis=1)
    hoop, mises, exact, error = stress_diagnostics(stress, centroids, volume, radius, remote)
    sample, sample_many = stress_sampler(points, tets, stress, tol)
    probes = radial_probes(sample, radius, remote, tol)
    strip = strip_probe(sample_many, radius, remote, hi[2])
    annulus = fixed_annulus(sample_many, radius, remote, hi[2])
    edges = unique_lengths(points, tets, TET_EDGES)
    hole_edges = unique_lengths(points, triangles[hole], ((0, 1), (1, 2), (2, 0)))
    hole_xy = unique_lengths(points[:, :2], triangles[hole], ((0, 1), (1, 2), (2, 0)))
    hole_xy = hole_xy[hole_xy > tol]
    hole_owners = np.unique(owners[hole])
    expected_volume = (np.prod((hi-lo)[:2]) - np.pi*radius**2)*(hi-lo)[2]
    if sha256(path) != original_hash:
        raise ValueError("The input mesh changed during the solve")
    signature = output / f"{label}.npz"
    np.savez_compressed(signature, points=points, tetrahedra=tets, triangles=triangles,
                       triangle_tags=tags, displacement=displacement.reshape(-1, 3),
                       tet_centroids=centroids, stress=stress, von_mises=mises, hoop_stress=hoop,
                       exact_stress=exact, radius=radius, remote_stress=remote)
    return dict(
        label=label, archive=str(path), archive_sha256=original_hash, solution=signature.name,
        solution_sha256=sha256(signature), nodes=len(points), dofs=len(displacement),
        tet4=len(tets), surface_tri3=len(triangles), outer_face_tags=side_tags,
        hole_face_tag=hole_tag, hole_wall_tri3=int(hole.sum()), body_tag=next(iter(body_tags)),
        bounds_mm=bounds.tolist(), mesh_volume_mm3=float(volume.sum()),
        ideal_cad_volume_mm3=float(expected_volume),
        relative_volume_error=float(volume.sum()/expected_volume-1),
        edge_p50_mm=float(np.median(edges)), edge_p95_mm=float(np.quantile(edges, .95)),
        hole_edge_p50_mm=float(np.median(hole_edges)),
        hole_circumferential_edge_p50_mm=float(np.median(hole_xy)),
        maximum_hole_wall_cell_hoop_stress_mpa=float(hoop[hole_owners].max()),
        hole_probe=probes[0], radial_probes=probes, physical_probe_strip=strip,
        fixed_annulus=annulus, **error,
        force_balance_n=force_balance.tolist(), moment_balance_n_mm=moment_balance.tolist(),
        rigid_support_reaction_relative=float(support_reaction),
        free_residual_relative=float(free_residual), strain_energy_n_mm=energy,
        applied_work_n_mm=work, assembly_seconds=assembled-start, solve_seconds=solved-assembled)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", nargs=2, metavar=("LABEL", "ARCHIVE"), action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--radius", type=float, default=5.)
    parser.add_argument("--young", type=float, default=210000.)
    parser.add_argument("--poisson", type=float, default=.30)
    parser.add_argument("--remote", type=float, default=100.)
    args = parser.parse_args()
    if not all(np.isfinite([args.radius, args.young, args.poisson, args.remote])) or min(
            args.radius, args.young, args.remote) <= 0 or not -1 < args.poisson < .5:
        parser.error("Provide finite positive radius, modulus and stress and -1 < poisson < 0.5")
    labels = [label for label, _ in args.mesh]
    if len(set(labels)) != len(labels) or not all(re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_-]*", label)
                                                 for label in labels):
        parser.error("Choose distinct path-safe labels")
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for label, archive in args.mesh:
        row = solve(label, Path(archive), args.output, args.radius, args.young, args.poisson, args.remote)
        if rows:
            np.testing.assert_allclose(row["bounds_mm"], rows[0]["bounds_mm"], rtol=0, atol=1e-8)
        rows.append(row)
        local_error = row["relative_near_hole_kirsch_stress_l2"]
        error_display = f"{100*local_error:.2f}%" if local_error is not None else "unresolved"
        print(f"{label}  {row['tet4']:,} Tet4  rim Kt {row['hole_probe']['fem_kt']:.4f}"
              f"  local stress L2 {error_display}")
    result = dict(
        schema="paacs-plate-hole-elasticity-v1",
        geometry=dict(bounds_mm=rows[0]["bounds_mm"], radius_mm=args.radius, hole_center_xy_mm=[0, 0]),
        material=dict(young_mpa=args.young, poisson=args.poisson),
        load=dict(remote_stress_mpa=args.remote, kirsch_peak_mpa=3*args.remote, kirsch_kt=3.,
                  boundary="Kirsch tractions on x/y truncation faces; free hole and z faces",
                  reference="Exact for the ideal circular 3D domain at Poisson ratio zero"
                  if args.poisson == 0 else "Plane-stress reference for finite-thickness 3D elasticity"),
        meshes=rows)
    (args.output / "summary.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n",
                                               encoding="utf-8")


if __name__ == "__main__":
    main()
