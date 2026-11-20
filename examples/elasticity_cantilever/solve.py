# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Solve an OCCT STEP cantilever on PAACS Tet4 archives."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from mesh_visual_data import read_archive, sha256

LENGTH, WIDTH, HEIGHT = 100., 10., 10.  # mm
YOUNG, POISSON, FORCE = 210_000., .3, 100.  # N/mm², dimensionless, N


def constitutive_matrix():
    shear = YOUNG / (2*(1+POISSON))
    lame = YOUNG*POISSON / ((1+POISSON)*(1-2*POISSON))
    elasticity = np.diag([2*shear]*3+[shear]*3)
    elasticity[:3, :3] += lame
    return elasticity


def assemble(points, tets, elasticity):
    positions = points[tets]
    jacobian = (positions[:, 1:] - positions[:, :1]).transpose(0, 2, 1)
    det = np.linalg.det(jacobian)
    if np.any(det <= 0):
        raise ValueError("The archive contains a nonpositive Tet4")
    volume = det/6
    gradients = np.empty((len(tets), 4, 3))
    gradients[:, 1:] = np.linalg.inv(jacobian)
    gradients[:, 0] = -gradients[:, 1:].sum(axis=1)
    strain = np.zeros((len(tets), 6, 12))
    dx, dy, dz = gradients.transpose(2, 0, 1)
    # Engineering shear strain ordering is xy, yz, xz.
    strain[:, 0, 0::3], strain[:, 1, 1::3], strain[:, 2, 2::3] = dx, dy, dz
    strain[:, 3, 0::3], strain[:, 3, 1::3] = dy, dx
    strain[:, 4, 1::3], strain[:, 4, 2::3] = dz, dy
    strain[:, 5, 0::3], strain[:, 5, 2::3] = dz, dx
    stiffness = ((strain.transpose(0, 2, 1) @ elasticity) @ strain) * volume[:, None, None]
    dofs = (3*tets[:, :, None] + np.arange(3)).reshape(-1, 12)
    rows, cols = np.broadcast_arrays(dofs[:, :, None], dofs[:, None, :])
    ndofs = 3*len(points)
    matrix = coo_matrix((stiffness.ravel(), (rows.ravel(), cols.ravel())), shape=(ndofs, ndofs)).tocsc()
    matrix.sum_duplicates()
    return matrix, strain, dofs, volume


def apply_constraints(matrix, load, constrained, prescribed):
    free = np.ones(matrix.shape[0], dtype=bool)
    free[constrained] = False
    displacement = np.zeros(matrix.shape[0])
    displacement[constrained] = prescribed
    effective_load = load[free] - matrix[free][:, constrained] @ displacement[constrained]
    displacement[free] = spsolve(matrix[free][:, free], effective_load)
    if not np.isfinite(displacement).all():
        raise ValueError("The elasticity system did not produce a finite solution")
    return displacement, free


def patch_test(matrix, strain, dofs, points, surface):
    gradient = 1e-4*np.array([[.4, -.2, .8], [.6, -.4, .1], [-.3, .2, .7]])
    exact = (points @ gradient.T).ravel()
    boundary_dofs = (3*np.unique(surface)[:, None] + np.arange(3)).ravel()
    displacement, free = apply_constraints(matrix, np.zeros_like(exact), boundary_dofs, exact[boundary_dofs])
    if not free.any():
        raise ValueError("The patch test requires interior nodes")
    difference = np.linalg.norm(displacement-exact)/np.linalg.norm(exact)
    gradient_strain = np.array([gradient[0, 0], gradient[1, 1], gradient[2, 2],
                                gradient[0, 1]+gradient[1, 0], gradient[1, 2]+gradient[2, 1],
                                gradient[0, 2]+gradient[2, 0]])
    recovered = (strain @ displacement[dofs, None]).squeeze(-1)
    strain_error = np.max(np.abs(recovered-gradient_strain))
    if difference > 1e-9 or strain_error > 1e-10:
        raise ValueError("Affine displacement patch test failed")
    return dict(relative_displacement_l2=float(difference), max_absolute_strain_error=float(strain_error),
                interior_dofs=int(free.sum()))


def cap(points, surface, tags, x):
    selected = np.all(np.isclose(points[surface, 0], x, rtol=0, atol=1e-9), axis=1)
    cap_tags = np.unique(tags[selected])
    if len(cap_tags) != 1 or not np.array_equal(selected, tags == cap_tags[0]):
        raise ValueError(f"Expected one complete CAD-tagged cap at x={x:g}")
    tris = surface[selected]
    a, b, c = points[tris].transpose(1, 0, 2)
    areas = .5*np.linalg.norm(np.cross(b-a, c-a), axis=1)
    np.testing.assert_allclose(areas.sum(), WIDTH*HEIGHT, rtol=1e-12)
    return tris, int(cap_tags[0]), areas


def solve(path, output):
    mesh = read_archive(path)
    points = mesh.points
    np.testing.assert_allclose(points.min(axis=0), [0., 0., 0.], atol=1e-9)
    np.testing.assert_allclose(points.max(axis=0), [LENGTH, WIDTH, HEIGHT], atol=1e-9)
    volume_cells = [cell for cell in mesh.cells if cell.dim == 3]
    surface_cells = [cell for cell in mesh.cells if cell.dim == 2]
    if not volume_cells or not surface_cells or any(cell.kind != 5 for cell in volume_cells) or any(
            cell.kind != 2 for cell in surface_cells):
        raise ValueError("Expected Tri3 exterior and Tet4 volume elements")
    tets = np.array([cell.ids for cell in volume_cells])
    surface = np.array([cell.ids for cell in surface_cells])
    tags = np.array([cell.tag for cell in surface_cells])
    clamp, clamp_tag, clamp_areas = cap(points, surface, tags, 0.)
    tip, tip_tag, tip_areas = cap(points, surface, tags, LENGTH)
    matrix, strain, dofs, volumes = assemble(points, tets, constitutive_matrix())
    np.testing.assert_allclose(volumes.sum(), LENGTH*WIDTH*HEIGHT, rtol=1e-12)

    weights = np.zeros(len(points))
    np.add.at(weights, tip.ravel(), np.repeat(tip_areas/3, 3))
    load = np.zeros((len(points), 3))
    load[:, 2] = -FORCE/(WIDTH*HEIGHT)*weights
    fixed = (3*np.unique(clamp)[:, None] + np.arange(3)).ravel()
    flat_load = load.ravel()
    displacement, free = apply_constraints(matrix, flat_load, fixed, 0.)
    residual = matrix @ displacement - flat_load
    reaction = np.zeros_like(residual)
    reaction[fixed] = residual[fixed]
    reaction = reaction.reshape(-1, 3)
    force_balance = (load+reaction).sum(axis=0)
    moment_balance = np.cross(points, load+reaction).sum(axis=0)
    np.testing.assert_allclose(force_balance, 0, atol=1e-6)
    np.testing.assert_allclose(moment_balance, 0, atol=1e-4)
    free_residual = np.linalg.norm(residual[free])/np.linalg.norm(flat_load)
    if free_residual > 1e-9:
        raise ValueError("Cantilever free-DOF residual exceeded tolerance")
    displacement_xyz = displacement.reshape(-1, 3)
    deflection = float(-(weights @ displacement_xyz[:, 2])/weights.sum())
    work = float(flat_load @ displacement)
    energy = float(.5*(displacement @ (matrix @ displacement)))
    np.testing.assert_allclose(work, 2*energy, rtol=1e-9)

    moment_of_inertia = WIDTH*HEIGHT**3/12
    euler = FORCE*LENGTH**3/(3*YOUNG*moment_of_inertia)
    shear = YOUNG/(2*(1+POISSON))
    timoshenko = euler + FORCE*LENGTH/((5/6)*shear*WIDTH*HEIGHT)
    edges = np.unique(np.sort(tets[:, ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))].reshape(-1, 2),
                              axis=1), axis=0)
    lengths = np.linalg.norm(points[edges[:, 1]]-points[edges[:, 0]], axis=1)
    result = dict(
        archive=str(path), archive_sha256=sha256(path), nodes=len(points), tet4=len(tets), surface_tri3=len(surface),
        edge_p50_mm=float(np.median(lengths)), edge_p95_mm=float(np.quantile(lengths, .95)),
        mesh_volume_mm3=float(volumes.sum()), solid_tags=sorted({cell.tag for cell in volume_cells}),
        fixed_face_tag=clamp_tag, loaded_face_tag=tip_tag, fixed_area_mm2=float(clamp_areas.sum()),
        loaded_area_mm2=float(tip_areas.sum()), tip_deflection_mm=deflection,
        euler_bernoulli_mm=euler, timoshenko_mm=timoshenko,
        error_vs_euler_percent=100*(deflection/euler-1), error_vs_timoshenko_percent=100*(deflection/timoshenko-1),
        reaction_force_n=reaction.sum(axis=0).tolist(), free_residual_relative=float(free_residual),
        force_balance_n=force_balance.tolist(), moment_balance_n_mm=moment_balance.tolist(),
        strain_energy_n_mm=energy, applied_work_n_mm=work,
        affine_patch=patch_test(matrix, strain, dofs, points, surface),
    )
    output.mkdir(parents=True, exist_ok=True)
    (output / "solve.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    np.savez_compressed(output / "displacement.npz", displacement_mm=displacement_xyz)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(solve(args.archive.resolve(), args.output.resolve()), allow_nan=False))


if __name__ == "__main__":
    main()
