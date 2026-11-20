# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Mesh a saved STEP-plate sizing recipe and verify its C3D4 solver export."""

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from examples.elasticity_plate_hole.solver_input import read_solver_input, write_solver_input
from scripts.mesh_visual_data import TET_EDGES, TET_FACES, read_archive, sha256


VARIANTS = ("global", "bore", "box", "combined")


def save_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def groups(mesh):
    count = Counter((cell.dim, cell.tag) for cell in mesh.cells)
    return {f"{dimension}:{tag}": n for (dimension, tag), n in sorted(count.items())}


def assert_export_matches(native, solver):
    np.testing.assert_array_equal(native.points, solver.points)
    for dimension in (2, 3):
        # CAD surface faces may be written in a different order. Preserve
        # orientation under cyclic rotation. Compare each Tet4 as written.
        def records(mesh):
            for cell in mesh.cells:
                if cell.dim == dimension:
                    ids = cell.ids
                    if dimension == 2:
                        start = ids.index(min(ids))
                        ids = ids[start:] + ids[:start]
                    yield cell.tag, ids
        if Counter(records(native)) != Counter(records(solver)):
            raise ValueError(f"Solver mesh changed dimension-{dimension} tagged connectivity")


def stats(values):
    if len(values) == 0:
        return {"count": 0, "median": None, "p95": None, "maximum": None}
    median, p95 = np.quantile(values, [.5, .95])
    return dict(count=len(values), median=float(median), p95=float(p95), maximum=float(values.max()))


def measures(mesh, recipe, radius, bore_tag):
    box = recipe["box"]
    lo, hi = np.array(box["lower_mm"]), np.array(box["upper_mm"])
    tets = np.array([cell.ids for cell in mesh.cells if cell.kind == 5], dtype=np.int64)
    edges = np.unique(np.sort(tets[:, TET_EDGES].reshape(-1, 2), axis=1), axis=0)
    a, b = mesh.points[edges].transpose(1, 0, 2)
    mid = (a+b)/2
    lengths = np.linalg.norm(b-a, axis=1)
    planar = np.isclose(a[:, 2], b[:, 2], atol=1e-10, rtol=0)
    inside = np.all((lo <= mid) & (mid <= hi), axis=1)
    near = np.linalg.norm(mid[:, :2], axis=1)
    bore_faces = np.array([cell.ids for cell in mesh.cells if cell.dim == 2 and cell.tag == bore_tag])
    bore_edges = np.unique(np.sort(bore_faces[:, ((0, 1), (1, 2), (2, 0))].reshape(-1, 2), axis=1), axis=0)
    a, b = mesh.points[bore_edges].transpose(1, 0, 2)
    rim = np.isclose(a[:, 2], b[:, 2], atol=1e-10, rtol=0)
    rim_lengths = np.linalg.norm(b-a, axis=1)[rim]
    rim_mid = ((a+b)/2)[rim]
    rim_in_box = np.all((lo <= rim_mid) & (rim_mid <= hi), axis=1)

    points = mesh.points[tets]
    signed_volume = np.linalg.det(points[:, 1:] - points[:, :1])/6
    if np.any(signed_volume <= 0):
        raise ValueError("Plate archive contains an invalid cell")
    uses = np.sort(tets[:, TET_FACES].reshape(-1, 3), axis=1)
    _, owners, counts = np.unique(uses, axis=0, return_inverse=True, return_counts=True)
    order = np.argsort(owners.ravel(), kind="stable")
    starts = np.cumsum(counts)-counts
    adjacent = starts[counts == 2]
    pairs = signed_volume[order[adjacent[:, None] + [0, 1]]//4]
    growth = np.cbrt(pairs.max(axis=1)/pairs.min(axis=1))
    return dict(rim_overlapping_box_mm=stats(rim_lengths[rim_in_box]),
                rim_other_mm=stats(rim_lengths[~rim_in_box]),
                left_bore_band_in_plane_mm=stats(lengths[planar & (mid[:, 0] < 0) & (near > radius) &
                                                             (near < radius+2)]),
                box_interior_in_plane_mm=stats(lengths[planar & inside & (near >= radius+3)]),
                far_in_plane_mm=stats(lengths[planar & (abs(mid[:, 0]) >= 30)]),
                adjacent_equivalent_size_ratio=stats(growth))


def mesh_variant(step, executable, output, geometry, sizing, variant):
    directory = output/variant
    directory.mkdir()
    archive = directory/"mesh.paacsmesh"
    deck = directory/"mesh.inp"
    bore, box = sizing["bore"], sizing["box"]
    bore_enabled = variant in ("bore", "combined")
    box_enabled = variant in ("box", "combined")
    arguments = [geometry["body_tag"], bore["face_tag"] if bore_enabled else 0,
                 sizing["global_mm"], sizing["minimum_mm"], sizing["growth_ratio"],
                 bore["target_mm"] if bore_enabled else 0, box["target_mm"] if box_enabled else 0,
                 *box["lower_mm"], *box["upper_mm"], sizing["layers"]]
    command = [str(executable), str(step), str(archive), "--controls", *map(str, arguments)]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    run = subprocess.run(command, env=env, text=True, encoding="utf-8", errors="replace",
                         capture_output=True, timeout=90)
    (directory/"mesh.log").write_text(run.stdout+run.stderr, encoding="utf-8")
    if run.returncode:
        raise ValueError(f"PAACS refused {variant}: {run.stderr.strip()}")
    report = json.loads(run.stdout.splitlines()[-1])
    if (report["outside_centroids"] or report["body_tag"] != geometry["body_tag"] or
            report["bore_tag"] != bore["face_tag"]):
        raise ValueError("Mesher report disagrees with the saved CAD body or bore selection")
    native = read_archive(archive)
    expected_groups = {f"2:{tag}" for tag in report["original_face_tags"]} | {f"3:{report['body_tag']}"}
    if groups(native).keys() != expected_groups:
        raise ValueError("Native archive lost a STEP boundary or body group")
    write_solver_input(native, deck)
    imported = read_solver_input(deck)
    assert_export_matches(native, imported)
    return dict(variant=variant, archive_sha256=sha256(archive), solver_input_sha256=sha256(deck),
                cad_groups=groups(native), solver_groups=groups(imported), **report,
                **measures(native, sizing, report["plate_radius_mm"], report["bore_tag"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project", type=Path)
    parser.add_argument("--mesher", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--ablation", action="store_true", help="Also mesh global, bore-only and box-only controls")
    args = parser.parse_args()
    project = json.loads(args.project.read_text(encoding="utf-8-sig"))
    if project.get("schema") != "paacs-plate-sizing-v1" or project["geometry"]["unit"] != "mm":
        raise ValueError("Expected a PAACS STEP plate project in millimetres")
    geo = project["geometry"]
    step = args.project.resolve().parent/geo["step"]
    if sha256(step) != geo["sha256"]:
        raise ValueError("STEP SHA-256 does not match the saved project selections")
    output = args.output.resolve()
    output.mkdir(parents=True)
    saved_step = output/"source.step"
    shutil.copyfile(step, saved_step)
    if sha256(saved_step) != geo["sha256"]:
        raise ValueError("The STEP changed while saving the project")
    geo["step"] = saved_step.name
    save_json(output/"project.json", project)
    executable = args.mesher.resolve()
    mesher_hash = sha256(executable)
    rows = [mesh_variant(saved_step, executable, output, geo, project["sizing"], variant)
            for variant in (VARIANTS if args.ablation else ("combined",))]
    if sha256(executable) != mesher_hash or sha256(saved_step) != geo["sha256"]:
        raise ValueError("A mesher input changed during the plate study")
    save_json(output/"metrics.json", dict(step_sha256=geo["sha256"], mesher_sha256=mesher_hash,
                                         sizing=project["sizing"], meshes=rows))
    for row in rows:
        print(f"{row['variant']}: {row['tet4']:,} Tet4, {len(row['cad_groups'])} preserved CAD groups")


if __name__ == "__main__":
    main()
