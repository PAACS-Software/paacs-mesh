# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Save native STEP sizing selections and replay them through paacs-cad-mesh."""

import argparse
from collections import defaultdict
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from examples.elasticity_plate_hole.plate_controls import assert_export_matches, groups, save_json, stats
from examples.elasticity_plate_hole.solver_input import read_solver_input, write_solver_input
from scripts.mesh_visual_data import TET_EDGES, TET_FACES, read_archive, section, sha256


SCHEMA = "paacs-native-step-sizing-v1"
VARIANTS = ("global", "face", "box", "combined")


def save_project(args):
    args.project.mkdir(parents=True)
    step = args.project / "source.step"
    shutil.copyfile(args.step, step)
    save_json(args.project / "project.json", {
        "schema": SCHEMA, "source": step.name, "source_sha256": sha256(step),
        "global_mm": args.size, "minimum_mm": args.minimum, "growth": args.growth,
        "face_tag": args.face_tag, "face_mm": args.face_size,
        "box_mm": args.box, "chord_mm": args.chord,
    })


def load_project(project):
    recipe = json.loads((project / "project.json").read_text(encoding="utf-8"))
    if recipe.get("schema") != SCHEMA or recipe.get("source") != "source.step":
        raise ValueError("Invalid native STEP sizing project")
    step = project / recipe["source"]
    if sha256(step) != recipe["source_sha256"]:
        raise ValueError("STEP source hash changed. Face tags need to be selected again.")
    return step, recipe


def mesh_variant(executable, step, recipe, output, variant):
    output.mkdir(parents=True)
    archive = output / "mesh.paacsmesh"
    command = [str(executable.resolve()), str(step.resolve()), str(recipe["global_mm"]),
               str(archive.resolve()), "--native-cad", "--native-sizing", "--report-json",
               "--trace-native-stages", "--threads", "1", "--native-validation-threads", "1",
               "--recovery-mode", "advanced", "--reuse-native-queries", "--max-work", "100000000",
               "--min-size", str(recipe["minimum_mm"]), "--growth", str(recipe["growth"])]
    if recipe["chord_mm"] is not None:
        command += ["--max-chord-error", str(recipe["chord_mm"])]
    if variant in ("face", "combined"):
        command += ["--face-size", str(recipe["face_tag"]), str(recipe["face_mm"])]
    if variant in ("box", "combined"):
        command += ["--size-box", *map(str, recipe["box_mm"])]
    save_json(output / "command.json", command)
    run = subprocess.run(command, capture_output=True, text=True, timeout=90)
    (output / "stdout.log").write_text(run.stdout, encoding="utf-8")
    (output / "stages.log").write_text(run.stderr, encoding="utf-8")
    report = next((line for line in reversed(run.stdout.splitlines()) if line.startswith("{")), None)
    if report is None:
        raise RuntimeError(f"{variant} produced no native receipt. See {output / 'stages.log'}")
    receipt = json.loads(report)
    save_json(output / "receipt.json", receipt)
    if (run.returncode != 0 or not receipt["archive_published"] or not receipt["audit_valid"]
            or receipt["bodies"][0]["route"] != "native_spatial_sizing"):
        raise RuntimeError(f"{variant} failed native sizing or audit. See {output / 'receipt.json'}")
    mesh = read_archive(archive)
    measures, _ = measure(mesh, recipe)
    body = receipt["bodies"][0]
    if (not body["boundary_embedding_certified"] or body["centroid_outside_distance_limit"] != 0
            or body["native_centroids_checked"] != measures["tet4"]):
        raise RuntimeError(f"{variant} did not complete strict CAD checks")
    include = output / "mesh.inp"
    write_solver_input(mesh, include)
    assert_export_matches(mesh, read_solver_input(include))
    measures.update(tags=groups(mesh), mesh_sha256=sha256(archive), solver_roundtrip=True)
    save_json(output / "measurements.json", measures)


def measure(mesh, recipe):
    volume_ids = [i for i, cell in enumerate(mesh.cells) if cell.dim == 3]
    tets = np.array([mesh.cells[i].ids for i in volume_ids])
    xyz = mesh.points[tets]
    centroids = xyz.mean(axis=1)
    edges = np.array([np.linalg.norm(xyz[:, a] - xyz[:, b], axis=1) for a, b in TET_EDGES]).T
    size = np.median(edges, axis=1)
    lo, hi = np.array(recipe["box_mm"][:3]), np.array(recipe["box_mm"][3:6])
    in_box = np.all((centroids >= lo) & (centroids <= hi), axis=1)
    box_distance = np.linalg.norm(np.maximum(0., np.maximum(lo - centroids, centroids - hi)), axis=1)
    shell_edges = set()
    for cell in mesh.cells:
        if cell.dim == 2 and cell.tag == recipe["face_tag"]:
            shell_edges.update(tuple(sorted((cell.ids[i], cell.ids[(i+1) % 3]))) for i in range(3))
    if not shell_edges:
        raise ValueError("Selected native face tag is missing from the mesh")
    a, b = np.array(sorted(shell_edges)).T
    selected_edges = np.linalg.norm(mesh.points[a] - mesh.points[b], axis=1)
    use = defaultdict(list)
    for i, tet in enumerate(tets):
        for face in TET_FACES:
            use[tuple(sorted(tet[list(face)]))].append(i)
    pairs = np.array([p for p in use.values() if len(p) == 2], dtype=int).reshape(-1, 2)
    pair_sizes = size[pairs]
    ratios = pair_sizes.max(axis=1) / pair_sizes.min(axis=1)
    near_box = box_distance <= recipe["global_mm"]*.5
    return {
        "nodes": len(mesh.points), "tet4": len(tets), "selected_face_edges_mm": stats(selected_edges),
        "box_tet_median_edges_mm": stats(size[in_box]),
        "outside_box_tet_median_edges_mm": stats(size[~in_box]),
        "volume_tet_median_edges_mm": stats(size),
        "adjacent_tet_size_ratio": stats(ratios),
        "near_box_adjacent_tet_size_ratio": stats(ratios[np.any(near_box[pairs], axis=1)]),
    }, (volume_ids, size)


def plot_sections(project, output, axis, coordinate):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    from matplotlib.patches import Rectangle

    _, recipe = load_project(project)
    variants = {}
    for variant in VARIANTS:
        mesh = read_archive(output / variant / "mesh.paacsmesh")
        metrics, (tets, sizes) = measure(mesh, recipe)
        polygons, owners, _ = section(mesh, tets, None, axis, coordinate)
        if not len(owners):
            raise ValueError(f"The requested plane misses the {variant} mesh")
        variants[variant] = (metrics, sizes[owners], polygons)
    scale = Normalize(min(row[1].min() for row in variants.values()),
                      np.quantile(np.concatenate([row[1] for row in variants.values()]), .98))
    fig, axes = plt.subplots(2, 2, figsize=(12, 11), sharex=True, sharey=True,
                             layout="constrained", facecolor="#f2f5f8")
    xy = [i for i in range(3) if i != axis]
    lo = np.array(recipe["box_mm"][:3])
    hi = np.array(recipe["box_mm"][3:6])
    cut = None
    for panel, (variant, (metrics, sizes, polygons)) in zip(axes.flat, variants.items()):
        cut = PolyCollection([p[:, xy] for p in polygons], array=sizes, norm=scale,
                             cmap="viridis", edgecolors="#173343", linewidths=.27)
        panel.add_collection(cut)
        if lo[axis] <= coordinate <= hi[axis]:
            panel.add_patch(Rectangle(lo[xy], *(hi-lo)[xy], fill=False, edgecolor="#ff6b4a",
                                      linestyle="--", linewidth=1.7, zorder=4))
        panel.autoscale()
        panel.set_aspect("equal")
        panel.set_xlabel(f'{"xyz"[xy[0]]} (mm)')
        panel.set_ylabel(f'{"xyz"[xy[1]]} (mm)')
        title = {"global": "Global only", "face": "Curved face only",
                 "box": "Refinement box only", "combined": "Face + box"}[variant]
        panel.set_title(title, loc="left", fontweight="bold", color="#173343", pad=19)
        face = metrics["selected_face_edges_mm"]["median"]
        box = metrics["box_tet_median_edges_mm"]["median"]
        box_label = f"{box:.3f} mm" if box is not None else "no centroids"
        growth = metrics["near_box_adjacent_tet_size_ratio"]["p95"]
        growth_label = f"{growth:.2f}×" if growth is not None else "n/a"
        panel.text(0, 1.013, f'{metrics["tet4"]:,} Tet4  ·  Face {face:.3f} mm  ·  Box {box_label}'
                   f'  ·  Near-box growth p95 {growth_label}',
                   transform=panel.transAxes, fontsize=8.5, color="#43596d")
    fig.colorbar(cut, ax=axes, label="Actual Tet4 median edge (mm)", extend="max", shrink=.85)
    fig.suptitle(f'Native STEP sizing at {"xyz"[axis]} = {coordinate:g} mm\n'
                 f'Exact mesh sections · dashed coral outline marks the refinement box · face {recipe["face_tag"]}',
                 fontsize=15, fontweight="bold", color="#173343")
    fig.savefig(output / "sizing-sections.png", dpi=190)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    save = commands.add_parser("save")
    save.add_argument("step", type=Path)
    save.add_argument("project", type=Path)
    save.add_argument("--size", type=float, required=True)
    save.add_argument("--minimum", type=float, default=0.)
    save.add_argument("--growth", type=float, default=1.3)
    save.add_argument("--face-tag", type=int, required=True)
    save.add_argument("--face-size", type=float, required=True)
    save.add_argument("--box", nargs=7, type=float, required=True, metavar="N")
    save.add_argument("--chord", type=float)
    replay = commands.add_parser("mesh")
    replay.add_argument("project", type=Path)
    replay.add_argument("executable", type=Path)
    replay.add_argument("output", type=Path)
    replay.add_argument("--variant", choices=VARIANTS, default="combined")
    plot = commands.add_parser("plot")
    plot.add_argument("project", type=Path)
    plot.add_argument("output", type=Path)
    plot.add_argument("--axis", choices="xyz", required=True)
    plot.add_argument("--coordinate", type=float, required=True)
    args = parser.parse_args()
    if args.action == "save":
        save_project(args)
    elif args.action == "mesh":
        step, recipe = load_project(args.project)
        mesh_variant(args.executable, step, recipe, args.output, args.variant)
    else:
        plot_sections(args.project, args.output, "xyz".index(args.axis), args.coordinate)


if __name__ == "__main__":
    main()
