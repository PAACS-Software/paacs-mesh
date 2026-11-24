"""Render measured Tet4 whole-cell cutaways without importing the native engine.

Requires NumPy and Matplotlib. Archive inputs are explicit; stored cells and
their zero-based archive IDs are retained. This is not a CAD or solver audit.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import struct
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


ROOT = Path(__file__).resolve().parents[1]
BG, INK, SUB = "#f4f7fa", "#17364b", "#526b7c"
FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3))


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    import sys
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


ARCHIVES = module("qualification_archive_reader", ROOT / "scripts/mesh_visual_data.py")
QUALITY = module("qualification_quality", ROOT / "python/paacs_mesh/quality_metrics.py")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def artifact(path):
    return {"path": str(path), "sha256": ARCHIVES.sha256(path)}


def point_hashes(mesh):
    points = np.asarray(mesh.points, dtype="<f8")
    order = np.lexsort((points[:, 2], points[:, 1], points[:, 0]))
    return {"nodes": len(points), "ordered_xyz_sha256": sha(points.tobytes()),
            "sorted_xyz_sha256": sha(points[order].tobytes())}


def surface_hashes(mesh):
    unoriented, oriented = [], []
    tags = Counter()
    for cell in mesh.cells:
        if cell.kind != 2 or cell.dim != 2:
            continue
        vertices = [struct.pack("<3d", *mesh.points[node]) for node in cell.ids]
        tag = struct.pack("<i", cell.tag)
        unoriented.append(tag + b"".join(sorted(vertices)))
        oriented.append(tag + min(b"".join(vertices[i:] + vertices[:i]) for i in range(3)))
        tags[cell.tag] += 1
    return {"triangles": len(unoriented), "face_tag_counts": dict(sorted(tags.items())),
            "unoriented_geometry_sha256": sha(b"".join(sorted(unoriented))),
            "oriented_geometry_sha256": sha(b"".join(sorted(oriented))),
            "method": "Exact stored binary64 coordinates and face tags; cyclic rotation preserves orientation"}


def summarize(values):
    if not len(values) or not np.isfinite(values).all():
        raise ValueError("Missing or nonfinite element measurements")
    quantiles = np.quantile(values, [0, .01, .05, .5, .95, 1])
    return dict(zip(("minimum", "p01", "p05", "median", "p95", "maximum"), map(float, quantiles)))


def prepare(path, label, args, started):
    print(json.dumps({"event": "measuring_archive", "path": str(path)}), flush=True)
    before = ARCHIVES.sha256(path)
    mesh = ARCHIVES.read_archive(path, max_bytes=args.max_bytes,
                                 max_nodes=args.max_nodes, max_cells=args.max_cells)
    ids = np.asarray([i for i, cell in enumerate(mesh.cells) if cell.kind == 5], dtype=np.int64)
    if any(cell.dim == 3 and cell.kind != 5 for cell in mesh.cells) or not len(ids):
        raise ValueError("Figure supports stored Tet4 volume cells only")
    connectivity = np.asarray([mesh.cells[i].ids for i in ids], dtype=np.int64)
    points = mesh.points[connectivity]
    metric_names = ("mean_ratio", "minimum_dihedral_deg", "scaled_jacobian")
    measurements = {name: np.empty(len(ids)) for name in metric_names}
    edge_max = np.empty(len(ids))
    for index, corners in enumerate(points):
        if index % 512 == 0 and time.perf_counter() - started > args.max_seconds:
            raise ValueError("Renderer measurement deadline exhausted")
        metrics = QUALITY.tet4_metrics(corners)
        if metrics["orientation"] != "positive" or metrics["range_limited"]:
            raise ValueError(f"Archive cell {ids[index]} has nonpositive or range-limited measurements")
        for name in metric_names:
            measurements[name][index] = metrics[name]
        edge_max[index] = metrics["edge_max"]
    # Drop complete cells whose centroid lies in the near y half. No point,
    # triangle or tetrahedron is clipped, moved, retriangulated or fabricated.
    center_y = float((mesh.points[:, 1].min() + mesh.points[:, 1].max()) / 2)
    keep = points.mean(axis=1)[:, 1] >= center_y
    exterior = {}
    for local in np.flatnonzero(keep):
        for face in FACES:
            face_ids = tuple(int(connectivity[local, node]) for node in face)
            key = tuple(sorted(face_ids))
            if key in exterior:
                exterior[key][2] += 1
            else:
                exterior[key] = [face_ids, int(local), 1]
    if any(value[2] > 2 for value in exterior.values()):
        raise ValueError("Retained cells have nonmanifold face incidence")
    faces = sorted((key, value) for key, value in exterior.items() if value[2] == 1)
    if not faces or len(faces) > args.max_faces:
        raise ValueError("Exposed-face output limit exceeded or cutaway empty")
    polygons = [mesh.points[list(value[0])] for _, value in faces]
    owner = np.asarray([value[1] for _, value in faces])
    quality = measurements["mean_ratio"]
    ranking = np.argsort(quality, kind="stable")
    worst = [{"cell_id": int(ids[i]), "mean_ratio": float(quality[i]),
              "minimum_dihedral_deg": float(measurements["minimum_dihedral_deg"][i]),
              "scaled_jacobian": float(measurements["scaled_jacobian"][i]),
              "connectivity": connectivity[i].tolist(), "centroid": points[i].mean(axis=0).tolist(),
              "retained_in_cutaway": bool(keep[i])} for i in ranking[:8]]
    mapping = [(int(ids[value[1]]), *key) for key, value in faces]
    receipt = {"label": label, "archive": artifact(path), "point_hashes": point_hashes(mesh),
        "surface_hashes": surface_hashes(mesh), "tetrahedra": len(ids),
        "body_tet4_counts": dict(Counter(mesh.cells[i].tag for i in ids)),
        "quality": {name: summarize(values) for name, values in measurements.items()},
        "cell_max_edge": summarize(edge_max), "mean_ratio_below_0p05": int(np.count_nonzero(quality < .05)),
        "mean_ratio_below_1e_8": int(np.count_nonzero(quality < 1e-8)),
        "cell_ids": "Zero-based index in the full archive cell sequence, including boundary cells",
        "worst_cells": worst, "cutaway": {"axis": "y", "value": center_y, "retained": "cell centroid >= value",
            "retained_tetrahedra": int(np.count_nonzero(keep)), "exposed_triangles": len(faces),
            "retained_cell_ids_sha256": sha(np.asarray(ids[keep], dtype="<i8").tobytes()),
            "face_owner_ids_sha256": sha(np.asarray(mapping, dtype="<i8").tobytes()),
            "face_mapping_rule": "Sorted node-ID triplets of exposed original faces paired with owning archive cell ID"}}
    if before != ARCHIVES.sha256(path):
        raise ValueError("Archive changed during figure preparation")
    return {"mesh": mesh, "polygons": polygons, "face_quality": quality[owner], "receipt": receipt}


def draw(ax, item, bounds, title):
    ax.set_facecolor(BG)
    colors = plt.get_cmap("RdYlBu")(Normalize(0, 1)(item["face_quality"]))
    ax.add_collection3d(Poly3DCollection(item["polygons"], facecolors=colors,
        edgecolors="#28444f88", linewidths=.21, zsort="average"))
    lo, hi = bounds
    padding = np.max(hi - lo) * .025
    ax.set(xlim=(lo[0]-padding, hi[0]+padding), ylim=(lo[1]-padding, hi[1]+padding),
           zlim=(lo[2]-padding, hi[2]+padding))
    ax.set_box_aspect(hi - lo)
    ax.view_init(elev=24, azim=-60)
    ax.set_proj_type("ortho")
    ax.set_axis_off()
    ax.set_title(title, fontsize=15, color=INK, weight="bold", pad=3)
    receipt = item["receipt"]
    q = receipt["quality"]["mean_ratio"]
    ax.text2D(.02, -.035,
        f"{receipt['tetrahedra']:,} Tet4  |  MR min {q['minimum']:.3g}  ·  p05 {q['p05']:.3f}  ·  median {q['median']:.3f}\n"
        f"{receipt['mean_ratio_below_0p05']:,} below 0.05  |  global minimum cell #{receipt['worst_cells'][0]['cell_id']}",
        transform=ax.transAxes, fontsize=10.5, color=SUB, linespacing=1.6)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("sphere-baseline", "sphere-coarse", "sphere-fine"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("cylinder-before", "cylinder-after", "initial-before", "initial-after", "cylinder-identity", "cylinder-result"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--max-bytes", type=int, default=64*1024**2)
    parser.add_argument("--max-nodes", type=int, default=100_000)
    parser.add_argument("--max-cells", type=int, default=100_000)
    parser.add_argument("--max-faces", type=int, default=100_000)
    parser.add_argument("--max-seconds", type=float, default=120.)
    args = parser.parse_args()
    if bool(args.cylinder_before) != bool(args.cylinder_after) or bool(args.initial_before) != bool(args.initial_after):
        parser.error("Both archives are required for each before/after pair")
    if min(args.max_bytes, args.max_nodes, args.max_cells, args.max_faces, args.max_seconds) <= 0:
        parser.error("Resource limits must be positive")
    started = time.perf_counter()
    inputs = [(args.sphere_baseline, "Sphere: previous fan"),
              (args.sphere_coarse, "Sphere: source retry, h = 0.5"),
              (args.sphere_fine, "Sphere: source retry, h = 0.25")]
    if args.cylinder_before:
        inputs += [(args.cylinder_before, "Cylinder: rebuilt topology"),
                   (args.cylinder_after, "Cylinder: retained topology")]
    items = [prepare(path, label, args, started) for path, label in inputs]
    provenance = {"schema": "paacs-quality-qualification-figure-v1",
        "scope": "Actual archived Tet4 faces and independent shape measurements; no native engine or CAD audit called",
        "metric_conventions": QUALITY.METRIC_CONVENTIONS,
        "camera": {"projection": "orthographic", "elevation": 24, "azimuth": -60},
        "color_scale": {"metric": "mean_ratio", "minimum": 0, "maximum": 1, "colormap": "RdYlBu"},
        "cutaway": "Whole Tet4 cells retained by a centroid half-space; original exposed faces only",
        "quality_scope": "All Tet4 cells, including removed and hidden cells; minimum-cell labels need not be visible faces",
        "timing_comparison_valid": False, "panels": [item["receipt"] for item in items]}
    if args.initial_before:
        before = ARCHIVES.read_archive(args.initial_before, max_bytes=args.max_bytes, max_nodes=args.max_nodes, max_cells=args.max_cells)
        after = ARCHIVES.read_archive(args.initial_after, max_bytes=args.max_bytes, max_nodes=args.max_nodes, max_cells=args.max_cells)
        bh, ah = point_hashes(before), point_hashes(after)
        provenance["initial_point_comparison"] = {"before": artifact(args.initial_before),
            "after": artifact(args.initial_after), "before_hashes": bh, "after_hashes": ah,
            "exact_ordered_xyz_match": bh == ah}
    if args.cylinder_before:
        before, after = items[-2:]
        provenance["final_boundary_comparison"] = {
            "exact_oriented_coordinate_and_tag_match": before["receipt"]["surface_hashes"] == after["receipt"]["surface_hashes"],
            "scope": "Stored Tri3 CAD-boundary cells; no CAD reevaluation"}
    if args.cylinder_identity:
        identity = read_json(args.cylinder_identity)
        if identity.get("reference_archive_sha256") != items[-2]["receipt"]["archive"]["sha256"]:
            raise ValueError("Cylinder reference archive differs from probe identity")
        provenance["cylinder_run_identity"] = {"artifact": artifact(args.cylinder_identity), "record": identity}
    if args.cylinder_result:
        result = read_json(args.cylinder_result)
        if result.get("published") is not True or result.get("reason") != "none" or result.get("probe_failed"):
            raise ValueError("Cylinder probe did not record successful native publication")
        provenance["cylinder_native_result"] = {"artifact": artifact(args.cylinder_result), "record": result,
            "basis": "Saved publication outcome; this figure does not independently repeat native CAD gates"}
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": INK})
    height = 11.5 if len(items) == 5 else 6.6
    fig = plt.figure(figsize=(16.5, height), facecolor=BG)
    grid = fig.add_gridspec(2 if len(items) == 5 else 1, 6, left=.025, right=.975,
                            top=.89 if len(items) == 5 else .82, bottom=.17, hspace=.34, wspace=.12)
    for index, item in enumerate(items):
        group = items[:3] if index < 3 else items[3:]
        lo = np.min([entry["mesh"].points.min(axis=0) for entry in group], axis=0)
        hi = np.max([entry["mesh"].points.max(axis=0) for entry in group], axis=0)
        slot = grid[0, 2*index:2*index+2] if index < 3 else grid[1, 3*(index-3):3*(index-2)]
        draw(fig.add_subplot(slot, projection="3d"), item, (lo, hi), item["receipt"]["label"])
    fig.text(.04, .957 if len(items) == 5 else .94, "PAACS  /  Actual Tet4 cutaways", fontsize=27, weight="bold")
    fig.text(.04, .923 if len(items) == 5 else .885,
             "Stored whole elements · same cut and camera within each geometry · one mean-ratio scale", fontsize=13, color=SUB)
    color_ax = fig.add_axes([.29, .086, .42, .016 if len(items) == 5 else .024])
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, 1), cmap="RdYlBu"), cax=color_ax, orientation="horizontal")
    bar.set_label("Tet4 mean ratio  ·  regular tetrahedron = 1", color=SUB, fontsize=11)
    fig.text(.04, .036, "No CAD geometry was rendered or changed. Statistics cover every Tet4; cell IDs and source hashes are recorded in the figure provenance.", fontsize=10.3, color=SUB)
    if args.cylinder_before:
        fig.text(.04, .015, "Cylinder: same requested h = 0.5 and wall h = 0.1. Its instrumented correctness run is not a timing comparison.", fontsize=10.3, color=SUB)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, facecolor=BG, metadata={"Software": "PAACS measured archive renderer"})
    plt.close(fig)
    provenance["figure"] = artifact(args.output)
    provenance["render_seconds"] = time.perf_counter() - started
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"event": "complete", "output": str(args.output), "seconds": provenance["render_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
