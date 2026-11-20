# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Plot measured in-plane edge sizes directly from the archived plate caps."""

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from mesh_visual_data import read_archive, sha256


def cap(run, row):
    path = run/row["variant"]/"mesh.paacsmesh"
    if sha256(path) != row["archive_sha256"]:
        raise ValueError("Plate mesh differs from the recorded archive")
    mesh = read_archive(path)
    tag = row["original_face_tags"][1]
    triangles = mesh.points[[cell.ids for cell in mesh.cells if cell.dim == 2 and cell.tag == tag]]
    if not len(triangles) or not np.all(triangles[:, :, 2] == triangles[0, 0, 2]):
        raise ValueError("Expected the original planar upper cap")
    return triangles[:, :, :2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="Initial study directory with all four meshes")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    receipt = json.loads((args.run/"metrics.json").read_text(encoding="utf-8"))
    box = receipt["sizing"]["box"]
    x, y, _ = box["lower_mm"]
    X, Y, _ = box["upper_mm"]
    names = {"global": "Global size", "bore": "Bore wall", "box": "Refinement box", "combined": "Bore wall + box"}
    background = "#f1f5f8"
    ink, muted = "#1a3349", "#476175"
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 7.9), facecolor=background)
    fig.subplots_adjust(left=.08, right=.90, bottom=.14, top=.82, wspace=.22, hspace=.44)
    fig.text(.08, .94, "STEP plate with combined sizing controls", color=ink, fontsize=24, weight="bold")
    fig.text(.08, .884, "Actual upper-cap triangles · Bore face 7 · Dashed rectangle marks the box footprint",
             color=muted, fontsize=10.5)
    norm = LogNorm(.22, 3.5)
    for ax, row in zip(axes.ravel(), receipt["meshes"], strict=True):
        triangles = cap(args.run, row)
        edges = np.roll(triangles, -1, axis=1)-triangles
        colors = np.median(np.linalg.norm(edges, axis=2), axis=1)
        collection = PolyCollection(triangles, cmap="viridis_r", norm=norm, edgecolors="#1d4056",
                                    linewidths=.22)
        collection.set_array(colors)
        ax.add_collection(collection)
        ax.add_patch(Rectangle((x, y), X-x, Y-y, fill=False, edgecolor="#d34c21",
                               linestyle="--", linewidth=1.25))
        ax.set(xlim=(-13, 21), ylim=(-9, 9), aspect="equal", facecolor="white")
        ax.set_xticks([-10, 0, 10, 20])
        ax.set_yticks([-8, 0, 8])
        ax.tick_params(colors=muted, labelsize=8)
        ax.set_title(f"{names[row['variant']]} · {row['tet4']:,} Tet4", loc="left", color=ink,
                     fontsize=11, weight="bold")
        ax.set_xlabel("x (mm)", color=muted, fontsize=9)
        ax.set_ylabel("y (mm)", color=muted, fontsize=9)
    bar = fig.colorbar(collection, cax=fig.add_axes([.91, .25, .016, .43]),
                       ticks=[.25, .5, 1, 2, 3.5], extend="both")
    bar.ax.set_yticklabels(["0.25", "0.50", "1.00", "2.00", "3.50"])
    bar.set_label("Median triangle edge (mm)", color=muted, fontsize=9)
    bar.ax.tick_params(colors=muted, labelsize=8)
    combined = next(row for row in receipt["meshes"] if row["variant"] == "combined")
    fig.text(.08, .055, f"Combined mesh · Median edge {combined['box_interior_in_plane_mm']['median']:.2f} mm "
             f"inside the box away from the bore and {combined['far_in_plane_mm']['median']:.2f} mm in the far field.",
             color=muted, fontsize=9)
    fig.savefig(args.output, dpi=160, facecolor=background)
    plt.close(fig)


if __name__ == "__main__":
    main()
