# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Plot actual Tet4 sections and stress convergence for the STEP plate."""

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from mesh_visual_data import Cell, Mesh, section, sha256


BACKGROUND = "#f1f5f8"
PANEL = "#ffffff"
INK = "#173047"
SUBTLE = "#425e72"
TEAL = "#167a83"
RUST = "#bb542e"


def load_section(run, row, plane):
    path = run / Path(row["solution"]).name
    if sha256(path) != row["solution_sha256"]:
        raise ValueError(f"Solution hash in summary does not match {path}")
    with np.load(path, allow_pickle=False) as data:
        points = data["points"]
        tets = data["tetrahedra"]
        stress = data["hoop_stress"]
    if (points.ndim != 2 or points.shape[1] != 3 or tets.ndim != 2 or
            tets.shape[1] != 4 or len(tets) != row["tet4"] or len(tets) == 0 or
            stress.shape != (len(tets),)):
        raise ValueError(f"Expected Tet4 coordinates and one hoop stress per cell in {path}")
    if (not np.isfinite(points).all() or not np.isfinite(stress).all() or
            not np.issubdtype(tets.dtype, np.integer) or tets.min() < 0 or tets.max() >= len(points)):
        raise ValueError(f"Invalid Tet4 plot data in {path}")
    if not points[:, 2].min() < plane < points[:, 2].max():
        raise ValueError(f"Section plane does not cut the plate in {path}")
    cells = [Cell(5, 3, 0, tuple(node_ids)) for node_ids in tets]
    polygons, owners, _ = section(Mesh(points, cells), range(len(cells)), None, axis=2, value=plane)
    if not len(polygons):
        raise ValueError(f"Section contains no tetrahedral intersections in {path}")
    return [poly[:, :2] for poly in polygons], stress[owners], len(tets)


def add_section(ax, polygons, values=None, norm=None, cmap=None):
    if values is None:
        collection = PolyCollection(polygons, facecolors="#c4e1e4", edgecolors="#286575", linewidths=.25)
    else:
        collection = PolyCollection(polygons, cmap=cmap, norm=norm, edgecolors="#233348",
                                    linewidths=.14)
        collection.set_array(values)
    ax.add_collection(collection)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()


def plate_panel(ax, case, title, bounds, closeup):
    polygons, _, count = case
    add_section(ax, polygons)
    lo, hi = bounds
    span = hi - lo
    ax.set(xlim=(lo[0] - .01*span[0], hi[0] + .01*span[0]),
           ylim=(lo[1] - .02*span[1], hi[1] + .02*span[1]))
    ax.set_title(f"{title} · {count:,} Tet4", loc="left", fontsize=13, pad=8,
                 color=INK, fontweight="bold")
    ax.add_patch(Rectangle((-closeup, -closeup), 2*closeup, 2*closeup,
                           fill=False, edgecolor=RUST, linewidth=1, linestyle="--"))


def stress_panel(ax, case, title, norm, cmap, closeup):
    polygons, stresses, _ = case
    add_section(ax, polygons, stresses, norm, cmap)
    ax.set(xlim=(-closeup, closeup), ylim=(-closeup, closeup))
    ax.set_title(title, loc="left", fontsize=11, pad=7, color=INK, fontweight="bold")


def family(meshes, label):
    rows = sorted((row for row in meshes if row["label"].startswith(f"{label}-")),
                  key=lambda row: row["hole_edge_p50_mm"], reverse=True)
    if not rows:
        raise ValueError(f"No {label} mesh in the solve summary")
    return rows


def style_chart(ax, title, xlabel, ylabel):
    ax.set_title(title, loc="left", fontsize=12, pad=6, color=INK, fontweight="bold")
    ax.set_xlabel(xlabel, color=SUBTLE, fontsize=9)
    ax.set_ylabel(ylabel, color=SUBTLE, fontsize=9)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.tick_params(colors=SUBTLE, labelsize=8)
    ax.grid(axis="y", alpha=.28)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color("#b8cad6")


def draw_convergence(ax, efficiency, uniform, local):
    # All solves use the same physical quadrature positions near the hole.
    for rows, color, marker in [(uniform, TEAL, "o"), (local, RUST, "s")]:
        ax.plot([row["hole_edge_p50_mm"] for row in rows],
                [row["physical_probe_strip"]["fem_sigma_tt_mpa"] for row in rows], marker=marker,
                color=color, linewidth=1.9, markersize=4.7)
        by_count = sorted(rows, key=lambda row: row["tet4"])
        efficiency.plot([row["tet4"] for row in by_count],
                        [100 * row["fixed_annulus"]["relative_kirsch_stress_l2"] for row in by_count],
                        marker=marker, color=color, linewidth=1.9, markersize=4.7)
    ax.axhline(uniform[0]["physical_probe_strip"]["kirsch_sigma_tt_mpa"], color="#586d87",
                linewidth=1.15, linestyle="--")
    style_chart(ax, "Mean hoop stress near the rim", "Median hole-wall edge (mm) · finer →", "Stress (MPa)")
    ax.invert_xaxis()
    style_chart(efficiency, "Stress error versus Tet4 count", "Total Tet4 cells", "Near-hole error (%)")
    efficiency.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    efficiency.xaxis.set_major_formatter(FuncFormatter(lambda n, _: f"{n:,.0f}"))


def draw(run, output):
    summary = json.loads((run / "summary.json").read_text(encoding="utf-8-sig"))
    geometry = summary["geometry"]
    bounds = np.asarray(geometry["bounds_mm"])
    lo, hi = bounds
    radius = geometry["radius_mm"]
    closeup = min(2.6*radius, .95*min(hi[:2]))
    plane = (lo[2] + hi[2])*.5 + .035*(hi[2] - lo[2])
    uniform = family(summary["meshes"], "uniform")
    local = family(summary["meshes"], "local")
    unrefined = load_section(run, uniform[-1], plane)
    refined = load_section(run, local[-1], plane)
    cmap = plt.get_cmap("RdYlBu_r")
    cropped_values = []
    for polygons, stresses, _ in (unrefined, refined):
        near = np.array([np.all(poly.min(axis=0) < closeup) and np.all(poly.max(axis=0) > -closeup)
                         for poly in polygons])
        cropped_values.extend(stresses[near])
    low, high = min(cropped_values), max(cropped_values)
    norm = Normalize(low, high)

    fig = plt.figure(figsize=(13, 9.2), facecolor=BACKGROUND)
    fig.text(.048, .948, "STEP plate with a hole", fontsize=27, color=INK, fontweight="bold")
    dimensions = " × ".join(f"{value:g}" for value in hi - lo)
    fig.text(.049, .907, f"{dimensions} mm plate · Ø{2*radius:g} mm hole · Three-dimensional Tet4 elasticity · Tension along x",
             fontsize=11, color=SUBTLE)
    fig.text(.049, .867, "MESH SECTIONS", fontsize=9, color=TEAL, fontweight="bold")
    plate_panel(fig.add_axes([.047, .524, .435, .306]), unrefined, "Finest uniform", bounds, closeup)
    plate_panel(fig.add_axes([.518, .524, .435, .306]), refined, "Finest local", bounds, closeup)

    fig.text(.049, .478, "HOOP STRESS NEAR THE HOLE", fontsize=9, color=TEAL, fontweight="bold")
    stress_panel(fig.add_axes([.047, .167, .205, .270]), unrefined, "Finest uniform", norm, cmap, closeup)
    stress_panel(fig.add_axes([.267, .167, .205, .270]), refined, "Finest local", norm, cmap, closeup)
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap),
                      cax=fig.add_axes([.049, .134, .421, .016]), orientation="horizontal")
    cb.set_label("Element hoop stress (MPa) · common scale", color=SUBTLE, fontsize=9)
    cb.ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    cb.ax.tick_params(colors=SUBTLE, labelsize=8)
    cb.outline.set_visible(False)

    ax = fig.add_axes([.565, .352, .387, .124], facecolor=PANEL)
    efficiency = fig.add_axes([.565, .149, .387, .124], facecolor=PANEL)
    draw_convergence(ax, efficiency, uniform, local)
    fig.legend(handles=[
        Line2D([], [], color=TEAL, marker="o", label="Uniform"),
        Line2D([], [], color=RUST, marker="s", label="Local"),
        Line2D([], [], color="#586d87", linestyle="--", label="Kirsch reference"),
    ], loc="lower left", bbox_to_anchor=(.555, .496), ncol=3, frameon=False, fontsize=8.6,
       labelcolor=SUBTLE, markerscale=.75)
    fig.text(.565, .073, "Stress uses a fixed near-rim band. Error uses the same annulus on all meshes.",
             fontsize=8.1, color=SUBTLE)
    fig.text(.049, .046, f"Sections cut the original tetrahedra at z = {plane:.2f} mm. "
             "Colors show hoop stress at each Tet4 centroid.", fontsize=9, color=SUBTLE)
    fig.savefig(output, dpi=160, facecolor=BACKGROUND)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="Directory with the solve summary and per-case NPZ files")
    parser.add_argument("output", type=Path, help="Output figure path")
    args = parser.parse_args()
    draw(args.run, args.output)


if __name__ == "__main__":
    main()
