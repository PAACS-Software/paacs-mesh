# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Plot the solved cantilever and displacement convergence."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import numpy as np

from solve import read_archive, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = json.loads((args.run / "convergence.json").read_text(encoding="utf-8"))
    last = rows[-1]
    folder = args.run / f"level-{last['level']}"
    archive = folder / "mesh.paacsmesh"
    if sha256(archive) != last["archive_sha256"]:
        raise ValueError("Convergence receipt and mesh archive hashes differ")
    mesh = read_archive(archive)
    tris = np.array([cell.ids for cell in mesh.cells if cell.kind == 2])
    side = tris[np.all(np.isclose(mesh.points[tris, 1], 0, atol=1e-9, rtol=0), axis=1)]
    with np.load(folder / "displacement.npz") as data:
        displacement = data["displacement_mm"]
    factor = 60
    warped = mesh.points + factor*displacement
    face_values = -displacement[side, 2].mean(axis=1)
    norm = Normalize(0, max(-displacement[:, 2]))
    palette = plt.get_cmap("viridis")
    bg, ink, sub = "#f3f6fa", "#183648", "#4d6477"
    fig = plt.figure(figsize=(12.1, 7.3), facecolor=bg)
    fig.text(.059, .932, "STEP cantilever  |  PAACS Tet4 elasticity", fontsize=25, color=ink, weight="bold")
    fig.text(.060, .887, "100 × 10 × 10 mm · Steel, E = 210 GPa, ν = 0.30 · Uniform 100 N transverse end load",
             fontsize=11, color=sub)

    ax = fig.add_axes([.075,.521,.75,.286], facecolor=bg)
    ax.add_collection(PolyCollection(warped[side][:, :, (0, 2)], facecolors=palette(norm(face_values)),
                                     edgecolors="#294654", linewidths=.12))
    ax.plot([0, 100, 100, 0, 0], [0, 0, 10, 10, 0], color=sub, linewidth=.8, linestyle="--", alpha=.65)
    ax.set(xlim=(-4,104), ylim=(-14,13))
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    cax = fig.add_axes([.846,.546,.016,.231])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=palette), cax=cax)
    cb.set_label("Downward displacement (mm)", color=sub)
    cb.ax.tick_params(colors=sub, labelsize=8)
    fig.text(.059,.474,f"Finest mesh · {last['tet4']:,} tetrahedra · Side view at {factor}× deformation",
             color=sub, fontsize=10)

    chart = fig.add_axes([.09,.106,.80,.301], facecolor=bg)
    edges = [row["edge_p50_mm"] for row in rows]
    deflections = [row["tip_deflection_mm"] for row in rows]
    chart.plot(edges, deflections, "o-", color="#087c89", linewidth=2.2, label="3D Tet4 finite elements")
    chart.axhline(last["euler_bernoulli_mm"], color="#a55d22", linestyle="--", linewidth=1.4,
                  label="Euler–Bernoulli")
    chart.axhline(last["timoshenko_mm"], color="#786195", linestyle=":", linewidth=1.4,
                  label="Timoshenko")
    for i, row in enumerate(rows):
        offset = 7 if i == 0 else -12
        chart.annotate(f"{row['tet4']:,}", (edges[i], deflections[i]), xytext=(0,offset),
                       textcoords="offset points", ha="center", color=sub, fontsize=8)
    chart.set_xlabel("Achieved median tetrahedron edge (mm)", color=sub)
    chart.set_ylabel("Loaded-face mean deflection (mm)", color=sub)
    chart.invert_xaxis()
    chart.tick_params(colors=sub)
    chart.grid(alpha=.25)
    chart.spines[["top","right"]].set_visible(False)
    chart.legend(frameon=False, loc="lower right", fontsize=8)
    fig.text(.06,.023,"Fixed face at x = 0. Labels along the convergence curve are tetrahedron counts.",
             color=sub, fontsize=9)
    fig.savefig(args.output, dpi=160, facecolor=bg)
    plt.close(fig)


if __name__ == "__main__":
    main()
