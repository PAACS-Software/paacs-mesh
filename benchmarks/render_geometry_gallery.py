#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Render actual CAD-CLI output. Optional dependency: matplotlib."""
import argparse
import json
import struct
import subprocess
import tempfile
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

parser = argparse.ArgumentParser()
parser.add_argument("cad_cli", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
fixtures = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "cad"
cases = [("beam_i", "I-beam / Prism6", []),
         ("perforated_plate", "Nine through-holes / Prism6", []),
         ("oblique_l", "Oblique L-profile / 12 graded layers", ["--layers", "12", "--layer-ratio", "4"]),
         ("box", "Mapped Hex8 / 12 graded layers", ["--layers", "12", "--layer-ratio", ".25"])]
fig = plt.figure(figsize=(12, 9), facecolor="#f6f8fa")
fig.suptitle("PAACS 1.0.0-rc.3  |  Engineering profiles and grading", x=.055, y=.975,
             ha="left", fontsize=19, fontweight="bold", color="#152d40")
fig.text(.055, .932, "Actual STEP imports and generated boundary meshes. C++26 / GCC 16.",
         fontsize=11, color="#536477")
with tempfile.TemporaryDirectory() as directory:
    for index, (name, title, extra) in enumerate(cases):
        archive = Path(directory) / f"{name}.paacsmesh"
        result = subprocess.run([str(args.cad_cli.resolve()), str(fixtures / f"{name}.step"),
                                 ".3", str(archive), *extra], capture_output=True, text=True, check=True)
        row = next(json.loads(line) for line in result.stdout.splitlines() if line.startswith("{"))
        data = archive.read_bytes()
        assert data[:8] == b"PAACSM8\0"
        nodes, cells, _ = struct.unpack_from("<QQQ", data, 16)
        points = [struct.unpack_from("<ddd", data, 40 + 32*i) for i in range(nodes)]
        offset = 40 + 32 * nodes
        polygons, colors = [], []
        palette = ["#75c8d3", "#368caa", "#e8ba74", "#b5dce1", "#5eb2c1", "#dae5ed"]
        for _ in range(cells):
            kind, dim, tag = struct.unpack_from("<B3xii", data, offset)
            count = [2, 3, 3, 6, 4, 4, 10, 6, 8, 5][kind]
            conn = struct.unpack_from("<" + "q" * count, data, offset + 12)
            offset += 12 + 8 * count
            if kind in (2, 4):
                polygons.append([points[i] for i in conn])
                colors.append(palette[tag % len(palette)])
        assert offset == len(data)
        ax = fig.add_subplot(2, 2, index + 1, projection="3d", facecolor="#f6f8fa")
        ax.add_collection3d(Poly3DCollection(polygons, facecolors=colors,
                                             edgecolors="#264553", linewidths=.32,
                                             alpha=1, zsort="average"))
        ranges = []
        for axis in range(3):
            lo, hi = min(p[axis] for p in points), max(p[axis] for p in points)
            ranges.append(hi-lo)
            [ax.set_xlim, ax.set_ylim, ax.set_zlim][axis](lo, hi)
        ax.set_box_aspect(ranges)
        ax.view_init(elev=32, azim=-58)
        ax.set_axis_off()
        ax.set_title(title, loc="left", fontsize=13, color="#152d40", pad=0)
        ax.text2D(0, -.035, f'{nodes:,} nodes  |  {cells:,} total cells  |  boundary audit passed',
                  transform=ax.transAxes, fontsize=10, color="#536477")
fig.text(.055, .025, "Controlled planar fixtures; curved CAD and general automatic decomposition remain outside this qualification.",
         fontsize=9, color="#536477")
fig.subplots_adjust(left=.055, right=.96, bottom=.095, top=.89, wspace=.18, hspace=.25)
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, dpi=170, facecolor=fig.get_facecolor())
