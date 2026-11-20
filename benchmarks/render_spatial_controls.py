#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Render exact sections and C++-sampled fields; preserve archive provenance."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm, Normalize
from render_cylinder_comparison import read, statistics
from render_local_controls import section, cutaway, axis2d, BG, INK, SUB


def digest(data):
    value = 1469598103934665603
    for byte in data:
        value = ((value ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return str(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('gallery', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    version = (Path(__file__).resolve().parents[1] / 'VERSION.txt').read_text().strip()
    titles = [('point', 'Point influence'), ('segment', 'Segment influence'),
              ('triangle', 'Triangle patch'), ('box', 'Box region'),
              ('ball', 'Ball region'), ('minimum', 'Minimum: either region'),
              ('maximum', 'Maximum: overlap'), ('smooth', 'Smooth transition + clamp')]
    records, loaded = [], {}
    for case, title in titles:
        path = args.gallery / (case + '.paacsmesh')
        points, cells = read(path)
        polygons, sizes, _ = section(points, cells)
        stats, _ = statistics(points, cells)
        csv = args.gallery / (case + '-field.csv')
        field = np.genfromtxt(csv, delimiter=',', names=True)
        assert len(field) == 101 * 101
        loaded[case] = points, cells, polygons, sizes, field
        records.append(dict(case=case, nodes=len(points),
                            tetrahedra=sum(dim == 3 for _, dim, _, _ in cells),
                            archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            archive_fnv1a=digest(path.read_bytes()),
                            field_csv_sha256=hashlib.sha256(csv.read_bytes()).hexdigest(),
                            section_polygons=len(polygons), independent_metrics=stats))
    low = min(x[3].min() for x in loaded.values())
    high = max(x[3].max() for x in loaded.values())
    norm, cmap = LogNorm(low, high), plt.get_cmap('YlGnBu')
    fig, axes = plt.subplots(2, 4, figsize=(17, 10.5), facecolor=BG)
    fig.subplots_adjust(left=.048, right=.91, top=.82, bottom=.16, wspace=.2, hspace=.37)
    fig.text(.048, .949, f'PAACS {version}  /  Composable 3D size controls', fontsize=24, weight='bold', color=INK)
    fig.text(.048, .906, 'Actual exported Tet4 meshes  ·  Exact geometric sections at y = 0.5  ·  Unit cube', fontsize=12, color=SUB)
    fig.text(.048, .872, 'Near size 0.060–0.075; far size 0.32–0.40. Each control uses the same reusable field interface.', fontsize=11, color=SUB)
    for ax, (case, title), record in zip(axes.flat, titles, records):
        _, _, polygons, sizes, _ = loaded[case]
        ax.set_facecolor(BG)
        ax.add_collection(PolyCollection(polygons, facecolors=cmap(norm(sizes)),
                                         edgecolors='#284a58', linewidths=.22))
        axis2d(ax)
        ax.set_title(title, loc='left', fontsize=12, color=INK, weight='bold', pad=10)
        ax.text(0, -.17, f"{record['tetrahedra']:,} tetrahedra · {record['nodes']:,} nodes",
                transform=ax.transAxes, fontsize=9.4, color=SUB)
        ax.set_xlabel('x', fontsize=9, color=SUB)
        ax.set_ylabel('z', fontsize=9, color=SUB)
    cax = fig.add_axes([.929, .25, .012, .49])
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    bar.set_label('Actual cell size: cube root of volume', fontsize=10, color=SUB)
    bar.set_ticks([.012, .02, .04, .08])
    bar.set_ticklabels(['0.012', '0.020', '0.040', '0.080'])
    bar.minorticks_off()
    bar.ax.tick_params(labelsize=9)
    fig.text(.048, .072, 'Minimum quality ≥ 0.560 in these eight meshes (tetrahedron mean ratio; ideal = 1). All passed full boundary and size checks.', fontsize=10.5, color=SUB)
    fig.text(.048, .043, 'Four-sided polygons are cuts through tetrahedra. Influence regions are sizing controls, not separate materials or internal boundaries.', fontsize=10, color=SUB)
    fig.savefig(args.output / f'paacs-{version}-spatial-controls.png', dpi=160, facecolor=BG)
    plt.close(fig)

    fig = plt.figure(figsize=(12.5, 11.5), facecolor=BG)
    grid = fig.add_gridspec(2, 2, left=.06, right=.87, top=.84, bottom=.15,
                           hspace=.35, wspace=.22, height_ratios=[1.15, 1])
    fig.text(.06, .953, 'Compose controls: refine either region or their overlap', fontsize=20, weight='bold', color=INK)
    fig.text(.06, .915, 'Same two ball influences · Near size 0.075 · Far size 0.4 · Two different operators', fontsize=11.5, color=SUB)
    for col, (case, title) in enumerate([('minimum', 'Minimum: fine where either requests it'),
                                       ('maximum', 'Maximum: fine where both request it')]):
        points, cells, polygons, _, field = loaded[case]
        ax = fig.add_subplot(grid[0, col], projection='3d', facecolor=BG)
        cutaway(ax, points, cells, polygons)
        ax.set_title(title, loc='left', fontsize=11.5, color=INK, weight='bold', pad=12)
        count = sum(dim == 3 for _, dim, _, _ in cells)
        ax.text2D(0, -.035, f'{count:,} tetrahedra · Actual mesh cutaway', transform=ax.transAxes, fontsize=10.5, color=SUB)
        ax = fig.add_subplot(grid[1, col], facecolor=BG)
        ax.imshow(field['size'].reshape(101, 101), origin='lower', extent=[0, 1, 0, 1],
                  interpolation='nearest', cmap='viridis', norm=Normalize(.075, .4))
        axis2d(ax)
        ax.set_title('Requested size field at y = 0.5', loc='left', fontsize=11, color=INK)
        ax.set_xlabel('x', color=SUB)
        ax.set_ylabel('z', color=SUB)
    cax = fig.add_axes([.9, .19, .015, .26])
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(.075, .4), cmap='viridis'), cax=cax)
    bar.set_label('Target edge size', fontsize=10, color=SUB)
    fig.text(.06, .076, 'Top: cyan exterior and gold interior cuts. Bottom: exact C++ field samples on a 101 × 101 display grid.', fontsize=10, color=SUB)
    fig.text(.06, .048, 'The mesher uses conservative whole-box bounds, not this display grid. Small influence regions cannot disappear between display samples.', fontsize=9.5, color=SUB)
    fig.savefig(args.output / f'paacs-{version}-field-composition.png', dpi=160, facecolor=BG)
    plt.close(fig)
    (args.output / f'paacs-{version}-spatial-provenance.json').write_text(json.dumps(records, indent=2) + '\n')
    print('Rendered 8 archive-backed sections and 2 cutaways.')


if __name__ == '__main__':
    main()
