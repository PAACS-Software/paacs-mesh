# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Render the actual hybrid duct interfaces, volume cells and axial flow."""

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from mesh_visual_data import Cell, EDGES, FACES, Mesh, clip_polygon, section, sha256


BG, INK, SUB = "#f1f5f8", "#173047", "#456176"
PRISM, TET, EXTERIOR = "#e9c183", "#7eb9c4", "#acced6"
INTERFACE, TEAL, RUST = "#244965", "#187681", "#bb542f"


def load_solution(run, row):
    path = run / Path(row["solution"]).name
    if sha256(path) != row["solution_sha256"]:
        raise ValueError(f"Solution checksum does not match the receipt: {path}")
    fields = ("points", "prisms", "tetrahedra", "velocity", "profile_coordinate_m",
              "profile_velocity", "profile_exact_velocity")
    with np.load(path, allow_pickle=False) as archive:
        data = {key: archive[key] for key in fields}
    points = data["points"]
    if points.ndim != 2 or points.shape[1] != 3 or data["velocity"].shape != (len(points),):
        raise ValueError("Expected 3D coordinates and one axial velocity per node")
    for key, count, width in [("prisms", row["prism6"], 6), ("tetrahedra", row["tet4"], 4)]:
        ids = data[key]
        if (ids.shape != (count, width) or not count or not np.issubdtype(ids.dtype, np.integer)
                or ids.min() < 0 or ids.max() >= len(points)):
            raise ValueError(f"Invalid {key} in {path}")
    profile = data["profile_coordinate_m"]
    if (profile.ndim != 1 or len(profile) < 2
            or any(data[key].shape != profile.shape for key in ("profile_velocity", "profile_exact_velocity"))
            or any(not np.isfinite(array).all() for array in data.values())):
        raise ValueError(f"Invalid coordinates, velocity or profile in {path}")
    return data


def face_inventory(blocks):
    exterior = []
    interface = np.empty((0, 3), dtype=int)
    for width in (3, 4):
        rows, kinds = [], []
        for kind, cells in blocks:
            local_faces = [face for face in FACES[kind] if len(face) == width]
            if local_faces:
                ids = cells[:, local_faces].reshape(-1, width)
                rows.append(ids)
                kinds.append(np.full(len(ids), kind, dtype=int))
        uses, kinds = np.concatenate(rows), np.concatenate(kinds)
        _, first, inverse, count = np.unique(np.sort(uses, axis=1), axis=0,
                                              return_index=True, return_inverse=True, return_counts=True)
        if np.any(count > 2):
            raise ValueError("Nonmanifold face incidence in the supplied geometry")
        exterior.extend(uses[first[count == 1]])
        if width == 3:
            prism_uses = np.bincount(inverse, weights=(kinds == 7))
            interface = uses[first[(count == 2) & (prism_uses == 1)]]
    if not len(interface):
        raise ValueError("No actual shared Prism6–Tet4 triangular facets")
    return exterior, interface


def cut_cells(points, blocks, axis, plane, keep_positive=False):
    cells = []
    for kind, ids in blocks:
        distances = points[ids, axis] - plane
        low, high = distances.min(axis=1), distances.max(axis=1)
        crosses = (low <= 0) & (high > 0 if keep_positive else high >= 0)
        cuts = ids[crosses]
        cells.extend(Cell(kind, 3, 1, tuple(cell)) for cell in cuts)
    polygons, owners, _ = section(Mesh(points, cells), range(len(cells)), None,
                                  axis=axis, value=plane)
    if not polygons:
        raise ValueError("The chosen section intersects no volume cells")
    return polygons, [cells[i] for i in owners]


def face_segments(points, faces, axis, plane):
    segments = []
    for ids in faces:
        p = points[ids]
        d = p[:, axis] - plane
        if min(d) > 0 or max(d) < 0:
            continue
        vertices = list(p[d == 0])
        for a, b in ((0, 1), (1, 2), (2, 0)):
            if d[a] * d[b] < 0:
                vertices.append(p[a] + (p[b]-p[a])*(-d[a])/(d[b]-d[a]))
        if len(vertices) == 2:
            segments.append(np.asarray(vertices))
        elif len(vertices) == 3:  # A genuinely coplanar shared triangle.
            segments.extend(np.asarray([vertices[a], vertices[b]]) for a, b in ((0, 1), (1, 2), (2, 0)))
    return segments


def color(kind):
    return PRISM if kind == 7 else TET


def cutaway(ax, points, blocks, exterior, interface, plane):
    polygons, owners = cut_cells(points, blocks, 0, plane, keep_positive=True)
    faces = [(poly, color(cell.kind)) for poly, cell in zip(polygons, owners)]
    view = np.array([-.87, .42, .37])
    view /= np.linalg.norm(view)
    right = np.cross([0., 1., 0.], view)
    right /= np.linalg.norm(right)
    screen = np.column_stack((right, np.cross(view, right)))
    for ids in exterior:
        poly = points[ids]
        if np.cross(poly[1]-poly[0], poly[2]-poly[0]) @ view <= 0:
            continue
        clipped = clip_polygon(poly, 0, plane, keep=1)
        if len(clipped) >= 3:
            faces.append((clipped, EXTERIOR))
    faces.sort(key=lambda item: item[0].mean(axis=0) @ view)
    projected = [poly @ screen for poly, _ in faces]
    ax.add_collection(PolyCollection(projected, facecolors=[c for _, c in faces],
                                     edgecolor=INTERFACE, linewidths=.19))
    ax.add_collection(LineCollection([edge @ screen for edge in face_segments(points, interface, 0, plane)],
                                     colors=INTERFACE, linewidths=.9))
    fit_mesh(ax, projected, .03)


def fit_mesh(ax, polygons, margin):
    locations = np.concatenate(polygons)
    lo, hi = locations.min(axis=0), locations.max(axis=0)
    pad = margin * min(hi-lo)
    ax.set(xlim=(lo[0]-pad, hi[0]+pad), ylim=(lo[1]-pad, hi[1]+pad), aspect="equal")
    ax.set_axis_off()
    ax.set_anchor("W")


def section_mesh(ax, points, cuts, owners, interface, plane):
    projected = [poly[:, :2] for poly in cuts]
    ax.add_collection(PolyCollection(projected, facecolors=[color(cell.kind) for cell in owners],
                                     edgecolor=INTERFACE, linewidths=.25))
    ax.add_collection(LineCollection([edge[:, :2] for edge in face_segments(points, interface, 2, plane)],
                                     colors=INTERFACE, linewidths=1.15))
    fit_mesh(ax, projected, .012)


def exposed_triangles(ax, points, interface):
    faces = points[interface]
    maximum_y = faces[:, :, 1].max()
    wall = faces[np.all(np.isclose(faces[:, :, 1], maximum_y, rtol=0, atol=1e-8), axis=1)]
    if not len(wall):
        raise ValueError("Expected a planar inner wall of shared triangular facets")
    projected = [poly[:, [2, 0]] for poly in wall]
    ax.add_collection(PolyCollection(projected, facecolor="#b9d8de",
                                     edgecolor=INTERFACE, linewidths=.34))
    fit_mesh(ax, projected, .045)
    return len(wall)


def plot_mesh(output, side, length, row, points, blocks, cuts, owners, interface, exterior, plane):
    fig = plt.figure(figsize=(13.2, 7.4), facecolor=BG)
    fig.text(.046, .930, "A prism shell joined to a tetrahedral core", fontsize=25, color=INK, weight="bold")
    fig.text(.047, .880, (f"{side:g} × {side:g} mm square duct · {length:g} mm length · "
                          f"{row['prism6']:,} Prism6 · {row['tet4']:,} Tet4"), fontsize=11, color=SUB)
    heading(fig, .047, .815, "Longitudinal cutaway at x = 0")
    heading(fig, .689, .815, "End-on mesh section")
    cutaway(fig.add_axes([.047, .392, .594, .397], facecolor=BG), points, blocks,
            exterior, interface, 0.)
    section_mesh(fig.add_axes([.688, .397, .257, .385], facecolor=BG), points, cuts, owners, interface, plane)

    heading(fig, .047, .337, "One extracted inner wall")
    count = exposed_triangles(fig.add_axes([.047, .106, .578, .209], facecolor=BG),
                              points, interface)
    fig.text(.047, .058, f"{count:,} triangular facets on this plane. Each facet has one Prism6 and one Tet4 owner.",
             fontsize=9, color=SUB)
    fig.legend(handles=[Patch(facecolor=PRISM, label="Wall-normal Prism6"),
                        Patch(facecolor=TET, label="Tet4 core"),
                        Line2D([], [], color=INTERFACE, label="Shared interface at the cut")],
               loc="upper left", bbox_to_anchor=(.685, .346), frameon=False,
               fontsize=9.5, labelcolor=SUB)
    fig.text(.694, .158, (f"{len(interface):,} shared triangles across four inner walls.\n"
                         f"First wall spacing {1000*row['first_wall_spacing_m']:.2f} mm."),
             fontsize=9, color=SUB, linespacing=1.55)
    fig.savefig(output, dpi=160, facecolor=BG)
    plt.close(fig)


def cut_vertex_values(points, velocity, cell, polygon, plane):
    xyz = points[list(cell.ids)]
    values = velocity[list(cell.ids)]
    distance = xyz[:, 2] - plane
    locations, samples = list(xyz[distance == 0]), list(values[distance == 0])
    for a, b in EDGES[cell.kind]:
        if distance[a] * distance[b] < 0:
            fraction = -distance[a] / (distance[b] - distance[a])
            locations.append(xyz[a] + fraction*(xyz[b]-xyz[a]))
            samples.append(values[a] + fraction*(values[b]-values[a]))
    delta = ((polygon[:, None] - np.asarray(locations))**2).sum(axis=2)
    closest = delta.argmin(axis=1)
    if np.any(delta[np.arange(len(polygon)), closest] > 1e-18):
        raise ValueError("A section vertex was not found on an original cell edge")
    return np.asarray(samples)[closest]


def flow_section(ax, points, velocity, cuts, owners, plane, norm):
    positions, values, triangles = [], [], []
    for poly, cell in zip(cuts, owners):
        first = len(positions)
        positions.extend(poly[:, :2])
        values.extend(cut_vertex_values(points, velocity, cell, poly, plane))
        triangles.extend((first, first+i, first+i+1) for i in range(1, len(poly)-1))
    xy = np.asarray(positions)
    return ax.tripcolor(xy[:, 0], xy[:, 1], np.asarray(triangles), np.asarray(values),
                        shading="gouraud", cmap="viridis", norm=norm)


def heading(fig, x, y, label):
    fig.text(x, y, label, fontsize=13, color=INK, weight="bold")


def chart(ax, xlabel, ylabel):
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.tick_params(colors=SUB, labelsize=8)
    ax.xaxis.label.set_color(SUB)
    ax.yaxis.label.set_color(SUB)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color("#bfccd5")
    ax.grid(axis="y", alpha=.22)
    ax.yaxis.set_major_locator(MaxNLocator(5))


def plot_flow(output, row, rows, pressure_reference, fluid, data, points, cuts, owners, plane):
    fig = plt.figure(figsize=(13.2, 6.25), facecolor=BG)
    fig.text(.047, .922, "Laminar flow in the hybrid duct", fontsize=25, color=INK, weight="bold")
    fig.text(.048, .862, (f"Reynolds number {fluid['reynolds']:g} · Bulk speed {fluid['bulk_speed_m_s']:g} m/s"
                          f" · Finest mesh {row['prism6'] + row['tet4']:,} mixed cells"), fontsize=11, color=SUB)
    heading(fig, .050, .778, "Axial velocity section")
    heading(fig, .408, .778, "Centerline velocity")
    heading(fig, .735, .778, "Pressure drop")

    ax = fig.add_axes([.057, .235, .254, .484], facecolor=BG)
    norm = Normalize(0., max(data["velocity"].max(), data["profile_exact_velocity"].max()))
    field = flow_section(ax, points, data["velocity"], cuts, owners, plane, norm)
    all_cuts = np.concatenate(cuts)
    lo, hi = all_cuts.min(axis=0), all_cuts.max(axis=0)
    ax.set(xlim=(lo[0], hi[0]), ylim=(lo[1], hi[1]), aspect="equal",
           xlabel="x (mm)", ylabel="y (mm)")
    ax.tick_params(colors=SUB, labelsize=8)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    colorbar = fig.colorbar(field, cax=fig.add_axes([.057, .146, .254, .018]), orientation="horizontal")
    colorbar.ax.tick_params(colors=SUB, labelsize=8)
    colorbar.set_label("Axial speed (m/s)", fontsize=8.5, color=SUB)
    colorbar.ax.xaxis.set_major_locator(MaxNLocator(5))

    ax = fig.add_axes([.409, .235, .254, .484], facecolor="#ffffff")
    x = data["profile_coordinate_m"]*1000
    ax.plot(x, data["profile_velocity"], color=TEAL, linewidth=2.45, label="Finest hybrid")
    ax.plot(x, data["profile_exact_velocity"], color=RUST, linewidth=1.5,
            linestyle="--", label="Square-duct series")
    ax.set(xlim=(x.min(), x.max()), ylim=(0., None))
    chart(ax, "x at y = 0, z = L/2 (mm)", "Axial speed (m/s)")
    ax.legend(frameon=False, fontsize=8, loc="lower center")

    ax = fig.add_axes([.737, .235, .217, .484], facecolor="#ffffff")
    spacing = [row["first_wall_spacing_m"]*1000 for row in rows]
    pressure = [row["pressure_drop_pa"] for row in rows]
    ax.plot(spacing, pressure, "o-", color=TEAL, linewidth=1.8, label="Hybrid solve")
    ax.axhline(pressure_reference, linestyle="--", color=RUST, linewidth=1.25, label="Square-duct series")
    for sample, x, y in zip(rows, spacing, pressure):
        ax.annotate(f"{sample['prism6'] + sample['tet4']:,}", (x, y), xytext=(0, 7),
                    textcoords="offset points", ha="center", color=SUB, fontsize=7.7)
    ax.margins(x=.14, y=.22)
    ax.invert_xaxis()
    chart(ax, "First wall spacing (mm) · finer →", "Pressure drop (Pa)")
    ax.legend(frameon=False, fontsize=7.6, loc="best")
    fig.text(.048, .042, (f"Section at z = {plane:.2f} mm. Shading interpolates solved edge values. "
                          "Pressure labels count all volume cells."),
             fontsize=9, color=SUB)
    fig.savefig(output, dpi=160, facecolor=BG)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="Directory holding summary.json and the solved NPZ files")
    parser.add_argument("output", type=Path, help="Directory for hybrid-mesh.png and flow-convergence.png")
    args = parser.parse_args()
    summary = json.loads((args.run / "summary.json").read_text(encoding="utf-8-sig"))
    rows = sorted(summary["meshes"], key=lambda item: item["first_wall_spacing_m"], reverse=True)
    final = rows[-1]
    data = load_solution(args.run, final)
    points = data["points"]*1000
    blocks = [(7, data["prisms"]), (5, data["tetrahedra"])]
    exterior, interface = face_inventory(blocks)
    geometry = summary["geometry"]
    side, length = geometry["side_m"]*1000, geometry["length_m"]*1000
    plane = points[:, 2].min() + .503*length
    cuts, owners = cut_cells(points, blocks, 2, plane)
    args.output.mkdir(parents=True, exist_ok=True)
    plot_mesh(args.output / "hybrid-mesh.png", side, length, final, points, blocks,
              cuts, owners, interface, exterior, plane)
    plot_flow(args.output / "flow-convergence.png", final, rows,
              summary["reference"]["pressure_drop_pa"], summary["fluid"], data, points, cuts, owners, plane)


if __name__ == "__main__":
    main()
