# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Render Prism6 cutaways and fully developed square-duct flow results."""

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib.tri import Triangulation
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from mesh_visual_data import Cell, FACES, Mesh, clip_polygon, read_archive, section, sha256


BG = "#f1f5f8"
INK = "#173047"
SUB = "#456176"
TEAL = "#177880"
RUST = "#ba552e"
GOLD = "#e8bf79"
OUTSIDE = "#8dc9cb"


def load_solution(run, row):
    path = run / Path(row["solution"]).name
    if sha256(path) != row["solution_sha256"]:
        raise ValueError(f"Solution checksum does not match its solve receipt: {path}")
    keys = ("points", "prisms", "section_points", "section_triangles", "section_velocity",
            "profile_coordinate_m", "profile_velocity", "profile_exact_velocity")
    with np.load(path, allow_pickle=False) as arrays:
        data = {key: arrays[key] for key in keys}
    points, prisms = data["points"], data["prisms"]
    xy, triangles = data["section_points"], data["section_triangles"]
    if points.ndim != 2 or points.shape[1] != 3 or xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError(f"Unexpected coordinate arrays in {path}")
    for name, ids, count, width, nodes in [
        ("prisms", prisms, row["prism6"], 6, len(points)),
        ("section triangles", triangles, row["section_tri3"], 3, len(xy)),
    ]:
        if (not len(ids) or ids.shape != (count, width) or
                not np.issubdtype(ids.dtype, np.integer) or ids.min() < 0 or ids.max() >= nodes):
            raise ValueError(f"Invalid {name} in {path}")
    if data["section_velocity"].shape != (len(xy),):
        raise ValueError(f"Expected one axial velocity per cross-section node in {path}")
    profile = data["profile_coordinate_m"]
    if (profile.ndim != 1 or len(profile) < 2 or
            any(data[key].shape != profile.shape for key in ("profile_velocity", "profile_exact_velocity"))):
        raise ValueError(f"Invalid line profile in {path}")
    if any(not np.isfinite(array).all() for array in data.values()):
        raise ValueError(f"Nonfinite coordinates or flow values in {path}")
    return data


def exterior_faces(prisms):
    # Group by sorted identity within each arity and retain the original
    # outward order. Shared interior faces occur twice.
    for width in (3, 4):
        local_ids = [face for face in FACES[7] if len(face) == width]
        uses = prisms[:, local_ids].reshape(-1, width)
        _, first, count = np.unique(np.sort(uses, axis=1), axis=0,
                                    return_index=True, return_counts=True)
        yield from uses[first[count == 1]]


def prism_sections(points, prisms, axis, plane):
    extent = points[prisms, axis]
    intersecting = prisms[(extent.min(axis=1) <= plane) & (extent.max(axis=1) >= plane)]
    cells = [Cell(7, 3, 1, tuple(ids)) for ids in intersecting]
    cuts, _, _ = section(Mesh(points, cells), range(len(cells)), None, axis=axis, value=plane)
    if not cuts:
        raise ValueError("The section does not intersect any Prism6")
    return cuts


def cutaway(ax, points, prisms, color_section=lambda _: GOLD):
    plane = (points[:, 0].min() + points[:, 0].max()) * .5
    cuts = prism_sections(points, prisms, 0, plane)

    # Orthographic camera facing the open cut, top wall and outlet.
    view = np.array([-.87, .42, .37])
    view /= np.linalg.norm(view)
    right = np.cross([0., 1., 0.], view)
    right /= np.linalg.norm(right)
    up = np.cross(view, right)
    screen = np.column_stack((right, up))
    faces = [(poly, color_section(poly)) for poly in cuts]
    for ids in exterior_faces(prisms):
        poly = points[ids]
        normal = np.cross(poly[1] - poly[0], poly[2] - poly[0])
        if normal @ view <= 0:
            continue
        visible = clip_polygon(poly, 0, plane, keep=1)
        if len(visible) >= 3:
            faces.append((visible, OUTSIDE))

    faces.sort(key=lambda pair: pair[0].mean(axis=0) @ view)
    projected = [poly @ screen for poly, _ in faces]
    ax.add_collection(PolyCollection(projected, facecolors=[color for _, color in faces],
                                     edgecolors="#355969", linewidths=.24))
    bounds = np.concatenate(projected)
    low, high = bounds.min(axis=0), bounds.max(axis=0)
    pad = .03 * (high[0] - low[0])
    ax.set(xlim=(low[0] - pad, high[0] + pad), ylim=(low[1] - pad, high[1] + pad),
           aspect="equal")
    ax.set_axis_off()
    ax.set_anchor("W")


def cross_section(ax, data, norm):
    xy = data["section_points"] * 1000
    triangles = Triangulation(xy[:, 0], xy[:, 1], data["section_triangles"])
    field = ax.tripcolor(triangles, data["section_velocity"], shading="gouraud",
                         cmap="viridis", norm=norm)
    ax.triplot(triangles, color="#172f40", linewidth=.30, alpha=.62)
    low, high = xy.min(axis=0), xy.max(axis=0)
    ax.set(xlim=(low[0], high[0]), ylim=(low[1], high[1]), aspect="equal",
           xlabel="x (mm)", ylabel="y (mm)")
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.tick_params(labelsize=8, colors=SUB)
    return field


def chart(ax, title, xlabel, ylabel):
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.set_title(title, loc="left", color=INK, fontweight="bold", fontsize=13, pad=9)
    ax.tick_params(labelsize=8.5, colors=SUB)
    ax.xaxis.label.set_color(SUB)
    ax.yaxis.label.set_color(SUB)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color("#bdcbd3")
    ax.grid(axis="y", alpha=.25)


def velocity_profile(ax, data):
    coordinate = data["profile_coordinate_m"] * 1000
    ax.plot(coordinate, data["profile_velocity"], color=TEAL, linewidth=2.6,
            label="Finest PAACS section")
    ax.plot(coordinate, data["profile_exact_velocity"], color=RUST, linewidth=1.45,
            dashes=(3, 3), label="Square-duct series")
    chart(ax, "Velocity through the centerline", "x at y = 0 (mm)", "Axial speed (m/s)")
    ax.set_xlim(coordinate.min(), coordinate.max())
    ax.set_ylim(bottom=0.)
    ax.legend(frameon=False, fontsize=8.5)


def pressure_drop(ax, rows, reference):
    wall_steps = [row["first_wall_spacing_m"] * 1000 for row in rows]
    drops = [row["pressure_drop_pa"] for row in rows]
    ax.plot(wall_steps, drops, "o-", color=TEAL, linewidth=1.9, label="PAACS section")
    ax.axhline(reference, color=RUST, dashes=(3, 3), linewidth=1.35, label="Square-duct series")
    for index, (row, spacing, drop) in enumerate(zip(rows, wall_steps, drops)):
        offset = 7 if drop < reference or index < len(rows) - 1 else -12
        ax.annotate(f"{row['prism6']:,}", (spacing, drop), xytext=(0, offset),
                    textcoords="offset points", ha="center", fontsize=8, color=SUB)
    chart(ax, "Pressure drop at matched bulk speed", "First wall spacing (mm) · finer →",
          "Pressure drop (Pa)")
    ax.margins(x=.12, y=.22)
    ax.invert_xaxis()
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.legend(frameon=False, fontsize=8.5, loc="best")


def render(run, output):
    summary = json.loads((run / "summary.json").read_text(encoding="utf-8-sig"))
    rows = sorted(summary["meshes"], key=lambda row: row["first_wall_spacing_m"], reverse=True)
    final = rows[-1]
    data = load_solution(run, final)
    geometry, fluid = summary["geometry"], summary["fluid"]
    side, length = geometry["side_m"] * 1000, geometry["length_m"] * 1000
    reference = summary["reference"]["pressure_drop_pa"]
    norm = Normalize(0., max(data["section_velocity"].max(), data["profile_exact_velocity"].max()))

    fig = plt.figure(figsize=(13.1, 8.9), facecolor=BG)
    fig.text(.049, .936, "Square duct with wall-graded prisms", fontsize=26, weight="bold", color=INK)
    fig.text(.050, .890, (f"{side:g} × {side:g} mm section · {length:g} mm length · "
             f"Bulk speed {fluid['bulk_speed_m_s']:g} m/s · Reynolds number {fluid['reynolds']:g}"),
             fontsize=11, color=SUB)
    fig.text(.050, .841, f"Prism6 cutaway  ·  {final['prism6']:,} cells",
             fontsize=13, weight="bold", color=INK)
    fig.text(.619, .841, "Wall layers and axial velocity", fontsize=13, weight="bold", color=INK)
    cutaway(fig.add_axes([.050, .503, .523, .306], facecolor=BG), data["points"] * 1000, data["prisms"])
    ax = fig.add_axes([.619, .503, .265, .300], facecolor=BG)
    field = cross_section(ax, data, norm)
    cax = fig.add_axes([.908, .526, .016, .255])
    cb = fig.colorbar(field, cax=cax)
    cb.set_label("Axial speed (m/s)", color=SUB, fontsize=9)
    cb.ax.tick_params(colors=SUB, labelsize=8)
    cb.outline.set_visible(False)
    fig.text(.050, .478, "Cyan exterior · Gold streamwise section · Flow along +z",
             fontsize=9, color=SUB)
    fig.text(.619, .436, (f"First wall interval {final['first_wall_spacing_m'] * 1000:.2f} mm"
                          f" · {final['section_tri3']:,} cross-section triangles"),
             fontsize=9, color=SUB)

    velocity_profile(fig.add_axes([.080, .118, .414, .259]), data)
    pressure_drop(fig.add_axes([.580, .118, .362, .259]), rows, reference)
    fig.text(.050, .034, ("Fully developed laminar solution on the extracted Prism6 cross-section. "
                         "Pressure labels show Prism6 counts."),
             fontsize=9, color=SUB)
    fig.savefig(output, dpi=160, facecolor=BG)
    plt.close(fig)


def render_wall_normal(archive, output):
    mesh = read_archive(archive)
    cells = [cell for cell in mesh.cells if cell.dim == 3]
    if not cells or any(cell.kind != 7 for cell in cells):
        raise ValueError("The wall-normal companion needs an all-Prism6 archive")
    prisms = np.asarray([cell.ids for cell in cells])
    # Rotate the x-axial archive into the z-axial projection coordinates.
    points = mesh.points[:, [2, 1, 0]].copy()
    points[:, 0] *= -1
    wall = points[prisms[:, :3], 1]
    across = points[prisms[:, 3:], 1] - wall
    span = points[prisms[:, 3:]][:, :, [0, 2]] - points[prisms[:, :3]][:, :, [0, 2]]
    if not np.allclose(span, 0., rtol=0., atol=1e-9) or np.any(across <= 0):
        raise ValueError("Prism sweep directions do not match the opposing y walls")
    schedule = np.unique(points[:, 1])
    gaps = np.diff(schedule)
    wall_rows = 5
    if len(gaps) <= wall_rows*2 or not np.allclose(gaps, gaps[::-1]):
        raise ValueError("Expected two mirrored wall-normal layer schedules")
    lo, hi = points.min(axis=0), points.max(axis=0)
    close_to_wall = schedule[wall_rows] - schedule[0]
    def color_section(poly):
        center = poly[:, 1].mean()
        return GOLD if min(center - lo[1], hi[1] - center) <= close_to_wall else "#c5e2e4"

    plane = lo[2] + .58*(hi[2]-lo[2])
    cuts = prism_sections(points, prisms, 2, plane)
    fig = plt.figure(figsize=(11.6, 4.7), facecolor=BG)
    fig.text(.049, .899, "Wall-normal Prism6 layers on two opposing walls",
             fontsize=22, color=INK, weight="bold")
    fig.text(.050, .824, (f"{len(prisms):,} prisms · First layer {gaps[0]:.2f} mm · "
                         f"{wall_rows} graded rows per wall · Maximum growth "
                         f"{np.maximum(gaps[1:]/gaps[:-1], gaps[:-1]/gaps[1:]).max():.2f}"),
             fontsize=10.1, color=SUB)
    cutaway(fig.add_axes([.047, .160, .577, .588], facecolor=BG), points, prisms, color_section)

    ax = fig.add_axes([.700, .153, .240, .602], facecolor=BG)
    ax.add_collection(PolyCollection([poly[:, :2] for poly in cuts],
                                     facecolors=[color_section(poly) for poly in cuts],
                                     edgecolors="#325362", linewidths=.36))
    span_xy = hi[:2]-lo[:2]
    ax.set(xlim=(lo[0] - .08*span_xy[0], hi[0] + .08*span_xy[0]),
           ylim=(lo[1] - .12*span_xy[1], hi[1] + .12*span_xy[1]), aspect="equal")
    ax.set_axis_off()
    ax.set_title("End-on section", loc="left", pad=0, fontsize=11, color=INK, weight="bold")
    for y, direction in [(lo[1], 1), (hi[1], -1)]:
        ax.annotate("", xy=(hi[0] + .044*span_xy[0], y + direction*close_to_wall),
                    xytext=(hi[0] + .044*span_xy[0], y),
                    arrowprops=dict(arrowstyle="-|>", color=RUST, linewidth=1.4))
    fig.text(.050, .067, "Gold marks the graded rows at the two y walls. Arrows show their inward wall normals.",
             fontsize=9, color=SUB)
    fig.savefig(output, dpi=160, facecolor=BG)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Solve directory or a wall-normal mesh archive")
    parser.add_argument("output", type=Path, help="Output figure path")
    parser.add_argument("--wall-normal", action="store_true", help="Plot the opposed-wall Prism6 mesh")
    args = parser.parse_args()
    if args.wall_normal:
        render_wall_normal(args.source, args.output)
    else:
        render(args.source, args.output)


if __name__ == "__main__":
    main()
