"""Deterministic, bounded SVG views of actual archived mesh elements."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from html import escape
from itertools import product
import hashlib
import json
import math
from pathlib import Path
import time

from ._cad_common import bounds, integer, number, sequence
from .mesh_archive import InspectionError, atomic_stream
from .mesh_inspection import METRICS, MeshInspection, TET_FACES, _centroid
from .quality_metrics import METRIC_CONVENTIONS


_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
_LABELS = {
    "mean_ratio": ("Mean ratio", "higher is better; regular Tet4 = 1", 1.),
    "scaled_jacobian": ("Signed scaled Jacobian", "higher is better; regular Tet4 = 1", 1.),
    "minimum_dihedral_deg": ("Minimum interior dihedral angle (degrees)", "regular Tet4 = 70.529 degrees", 70.52877936550931),
}


def _colour(value, ideal, metric):
    fraction = max(0., min(1., value / ideal))
    # Fixed interpolation; the colour map never replaces the measured value.
    low, high = (195, 45, 49), (235, 163, 57)
    if metric == "scaled_jacobian" and value < 0:
        low, high = (97, 47, 128), low
        fraction = max(0., min(1., value + 1))
    return "#" + "".join(f"{round(a + (b - a) * fraction):02x}" for a, b in zip(low, high))


def _camera(frame, azimuth, elevation, width, height):
    magnitude = max(abs(value) for point in frame for value in point) or 1.
    lower, upper = [[value / magnitude for value in point] for point in frame]
    center = [a / 2 + b / 2 for a, b in zip(lower, upper)]
    extent = max(b - a for a, b in zip(lower, upper)) or 1.
    az, el = math.radians(azimuth), math.radians(elevation)
    right = (-math.sin(az), math.cos(az), 0.)
    up = (-math.sin(el) * math.cos(az), -math.sin(el) * math.sin(az), math.cos(el))
    toward = (math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el))

    def project(point):
        local = [(point[i] / magnitude - center[i]) / extent for i in range(3)]
        result = [sum(local[i] * vector[i] for i in range(3)) for vector in (right, up, toward)]
        if not all(math.isfinite(value) for value in result):
            raise InspectionError("projection_range", "Coordinates cannot be represented in this camera frame")
        return result

    corners = [project(point) for point in product(*zip(*frame))]
    low = [min(point[i] for point in corners) for i in range(2)]
    high = [max(point[i] for point in corners) for i in range(2)]
    available_width, available_height = width - 64, height - 208
    scale = min(available_width / (high[0] - low[0] or 1.), available_height / (high[1] - low[1] or 1.))
    middle = [(a / 2 + b / 2) for a, b in zip(low, high)]

    def pixel(point):
        x, y, depth = project(point)
        return ((x - middle[0]) * scale + width / 2,
                -(y - middle[1]) * scale + 112 + available_height / 2, depth)

    return pixel


def visualize_mesh(inspection, path, *, metric="mean_ratio", highlight_below=.2, cad=None,
                   cut_axis=None, cut_value=None, view_bounds=None, azimuth=35, elevation=25,
                   width=1200, height=900, max_triangles=100000, max_highlights=10000,
                   max_output_bytes=32 * 1024 * 1024, highlight_cell_ids=None):
    """Write actual cells, with all selected poor cells overlaid as an x-ray.

    The translucent context is the native boundary (derived from Tet4 faces if
    absent) or a source-matched CAD display tessellation. A cut retains whole
    Tet4 cells whose centroid coordinate is <= cut_value; it never cuts a cell
    with an artificial plane. Highlights remain an explicit full-mesh x-ray,
    including cells removed from the context by that cut. A close-up uses the
    supplied world frame and clips only drawing outside the SVG viewport.
    Reuse view_bounds, angles and dimensions for comparable before/after views.
    No native validity or simulation certification is performed.
    """
    if not isinstance(inspection, MeshInspection):
        raise ValueError("inspection must be a MeshInspection")
    if metric not in METRICS:
        raise ValueError("Unknown Tet4 quality metric")
    threshold = number(highlight_below, "highlight_below")
    width, height = integer(width, "width"), integer(height, "height")
    if not 320 <= width <= 8192 or not 320 <= height <= 8192:
        raise ValueError("SVG dimensions must be 320..8192")
    max_triangles = integer(max_triangles, "max_triangles")
    max_highlights = integer(max_highlights, "max_highlights", minimum=0)
    max_output_bytes = integer(max_output_bytes, "max_output_bytes")
    units = inspection.units
    if units is not None and (not isinstance(units, str) or any(
            not (character in "\t\n\r" or 0x20 <= ord(character) <= 0xD7FF or
                 0xE000 <= ord(character) <= 0xFFFD or 0x10000 <= ord(character) <= 0x10FFFF)
            for character in units)):
        raise InspectionError("display_metadata", "Stored units are not valid XML text")
    azimuth, elevation = number(azimuth, "azimuth"), number(elevation, "elevation")
    if not -90 <= elevation <= 90:
        raise ValueError("elevation must be between -90 and 90 degrees")
    azimuth %= 360
    if cut_axis is not None:
        if isinstance(cut_axis, bool) or cut_axis not in ("x", "y", "z", 0, 1, 2):
            raise ValueError("cut_axis must be x, y, z or 0, 1, 2")
        cut_axis = ("x", "y", "z").index(cut_axis) if isinstance(cut_axis, str) else int(cut_axis)
        cut_value = number(cut_value, "cut_value")
    elif cut_value is not None:
        raise ValueError("cut_value requires cut_axis")
    frame = bounds(view_bounds, "view_bounds") if view_bounds is not None else None
    selected = None
    if highlight_cell_ids is not None:
        values = sequence(highlight_cell_ids, "highlight_cell_ids", maximum=max_highlights)
        selected = {integer(value, "cell_id", minimum=0) for value in values}
        if len(selected) != len(values) or any(index >= inspection.cell_count for index in selected):
            raise ValueError("highlight_cell_ids must contain distinct IDs in this mesh")

    archive = inspection._archive
    protected = list(archive.protected_paths)
    tet_count, native_triangles = archive.type_counts[5], archive.type_counts[2]
    graph_needed = cut_axis is not None or (not native_triangles and cad is None)
    graph_faces = 4 * tet_count if graph_needed else 0
    face_bound = native_triangles + graph_faces
    admitted_faces = min(max_triangles, face_bound)
    highlight_bound = min(tet_count, max_highlights)
    # Includes dictionary/tuple graph entries, depth sorting, independent world
    # coordinates, metadata encoding/escaping, and streamed polygon records.
    scratch = 65536 + 768 * graph_faces + 1024 * admitted_faces + 4096 * highlight_bound
    cad_points = min(250000, max_triangles * 3) if cad is not None else 0
    if cad is not None:
        scratch += 384 * (cad_points + max_triangles)
        if inspection.source_sha256 is None or getattr(cad, "source_sha256", None) != inspection.source_sha256:
            raise InspectionError("source_mismatch", "CAD context requires the same known source SHA-256 as the inspection")
        source_path = cad.info.get("source_path")
        if source_path:
            protected.append(Path(source_path))
    inspection._reserve(scratch, "Mesh SVG visualization")
    started, work = time.perf_counter(), 0

    def charge(amount=1):
        nonlocal work
        work += amount
        inspection._deadline(started, work)

    context, graph, highlights = [], {}, []
    missing_metrics = unsupported_cells = retained = 0
    observed_lower, observed_upper = [math.inf] * 3, [-math.inf] * 3
    for node_id in range(inspection.node_count):
        charge()
        point, _ = archive.node(node_id)
        for axis in range(3):
            observed_lower[axis] = min(observed_lower[axis], point[axis])
            observed_upper[axis] = max(observed_upper[axis], point[axis])
    if frame is None:
        if not inspection.node_count:
            raise InspectionError("display_geometry", "Mesh has no coordinates to display")
        frame = [observed_lower, observed_upper]

    for cell_id in range(inspection.cell_count):
        charge()
        kind, nodes, classification = archive.cell(cell_id)
        if kind == 2 and classification["dimension"] == 2 and cad is None:
            if len(context) >= max_triangles:
                raise InspectionError("triangle_budget", "Native boundary exceeds max_triangles; no triangles were omitted")
            context.append((tuple(nodes), cell_id, classification["tag"]))
        if kind != 5:
            unsupported_cells += kind not in (0, 1, 2, 3, 4)
            if selected is not None and cell_id in selected:
                raise InspectionError("unsupported_metric", f"Selected cell {cell_id} is not Tet4")
            continue
        points = [archive.node(index)[0] for index in nodes]
        value = inspection._metric(cell_id, metric)
        if value is None:
            missing_metrics += 1
            if selected is not None and cell_id in selected:
                raise InspectionError("unsupported_metric", f"Selected cell {cell_id} has no finite {metric}")
        is_highlight = cell_id in selected if selected is not None else value is not None and value < threshold
        if is_highlight:
            if len(highlights) >= max_highlights:
                raise InspectionError("highlight_budget", "All matching cells exceed max_highlights; no cells were silently omitted")
            highlights.append({"cell_id": cell_id, "node_ids": list(nodes), "points": points,
                               "body_tag": classification["tag"] if classification["dimension"] == 3 else None,
                               "metric": metric, "value": value})
        keep = cut_axis is None or _centroid(points)[cut_axis] <= cut_value
        retained += keep
        if graph_needed:
            for face in TET_FACES:
                charge()
                face_nodes = tuple(nodes[index] for index in face)
                key = tuple(sorted(face_nodes))
                if key not in graph:
                    graph[key] = [face_nodes, cell_id, 0, 0, None, None]
                record = graph[key]
                record[2] += 1
                if keep:
                    record[3] += 1
                    record[4], record[5] = face_nodes, cell_id

    cut_faces = []
    for record in graph.values():
        charge()
        nodes, cell_id, uses, kept_uses, kept_nodes, kept_cell_id = record
        if uses == 1 and not native_triangles and cad is None:
            if len(context) + len(cut_faces) >= max_triangles:
                raise InspectionError("triangle_budget", "Whole-cell boundary exceeds max_triangles; no triangles were omitted")
            context.append((nodes, cell_id, None))
        if cut_axis is not None and uses > 1 and kept_uses == 1:
            if len(context) + len(cut_faces) >= max_triangles:
                raise InspectionError("triangle_budget", "Whole-cell boundary exceeds max_triangles; no triangles were omitted")
            cut_faces.append((kept_nodes, kept_cell_id, None))

    if cad is None and (archive.type_counts[3] or archive.type_counts[4] or
                        (unsupported_cells and not native_triangles)):
        raise InspectionError("unsupported_context", "A complete context for these cell families requires source-matched CAD display geometry")

    cad_metadata = None
    if cad is not None:
        charge()
        display = cad.preview(linear_deflection=.02, max_points=cad_points, max_triangles=max_triangles).to_dict()
        charge()
        if display.get("source_sha256") != inspection.source_sha256:
            raise InspectionError("source_mismatch", "CAD preview source identity differs from the inspection")
        preview_points, triangles = display["points"], display["triangles"]
        if len(preview_points) > cad_points or len(triangles) + len(cut_faces) > max_triangles:
            raise InspectionError("triangle_budget", "CAD context plus cut surface exceeds max_triangles")
        cad_metadata = {"display_geometry_only": True, "source_sha256": display["source_sha256"],
                        "linear_deflection": display.get("linear_deflection"),
                        "angular_deflection_radians": display.get("angular_deflection_radians"),
                        "triangle_count": len(triangles), "units": display.get("units"),
                        "scope": "CAD display tessellation; deflection is not a certified mesh error bound"}
        for index, nodes in enumerate(triangles):
            charge()
            if len(nodes) != 3 or any(type(node) is not int or not 0 <= node < len(preview_points) for node in nodes):
                raise InspectionError("display_geometry", "CAD preview has invalid triangle connectivity")
            context.append((tuple(tuple(number(value, "CAD display coordinate") for value in preview_points[node]) for node in nodes),
                            None, display["face_tags"][index]))
    if not context and not cut_faces:
        raise InspectionError("display_geometry", "No supported boundary or Tet4 context is available")

    pixel = _camera(frame, azimuth, elevation, width, height)

    def coordinates(record, is_cad=False):
        return record[0] if is_cad else [archive.node(index)[0] for index in record[0]]

    # O(n log n) depth ordering is pre-admitted; no geometry is synthesized.
    charge(len(context) * max(1, len(context).bit_length()) +
           len(cut_faces) * max(1, len(cut_faces).bit_length()))
    context.sort(key=lambda record: sum(pixel(point)[2] for point in coordinates(record, cad is not None)))
    cut_faces.sort(key=lambda record: sum(pixel(point)[2] for point in coordinates(record)))
    label, convention, ideal = _LABELS[metric]
    charge(len(highlights))
    selected_minimum = min((cell["value"] for cell in highlights), default=None)
    tick_values = ([-1., -.5, 0., .5, 1.] if metric == "scaled_jacobian" else
                   [0., 15., 30., 45., 60., ideal] if metric == "minimum_dihedral_deg" else
                   [0., .25, .5, .75, 1.])
    legend = {"minimum": tick_values[0], "maximum": tick_values[-1], "ticks": tick_values,
              "units": "degrees" if metric == "minimum_dihedral_deg" else "dimensionless",
              "zero_colour": _colour(0., ideal, metric),
              "negative_colour": _colour(-1., ideal, metric) if metric == "scaled_jacobian" else None}
    cropped_frame = any(frame[0][axis] > observed_lower[axis] or frame[1][axis] < observed_upper[axis]
                        for axis in range(3))
    metadata = {"schema": "paacs-mesh-visualization-v1", **inspection.identity(),
                "metric": metric, "metric_convention": METRIC_CONVENTIONS[metric],
                "highlight_selection": "explicit_cell_ids" if selected is not None else "metric_below_threshold",
                "highlight_below": threshold, "highlight_count": len(highlights), "highlights": highlights,
                "selected_metric_minimum": selected_minimum,
                "legend": legend,
                "highlight_scope": "all selected Tet4 cells, full-mesh x-ray overlay; independent of context cut; viewport may crop off-screen drawing",
                "unavailable_tet4_metric_count": missing_metrics, "unsupported_volume_cell_count": unsupported_cells,
                "context": "cad_display_tessellation" if cad is not None else "native_boundary_triangles" if native_triangles else "boundary_derived_from_actual_tet4_faces",
                "context_triangles": len(context), "cut_surface_triangles": len(cut_faces),
                "cut": None if cut_axis is None else {"axis": "xyz"[cut_axis], "value": cut_value,
                    "rule": "keep whole Tet4 cells with centroid[axis] <= value; expose original cell faces; no plane clipping",
                    "retained_tet4": retained, "total_tet4": tet_count},
                "camera": {"projection": "orthographic", "azimuth_degrees": azimuth,
                           "elevation_degrees": elevation, "view_bounds": frame,
                           "width": width, "height": height, "world_frame_explicit": view_bounds is not None,
                           "world_frame_crops_mesh": cropped_frame},
                "cad": cad_metadata, "native_validity_rechecked": False}
    charge()
    encoded_metadata = escape(json.dumps(metadata, sort_keys=True, separators=(",", ":"), allow_nan=False))
    total_bytes = 0
    digest = hashlib.sha256()
    destination = Path(path).absolute()
    with atomic_stream(destination, protected=protected) as stream:
        def emit(text):
            nonlocal total_bytes
            charge()
            data = (text + "\n").encode("utf-8")
            if len(data) > max_output_bytes - total_bytes:
                raise InspectionError("output_budget", "SVG exceeds max_output_bytes; destination was not replaced")
            stream.write(data)
            digest.update(data)
            total_bytes += len(data)

        emit(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">')
        emit('<title>PAACS Mesh quality inspection</title>')
        emit('<desc>Actual mesh elements. All highlighted cells are drawn as an x-ray overlay; CAD context, when provided, is display tessellation only.</desc>')
        emit(f'<metadata id="paacs-mesh-visualization">{encoded_metadata}</metadata>')
        emit('<rect width="100%" height="100%" fill="#f3f6fa"/>')
        emit(f'<defs><clipPath id="mesh-viewport"><rect x="24" y="112" width="{width - 48}" height="{height - 208}"/></clipPath>')
        emit('<linearGradient id="quality-colours" x1="0%" y1="0%" x2="100%" y2="0%">')
        for value in ([legend["minimum"], 0., legend["maximum"]] if metric == "scaled_jacobian" else [0., ideal]):
            offset = 100 * (value - legend["minimum"]) / (legend["maximum"] - legend["minimum"])
            emit(f'<stop offset="{offset:g}%" stop-color="{_colour(value, ideal, metric)}"/>')
        emit('</linearGradient></defs>')
        emit('<g font-family="sans-serif" fill="#17364b">')
        emit('<text x="28" y="32" font-size="22">PAACS Mesh · actual Tet4 inspection</text>')
        selection_note = f'threshold &lt; {threshold:.6g}' if selected is None else 'explicit native cell IDs'
        emit(f'<text x="28" y="56" font-size="13">{escape(label)} · {selection_note}</text>')
        emit(f'<text x="28" y="78" font-size="13">X-ray overlay: all {len(highlights):,} selected cells, including interiors; no occlusion filtering</text>')
        cut_label = ("Close-up viewport; full-object context clipped to frame" if cropped_frame else "Whole-object context")
        if cut_axis is not None:
            cut_label = ("Close-up viewport; " if cropped_frame else "") + f'Whole-cell cut: centroid {"xyz"[cut_axis]} <= {cut_value:.6g}; no plane clipping'
            if cropped_frame:
                cut_label += '; context clipped to frame'
        emit(f'<text x="28" y="98" font-size="12">{escape(cut_label)}</text></g>')
        legend_width = min(342, width - 56)
        legend_x = width - legend_width - 28
        legend_y = 16 if width >= 900 else height - 190
        bar_x, bar_y, bar_width = legend_x + 12, legend_y + 26, legend_width - 24
        emit(f'<g id="quality-legend" font-family="sans-serif" fill="#17364b"><rect x="{legend_x}" y="{legend_y}" width="{legend_width}" height="88" rx="4" fill="#ffffff" fill-opacity=".96" stroke="#c4d1db"/>')
        emit(f'<text x="{bar_x}" y="{legend_y + 17}" font-size="11">Numeric colour scale · higher is better</text>')
        emit(f'<rect x="{bar_x}" y="{bar_y}" width="{bar_width}" height="12" fill="url(#quality-colours)"/>')
        for index, value in enumerate(tick_values):
            x = bar_x + bar_width * (value - legend["minimum"]) / (legend["maximum"] - legend["minimum"])
            tick = f'{value:.3f}'.rstrip('0').rstrip('.')
            if metric == "minimum_dihedral_deg":
                tick += '°'
            anchor = "start" if index == 0 else "end" if index + 1 == len(tick_values) else "middle"
            emit(f'<line x1="{x:.4f}" y1="{bar_y + 12}" x2="{x:.4f}" y2="{bar_y + 17}" stroke="#17364b"/><text data-metric-tick="{value:.17g}" x="{x:.4f}" y="{bar_y + 29}" text-anchor="{anchor}" font-size="10">{tick}</text>')
        emit(f'<text x="{bar_x}" y="{legend_y + 78}" font-size="10">{escape(convention)}</text></g>')
        emit('<g clip-path="url(#mesh-viewport)">')

        def polygon(record, kind):
            points = coordinates(record, cad is not None and kind == "context")
            projected = [pixel(point) for point in points]
            encoded = " ".join(f"{x:.6f},{y:.6f}" for x, y, _ in projected)
            cell_id, face_tag = record[1], record[2]
            identity = f' data-cell-id="{cell_id}"' if cell_id is not None else ""
            node_identity = (f' data-node-ids="{",".join(str(index) for index in record[0])}"'
                             if not (cad is not None and kind == "context") else "")
            fill = "#69adb9" if kind == "context" else "#e9bd62"
            opacity = ".12" if kind == "context" else ".55"
            emit(f'<polygon data-layer="{kind}"{identity}{node_identity} points="{encoded}" fill="{fill}" fill-opacity="{opacity}" stroke="#527583" stroke-opacity=".35" stroke-width=".45"><title>{kind}; cell {cell_id}; CAD face {face_tag}</title></polygon>')

        for record in context:
            polygon(record, "context")
        for record in cut_faces:
            polygon(record, "whole-cell-cut")
        # Deliberately last: every selected cell remains visible through context.
        emit('<g id="xray-highlights" fill="none" stroke-width="1.8" stroke-linejoin="round">')
        for cell in highlights:
            charge(6)
            projected = [pixel(point) for point in cell["points"]]
            drawing = " ".join(f'M {projected[a][0]:.6f},{projected[a][1]:.6f} L {projected[b][0]:.6f},{projected[b][1]:.6f}' for a, b in _EDGES)
            colour = _colour(cell["value"], ideal, metric)
            emit(f'<path data-cell-id="{cell["cell_id"]}" data-node-ids="{",".join(str(index) for index in cell["node_ids"])}" data-metric="{metric}" data-value="{cell["value"]:.17g}" d="{drawing}" stroke="{colour}"><title>Cell {cell["cell_id"]}; {metric}={cell["value"]:.17g}; body {cell["body_tag"]}; node IDs {cell["node_ids"]}</title></path>')
        emit('</g></g>')
        if len(highlights) == 1:
            cell = highlights[0]
            detail = f'Native cell {cell["cell_id"]} · {metric} = {cell["value"]:.17g}'
        elif highlights:
            detail = f'Selected {metric} minimum: {selected_minimum:.6g}'
            if len(highlights) <= 4:
                detail += ' · native cell IDs: ' + ', '.join(str(cell["cell_id"]) for cell in highlights)
        else:
            detail = 'No finite Tet4 metric values selected'
        detail += f' · unavailable Tet4 metrics: {missing_metrics}; unsupported volume cells: {unsupported_cells}'
        emit(f'<text x="28" y="{height - 90}" font-family="sans-serif" font-size="11" fill="#17364b">{escape(detail)}</text>')
        units_label = units if units is not None else "unknown (bare archive has no stored units)"
        emit(f'<text x="28" y="{height - 70}" font-family="sans-serif" font-size="12" fill="#17364b">Units: {escape(units_label)} · native IDs · numeric quality legend · no simulation certificate</text>')
        context_note = (f'CAD display approximation; linear deflection {cad_metadata["linear_deflection"]}; not a certified error bound'
                        if cad is not None else 'Translucent context: actual native boundary or original Tet4 boundary faces')
        emit(f'<text x="28" y="{height - 49}" font-family="sans-serif" font-size="11" fill="#526d80">{escape(context_note)}</text>')
        emit(f'<text x="28" y="{height - 28}" font-family="monospace" font-size="10" fill="#526d80">Mesh SHA-256 {inspection.mesh_sha256}</text>')
        emit('</svg>')
        charge()
    return {"schema": "paacs-mesh-visualization-receipt-v1", "path": str(destination),
            **inspection.identity(), "metric": metric, "highlight_count": len(highlights),
            "highlighted_cell_ids": [cell["cell_id"] for cell in highlights],
            "highlight_selection": metadata["highlight_selection"], "highlight_below": threshold,
            "highlight_scope": metadata["highlight_scope"],
            "selected_metric_minimum": selected_minimum,
            "legend": legend,
            "unavailable_tet4_metric_count": missing_metrics, "unsupported_volume_cell_count": unsupported_cells,
            "context_triangles": len(context), "cut_surface_triangles": len(cut_faces),
            "camera": metadata["camera"], "cut": metadata["cut"], "cad": cad_metadata,
            "native_validity_rechecked": False, "bytes": total_bytes, "svg_sha256": digest.hexdigest()}
