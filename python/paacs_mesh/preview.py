"""Headless, tagged CAD previews for checking selections before meshing."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from dataclasses import dataclass
from html import escape
import json
import math
import os
from pathlib import Path
import tempfile

from ._cad_common import SelectionError, boolean, integer, json_text, number, protected_destination


def _project_points(points):
    """Normalize before projecting, so finite large placements cannot overflow."""
    if not points:
        raise ValueError("A CAD preview must contain display points")
    magnitude = max((abs(value) for point in points for value in point), default=1) or 1
    scaled = [tuple(value / magnitude for value in point) for point in points]
    centre = [min(point[i] for point in scaled) / 2 + max(point[i] for point in scaled) / 2
              for i in range(3)]
    local = [tuple(point[i] - centre[i] for i in range(3)) for point in scaled]
    extent = max((abs(value) for point in local for value in point), default=1) or 1
    local = [tuple(value / extent for value in point) for point in local]
    projected = [(0.8660254037844386 * (x - y), .5 * (x + y) - z, x + y + z)
                 for x, y, z in local]
    return local, projected


def _visible_labels(projected, triangles, tags, selected):
    """Find visible triangle-centroid anchors with bounded spatial/depth tests.

    Labels are an aid to selection, not a visibility or geometry certificate.
    A dense scene may omit labels rather than enter an unbounded all-pairs scan.
    """
    cells, max_references, max_tests, max_faces = 32, 1000000, 1000000, 1000
    low = [min(point[i] for point in projected) for i in range(2)]
    span = [max(max(point[i] for point in projected) - low[i], 1e-30) for i in range(2)]
    grid = [[] for _ in range(cells * cells)]
    candidates = {}
    references = 0

    def cell(value, axis):
        return min(cells - 1, max(0, int((value - low[axis]) / span[axis] * cells)))

    for index, tri in enumerate(triangles):
        a, b, c = [projected[n] for n in tri]
        area = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
        if area <= 1e-30:
            continue
        x0, x1 = cell(min(a[0], b[0], c[0]), 0), cell(max(a[0], b[0], c[0]), 0)
        y0, y1 = cell(min(a[1], b[1], c[1]), 1), cell(max(a[1], b[1], c[1]), 1)
        references += (x1 - x0 + 1) * (y1 - y0 + 1)
        if references > max_references:
            return {}, True
        for y in range(y0, y1 + 1):
            for x in range(x0, x1 + 1):
                grid[y * cells + x].append(index)
        tag = tags[index]
        if selected and tag not in selected:
            continue
        if tag not in candidates:
            if len(candidates) == max_faces:
                continue
            candidates[tag] = []
        choices = candidates[tag]
        choices.append((area, index))
        choices.sort(reverse=True)
        del choices[8:]

    depth_tolerance = max(max(p[2] for p in projected) - min(p[2] for p in projected), 1) * 1e-10
    visible, tests = {}, 0
    for tag, choices in candidates.items():
        for _, index in choices:
            tri = triangles[index]
            x, y, depth = [sum(projected[n][i] for n in tri) / 3 for i in range(3)]
            hidden = False
            for other in grid[cell(y, 1) * cells + cell(x, 0)]:
                if other == index:
                    continue
                tests += 1
                if tests > max_tests:
                    return visible, True
                a, b, c = [projected[n] for n in triangles[other]]
                determinant = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1])
                if abs(determinant) <= 1e-30:
                    continue
                u = ((b[1] - c[1]) * (x - c[0]) + (c[0] - b[0]) * (y - c[1])) / determinant
                v = ((c[1] - a[1]) * (x - c[0]) + (a[0] - c[0]) * (y - c[1])) / determinant
                w = 1 - u - v
                if min(u, v, w) >= -1e-10 and u * a[2] + v * b[2] + w * c[2] > depth + depth_tolerance:
                    hidden = True
                    break
            if not hidden:
                visible[tag] = index
                break
    return visible, len(candidates) == max_faces


@dataclass(frozen=True)
class CadPreview:
    """Owned display triangles, explicitly separate from a simulation mesh."""

    __slots__ = ("_owner", "_json")

    _owner: object
    _json: str

    def to_dict(self):
        return json.loads(self._json)

    def write_svg(self, path, *, width=960, height=720, labels=True, title="PAACS CAD selection"):
        """Write a static isometric selection preview without GUI dependencies."""
        width, height = integer(width, "width"), integer(height, "height")
        if width > 4096 or height > 4096 or width < 200 or height < 200:
            raise ValueError("SVG dimensions must be 200..4096")
        labels = boolean(labels, "labels")
        if (not isinstance(title, str) or len(title) > 1024 or
                any(ord(c) < 32 and c not in "\t\n\r" for c in title)):
            raise ValueError("title must contain at most 1024 characters")
        data = self.to_dict()
        points, triangles, tags = data["points"], data["triangles"], data["face_tags"]
        if len(triangles) > 100000:
            raise ValueError("SVG output is limited to 100000 display triangles")
        selected = set(data["selected_face_tags"])
        local, projected = _project_points(points)
        low = [min(p[i] for p in projected) for i in range(2)]
        high = [max(p[i] for p in projected) for i in range(2)]
        span = [max(high[i] - low[i], 1e-30) for i in range(2)]
        scale = min((width - 64) / span[0], (height - 128) / span[1])
        centre = [(low[i] + high[i]) / 2 for i in range(2)]
        coords = [((p[0] - centre[0]) * scale + width / 2,
                   (p[1] - centre[1]) * scale + height / 2 + 16) for p in projected]
        order = sorted(range(len(triangles)), key=lambda i: sum(projected[n][2] for n in triangles[i]))
        xml = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img">',
               f'<title>{escape(title)}</title>',
               '<rect width="100%" height="100%" fill="#f3f6fa"/>',
               f'<text x="24" y="32" font-family="sans-serif" font-size="20" fill="#17364b">{escape(title)}</text>',
               '<text x="24" y="54" font-family="sans-serif" font-size="12" fill="#526d80">Display tessellation · gold = selected faces · units mm</text>']
        for index in order:
            tri = triangles[index]
            face_tag = tags[index]
            a, b, c = [local[n] for n in tri]
            u, v = [b[i] - a[i] for i in range(3)], [c[i] - a[i] for i in range(3)]
            normal = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
            norm = math.hypot(*normal)
            light = .76 + .2 * abs(sum(normal) / (math.sqrt(3) * norm)) if norm else .76
            base = (242, 196, 101) if face_tag in selected else (115, 188, 197)
            colour = "#" + "".join(f"{round(channel * light):02x}" for channel in base)
            polygon = " ".join(f"{coords[n][0]:.3f},{coords[n][1]:.3f}" for n in tri)
            xml.append(f'<polygon points="{polygon}" fill="{colour}" stroke="#365365" stroke-opacity="0.22" stroke-width="0.35"/>')
        label_limited = False
        if labels:
            visible, label_limited = _visible_labels(projected, triangles, tags, selected)
            for tag, index in visible.items():
                tri = triangles[index]
                x, y = [sum(coords[n][i] for n in tri) / 3 for i in range(2)]
                xml.append(f'<text x="{x:.3f}" y="{y:.3f}" font-family="sans-serif" font-size="13" '
                           f'text-anchor="middle" fill="#17364b" stroke="#ffffff" stroke-width="3" '
                           f'paint-order="stroke">face {tag}</text>')
        xml.append(f'<text x="24" y="{height - 20}" font-family="sans-serif" font-size="11" fill="#526d80">'
                   f'Source {data["source_sha256"][:16]} · {len(triangles):,} display triangles · not a simulation mesh</text>')
        label_note = "label budget reached; some labels omitted" if label_limited else "hidden or very small face labels may be omitted"
        xml.append(f'<text x="24" y="{height - 5}" font-family="sans-serif" font-size="10" fill="#526d80">'
                   f'Isometric display projection · {label_note}</text></svg>')
        payload = "\n".join(xml).encode("utf-8")
        if len(payload) > 32 * 1024 * 1024:
            raise ValueError("SVG exceeds the 32 MiB display budget")
        destination = protected_destination(path, self._owner.info["source_path"])
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode="wb", dir=destination.parent,
                                             prefix=f".{destination.name}.", delete=False) as stream:
                temporary = Path(stream.name)
                stream.write(payload)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, destination)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        return destination


def preview_model(model, selection=None, *, linear_deflection=.02, max_points=50000, max_triangles=100000):
    """Copy-tessellate CAD for display; serialized preview data is limited to 8 MiB.

    Output limits do not impose a hard memory/time bound on the CAD tessellator.
    The preview does not certify mesh validity or geometric approximation error.
    """
    model._assert_identity()
    if selection is not None:
        model._require_selection(selection)
        selection.validate()
    deflection = number(linear_deflection, "linear_deflection", minimum=0)
    if deflection <= 0:
        raise ValueError("linear_deflection must be positive")
    points, triangles = integer(max_points, "max_points"), integer(max_triangles, "max_triangles")
    if points > 250000 or triangles > 500000:
        raise ValueError("Preview limits must not exceed 250000 points or 500000 triangles")
    try:
        data = model._native._preview(linear_deflection=deflection, max_points=points, max_triangles=triangles)
        data.update(schema="paacs-cad-preview-v1", display_geometry_only=True, source_model_unmodified=True,
                    source_sha256=model.source_sha256, kernel_version=model.kernel_version,
                    units="mm", selected_face_tags=list(selection.tags) if selection is not None else [])
        return CadPreview(model, json_text(data))
    except (RuntimeError, ValueError) as error:
        raise SelectionError("preview_failed", str(error), stage="preview") from error
