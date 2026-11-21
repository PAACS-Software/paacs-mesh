"""Tagged display previews preserve the CAD used by the native mesher."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
import xml.etree.ElementTree as ET

import paacs_mesh as pm
from test_native import fixture_path


SVG = "{http://www.w3.org/2000/svg}"


def triangle_normal(points, indices):
    a, b, c = [points[n] for n in indices]
    u, v = [b[i] - a[i] for i in range(3)], [c[i] - a[i] for i in range(3)]
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]


class CadPreviewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = fixture_path()
        cls.model = pm.read_step(cls.fixture)

    def test_owned_tagged_display_and_reference_orientation(self):
        selection = self.model.select_faces(surface="cylinder")
        preview = self.model.preview(selection, linear_deflection=.04)
        self.assertIsInstance(preview, pm.CadPreview)
        with self.assertRaises(AttributeError):
            preview.extra = True
        with self.assertRaises(AttributeError):
            selection.tags = (999,)
        data = preview.to_dict()
        self.assertTrue(data["display_geometry_only"])
        self.assertTrue(data["source_model_unmodified"])
        self.assertEqual(data["source_sha256"], self.model.source_sha256)
        self.assertEqual(data["selected_face_tags"], list(selection.tags))
        self.assertEqual(data["units"], "mm")
        self.assertEqual(data["linear_deflection"], .04)
        self.assertNotIn("valid", data)
        points, triangles, tags = data["points"], data["triangles"], data["face_tags"]
        self.assertGreater(len(points), 0)
        self.assertEqual(len(triangles), len(tags))
        self.assertEqual(set(tags), set(self.model.face_tags))
        self.assertTrue(all(len(point) == 3 and all(math.isfinite(v) for v in point) for point in points))
        records = {face["tag"]: face for face in self.model.inspect()["entities"]}
        for triangle, tag in zip(triangles, tags):
            self.assertEqual(len(set(triangle)), 3)
            self.assertTrue(all(0 <= n < len(points) for n in triangle))
            normal = triangle_normal(points, triangle)
            self.assertGreater(math.hypot(*normal), 0)
            if records[tag]["normal"] is not None:
                self.assertGreater(sum(a * b for a, b in zip(normal, records[tag]["normal"])), 0)
        data["face_tags"].clear()
        self.assertEqual(set(preview.to_dict()["face_tags"]), set(self.model.face_tags))

    def test_success_and_refusal_leave_source_and_meshing_unchanged(self):
        model = pm.read_step(self.fixture)
        original_bytes, original_info = self.fixture.read_bytes(), model.info
        with tempfile.TemporaryDirectory() as directory:
            before, after = Path(directory) / "before.paacsmesh", Path(directory) / "after.paacsmesh"
            model.mesh(.7).write(before)
            model.preview(linear_deflection=.02)
            for options in ({"max_points": 1}, {"max_triangles": 1}):
                with self.subTest(options=options), self.assertRaises(pm.SelectionError) as caught:
                    model.preview(**options)
                self.assertEqual(caught.exception.reason, "preview_failed")
            self.assertEqual(model.info, original_info)
            self.assertEqual(self.fixture.read_bytes(), original_bytes)
            model.mesh(.7).write(after)
            self.assertEqual(before.read_bytes(), after.read_bytes())

    def test_assembly_placements_and_face_identity(self):
        model = pm.read_step(self.fixture.with_name("assembly.step"))
        data = model.preview(linear_deflection=.05).to_dict()
        records = {face["tag"]: face for face in model.inspect()["entities"]}
        self.assertEqual(set(data["face_tags"]), set(model.face_tags))
        self.assertEqual({owner["body_tag"] for face in records.values() for owner in face["owners"]}, {1, 2})
        self.assertAlmostEqual(min(p[0] for p in data["points"]), 0, places=8)
        self.assertAlmostEqual(max(p[0] for p in data["points"]), 5, places=8)
        for triangle, tag in zip(data["triangles"], data["face_tags"]):
            face = records[tag]
            self.assertEqual(len(face["owners"]), 1)
            for n in triangle:
                p = data["points"][n]
                self.assertTrue(all(lo - 1e-7 <= v <= hi + 1e-7
                                    for lo, v, hi in zip(face["bounds"]["lower"], p, face["bounds"]["upper"])))
                self.assertTrue(p[0] <= 2 + 1e-7 or p[0] >= 4 - 1e-7)
            normal = triangle_normal(data["points"], triangle)
            self.assertGreater(sum(a * b for a, b in zip(normal, face["normal"])), 0)

    def test_selection_and_request_admission(self):
        foreign = pm.read_step(self.fixture).select_faces(surface="cylinder")
        with self.assertRaises(pm.SelectionError) as caught:
            self.model.preview(foreign)
        self.assertEqual(caught.exception.reason, "foreign_selection")
        for options in ({"linear_deflection": 0}, {"linear_deflection": float("nan")},
                        {"max_points": 0}, {"max_triangles": True},
                        {"max_points": 250001}, {"max_triangles": 500001}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.model.preview(**options)

    def test_svg_export_and_source_protection(self):
        selected = self.model.select_faces(surface="cylinder")
        preview = self.model.preview(selected)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "preview.svg"
            self.assertEqual(preview.write_svg(path, title='Cylinder <wall> & "selection"'), path)
            root = ET.parse(path).getroot()
            self.assertEqual(root.find(SVG + "title").text, 'Cylinder <wall> & "selection"')
            self.assertEqual(len(root.findall(SVG + "polygon")), len(preview.to_dict()["triangles"]))
            self.assertIn("not a simulation mesh", path.read_text(encoding="utf-8"))
            old = path.read_bytes()
            with self.assertRaises(ValueError):
                preview.write_svg(path, width=10)
            self.assertEqual(path.read_bytes(), old)
        with self.assertRaises(ValueError):
            preview.write_svg(self.fixture)


class SvgProjectionTests(unittest.TestCase):
    """Synthetic display triangles isolate rendering behavior, without CAD claims."""

    def preview(self, points, triangles, tags, selected=()):
        owner = SimpleNamespace(info={"source_path": str(fixture_path())})
        data = {"points": points, "triangles": triangles, "face_tags": tags,
                "selected_face_tags": list(selected), "source_sha256": "a" * 64}
        return pm.CadPreview(owner, json.dumps(data))

    def test_hidden_faces_do_not_receive_labels_over_front_faces(self):
        # The second triangle has the same projection, nearer the camera.
        points = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 1), (2, 1, 1), (1, 2, 1)]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "occlusion.svg"
            self.preview(points, [(0, 1, 2), (3, 4, 5)], [1, 2]).write_svg(path)
            texts = [node.text for node in ET.parse(path).getroot().findall(SVG + "text")]
            self.assertIn("face 2", texts)
            self.assertNotIn("face 1", texts)
            self.preview(points, [(0, 1, 2), (3, 4, 5)], [1, 2], [1]).write_svg(path)
            texts = [node.text for node in ET.parse(path).getroot().findall(SVG + "text")]
            self.assertFalse(any(text.startswith("face ") for text in texts))

    def test_large_finite_coordinates_project_without_overflow(self):
        points = [(1e308, 0, 0), (0, 1e308, 0), (0, 0, 1e308)]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "large.svg"
            self.preview(points, [(0, 1, 2)], [1]).write_svg(path)
            root = ET.parse(path).getroot()
            for node in root.findall(SVG + "polygon"):
                coordinates = node.attrib["points"].replace(",", " ").split()
                self.assertTrue(all(math.isfinite(float(value)) for value in coordinates))


if __name__ == "__main__":
    unittest.main()
