"""Verify SVGs against original archive coordinates, IDs and measured metrics."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import json
import math
import os
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import Mock, patch
import xml.etree.ElementTree as ET

import paacs_mesh as pm
from paacs_mesh.mesh_archive import InspectionError
from paacs_mesh.mesh_inspection import inspect_mesh
from paacs_mesh.mesh_visualization import visualize_mesh


SVG = {"s": "http://www.w3.org/2000/svg"}


def write_archive(path, points, cells):
    payload = bytearray(struct.pack("<8sIIQQQ", b"PAACSM8\0", 1, 0, len(points), len(cells), sum(len(cell[1]) for cell in cells)))
    for point in points:
        payload.extend(struct.pack("<3dii", *point, 3, 11))
    for kind, nodes, dimension, tag in cells:
        payload.extend(struct.pack("<B3xii", kind, dimension, tag))
        payload.extend(struct.pack(f"<{len(nodes)}q", *nodes))
    path.write_bytes(payload)
    return path


def example(path):
    regular = [(3, 0, 0), (4, 0, 0), (3.5, math.sqrt(3) / 2, 0),
               (3.5, math.sqrt(3) / 6, math.sqrt(2 / 3))]
    points = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (.3, .3, .02), (.3, .3, -.02)] + regular
    cells = [(0, (0, 1), 1, 7), (5, (0, 1, 2, 3), 3, 11),
             (5, (0, 2, 1, 4), 3, 11), (5, (5, 6, 7, 8), 3, 22)]
    return write_archive(path, points, cells)


def svg_data(path):
    root = ET.parse(path).getroot()
    metadata = json.loads(root.find("s:metadata", SVG).text)
    paths = root.findall("s:g/s:g[@id='xray-highlights']/s:path", SVG)
    return root, metadata, paths


class MeshVisualizationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.source = example(self.root / "source.paacsmesh")
        self.inspection = inspect_mesh(self.source)

    def tearDown(self):
        self.inspection.close()
        self.directory.cleanup()

    def test_deterministic_svg_matches_python_query_and_original_vertices(self):
        first, second = self.root / "first.svg", self.root / "second.svg"
        receipt = visualize_mesh(self.inspection, first, highlight_below=.2)
        other = visualize_mesh(self.inspection, second, highlight_below=.2)
        self.assertEqual(first.read_bytes(), second.read_bytes())
        self.assertEqual(receipt["svg_sha256"], hashlib.sha256(first.read_bytes()).hexdigest())
        self.assertEqual(receipt["svg_sha256"], other["svg_sha256"])
        _, data, paths = svg_data(first)
        query = self.inspection.query(metric="mean_ratio", quality_below=.2, associations=False)
        self.assertEqual(data["highlight_count"], 2)
        expected = {row["cell_id"]: row for row in query["cells"]}
        self.assertEqual(set(expected), {1, 2})
        self.assertEqual([int(path.attrib["data-cell-id"]) for path in paths], [1, 2])
        for cell in data["highlights"]:
            self.assertEqual(cell["node_ids"], expected[cell["cell_id"]]["node_ids"])
            self.assertEqual(cell["points"], expected[cell["cell_id"]]["points"])
            self.assertEqual(cell["value"], expected[cell["cell_id"]]["metrics"]["mean_ratio"])
        self.assertEqual(data["mesh_sha256"], self.inspection.mesh_sha256)
        self.assertIsNone(data["source_sha256"])
        self.assertIsNone(data["units"])
        self.assertIn("x-ray", data["highlight_scope"])
        self.assertEqual(receipt["bytes"], first.stat().st_size)
        self.assertTrue(all(path.find("s:title", SVG) is not None for path in paths))

    def test_centroid_cut_uses_whole_original_faces_and_keeps_all_xray_cells(self):
        output = self.root / "cut.svg"
        receipt = visualize_mesh(self.inspection, output, cut_axis="z", cut_value=.003)
        root, data, paths = svg_data(output)
        self.assertEqual(data["cut"]["retained_tet4"], 1)
        self.assertEqual(data["cut"]["total_tet4"], 3)
        self.assertEqual(receipt["cut_surface_triangles"], 1)
        faces = root.findall("s:g/s:polygon[@data-layer='whole-cell-cut']", SVG)
        self.assertEqual(len(faces), 1)
        self.assertEqual(int(faces[0].attrib["data-cell-id"]), 2)
        nodes = [int(value) for value in faces[0].attrib["data-node-ids"].split(",")]
        self.assertEqual(set(nodes), {0, 1, 2})
        self.assertTrue(all(self.inspection._archive.node(index)[0][2] == 0 for index in nodes))
        self.assertEqual([int(path.attrib["data-cell-id"]) for path in paths], [1, 2])

    def test_explicit_ids_and_fixed_world_frame(self):
        frame = [[-1, -1, -1], [5, 2, 2]]
        output = self.root / "selected.svg"
        receipt = visualize_mesh(self.inspection, output, highlight_cell_ids=[3], view_bounds=frame,
                                 azimuth=15, elevation=10, metric="scaled_jacobian")
        _, data, paths = svg_data(output)
        self.assertEqual(receipt["highlighted_cell_ids"], [3])
        self.assertEqual(data["highlight_selection"], "explicit_cell_ids")
        self.assertEqual(data["camera"]["view_bounds"], frame)
        self.assertTrue(data["camera"]["world_frame_explicit"])
        self.assertEqual(len(paths), 1)
        self.assertAlmostEqual(float(paths[0].attrib["data-value"]), 1, places=14)
        self.assertIn("maximum of four incident", data["metric_convention"])
        again = self.root / "selected-again.svg"
        repeated = visualize_mesh(self.inspection, again, highlight_cell_ids=[3], view_bounds=frame,
                                  azimuth=15, elevation=10, metric="scaled_jacobian")
        self.assertEqual(receipt["camera"], repeated["camera"])
        self.assertEqual(first := output.read_bytes(), again.read_bytes())
        self.assertGreater(len(first), 0)

    def test_numeric_legends_include_signed_zero_and_degree_units(self):
        for metric, ticks, units in (("mean_ratio", [0, .25, .5, .75, 1], "dimensionless"),
                                     ("scaled_jacobian", [-1, -.5, 0, .5, 1], "dimensionless"),
                                     ("minimum_dihedral_deg", [0, 15, 30, 45, 60, 70.52877936550931], "degrees")):
            with self.subTest(metric=metric):
                output = self.root / (metric + ".svg")
                receipt = visualize_mesh(self.inspection, output, metric=metric)
                root, data, _ = svg_data(output)
                values = root.findall("s:g[@id='quality-legend']/s:text[@data-metric-tick]", SVG)
                self.assertEqual([float(node.attrib["data-metric-tick"]) for node in values], ticks)
                self.assertEqual(data["legend"]["units"], units)
                self.assertEqual(data["legend"], receipt["legend"])
                if metric == "scaled_jacobian":
                    self.assertNotEqual(data["legend"]["negative_colour"], data["legend"]["zero_colour"])
                    stops = root.findall("s:defs/s:linearGradient/s:stop", SVG)
                    self.assertEqual([node.attrib["offset"] for node in stops], ["0%", "50%", "100%"])
                if metric == "minimum_dihedral_deg":
                    self.assertTrue(all(node.text.endswith("°") for node in values))
        output = self.root / "cell-caption.svg"
        visualize_mesh(self.inspection, output, highlight_cell_ids=[1], view_bounds=[[.1, .1, -.01], [.4, .4, .03]])
        root, data, _ = svg_data(output)
        text = " ".join(root.itertext())
        self.assertIn("Native cell 1", text)
        self.assertIn(format(data["highlights"][0]["value"], ".17g"), text)
        self.assertIn("Close-up viewport; full-object context clipped to frame", text)
        points = [(0, 0, 0), (1, 0, 0), (.5, math.sqrt(3) / 2, 0),
                  (.5, math.sqrt(3) / 6, math.sqrt(2 / 3)),
                  (3, 0, 0), (4, 0, 0), (3, 1, 0), (3.2, .2, 0)]
        source = write_archive(self.root / "signed.paacsmesh", points,
                               [(5, (0, 2, 1, 3), 3, 11), (5, (4, 5, 6, 7), 3, 22)])
        with inspect_mesh(source) as inspection:
            output = self.root / "signed.svg"
            visualize_mesh(inspection, output, metric="scaled_jacobian")
            _, data, paths = svg_data(output)
            self.assertLess(data["highlights"][0]["value"], 0)
            self.assertEqual(data["highlights"][1]["value"], 0)
            self.assertEqual(paths[0].attrib["stroke"], data["legend"]["negative_colour"])
            self.assertEqual(paths[1].attrib["stroke"], data["legend"]["zero_colour"])

    def test_all_150_poor_cells_are_included(self):
        points, cells = [], []
        for index in range(150):
            first = len(points)
            points.extend([(index * 2, 0, 0), (index * 2 + 1, 0, 0), (index * 2, 1, 0), (index * 2 + .3, .3, .01)])
            cells.append((5, tuple(range(first, first + 4)), 3, index + 1))
        source = write_archive(self.root / "many.paacsmesh", points, cells)
        with inspect_mesh(source) as inspection:
            output = self.root / "many.svg"
            receipt = visualize_mesh(inspection, output)
            _, data, paths = svg_data(output)
            self.assertEqual(receipt["highlight_count"], 150)
            self.assertEqual([cell["cell_id"] for cell in data["highlights"]], list(range(150)))
            self.assertEqual([int(path.attrib["data-cell-id"]) for path in paths], list(range(150)))

    def test_budget_refusals_preserve_existing_output_and_cleanup(self):
        output = self.root / "protected.svg"
        output.write_bytes(b"existing output")
        for options, reason in (({"max_highlights": 1}, "highlight_budget"),
                                ({"max_triangles": 1}, "triangle_budget"),
                                ({"max_output_bytes": 100}, "output_budget")):
            with self.subTest(reason=reason):
                with self.assertRaises(InspectionError) as caught:
                    visualize_mesh(self.inspection, output, **options)
                self.assertEqual(caught.exception.reason, reason)
                self.assertEqual(output.read_bytes(), b"existing output")
                self.assertFalse(list(self.root.glob(".paacs-inspection-*")))
        original_work = self.inspection.max_work
        self.inspection.max_work = 1
        with self.assertRaises(InspectionError) as caught:
            visualize_mesh(self.inspection, output)
        self.assertEqual(caught.exception.reason, "work_budget")
        self.inspection.max_work = original_work
        with patch.object(self.inspection, "_reserve", side_effect=InspectionError("memory_budget", "test reservation")):
            with self.assertRaises(InspectionError) as caught:
                visualize_mesh(self.inspection, output)
        self.assertEqual(caught.exception.reason, "memory_budget")
        with patch.object(self.inspection, "_deadline", side_effect=InspectionError("time_budget", "test deadline")):
            with self.assertRaises(InspectionError) as caught:
                visualize_mesh(self.inspection, output)
        self.assertEqual(caught.exception.reason, "time_budget")
        self.assertEqual(output.read_bytes(), b"existing output")

    def test_protects_every_archive_input_path_and_hardlink(self):
        extra = self.root / "report.json"
        extra.write_text("input provenance", encoding="utf-8")
        self.inspection._archive.protected_paths.append(extra)
        for original in (self.source, extra, self.inspection._archive.path):
            before = original.read_bytes()
            with self.assertRaisesRegex(ValueError, "input"):
                visualize_mesh(self.inspection, original)
            self.assertEqual(original.read_bytes(), before)
        linked = self.root / "source-link.svg"
        os.link(self.source, linked)
        with self.assertRaisesRegex(ValueError, "input"):
            visualize_mesh(self.inspection, linked)
        self.assertEqual(linked.read_bytes(), self.source.read_bytes())

    def test_source_mismatch_refuses_before_cad_preview(self):
        cad = Mock(source_sha256="a" * 64)
        with self.assertRaises(InspectionError) as caught:
            visualize_mesh(self.inspection, self.root / "mismatch.svg", cad=cad)
        self.assertEqual(caught.exception.reason, "source_mismatch")
        cad.preview.assert_not_called()
        self.assertFalse((self.root / "mismatch.svg").exists())

    def test_unknown_and_unsupported_metrics_are_not_zero_quality(self):
        with self.assertRaises(ValueError):
            visualize_mesh(self.inspection, self.root / "unknown.svg", metric="aspect_ratio")
        with self.assertRaises(InspectionError) as caught:
            visualize_mesh(self.inspection, self.root / "line.svg", highlight_cell_ids=[0])
        self.assertEqual(caught.exception.reason, "unsupported_metric")
        source = write_archive(self.root / "degenerate.paacsmesh", [(0, 0, 0), (1, 0, 0), (0, 1, 0), (.2, .2, 0)],
                               [(5, (0, 1, 2, 3), 3, 11)])
        with inspect_mesh(source) as inspection:
            output = self.root / "degenerate.svg"
            receipt = visualize_mesh(inspection, output, metric="minimum_dihedral_deg", highlight_below=10)
            _, data, paths = svg_data(output)
            self.assertEqual(receipt["highlight_count"], 0)
            self.assertEqual(data["unavailable_tet4_metric_count"], 1)
            self.assertEqual(paths, [])

    def test_atomic_write_failure_leaves_destination_unchanged(self):
        output = self.root / "unchanged.svg"
        output.write_bytes(b"keep")
        with patch("paacs_mesh.mesh_archive.os.replace", side_effect=OSError("simulated rename failure")):
            with self.assertRaisesRegex(OSError, "rename"):
                visualize_mesh(self.inspection, output)
        self.assertEqual(output.read_bytes(), b"keep")
        self.assertFalse(list(self.root.glob(".paacs-inspection-*")))

    def test_genuine_source_matched_cad_context(self):
        override = os.environ.get("PAACS_TEST_STEP_FIXTURE")
        fixture = Path(override) if override else Path(__file__).resolve().parents[2] / "tests/fixtures/cad/cylinder.step"
        model = pm.read_step(fixture)
        mesh = model.mesh(.5, assembly="conforming", workers=1, validation_workers=1, max_seconds=90)
        with inspect_mesh(mesh) as inspection:
            output = self.root / "cad.svg"
            receipt = visualize_mesh(inspection, output, cad=model)
            _, data, _ = svg_data(output)
            self.assertEqual(data["source_sha256"], model.source_sha256)
            self.assertEqual(data["units"], "mm")
            self.assertEqual(data["context"], "cad_display_tessellation")
            self.assertTrue(data["cad"]["display_geometry_only"])
            self.assertEqual(data["cad"]["linear_deflection"], .02)
            self.assertGreater(receipt["context_triangles"], 0)


if __name__ == "__main__":
    unittest.main()
