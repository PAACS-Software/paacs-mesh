"""Installed CLI checks for source-bound mesh inspection, SVG and comparison."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import paacs_mesh as pm
from test_native import fixture_path


SVG = {"s": "http://www.w3.org/2000/svg"}


def write_archive(path):
    # One right tetrahedron and one thin tetrahedron, with stable original IDs.
    points = [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.),
              (3., 0., 0.), (4., 0., 0.), (3., 1., 0.), (3.3, .3, .0001)]
    cells = [(2, 2, 11, (0, 2, 1)), (5, 3, 7, (0, 1, 2, 3)),
             (5, 3, 9, (4, 5, 6, 7)), (0, 1, 21, (0, 1))]
    data = bytearray(struct.pack("<8sIIQQQ", b"PAACSM8\0", 1, 0, len(points), len(cells),
                                 sum(len(row[3]) for row in cells)))
    for point in points:
        data.extend(struct.pack("<3dii", *point, 3, 0))
    for kind, dimension, tag, nodes in cells:
        data.extend(struct.pack("<B3sii", kind, b"\0\0\0", dimension, tag))
        data.extend(struct.pack(f"<{len(nodes)}q", *nodes))
    path.write_bytes(data)
    return path


class MeshCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="paacs-mesh-cli-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name)
        cls.archive = write_archive(cls.root / "mesh.paacsmesh")
        cls.digest = hashlib.sha256(cls.archive.read_bytes()).hexdigest()

    def command(self, *arguments, expected=0, machine=True, timeout=90):
        child = subprocess.run([sys.executable, "-I", "-m", "paacs_mesh", *map(str, arguments)],
                               capture_output=True, text=True, timeout=timeout)
        self.assertEqual(child.returncode, expected, child.stdout + child.stderr)
        self.assertNotIn("Traceback", child.stderr)
        if machine:
            # Parse the entire stream. Native messages must never prefix JSON.
            return json.loads(child.stdout)
        return child

    def test_inspection_json_pagination_original_ids_and_unknown_provenance(self):
        first = self.command("inspect-mesh", self.archive, "--limit", 1, "--json")
        self.assertEqual(first["schema"], "paacs-mesh-cell-page-v1")
        self.assertEqual(first["mesh_sha256"], self.digest)
        self.assertIsNone(first["source_sha256"])
        self.assertIsNone(first["units"])
        self.assertEqual(first["total"], 2)
        self.assertEqual(first["cells"][0]["cell_id"], 1)
        self.assertEqual(first["cells"][0]["node_ids"], [0, 1, 2, 3])
        self.assertEqual(first["cells"][0]["body_tag"], 7)
        second = self.command("inspect-mesh", self.archive, "--offset", first["next_offset"],
                              "--limit", 1, "--json")
        self.assertEqual(second["cells"][0]["cell_id"], 2)
        self.assertIsNone(second["next_offset"])
        boundary = self.command("inspect-mesh", self.archive, "--cell", 0, "--json")
        self.assertEqual(boundary["cells"][0]["cell_type"], "tri3")
        self.assertIsNone(boundary["cells"][0]["metrics"]["mean_ratio"])

    def test_quality_body_region_and_sort_filters(self):
        poor = self.command("inspect-mesh", self.archive, "--quality-below", .2,
                            "--sort-by", "mean_ratio", "--no-associations", "--json")
        self.assertEqual([row["cell_id"] for row in poor["cells"]], [2])
        body = self.command("inspect-mesh", self.archive, "--body", 7, "--json")
        self.assertEqual([row["cell_id"] for row in body["cells"]], [1])
        region = ["--region", .8, 0, 0, .9, .1, .1]
        centroid = self.command("inspect-mesh", self.archive, *region, "--json")
        overlap = self.command("inspect-mesh", self.archive, *region,
                               "--spatial-mode", "aabb_overlap", "--json")
        self.assertEqual(centroid["total"], 0)
        self.assertEqual([row["cell_id"] for row in overlap["cells"]], [1])
        descending = self.command("inspect-mesh", self.archive, "--sort-by", "mean_ratio",
                                  "--descending", "--json")
        self.assertEqual([row["cell_id"] for row in descending["cells"]], [1, 2])

    def test_summary_and_human_output_do_not_claim_bare_archive_certification(self):
        report = self.command("inspect-mesh", self.archive, "--summary", "--spatial-bins", 2, "--json")
        self.assertEqual(report["schema"], "paacs-mesh-quality-summary-v1")
        self.assertEqual(report["scope"], "entire_snapshot")
        self.assertEqual(report["cell_count"], 4)
        self.assertEqual(report["quality"]["mean_ratio"]["count"], 2)
        self.assertEqual(report["validation"]["reported_native_status"], "unavailable")
        self.assertFalse(report["validation"]["independent_native_reaudit_performed"])
        self.assertEqual(report["spatial_distribution"]["bins_per_axis"], 2)
        child = self.command("inspect-mesh", self.archive, "--limit", 1, machine=False)
        self.assertIn("Cell 1 (tet4), body 7", child.stdout)
        self.assertIn("Next page: --offset 1", child.stdout)
        self.assertIn(self.digest, child.stdout)

    def test_summary_rejects_filters_and_proximity_options_instead_of_ignoring_them(self):
        alternatives = [("--cell", 1), ("--quality-below", .2),
                        ("--region", 0, 0, 0, 1, 1, 1), ("--body", 7),
                        ("--orientation", "positive"), ("--cell-type", "all"),
                        ("--source", fixture_path()), ("--cad-distance", .2),
                        ("--cad-face", 1), ("--cad-all-faces",), ("--offset", 1),
                        ("--limit", 10), ("--sort-by", "mean_ratio"), ("--descending",),
                        ("--metric", "scaled_jacobian"), ("--spatial-mode", "aabb_overlap"),
                        ("--no-associations",)]
        for option in alternatives:
            with self.subTest(option=option):
                report = self.command("inspect-mesh", self.archive, "--summary", *option,
                                      "--json", expected=2)
                self.assertEqual(report["schema"], "paacs-mesh-cli-error-v1")
                self.assertEqual(report["reason"], "invalid_argument")
                self.assertIn("entire snapshot", report["message"])

    def test_query_errors_are_single_json_documents_with_stable_exit_codes(self):
        for option in (("--cell", 99), ("--limit", 1001), ("--limit", 0),
                       ("--cad-distance", .1), ("--cad-face", 1), ("--nonsense",)):
            with self.subTest(option=option):
                report = self.command("inspect-mesh", self.archive, *option, "--json", expected=2)
                self.assertEqual(report["status"], "failed")
                self.assertEqual(report["stage"], "input")
                self.assertEqual(report["reason"], "invalid_argument")
        malformed = self.root / "truncated.paacsmesh"
        malformed.write_bytes(b"short")
        report = self.command("inspect-mesh", malformed, "--json", expected=1)
        self.assertEqual(report["stage"], "mesh_inspection")
        self.assertEqual(report["reason"], "archive_header")
        limited = self.command("inspect-mesh", self.archive, "--max-cells", 1, "--json", expected=1)
        self.assertEqual(limited["reason"], "count_budget")

    def test_vtk_export_preserves_original_ids_and_protects_the_input(self):
        output = self.root / "inspected.vtk"
        report = self.command("inspect-mesh", self.archive, "--cell", 2, "--vtk", output, "--json")
        self.assertEqual(report["exports"]["vtk"], str(output))
        vtk = output.read_text(encoding="utf-8")
        self.assertIn("SCALARS native_cell_id long 1\nLOOKUP_TABLE default\n0\n1\n2\n3", vtk)
        self.assertIn("SCALARS mean_ratio double 1", vtk)
        self.assertIn("SCALARS minimum_dihedral_deg double 1", vtk)
        self.assertIn("SCALARS scaled_jacobian double 1", vtk)
        self.assertIn(self.digest, vtk)
        failure = self.command("inspect-mesh", self.archive, "--vtk", self.archive, "--json", expected=2)
        self.assertIn("input", failure["message"])
        self.assertEqual(hashlib.sha256(self.archive.read_bytes()).hexdigest(), self.digest)

    def test_svg_receipt_exact_ids_metric_and_whole_element_cut(self):
        output = self.root / "quality.svg"
        report = self.command("visualize", self.archive, "--output", output, "--cell", 2,
                              "--cut-axis", "x", "--cut-value", 1, "--color-by", "scaled_jacobian",
                              "--width", 640, "--height", 480, "--json")
        self.assertEqual(report["schema"], "paacs-mesh-visualization-receipt-v1")
        self.assertEqual(report["highlighted_cell_ids"], [2])
        self.assertEqual(report["mesh_sha256"], self.digest)
        self.assertEqual(report["bytes"], output.stat().st_size)
        self.assertEqual(report["svg_sha256"], hashlib.sha256(output.read_bytes()).hexdigest())
        self.assertFalse(report["native_validity_rechecked"])
        root = ET.parse(output).getroot()
        metadata = json.loads(root.find("s:metadata", SVG).text)
        self.assertEqual(metadata["mesh_sha256"], self.digest)
        self.assertEqual(metadata["cut"]["retained_tet4"], 1)
        self.assertEqual(metadata["cut"]["total_tet4"], 2)
        self.assertIn("whole Tet4 cells", metadata["cut"]["rule"])
        self.assertEqual(metadata["highlights"][0]["cell_id"], 2)
        self.assertEqual(metadata["highlights"][0]["node_ids"], [4, 5, 6, 7])
        self.assertEqual(metadata["highlight_scope"], report["highlight_scope"])
        paths = root.findall(".//s:g[@id='xray-highlights']/s:path", SVG)
        self.assertEqual([path.get("data-cell-id") for path in paths], ["2"])
        self.assertEqual(paths[0].get("data-node-ids"), "4,5,6,7")
        self.assertEqual(paths[0].get("data-metric"), "scaled_jacobian")
        self.assertAlmostEqual(float(paths[0].get("data-value")), metadata["selected_metric_minimum"])

    def test_svg_budget_and_invalid_selection_preserve_existing_outputs(self):
        output = self.root / "protected.svg"
        output.write_text("original", encoding="utf-8")
        limited = self.command("visualize", self.archive, "--output", output,
                               "--max-output-bytes", 20, "--json", expected=1)
        self.assertEqual(limited["reason"], "output_budget")
        self.assertEqual(output.read_text(encoding="utf-8"), "original")
        for options, code in ((("--cell", 99), 2), (("--cell", 0), 1),
                              (("--cut-value", 0), 2), (("--cut-axis", "z"), 2)):
            with self.subTest(options=options):
                failure = self.command("visualize", self.archive, "--output", output, *options,
                                       "--json", expected=code)
                self.assertEqual(failure["status"], "failed")
                self.assertEqual(output.read_text(encoding="utf-8"), "original")
        self.command("visualize", self.archive, "--output", self.archive, "--json", expected=2)
        self.assertEqual(hashlib.sha256(self.archive.read_bytes()).hexdigest(), self.digest)

    def test_bare_archive_comparison_is_refused_without_inventing_CAD_identity(self):
        report = self.command("compare", self.archive, self.archive, "--json", expected=1)
        self.assertEqual(report["stage"], "mesh_inspection")
        self.assertEqual(report["reason"], "comparison_source")

    def test_many_distinct_boundary_face_distributions_respect_memory_budget(self):
        archive = self.root / "many-boundary-faces.paacsmesh"
        count = 20_000
        with archive.open("wb") as stream:
            stream.write(struct.pack("<8sIIQQQ", b"PAACSM8\0", 1, 0, 3, count, count * 3))
            for point in ((0., 0., 0.), (1., 0., 0.), (0., 1., 0.)):
                stream.write(struct.pack("<3dii", *point, 2, 0))
            for face_tag in range(1, count + 1):
                stream.write(struct.pack("<B3xii3q", 2, 2, face_tag, 0, 1, 2))
        # The archive itself and compact measurement columns fit. Materializing
        # 20,000 per-face distributions and owned JSON copies exceeds 12 MiB.
        failure = self.command("inspect-mesh", archive, "--summary", "--max-bytes", 12 * 1024 ** 2,
                               "--json", expected=1)
        self.assertEqual(failure["stage"], "mesh_inspection")
        self.assertEqual(failure["reason"], "memory_budget")

    @unittest.skipUnless(os.name == "posix" and hasattr(os, "mkfifo"), "POSIX named-pipe fixture")
    def test_named_pipe_bundle_metadata_is_refused_without_waiting_for_a_writer(self):
        for filename in ("inspection.json", "report.json"):
            with self.subTest(filename=filename):
                bundle = self.root / ("fifo-" + filename)
                bundle.mkdir()
                shutil.copyfile(self.archive, bundle / "mesh.paacsmesh")
                if filename == "report.json":
                    (bundle / "inspection.json").write_text(json.dumps({
                        "schema": "paacs-mesh-inspection-bundle-v1", "mesh_sha256": self.digest,
                        "report_sha256": "a" * 64, "source_sha256": None, "units": None,
                    }), encoding="utf-8")
                os.mkfifo(bundle / filename)
                # A blocking metadata open would wait indefinitely; the child
                # deadline makes that regression fail without hanging the suite.
                failure = self.command("inspect-mesh", bundle, "--summary", "--json", expected=1, timeout=5)
                self.assertEqual(failure["reason"], "bundle_metadata")


class NativeMeshComparisonCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="paacs-native-compare-cli-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name)
        cls.model = pm.read_step(fixture_path())
        cls.first, cls.second = cls.root / "coarse", cls.root / "finer"
        for size, destination in ((.7, cls.first), (.5, cls.second)):
            mesh = cls.model.mesh(size, workers=1, validation_workers=1)
            with pm.inspect_mesh(mesh) as inspection:
                inspection.save(destination)

    def command(self, *arguments, expected=0):
        child = subprocess.run([sys.executable, "-I", "-m", "paacs_mesh", *map(str, arguments)],
                               capture_output=True, text=True, timeout=90)
        self.assertEqual(child.returncode, expected, child.stdout + child.stderr)
        self.assertNotIn("Traceback", child.stderr)
        return json.loads(child.stdout)

    def test_native_valid_comparison_outputs_json_and_atomic_file(self):
        output = self.root / "comparison.json"
        comparison = self.command("compare", self.first, self.second, "--output", output,
                                  "--spatial-bins", 2, "--json")
        self.assertEqual(comparison["schema"], "paacs-mesh-comparison-v1")
        self.assertEqual(comparison, json.loads(output.read_text(encoding="utf-8")))
        self.assertEqual(comparison["source_sha256"], self.model.source_sha256)
        self.assertEqual(comparison["units"], "mm")
        for name in ("first", "second"):
            self.assertIn(comparison[name]["validation"]["reported_native_status"], ("valid", "valid_with_warnings"))
            self.assertFalse(comparison[name]["validation"]["independent_native_reaudit_performed"])
            self.assertEqual(set(comparison[name]["quality"]), {"mean_ratio", "minimum_dihedral_deg", "scaled_jacobian"})
            self.assertEqual(comparison[name]["quality"]["mean_ratio"]["unavailable"], 0)
        self.assertEqual(comparison["first"]["sizing"]["requested"]["global_size"], .7)
        self.assertEqual(comparison["second"]["sizing"]["requested"]["global_size"], .5)
        self.assertFalse(comparison["timings"]["end_to_end_runtime_comparison_available"])
        self.assertEqual(comparison["spatial_comparison"]["bins_per_axis"], 2)

    def test_comparison_cannot_replace_either_bound_report_and_refuses_work_limit(self):
        protected = self.second / "report.json"
        original = protected.read_bytes()
        failure = self.command("compare", self.first, self.second, "--output", protected,
                               "--json", expected=2)
        self.assertIn("input", failure["message"])
        self.assertEqual(protected.read_bytes(), original)
        limited = self.command("compare", self.first, self.second, "--max-work", 1, "--json", expected=1)
        self.assertEqual(limited["reason"], "work_budget")

    def test_malformed_hash_paired_sizing_metadata_has_structured_refusal(self):
        malformed = self.root / "malformed-metadata"
        shutil.copytree(self.first, malformed)
        report_path, manifest_path = malformed / "report.json", malformed / "inspection.json"
        original = json.loads(report_path.read_text(encoding="utf-8"))
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for sizing in ([1], {"face_sizes": [1]}, {"volume_controls": [1]}):
            with self.subTest(sizing=sizing):
                report = dict(original, sizing=sizing)
                payload = json.dumps(report, allow_nan=False).encode("utf-8")
                report_path.write_bytes(payload)
                manifest["report_sha256"] = hashlib.sha256(payload).hexdigest()
                manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
                failure = self.command("inspect-mesh", malformed, "--summary", "--json", expected=1)
                self.assertEqual(failure["stage"], "mesh_inspection")
                self.assertEqual(failure["reason"], "inspection_report")


if __name__ == "__main__":
    unittest.main()
