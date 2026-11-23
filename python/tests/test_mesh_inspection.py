"""Independent archive/query/field checks; native validity remains a separate gate."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import json
import math
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from unittest.mock import patch

import paacs_mesh as pm
from paacs_mesh.mesh_inspection import distribution
from test_native import fixture_path


def write_archive(path, points=None, cells=None):
    points = points or [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.)]
    if cells is None:
        cells = [(2, 2, 11, (0, 2, 1)), (5, 3, 7, (0, 1, 2, 3)), (5, 3, 0, (0, 2, 1, 3))]
    data = bytearray(struct.pack("<8sIIQQQ", b"PAACSM8\0", 1, 0, len(points), len(cells), sum(len(row[3]) for row in cells)))
    for point in points:
        data.extend(struct.pack("<3dii", *point, 3, 0))
    for kind, dimension, tag, nodes in cells:
        data.extend(struct.pack("<B3sii", kind, b"\0\0\0", dimension, tag))
        data.extend(struct.pack("<" + "q" * len(nodes), *nodes))
    path.write_bytes(data)
    return path


class MeshInspectionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.archive = write_archive(self.root / "mesh.paacsmesh")

    def tearDown(self):
        self.temporary.cleanup()

    def test_native_ids_pagination_lookup_and_unknown_archive_provenance(self):
        with pm.inspect_mesh(self.archive) as inspection:
            first = inspection.query(limit=1)
            second = inspection.query(offset=first["next_offset"], limit=1)
            self.assertEqual([first["cells"][0]["cell_id"], second["cells"][0]["cell_id"]], [1, 2])
            self.assertEqual(first["total"], 2)
            self.assertIsNone(second["next_offset"])
            self.assertEqual(inspection.mesh_sha256, hashlib.sha256(self.archive.read_bytes()).hexdigest())
            self.assertIsNone(inspection.source_sha256)
            self.assertIsNone(inspection.units)
            self.assertEqual(inspection.validation["reported_native_status"], "unavailable")
            boundary = inspection.query(cell_id=0)["cells"][0]
            self.assertEqual(boundary["cell_type"], "tri3")
            self.assertIsNone(boundary["metrics"]["mean_ratio"])
            record = first["cells"][0]
            self.assertEqual(record["node_ids"], [0, 1, 2, 3])
            self.assertEqual(record["centroid"], [.25, .25, .25])
            self.assertEqual(record["body_tag"], 7)
            self.assertEqual(record["cad_boundary_associations"][0]["boundary_cell_id"], 0)
            self.assertEqual(record["cad_boundary_associations"][0]["tag"], 11)
            self.assertFalse(record["nearby_cad_surfaces"]["available"])

    def test_filters_sorting_and_tag_zero(self):
        with pm.open_mesh(self.archive) as inspection:
            inverted = inspection.query(quality_below=.01)
            self.assertEqual([r["cell_id"] for r in inverted["cells"]], [2])
            self.assertEqual(inverted["cells"][0]["orientation"], "inverted")
            self.assertEqual(inverted["cells"][0]["body_tag"], 0)
            self.assertFalse(inverted["cells"][0]["cad_body_identity_known"])
            self.assertEqual(inspection.query(body_tags=[0])["total"], 1)
            self.assertEqual(inspection.query(region=[[.2, .2, .2], [.3, .3, .3]])["total"], 2)
            self.assertEqual(inspection.query(region=[[.8, .8, .8], [1, 1, 1]])["total"], 0)
            self.assertEqual(inspection.query(region=[[.8, .8, .8], [1, 1, 1]], spatial_mode="aabb_overlap")["total"], 2)
            self.assertEqual(inspection.query(sort_by="mean_ratio")["cells"][0]["cell_id"], 2)
            self.assertEqual(inspection.query(sort_by="mean_ratio", descending=True)["cells"][0]["cell_id"], 1)
            with self.assertRaises(ValueError):
                inspection.query(body_tags=iter(range(1000000)))
            for value in (True, -1, 3):
                with self.assertRaises(ValueError):
                    inspection.query(cell_id=value)

    def test_full_bad_population_is_available_not_only_top_twenty(self):
        cells = [(5, 3, 1, (0, 2, 1, 3)) for _ in range(150)]
        write_archive(self.archive, cells=cells)
        with pm.open_mesh(self.archive) as inspection:
            ids = []
            offset = 0
            while offset is not None:
                page = inspection.query(offset=offset, limit=17, quality_below=.2, associations=False)
                ids.extend(row["cell_id"] for row in page["cells"])
                offset = page["next_offset"]
            self.assertEqual(ids, list(range(150)))

    def test_summary_definitions_and_spatial_counts(self):
        with pm.open_mesh(self.archive) as inspection:
            summary = inspection.summary(spatial_bins=2)
            self.assertEqual(summary["orientation_counts"], {"positive": 1, "inverted": 1, "degenerate": 0, "unsupported": 1})
            self.assertEqual(summary["quality"]["mean_ratio"]["count"], 2)
            self.assertEqual(summary["quality"]["mean_ratio"]["minimum"], 0.)
            self.assertAlmostEqual(summary["quality"]["mean_ratio"]["maximum"], 12 * (.5)**(2/3) / 9)
            self.assertEqual(sum(row["count"] for row in summary["spatial_distribution"]["bins"]), 2)
            for metric, global_population in summary["quality"].items():
                populations = [row["quality"][metric] for row in summary["spatial_distribution"]["bins"]]
                self.assertEqual(sum(row["count"] for row in populations), global_population["count"])
                self.assertEqual(sum(row["unavailable"] for row in populations), global_population["unavailable"])
                self.assertEqual([sum(row["histogram"]["counts"][index] for row in populations)
                                  for index in range(len(global_population["histogram"]["counts"]))],
                                 global_population["histogram"]["counts"])
            self.assertEqual(summary["body_tet4_counts"], {"0": 1, "7": 1})
            self.assertEqual(summary["sizing"]["achieved"]["boundary_faces"]["11"]["count"], 3)
            self.assertIsNone(summary["sizing"]["requested"]["global_size"])
            json.dumps(summary, allow_nan=False)

    def test_extreme_finite_statistics_and_centroid_bins(self):
        maximum = sys.float_info.max
        self.assertEqual(distribution([maximum] * 7)["mean"], maximum)
        self.assertEqual(distribution([-maximum, maximum])["median"], 0.)
        tiny = math.ulp(0.)
        write_archive(self.archive, points=[(tiny, 0, 0), (tiny * 2, tiny, 0), (tiny, 0, tiny), (tiny, tiny, tiny)],
                      cells=[(5, 3, 0, (0, 1, 2, 3))])
        with pm.open_mesh(self.archive) as inspection:
            self.assertGreaterEqual(inspection.query(cell_id=0)["cells"][0]["centroid"][0], tiny)
            self.assertTrue(all(0 <= index < 4 for row in inspection.summary()["spatial_distribution"]["bins"] for index in row["index"]))
        write_archive(self.archive, points=[(-maximum, 0, 0), (maximum, 0, 0), (0, maximum, 0), (0, 0, maximum)],
                      cells=[(5, 3, 0, (0, 1, 2, 3))])
        with pm.open_mesh(self.archive) as inspection:
            json.dumps(inspection.summary(), allow_nan=False)

    def test_input_snapshot_is_immutable(self):
        with pm.open_mesh(self.archive) as inspection:
            digest = inspection.mesh_sha256
            self.archive.write_bytes(b"changed")
            self.assertEqual(inspection.mesh_sha256, digest)
            self.assertEqual(inspection.query(cell_id=1)["cells"][0]["node_ids"], [0, 1, 2, 3])

    def test_parser_and_operation_budgets_and_atomic_export(self):
        with self.assertRaises(pm.InspectionError):
            pm.open_mesh(self.archive, max_cells=1)
        with self.assertRaises(pm.InspectionError):
            pm.open_mesh(self.archive, max_seconds=1e-12)
        self.archive.write_bytes(b"short")
        with self.assertRaises(pm.InspectionError):
            pm.open_mesh(self.archive)
        write_archive(self.archive)
        with pm.open_mesh(self.archive, max_bytes=65536) as inspection:
            with self.assertRaises(pm.InspectionError):
                inspection.query(limit=1000)
        with pm.open_mesh(self.archive) as inspection:
            output = self.root / "output.vtk"
            output.write_text("original")
            with self.assertRaises(pm.InspectionError):
                inspection.write_vtk(output, max_output_bytes=20)
            self.assertEqual(output.read_text(), "original")
            with self.assertRaises(ValueError):
                inspection.write_vtk(self.archive)
            inspection.write_vtk(output)
            text = output.read_text()
            self.assertIn("SCALARS native_cell_id long 1", text)
            self.assertIn("SCALARS minimum_dihedral_deg double 1", text)
            self.assertIn("SCALARS scaled_jacobian double 1", text)
            self.assertIn(inspection.mesh_sha256, text)
        with pm.open_mesh(self.archive, max_work=30) as inspection:
            output = self.root / "limited.vtk"
            output.write_text("original")
            with self.assertRaises(pm.InspectionError):
                inspection.write_vtk(output)
            self.assertEqual(output.read_text(), "original")


class NativeInspectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = pm.read_step(fixture_path())
        cls.mesh = cls.model.mesh(.5, workers=1, validation_workers=1)

    def test_native_metrics_and_saved_bundle_hashes(self):
        with tempfile.TemporaryDirectory() as directory, pm.inspect_mesh(self.mesh) as inspection:
            summary = inspection.summary()
            self.assertEqual(inspection.validation["reported_native_status"] in ("valid", "valid_with_warnings"), True)
            self.assertEqual(inspection.units, "mm")
            self.assertEqual(inspection.source_sha256, self.model.source_sha256)
            self.assertAlmostEqual(summary["quality"]["mean_ratio"]["minimum"], self.mesh.report["quality"]["minimum"], places=10)
            self.assertEqual(summary["cell_counts"]["tet4"], self.mesh.tetrahedron_count)
            bundle = inspection.save(Path(directory) / "run")
            with pm.open_mesh(bundle) as reopened:
                self.assertEqual(reopened.mesh_sha256, inspection.mesh_sha256)
                self.assertEqual(reopened.query(limit=3)["cells"], inspection.query(limit=3)["cells"])
                self.assertIn("not_recertified", reopened.validation["basis"])
                for name in ("mesh.paacsmesh", "report.json", "inspection.json"):
                    with self.assertRaises(ValueError):
                        reopened.write_vtk(bundle / name)
            (bundle / "mesh.paacsmesh").write_bytes(b"tampered")
            with self.assertRaises(pm.InspectionError):
                pm.open_mesh(bundle)

    def test_bundle_missing_hash_metadata_bomb_and_duplicate_keys_refused(self):
        with tempfile.TemporaryDirectory() as directory, pm.inspect_mesh(self.mesh) as inspection:
            bundle = inspection.save(Path(directory) / "run")
            manifest_path = bundle / "inspection.json"
            original = json.loads(manifest_path.read_text())
            missing = dict(original)
            missing.pop("mesh_sha256")
            manifest_path.write_text(json.dumps(missing))
            with self.assertRaises(pm.InspectionError):
                pm.open_mesh(bundle)
            report_path = bundle / "report.json"
            report = dict(self.mesh.report, padding=[{} for _ in range(20000)])
            payload = json.dumps(report).encode()
            report_path.write_bytes(payload)
            manifest = dict(original, report_sha256=hashlib.sha256(payload).hexdigest())
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaises(pm.InspectionError):
                pm.open_mesh(bundle, max_bytes=768 * 1024)
            payload = b'{"audit": {}, "audit": {}}'
            report_path.write_bytes(payload)
            manifest["report_sha256"] = hashlib.sha256(payload).hexdigest()
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaises(pm.InspectionError):
                pm.open_mesh(bundle)


if __name__ == "__main__":
    unittest.main()
