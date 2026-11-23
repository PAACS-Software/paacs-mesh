"""Guarded, bounded CAD samples; native tests distinguish trims from planes."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import copy
import json
import math
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch

import paacs_mesh as pm
from paacs_mesh.cad import CadModel
from paacs_mesh.cad_proximity import proximity_scratch_bytes, query_facial_proximity
from paacs_mesh.mesh_inspection import MeshInspection
from test_native import fixture_path


class FakeNative:
    """Receipt/budget fixture only; not evidence of native CAD capability."""
    source_sha256 = "a" * 64
    kernel_version = "test-kernel"
    model_identity = 17
    face_tags, body_tags, edge_tags, vertex_tags = (1, 2, 3), (1,), (), ()
    info = {}

    def __init__(self):
        self.calls = []
        self.fail_call = None
        self.transform = lambda value: value

    def _face_distances(self, points, tags, **options):
        self.calls.append((copy.deepcopy(points), list(tags), dict(options)))
        if len(self.calls) == self.fail_call:
            raise RuntimeError("synthetic kernel failure")
        rows = []
        for point in points:
            faces = [{"face_tag": tag, "distance": abs(point[0] - tag),
                      "closest_point": [float(tag), point[1], point[2]], "tolerance": 1e-7}
                     for tag in tags]
            best = min(faces, key=lambda face: face["distance"])
            row = {"point": list(point), "minimum_distance": best["distance"],
                   "closest_face_tag": best["face_tag"], "closest_point": best["closest_point"],
                   "face_tolerance": best["tolerance"], "measured_face_count": len(tags),
                   "nearest_face_tags": [face["face_tag"] for face in faces
                                         if face["distance"] - best["distance"] <= 1e-7]}
            if options["include_faces"]:
                row["faces"] = faces
            rows.append(row)
        return self.transform({"schema": "paacs-cad-face-distances-v1", "source_sha256": self.source_sha256,
                               "kernel_version": self.kernel_version,
                               "evaluations": len(points) * len(tags), "rows": rows})


class ProximityContractTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        path = Path(self.temporary.name) / "unit.paacsmesh"
        count = 130
        payload = bytearray(struct.pack("<8sIIQQQ", b"PAACSM8\0", 1, 0, 4, count, 4 * count))
        for xyz in ((0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.)):
            payload.extend(struct.pack("<3dii", *xyz, 3, 1))
        payload.extend(struct.pack("<B3xii4q", 5, 3, 1, 0, 1, 2, 3) * count)
        path.write_bytes(payload)
        self.inspection = MeshInspection(path)
        self.addCleanup(self.inspection.close)
        self.native = FakeNative()
        self.model = CadModel(self.native)
        # Deliberately synthetic provenance for helper-contract tests. Native
        # integration below separately exercises real reports and OCCT geometry.
        archive = self.inspection._archive
        archive.report = {"kernel_version": self.native.kernel_version}
        archive.source_sha256, archive.units = self.native.source_sha256, "mm"

    def test_ids_centroid_samples_complete_faces_and_chunk_budgets(self):
        result = query_facial_proximity(self.inspection, list(range(130)), self.model)
        self.assertEqual([row["cell_id"] for row in result["cells"]], list(range(130)))
        self.assertEqual(len(self.native.calls), 3)
        self.assertEqual(result["evaluations"], 390)
        self.assertTrue(result["complete"])
        self.assertEqual(result["mesh_sha256"], self.inspection.mesh_sha256)
        for row in result["cells"]:
            self.assertEqual(row["sample_point"], [.25, .25, .25])
            self.assertEqual(row["nearest_face_tags"], [1])
            self.assertEqual(row["minimum_distance"], .75)
            self.assertNotIn("faces", row)
        for points, tags, options in self.native.calls:
            self.assertLessEqual(len(points), 64)
            self.assertEqual(options["max_work"], len(points) * len(tags))
            self.assertGreater(options["max_seconds"], 0)
        self.assertIn("boundary_ownership", result["does_not_establish"])
        json.dumps(result, allow_nan=False)

    def test_full_face_rows_are_opt_in_and_candidate_order_is_not_renumbered(self):
        result = query_facial_proximity(self.inspection, [129, 7, 0], self.model,
                                       face_tags=[3, 1], include_faces=True)
        self.assertEqual([row["cell_id"] for row in result["cells"]], [129, 7, 0])
        self.assertEqual(result["selected_face_tags"], [1, 3])
        self.assertEqual([face["face_tag"] for face in result["cells"][0]["faces"]], [1, 3])

    def test_empty_candidates_return_empty_complete_evidence_without_native_queries(self):
        result = query_facial_proximity(self.inspection, [], self.model)
        self.assertEqual(result["cells"], [])
        self.assertEqual(result["evaluations"], 0)
        self.assertTrue(result["complete"])
        self.assertFalse(self.native.calls)

    def test_query_distance_filter_keeps_exact_total_before_pagination_and_native_ids(self):
        page = self.inspection.query(cad=self.model, cad_distance=.75, cad_face_tags=[1],
                                     offset=127, limit=2, sort_by="mean_ratio", associations=False)
        self.assertEqual(page["total"], 130)
        self.assertEqual(page["next_offset"], 129)
        self.assertEqual([row["cell_id"] for row in page["cells"]], [127, 128])
        self.assertEqual(page["cad_proximity"]["evaluations"], 130)
        for row in page["cells"]:
            self.assertTrue(row["nearby_cad_surfaces"]["available"])
            self.assertEqual(row["nearby_cad_surfaces"]["cell_id"], row["cell_id"])
            self.assertEqual(row["nearby_cad_surfaces"]["minimum_distance"], .75)
        empty = self.inspection.query(cad=self.model, cad_distance=.74, cad_face_tags=[1], associations=False)
        self.assertEqual(empty["total"], 0)
        self.assertEqual(empty["cells"], [])
        self.native.calls.clear()
        filtered = self.inspection.query(cad=self.model, quality_below=0, associations=False)
        self.assertEqual(filtered["cad_proximity"]["evaluations"], 0)
        self.assertFalse(self.native.calls)

    def test_missing_or_different_source_kernel_and_units_refuse_before_queries(self):
        archive = self.inspection._archive
        for attribute, bad in (("source_sha256", None), ("source_sha256", "b" * 64), ("units", None)):
            original = getattr(archive, attribute)
            setattr(archive, attribute, bad)
            with self.assertRaises(pm.InspectionError):
                query_facial_proximity(self.inspection, [0], self.model)
            setattr(archive, attribute, original)
        for bad in (None, "other-kernel"):
            archive.report["kernel_version"] = bad
            with self.assertRaises(pm.InspectionError):
                query_facial_proximity(self.inspection, [0], self.model)
        self.assertFalse(self.native.calls)

    def test_work_memory_and_enclosing_scratch_admit_before_native_queries(self):
        with self.assertRaises(pm.InspectionError):
            query_facial_proximity(self.inspection, list(range(130)), self.model, max_work=389)
        with self.assertRaises(pm.InspectionError):
            query_facial_proximity(self.inspection, [0], self.model,
                                   _reserved_bytes=self.inspection.max_bytes)
        original_work = self.inspection.max_work
        self.inspection.max_work = 2
        with self.assertRaises(pm.InspectionError):
            query_facial_proximity(self.inspection, [0], self.model)
        self.inspection.max_work = original_work
        previous = self.inspection.max_bytes
        self.inspection.max_bytes = proximity_scratch_bytes(1, 3) - 1
        with self.assertRaises(pm.InspectionError):
            query_facial_proximity(self.inspection, [0], self.model)
        self.inspection.max_bytes = previous
        self.assertFalse(self.native.calls)

    def test_invalid_ids_tags_and_unbounded_generators_are_rejected(self):
        for ids, tags in (([-1], [1]), ([130], [1]), ([0, 0], [1]), ([0], []),
                          ([0], [1, 1]), ([0], [999]), ((i for i in range(2)), [1])):
            with self.assertRaises((ValueError, pm.InspectionError)):
                query_facial_proximity(self.inspection, ids, self.model, face_tags=tags)
        self.assertFalse(self.native.calls)

    def test_second_chunk_failure_and_partial_receipts_never_return_partial_minima(self):
        self.native.fail_call = 2
        with self.assertRaises(pm.InspectionError) as caught:
            query_facial_proximity(self.inspection, list(range(130)), self.model)
        self.assertEqual(caught.exception.reason, "cad_query_failed")
        self.native.fail_call = None
        self.native.transform = lambda receipt: dict(receipt, rows=[])
        with self.assertRaises(pm.InspectionError) as caught:
            query_facial_proximity(self.inspection, [0], self.model)
        self.assertEqual(caught.exception.reason, "cad_query_identity")
        self.native.transform = lambda receipt: dict(receipt, rows=[dict(receipt["rows"][0], faces=[])])
        with self.assertRaises(pm.InspectionError) as caught:
            query_facial_proximity(self.inspection, [0], self.model, include_faces=True)
        self.assertEqual(caught.exception.reason, "cad_query_receipt")

    def test_overdue_native_call_is_distinguished_from_kernel_failure(self):
        clock = [0.]

        def overdue(_):
            clock[0] = 2.
            raise RuntimeError("synthetic cooperative deadline")

        self.native.transform = overdue
        with patch("paacs_mesh.cad_proximity.time.perf_counter", side_effect=lambda: clock[0]):
            with self.assertRaises(pm.InspectionError) as caught:
                query_facial_proximity(self.inspection, [0], self.model, max_seconds=1)
        self.assertEqual(caught.exception.reason, "cad_time_budget")


class NativeTrimmedProximityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(pm._native.CadModel, "_face_distances"):
            raise AssertionError("Installed native module must provide the CAD proximity binding")
        cls.model = pm.read_step(fixture_path())
        cls.mesh = cls.model.mesh(.7, workers=1, validation_workers=1)

    def test_outside_disk_projection_measures_trim_not_infinite_plane(self):
        top = self.model.select_faces(normal=[0, 0, 1]).tags
        face = self.model.inspect(tags=top)["entities"][0]
        center = face["area_centroid"]
        sample = [center[0] + 2, center[1], center[2]]
        result = self.model._native._face_distances([sample], list(top), include_faces=True)
        row = result["rows"][0]
        self.assertAlmostEqual(row["minimum_distance"], 1, places=9)
        self.assertAlmostEqual(math.dist(sample, row["closest_point"]), 1, places=9)
        self.assertAlmostEqual(math.dist(center, row["closest_point"]), 1, places=9)
        self.assertEqual(row["closest_face_tag"], top[0])
        self.assertEqual(row["measured_face_count"], 1)
        self.assertGreater(row["face_tolerance"], 0)

    def test_trimmed_cylinder_end_is_not_an_infinite_cylinder(self):
        side = self.model.select_faces(surface="cylinder").tags
        row = self.model._native._face_distances([[0, 0, 3]], list(side))["rows"][0]
        self.assertAlmostEqual(row["minimum_distance"], math.sqrt(2), places=9)
        self.assertAlmostEqual(row["closest_point"][2], 2, places=9)
        self.assertAlmostEqual(math.hypot(*row["closest_point"][:2]), 1, places=9)

    def test_live_and_reopened_identity_guarded_centroid_evidence(self):
        top = self.model.select_faces(normal=[0, 0, 1]).tags
        with tempfile.TemporaryDirectory() as temporary, pm.inspect_mesh(self.mesh) as inspection:
            ids = [row["cell_id"] for row in inspection.query(limit=3, associations=False)["cells"]]
            result = query_facial_proximity(inspection, ids, self.model, face_tags=top)
            for row in result["cells"]:
                self.assertAlmostEqual(row["minimum_distance"], 2 - row["sample_point"][2], places=8)
            bundle = inspection.save(Path(temporary) / "bundle")
            with pm.open_mesh(bundle) as reopened:
                other = query_facial_proximity(reopened, ids, self.model, face_tags=top)
                self.assertEqual(other["cells"], result["cells"])
                self.assertEqual(other["mesh_sha256"], result["mesh_sha256"])

    def test_native_private_admission_and_unknown_tags_are_refused(self):
        for kwargs in ({"max_work": 1}, {"max_result_bytes": 1}, {"max_work": 100001}):
            with self.assertRaises(ValueError):
                self.model._native._face_distances([[0, 0, 1]], self.model.face_tags, **kwargs)
        with self.assertRaises(ValueError):
            self.model._native._face_distances([[0, 0, 1]], [999999])


if __name__ == "__main__":
    unittest.main()
