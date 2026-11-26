"""Installed API checks for the opt-in, owned-fixture boundary research path."""
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
import hashlib
import json
import os
from pathlib import Path
import struct
import tempfile
import unittest

import paacs_mesh as pm
from paacs_mesh import experimental


class ExperimentalBoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixture = Path(os.environ.get("PAACS_TEST_STEP_FIXTURE", Path(__file__).resolve().parents[2] / "tests/fixtures/cad/cylinder.step"))
        cls.model = pm.read_step(fixture, retain_source=True)
        cls.boundary = experimental.prepare_boundary(cls.model, size=.8, max_seconds=60)

    def assert_qualified(self, mesh):
        report = mesh.report
        self.assertTrue(report["mesh_published"])
        self.assertEqual(report["reason"], "none")
        for section in ("audit", "quality"):
            self.assertTrue(report[section]["complete"])
            self.assertTrue(report[section]["valid"])
        self.assertEqual(set(mesh.body_tags), set(self.model.body_tags))
        self.assertEqual(set(mesh.face_tags), set(self.model.face_tags))
        self.assertEqual(report["experimental_boundary_volume"]["boundary_sha256"], self.boundary.boundary_sha256)
        self.assertTrue(report["experimental_boundary_volume"]["boundary_unchanged"])

    def test_preparation_save_and_fresh_import_precede_any_volume_request(self):
        report = self.boundary.report
        self.assertTrue(report["prepared"])
        self.assertGreater(report["triangles"], 0)
        self.assertNotIn("tetrahedra", report)
        with tempfile.TemporaryDirectory() as temporary:
            path = self.boundary.save(Path(temporary) / "boundary")
            loaded = experimental.load_boundary(path)
        self.assertEqual(loaded.source_sha256, self.model.source_sha256)
        self.assertEqual(loaded.boundary_sha256, self.boundary.boundary_sha256)
        self.assertTrue(loaded.report["rebound"])
        self.assertNotEqual(loaded._model.model_identity, self.model.model_identity)
        self.assert_qualified(loaded.mesh())

    def test_classical_cloud_replays_as_strict_supplied_points(self):
        classical = self.boundary.mesh()
        self.assert_qualified(classical)
        cloud = classical.initial_points(max_bytes=32*1024**2)
        packed = b"".join(struct.pack("<3d", *point) for point in cloud)
        self.assertEqual(hashlib.sha256(packed).hexdigest(), classical.report["experimental_boundary_volume"]["initial_points_sha256"])
        supplied = self.boundary.mesh(cloud)
        self.assert_qualified(supplied)
        report = supplied.report["experimental_boundary_volume"]
        self.assertEqual(report["mode"], "supplied")
        self.assertTrue(report["all_supplied_points_retained"])
        self.assertEqual(report["supplied_points"], len(cloud))
        self.assertEqual(report["retained_points"], len(cloud))
        self.assertEqual(supplied.initial_points(max_bytes=32*1024**2), cloud)

    def test_outside_supplied_point_refuses_without_replacement(self):
        before = self.boundary.boundary_sha256
        with self.assertRaises(pm.MeshingError) as refusal:
            self.boundary.mesh([(1e9, 1e9, 1e9)])
        report = refusal.exception.report
        self.assertFalse(report["mesh_published"])
        self.assertEqual(report["experimental_boundary_volume"]["mode"], "supplied")
        self.assertEqual(report["experimental_boundary_volume"]["supplied_points"], 1)
        self.assertEqual(report["experimental_boundary_volume"]["supplied_input_sha256"], hashlib.sha256(struct.pack("<3d",1e9,1e9,1e9)).hexdigest())
        self.assertEqual(self.boundary.boundary_sha256, before)

    def test_binding_replays_effective_chord_and_preserves_requested_metadata(self):
        # Exercise the actual C++ admission rule rather than mocking replay:
        # a tightened snapshot must not be passed the original coarser request.
        options = json.loads(self.boundary._options)
        requested = self.boundary.report["effective_chord_error"] * 2
        options["max_chord_error"] = requested
        mesh = pm._native._experimental_mesh_boundary(self.model._native, self.boundary._boundary, options, None)
        self.assert_qualified(mesh)
        self.assertEqual(mesh.report["fidelity"]["requested_max_chord_error"], requested)

    def test_rehashed_archive_still_requires_original_cad_rebind(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = self.boundary.save(Path(temporary) / "boundary")
            archive = bytearray((path / "boundary.paacsmesh").read_bytes())
            struct.pack_into("<d", archive, 40, struct.unpack_from("<d", archive, 40)[0] + 1000.)
            (path / "boundary.paacsmesh").write_bytes(archive)
            manifest = json.loads((path / "manifest.json").read_text())
            digest = hashlib.sha256(archive).hexdigest()
            manifest["boundary_sha256"] = manifest["files"]["boundary.paacsmesh"]["sha256"] = digest
            (path / "manifest.json").write_text(json.dumps(manifest))
            with self.assertRaises(pm.MeshingError) as refusal:
                experimental.load_boundary(path)
            self.assertEqual(refusal.exception.stage, "boundary_rebind")


if __name__ == "__main__":
    unittest.main()
