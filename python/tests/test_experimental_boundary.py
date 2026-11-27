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
from paacs_mesh import pipeline
from sizing_data import read_archive


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

    def test_public_facade_and_standalone_generation_precede_volume(self):
        self.assertIs(pipeline.PreparedBoundary, experimental.PreparedBoundary)
        self.assertIs(pipeline.prepare_boundary, experimental.prepare_boundary)
        before = self.boundary.boundary_sha256
        cloud = self.boundary.generate_initial_points()
        self.assertIsInstance(cloud, pipeline.InitialPointCloud)
        report = cloud.report
        self.assertTrue(report["generated"])
        self.assertTrue(report["initial_points_available"])
        self.assertTrue(report["boundary_unchanged"])
        self.assertFalse(report["recovery_performed"])
        self.assertFalse(report["volume_construction_performed"])
        self.assertFalse(report.get("mesh_published", False))
        self.assertNotIn("audit", report)
        self.assertEqual(cloud.source_sha256, self.boundary.source_sha256)
        self.assertEqual(cloud.boundary_sha256, before)
        self.assertEqual(cloud.sha256, report["initial_points_sha256"])
        points = cloud.points(max_bytes=32*1024**2)
        self.assertEqual(len(points), cloud.point_count)
        classical = self.boundary.mesh()
        self.assertEqual(points, classical.initial_points(max_bytes=32*1024**2))
        supplied = self.boundary.mesh(cloud)
        self.assert_qualified(supplied)
        self.assertEqual(supplied.report["experimental_boundary_volume"]["supplied_input_sha256"], cloud.sha256)
        self.assertEqual(self.boundary.boundary_sha256, before)

    def test_generated_point_bundle_roundtrip_and_immutable_copies(self):
        cloud = self.boundary.generate_initial_points()
        with tempfile.TemporaryDirectory() as temporary:
            path = cloud.save(Path(temporary)/"points")
            loaded = pipeline.load_initial_points(path)
            self.assertEqual(loaded.sha256, cloud.sha256)
            self.assertEqual(loaded.points(max_bytes=32*1024**2), cloud.points(max_bytes=32*1024**2))
            self.assert_qualified(self.boundary.mesh(loaded))
            with self.assertRaises(FileExistsError):cloud.save(path)
        with self.assertRaises(AttributeError):cloud._payload = b""
        copy = cloud.report;copy["generated"] = False
        self.assertTrue(cloud.report["generated"])
        with self.assertRaises(ValueError):cloud.points(max_bytes=1)

    def test_cloud_identity_and_malformed_requests_refuse_without_replacement(self):
        cloud = self.boundary.generate_initial_points()
        other = pipeline.prepare_boundary(self.model, size=.9, max_seconds=60)
        with self.assertRaisesRegex(ValueError, "source/boundary/build/options"):
            other.mesh(cloud)
        for points in ([(True, 0., 0.)], [(float("nan"), 0., 0.)], [(0., 0.)]):
            with self.assertRaises(ValueError):self.boundary.mesh(points)
        point = cloud.points(max_bytes=32*1024**2)
        if point:
            with self.assertRaises(pm.MeshingError) as refusal:
                self.boundary.mesh([point[0], point[0]])
            self.assertFalse(refusal.exception.report["mesh_published"])

    def test_box_and_through_hole_stage_replay_preserve_boundary_and_tags(self):
        fixtures = Path(__file__).resolve().parents[2]/"tests/fixtures/cad"
        for name in ("box.step", "through_hole.step"):
            with self.subTest(fixture=name):
                model = pm.read_step(fixtures/name, retain_source=True)
                boundary = pipeline.prepare_boundary(model, size=.4, max_seconds=60)
                cloud = boundary.generate_initial_points()
                mesh = boundary.mesh(cloud)
                report = mesh.report
                self.assertTrue(report["mesh_published"])
                self.assertTrue(report["audit"]["valid"])
                self.assertTrue(report["quality"]["valid"])
                self.assertTrue(report["experimental_boundary_volume"]["boundary_unchanged"])
                self.assertEqual(report["experimental_boundary_volume"]["boundary_sha256"], boundary.boundary_sha256)
                self.assertEqual(set(mesh.face_tags), set(model.face_tags))
                self.assertEqual(set(mesh.body_tags), set(model.body_tags))
                self.assertEqual(report["experimental_boundary_volume"]["admitted_points"], cloud.point_count)

    def test_unsupported_assembly_and_controls_do_not_enter_pipeline(self):
        fixture = Path(__file__).resolve().parents[2]/"tests/fixtures/cad/assembly.step"
        assembly = pm.read_step(fixture, retain_source=True)
        with self.assertRaisesRegex(ValueError, "one original solid"):
            pipeline.prepare_boundary(assembly, size=.8)
        for options in ({"face_sizes": []}, {"growth": 1.3}, {"workers": 2}, {"max_nodes": True}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                pipeline.prepare_boundary(self.model, size=.8, **options)

    def test_real_box_and_cylinder_surface_edge_exports_keep_native_ids_and_metadata(self):
        fixture = Path(__file__).resolve().parents[2]/"tests/fixtures/cad/box.step"
        for name, model in (("cylinder", self.model), ("box", pm.read_step(fixture, retain_source=True))):
            with self.subTest(fixture=name), tempfile.TemporaryDirectory() as temporary:
                boundary = pipeline.prepare_boundary(model, size=.8, max_seconds=60)
                surface, edges = boundary.surface(), boundary.edges()
                directory = Path(temporary)
                surface_path = surface.save(directory/"surface")
                edge_path = edges.save(directory/"edges")
                for kind, path, artifact in (("surface", surface_path, surface), ("cad_edges", edge_path, edges)):
                    manifest = json.loads((path/"manifest.json").read_text())
                    self.assertEqual(manifest["kind"], kind)
                    self.assertFalse(manifest["volume_qualified"])
                    self.assertEqual(manifest["source_sha256"], model.source_sha256)
                    self.assertEqual(manifest["boundary_sha256"], boundary.boundary_sha256)
                    for member, receipt in manifest["files"].items():
                        data = (path/member).read_bytes()
                        self.assertEqual(receipt["bytes"], len(data))
                        self.assertEqual(receipt["sha256"], hashlib.sha256(data).hexdigest())
                    metadata = json.loads((path/"metadata.json").read_text())
                    self.assertEqual(metadata["cad_topology"], boundary.cad_topology())
                    self.assertEqual(artifact.source_sha256, model.source_sha256)
                    self.assertTrue((path/"mesh.vtk").read_bytes().startswith(b"# vtk"))
                surface_archive = read_archive(surface_path/"mesh.paacsmesh")
                edge_archive = read_archive(edge_path/"mesh.paacsmesh")
                self.assertTrue(surface_archive.cells);self.assertTrue(edge_archive.cells)
                self.assertEqual({cell.kind for cell in surface_archive.cells}, {2})
                self.assertEqual({cell.kind for cell in edge_archive.cells}, {0})
                self.assertEqual(surface_archive.points, edge_archive.points)
                self.assertEqual(surface_archive.node_tags, edge_archive.node_tags)
                self.assertEqual({cell.tag for cell in surface_archive.cells}, set(model.face_tags))
                expected_edges = {edge["edge_tag"] for edge in edges.report["edges"] if edge["natural_line_order"]}
                self.assertEqual({cell.tag for cell in edge_archive.cells}, expected_edges)
                self.assertEqual(hashlib.sha256((surface_path/"mesh.paacsmesh").read_bytes()).hexdigest(), boundary.boundary_sha256)
                self.assertEqual(edges.report["node_id_space"], "qualified_boundary")


if __name__ == "__main__":
    unittest.main()
