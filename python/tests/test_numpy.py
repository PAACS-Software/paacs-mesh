"""Installed-wheel checks for independent NumPy mesh snapshots."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import gc
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import paacs_mesh as pm
from paacs_mesh import pipeline

from sizing_data import read_archive
from test_native import fixture_path

try:
    import numpy as np
except ImportError:
    if os.environ.get("PAACS_REQUIRE_NUMPY") == "1":
        raise RuntimeError("NumPy qualification was requested but NumPy is unavailable")
    np = None


CELL_NAMES = (
    "line2", "line3", "tri3", "tri6", "quad4", "tet4", "tet10", "prism6", "hex8", "pyramid5",
)
CELL_WIDTHS = (2, 3, 3, 6, 4, 4, 10, 6, 8, 5)


def numerical_arrays(snapshot):
    """Iterate every numeric payload included in the extraction budget."""
    yield snapshot["points"]
    yield snapshot["node_entity_dimensions"]
    yield snapshot["node_entity_tags"]
    for block in snapshot["cells"].values():
        yield block["connectivity"]
        yield block["entity_dimensions"]
        yield block["entity_tags"]
        yield block["cell_indices"]


class OptionalNumPyTests(unittest.TestCase):
    def test_base_environment_without_numpy_when_requested(self):
        if os.environ.get("PAACS_EXPECT_NO_NUMPY") == "1":
            self.assertIsNone(np, "The stdlib-only qualification environment contains NumPy")

    def test_import_meshing_and_export_without_numpy(self):
        script = r'''
import importlib.abc
from pathlib import Path
import sys

class NoNumPy(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "numpy" or fullname.startswith("numpy."):
            raise ModuleNotFoundError("NumPy intentionally unavailable", name=fullname)

sys.meta_path.insert(0, NoNumPy())
import paacs_mesh as pm
from paacs_mesh import pipeline
model = pm.read_step(Path(sys.argv[1]))
mesh = model.mesh(.3)
assert mesh.report["audit"]["valid"]
mesh.write_archive(Path(sys.argv[2]))
assert Path(sys.argv[2]).is_file()
assert not any(name == "numpy" or name.startswith("numpy.") for name in sys.modules)
try:
    mesh.to_numpy()
except ImportError as error:
    message = str(error).lower()
    assert "numpy" in message and "install" in message, str(error)
else:
    raise AssertionError("Array extraction succeeded with NumPy unavailable")
boundary = pipeline.prepare_boundary(model, size=.8, max_seconds=60)
cloud = boundary.generate_initial_points()
assert cloud.report["generated"] and not cloud.report["volume_construction_performed"]
assert isinstance(cloud.points(), tuple)
assert boundary.mesh(cloud).report["audit"]["valid"]
try:
    cloud.to_numpy()
except ImportError as error:
    assert "numpy" in str(error).lower() and "install" in str(error).lower()
else:
    raise AssertionError("Point arrays succeeded with NumPy unavailable")
assert not any(name == "numpy" or name.startswith("numpy.") for name in sys.modules)
'''
        with tempfile.TemporaryDirectory(prefix="paacs-without-numpy-") as directory:
            archive = Path(directory) / "without-numpy.paacsmesh"
            child = subprocess.run(
                [sys.executable, "-I", "-c", script, str(fixture_path()), str(archive)],
                capture_output=True, text=True, timeout=60,
            )
            self.assertEqual(child.returncode, 0, child.stdout + child.stderr)
            self.assertGreater(len(read_archive(archive).cells), 0)


@unittest.skipIf(np is None, "NumPy extra is not installed")
class NumPySnapshotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = fixture_path()
        cls.model = pm.read_step(cls.fixture)
        cls.mesh = cls.model.mesh(.3)
        cls.temporary = tempfile.TemporaryDirectory(prefix="paacs-numpy-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.directory = Path(cls.temporary.name)
        cls.archive_path = cls.directory / "original.paacsmesh"
        cls.mesh.write_archive(cls.archive_path)
        cls.archive = read_archive(cls.archive_path)
        cls.snapshot = cls.mesh.to_numpy()

    def assert_snapshot_matches_archive(self, snapshot, archive):
        self.assertEqual(snapshot["units"], "mm")
        points = snapshot["points"]
        self.assertEqual(points.dtype, np.dtype("float64"))
        self.assertEqual(points.shape, (len(archive.points), 3))
        np.testing.assert_array_equal(points, archive.points)
        for key, column in (("node_entity_dimensions", 0), ("node_entity_tags", 1)):
            array = snapshot[key]
            self.assertEqual(array.dtype, np.dtype("int32"))
            self.assertEqual(array.shape, (len(archive.points),))
            np.testing.assert_array_equal(array, [tag[column] for tag in archive.node_tags])

        expected_kinds = {cell.kind for cell in archive.cells}
        self.assertEqual(set(snapshot["cells"]), {CELL_NAMES[kind] for kind in expected_kinds})
        seen = []
        for kind in expected_kinds:
            expected = [(index, cell) for index, cell in enumerate(archive.cells) if cell.kind == kind]
            block = snapshot["cells"][CELL_NAMES[kind]]
            connectivity = block["connectivity"]
            self.assertEqual(connectivity.dtype, np.dtype("int64"))
            self.assertEqual(connectivity.shape, (len(expected), CELL_WIDTHS[kind]))
            np.testing.assert_array_equal(connectivity, [cell.nodes for _, cell in expected])
            self.assertTrue(np.all(connectivity >= 0))
            self.assertTrue(np.all(connectivity < len(archive.points)))
            for key, dtype, values in (
                ("entity_dimensions", "int32", [cell.dimension for _, cell in expected]),
                ("entity_tags", "int32", [cell.tag for _, cell in expected]),
                ("cell_indices", "int64", [index for index, _ in expected]),
            ):
                array = block[key]
                self.assertEqual(array.dtype, np.dtype(dtype))
                self.assertEqual(array.shape, (len(expected),))
                np.testing.assert_array_equal(array, values)
            seen.extend(block["cell_indices"].tolist())
        self.assertEqual(sorted(seen), list(range(len(archive.cells))))
        for array in numerical_arrays(snapshot):
            self.assertTrue(array.flags.c_contiguous)
            self.assertTrue(array.flags.owndata)
            self.assertFalse(array.flags.writeable)

    def test_archive_parity_and_array_contract(self):
        self.assert_snapshot_matches_archive(self.snapshot, self.archive)
        self.assertEqual(len(self.snapshot["points"]), self.mesh.node_count)
        self.assertEqual(len(self.snapshot["cells"]["tet4"]["connectivity"]),
                         self.mesh.tetrahedron_count)

    def test_tet_orientation_volume_and_cad_groups(self):
        block = self.snapshot["cells"]["tet4"]
        vertices = self.snapshot["points"][block["connectivity"]]
        edges = vertices[:, 1:] - vertices[:, :1]
        signed_volumes = np.einsum("ij,ij->i", edges[:, 0], np.cross(edges[:, 1], edges[:, 2])) / 6
        self.assertTrue(np.all(np.isfinite(signed_volumes)))
        self.assertTrue(np.all(signed_volumes > 0))
        self.assertAlmostEqual(float(signed_volumes.sum()),
                               self.mesh.report["quality"]["total_volume"], places=11)
        np.testing.assert_array_equal(np.unique(block["entity_dimensions"]), [3])
        np.testing.assert_array_equal(np.unique(block["entity_tags"]), self.mesh.body_tags)
        triangles = self.snapshot["cells"]["tri3"]
        np.testing.assert_array_equal(np.unique(triangles["entity_dimensions"]), [2])
        np.testing.assert_array_equal(np.unique(triangles["entity_tags"]), self.mesh.face_tags)

    def test_snapshot_edits_cannot_change_native_mesh(self):
        snapshot = self.mesh.to_numpy()
        original = self.archive_path.read_bytes()
        for edited, baseline in zip(numerical_arrays(snapshot), numerical_arrays(self.snapshot)):
            self.assertFalse(np.shares_memory(edited, baseline))
            with self.assertRaises(ValueError):
                edited.flat[0] = 0
            edited.setflags(write=True)
            edited.fill(-123)
        saved = self.directory / "after-snapshot-edits.paacsmesh"
        self.mesh.write_archive(saved)
        self.assertEqual(saved.read_bytes(), original)
        self.assert_snapshot_matches_archive(self.mesh.to_numpy(), self.archive)
        self.assertTrue(self.mesh.report["audit"]["valid"])

    def test_snapshot_and_slices_outlive_mesh_model_and_source(self):
        with tempfile.TemporaryDirectory(prefix="paacs-numpy-lifetime-") as directory:
            source = Path(directory) / "temporary.step"
            shutil.copyfile(self.fixture, source)
            model = pm.read_step(source)
            mesh = model.mesh(.3)
            snapshot = mesh.to_numpy()
            points = snapshot["points"]
            connectivity = snapshot["cells"]["tet4"]["connectivity"][::2]
            point_copy, cell_copy = points.copy(), connectivity.copy()
            del snapshot, mesh, model
            source.unlink()
            gc.collect()
            np.testing.assert_array_equal(points, point_copy)
            np.testing.assert_array_equal(connectivity, cell_copy)
            self.assertTrue(np.all(np.isfinite(points[connectivity])))

    def test_assembly_body_ids_and_node_identity_survive(self):
        mesh = pm.read_step(self.fixture.with_name("assembly.step")).mesh(.3)
        snapshot = mesh.to_numpy()
        archive_path = self.directory / "assembly.paacsmesh"
        mesh.write_archive(archive_path)
        self.assert_snapshot_matches_archive(snapshot, read_archive(archive_path))
        self.assertEqual(snapshot["points"].shape[0], mesh.node_count)
        block = snapshot["cells"]["tet4"]
        np.testing.assert_array_equal(np.unique(block["entity_tags"]), [1, 2])
        first = np.unique(block["connectivity"][block["entity_tags"] == 1])
        second = np.unique(block["connectivity"][block["entity_tags"] == 2])
        self.assertEqual(np.intersect1d(first, second).size, 0)
        self.assertEqual(mesh.report["mesh_scope"], "independent_parts")
        self.assertFalse(mesh.report["interbody_nonoverlap_certified"])

    def test_exact_payload_budget_and_refusal_leave_mesh_usable(self):
        required = sum(array.nbytes for array in numerical_arrays(self.snapshot))
        with self.assertRaises(MemoryError):
            self.mesh.to_numpy(max_bytes=required - 1)
        with self.assertRaises(MemoryError):
            self.mesh.to_numpy(max_bytes=1)
        self.assert_snapshot_matches_archive(self.mesh.to_numpy(max_bytes=required), self.archive)
        saved = self.directory / "after-budget-refusal.paacsmesh"
        self.mesh.write_archive(saved)
        self.assertEqual(saved.read_bytes(), self.archive_path.read_bytes())

    def test_invalid_budgets_are_rejected(self):
        for budget in (False, True, 0, -1, 1.5, None, "1024", 2 ** 65):
            with self.subTest(max_bytes=budget), self.assertRaises((TypeError, ValueError, OverflowError)):
                self.mesh.to_numpy(max_bytes=budget)


@unittest.skipIf(np is None, "NumPy extra is not installed")
class PipelineNumPyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = pm.read_step(fixture_path(), retain_source=True)
        cls.boundary = pipeline.prepare_boundary(cls.model, size=.8, max_seconds=60)
        cls.cloud = cls.boundary.generate_initial_points()
        cls.temporary = tempfile.TemporaryDirectory(prefix="paacs-pipeline-numpy-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.bundle = cls.boundary.save(Path(cls.temporary.name)/"boundary")
        cls.archive = read_archive(cls.bundle/"boundary.paacsmesh")

    def test_boundary_snapshot_is_exact_tagged_oriented_tri3_and_owned(self):
        snapshot = self.boundary.to_numpy()
        NumPySnapshotTests.assert_snapshot_matches_archive(self, snapshot, self.archive)
        self.assertEqual(set(snapshot["cells"]), {"tri3"})
        tri = snapshot["cells"]["tri3"]
        self.assertEqual(set(tri["entity_tags"].tolist()), set(self.model.face_tags))
        vertices = snapshot["points"][tri["connectivity"]]
        areas = np.linalg.norm(np.cross(vertices[:, 1]-vertices[:, 0], vertices[:, 2]-vertices[:, 0]), axis=1)/2
        self.assertTrue(np.all(np.isfinite(areas)));self.assertTrue(np.all(areas > 0))
        volume = np.einsum("ij,ij->i", vertices[:, 0], np.cross(vertices[:, 1], vertices[:, 2])).sum()/6
        self.assertGreater(volume, 0)
        self.assertLess(abs(float(volume)/self.model.info["cad_volume"]-1), .01)

    def test_boundary_copy_edits_and_copy_budget_cannot_change_frozen_native_data(self):
        before = (self.bundle/"boundary.paacsmesh").read_bytes();snapshot = self.boundary.to_numpy()
        for array in numerical_arrays(snapshot):
            with self.assertRaises(ValueError):array.flat[0] = 0
            array.setflags(write=True);array.fill(-999)
        NumPySnapshotTests.assert_snapshot_matches_archive(self, self.boundary.to_numpy(), self.archive)
        with self.assertRaises(MemoryError):self.boundary.to_numpy(max_bytes=1)
        path = self.boundary.save(Path(self.temporary.name)/"boundary-after-copy-edits")
        self.assertEqual((path/"boundary.paacsmesh").read_bytes(), before)

    def test_cloud_arrays_are_bounded_independent_readonly_float64_copies(self):
        cloud = self.cloud;array = cloud.to_numpy();other = cloud.to_numpy();digest = cloud.sha256
        self.assertEqual(array.shape, (cloud.point_count, 3));self.assertEqual(array.dtype, np.dtype("float64"))
        self.assertTrue(array.flags.owndata);self.assertTrue(array.flags.c_contiguous)
        self.assertFalse(array.flags.writeable);self.assertFalse(np.shares_memory(array, other))
        np.testing.assert_array_equal(array, np.asarray(cloud.points(max_bytes=32*1024**2)).reshape((-1, 3)))
        if array.size:
            with self.assertRaises(ValueError):array.flat[0] = 0
            array.setflags(write=True);array.fill(-999)
        self.assertEqual(cloud.sha256, digest)
        np.testing.assert_array_equal(cloud.to_numpy(), other)
        with self.assertRaises(ValueError):cloud.to_numpy(max_bytes=1)
        with self.assertRaises(ValueError):cloud.to_numpy(max_bytes=True)

    def test_generated_cloud_array_is_accepted_as_explicit_points_without_mutation(self):
        array = self.cloud.to_numpy();before = array.tobytes()
        mesh = self.boundary.mesh(array)
        self.assertTrue(mesh.report["mesh_published"])
        self.assertEqual(mesh.report["experimental_boundary_volume"]["supplied_input_sha256"], self.cloud.sha256)
        self.assertEqual(array.tobytes(), before)


if __name__ == "__main__":
    unittest.main()
