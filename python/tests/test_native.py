"""Acceptance checks for the installed native STEP wheel (stdlib only)."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from concurrent.futures import ThreadPoolExecutor
import gc
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

import paacs_mesh as pm


def fixture_path():
    standalone = os.environ.get("PAACS_TEST_STEP_FIXTURE")
    if standalone:
        return Path(standalone)
    return Path(__file__).resolve().parents[2] / "tests" / "fixtures" / "cad" / "cylinder.step"


def archive_data(path):
    """Read the public archive format independently of the extension."""
    data = path.read_bytes()
    magic, version, flags, nodes, cells, connectivity = struct.unpack_from("<8sIIQQQ", data)
    assert magic == b"PAACSM8\0" and version == 1 and flags == 0
    points, node_tags, cell_tags = [], [], []
    at = 40
    for _ in range(nodes):
        x, y, z, dim, tag = struct.unpack_from("<dddii", data, at)
        points.append((x, y, z))
        node_tags.append((dim, tag))
        at += 32
    widths = (2, 3, 3, 6, 4, 4, 10, 6, 8, 5)
    actual_connectivity = 0
    for _ in range(cells):
        cell_type, dim, tag = struct.unpack_from("<B3xii", data, at)
        cell_tags.append((dim, tag))
        width = widths[cell_type]
        at += 12
        indices = struct.unpack_from(f"<{width}q", data, at)
        assert all(0 <= index < nodes for index in indices)
        at += 8 * width
        actual_connectivity += width
    assert at == len(data) and actual_connectivity == connectivity
    return points, node_tags, cell_tags


class NativeWheelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = fixture_path()
        if not cls.fixture.is_file():
            raise RuntimeError(f"STEP fixture missing: {cls.fixture}")
        cls.model = pm.read_step(cls.fixture)
        cls.mesh = cls.model.mesh(0.3)

    def test_installed_package_when_requested(self):
        if os.environ.get("PAACS_EXPECT_INSTALLED") == "1":
            self.assertIn("site-packages", str(Path(pm.__file__).resolve()).lower())

    def test_cylinder_units_and_native_validity(self):
        info = self.model.info
        self.assertEqual(info["units"], "mm")
        self.assertEqual(info["millimetres_per_model_unit"], 1)
        self.assertEqual(info["body_count"], 1)
        self.assertEqual(info["curved_faces"], 1)
        self.assertAlmostEqual(info["cad_volume"], 2 * math.pi, places=9)
        report = self.mesh.report
        self.assertTrue(report["mesh_published"])
        self.assertEqual(report["reason"], "none")
        self.assertTrue(report["audit"]["complete"])
        self.assertTrue(report["audit"]["valid"])
        self.assertTrue(report["quality"]["complete"])
        self.assertTrue(report["quality"]["valid"])
        self.assertEqual(report["quality"]["inverted_cells"], 0)
        self.assertGreater(report["tetrahedra"], 0)
        self.assertEqual(report["workers"], 1)
        self.assertEqual(report["validation_workers"], 1)
        self.assertEqual(report["mesh_scope"], "independent_parts")
        self.assertFalse(report["interbody_nonoverlap_certified"])
        self.assertEqual(report["bodies"][0]["route"], "cylinder")
        self.assertEqual(self.mesh.face_tags, [1, 2, 3])
        self.assertEqual(self.mesh.body_tags, [1])
        json.dumps(report, allow_nan=False)
        self.assertLess(abs(report["quality"]["total_volume"] / (2 * math.pi) - 1), 0.01)

    def test_unicode_paths_and_tagged_exports(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "PAACS_é_网格"
            root.mkdir()
            source = root / "圆柱.step"
            shutil.copyfile(self.fixture, source)
            mesh = pm.read_step(source).mesh(0.3)
            archive = root / "网格.paacsmesh"
            vtk = root / "网格.vtk"
            mesh.write(archive)
            mesh.write(vtk)
            points, node_tags, cell_tags = archive_data(archive)
            self.assertEqual(len(points), mesh.node_count)
            self.assertEqual(len(cell_tags), mesh.cell_count)
            self.assertEqual(sorted({tag for dim, tag in cell_tags if dim == 3}), mesh.body_tags)
            self.assertEqual(sorted({tag for dim, tag in cell_tags if dim == 2}), mesh.face_tags)
            tokens = vtk.read_text(encoding="ascii").split()
            at = tokens.index("POINTS")
            self.assertEqual(int(tokens[at + 1]), len(points))
            written = [float(x) for x in tokens[at + 3:at + 3 + 3 * len(points)]]
            self.assertEqual(written, [x for point in points for x in point])
            # Both POINT_DATA and CELL_DATA classifications must survive.
            at = 0
            for expected in (cell_tags, node_tags):
                for name, column in (("entity_tag", 1), ("entity_dimension", 0)):
                    at = tokens.index(name, at)
                    start = at + 5  # name int 1 LOOKUP_TABLE default
                    self.assertEqual([int(x) for x in tokens[start:start + len(expected)]],
                                     [tag[column] for tag in expected])
                    at = start + len(expected)

    def test_assembly_preserves_original_bodies(self):
        model = pm.read_step(self.fixture.with_name("assembly.step"))
        self.assertEqual(model.info["body_count"], 2)
        mesh = model.mesh(0.3)
        report = mesh.report
        self.assertEqual(mesh.body_tags, [1, 2])
        self.assertEqual(report["validated_parts"], 2)
        self.assertEqual(sorted(body["solid_tag"] for body in report["bodies"]), [1, 2])
        self.assertTrue(report["audit"]["valid"])
        self.assertTrue(report["quality"]["valid"])
        self.assertEqual(report["mesh_scope"], "independent_parts")
        self.assertFalse(report["interbody_nonoverlap_certified"])
        self.assertAlmostEqual(report["quality"]["total_volume"], 3, places=10)
        self.assertTrue(all(pair["state"] == "no_material_overlap"
                            for pair in report["interferences"]))

    def test_failure_reasons_and_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            missing = Path(directory) / "missing.step"
            with self.assertRaises(pm.CadImportError) as caught:
                pm.read_step(missing)
            self.assertEqual(caught.exception.reason, "file_unavailable")
            self.assertEqual(caught.exception.stage, "import")
            self.assertFalse(caught.exception.report["imported"])
            malformed = Path(directory) / "malformed.step"
            malformed.write_text("this is not STEP", encoding="ascii")
            with self.assertRaises(pm.CadImportError) as caught:
                pm.read_step(malformed)
            self.assertEqual(caught.exception.reason, "parse_failed")
        with self.assertRaises(pm.CadImportError) as caught:
            pm.read_step(self.fixture, max_file_bytes=1)
        self.assertEqual(caught.exception.reason, "file_budget")
        with self.assertRaises(pm.MeshingError) as caught:
            self.model.mesh(0.3, max_nodes=1)
        self.assertFalse(caught.exception.report["mesh_published"])
        self.assertNotEqual(caught.exception.reason, "none")
        self.assertTrue(self.mesh.report["mesh_published"])

    def test_invalid_options_are_rejected(self):
        for size in (0, -1, math.nan, math.inf, True, None, "0.3"):
            with self.subTest(size=size), self.assertRaises((ValueError, TypeError)):
                self.model.mesh(size)
        for options in ({"recovery": "typo"}, {"assembly": "fuse"},
                        {"workers": 0}, {"workers": True}, {"validation_workers": 9},
                        {"max_nodes": -1}, {"max_cells": 2 ** 65}, {"max_work": True},
                        {"max_seconds": math.nan}, {"max_seconds": 1e300}):
            with self.subTest(options=options), self.assertRaises((ValueError, TypeError, OverflowError)):
                self.model.mesh(0.3, **options)
        with self.assertRaises(ValueError):
            pm.read_step("bad\0path.step")

    def test_export_failure_and_source_protection(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.step"
            shutil.copyfile(self.fixture, source)
            original = source.read_bytes()
            mesh = pm.read_step(source).mesh(0.3)
            for writer in (mesh.write_vtk, mesh.write_archive):
                with self.assertRaises(ValueError):
                    writer(source)
                with self.assertRaises(pm.ExportError) as caught:
                    writer(Path(directory) / "absent" / "result.vtk")
                self.assertEqual(caught.exception.stage, "export")
            self.assertEqual(source.read_bytes(), original)
            with self.assertRaises(ValueError):
                mesh.write(Path(directory) / "unsupported.stl")

    def test_mesh_outlives_model_and_input(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "temporary.step"
            shutil.copyfile(self.fixture, source)
            model = pm.read_step(source)
            mesh = model.mesh(0.3)
            del model
            gc.collect()
            source.unlink()
            archive = Path(directory) / "independent.paacsmesh"
            mesh.write_archive(archive)
            self.assertEqual(len(archive_data(archive)[0]), mesh.node_count)

    def test_shared_model_serializes_native_operations(self):
        with ThreadPoolExecutor(max_workers=2) as executor:
            results = list(executor.map(lambda _: self.model.mesh(0.3), range(2)))
        self.assertEqual(results[0].tetrahedron_count, results[1].tetrahedron_count)
        self.assertTrue(all(result.report["audit_valid"] for result in results))

    def test_distinct_models_serialize_internal_worker_dispatch(self):
        models = [pm.read_step(self.fixture), pm.read_step(self.fixture)]
        with ThreadPoolExecutor(max_workers=2) as executor:
            results = list(executor.map(lambda model: model.mesh(0.3, workers=2), models))
        self.assertEqual(results[0].tetrahedron_count, results[1].tetrahedron_count)
        for result in results:
            self.assertEqual(result.report["workers"], 2)
            self.assertTrue(result.report["audit"]["complete"])
            self.assertTrue(result.report["audit"]["valid"])
            self.assertTrue(result.report["quality"]["valid"])

    def test_cli_uses_the_native_engine(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "cli.paacsmesh"
            child = subprocess.run([sys.executable, "-m", "paacs_mesh", str(self.fixture),
                                    str(output), "--size", "0.3"],
                                   capture_output=True, text=True, timeout=60)
            self.assertEqual(child.returncode, 0, child.stderr)
            # OCCT may print transfer diagnostics before the report.
            start = child.stdout.find("{\n")
            report = json.loads(child.stdout[start:])
            self.assertTrue(report["mesh_published"])
            self.assertEqual(len(archive_data(output)[0]), self.mesh.node_count)


if __name__ == "__main__":
    unittest.main()
