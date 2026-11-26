"""Pure archive and saved-receipt integrity checks; no native engine imports."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import importlib.util
import json
from pathlib import Path
import struct
import tempfile
import time
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("initial_generation_comparison", ROOT / "benchmarks/compare_initial_generation.py")
COMPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARE)


def archive(order=(0, 1, 2, 3), *, face_tag=1, vertex_tag=1):
    points = ((0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.))
    ids = {old: new for new, old in enumerate(order)}
    cells = [(5, 3, 1, tuple(ids[i] for i in range(4)))]
    for face in ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)):
        cells.append((2, 2, face_tag, tuple(ids[i] for i in face)))
    result = bytearray(struct.pack("<8sIIQQQ", b"PAACSM8\0", 1, 0, 4, 5, 16))
    for old in order:
        result.extend(struct.pack("<3dii", *points[old], 0, old + vertex_tag))
    for kind, dim, tag, vertices in cells:
        result.extend(struct.pack("<B3sii", kind, b"\0\0\0", dim, tag))
        result.extend(struct.pack("<" + "q" * len(vertices), *vertices))
    return bytes(result)


def limits():
    return SimpleNamespace(max_bytes=1024**2, max_nodes=100, max_cells=100,
                           max_working_bytes=128 * 1024**2, max_seconds=5, max_plot_tets=100)


class InitialGenerationComparisonTests(unittest.TestCase):
    def test_all_public_owned_manifest_fixtures_are_supported(self):
        manifest = json.loads((ROOT / "benchmarks/manifests/initial-generation-owned.json").read_text())
        self.assertEqual(len(manifest["cases"]), 14)
        fixtures = {case["fixture"] for case in manifest["cases"]}
        self.assertEqual(len(fixtures), 8)
        self.assertEqual(COMPARE.OWNED, fixtures)
        for case in manifest["cases"]:
            self.assertEqual(case["kind"], "owned")
            self.assertTrue((ROOT / case["fixture"]).is_file())

    def test_point_order_independent_surface_identity_and_explicit_final_scope(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            a, b = directory / "a.paacsmesh", directory / "b.paacsmesh"
            a.write_bytes(archive()); b.write_bytes(archive((3, 1, 0, 2)))
            first, _ = COMPARE.signatures(a, limits(), time.perf_counter())
            second, _ = COMPARE.signatures(b, limits(), time.perf_counter())
            self.assertEqual(first["tetrahedra"], 1)
            self.assertNotEqual(first["final_all_points"]["ordered_xyz_sha256"], second["final_all_points"]["ordered_xyz_sha256"])
            self.assertEqual(first["final_all_points"]["multiset_xyz_and_classification_sha256"],
                             second["final_all_points"]["multiset_xyz_and_classification_sha256"])
            self.assertEqual(first["surface"]["oriented_coordinates_and_face_tags_sha256"],
                             second["surface"]["oriented_coordinates_and_face_tags_sha256"])
            self.assertFalse(first["initial_point_cloud"]["available"])

    def test_node_and_face_tags_change_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mesh.paacsmesh"
            path.write_bytes(archive()); initial, _ = COMPARE.signatures(path, limits(), time.perf_counter())
            path.write_bytes(archive(vertex_tag=9)); node, _ = COMPARE.signatures(path, limits(), time.perf_counter())
            self.assertEqual(initial["final_all_points"]["ordered_xyz_sha256"], node["final_all_points"]["ordered_xyz_sha256"])
            self.assertNotEqual(initial["final_all_points"]["ordered_xyz_and_classification_sha256"],
                                node["final_all_points"]["ordered_xyz_and_classification_sha256"])
            path.write_bytes(archive(face_tag=8)); face, _ = COMPARE.signatures(path, limits(), time.perf_counter())
            self.assertNotEqual(initial["surface"]["oriented_coordinates_and_face_tags_sha256"],
                                face["surface"]["oriented_coordinates_and_face_tags_sha256"])

    def test_count_and_memory_refusals(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mesh.paacsmesh"; path.write_bytes(archive())
            for key, value in (("max_cells", 4), ("max_bytes", 40), ("max_working_bytes", 1)):
                options = limits(); setattr(options, key, value)
                with self.assertRaises(ValueError):
                    COMPARE.signatures(path, options, time.perf_counter())

    def test_native_report_substitution_is_refused_before_archive_analysis(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            fixture = "tests/fixtures/cad/thin_plate.step"
            source = COMPARE.digest((ROOT / fixture).read_bytes())
            native = {"source_sha256": source, "status": "refused", "reason": "cancelled"}
            payload = json.dumps(native).encode(); (directory / "refusal-report.json").write_bytes(payload)
            case = {"case_id": "unit", "specification": {"fixture": fixture}, "source_sha256": source,
                    "status": "refused", "native_report_file": "refusal-report.json",
                    "native_report_sha256": COMPARE.digest(payload)}
            (directory / "case.json").write_text(json.dumps(case))
            retained = COMPARE.load_case(directory, limits(), time.perf_counter(), {})
            self.assertIsNone(retained["archive"])
            native["reason"] = "geometry_error"
            (directory / "refusal-report.json").write_text(json.dumps(native))
            with self.assertRaisesRegex(ValueError, "Native report hash"):
                COMPARE.load_case(directory, limits(), time.perf_counter(), {})


if __name__ == "__main__":
    unittest.main()
