"""Installed-wheel checks for native CAD face, volume and growth controls."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import paacs_mesh as pm

from sizing_data import measure, read_archive
from test_native import fixture_path


LOWER = (-.3, -.3, .7)
UPPER = (.3, .3, 1.3)
VARIANTS = ("global", "face", "box", "combined")


class NativeSizingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = fixture_path()
        cls.model = pm.read_step(cls.fixture)
        cls.temporary = tempfile.TemporaryDirectory(prefix="paacs-python-sizing-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.directory = Path(cls.temporary.name)
        cls.meshes, cls.archives, cls.measurements = {}, {}, {}
        for variant in VARIANTS:
            mesh = cls.model.mesh(.7, **cls.options(variant))
            archive = cls.directory / f"{variant}.paacsmesh"
            mesh.write_archive(archive)
            cls.meshes[variant] = mesh
            cls.archives[variant] = archive.read_bytes()
            cls.measurements[variant] = measure(read_archive(archive), 1, LOWER, UPPER, .7, .15)

    @staticmethod
    def options(variant="combined", **overrides):
        values = dict(
            face_sizes=[pm.FaceSize(1, .16)] if variant in ("face", "combined") else [],
            volume_controls=[pm.BoxSize(LOWER, UPPER, .15)] if variant in ("box", "combined") else [],
            min_size=.12, growth=1.55, native_sizing=True,
            assembly="conforming", recovery="standard", max_work=100_000_000,
            workers=1, validation_workers=1,
        )
        values.update(overrides)
        return values

    @classmethod
    def tearDownClass(cls):
        destination = os.environ.get("PAACS_SIZING_RESULTS")
        if destination:
            path = Path(destination)
            path.parent.mkdir(parents=True, exist_ok=True)
            payload = {
                "schema": "paacs-python-sizing-checks-v1",
                "source_sha256": hashlib.sha256(cls.fixture.read_bytes()).hexdigest(),
                "sizes_are_best_effort": True,
                "measurements": cls.measurements,
                "requested_controls": cls.meshes["combined"].report["sizing"],
            }
            path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    def assert_native_result(self, mesh):
        report = mesh.report
        self.assertTrue(report["mesh_published"])
        self.assertEqual(report["reason"], "none")
        self.assertTrue(report["audit"]["complete"])
        self.assertTrue(report["audit"]["valid"])
        self.assertTrue(report["quality"]["complete"])
        self.assertTrue(report["quality"]["valid"])
        self.assertEqual(report["quality"]["inverted_cells"], 0)
        self.assertEqual(mesh.face_tags, [1, 2, 3])
        self.assertEqual(mesh.body_tags, [1])
        self.assertEqual(report["mesh_scope"], "conforming_domain")
        self.assertEqual(len(report["bodies"]), 1)
        body = report["bodies"][0]
        self.assertEqual(body["route"], "native_spatial_sizing")
        self.assertEqual(body["reason"], "none")
        self.assertTrue(body["boundary_embedding_certified"])
        self.assertEqual(body["native_centroids_checked"], mesh.tetrahedron_count)
        self.assertEqual(body["centroid_outside_distance_limit"], 0)
        self.assertLessEqual(body["maximum_boundary_sample_error"], body["surface_chord_error_limit"] + 1e-12)
        self.assertLessEqual(body["maximum_boundary_node_error"], body["surface_chord_error_limit"] + 1e-12)
        self.assertLessEqual(body["relative_volume_error"], .01)
        json.dumps(report, allow_nan=False)

    def saved(self, mesh, name):
        path = self.directory / f"{name}.paacsmesh"
        mesh.write_archive(path)
        return path.read_bytes()

    def record(self, mesh, name, box_target=.15):
        path = self.directory / f"{name}.paacsmesh"
        mesh.write_archive(path)
        self.measurements[name] = measure(read_archive(path), 1, LOWER, UPPER, .7, box_target)
        return path.read_bytes()

    def test_imported_face_and_body_tags_are_copies(self):
        self.assertEqual(self.model.face_tags, [1, 2, 3])
        self.assertEqual(self.model.body_tags, [1])
        faces = self.model.face_tags
        faces.clear()
        bodies = self.model.body_tags
        bodies.append(99)
        self.assertEqual(self.model.face_tags, [1, 2, 3])
        self.assertEqual(self.model.body_tags, [1])
        self.assertEqual(len(self.model.face_tags), self.model.info["face_count"])
        self.assertEqual(len(self.model.body_tags), self.model.info["body_count"])

    def test_controlled_routes_keep_native_checks_and_exported_tags(self):
        for variant in VARIANTS:
            with self.subTest(variant=variant):
                mesh = self.meshes[variant]
                self.assert_native_result(mesh)
                metrics = self.measurements[variant]
                self.assertEqual(metrics["tet4"], mesh.tetrahedron_count)
                self.assertEqual(sorted(metrics["tags"]), ["2:1", "2:2", "2:3", "3:1"])
                self.assertGreater(metrics["selected_face_edges_mm"]["count"], 0)

    def test_face_and_volume_controls_both_refine_their_regions(self):
        global_result, face, box, combined = (self.measurements[name] for name in VARIANTS)
        self.assertLess(face["selected_face_edges_mm"]["median"],
                        .9 * global_result["selected_face_edges_mm"]["median"])
        self.assertLess(combined["selected_face_edges_mm"]["median"],
                        .9 * box["selected_face_edges_mm"]["median"])
        self.assertGreaterEqual(combined["box_tet_median_edges_mm"]["count"], 20)
        self.assertGreaterEqual(box["box_tet_median_edges_mm"]["count"], 20)
        self.assertLess(combined["box_tet_median_edges_mm"]["median"],
                        .5 * face["box_tet_median_edges_mm"]["median"])
        self.assertLess(combined["box_tet_median_edges_mm"]["p95"],
                        .6 * face["box_tet_median_edges_mm"]["p95"])
        # Retain the tail and adjacent-size diagnostics, rather than certifying
        # a uniform physical edge size or a hard adjacent growth ratio.
        for key in ("minimum", "median", "p95", "maximum"):
            self.assertTrue(math.isfinite(combined["box_tet_median_edges_mm"][key]))
            self.assertTrue(math.isfinite(combined["near_box_adjacent_tet_size_ratio"][key]))
        self.assertIn("box_tets_above_reference_target", combined)
        self.assertIn("box_tets_above_twice_reference_target", combined)

    def test_requested_sizing_report_is_owned_and_json_safe(self):
        for variant in VARIANTS:
            with self.subTest(variant=variant):
                sizing = self.meshes[variant].report["sizing"]
                self.assertEqual(sizing["global_size"], .7)
                self.assertEqual(sizing["min_size"], .12)
                self.assertEqual(sizing["growth"], 1.55)
                self.assertIs(sizing["native_sizing"], True)
                self.assertIs(sizing["requested_native_sizing"], True)
                self.assertEqual(sizing["mode"], "spatial")
                self.assertEqual(len(sizing["face_sizes"]), int(variant in ("face", "combined")))
                self.assertEqual(len(sizing["volume_controls"]), int(variant in ("box", "combined")))
        mesh = self.meshes["combined"]
        report = mesh.report
        face = report["sizing"]["face_sizes"][0]
        box = report["sizing"]["volume_controls"][0]
        self.assertEqual(face, dict(tag=1, size=.16, effective_size=.16, growth=1.55, requested_growth=None))
        self.assertEqual(box, dict(kind="box", lower=list(LOWER), upper=list(UPPER),
                                   size=.15, effective_size=.15, growth=1.55, requested_growth=None))
        face["tag"] = 999
        box["lower"][0] = -100
        self.assertEqual(mesh.report["sizing"]["face_sizes"][0]["tag"], 1)
        self.assertEqual(mesh.report["sizing"]["volume_controls"][0]["lower"], list(LOWER))
        json.dumps(mesh.report, allow_nan=False)

    def test_minimum_clamps_requests_to_the_same_effective_field(self):
        archives = []
        for name, requested in (("clamp_low", .10), ("clamp_effective", .15)):
            mesh = self.model.mesh(.7, **self.options(
                min_size=.15, volume_controls=[pm.BoxSize(LOWER, UPPER, requested)]))
            self.assert_native_result(mesh)
            sizing = mesh.report["sizing"]
            self.assertEqual(sizing["min_size"], .15)
            self.assertEqual(sizing["volume_controls"][0]["size"], requested)
            self.assertEqual(sizing["volume_controls"][0]["effective_size"], .15)
            archives.append(self.record(mesh, name))
        self.assertEqual(archives[0], archives[1])

    def test_controls_automatically_select_the_same_native_sizing_route(self):
        mesh = self.model.mesh(.7, **self.options(native_sizing=None))
        self.assert_native_result(mesh)
        sizing = mesh.report["sizing"]
        self.assertIs(sizing["native_sizing"], True)
        self.assertIsNone(sizing["requested_native_sizing"])
        self.assertEqual(sizing["mode"], "spatial")
        self.assertEqual(self.saved(mesh, "automatic_route"), self.archives["combined"])

    def test_growth_changes_the_native_mesh_without_a_false_size_guarantee(self):
        mesh = self.model.mesh(.7, **self.options(growth=1.25))
        self.assert_native_result(mesh)
        sizing = mesh.report["sizing"]
        self.assertEqual(sizing["growth"], 1.25)
        self.assertEqual(sizing["face_sizes"][0]["growth"], 1.25)
        self.assertEqual(sizing["volume_controls"][0]["growth"], 1.25)
        self.assertNotEqual(self.record(mesh, "growth_125"), self.archives["combined"])

    def test_control_growth_overrides_and_inheritance(self):
        inherited = pm.FaceSize(1, .16)
        explicit_face = pm.FaceSize(1, .16, growth=1.55)
        explicit_box = pm.BoxSize(LOWER, UPPER, .15, growth=1.25)
        self.assertIsNone(inherited.growth)
        self.assertEqual(explicit_face.growth, 1.55)
        self.assertEqual(explicit_box.growth, 1.25)
        mesh = self.model.mesh(.7, **self.options(
            face_sizes=[explicit_face], volume_controls=[explicit_box]))
        self.assert_native_result(mesh)
        sizing = mesh.report["sizing"]
        face = sizing["face_sizes"][0]
        box = sizing["volume_controls"][0]
        self.assertEqual((face["growth"], face["requested_growth"]), (1.55, 1.55))
        self.assertEqual((box["growth"], box["requested_growth"]), (1.25, 1.25))
        self.record(mesh, "per_control_growth")

    def test_legacy_cap_sizing_toggle_keeps_the_requested_policy(self):
        mesh = self.model.mesh(.7, face_sizes=[pm.FaceSize(2, .3)],
                               min_size=.12, growth=1.55, native_sizing=False,
                               max_work=100_000_000)
        report = mesh.report
        self.assertTrue(report["audit"]["valid"])
        self.assertTrue(report["quality"]["valid"])
        self.assertEqual(report["bodies"][0]["route"], "cylinder")
        self.assertEqual(report["mesh_scope"], "independent_parts")
        self.assertFalse(report["interbody_nonoverlap_certified"])
        self.assertFalse(report["sizing"]["native_sizing"])
        self.assertIs(report["sizing"]["requested_native_sizing"], False)
        self.assertEqual(report["sizing"]["mode"], "standard")
        output = self.directory / "legacy-cap-cli.paacsmesh"
        child = subprocess.run(
            [sys.executable, "-I", "-m", "paacs_mesh", str(self.fixture), str(output),
             "--size", ".7", "--face-size", "2", ".3", "--no-native-sizing",
             "--min-size", ".12", "--growth", "1.55", "--max-work", "100000000"],
            capture_output=True, text=True, timeout=90)
        self.assertEqual(child.returncode, 0, child.stderr)
        self.assertEqual(output.read_bytes(), self.saved(mesh, "legacy-cap-python"))

    def test_typed_controls_copy_coordinates_and_are_immutable(self):
        lower, upper = list(LOWER), list(UPPER)
        face = pm.FaceSize(1, .16)
        box = pm.BoxSize(lower, upper, .15)
        lower[0], upper[2] = -99, 99
        self.assertEqual(box.lower, LOWER)
        self.assertEqual(box.upper, UPPER)
        self.assertIsNone(box.growth)
        for value, attribute, replacement in ((face, "tag", 2), (face, "size", .2),
                                               (face, "growth", 2), (box, "lower", LOWER),
                                               (box, "upper", UPPER), (box, "size", .2),
                                               (box, "growth", 2)):
            with self.subTest(attribute=attribute), self.assertRaises(AttributeError):
                setattr(value, attribute, replacement)

    def test_invalid_control_values_are_rejected(self):
        for tag in (True, 0, -1, 2 ** 31, 1.5, "1", None):
            with self.subTest(tag=tag), self.assertRaises((ValueError, TypeError, OverflowError)):
                pm.FaceSize(tag, .16)
        for size in (True, 0, -.1, math.nan, math.inf, None, ".16"):
            with self.subTest(size=size), self.assertRaises((ValueError, TypeError)):
                pm.FaceSize(1, size)
            with self.subTest(box_size=size), self.assertRaises((ValueError, TypeError)):
                pm.BoxSize(LOWER, UPPER, size)
        for growth in (True, 0, .99, math.nan, math.inf, "1.3"):
            with self.subTest(face_growth=growth), self.assertRaises((ValueError, TypeError)):
                pm.FaceSize(1, .16, growth=growth)
            with self.subTest(box_growth=growth), self.assertRaises((ValueError, TypeError)):
                pm.BoxSize(LOWER, UPPER, .15, growth=growth)
        for lower, upper in (((0, 0), (1, 1, 1)), ((0, 0, 0, 0), (1, 1, 1)),
                             ((0, 0, 0), (0, 1, 1)), ((0, 0, 0), (-1, 1, 1)),
                             ((math.nan, 0, 0), (1, 1, 1)), ((0, 0, 0), (math.inf, 1, 1)),
                             ((True, 0, 0), (1, 1, 1)), ("abc", (1, 1, 1))):
            with self.subTest(lower=lower, upper=upper), self.assertRaises((ValueError, TypeError)):
                pm.BoxSize(lower, upper, .15)

    def test_invalid_mesh_scalars_and_untyped_payloads_are_rejected(self):
        invalid = [
            {"min_size": -.1}, {"min_size": .8}, {"min_size": math.nan},
            {"min_size": math.inf}, {"min_size": True}, {"min_size": ".12"},
            {"growth": False}, {"growth": .99}, {"growth": math.nan},
            {"growth": math.inf}, {"growth": "1.55"},
            {"native_sizing": 1}, {"native_sizing": "yes"},
            {"face_sizes": [(1, .16)]}, {"face_sizes": {1: .16}},
            {"face_sizes": "wall"}, {"face_sizes": pm.FaceSize(1, .16)},
            {"volume_controls": [{"lower": LOWER, "upper": UPPER, "size": .15}]},
            {"volume_controls": [pm.FaceSize(1, .16)]},
            {"face_sizes": [pm.BoxSize(LOWER, UPPER, .15)]},
            {"native_sizing": False},  # Boxes require the native spatial route.
        ]
        for overrides in invalid:
            with self.subTest(options=overrides), self.assertRaises((ValueError, TypeError, OverflowError)):
                self.model.mesh(.7, **self.options(**overrides))

    def test_native_refusals_preserve_the_previous_mesh(self):
        failed = (
            (dict(face_sizes=[pm.FaceSize(999999, .16)]), "invalid_request"),
            (dict(face_sizes=[pm.FaceSize(1, .16), pm.FaceSize(1, .17)]), "invalid_request"),
            (dict(face_sizes=[pm.FaceSize(1, .8)]), "invalid_request"),
            (dict(volume_controls=[pm.BoxSize(LOWER, UPPER, .8)]), "invalid_request"),
            (dict(max_work=1), "budget_exceeded"),
            (dict(max_nodes=1), "budget_exceeded"),
            (dict(assembly="preserve_parts"), "unsupported_control"),
        )
        for overrides, reason in failed:
            with self.subTest(options=overrides):
                with self.assertRaises(pm.MeshingError) as caught:
                    self.model.mesh(.7, **self.options(**overrides))
                self.assertEqual(caught.exception.reason, reason)
                self.assertFalse(caught.exception.report["mesh_published"])
                self.assertEqual(caught.exception.stage, "mesh")
                json.dumps(caught.exception.report, allow_nan=False)
                if overrides.get("assembly") == "preserve_parts":
                    report = caught.exception.report
                    self.assertEqual(report["requested_assembly"], "preserve_parts")
                    self.assertEqual(report["mesh_scope"], "independent_parts")
                    self.assertFalse(report["interbody_nonoverlap_certified"])
                    self.assertIn("conforming", report.get("hint", ""))
        assembly = pm.read_step(self.fixture.with_name("assembly.step"))
        with self.assertRaises(pm.MeshingError) as caught:
            assembly.mesh(.7, **self.options())
        self.assertEqual(caught.exception.reason, "unsupported_control")
        self.assertFalse(caught.exception.report["mesh_published"])
        self.assertEqual(self.saved(self.meshes["combined"], "after_refusals"), self.archives["combined"])

    def test_cli_combined_controls_match_the_python_archive(self):
        output = self.directory / "cli.paacsmesh"
        command = [sys.executable, "-I", "-m", "paacs_mesh", str(self.fixture), str(output),
                   "--size", ".7", "--assembly", "conforming", "--native-sizing",
                   "--face-size", "1", ".16",
                   "--size-box", "-.3", "-.3", ".7", ".3", ".3", "1.3", ".15",
                   "--min-size", ".12", "--growth", "1.55", "--max-work", "100000000"]
        child = subprocess.run(command, capture_output=True, text=True, timeout=90)
        self.assertEqual(child.returncode, 0, child.stderr)
        self.assertEqual(output.read_bytes(), self.archives["combined"])
        report = json.loads(child.stdout[child.stdout.find("{\n"):])
        self.assertEqual(report["sizing"], self.meshes["combined"].report["sizing"])
        previous = output.read_bytes()
        child = subprocess.run([sys.executable, "-I", "-m", "paacs_mesh", str(self.fixture), str(output),
                                "--size", ".7", "--assembly", "conforming", "--native-sizing",
                                "--face-size", "1", ".16", "--face-size", "1", ".17"],
                               capture_output=True, text=True, timeout=90)
        self.assertNotEqual(child.returncode, 0)
        self.assertEqual(output.read_bytes(), previous)
        failure = json.loads(child.stderr)
        self.assertEqual(failure["reason"], "invalid_request")
        self.assertFalse(failure["report"]["mesh_published"])

    def test_cli_integer_overflow_is_a_clean_refusal(self):
        output = self.directory / "held-cli.paacsmesh"
        previous = self.archives["combined"]
        output.write_bytes(previous)
        for options in (("--max-work", "-1"), ("--max-work", str(2 ** 65)),
                        ("--face-size", str(2 ** 65), ".16")):
            with self.subTest(options=options):
                child = subprocess.run(
                    [sys.executable, "-I", "-m", "paacs_mesh", str(self.fixture), str(output),
                     "--size", ".7", *options], capture_output=True, text=True, timeout=90)
                self.assertEqual(child.returncode, 2, child.stderr)
                self.assertNotIn("Traceback", child.stderr)
                self.assertEqual(output.read_bytes(), previous)


if __name__ == "__main__":
    unittest.main()
