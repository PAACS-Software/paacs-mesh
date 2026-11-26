"""Source-bound comparisons use independent cell populations and one spatial grid."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import json
import math
from pathlib import Path
import struct
import tempfile
import unittest

import paacs_mesh as pm
from paacs_mesh.mesh_archive import InspectionError
from paacs_mesh.mesh_comparison import compare_meshes
from paacs_mesh.mesh_inspection import METRICS, MeshInspection
from test_native import fixture_path


RIGHT = [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.)]
REGULAR = [(0., 0., 0.), (1., 0., 0.), (.5, math.sqrt(3) / 2, 0.),
           (.5, math.sqrt(3) / 6, math.sqrt(2 / 3))]


def write_bundle(directory, *, points=None, cells=None, source="a" * 64, units="mm", changes=None):
    """Synthetic report pairing exercises interpretation, not a native recertification."""
    directory.mkdir()
    points = RIGHT if points is None else points
    cells = [(2, 2, 11, (0, 2, 1)), (5, 3, 7, (0, 1, 2, 3))] if cells is None else cells
    payload = bytearray(struct.pack("<8sIIQQQ", b"PAACSM8\0", 1, 0, len(points), len(cells),
                                    sum(len(row[3]) for row in cells)))
    for point in points:
        payload.extend(struct.pack("<3dii", *point, 3, 0))
    for kind, dimension, tag, nodes in cells:
        payload.extend(struct.pack("<B3sii", kind, b"\0\0\0", dimension, tag))
        payload.extend(struct.pack("<" + "q" * len(nodes), *nodes))
    report = {
        "status": "succeeded", "reason": "none", "mesh_published": True,
        "source_sha256": source, "kernel_version": "synthetic-test-kernel-A", "units": units, "nodes": len(points), "cells": len(cells),
        "tetrahedra": sum(row[0] == 5 for row in cells), "mesh_scope": "conforming_domain",
        "interbody_nonoverlap_certified": True,
        "audit": {"complete": True, "valid": True, "missing_boundary": 0,
                  "unexpected_boundary": 0, "nonmanifold_faces": 0, "inconsistent_orientation": 0,
                  "invalid_connectivity": 0, "nonfinite_nodes": 0},
        "quality": {"complete": True, "valid": True, "inverted_cells": 0,
                    "minimum": .5, "warning_threshold": 0., "minimum_required": 0.},
        "bodies": [{"solid_tag": 7, "reason": "none", "cad_volume": 1., "mesh_volume": .99,
                    "relative_volume_error": .01, "maximum_boundary_node_error": .002,
                    "maximum_boundary_sample_error": .003, "surface_chord_error_limit": .1,
                    "centroid_outside_distance_limit": .05}], "interferences": [],
        "sizing": {"global_size": 1., "requested_min_size": 0., "min_size": .1, "growth": 1.4,
                   "requested_native_sizing": True,
                   "face_sizes": [{"tag": 11, "size": .2, "effective_size": .2,
                                   "requested_growth": 1.2, "growth": 1.2}],
                   "volume_controls": [{"kind": "box", "lower": [0., 0., 0.], "upper": [.5, .5, .5],
                                        "size": .15, "effective_size": .15,
                                        "requested_growth": 1.3, "growth": 1.3}]},
        "effective_policy": {"wrapper_native_sizing": True, "body_routes": []},
        "seconds": {"import": .01, "meshing_total": .2, "surface": .04, "volume": .12},
    }
    report.update(changes or {})
    report_bytes = json.dumps(report, allow_nan=False, separators=(",", ":")).encode()
    (directory / "mesh.paacsmesh").write_bytes(payload)
    (directory / "report.json").write_bytes(report_bytes)
    (directory / "inspection.json").write_text(json.dumps({
        "schema": "paacs-mesh-inspection-bundle-v1", "mesh_sha256": hashlib.sha256(payload).hexdigest(),
        "report_sha256": hashlib.sha256(report_bytes).hexdigest(), "source_sha256": source, "units": units,
    }), encoding="utf-8")
    return directory


class MeshComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.sequence = 0

    def inspected(self, *, max_bytes=256 * 1024 * 1024, **arguments):
        self.sequence += 1
        path = write_bundle(self.root / str(self.sequence), **arguments)
        snapshot = MeshInspection(path, max_bytes=max_bytes)
        self.addCleanup(snapshot.close)
        return snapshot

    def test_quality_distributions_counts_sizing_and_reported_fidelity(self):
        first = self.inspected()
        second = self.inspected(points=REGULAR, changes={
            "seconds": {"import": .01, "meshing_total": .1, "surface": .02, "volume": .06},
            "bodies": [{"solid_tag": 7, "reason": "none", "cad_volume": 1., "mesh_volume": .995,
                        "relative_volume_error": .005, "maximum_boundary_node_error": .001}],
        })
        comparison = compare_meshes(first, second)
        self.assertEqual(comparison["schema"], "paacs-mesh-comparison-v1")
        self.assertEqual(comparison["counts_delta"], {"node_count": 0, "cell_count": 0, "tetrahedron_count": 0})
        self.assertEqual(comparison["first"]["cell_counts"]["tri3"], 1)
        self.assertEqual(comparison["second"]["tetrahedron_count"], 1)
        for name in METRICS:
            self.assertEqual(comparison["first"]["quality"][name]["count"], 1)
            expected = second._metric(1, name) - first._metric(1, name)
            self.assertAlmostEqual(comparison["quality_delta"][name]["median"], expected)
            for key in ("minimum", "p01", "p05", "median", "p95", "maximum", "mean"):
                self.assertIsNotNone(comparison["second"]["quality"][name][key])
        sizing = comparison["first"]["sizing"]
        self.assertEqual(sizing["requested"]["face_sizes"][0]["tag"], 11)
        self.assertEqual(sizing["effective"]["resolved_routing"], first.native_report["effective_policy"])
        self.assertEqual(sizing["effective"]["sizing_completion_certificate"], "unavailable")
        self.assertFalse(sizing["achieved_growth"]["available"])
        self.assertEqual(sizing["achieved"]["boundary_faces"]["11"]["count"], 3)
        self.assertEqual(sizing["achieved"]["tet4"]["edge_max"]["count"], 1)
        body = comparison["cad_error"]["bodies"][0]
        self.assertEqual(body["body_tag"], 7)
        self.assertAlmostEqual(body["delta"]["maximum_boundary_node_error"], -.001)
        self.assertAlmostEqual(body["delta"]["relative_volume_error"], -.005)
        self.assertIsNone(body["delta"]["maximum_boundary_sample_error"])
        self.assertEqual(comparison["first"]["validation"]["basis"], "hash_paired_saved_native_report_not_recertified")
        self.assertFalse(comparison["second"]["validation"]["independent_native_reaudit_performed"])
        json.dumps(comparison, allow_nan=False)

    def test_shared_grid_prevents_false_changes_from_different_unused_extents(self):
        points = RIGHT + [(x + 4., y, z) for x, y, z in RIGHT]
        cells = [(5, 3, 7, (0, 1, 2, 3)), (5, 3, 7, (4, 5, 6, 7))]
        first = self.inspected(points=points, cells=cells)
        second = self.inspected(points=points + [(100., 0., 0.)], cells=cells)
        comparison = compare_meshes(first, second)
        self.assertNotEqual(first.mesh_sha256, second.mesh_sha256)
        grid = comparison["spatial_comparison"]
        self.assertEqual(grid["bounds"], [[0., 0., 0.], [100., 1., 1.]])
        self.assertEqual(grid["affected_region_indices"], [])
        self.assertEqual(len(grid["regions"]), 1)
        self.assertEqual(grid["regions"][0]["first"]["count"], 2)
        self.assertEqual(grid["regions"][0]["delta"]["count"], 0)

    def test_changed_regions_do_not_pair_identical_cell_ids_across_hashes(self):
        first = self.inspected(cells=[(5, 3, 7, (0, 1, 2, 3))])
        second = self.inspected(points=[(x + 4., y, z) for x, y, z in RIGHT],
                                cells=[(5, 3, 7, (0, 1, 2, 3))])
        comparison = compare_meshes(first, second)
        regions = comparison["spatial_comparison"]["regions"]
        self.assertEqual(len(regions), 2)
        self.assertEqual(sorted(row["delta"]["count"] for row in regions), [-1, 1])
        self.assertTrue(comparison["cell_correspondence"].startswith("not_assumed"))
        self.assertNotIn("paired_cells", comparison)
        for row in regions:
            self.assertEqual(len(row["index"]), 3)
            self.assertEqual(len(row["bounds"]), 2)

    def test_low_quality_counts_and_unknown_body_identity_are_measured(self):
        first = self.inspected(cells=[(5, 3, 7, (0, 1, 2, 3))])
        thin = [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (.3, .3, .00001)]
        second = self.inspected(points=thin, cells=[(5, 3, 0, (0, 1, 2, 3))])
        comparison = compare_meshes(first, second, spatial_bins=1)
        row = comparison["spatial_comparison"]["regions"][0]
        self.assertEqual(row["first"]["body_tags"], [7])
        self.assertEqual(row["second"]["body_tags"], [])
        self.assertEqual(row["second"]["unclassified_body_count"], 1)
        self.assertEqual(row["delta"]["below_0_2"], 1)
        self.assertLess(row["delta"]["minimum_mean_ratio"], 0)
        self.assertEqual(comparison["second"]["body_tet4_counts"], {})

    def test_native_timing_never_guesses_workflow_totals_or_retries(self):
        first = self.inspected(changes={"attempts": [{"seconds": 99.}], "elapsed_seconds": 123.})
        second = self.inspected()
        comparison = compare_meshes(first, second)
        timing = comparison["timings"]
        self.assertFalse(timing["end_to_end_runtime_comparison_available"])
        for profile in (timing["first"], timing["second"]):
            self.assertEqual(profile["meshing_scope"], "final_native_attempt_only")
            self.assertFalse(profile["retry_timing_available"])
            self.assertIsNone(profile["workflow_total_seconds"])
            self.assertIsNone(profile["failed_retry_seconds"])
            self.assertIsNone(profile["retry_count"])
        self.assertEqual(timing["reported_seconds_delta"]["meshing_total"], 0.)

    def test_source_units_validity_and_input_types_are_required(self):
        good = self.inspected()
        for other, reason in ((self.inspected(source="b" * 64), "comparison_source"),
                              (self.inspected(source=None), "comparison_source"),
                              (self.inspected(units="m"), "comparison_units"),
                              (self.inspected(changes={"mesh_published": False, "status": "failed"}),
                               "comparison_validation")):
            with self.subTest(reason=reason), self.assertRaises(InspectionError) as caught:
                compare_meshes(good, other)
            self.assertEqual(caught.exception.reason, reason)
        bare_path = good._archive.original_path
        with MeshInspection(bare_path) as bare, self.assertRaises(InspectionError) as caught:
            compare_meshes(good, bare)
        self.assertEqual(caught.exception.reason, "comparison_source")
        with self.assertRaises(TypeError):
            compare_meshes(good, {})

    def test_warning_status_is_accepted_and_duplicate_body_claims_refused(self):
        good = self.inspected()
        warning = self.inspected(changes={"quality_warning": True})
        self.assertEqual(compare_meshes(good, warning)["second"]["validation"]["reported_native_status"],
                         "valid_with_warnings")
        duplicate = self.inspected(changes={"bodies": [{"solid_tag": 7, "reason": "none"},
                                                     {"solid_tag": 7, "reason": "none"}]})
        with self.assertRaises(InspectionError) as caught:
            compare_meshes(good, duplicate)
        self.assertEqual(caught.exception.reason, "comparison_report")

    def test_preserved_assembly_interference_remains_explicit(self):
        separate = self.inspected(changes={"mesh_scope": "separate_bodies",
            "interbody_nonoverlap_certified": False,
            "interferences": [{"first_solid_tag": 7, "second_solid_tag": 9,
                               "state": "material_overlap", "overlap_volume": .01}]})
        comparison = compare_meshes(separate, separate)
        for name in ("first", "second"):
            self.assertEqual(comparison[name]["validation"]["reported_native_status"], "valid_with_warnings")
            self.assertFalse(comparison[name]["assembly"]["interbody_nonoverlap_certified"])
            self.assertEqual(comparison[name]["assembly"]["interferences"][0]["state"], "material_overlap")

    def test_work_time_memory_and_argument_budgets_are_enforced(self):
        first, second = self.inspected(), self.inspected()
        with self.assertRaises(InspectionError) as caught:
            compare_meshes(first, second, max_work=1)
        self.assertEqual(caught.exception.reason, "work_budget")
        with self.assertRaises(InspectionError) as caught:
            compare_meshes(first, second, max_seconds=1e-12)
        self.assertEqual(caught.exception.reason, "time_budget")
        smaller = self.inspected(max_bytes=4 * 1024 * 1024)
        with self.assertRaises(InspectionError) as caught:
            compare_meshes(first, smaller)
        self.assertEqual(caught.exception.reason, "memory_budget")
        for keyword, value in (("spatial_bins", 0), ("spatial_bins", 17), ("spatial_bins", True),
                               ("max_work", 0), ("max_seconds", 0), ("max_seconds", math.inf)):
            with self.subTest(keyword=keyword, value=value), self.assertRaises(ValueError):
                compare_meshes(first, second, **{keyword: value})

    def test_missing_or_changed_kernel_does_not_pair_CAD_body_tags(self):
        first = self.inspected()
        for kernel in (None, "synthetic-test-kernel-B"):
            with self.subTest(kernel=kernel):
                second = self.inspected(changes={"kernel_version": kernel})
                comparison = compare_meshes(first, second)
                self.assertEqual(comparison["counts_delta"]["tetrahedron_count"], 0)
                self.assertFalse(comparison["cad_error"]["body_comparison_available"])
                self.assertEqual(comparison["cad_error"]["bodies"], [])
                self.assertEqual(comparison["cad_error"]["first_body_measurements"][0]["body_tag"], 7)
                self.assertEqual(comparison["cad_error"]["second_body_measurements"][0]["relative_volume_error"], .01)

    def test_output_budget_rejects_large_valid_input_metadata(self):
        # Each source report is within 8 MiB. Their combined preserved policy
        # exceeds the output limit; explicitly admit their conservative scratch.
        huge = self.inspected(max_bytes=2 * 1024 ** 3,
                              changes={"effective_policy": {"actual_native_note": "x" * (4 * 1024 ** 2)}})
        with self.assertRaises(InspectionError) as caught:
            compare_meshes(huge, huge)
        self.assertEqual(caught.exception.reason, "output_budget")

    def test_result_is_owned_and_input_reports_are_unchanged(self):
        first, second = self.inspected(), self.inspected()
        original = first.native_report
        comparison = compare_meshes(first, second)
        comparison["first"]["sizing"]["effective"]["resolved_routing"]["wrapper_native_sizing"] = False
        comparison["first"]["identity"]["source_sha256"] = "changed"
        self.assertEqual(first.native_report, original)
        self.assertEqual(first.source_sha256, "a" * 64)

    def test_extreme_finite_bounds_produce_finite_json(self):
        scale = 1e300
        first = self.inspected(points=[(x * scale - scale / 2, y * scale, z * scale) for x, y, z in RIGHT],
                               cells=[(5, 3, 7, (0, 1, 2, 3))])
        second = self.inspected(points=[(x * scale + scale / 2, y * scale, z * scale) for x, y, z in RIGHT],
                                cells=[(5, 3, 7, (0, 1, 2, 3))])
        comparison = compare_meshes(first, second)
        self.assertTrue(math.isfinite(comparison["first"]["sizing"]["achieved"]["tet4"]
                                      ["volume_equivalent_edge"]["minimum"]))
        json.dumps(comparison, allow_nan=False)


class NativeMeshComparisonTests(unittest.TestCase):
    def test_owned_cylinder_native_meshes_and_saved_bundle(self):
        model = pm.read_step(fixture_path())
        coarse = model.mesh(.7, workers=1, validation_workers=1)
        finer = model.mesh(.5, workers=1, validation_workers=1)
        with MeshInspection(coarse) as first, MeshInspection(finer) as second:
            comparison = compare_meshes(first, second)
            self.assertEqual(comparison["source_sha256"], model.source_sha256)
            self.assertEqual(comparison["units"], "mm")
            self.assertEqual(comparison["first"]["node_count"], coarse.node_count)
            self.assertEqual(comparison["second"]["node_count"], finer.node_count)
            self.assertEqual(comparison["first"]["sizing"]["requested"]["global_size"], .7)
            self.assertEqual(comparison["second"]["sizing"]["requested"]["global_size"], .5)
            for name in METRICS:
                self.assertGreater(comparison["second"]["quality"][name]["count"], 0)
                self.assertAlmostEqual(comparison["quality_delta"][name]["median"],
                                       second.summary()["quality"][name]["median"] -
                                       first.summary()["quality"][name]["median"])
            self.assertFalse(comparison["timings"]["end_to_end_runtime_comparison_available"])
            self.assertTrue(all(row["body_tag"] > 0 for row in comparison["cad_error"]["bodies"]))
            with tempfile.TemporaryDirectory() as directory:
                saved = Path(directory) / "saved"
                second.save(saved)
                with MeshInspection(saved) as reopened:
                    same = compare_meshes(second, reopened)
                    self.assertEqual(same["spatial_comparison"]["affected_region_indices"], [])
                    self.assertEqual(same["quality_delta"]["mean_ratio"]["minimum"], 0.)


if __name__ == "__main__":
    unittest.main()
