"""Native regression: surface overrides reach insertion, not only metadata."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import unittest

import paacs_mesh as pm
from test_native import fixture_path


class NativeSurfaceConstructionLimitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Intentionally strict: an older native wheel must fail these tests,
        # rather than report success after skipping the new construction API.
        cls.model = pm.read_step(fixture_path())
        cls.wall = cls.model.select_faces(surface="cylinder", radius=1, expect=1)

    def _check_surface_limit(self, points, chord_points, candidates, sites, old_limit):
        mesh = self.model.mesh(
            .5, face_sizes=[pm.FaceSize(tag, .1) for tag in self.wall.tags],
            assembly="conforming", workers=1, validation_workers=1, growth=1.3,
            construction_effort="legacy", max_surface_points=points,
            max_surface_chord_points=chord_points,
            max_interior_candidates=candidates, max_interior_sites=sites,
            max_nodes=500_000, max_cells=3_000_000, max_work=100_000_000,
            max_output_bytes=256 * 1024**2, max_working_bytes=1024**3,
            max_audit_bytes=512 * 1024**2, max_seconds=90,
        )
        report = mesh.report
        self.assertTrue(report["mesh_published"])
        for name in ("audit", "quality"):
            self.assertTrue(report[name]["complete"])
            self.assertTrue(report[name]["valid"])
        self.assertEqual(report["quality"]["inverted_cells"], 0)
        self.assertEqual(report["construction"]["effort"], "legacy")
        walls = [(body["solid_tag"], face) for body in report["bodies"]
                 for face in body["surface_spacing"] if face["face_tag"] in self.wall.tags]
        self.assertEqual(len(walls), 1)
        body_tag, wall = walls[0]
        self.assertEqual(wall["point_limit"], points)
        self.assertEqual(wall["chord_point_limit"], chord_points)
        self.assertGreater(wall["added_points"], old_limit)
        self.assertLessEqual(wall["added_points"], wall["point_limit"])
        self.assertLessEqual(wall["chord_added_points"], wall["chord_point_limit"])
        if wall["pending_entries"]:
            self.assertFalse(wall["queue_exhausted"])
            if wall["termination"] == "point_limit":
                self.assertTrue(wall["budget_limited"])
                warnings = [issue for issue in pm.diagnose(report)["issues"]
                            if issue["code"] == "construction_budget_limited"]
                self.assertTrue(any(issue["severity"] == "warning" and body_tag in issue["body_tags"]
                                    for issue in warnings))

    def test_double_surface_allowance_inserts_beyond_the_old_256_point_cap(self):
        self._check_surface_limit(512, 192, 32768, 1024, 256)

    def test_quadruple_surface_allowance_inserts_beyond_the_512_point_cap(self):
        self._check_surface_limit(1024, 384, 65536, 2048, 512)


if __name__ == "__main__":
    unittest.main()
