"""Installed, bounded body-relation queries on project-owned geometry."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import json
import unittest

import paacs_mesh as pm
from test_native import fixture_path


class BodyRelationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = pm.read_step(fixture_path().with_name("assembly.step"))

    def test_disconnected_parts_are_query_only(self):
        before = self.model.info
        report = self.model.body_relations([(2, 1)])
        self.assertEqual(report["schema"], "paacs-cad-body-relations-v1")
        self.assertTrue(report["completed"])
        self.assertTrue(report["source_model_verified_unchanged"])
        self.assertEqual(report["source_sha256"], self.model.source_sha256)
        self.assertEqual(report["model_identity"], self.model.model_identity)
        self.assertEqual(report["units"], "mm")
        self.assertEqual(len(report["pairs"]), 1)
        pair = report["pairs"][0]
        self.assertEqual([pair["first_body_tag"], pair["second_body_tag"]], [1, 2])
        self.assertEqual(pair["relation"], "disconnected")
        self.assertFalse(report["mesh_interface_certified"])
        self.assertEqual(self.model.info, before)
        json.dumps(report, allow_nan=False)
        pair["relation"] = "changed by caller"
        self.assertEqual(self.model.body_relations([(1, 2)])["pairs"][0]["relation"], "disconnected")

    def test_invalid_pairs_and_limits(self):
        bad = [([], {}), ([(1, 1)], {}), ([(1, 2), (2, 1)], {}),
               ([(True, 2)], {}), ([(1, 2, 3)], {}), ("12", {}),
               ([(1, 2)] * 65, {}), ([(1, 2)], {"max_seconds": 0}),
               ([(1, 2)], {"max_work": True}), ([(1, 2)], {"max_seconds": 1e308}),
               ([(1, 2)], {"max_working_bytes": 2 ** 31}),
               ([(1, 2)], {"distance_tolerance": float("nan")})]
        for pairs, kwargs in bad:
            with self.subTest(pairs=str(pairs)[:50], kwargs=kwargs):
                with self.assertRaises(pm.SelectionError) as caught:
                    self.model.body_relations(pairs, **kwargs)
                self.assertEqual(caught.exception.reason, "invalid_argument")
                self.assertEqual(caught.exception.stage, "body_relations")
        with self.assertRaises(pm.SelectionError) as caught:
            self.model.body_relations([(1, 999999)])
        self.assertEqual(caught.exception.reason, "unknown_body")

    def test_resource_failure_keeps_receipt_and_no_partial_result(self):
        with self.assertRaises(pm.SelectionError) as caught:
            self.model.body_relations([(1, 2)], max_working_bytes=1)
        self.assertEqual(caught.exception.reason, "budget_exceeded")
        self.assertFalse(caught.exception.report["completed"])
        self.assertEqual(caught.exception.report["pairs"], [])
        self.assertEqual(caught.exception.report["source_sha256"], self.model.source_sha256)
        self.assertFalse(caught.exception.report["mesh_interface_certified"])

    def test_direct_native_caps_are_checked_before_query(self):
        native = self.model._native
        for options in ({"max_work": 10_000_001},
                        {"max_working_bytes": 1024 ** 3 + 1},
                        {"max_work": True}, {"max_seconds": 0}, {"max_seconds": 1e308},
                        {"distance_tolerance": float("nan")}):
            with self.subTest(options=options):
                with self.assertRaises((TypeError, ValueError)):
                    native._body_relations([(1, 2)], **options)
        with self.assertRaises(ValueError):
            native._body_relations([(1, 2)] * 65)

    def test_direct_native_captures_selected_sequence_length_once(self):
        class ChangingLength:
            def __init__(self):
                self.length_calls = 0
                self.indices = []

            def __len__(self):
                self.length_calls += 1
                # A repeated length check or iterator end must not enlarge the
                # native conversion after its original input admission.
                return 1 if self.length_calls <= 3 else 1_000_000

            def __getitem__(self, index):
                self.indices.append(index)
                if index != 0:
                    raise AssertionError("native conversion exceeded the admitted pair count")
                return (1, 2)

        selected = ChangingLength()
        report = self.model._native._body_relations(selected)
        self.assertTrue(report["completed"])
        self.assertEqual(report["requested_pairs"], 1)
        self.assertEqual(selected.length_calls, 1)
        self.assertEqual(selected.indices, [0])


if __name__ == "__main__":
    unittest.main()
