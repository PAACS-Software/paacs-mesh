"""Body-relation receipt semantics without native geometry work."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import unittest

from paacs_mesh._cad_common import SelectionError
from paacs_mesh.body_relations import inspect_body_relations


class _Native:
    def __init__(self, result):
        self.result = result
        self.calls = []

    def _body_relations(self, pairs, **options):
        self.calls.append((pairs, options))
        return self.result


class _Model:
    source_sha256 = "a" * 64
    model_identity = 73
    body_tags = [1, 2, 3]

    def __init__(self, result):
        self._native = _Native(result)

    def _assert_identity(self):
        pass


def _receipt(*, completed=True):
    return {
        "schema": "paacs-cad-body-relations-v1",
        "source_sha256": _Model.source_sha256,
        "model_identity": _Model.model_identity,
        "units": "mm",
        "completed": completed,
        "reason": "none" if completed else "budget_exceeded",
        "diagnostic": "body_relations_measured" if completed else "body_relations_work_budget",
        "source_audit_available": True,
        "source_model_verified_unchanged": True,
        "mesh_interface_certified": False,
        "requested_pairs": 1,
        "pairs": [{"first_body_tag": 1, "second_body_tag": 2,
                   "relation": "unresolved", "evidence": "native_boundary_common",
                   "issue": "near_but_unproved", "minimum_distance": .001}],
    }


class BodyRelationReceiptTests(unittest.TestCase):
    def test_completed_unresolved_pair_is_a_measured_query(self):
        native = _receipt()
        model = _Model(native)
        result = inspect_body_relations(model, [(2, 1)])
        self.assertTrue(result["completed"])
        self.assertEqual(result["pairs"][0]["relation"], "unresolved")
        self.assertFalse(result["mesh_interface_certified"])
        self.assertEqual(model._native.calls[0][0], [[1, 2]])
        result["pairs"][0]["relation"] = "caller changed its copy"
        self.assertEqual(native["pairs"][0]["relation"], "unresolved")

    def test_failed_batch_preserves_partial_measurements_as_diagnostics(self):
        native = _receipt(completed=False)
        native["pairs"][0].update(relation="face_contact", issue="none", native_contact_area=1.)
        model = _Model(native)
        with self.assertRaises(SelectionError) as caught:
            inspect_body_relations(model, [(1, 2)])
        error = caught.exception
        self.assertEqual((error.stage, error.reason), ("body_relations", "budget_exceeded"))
        self.assertFalse(error.report["completed"])
        self.assertEqual(error.report["pairs"][0]["native_contact_area"], 1.)
        self.assertFalse(error.report["mesh_interface_certified"])
        error.report["pairs"][0]["native_contact_area"] = 999.
        self.assertEqual(native["pairs"][0]["native_contact_area"], 1.)

    def test_mismatched_source_or_model_is_never_returned(self):
        for field, value in [("source_sha256", "b" * 64), ("model_identity", 74)]:
            with self.subTest(field=field):
                native = _receipt()
                native[field] = value
                with self.assertRaises(SelectionError) as caught:
                    inspect_body_relations(_Model(native), [(1, 2)])
                self.assertEqual((caught.exception.stage, caught.exception.reason),
                                 ("body_relations", "model_changed"))

    def test_invalid_selection_never_enters_native(self):
        model = _Model(_receipt())
        with self.assertRaises(SelectionError) as caught:
            inspect_body_relations(model, [(1, 99)])
        self.assertEqual(caught.exception.reason, "unknown_body")
        self.assertEqual(model._native.calls, [])

    def test_identity_failure_is_scoped_to_body_relations(self):
        model = _Model(_receipt())

        def changed():
            raise SelectionError("model_changed", "changed before query", {"identity": 74})

        model._assert_identity = changed
        with self.assertRaises(SelectionError) as caught:
            inspect_body_relations(model, [(1, 2)])
        self.assertEqual((caught.exception.stage, caught.exception.reason),
                         ("body_relations", "model_changed"))
        self.assertEqual(caught.exception.report, {"identity": 74})
        self.assertEqual(model._native.calls, [])

    def test_scalar_bounds_refuse_before_native(self):
        bad = [
            {"max_work": 0}, {"max_work": 10_000_001}, {"max_work": True},
            {"max_working_bytes": 0}, {"max_working_bytes": 1024 ** 3 + 1},
            {"max_seconds": 0}, {"max_seconds": float("inf")},
            {"distance_tolerance": -1}, {"distance_tolerance": True},
        ]
        for options in bad:
            with self.subTest(options=options):
                model = _Model(_receipt())
                with self.assertRaises(SelectionError) as caught:
                    inspect_body_relations(model, [(1, 2)], **options)
                self.assertEqual((caught.exception.stage, caught.exception.reason),
                                 ("body_relations", "invalid_argument"))
                self.assertEqual(model._native.calls, [])

    def test_explicit_maxima_and_zero_separation_threshold_are_forwarded(self):
        model = _Model(_receipt())
        inspect_body_relations(model, [(1, 2)], distance_tolerance=0,
                               max_work=10_000_000, max_working_bytes=1024 ** 3,
                               max_seconds=.01)
        pairs, options = model._native.calls[0]
        self.assertEqual(pairs, [[1, 2]])
        self.assertEqual(options, {"distance_tolerance": 0., "max_work": 10_000_000,
                                  "max_working_bytes": 1024 ** 3, "max_seconds": .01})


if __name__ == "__main__":
    unittest.main()
