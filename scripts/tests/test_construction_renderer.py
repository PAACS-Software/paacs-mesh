"""Metadata/file guards for retained-run figures; no CAD or native meshing."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
BENCHMARKS = ROOT / "benchmarks"
with patch.object(sys, "path", [str(BENCHMARKS), *sys.path]):
    SPEC = importlib.util.spec_from_file_location("construction_renderer_test_module", BENCHMARKS / "render_construction_budgets.py")
    RENDERER = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(RENDERER)


class ReachedRendering(RuntimeError):
    pass


class FakeInspection:
    """Isolate pre-render guards; this fixture makes no native validity claim."""
    node_count = 4
    validation = {"basis": RENDERER.SAVED_BASIS, "native_reaudit_performed": False}

    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass

    def summary(self, **_):
        metric = {"minimum": .5, "p05": .5, "median": .5}
        return {"cell_counts": {"tet4": 1}, "bounds": [[0., 0., 0.], [1., 1., 1.]],
                "quality": {name: dict(metric) for name in ("mean_ratio", "minimum_dihedral_deg", "scaled_jacobian")}}

    def query(self, **_):
        return {"total": 0}

    def visualize(self, *_args, **_options):
        raise ReachedRendering("Explicit comparison guards passed")


class ConstructionRendererGuardTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.cases = []
        for index in range(2):
            case_dir = self.directory / str(index)
            case_dir.mkdir()
            payload = f"unparsed archive payload {index}".encode()
            (case_dir / "mesh.paacsmesh").write_bytes(payload)
            case = {"schema": "paacs-construction-budget-case-v1", "status": "succeeded", "native_status": "succeeded",
                    "retry_count": 0, "workers": 1, "validation_workers": 1, "source_sha256": "a" * 64,
                    "wall_size_mm": .1, "global_size_mm": .5, "native_seconds": 1., "nodes": 4, "tetrahedra": 1,
                    "native_extension_sha256": "b" * 64, "mesh_sha256": hashlib.sha256(payload).hexdigest(),
                    "global_max_work": 100000000, "global_max_working_bytes": 1024**3,
                    "options": {"construction_effort": "legacy"}, "profile": "legacy"}
            report = {"mesh_published": True, "reason": "none", "audit": {"complete": True, "valid": True},
                      "quality": {"complete": True, "valid": True}, "source_sha256": "a" * 64,
                      "units": "mm", "kernel_version": "8.0.1", "requested_size": .5,
                      "sizing": {"face_sizes": [{"tag": 1, "size": .1}]},
                      "seconds": {"meshing_total": 1.}, "construction": {"effort": "legacy"}, "nodes": 4, "cells": 1}
            (case_dir / "case.json").write_text(json.dumps(case), encoding="utf-8")
            (case_dir / "native-report.json").write_text(json.dumps(report), encoding="utf-8")
            self.cases.append(case_dir)
        self.args = SimpleNamespace(case=[["A", str(self.cases[0])], ["B", str(self.cases[1])]],
            output=self.directory / "new-output", development_native=None,
            max_input_bytes=8 * RENDERER.MIB, max_output_bytes=8 * RENDERER.MIB,
            max_bytes=64 * RENDERER.MIB, max_cells=1000, max_nodes=1000, max_work=100000,
            max_seconds=10, total_seconds=30, allow_mixed_native=False, allow_mixed_budgets=False,
            view_bounds=None, cut_axis="y", cut_value=0., highlight_below=.1, azimuth=35., elevation=25.,
            max_triangles=1000, max_highlights=1000, columns=2, png="none")
        package = ModuleType("paacs_mesh")
        package.__path__ = []
        package.inspect_mesh = lambda *_args, **_kwargs: FakeInspection()
        workflow = ModuleType("paacs_mesh.workflow")
        workflow._rename_new_directory = lambda *_: self.fail("No test should publish an output directory")
        self.modules = {"paacs_mesh": package, "paacs_mesh.workflow": workflow}

    def change_case(self, index, **fields):
        path = self.cases[index] / "case.json"
        case = json.loads(path.read_text())
        path.write_text(json.dumps(dict(case, **fields)), encoding="utf-8")

    def render(self):
        with patch.dict(sys.modules, self.modules):
            return RENDERER.render(self.args)

    def change_report(self, index, **fields):
        path = self.cases[index] / "native-report.json"
        report = json.loads(path.read_text())
        path.write_text(json.dumps(dict(report, **fields)), encoding="utf-8")

    def add_resources(self, index, policy="normal", work=None, memory=512 * RENDERER.MIB):
        self.change_case(index, options={"construction_effort": "legacy", "resource_policy": policy},
                         global_max_work=work, global_max_working_bytes=None)
        self.change_report(index, resources={"policy": policy, "requested_max_work": work,
            "effective_max_work": work if work is not None else (2**62 - 1 if policy == "normal" else 50_000_000),
            "logical_work_cap_disabled": policy == "normal" and work is None,
            "max_working_bytes": memory})

    def pair_report(self, index):
        digest = hashlib.sha256((self.cases[index] / "native-report.json").read_bytes()).hexdigest()
        self.change_case(index, native_report_file="native-report.json", native_report_sha256=digest)

    def caption(self, index):
        case = json.loads((self.cases[index] / "case.json").read_text())
        report = json.loads((self.cases[index] / "native-report.json").read_text())
        stats = {"tet4": 1, "nodes": 4, "native_seconds": 1., "below_threshold": 0,
                 "threshold": .1, "native_sha256": "b" * 64,
                 "effective_resources": RENDERER._resource_receipt(case, report)}
        stats.update({name: {"minimum": .5, "p05": .5, "median": .5}
                      for name in ("mean_ratio", "dihedral", "jacobian")})
        return "\n".join(RENDERER._caption(stats))

    def test_changed_archive_hash_is_rejected_without_publication(self):
        before = (self.cases[0] / "mesh.paacsmesh").read_bytes()
        self.change_case(0, mesh_sha256="0" * 64)
        with self.assertRaisesRegex(ValueError, "Archive differs"):
            self.render()
        self.assertFalse(self.args.output.exists())
        self.assertEqual((self.cases[0] / "mesh.paacsmesh").read_bytes(), before)

    def test_mixed_global_budgets_require_explicit_flag(self):
        self.change_case(1, global_max_work=50000000, global_max_working_bytes=512 * RENDERER.MIB)
        with self.assertRaisesRegex(ValueError, "allow-mixed-budgets"):
            self.render()
        self.assertFalse(self.args.output.exists())
        self.args.allow_mixed_budgets = True
        with self.assertRaises(ReachedRendering):
            self.render()
        self.assertFalse(self.args.output.exists())

    def test_mixed_recorded_native_binaries_require_explicit_flag(self):
        self.change_case(1, native_extension_sha256="c" * 64)
        with self.assertRaisesRegex(ValueError, "allow-mixed-native"):
            self.render()
        self.assertFalse(self.args.output.exists())
        self.args.allow_mixed_native = True
        with self.assertRaises(ReachedRendering):
            self.render()


    def test_normal_policy_is_verified_in_native_resources(self):
        for index in range(2):
            self.add_resources(index)
        with self.assertRaises(ReachedRendering):
            self.render()
        caption = self.caption(0)
        self.assertIn("no practical work ceiling (normal)", caption)
        self.assertIn("512 MiB modeled", caption)
        self.assertNotIn("not recorded", caption)
        self.assertFalse(self.args.output.exists())

    def test_requested_native_resource_policy_mismatch_is_rejected(self):
        self.add_resources(0)
        report = json.loads((self.cases[0] / "native-report.json").read_text())
        report["resources"]["policy"] = "bounded"
        self.change_report(0, resources=report["resources"])
        with self.assertRaisesRegex(ValueError, "differ on resource policy"):
            self.render()
        self.assertFalse(self.args.output.exists())

    def test_requested_construction_cap_mismatch_is_rejected(self):
        self.change_case(0, options={"construction_effort": "legacy", "max_surface_points": 32})
        self.change_report(0, construction={"effort": "legacy", "max_surface_points": 16})
        with self.assertRaisesRegex(ValueError, "construction settings"):
            self.render()
        self.assertFalse(self.args.output.exists())

    def test_requested_global_work_mismatch_is_rejected(self):
        self.add_resources(0, work=100_000_000)
        self.change_case(0, global_max_work=50_000_000)
        with self.assertRaisesRegex(ValueError, "requested global work"):
            self.render()
        self.assertFalse(self.args.output.exists())

    def test_same_omitted_limits_with_different_effective_policies_are_mixed(self):
        self.add_resources(0, "normal")
        self.add_resources(1, "bounded")
        with self.assertRaisesRegex(ValueError, "allow-mixed-budgets"):
            self.render()
        self.args.allow_mixed_budgets = True
        with self.assertRaises(ReachedRendering):
            self.render()
        self.assertFalse(self.args.output.exists())

    def test_same_omitted_limits_with_different_effective_memory_are_mixed(self):
        self.add_resources(0, memory=512 * RENDERER.MIB)
        self.add_resources(1, memory=1024 * RENDERER.MIB)
        with self.assertRaisesRegex(ValueError, "allow-mixed-budgets"):
            self.render()
        self.assertFalse(self.args.output.exists())

    def test_legacy_receipt_does_not_invent_effective_limits(self):
        caption = self.caption(0)
        self.assertIn("Global work effective unavailable (requested 1e+08)", caption)
        self.assertIn("working allowance effective unavailable (requested 1024 MiB)", caption)
        self.assertNotIn("no practical work ceiling", caption)
        with self.assertRaises(ReachedRendering):
            self.render()
        self.assertFalse(self.args.output.exists())

    def test_global_only_runs_accept_no_wall_request_or_face_controls(self):
        for index in range(2):
            self.change_case(index, wall_size_mm=None)
            self.change_report(index, sizing={"face_sizes": [], "volume_controls": []})
        with self.assertRaises(ReachedRendering):
            self.render()
        self.assertFalse(self.args.output.exists())

    def test_global_only_runs_refuse_hidden_face_controls(self):
        self.change_case(0, wall_size_mm=None)
        with self.assertRaisesRegex(ValueError, "global-only.*face-sizing"):
            self.render()
        self.assertFalse(self.args.output.exists())

    def test_global_only_size_must_remain_positive_and_match_native_request(self):
        case = {"wall_size_mm": None, "global_size_mm": .5}
        report = {"requested_size": .5, "sizing": {"face_sizes": []}}
        self.assertEqual(RENDERER._sizing_targets(case, report), (None, .5))
        for value in (0, -1, None, True):
            with self.subTest(value=value), self.assertRaises(ValueError):
                RENDERER._sizing_targets(dict(case, global_size_mm=value), report)
        with self.assertRaisesRegex(ValueError, "agree on the positive global size"):
            RENDERER._sizing_targets(case, dict(report, requested_size=.4))
        self.assertFalse(self.args.output.exists())

    def test_matching_native_report_file_hash_passes_without_native_calls(self):
        for index in range(2):
            self.pair_report(index)
        with self.assertRaises(ReachedRendering):
            self.render()
        self.assertFalse(self.args.output.exists())

    def test_substituted_native_report_is_rejected_before_publication(self):
        self.pair_report(0)
        archive = (self.cases[0] / "mesh.paacsmesh").read_bytes()
        self.change_report(0, nodes=400)
        with self.assertRaisesRegex(ValueError, "Native report differs.*file SHA-256"):
            self.render()
        self.assertFalse(self.args.output.exists())
        self.assertEqual((self.cases[0] / "mesh.paacsmesh").read_bytes(), archive)

    def test_native_report_pairing_requires_expected_filename(self):
        self.pair_report(0)
        self.change_case(0, native_report_file="other-report.json")
        with self.assertRaisesRegex(ValueError, "filename must be native-report.json"):
            self.render()
        self.assertFalse(self.args.output.exists())

    def test_available_native_node_and_tet_counts_must_match_archive(self):
        for fields in ({"nodes": 5}, {"tetrahedra": 2}, {"tetrahedra": True}):
            with self.subTest(fields=fields):
                self.change_report(0, nodes=4, tetrahedra=1)
                self.change_report(0, **fields)
                with self.assertRaisesRegex(ValueError, "native report counts differ"):
                    self.render()
                self.assertFalse(self.args.output.exists())

    def test_effective_work_must_match_declared_global_policy(self):
        for policy, requested, changed in (("normal", 100_000_000, 50_000_000),
                ("normal", 100_000_000, 100_000_001),
                ("bounded", None, 100_000_000), ("normal", None, 50_000_000)):
            with self.subTest(policy=policy, requested=requested, effective=changed):
                self.add_resources(0, policy, work=requested)
                case = json.loads((self.cases[0] / "case.json").read_text())
                report = json.loads((self.cases[0] / "native-report.json").read_text())
                report["resources"]["effective_max_work"] = changed
                with self.assertRaisesRegex(ValueError, "effective work|bounded default work|inconsistent practical ceiling"):
                    RENDERER._resource_receipt(case, report)


if __name__ == "__main__":
    unittest.main()
