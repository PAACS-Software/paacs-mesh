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


if __name__ == "__main__":
    unittest.main()
