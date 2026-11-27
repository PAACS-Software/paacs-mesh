"""Pure artifact/export tests; fake handles do not certify Native geometry."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import json
import os
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch

from paacs_mesh.pipeline_export import GeometryArtifact


def archive(kind="surface", *, invalid_index=False):
    nodes = [(0., 0., 0., 0, 1), (1., 0., 0., 0, 2), (0., 1., 0., 2, 3)]
    ids = (0, 1, 2) if kind == "surface" else (0, 1)
    if invalid_index:
        ids = (*ids[:-1], 99)
    header = struct.pack("<8sIIQQQ", b"PAACSM8\0", 1, 0, len(nodes), 1, len(ids))
    data = b"".join(struct.pack("<3dii", *node) for node in nodes)
    cell = struct.pack("<B3sii", 2 if kind == "surface" else 0, b"\0\0\0", 2 if kind == "surface" else 1, 3)
    return header + data + cell + struct.pack("<" + "Q" * len(ids), *ids)


class Handle:
    def __init__(self, kind="surface", payload=None):
        self.payload = archive(kind) if payload is None else payload
        self.report = {"prepared": True, "nodes": 3, "triangles": 1} if kind == "surface" else {
            "complete": True, "nodes": 3, "segments": 1, "unused_nodes": 1,
            "parameter_order": "increasing native source parameter", "edges": [
                {"edge_tag": 3, "natural_node_order": [0, 1], "natural_line_order": [0],
                 "natural_parameters": [0., 1.], "body_tags": [1, 2], "face_tags": [3, 4],
                 "uses": [{"face_tag": 3, "sense": -1, "seam": True}]}]}
        self.writes = []
        self.arrays = {"points": object(), "cells": object()}

    def archive(self, maximum):
        if len(self.payload) > maximum:
            raise ValueError("native archive byte allowance")
        return self.payload

    def write(self, path, maximum):
        data = self.payload if Path(path).suffix == ".paacsmesh" else b"# vtk DataFile Version 3.0\nfixture\nASCII\n"
        if len(data) > maximum:
            raise ValueError("native output byte allowance")
        self.writes.append((Path(path), maximum))
        Path(path).write_bytes(data)

    def to_numpy(self, *, max_bytes):
        self.array_limit = max_bytes
        return dict(self.arrays)


def artifact(kind="surface", *, payload=None, protected=()):
    native = Handle(kind, payload)
    metadata = {"kind": kind, "units": "mm", "source_sha256": "a" * 64,
                "boundary_sha256": hashlib.sha256(archive()).hexdigest(),
                "cad_topology": {"faces": [{"tag": 3, "owners": [{"body_tag": 1}, {"body_tag": 2}]}]}}
    return GeometryArtifact(native, metadata, protected_paths=protected), native, metadata


class PipelineExportTests(unittest.TestCase):
    def test_identity_metadata_and_native_report_are_copied(self):
        view, native, metadata = artifact("cad_edges")
        metadata["units"] = "metres"
        native.report["edges"].clear()
        copy = view.report
        copy["edges"].clear()
        self.assertEqual(view.units, "mm")
        self.assertEqual(view.source_sha256, "a" * 64)
        self.assertEqual(view.report["edges"][0]["body_tags"], [1, 2])
        self.assertTrue(view.report["edges"][0]["uses"][0]["seam"])

    def test_archives_keep_unused_boundary_nodes_and_distinguish_kind(self):
        for kind in ("surface", "cad_edges"):
            view, native, _ = artifact(kind)
            self.assertEqual(view.archive(), native.payload)
            self.assertEqual(struct.unpack_from("<Q", view.archive(), 16)[0], 3)
        view, _, _ = artifact("cad_edges", payload=archive())
        with self.assertRaisesRegex(ValueError, "header/type/count|unexpected cells"):
            view.archive()

    def test_connectivity_and_node_finiteness_refuse(self):
        for payload in (archive(invalid_index=True),
                        archive()[:40] + struct.pack("<d", float("nan")) + archive()[48:]):
            view, _, _ = artifact(payload=payload)
            with self.assertRaises(ValueError):
                view.archive()

    def test_boundary_identity_is_verified_on_surface_export(self):
        view, native, _ = artifact()
        changed = bytearray(native.payload)
        struct.pack_into("<d", changed, 40, 0.1)
        native.payload = bytes(changed)
        with self.assertRaisesRegex(ValueError, "frozen boundary"):
            view.archive()

    def test_invalid_metadata_and_incomplete_native_stage_refuse(self):
        _, native, metadata = artifact()
        for key, value in (("kind", "tet4"), ("units", "m"), ("source_sha256", "A" * 64)):
            with self.subTest(key=key), self.assertRaises(ValueError):
                GeometryArtifact(native, {**metadata, key: value})
        native.report["prepared"] = False
        with self.assertRaisesRegex(ValueError, "completed native"):
            GeometryArtifact(native, metadata)

    def test_array_extraction_delegates_with_original_ids_and_copied_metadata(self):
        view, native, _ = artifact("cad_edges")
        data = view.to_numpy(max_bytes=256)
        self.assertIs(data["points"], native.arrays["points"])
        self.assertEqual(native.array_limit, 256)
        self.assertIs(data["volume_qualified"], False)
        data["artifact_metadata"]["units"] = "bad"
        self.assertEqual(view.units, "mm")

    def test_single_file_writes_delegate_without_unannounced_sidecars(self):
        view, native, _ = artifact()
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "surface.paacsmesh"
            self.assertEqual(view.write(target), target.absolute())
            self.assertEqual(target.read_bytes(), native.payload)
            self.assertEqual([item.name for item in Path(directory).iterdir()], [target.name])
            with self.assertRaisesRegex(ValueError, "support"):
                view.write(Path(directory) / "surface.msh")

    def test_source_and_hardlink_are_protected(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.vtk"
            source.write_bytes(b"original source")
            view, native, _ = artifact(protected=(source,))
            with self.assertRaises(ValueError):
                view.write(source)
            self.assertEqual(source.read_bytes(), b"original source")
            link = Path(directory) / "link.vtk"
            try:
                os.link(source, link)
            except OSError:
                return
            with self.assertRaises(ValueError):
                view.write(link)
            self.assertEqual(native.writes, [])

    def test_bundle_hashes_metadata_incidence_and_no_volume_claim(self):
        view, _, _ = artifact("cad_edges")
        with tempfile.TemporaryDirectory() as directory:
            saved = view.save(Path(directory) / "edges")
            manifest = json.loads((saved / "manifest.json").read_text())
            self.assertFalse(manifest["volume_qualified"])
            self.assertEqual(set(manifest["files"]), {"mesh.paacsmesh", "mesh.vtk", "metadata.json", "report.json"})
            for name, record in manifest["files"].items():
                data = (saved / name).read_bytes()
                self.assertEqual(record, {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
            metadata = json.loads((saved / "metadata.json").read_text())
            self.assertEqual(metadata["node_id_space"], "original_boundary_zero_based")
            self.assertEqual(metadata["cad_topology"]["faces"][0]["owners"], [{"body_tag": 1}, {"body_tag": 2}])
            report = json.loads((saved / "report.json").read_text())
            self.assertEqual(report["edges"][0]["natural_parameters"], [0., 1.])
            self.assertEqual(report["unused_nodes"], 1)
            self.assertFalse((saved / "source.step").exists())

    def test_bundle_cannot_replace_existing_directory(self):
        view, native, _ = artifact()
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                view.save(directory)
            self.assertEqual(native.writes, [])

    def test_refused_aggregate_budget_publishes_nothing(self):
        view, _, _ = artifact()
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "output"
            with self.assertRaises(ValueError):
                view.save(target, max_output_bytes=200)
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_refused_native_writer_removes_staging(self):
        view, native, _ = artifact()
        with tempfile.TemporaryDirectory() as directory, patch.object(native, "write", side_effect=ValueError("refused")):
            with self.assertRaisesRegex(ValueError, "refused"):
                view.save(Path(directory) / "output")
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_work_and_clock_refuse_before_writing(self):
        view, native, _ = artifact()
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "work"):
                view.write(Path(directory) / "out.vtk", max_work=1)
            with patch("paacs_mesh.pipeline_export.time.perf_counter", side_effect=[0., 2.]):
                with self.assertRaisesRegex(ValueError, "time"):
                    view.write(Path(directory) / "out.vtk", max_seconds=1)
            self.assertEqual(native.writes, [])

    def test_invalid_limits_refuse_without_native_write(self):
        view, native, _ = artifact()
        for options in ({"max_bytes": True}, {"max_work": 0}, {"max_seconds": 0}, {"max_seconds": float("inf")}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                view.archive(**options)
        self.assertEqual(native.writes, [])


if __name__ == "__main__":
    unittest.main()
