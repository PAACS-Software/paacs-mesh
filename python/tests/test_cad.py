"""Installed CAD inspection, identity and guarded selection checks."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest

import paacs_mesh as pm
from test_native import fixture_path


class CadInspectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = fixture_path()
        cls.model = pm.read_step(cls.fixture)

    def test_source_identity_and_copied_info(self):
        expected = hashlib.sha256(self.fixture.read_bytes()).hexdigest()
        self.assertEqual(self.model.source_sha256, expected)
        self.assertEqual(self.model.info["source_sha256"], expected)
        self.assertEqual(self.model.kernel_version, "8.0.1")
        self.assertNotEqual(self.model.model_identity, pm.read_step(self.fixture).model_identity)
        altered = self.model.info
        altered["source_sha256"] = "0" * 64
        self.assertEqual(self.model.info["source_sha256"], expected)
        mesh = self.model.mesh(0.7)
        self.assertEqual(mesh.report["source_sha256"], expected)

    def test_cylinder_metrics_and_orientation(self):
        page = self.model.inspect(topology=True)
        self.assertEqual(page["total"], 3)
        self.assertIsNone(page["next_offset"])
        faces = page["entities"]
        self.assertEqual(sorted(face["surface"] for face in faces), ["cylinder", "plane", "plane"])
        side = next(face for face in faces if face["surface"] == "cylinder")
        self.assertAlmostEqual(side["area"], 4 * math.pi, places=7)
        self.assertAlmostEqual(side["geometry"]["radius"], 1, places=12)
        self.assertIsNone(side["normal"])
        self.assertTrue(side["edge_tags"])
        for face in faces:
            self.assertEqual(face["owners"][0]["body_tag"], 1)
            self.assertGreater(face["tolerance"], 0)
            self.assertEqual(len(face["bounds"]["lower"]), 3)
        top = self.model.select_faces(normal=[0, 0, 1])
        bottom = self.model.select_faces(normal=[0, 0, -1])
        self.assertNotEqual(top.tags, bottom.tags)
        self.assertAlmostEqual(top.guards[0]["area"], math.pi, places=7)
        self.assertTrue(top.validate())
        body = self.model.inspect("bodies")["entities"][0]
        self.assertAlmostEqual(body["volume"], 2 * math.pi, places=9)
        self.assertEqual(len(body["face_uses"]), 3)
        self.assertEqual(len(self.model.inspect("edges")["entities"]), len(self.model.edge_tags))
        self.assertEqual(len(self.model.inspect("vertices")["entities"]), len(self.model.vertex_tags))

    def test_selection_predicates_cardinality_and_ownership(self):
        side = self.model.select_faces(radius=1)
        self.assertEqual(len(side), 1)
        caps = self.model.select_faces(surface="plane", area=math.pi, expect=2)
        self.assertEqual(len(caps), 2)
        self.assertEqual(self.model.select_faces(tags=list(caps)).tags, caps.tags)
        with self.assertRaises(pm.SelectionError) as caught:
            self.model.select_faces(surface="plane")
        self.assertEqual(caught.exception.reason, "cardinality")
        with self.assertRaises(pm.SelectionError):
            self.model.select_faces(tags=[99999])
        with self.assertRaises(pm.SelectionError) as caught:
            pm.read_step(self.fixture).define_group("foreign", caps)
        self.assertEqual(caught.exception.reason, "foreign_selection")
        changed = caps.query
        changed["surface"] = "sphere"
        self.assertEqual(caps.query["surface"], "plane")
        neighbours = self.model.select_faces(adjacent_to=list(side), expect=2)
        self.assertEqual(neighbours.tags, caps.tags)

    def test_inspection_pagination_and_budgets(self):
        page = self.model.inspect(limit=1)
        self.assertEqual(len(page["entities"]), 1)
        self.assertEqual(page["next_offset"], 1)
        self.assertEqual(self.model.inspect(offset=99)["entities"], [])
        for kwargs in ({"limit": 1001}, {"offset": -1}, {"limit": True}, {"kind": "anything"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(pm.SelectionError):
                self.model.inspect(**kwargs)
        with self.assertRaises(pm.SelectionError) as caught:
            self.model.select_faces(surface="plane", expect=2, max_faces=1)
        self.assertEqual(caught.exception.reason, "scan_budget")
        json.dumps(self.model.inspect(), allow_nan=False)

    def test_embedded_representation_names_and_disabled_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "named.step"
            text = self.fixture.read_text(encoding="ascii")
            text = text.replace("#17 = ADVANCED_FACE('',", "#17 = ADVANCED_FACE('bore wall',")
            text = text.replace("#105 = ADVANCED_FACE('',", "#105 = ADVANCED_FACE('cap',")
            text = text.replace("#109 = ADVANCED_FACE('',", "#109 = ADVANCED_FACE('cap',")
            text = text.replace("MANIFOLD_SOLID_BREP('',", "MANIFOLD_SOLID_BREP('part A',")
            source.write_text(text, encoding="ascii")
            model = pm.read_step(source)
            self.assertTrue(model.info["representation_names"]["complete"])
            self.assertEqual(model.select_faces(name="bore wall").tags, model.select_faces(surface="cylinder").tags)
            self.assertEqual(len(model.select_faces(name="cap", expect=2)), 2)
            with self.assertRaises(pm.SelectionError) as caught:
                model.select_faces(name="cap")
            self.assertEqual(caught.exception.reason, "cardinality")
            disabled = pm.read_step(source, names=False)
            with self.assertRaises(pm.SelectionError) as caught:
                disabled.select_faces(name="bore wall")
            self.assertEqual(caught.exception.reason, "names_disabled")
            a, b = Path(directory) / "a.paacsmesh", Path(directory) / "b.paacsmesh"
            model.mesh(0.7).write(a)
            disabled.mesh(0.7).write(b)
            self.assertEqual(a.read_bytes(), b.read_bytes())

    def test_intent_round_trip_and_atomic_refusal(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "cylinder.paacs.json"
            model = pm.read_step(self.fixture)
            caps = model.select_faces(surface="plane", expect=2)
            group = model.define_group("caps", caps, group_id="end-caps-v1")
            self.assertEqual(group["id"], "end-caps-v1")
            model.save_intent(manifest, revision="test-1")
            reopened = pm.read_step(self.fixture, intent=manifest)
            self.assertEqual(reopened.select_faces(group="caps").tags, caps.tags)
            self.assertEqual(reopened.groups["caps"]["id"], "end-caps-v1")
            broken = json.loads(manifest.read_text(encoding="utf-8"))
            broken["groups"][0]["guards"][0]["area"] += 1
            manifest.write_text(json.dumps(broken), encoding="utf-8")
            before = reopened.groups
            with self.assertRaises(pm.IntentError):
                reopened.load_intent(manifest)
            self.assertEqual(reopened.groups, before)
            manifest.write_text('{"schema":"x","schema":"y"}', encoding="utf-8")
            with self.assertRaises(pm.IntentError):
                reopened.load_intent(manifest)

    def test_unicode_and_quoted_representation_name(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "unicode.step"
            value = self.fixture.read_text(encoding="ascii").replace(
                "#17 = ADVANCED_FACE('',", r"#17 = ADVANCED_FACE('Paul''s \X2\00E9\X0\ wall',")
            source.write_text(value, encoding="ascii")
            model = pm.read_step(source)
            named = model.select_faces(name="Paul's é wall")
            self.assertEqual(named.tags, model.select_faces(surface="cylinder").tags)

    def test_source_changed_and_retained_snapshot(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "copy.step"
            payload = self.fixture.read_bytes()
            path.write_bytes(payload)
            retained = pm.read_step(path, retain_source=True)
            ordinary = pm.read_step(path)
            expected = retained.source_sha256
            path.write_bytes(payload + b"\n")
            self.assertEqual(ordinary.source_sha256, expected)
            self.assertEqual(retained._source_snapshot(), payload)
            with self.assertRaises(pm.IntentError) as caught:
                ordinary._source_snapshot()
            self.assertEqual(caught.exception.reason, "source_changed")
            path.unlink()
            self.assertEqual(retained._source_snapshot(), payload)
            self.assertEqual(retained.inspect()["total"], 3)

    def test_explicit_semantic_rebind_preserves_group_id(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, changed = root / "cylinder.paacs.json", root / "edited.step"
            model = pm.read_step(self.fixture)
            model.define_group("side", model.select_faces(surface="cylinder"), group_id="side-1")
            model.save_intent(manifest)
            changed.write_bytes(self.fixture.read_bytes() + b"\n")
            edited = pm.read_step(changed)
            with self.assertRaises(pm.IntentError):
                edited.load_intent(manifest)
            edited.load_intent(manifest, rebind=True)
            self.assertEqual(edited.groups["side"]["id"], "side-1")
            self.assertEqual(edited.select_faces(group="side").tags, model.select_faces(group="side").tags)


if __name__ == "__main__":
    unittest.main()
