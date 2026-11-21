"""Regression checks for companion-group body scope and atomic ownership guards."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import json
from pathlib import Path
import tempfile
import unittest

import paacs_mesh as pm
from test_native import fixture_path


class CadGroupScopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = fixture_path().with_name("assembly.step")

    def test_saved_group_retains_its_body_scope_for_further_selection(self):
        with tempfile.TemporaryDirectory(prefix="paacs-group-scope-") as directory:
            model = pm.read_step(self.fixture)
            self.assertEqual(model.body_tags, [1, 2])
            selected = model.select_faces(body=1, normal=[0, 0, 1])
            self.assertEqual(selected.body, 1)
            model.define_group("body one top", selected, "body-one-top")
            manifest = Path(directory) / "assembly.intent.json"
            model.save_intent(manifest)
            reopened = pm.read_step(self.fixture, intent=manifest)
            restored = reopened.select_faces(group="body one top")
            self.assertEqual(restored.tags, selected.tags)
            self.assertEqual(restored.body, 1)
            self.assertTrue(restored.validate())
            filtered = reopened.select_faces(group="body one top", normal=[0, 0, 1])
            self.assertEqual(filtered.tags, selected.tags)
            self.assertEqual(filtered.body, 1)
            with self.assertRaises(pm.SelectionError):
                reopened.select_faces(group="body one top", body=2)

    def test_existing_but_wrong_owner_in_second_group_refuses_atomically(self):
        with tempfile.TemporaryDirectory(prefix="paacs-group-owner-") as directory:
            root = Path(directory)
            model = pm.read_step(self.fixture)
            top_one = model.select_faces(body=1, normal=[0, 0, 1])
            model.define_group("a valid", top_one, "valid-id")
            model.define_group("z invalid", top_one, "invalid-id")
            manifest = root / "assembly.intent.json"
            model.save_intent(manifest, revision="incoming")
            payload = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertEqual(payload["groups"][1]["name"], "z invalid")
            broken = payload["groups"][1]
            self.assertEqual(broken["guards"][0]["owners"], [1])
            # Body 2 is valid in the model, but it does not own this face.
            broken["body"] = broken["query"]["body"] = 2
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            target = pm.read_step(self.fixture)
            target.define_group("existing", target.select_faces(body=2, normal=[0, 0, 1]), "existing-id")
            target.save_intent(root / "previous.json", revision="stable")
            before = target.groups
            with self.assertRaises(pm.IntentError) as caught:
                target.load_intent(manifest)
            self.assertEqual(caught.exception.stage, "intent")
            self.assertEqual(target.groups, before)
            receipt = target.save_intent(root / "after-refusal.json")
            self.assertEqual(receipt["revision"], "stable")
            self.assertEqual(target.select_faces(group="existing").body, 2)


if __name__ == "__main__":
    unittest.main()
