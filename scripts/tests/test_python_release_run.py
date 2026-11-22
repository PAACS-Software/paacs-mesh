"""Cross-run publication admits only immutable, fully qualified release artifacts."""
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from copy import deepcopy
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import check_python_release_run as release


class QualifiedRunTests(unittest.TestCase):
    def setUp(self):
        self.sha, self.workflow_sha = "a" * 40, "b" * 40
        self.current, self.previous = "101", "202"
        self.tag = "v1.2.0-rc.2"
        self.base = {"repository": {"full_name": release.REPOSITORY, "id": 7},
            "head_repository": {"full_name": release.REPOSITORY, "id": 7},
            "event": "workflow_dispatch", "head_branch": "main", "path": release.WORKFLOW,
            "workflow_id": 17, "head_sha": self.workflow_sha, "status": "in_progress"}
        self.old = dict(deepcopy(self.base), head_sha=self.sha, status="completed", conclusion="failure")
        self.jobs = [{"name": name, "status": "completed", "conclusion": "success"}
                     for name in sorted(release.QUALIFY | {"prepare"})]
        self.jobs += [{"name": "assemble", "status": "completed", "conclusion": "failure"}]
        self.jobs += [{"name": name, "status": "completed", "conclusion": "skipped"}
                      for name in ("source-release", "publish", "verify-index")]
        self.artifacts = [{"name": name, "id": index + 1, "expired": False,
            "expires_at": "2099-01-01T00:00:00Z", "digest": "sha256:" + "c" * 64,
            "workflow_run": {"id": int(self.previous), "head_sha": self.sha,
                             "head_branch": "main", "repository_id": 7, "head_repository_id": 7}}
            for index, name in enumerate(release.BUILDERS + release.CONSUMERS)]
        self.calls = []
        self.addCleanup(patch.stopall)
        patch.object(release, "release_source", return_value=(self.sha, "1.2.0-rc.2")).start()
        patch.object(release, "published", return_value=False).start()
        patch.object(release, "github", side_effect=self.api).start()

    def api(self, path):
        self.calls.append(path)
        if path == "actions/runs/" + self.current:
            return deepcopy(self.base)
        if path == "actions/runs/" + self.previous:
            return deepcopy(self.old)
        if path == "releases/tags/" + self.tag:
            return {"draft": False, "tag_name": self.tag}
        if path.startswith("actions/runs/" + self.previous + "/jobs?"):
            return {"total_count": len(self.jobs), "jobs": deepcopy(self.jobs)}
        if path.startswith("actions/runs/" + self.previous + "/artifacts?"):
            return {"total_count": len(self.artifacts), "artifacts": deepcopy(self.artifacts)}
        self.fail("Unexpected API request: " + path)

    def prepare(self, resume=True):
        return release.prepare(Path("."), self.tag, self.workflow_sha, self.current,
                               self.previous if resume else "")

    def test_normal_publication_qualifies_the_tagged_source(self):
        result = self.prepare(False)
        self.assertEqual(result["source-sha"], self.sha)
        self.assertEqual(result["artifact-run-id"], self.current)
        self.assertEqual(result["reuse-qualified"], "false")
        self.assertEqual(result["builder-artifact-ids"], "")
        self.assertFalse(any(self.previous in path for path in self.calls))

    def test_resume_binds_original_source_and_all_six_immutable_artifact_ids(self):
        result = self.prepare()
        self.assertEqual(result, {"source-sha": self.sha, "artifact-run-id": self.previous,
            "reuse-qualified": "true", "builder-artifact-ids": "1,2,3",
            "consumer-artifact-ids": "4,5,6"})
        self.assertNotEqual(result["source-sha"], self.workflow_sha)

    def test_run_identifiers_are_bounded_and_cannot_select_the_current_run(self):
        for value in ("0", "-1", "1e3", "123/456", "1\n2", "1" * 21, self.current):
            with self.subTest(value=value), self.assertRaises(RuntimeError):
                release.prepare(Path("."), self.tag, self.workflow_sha, self.current, value)

    def test_source_run_must_be_the_exact_trusted_main_release_workflow(self):
        original = deepcopy(self.old)
        for key, value in (("head_sha", "d" * 40), ("event", "pull_request"),
                           ("head_branch", "feature"), ("path", ".github/workflows/other.yml"),
                           ("workflow_id", 99), ("status", "in_progress"),
                           ("repository", {"full_name": "other/project"}),
                           ("head_repository", {"full_name": "fork/paacs-mesh"})):
            self.old = dict(deepcopy(original), **{key: value})
            with self.subTest(key=key), self.assertRaises(RuntimeError):
                self.prepare()
        self.old = original

    def test_missing_failed_or_duplicated_qualification_job_refuses_resume(self):
        original = deepcopy(self.jobs)
        mutations = [original[:-1], original + [deepcopy(original[0])]]
        for conclusion in ("failure", "skipped", "cancelled"):
            changed = deepcopy(original)
            changed[0]["conclusion"] = conclusion
            mutations.append(changed)
        for jobs in mutations:
            self.jobs = jobs
            with self.assertRaises(RuntimeError):
                self.prepare()
        self.jobs = original

    def test_a_previous_upload_attempt_is_never_resumed(self):
        original = deepcopy(self.jobs)
        for name in ("publish", "source-release"):
            for conclusion in ("success", "failure"):
                self.jobs = deepcopy(original)
                next(j for j in self.jobs if j["name"] == name)["conclusion"] = conclusion
                with self.subTest(name=name, conclusion=conclusion), self.assertRaises(RuntimeError):
                    self.prepare()
        self.jobs = original

    def test_missing_or_duplicate_consumer_artifacts_refuse_resume(self):
        original = deepcopy(self.artifacts)
        for rows in (original[:-1], original + [deepcopy(original[-1])]):
            self.artifacts = rows
            with self.assertRaises(RuntimeError):
                self.prepare()
        self.artifacts = original

    def test_expired_or_unidentified_artifact_refuses_resume(self):
        original = deepcopy(self.artifacts)
        for change in ({"expired": True}, {"expires_at": "2020-01-01T00:00:00Z"},
                       {"digest": None}, {"digest": "sha256:wrong"}, {"id": 0}):
            self.artifacts = deepcopy(original)
            self.artifacts[0].update(change)
            with self.subTest(change=change), self.assertRaises((RuntimeError, TypeError)):
                self.prepare()
        self.artifacts = original

    def test_artifact_source_run_commit_and_repository_must_match(self):
        original = deepcopy(self.artifacts)
        for change in ({"id": 999}, {"head_sha": "d" * 40}, {"head_branch": "feature"},
                       {"repository_id": 99}, {"head_repository_id": 99}):
            self.artifacts = deepcopy(original)
            self.artifacts[0]["workflow_run"].update(change)
            with self.subTest(change=change), self.assertRaises(RuntimeError):
                self.prepare()
        self.artifacts = original

    def test_current_workflow_and_existing_pypi_version_are_checked(self):
        self.base["head_sha"] = "d" * 40
        with self.assertRaises(RuntimeError):
            self.prepare()
        self.base["head_sha"] = self.workflow_sha
        release.published.return_value = True
        with self.assertRaises(RuntimeError):
            self.prepare()

    def test_metadata_pagination_is_bounded(self):
        with patch.object(release, "github", return_value={"total_count": 501, "jobs": [{}] * 100}):
            with self.assertRaises(RuntimeError):
                release.pages("actions/runs/202/jobs?filter=latest", "jobs")


class TaggedSourceTests(unittest.TestCase):
    def test_tagged_source_is_resolved_separately_from_clean_workflow_revision(self):
        source, workflow = "a" * 40, "b" * 40
        with patch.object(release, "source_revision", return_value={"commit": workflow, "clean": True}), \
                patch.object(release.subprocess, "check_output", side_effect=["tag", source, "1.2.0-rc.2\n"]), \
                patch.object(release.subprocess, "run") as ancestor:
            self.assertEqual(release.release_source(Path("."), "v1.2.0-rc.2", workflow),
                             (source, "1.2.0-rc.2"))
            self.assertEqual(ancestor.call_args.args[0][-2:], [source, workflow])

    def test_dirty_workflow_lightweight_tag_mismatched_version_and_nonancestor_refuse(self):
        source, workflow = "a" * 40, "b" * 40
        with patch.object(release, "source_revision", return_value={"commit": workflow, "clean": False}):
            with self.assertRaises(RuntimeError):
                release.release_source(Path("."), "v1.2.0-rc.2", workflow)
        for responses in (["commit"], ["tag", source, "1.2.0-rc.1"]):
            with patch.object(release, "source_revision", return_value={"commit": workflow, "clean": True}), \
                    patch.object(release.subprocess, "check_output", side_effect=responses):
                with self.assertRaises(RuntimeError):
                    release.release_source(Path("."), "v1.2.0-rc.2", workflow)
        with patch.object(release, "source_revision", return_value={"commit": workflow, "clean": True}), \
                patch.object(release.subprocess, "check_output", side_effect=["tag", source, "1.2.0-rc.2"]), \
                patch.object(release.subprocess, "run", side_effect=subprocess.CalledProcessError(1, "git")):
            with self.assertRaises(subprocess.CalledProcessError):
                release.release_source(Path("."), "v1.2.0-rc.2", workflow)


if __name__ == "__main__":
    unittest.main()
