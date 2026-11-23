#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Resolve release source and admit only fully qualified, unpublished artifacts."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

from packaging.version import Version

from build_python_wheel import source_revision
from python_compatibility import PLATFORM_TAGS, consumers, platforms, targets


REPOSITORY = "PAACS-Software/paacs-mesh"
WORKFLOW = ".github/workflows/publish-python.yml"
BUILDERS = tuple(row["artifact"] for row in targets())
CONSUMERS = tuple(row["consumer_artifact"] for row in consumers())
SDK_ARTIFACTS = {row["sdk_artifact"] for row in platforms()}
QUALIFY = {"qualify / " + name for name in (
    ["release-assembly-checks", "compatibility-matrix"]
    + ["OCCT SDK (" + platform + ")" for platform in PLATFORM_TAGS]
    + [row["job_name"] for row in targets()]
    + [row["job_name"] for row in consumers()])}


def run_id(value: str) -> str:
    if not re.fullmatch(r"[1-9][0-9]{0,19}", value):
        raise RuntimeError("Expected a bounded positive GitHub run ID")
    return value


def github(path: str) -> dict:
    token = os.environ.get("GH_TOKEN")
    if not token:
        raise RuntimeError("GitHub Actions read access is required")
    request = Request("https://api.github.com/repos/" + REPOSITORY + "/" + path,
        headers={"Authorization": "Bearer " + token, "Accept": "application/vnd.github+json",
                 "X-GitHub-Api-Version": "2022-11-28"})
    try:
        with urlopen(request, timeout=30) as response:
            data = response.read(5_000_001)
        if len(data) > 5_000_000:
            raise RuntimeError("GitHub release response exceeds its budget")
        return json.loads(data)
    except (HTTPError, URLError) as error:
        raise RuntimeError("GitHub release metadata could not be verified") from error


def pages(path: str, key: str) -> list[dict]:
    records = []
    separator = "&" if "?" in path else "?"
    for page in range(1, 6):
        result = github(f"{path}{separator}per_page=100&page={page}")
        rows = result[key]
        records.extend(rows)
        if len(records) == result["total_count"]:
            return records
        if len(rows) != 100:
            break
    raise RuntimeError("GitHub release metadata is incomplete or exceeds its page budget")


def published(version: str) -> bool:
    try:
        with urlopen("https://pypi.org/pypi/paacs-mesh/" + quote(version, safe="") + "/json",
                     timeout=15):
            return True
    except HTTPError as error:
        if error.code == 404:
            return False
        raise RuntimeError("PyPI publication status could not be verified") from error
    except URLError as error:
        raise RuntimeError("PyPI publication status could not be verified") from error


def release_source(source: Path, tag: str, workflow_sha: str) -> tuple[str, str]:
    if (not re.fullmatch(r"[0-9a-f]{40}", workflow_sha)
            or not re.fullmatch(r"v[A-Za-z0-9.!+_-]{1,80}", tag)):
        raise RuntimeError("Invalid workflow revision or release tag")
    revision = source_revision(source)
    if revision.get("commit") != workflow_sha or revision.get("clean") is not True:
        raise RuntimeError("The release workflow requires its exact clean source revision")
    def git(*args):
        return subprocess.check_output(["git", "-C", str(source), *args], text=True).strip()
    reference = "refs/tags/" + tag
    if git("cat-file", "-t", reference) != "tag":
        raise RuntimeError("The release requires an annotated immutable tag")
    commit = git("rev-parse", "--verify", reference + "^{commit}")
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise RuntimeError("Release tag did not resolve to a commit")
    version = git("show", commit + ":VERSION.txt")
    if tag != "v" + version:
        raise RuntimeError("Release tag and tagged VERSION.txt differ")
    Version(version)
    subprocess.run(["git", "-C", str(source), "merge-base", "--is-ancestor", commit, workflow_sha],
                   check=True, capture_output=True)
    return commit, version


def trusted_run(run: dict, sha: str, *, completed: bool = False):
    if (run.get("repository", {}).get("full_name") != REPOSITORY
            or run.get("head_repository", {}).get("full_name") != REPOSITORY
            or run.get("event") != "workflow_dispatch" or run.get("head_branch") != "main"
            or run.get("head_sha") != sha or run.get("path") != WORKFLOW
            or type(run.get("workflow_id")) is not int or run["workflow_id"] <= 0
            or (completed and run.get("status") != "completed")):
        raise RuntimeError("Artifacts require a trusted main release dispatch at the tagged source")


def qualification(identifier: str, sha: str, workflow_id: int) -> tuple[str, str]:
    previous = github("actions/runs/" + identifier)
    trusted_run(previous, sha, completed=True)
    if previous["workflow_id"] != workflow_id or previous.get("conclusion") != "failure":
        raise RuntimeError("Resume requires the same publication workflow stopped before upload")
    jobs = pages("actions/runs/" + identifier + "/jobs?filter=latest", "jobs")
    expected = QUALIFY | {"prepare", "assemble", "source-release", "publish", "verify-index"}
    if len(jobs) != len(expected) or {j["name"] for j in jobs} != expected:
        raise RuntimeError("The source run lacks the complete publication qualification job set")
    for job in jobs:
        wanted = ("success" if job["name"] in QUALIFY | {"prepare"} else
                  "failure" if job["name"] == "assemble" else "skipped")
        if job.get("status") != "completed" or job.get("conclusion") != wanted:
            raise RuntimeError("Resume requires all qualification passes and no attempted publication")
    artifacts = pages("actions/runs/" + identifier + "/artifacts", "artifacts")
    required = set(BUILDERS + CONSUMERS)
    selected = {a["name"]: a for a in artifacts if a["name"] in required}
    if (set(selected) != required or sum(a["name"] in required for a in artifacts) != len(required)
            or len({a["name"] for a in artifacts}) != len(artifacts)
            or len({a.get("id") for a in artifacts}) != len(artifacts)
            or any(a["name"] not in required | SDK_ARTIFACTS for a in artifacts)):
        raise RuntimeError("The source run lacks the complete distinct builder/consumer artifact matrix")
    now = datetime.now(timezone.utc)
    for artifact in selected.values():
        run = artifact.get("workflow_run", {})
        if (type(artifact.get("id")) is not int or artifact["id"] <= 0
                or artifact.get("expired") is not False
                or datetime.fromisoformat(artifact["expires_at"].replace("Z", "+00:00")) <= now
                or not isinstance(artifact.get("digest"), str)
                or not re.fullmatch(r"sha256:[0-9a-f]{64}", artifact["digest"])
                or str(run.get("id")) != identifier or run.get("head_sha") != sha
                or run.get("head_branch") != "main"
                or type(run.get("repository_id")) is not int or run["repository_id"] <= 0
                or run["repository_id"] != previous["repository"].get("id")
                or run.get("repository_id") != run.get("head_repository_id")):
            raise RuntimeError("A qualified artifact is expired, mutable or from another source")
    return tuple(",".join(str(selected[name]["id"]) for name in names)
                 for names in (BUILDERS, CONSUMERS))


def prepare(source: Path, tag: str, workflow_sha: str, current_id: str,
            qualified_id: str = "") -> dict[str, str]:
    current_id = run_id(current_id)
    sha, version = release_source(source, tag, workflow_sha)
    current = github("actions/runs/" + current_id)
    trusted_run(current, workflow_sha)
    release = github("releases/tags/" + quote(tag, safe=""))
    if release.get("draft") is not False or release.get("tag_name") != tag:
        raise RuntimeError("The tagged GitHub release must already be published")
    if published(str(Version(version))):
        raise RuntimeError("This version is already published to PyPI; do not repeat upload")
    builders = consumers = ""
    if qualified_id:
        qualified_id = run_id(qualified_id)
        if qualified_id == current_id:
            raise RuntimeError("A release cannot resume its own unfinished qualification")
        builders, consumers = qualification(qualified_id, sha, current["workflow_id"])
    return {"source-sha": sha, "artifact-run-id": qualified_id or current_id,
            "reuse-qualified": "true" if qualified_id else "false",
            "builder-artifact-ids": builders, "consumer-artifact-ids": consumers}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--tag", required=True)
    parser.add_argument("--workflow-sha", required=True)
    parser.add_argument("--qualified-run-id", default="")
    args = parser.parse_args()
    if os.environ.get("GITHUB_REPOSITORY") != REPOSITORY:
        raise RuntimeError("Publication is restricted to the PAACS repository")
    result = prepare(args.source.resolve(), args.tag, args.workflow_sha,
                     os.environ.get("GITHUB_RUN_ID", ""), args.qualified_run_id)
    with Path(os.environ["GITHUB_OUTPUT"]).open("a", encoding="utf-8") as stream:
        for key, value in result.items():
            stream.write(f"{key}={value}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
