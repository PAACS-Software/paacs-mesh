#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Prepare qualified wheels and matching source; PyPI upload remains separate."""
from __future__ import annotations

import argparse
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import zipfile

from packaging.specifiers import SpecifierSet
from packaging.utils import canonicalize_name, parse_wheel_filename
from packaging.version import Version

from build_python_wheel import NUMPY_TEST_VERSIONS, PLATFORM_TAGS, occt_pin, sha256


ARTIFACTS = {
    "windows-x86_64": "paacs-mesh-cp312-win_amd64",
    "macos-arm64": "paacs-mesh-cp312-macosx_15_0_arm64",
    "linux-x86_64": "paacs-mesh-cp312-manylinux_2_28_x86_64",
}


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path: Path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def revision(source: Path):
    def git(*args):
        return subprocess.check_output(["git", "-c", "safe.directory=" + source.resolve().as_posix(),
            "-C", str(source), *args], text=True).strip()
    return {"commit": git("rev-parse", "HEAD"), "tree": git("rev-parse", "HEAD^{tree}"),
            "clean": not git("status", "--porcelain", "--untracked-files=normal")}


def expected_revision(value: str):
    if not re.fullmatch(r"[0-9a-f]{40}", value):
        raise RuntimeError("Release source must be a full immutable Git commit")


def preflight(source: Path, tag: str, expected: str):
    expected_revision(expected)
    version = (source / "VERSION.txt").read_text(encoding="utf-8").strip()
    if tag != "v" + version:
        raise RuntimeError("Release tag must exactly match v plus VERSION.txt")
    Version(version)
    actual = revision(source)
    tagged = subprocess.check_output(["git", "-C", str(source), "rev-parse", "--verify",
        f"refs/tags/{tag}^{{commit}}"], text=True).strip()
    if actual["commit"] != expected or tagged != expected or not actual["clean"]:
        raise RuntimeError("Release tag, workflow source and clean checkout must identify the same commit")
    return {"version": str(Version(version)), "tag": tag, **actual}


def artifact_file(directory: Path, name: str):
    if not isinstance(name, str) or Path(name).name != name or name in {".", ".."}:
        raise RuntimeError("Artifact receipts must name a plain file")
    path = directory / name
    if not path.is_file():
        raise RuntimeError(f"Missing release artifact: {name}")
    return path


def checked_file(directory: Path, record: dict):
    path = artifact_file(directory, record["file"])
    if sha256(path) != record["sha256"]:
        raise RuntimeError(f"Artifact checksum differs from qualification: {path.name}")
    return path


def tests_passed(receipt: dict):
    fields = [receipt.get(key) for key in ("discovered", "passed", "skipped")]
    if (any(type(v) is not int or v < 0 for v in fields) or fields[1] == 0
            or fields[0] != fields[1] + fields[2]):
        raise RuntimeError("The artifact lacks a nonempty successful installed-test result")


def numpy_passed(rows: list):
    if (not isinstance(rows, list) or len(rows) != len(NUMPY_TEST_VERSIONS)
            or {row.get("version") for row in rows} != set(NUMPY_TEST_VERSIONS)
            or any(row.get("status") != "passed" or type(row.get("tests")) is not int
                   or row["tests"] <= 0 for row in rows)):
        raise RuntimeError("Both supported NumPy compatibility tests must pass")


def check_source(receipt: dict, expected: str):
    source = receipt.get("source_revision", {})
    if source.get("commit") != expected or source.get("clean") is not True:
        raise RuntimeError("Artifact qualification came from another or dirty source revision")


def wheel_metadata(path: Path, version: Version, platform: str):
    name, parsed, build, tags = parse_wheel_filename(path.name)
    if (canonicalize_name(name) != "paacs-mesh" or parsed != version or build
            or not tags or any(t.interpreter != "cp312" or t.abi != "cp312"
                               or t.platform != PLATFORM_TAGS[platform] for t in tags)):
        raise RuntimeError("Wheel name, version or CPython/platform tags differ from the qualified release")
    with zipfile.ZipFile(path) as archive:
        metadata = [n for n in archive.namelist() if n.endswith(".dist-info/METADATA")]
        if len(metadata) != 1 or archive.getinfo(metadata[0]).file_size > 1_000_000:
            raise RuntimeError("Wheel must contain one bounded metadata record")
        record = BytesParser().parsebytes(archive.read(metadata[0]))
        if (canonicalize_name(record["Name"] or "") != "paacs-mesh"
                or Version(record["Version"]) != version
                or SpecifierSet(record["Requires-Python"] or "") != SpecifierSet(">=3.12,<3.13")):
            raise RuntimeError("Wheel metadata must identify the CPython 3.12 PAACS release")


def assemble(source: Path, builders: Path, consumers: Path, expected: str, output: Path):
    expected_revision(expected)
    if output.exists():
        raise RuntimeError("Refusing to overwrite an existing release staging directory")
    version = Version((source / "VERSION.txt").read_text(encoding="utf-8").strip())
    pin = occt_pin(source)
    selected, evidence, archives = [], {}, {}
    for platform, artifact in ARTIFACTS.items():
        directory = builders / artifact
        receipt_path = artifact_file(directory, "wheel-validation.json")
        receipt = load(receipt_path)
        check_source(receipt, expected)
        if receipt.get("platform") != platform or not str(receipt.get("installed_smoke", "")).startswith("passed "):
            raise RuntimeError("Builder did not qualify this platform")
        tests_passed(receipt.get("installed_tests", {}))
        numpy_passed(receipt.get("numpy_validation"))
        wheel = checked_file(directory, {"file": receipt["wheel"], "sha256": receipt["sha256"]})
        if list(directory.glob("*.whl")) != [wheel]:
            raise RuntimeError("Builder artifact must contain exactly its qualified wheel")
        wheel_metadata(wheel, version, platform)
        consumer_path = artifact_file(consumers / ("paacs-mesh-installed-" + platform), "consumer-validation.json")
        consumer = load(consumer_path)
        check_source(consumer, expected)
        if (consumer.get("platform") != platform or consumer.get("wheel") != wheel.name
                or consumer.get("sha256") != receipt["sha256"]
                or not str(consumer.get("qualification", "")).startswith("passed ")):
            raise RuntimeError("The independent consumer did not qualify this exact wheel")
        tests_passed(consumer.get("installed_tests", {}))
        numpy_passed(consumer.get("numpy_validation"))
        archive = checked_file(directory, receipt["occt_source"])
        if archive.name != f"occt-{pin['VERSION']}.tar.gz":
            raise RuntimeError("OCCT source version does not match the build recipe")
        with archive.open("rb") as stream:
            if hashlib.file_digest(stream, "sha512").hexdigest() != pin["SHA512"].lower():
                raise RuntimeError("OCCT source differs from the pinned upstream source")
        archives[archive.name] = archive
        audit = checked_file(directory, {"file": receipt["repair"]["audit"], "sha256": receipt["repair"]["audit_sha256"]})
        parity = checked_file(directory, receipt["native_parity"])
        parity_result = load(parity)
        if parity_result.get("status") != "passed" or len(parity_result.get("cases", [])) != 3:
            raise RuntimeError("Native/Python parity qualification is missing")
        evidence[platform] = {"builder": receipt_path, "consumer": consumer_path, "parity": parity, "dependencies": audit}
        selected.append(wheel)
    # Validation finishes before staging; a failed gate cannot leave an upload set.
    dist, assets = output / "dist", output / "release-assets"
    dist.mkdir(parents=True)
    assets.mkdir()
    for wheel in selected:
        shutil.copyfile(wheel, dist / wheel.name)
    for name, path in archives.items():
        shutil.copyfile(path, assets / name)
    for platform, records in evidence.items():
        for kind, path in records.items():
            shutil.copyfile(path, assets / f"{platform}-{kind}{path.suffix}")
    manifest = {"schema": "paacs-python-release-v1", "source_commit": expected, "version": str(version),
        "wheels": {p.name: sha256(p) for p in selected},
        "source_and_receipts": {p.name: sha256(p) for p in assets.iterdir()},
        "qualification": "Three repaired platform wheels and their separate installed consumers; same clean immutable source"}
    write(output / "release-manifest.json", manifest)
    shutil.copyfile(output / "release-manifest.json", assets / "release-manifest.json")
    return manifest


def verify(manifest: dict, directory: Path, expected: str):
    expected_revision(expected)
    if (manifest.get("schema") != "paacs-python-release-v1" or manifest.get("source_commit") != expected
            or len(manifest.get("wheels", {})) != 3):
        raise RuntimeError("Publication manifest has the wrong source or platform set")
    files = {p.name for p in directory.iterdir() if p.is_file()}
    if files != set(manifest["wheels"]):
        raise RuntimeError("Publication directory must contain exactly the qualified wheels")
    for name, digest in manifest["wheels"].items():
        if sha256(artifact_file(directory, name)) != digest:
            raise RuntimeError("Publication wheel changed after assembly")


def publish_sources(directory: Path, tag: str, expected: str):
    """Keep exact matching source assets available before any PyPI upload."""
    expected_revision(expected)
    manifest = load(artifact_file(directory, "release-manifest.json"))
    if manifest.get("source_commit") != expected:
        raise RuntimeError("Source assets do not identify the qualified revision")
    records = manifest["source_and_receipts"]
    for name, digest in records.items():
        if sha256(artifact_file(directory, name)) != digest:
            raise RuntimeError("Source assets changed after assembly")
    names = set(records) | {"release-manifest.json"}
    if {p.name for p in directory.iterdir() if p.is_file()} != names:
        raise RuntimeError("Source staging contains unexpected release assets")
    repository = "PAACS-Software/paacs-mesh"
    view = json.loads(subprocess.check_output(["gh", "release", "view", tag, "--repo", repository,
        "--json", "tagName,isDraft,assets"], text=True))
    if view["isDraft"] or view["tagName"] != tag:
        raise RuntimeError("The matching GitHub release must be public before publication")
    existing = {asset["name"] for asset in view["assets"]}
    for name in sorted(names):
        path = artifact_file(directory, name)
        if name in existing:
            with tempfile.TemporaryDirectory(prefix="paacs-release-source-") as temporary:
                subprocess.run(["gh", "release", "download", tag, "--repo", repository,
                    "--pattern", name, "--dir", temporary], check=True)
                if sha256(Path(temporary) / name) != sha256(path):
                    raise RuntimeError(f"Published asset differs; refusing to replace {name}")
        else:
            subprocess.run(["gh", "release", "upload", tag, "--repo", repository, str(path)], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("preflight", "assemble", "verify", "source-assets"))
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--expected-sha", required=True)
    parser.add_argument("--tag")
    parser.add_argument("--builders", type=Path)
    parser.add_argument("--consumers", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--dist", type=Path)
    args = parser.parse_args()
    if args.command == "preflight":
        if args.tag is None:
            parser.error("preflight requires --tag")
        print(json.dumps(preflight(args.source, args.tag, args.expected_sha)))
    elif args.command == "assemble":
        if None in (args.builders, args.consumers, args.output):
            parser.error("assemble requires --builders, --consumers and --output")
        print(json.dumps(assemble(args.source, args.builders, args.consumers, args.expected_sha, args.output), indent=2))
    elif args.command == "verify":
        if args.manifest is None or args.dist is None:
            parser.error("verify requires --manifest and --dist")
        verify(load(args.manifest), args.dist, args.expected_sha)
        print("Release wheel checksums verified")
    else:
        if args.dist is None or args.tag is None:
            parser.error("source-assets requires --tag and --dist")
        publish_sources(args.dist, args.tag, args.expected_sha)


if __name__ == "__main__":
    main()
