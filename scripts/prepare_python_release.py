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

from build_python_wheel import occt_pin, sha256
from python_compatibility import PLATFORM_TAGS, REQUIRES_PYTHON, consumers as compatibility_consumers, targets


TARGETS = {row["id"]: row for row in targets()}
CONSUMERS = {row["id"]: row for row in compatibility_consumers()}
ARTIFACTS = {name: row["artifact"] for name, row in TARGETS.items()}


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


def numpy_passed(rows: list, expected):
    if (not isinstance(rows, list) or len(rows) != len(expected)
            or any(not isinstance(row, dict) for row in rows)
            or {row.get("version") for row in rows} != set(expected)
            or any(row.get("status") != "passed" or type(row.get("tests")) is not int
                   or row["tests"] <= 0 for row in rows)):
        raise RuntimeError("Every NumPy version required for this Python/platform must pass")


def check_source(receipt: dict, expected: str):
    source = receipt.get("source_revision", {})
    if source.get("commit") != expected or source.get("clean") is not True:
        raise RuntimeError("Artifact qualification came from another or dirty source revision")


def wheel_metadata(path: Path, version: Version, platform: str, python_version="3.12"):
    name, parsed, build, tags = parse_wheel_filename(path.name)
    expected_platform = PLATFORM_TAGS[platform]
    allowed_platforms = {expected_platform}
    if platform in ("linux-x86_64", "linux-arm64"):
        architecture = "x86_64" if platform == "linux-x86_64" else "aarch64"
        allowed_platforms.add("manylinux_2_27_" + architecture)
    python_tag = "cp" + python_version.replace(".", "")
    platforms = {t.platform for t in tags}
    if (canonicalize_name(name) != "paacs-mesh" or parsed != version or build
            or not tags or any(t.interpreter != python_tag or t.abi != python_tag
                               for t in tags)
            or expected_platform not in platforms or not platforms <= allowed_platforms):
        raise RuntimeError("Wheel name, version or CPython/platform tags differ from the qualified release")
    with zipfile.ZipFile(path) as archive:
        metadata = [n for n in archive.namelist() if n.endswith(".dist-info/METADATA")]
        if len(metadata) != 1 or archive.getinfo(metadata[0]).file_size > 1_000_000:
            raise RuntimeError("Wheel must contain one bounded metadata record")
        record = BytesParser().parsebytes(archive.read(metadata[0]))
        if (canonicalize_name(record["Name"] or "") != "paacs-mesh"
                or Version(record["Version"]) != version
                or SpecifierSet(record["Requires-Python"] or "") != SpecifierSet(REQUIRES_PYTHON)):
            raise RuntimeError("Wheel metadata must identify the qualified PAACS Python release")


def receipt_target(receipt: dict, row: dict, *, consumer=False):
    expected = {"target_id": row.get("target_id", row["id"]), "platform": row["platform"],
                "python_version": row["python"], "python_tag": row["python_tag"]}
    if consumer:
        expected.update(consumer_id=row["id"], consumer_os=row["consumer_os"])
    if any(receipt.get(key) != value for key, value in expected.items()):
        raise RuntimeError("Qualification receipt does not identify its exact Python/platform consumer")
    host = receipt.get("host", {})
    if not isinstance(host, dict):
        raise RuntimeError("Qualification receipt lacks actual host evidence")
    platform = row["platform"]
    system = "Windows" if platform.startswith("windows-") else "Darwin" if platform == "macos-arm64" else "Linux"
    machines = {"AMD64", "x86_64", "x64"} if platform.endswith("x86_64") else {"ARM64", "arm64", "aarch64"}
    if host.get("system") != system or host.get("machine") not in machines:
        raise RuntimeError("Qualification host operating system or architecture differs from its target")
    if system == "Darwin":
        actual = str(host.get("macos_version", "")).split(".")[0]
        required = row.get("required_os_major") if consumer else 15
        if not actual.isdigit() or int(actual) < 15 or (consumer and int(actual) != required):
            raise RuntimeError("macOS consumer did not execute on its required OS major")
    elif platform == "windows-arm64":
        build = str(host.get("version", "")).split(".")
        windows11 = host.get("release") == "11" or (
            len(build) >= 3 and all(value.isdigit() for value in build[:3])
            and int(build[0]) >= 10 and int(build[2]) >= 22000)
        if not windows11:
            raise RuntimeError("Windows ARM64 qualification requires Windows 11 host evidence")
    elif system == "Linux":
        libc = str(host.get("libc", "")).lower()
        version = str(host.get("libc_version", ""))
        if libc not in ("glibc", "gnu libc") or not re.fullmatch(r"[0-9]+\.[0-9]+(?:\.[0-9]+)*", version):
            raise RuntimeError("Linux qualification requires glibc host evidence")
        if tuple(int(value) for value in version.split(".")[:2]) < (2, 28):
            raise RuntimeError("Linux qualification requires glibc 2.28 or newer")


def artifact_directories(directory: Path, expected: set[str]):
    records = list(directory.iterdir())
    if any(not path.is_dir() for path in records) or {path.name for path in records} != expected:
        raise RuntimeError("Release artifact directories do not match the complete compatibility matrix")


def check_build_info(directory: Path, receipt: dict, row: dict, version: str):
    """Bind release-only compiler/configuration evidence to the qualified wheel."""
    record = receipt.get("build_info")
    if (not isinstance(record, dict) or record.get("file") != "native-build-info.txt"
            or not re.fullmatch(r"[0-9a-f]{64}", str(record.get("sha256")))):
        raise RuntimeError("The wheel lacks hash-bound native compiler and configuration evidence")
    path = checked_file(directory, record)
    if path.stat().st_size > 65_536:
        raise RuntimeError("Native build information exceeds its metadata budget")
    fields = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key in fields:
            raise RuntimeError("Native build information contains duplicate evidence fields")
        fields[key] = value
    expected = {"release": version, "cmake": "4.4.3", "standard": "23", "configuration": "Release",
                "strict_fp": "ON", "ipo": "OFF", "native": "OFF", "sanitizer": "none"}
    if any(fields.get(key) != value for key, value in expected.items()) or "flags" not in fields:
        raise RuntimeError("The native wheel build does not match its qualified release configuration")
    if any(flag in fields["flags"] for flag in ("/fp:fast", "-ffast-math", "-Ofast", "-march=native",
                                               "-funsafe-math-optimizations")):
        raise RuntimeError("The native wheel build enables unqualified floating-point or CPU flags")
    compiler = re.fullmatch(r"(MSVC|GNU|Clang|AppleClang) ([0-9]+(?:\.[0-9]+)+)", fields.get("compiler", ""))
    platform = row["platform"]
    permitted = {"MSVC"} if platform.startswith("windows-") else {"AppleClang"} if platform == "macos-arm64" else {"GNU", "Clang"}
    if compiler is None or compiler[1] not in permitted:
        raise RuntimeError("The actual native compiler does not match its wheel platform")
    if platform.startswith("windows-"):
        toolchain = receipt.get("toolchain", {})
        if not isinstance(toolchain, dict):
            raise RuntimeError("The Windows wheel omits its selected toolchain evidence")
        tools = re.fullmatch(r"14\.([0-9]+)(?:\.[0-9]+)+", str(toolchain.get("tools_version", "")))
        native = compiler[2].split(".")
        if tools is None or native[:2] != ["19", tools[1]]:
            raise RuntimeError("Selected Windows tools and actual native compiler versions differ")
        requested = row.get("msvc_toolset")
        if requested and (toolchain.get("requested_toolset") != requested or "14." + tools[1] != requested):
            raise RuntimeError("The actual Windows ARM64 compiler does not match its pinned toolset")
    return path


def assemble(source: Path, builders: Path, consumers: Path, expected: str, output: Path):
    expected_revision(expected)
    if output.exists():
        raise RuntimeError("Refusing to overwrite an existing release staging directory")
    source_version = (source / "VERSION.txt").read_text(encoding="utf-8").strip()
    version = Version(source_version)
    pin = occt_pin(source)
    artifact_directories(builders, set(ARTIFACTS.values()))
    artifact_directories(consumers, {row["consumer_artifact"] for row in CONSUMERS.values()})
    selected, evidence, archives, qualified = [], {}, {}, {}
    for identifier, row in TARGETS.items():
        platform, artifact = row["platform"], row["artifact"]
        directory = builders / artifact
        receipt_path = artifact_file(directory, "wheel-validation.json")
        receipt = load(receipt_path)
        check_source(receipt, expected)
        receipt_target(receipt, row)
        if receipt.get("platform") != platform or not str(receipt.get("installed_smoke", "")).startswith("passed "):
            raise RuntimeError("Builder did not qualify this platform")
        tests_passed(receipt.get("installed_tests", {}))
        numpy_passed(receipt.get("numpy_validation"), row["numpy_versions"])
        wheel = checked_file(directory, {"file": receipt["wheel"], "sha256": receipt["sha256"]})
        if list(directory.glob("*.whl")) != [wheel]:
            raise RuntimeError("Builder artifact must contain exactly its qualified wheel")
        wheel_metadata(wheel, version, platform, row["python"])
        build_info = check_build_info(directory, receipt, row, source_version)
        target_consumers = [item for item in CONSUMERS.values() if item["target_id"] == identifier]
        consumer_evidence = {}
        for consumer_row in target_consumers:
            consumer_path = artifact_file(consumers / consumer_row["consumer_artifact"], "consumer-validation.json")
            consumer = load(consumer_path)
            check_source(consumer, expected)
            receipt_target(consumer, consumer_row, consumer=True)
            if (consumer.get("wheel") != wheel.name or consumer.get("sha256") != receipt["sha256"]
                    or consumer.get("build_info") != receipt["build_info"]
                    or consumer.get("toolchain") != receipt.get("toolchain")
                    or not str(consumer.get("qualification", "")).startswith("passed ")):
                raise RuntimeError("The independent consumer did not qualify this exact wheel")
            tests_passed(consumer.get("installed_tests", {}))
            numpy_passed(consumer.get("numpy_validation"), row["numpy_versions"])
            consumer_evidence[consumer_row["id"]] = consumer_path
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
        evidence[identifier] = {"builder": receipt_path, "parity": parity, "dependencies": audit,
                                "native-build-info": build_info}
        for consumer_id, consumer_path in consumer_evidence.items():
            evidence[consumer_id] = {"consumer": consumer_path}
        qualified[identifier] = {"wheel": wheel.name, "platform": platform, "python": row["python"],
                                "python_tag": row["python_tag"], "consumer_ids": sorted(consumer_evidence)}
        selected.append(wheel)
    if len({wheel.name for wheel in selected}) != len(TARGETS):
        raise RuntimeError("Release targets must identify distinct Python/platform wheels")
    # Validation finishes before staging; a failed gate cannot leave an upload set.
    dist, assets = output / "dist", output / "release-assets"
    dist.mkdir(parents=True)
    assets.mkdir()
    for wheel in selected:
        shutil.copyfile(wheel, dist / wheel.name)
    for name, path in archives.items():
        shutil.copyfile(path, assets / name)
    for identifier, records in evidence.items():
        for kind, path in records.items():
            shutil.copyfile(path, assets / f"{identifier}-{kind}{path.suffix}")
    manifest = {"schema": "paacs-python-release-v1", "source_commit": expected, "version": str(version),
        "wheels": {p.name: sha256(p) for p in selected},
        "source_and_receipts": {p.name: sha256(p) for p in assets.iterdir()},
        "targets": qualified, "qualified_consumers": sorted(CONSUMERS),
        "qualification": "Complete Python/platform wheels and independent OS consumers; same clean immutable source"}
    write(output / "release-manifest.json", manifest)
    shutil.copyfile(output / "release-manifest.json", assets / "release-manifest.json")
    return manifest


def verify(manifest: dict, directory: Path, expected: str):
    expected_revision(expected)
    if (manifest.get("schema") != "paacs-python-release-v1" or manifest.get("source_commit") != expected
            or len(manifest.get("wheels", {})) != len(TARGETS)):
        raise RuntimeError("Publication manifest has the wrong source or platform set")
    verify_matrix(manifest)
    files = {p.name for p in directory.iterdir() if p.is_file()}
    if files != set(manifest["wheels"]):
        raise RuntimeError("Publication directory must contain exactly the qualified wheels")
    for name, digest in manifest["wheels"].items():
        if sha256(artifact_file(directory, name)) != digest:
            raise RuntimeError("Publication wheel changed after assembly")


def verify_matrix(manifest: dict):
    """Bind every wheel and OS consumer to the same exact compatibility contract."""
    if (set(manifest.get("targets", {})) != set(TARGETS)
            or manifest.get("qualified_consumers") != sorted(CONSUMERS)):
        raise RuntimeError("Publication manifest lacks the complete Python/platform/OS matrix")
    wheels = manifest.get("wheels", {})
    if not isinstance(wheels, dict) or len(wheels) != len(TARGETS):
        raise RuntimeError("Publication manifest lacks the exact qualified wheel set")
    selected = set()
    version = Version(manifest["version"])
    for identifier, row in TARGETS.items():
        record = manifest["targets"][identifier]
        expected_consumers = sorted(name for name, item in CONSUMERS.items() if item["target_id"] == identifier)
        if (record.get("platform") != row["platform"] or record.get("python") != row["python"]
                or record.get("python_tag") != row["python_tag"] or record.get("consumer_ids") != expected_consumers):
            raise RuntimeError("Publication target differs from its required Python or OS consumers")
        name = record.get("wheel")
        if not isinstance(name, str) or Path(name).name != name or name not in wheels or name in selected:
            raise RuntimeError("Publication targets do not map to distinct qualified wheels")
        _, parsed, build, tags = parse_wheel_filename(name)
        required_platform = row["platform_tag"]
        allowed = {required_platform}
        if row["platform"].startswith("linux-"):
            allowed.add(required_platform.replace("manylinux_2_28_", "manylinux_2_27_"))
        if (parsed != version or build or not tags
                or any(tag.interpreter != row["python_tag"] or tag.abi != row["python_tag"] for tag in tags)
                or required_platform not in {tag.platform for tag in tags}
                or not {tag.platform for tag in tags} <= allowed
                or canonicalize_name(parse_wheel_filename(name)[0]) != "paacs-mesh"
                or not re.fullmatch(r"[0-9a-f]{64}", str(wheels[name]))):
            raise RuntimeError("Publication wheel tags or checksum differ from the target contract")
        selected.add(name)
    if selected != set(wheels):
        raise RuntimeError("Publication manifest has unqualified wheels")


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
