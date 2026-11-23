#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Shared, fail-closed build and qualification matrix for native Python wheels.

This is a qualification plan, not a claim that an unreleased wheel is supported.
The complete builder and independent consumer sets must pass before publication.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path
import re
import sys


MANIFEST = Path(__file__).resolve().parents[1] / "python/compatibility.json"
_VERSIONS = ("3.12", "3.13", "3.14", "3.15")
_PLATFORMS = {
    "windows-x86_64": ("win_amd64", "windows-2022", "x64"),
    "macos-arm64": ("macosx_15_0_arm64", "macos-15", "arm64"),
    "linux-x86_64": ("manylinux_2_28_x86_64", "ubuntu-24.04", "x64"),
    "linux-arm64": ("manylinux_2_28_aarch64", "ubuntu-24.04-arm", "arm64"),
    "windows-arm64": ("win_arm64", "windows-11-arm", "arm64"),
}
_CONSUMERS = {
    "windows-x86_64": {"windows2022": ("windows-2022", None, False)},
    "macos-arm64": {"macos15": ("macos-15", 15, False),
                    "macos26": ("macos-26", 26, False),
                    "macos27": ("xcode-27", 27, True)},
    "linux-x86_64": {"linuxglibc228": ("ubuntu-24.04", None, False)},
    "linux-arm64": {"linuxglibc228": ("ubuntu-24.04-arm", None, False)},
    "windows-arm64": {"windows11arm": ("windows-11-arm", None, False)},
}


def validate_manifest(value: dict) -> dict:
    """Reject incomplete plans and accidental platform/ABI expansion."""
    root_fields = {"schema", "python_versions", "requires_python", "numpy_current",
                   "numpy_legacy", "platforms"}
    if (not isinstance(value, dict) or set(value) != root_fields
            or value.get("schema") != "paacs-python-compatibility-v1"
            or value.get("python_versions") != list(_VERSIONS)
            or value.get("requires_python") != ">=3.12,<3.16"):
        raise RuntimeError("Compatibility manifest must describe the complete CPython 3.12–3.15 plan")
    if any(not isinstance(value.get(key), str)
           or not re.fullmatch(r"\d+\.\d+\.\d+", value[key])
           for key in ("numpy_current", "numpy_legacy")):
        raise RuntimeError("Compatibility NumPy tests require exact release versions")
    if value["numpy_current"] == value["numpy_legacy"]:
        raise RuntimeError("Legacy and current NumPy tests must be distinct")
    entries = value.get("platforms")
    if not isinstance(entries, list) or len(entries) != len(_PLATFORMS):
        raise RuntimeError("Compatibility manifest must contain all five platform architectures")
    found = set()
    common = {"id", "platform_tag", "runner", "architecture", "jobs", "timeout_minutes",
              "legacy_numpy", "consumers"}
    for row in entries:
        if not isinstance(row, dict) or not isinstance(row.get("id"), str):
            raise RuntimeError("Every compatibility platform must be a named object")
        name = row["id"]
        if name not in _PLATFORMS or name in found:
            raise RuntimeError("Unknown or duplicate compatibility platform: " + name)
        found.add(name)
        tag, runner, architecture = _PLATFORMS[name]
        extra = ({"image", "image_architecture", "minimum_glibc"} if name.startswith("linux-")
                 else {"minimum_macos"} if name == "macos-arm64" else set())
        allowed = common | extra | ({"image_config_sha256"} if name.startswith("linux-") else set())
        if (not common | extra <= row.keys() or not row.keys() <= allowed
                or (row["platform_tag"], row["runner"], row["architecture"]) != (tag, runner, architecture)
                or type(row["jobs"]) is not int or not 1 <= row["jobs"] <= 8
                or type(row["timeout_minutes"]) is not int or not 1 <= row["timeout_minutes"] <= 360
                or type(row["legacy_numpy"]) is not bool
                or row["legacy_numpy"] != (name != "windows-arm64")):
            raise RuntimeError("Unexpected platform build settings: " + name)
        if name.startswith("linux-"):
            arch = "aarch64" if name == "linux-arm64" else "x86_64"
            if (row["minimum_glibc"] != "2.28"
                    or row["image_architecture"] != ("arm64" if arch == "aarch64" else "amd64")
                    or not isinstance(row["image"], str)
                    or not re.fullmatch(r"quay\.io/pypa/manylinux_2_28_" + arch
                                        + r"@sha256:[0-9a-f]{64}", row["image"])
                    or ("image_config_sha256" in row
                        and not re.fullmatch(r"[0-9a-f]{64}", str(row["image_config_sha256"])))):
                raise RuntimeError("Linux builds require a pinned, matching native manylinux 2.28 image")
        if name == "macos-arm64" and row["minimum_macos"] != "15.0":
            raise RuntimeError("macOS ARM64 wheels must retain deployment target 15.0")
        environments = row["consumers"]
        expected = _CONSUMERS[name]
        if not isinstance(environments, list) or len(environments) != len(expected):
            raise RuntimeError("Missing or extra independent consumer environments: " + name)
        seen = set()
        for environment in environments:
            if (not isinstance(environment, dict)
                    or set(environment) != {"id", "runner", "required_os_major", "preview"}
                    or not isinstance(environment.get("id"), str)):
                raise RuntimeError("Malformed independent consumer environment: " + name)
            identity = environment["id"]
            if identity in seen or identity not in expected:
                raise RuntimeError("Unknown or duplicate independent consumer: " + identity)
            seen.add(identity)
            if ((environment["runner"], environment["required_os_major"], environment["preview"])
                    != expected[identity] or type(environment["preview"]) is not bool
                    or (environment["required_os_major"] is not None
                        and type(environment["required_os_major"]) is not int)):
                raise RuntimeError("Independent consumer OS settings differ: " + identity)
    return deepcopy(value)


def load_manifest(path: Path = MANIFEST) -> dict:
    path = Path(path)
    if path.stat().st_size > 64_000:
        raise RuntimeError("Compatibility manifest exceeds its bounded size")
    def unique_object(pairs):
        result = {}
        for key, item in pairs:
            if key in result:
                raise RuntimeError("Duplicate compatibility manifest key: " + key)
            result[key] = item
        return result
    try:
        return validate_manifest(json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique_object))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise RuntimeError("Compatibility manifest must contain UTF-8 JSON") from error


_CONFIG = load_manifest()
PYTHON_VERSIONS = tuple(_CONFIG["python_versions"])
REQUIRES_PYTHON = _CONFIG["requires_python"]
PLATFORM_TAGS = {row["id"]: row["platform_tag"] for row in _CONFIG["platforms"]}


def _version(version: str | None) -> str:
    result = version if version is not None else f"{sys.version_info.major}.{sys.version_info.minor}"
    if result not in PYTHON_VERSIONS:
        raise RuntimeError("Unsupported qualification interpreter: " + str(result))
    return result


def python_tag(version: str | None = None) -> str:
    return "cp" + _version(version).replace(".", "")


def platforms() -> list[dict]:
    """Return independent copies of the five native SDK build configurations."""
    result = []
    for configured in _CONFIG["platforms"]:
        row = deepcopy(configured)
        row["platform"] = row["id"]
        row.setdefault("image", "")
        row["sdk_artifact"] = "paacs-mesh-occt-" + row["id"]
        result.append(row)
    return result


def numpy_versions(platform: str, python_version: str | None = None) -> tuple[str, ...]:
    version = _version(python_version)
    row = next((p for p in _CONFIG["platforms"] if p["id"] == platform), None)
    if row is None:
        raise RuntimeError("Unknown qualification platform: " + str(platform))
    if row["legacy_numpy"] and version == "3.12":
        return (_CONFIG["numpy_legacy"], _CONFIG["numpy_current"])
    return (_CONFIG["numpy_current"],)


def targets() -> list[dict]:
    result = []
    for platform in platforms():
        for version in PYTHON_VERSIONS:
            row = {key: deepcopy(value) for key, value in platform.items() if key != "consumers"}
            tag = python_tag(version)
            row.update(id=f"{tag}-{platform['id']}", python=version, python_tag=tag,
                       artifact=f"paacs-mesh-{tag}-{platform['platform_tag']}",
                       numpy_versions=list(numpy_versions(platform["id"], version)))
            row["job_name"] = f"wheel ({row['id']})"
            result.append(row)
    return result


def target(platform: str, python_version: str | None = None) -> dict:
    version = _version(python_version)
    for row in targets():
        if row["platform"] == platform and row["python"] == version:
            return row
    raise RuntimeError("Unknown qualification platform: " + str(platform))


def target_by_id(identity: str) -> dict:
    for row in targets():
        if row["id"] == identity:
            return row
    raise RuntimeError("Unknown qualification target: " + str(identity))


def consumers() -> list[dict]:
    result = []
    configurations = {p["id"]: p["consumers"] for p in _CONFIG["platforms"]}
    for build in targets():
        for environment in configurations[build["platform"]]:
            row = deepcopy(build)
            row.update(id=build["id"] + "-" + environment["id"], target_id=build["id"],
                       consumer_os=environment["id"], runner=environment["runner"],
                       required_os_major=environment["required_os_major"], preview=environment["preview"])
            row["consumer_artifact"] = "paacs-mesh-consumer-" + row["id"]
            row["job_name"] = f"installed wheel ({row['id']})"
            result.append(row)
    return result


def consumer(identity: str) -> dict:
    for row in consumers():
        if row["id"] == identity:
            return row
    raise RuntimeError("Unknown qualification consumer: " + str(identity))


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    for kind in ("targets", "consumers", "platforms"):
        selection.add_argument("--" + kind, action="store_const", const=kind, dest="kind")
    parser.add_argument("--github-output", type=Path,
                        help="Also append the selected named matrix to this GitHub Actions output file")
    args = parser.parse_args(argv)
    rows = {"targets": targets, "consumers": consumers, "platforms": platforms}[args.kind]()
    encoded = json.dumps({"include": rows}, separators=(",", ":"), allow_nan=False)
    if args.github_output is not None:
        with args.github_output.open("a", encoding="utf-8") as stream:
            stream.write(args.kind + "=" + encoded + "\n")
    print(encoded)


if __name__ == "__main__":
    main()
