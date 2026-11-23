#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Qualify a downloaded wheel artifact on a fresh machine without its build SDK."""

import argparse
import hashlib
import json
from pathlib import Path

from build_python_wheel import installed_smoke, occt_pin, sha256, source_revision, target_identity, validate_contents
from python_compatibility import consumer


def consumer_identity(identifier: str) -> dict:
    """Bind independent qualification to its declared interpreter, CPU and OS."""
    selected = consumer(identifier)
    actual = target_identity()
    if actual["target_id"] != selected["target_id"]:
        raise RuntimeError("The consumer does not match this native platform and CPython interpreter")
    required_major = selected.get("required_os_major")
    if required_major is not None:
        try:
            major = int(actual["host"]["macos_version"].split(".")[0])
        except (KeyError, ValueError) as error:
            raise RuntimeError("The macOS consumer must identify its actual operating-system version") from error
        if major != required_major:
            raise RuntimeError(f"The consumer requires macOS {required_major}, found macOS {major}")
    if selected["platform"] == "windows-arm64":
        host = actual["host"]
        version = host["version"].split(".")
        try:
            build = int(version[2])
        except (IndexError, ValueError) as error:
            raise RuntimeError("Windows ARM64 qualification requires an identifiable Windows 11 build") from error
        if build < 22000:
            raise RuntimeError("Windows ARM64 qualification requires Windows 11 or newer")
    return {**actual, "consumer_id": identifier, "consumer_os": selected["consumer_os"]}


def artifact_file(directory: Path, name: str) -> Path:
    if not isinstance(name, str) or Path(name).name != name or name in {".", ".."}:
        raise RuntimeError("Artifact receipts must use plain file names")
    result = directory / name
    if not result.is_file():
        raise RuntimeError(f"The artifact omits {name}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--consumer-id", required=True, help="exact independent consumer ID from the compatibility manifest")
    args = parser.parse_args()
    source, directory = args.source.resolve(), args.artifact_dir.resolve()
    build_dir = (args.build_dir or source / "build/python-wheel-consumer").resolve()
    identity = consumer_identity(args.consumer_id)
    current_platform = identity["platform"]
    receipt = json.loads((directory / "wheel-validation.json").read_text(encoding="utf-8"))
    source_record = source_revision(source)
    if (source_record["clean"] is not True or receipt.get("source_revision") != source_record):
        raise RuntimeError("Independent qualification requires the same clean source revision as its builder")
    if receipt["platform"] != current_platform:
        raise RuntimeError("The downloaded wheel was qualified for a different platform")
    for field in ("target_id", "python_version", "python_tag"):
        if receipt.get(field) != identity[field]:
            raise RuntimeError(f"The downloaded wheel's {field} does not match this consumer")
    wheel = artifact_file(directory, receipt["wheel"])
    if sha256(wheel) != receipt["sha256"]:
        raise RuntimeError("The wheel does not match the builder's receipt")
    validate_contents(wheel)
    pin = occt_pin(source)
    archive_record = receipt["occt_source"]
    if archive_record["file"] != f"occt-{pin['VERSION']}.tar.gz":
        raise RuntimeError("The artifact does not identify the pinned OCCT source")
    archive = artifact_file(directory, archive_record["file"])
    with archive.open("rb") as stream:
        actual_source = hashlib.file_digest(stream, "sha512").hexdigest()
    if actual_source.lower() != pin["SHA512"].lower() or sha256(archive) != archive_record["sha256"]:
        raise RuntimeError("The OCCT source archive does not match its pin and receipt")
    audit = artifact_file(directory, receipt["repair"]["audit"])
    if sha256(audit) != receipt["repair"]["audit_sha256"]:
        raise RuntimeError("The dependency audit does not match the builder's receipt")
    parity_record = receipt["native_parity"]
    parity = artifact_file(directory, parity_record["file"])
    if sha256(parity) != parity_record["sha256"]:
        raise RuntimeError("The native CLI parity result does not match the builder's receipt")

    # This job has checked out source and downloaded artifacts only. No SDK is
    # configured or built; on Linux its mount also differs from the builder.
    numpy_results = installed_smoke(source, wheel, build_dir)
    result = {
        "source_revision": source_record,
        **identity,
        "wheel": wheel.name,
        "sha256": sha256(wheel),
        "installed_tests": json.loads((build_dir / "installed-validation.json").read_text(encoding="utf-8")),
        "numpy_validation": numpy_results,
        "qualification": "passed with isolated Python and cleared SDK, loader and CAD paths",
    }
    (build_dir / "consumer-validation.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
