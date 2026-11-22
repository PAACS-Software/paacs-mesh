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

from build_python_wheel import host_platform, installed_smoke, occt_pin, sha256, source_revision, validate_contents


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
    args = parser.parse_args()
    source, directory = args.source.resolve(), args.artifact_dir.resolve()
    build_dir = (args.build_dir or source / "build/python-wheel-consumer").resolve()
    current_platform = host_platform()
    receipt = json.loads((directory / "wheel-validation.json").read_text(encoding="utf-8"))
    source_record = source_revision(source)
    if (source_record["clean"] is not True or receipt.get("source_revision") != source_record):
        raise RuntimeError("Independent qualification requires the same clean source revision as its builder")
    if receipt["platform"] != current_platform:
        raise RuntimeError("The downloaded wheel was qualified for a different platform")
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
        "platform": current_platform,
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
