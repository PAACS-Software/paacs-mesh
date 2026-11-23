#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Compare an installed Python STEP wheel with the same platform's native CLI."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import paacs_mesh as pm


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-cli", type=Path, required=True)
    parser.add_argument("--sdk-lib-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    native = args.native_cli.resolve(strict=True)
    libraries = args.sdk_lib_dir.resolve(strict=True)
    source = Path(__file__).resolve().parents[1]
    fixtures = source / "tests/fixtures/cad"
    if "site-packages" not in str(Path(pm.__file__).resolve()).lower():
        raise RuntimeError("Parity must exercise the installed wheel")

    # Only the native executable receives the SDK loader path. Python has
    # already loaded the repaired wheel under the caller's isolated settings.
    native_env = dict(os.environ)
    variable = "PATH" if sys.platform == "win32" else (
        "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH")
    native_env[variable] = str(libraries) + os.pathsep + native_env.get(variable, "")
    common = ["--native-cad", "--report-json", "--recovery-mode", "standard",
              "--threads", "1", "--native-validation-threads", "1"]
    cases = [
        ("cylinder", "cylinder.step", .3, {}, ["--independent-parts"]),
        ("assembly", "assembly.step", .3, {}, ["--independent-parts"]),
        ("combined-sizing", "cylinder.step", .7, {
            "assembly": "conforming", "native_sizing": True,
            "face_sizes": [pm.FaceSize(1, .16)],
            "volume_controls": [pm.BoxSize((-.3, -.3, .7), (.3, .3, 1.3), .15)],
            "min_size": .12, "growth": 1.55, "max_work": 100_000_000,
        }, ["--native-sizing", "--growth", "1.55", "--face-size", "1", ".16",
            "--size-box", "-.3", "-.3", ".7", ".3", ".3", "1.3", ".15",
            "--min-size", ".12", "--max-work", "100000000"]),
        ("construction-caps", "cylinder.step", .5, {
            "assembly": "conforming", "construction_effort": "thorough",
            "max_surface_points": 400, "max_surface_chord_points": 128,
            "max_interior_candidates": 20000, "max_interior_sites": 256,
            "max_work": 100_000_000,
        }, ["--construction-effort", "thorough", "--max-surface-points", "400",
            "--max-surface-chord-points", "128", "--max-interior-candidates", "20000",
            "--max-interior-sites", "256", "--max-work", "100000000"]),
    ]
    results = []
    with tempfile.TemporaryDirectory(prefix="paacs-native-parity-") as directory:
        root = Path(directory)
        for name, filename, size, options, flags in cases:
            step = fixtures / filename
            mesh = pm.read_step(step).mesh(size, workers=1, validation_workers=1, **options)
            report = mesh.report
            if not report["audit"]["valid"] or not report["quality"]["valid"]:
                raise RuntimeError(f"Python native checks failed for {name}")
            python_archive = root / f"{name}-python.paacsmesh"
            native_archive = root / f"{name}-native.paacsmesh"
            mesh.write_archive(python_archive)
            child = subprocess.run(
                [str(native), str(step), str(size), str(native_archive), *common, *flags],
                env=native_env, capture_output=True, text=True,
                encoding="utf-8", errors="replace", timeout=120,
            )
            if child.returncode:
                raise RuntimeError(f"Native CLI failed for {name} (exit {child.returncode}, "
                                   f"0x{child.returncode & 0xffffffff:08x}):\n"
                                   f"{child.stdout}\n{child.stderr}")
            python_hash, native_hash = digest(python_archive), digest(native_archive)
            if python_hash != native_hash:
                raise RuntimeError(f"Python and native CLI archives differ for {name}")
            if name == "construction-caps":
                native_reports = []
                for line in child.stdout.splitlines():
                    try:
                        value = json.loads(line)
                    except ValueError:
                        continue
                    if isinstance(value, dict) and value.get("schema") == "paacs-native-cad-report-v1":
                        native_reports.append(value)
                requested = {key: value for key, value in report["construction"].items() if key != "scope"}
                if len(native_reports) != 1 or native_reports[0].get("construction") != requested:
                    raise RuntimeError("Python and native CLI construction receipts differ")
            results.append({
                "case": name, "source": filename, "source_sha256": digest(step),
                "archives_identical": True, "archive_sha256": python_hash,
                "nodes": mesh.node_count, "cells": mesh.cell_count,
                "tetrahedra": mesh.tetrahedron_count,
                "face_tags": mesh.face_tags, "body_tags": mesh.body_tags,
                "volume": report["quality"]["total_volume"],
                "mesh_scope": report["mesh_scope"],
            })
    receipt = {"schema": "paacs-python-native-parity-v1", "status": "passed",
               "platform": sys.platform, "engine_version": pm.native_version,
               "workers": 1, "validation_workers": 1, "cases": results}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
