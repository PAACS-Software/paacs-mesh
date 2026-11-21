#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Compare Python STEP sizing controls; optionally plot actual mesh sections."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

import paacs_mesh as pm


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="new output directory")
    parser.add_argument("--step", type=Path, default=root / "tests/fixtures/cad/cylinder.step")
    parser.add_argument("--size", type=float, default=.7)
    parser.add_argument("--face-tag", type=int, default=1)
    parser.add_argument("--face-size", type=float, default=.16)
    parser.add_argument("--box", nargs=7, type=float, default=[-.3, -.3, .7, .3, .3, 1.3, .15])
    parser.add_argument("--min-size", type=float, default=.12)
    parser.add_argument("--growth", type=float, default=1.55)
    parser.add_argument("--max-work", type=int, default=100_000_000)
    parser.add_argument("--plot", action="store_true", help="requires NumPy and Matplotlib")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists; choose a new directory")
    model = pm.read_step(args.step)
    face = pm.FaceSize(tag=args.face_tag, size=args.face_size)
    box = pm.BoxSize(lower=tuple(args.box[:3]), upper=tuple(args.box[3:6]), size=args.box[6])
    args.output.mkdir(parents=True)
    source = args.output / "source.step"
    shutil.copyfile(args.step, source)
    recipe = {
        "schema": "paacs-native-step-sizing-v1", "source": source.name,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "global_mm": args.size, "minimum_mm": args.min_size, "growth": args.growth,
        "face_tag": args.face_tag, "face_mm": args.face_size, "box_mm": args.box,
        "chord_mm": None,
    }
    save_json(args.output / "project.json", recipe)
    rows = []
    for variant in ("global", "face", "box", "combined"):
        directory = args.output / variant
        directory.mkdir()
        mesh = model.mesh(args.size, native_sizing=True, assembly="conforming",
                          min_size=args.min_size, growth=args.growth,
                          face_sizes=[face] if variant in ("face", "combined") else [],
                          volume_controls=[box] if variant in ("box", "combined") else [],
                          max_work=args.max_work, workers=1, validation_workers=1)
        mesh.write(directory / "mesh.paacsmesh")
        mesh.write(directory / "mesh.vtk")
        save_json(directory / "report.json", mesh.report)
        rows.append({"variant": variant, "nodes": mesh.node_count,
                     "tetrahedra": mesh.tetrahedron_count,
                     "audit_valid": mesh.report["audit"]["valid"]})
    if args.plot:
        # Plotting is a repository tool, not a dependency of the installed mesher.
        sys.path.insert(0, str(root))
        from scripts.native_step_sizing import measure, plot_sections
        from scripts.mesh_visual_data import read_archive

        for row in rows:
            directory = args.output / row["variant"]
            metrics, _ = measure(read_archive(directory / "mesh.paacsmesh"), recipe)
            save_json(directory / "measurements.json", metrics)
            row["measurements"] = metrics
        plot_sections(args.output, args.output, 1, 0.)
    save_json(args.output / "results.json", rows)
    print(json.dumps(rows, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
