# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Mesh a cantilever STEP with PAACS and solve three or more size levels."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from solve import LENGTH, WIDTH, HEIGHT, sha256, solve


def read_report(stdout):
    for line in reversed(stdout.splitlines()):
        if line.startswith("{"):
            result = json.loads(line)
            if result.get("schema") == "paacs-native-cad-report-v1":
                return result
    raise ValueError("PAACS did not emit a native CAD receipt")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", required=True, type=Path)
    parser.add_argument("--mesher", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--size", nargs="+", type=float, default=[5., 2.5, 1.25, .8],
                        help="Requested target edge lengths in mm, coarse to fine")
    args = parser.parse_args()
    if len(args.size) < 3 or any(a <= b for a, b in zip(args.size, args.size[1:])) or args.size[-1] <= 0:
        parser.error("Use at least three strictly decreasing positive target sizes")
    output = args.output.resolve()
    output.mkdir(parents=True)
    step, mesher = args.step.resolve(), args.mesher.resolve()
    step_hash, executable_hash = sha256(step), sha256(mesher)
    diagonal = (LENGTH**2+WIDTH**2+HEIGHT**2)**.5
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    rows = []
    for i, size in enumerate(args.size, start=1):
        folder = output / f"level-{i}"
        folder.mkdir()
        archive = folder / "mesh.paacsmesh"
        command = [str(mesher), str(step), str(size), str(archive), "--native-cad", "--elements", "tet4",
                   "--threads", "1", "--native-validation-threads", "1", "--report-json",
                   "--max-work", "300000000", "--max-cells", "1000000"]
        started = time.perf_counter()
        process = subprocess.run(command, env=env, capture_output=True, text=True, errors="replace", timeout=90)
        elapsed = time.perf_counter()-started
        (folder / "mesh.stdout.log").write_text(process.stdout, encoding="utf-8")
        (folder / "mesh.stderr.log").write_text(process.stderr, encoding="utf-8")
        result = read_report(process.stdout)
        if (process.returncode or not result.get("audit_valid") or not result.get("archive_published") or
                result.get("solids") != 1 or result.get("millimetres_per_model_unit") != 1):
            raise ValueError(f"PAACS refused level {i}: {result.get('reason')}")
        receipt = dict(command=command, wall_seconds=elapsed, size_mm=size, size_factor=size/diagonal,
                       step_sha256=step_hash, mesher_sha256=executable_hash, result=result)
        (folder / "mesh-receipt.json").write_text(json.dumps(receipt, indent=2)+"\n", encoding="utf-8")
        solved = solve(archive, folder)
        row = dict(level=i, target_size_mm=size, mesh_wall_seconds=elapsed,
                   native_route=result["bodies"][0]["route"], step_sha256=step_hash,
                   mesher_sha256=executable_hash, **solved)
        rows.append(row)
        (output / "convergence.json").write_text(json.dumps(rows, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        print(f"{size:g} mm · {row['tet4']:,} Tet4 · {row['tip_deflection_mm']:.8f} mm", flush=True)
    assert sha256(step) == step_hash and sha256(mesher) == executable_hash


if __name__ == "__main__":
    main()
