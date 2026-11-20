# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Generate the paired uniform and hole-local plate meshes."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from mesh_visual_data import TET_EDGES, read_archive, sha256

LEVELS = ((4., 2), (2.7, 3), (1.8, 4))


def mesh_metrics(archive, report):
    mesh = read_archive(archive)
    tets = np.array([cell.ids for cell in mesh.cells if cell.kind == 5])
    if len(tets) != report["tet4"]:
        raise ValueError("Archive and fixture receipt Tet4 counts differ")
    edges = np.unique(np.sort(tets[:, TET_EDGES].reshape(-1, 2), axis=1), axis=0)
    a, b = mesh.points[edges].transpose(1,0,2)
    planar_size = np.linalg.norm((b-a)[:, :2], axis=1)
    radii = np.linalg.norm(((a+b)/2)[:, :2], axis=1)
    planar = np.isclose(a[:,2],b[:,2],atol=1e-10,rtol=0)
    near = planar & (radii <= 8.)
    far = planar & (radii >= 18.)
    if not np.any(near) or not np.any(far):
        raise ValueError("Expected both near-hole and far-field sample edges")
    low_rim = (np.isclose(mesh.points[:,2],-1,atol=1e-10,rtol=0) &
               np.isclose(np.linalg.norm(mesh.points[:, :2],axis=1),report["plate_radius_mm"],atol=1e-8,rtol=0))
    rim = low_rim[edges].all(axis=1)
    far_points = mesh.points[np.linalg.norm(mesh.points[:, :2],axis=1)>=18.]
    far_points = far_points[np.lexsort(far_points.T[::-1])].astype("<f8",copy=False)
    return dict(
        archive_sha256=sha256(archive), near_median_inplane_edge_mm=float(np.median(planar_size[near])),
        far_median_inplane_edge_mm=float(np.median(planar_size[far])),
        rim_median_edge_mm=float(np.median(planar_size[rim])),
        edges_near=int(near.sum()), edges_far=int(far.sum()), rim_edges=int(rim.sum()),
        far_points_count=len(far_points), far_points_sha256=hashlib.sha256(far_points.tobytes()).hexdigest(),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=Path, required=True)
    parser.add_argument("--mesher", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True)
    step, mesher = args.step.resolve(), args.mesher.resolve()
    step_hash, mesher_hash = sha256(step), sha256(mesher)
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    results = []
    for i, (far, layers) in enumerate(LEVELS,1):
        for kind, near in (("uniform", far), ("local", .45*far)):
            label = f"{kind}-{i}"
            folder = output/label
            folder.mkdir()
            archive = folder/"mesh.paacsmesh"
            command = [str(mesher),str(step),str(archive),str(far),str(near),str(layers)]
            start = time.perf_counter()
            process = subprocess.run(command, env=environment, capture_output=True, text=True,
                                     errors="replace", timeout=90)
            elapsed = time.perf_counter()-start
            (folder/"mesh.stdout.log").write_text(process.stdout, encoding="utf-8")
            (folder/"mesh.stderr.log").write_text(process.stderr, encoding="utf-8")
            if process.returncode:
                raise ValueError(f"PAACS fixture refused {label}: {process.stderr.strip()}")
            report = json.loads(process.stdout.splitlines()[-1])
            row = dict(label=label, variant=kind, level=i, far_target_mm=far, rim_target_mm=near,
                       command=command, wall_seconds=elapsed, step_sha256=step_hash,
                       mesher_sha256=mesher_hash, archive=str(archive), **report,
                       **mesh_metrics(archive,report))
            if kind == "local" and (row["far_points_sha256"] != results[-1]["far_points_sha256"] or
                    row["far_median_inplane_edge_mm"] != results[-1]["far_median_inplane_edge_mm"]):
                raise ValueError(f"Far-field sampling changed for {label}")
            results.append(row)
            (folder/"mesh-receipt.json").write_text(json.dumps(row,indent=2)+"\n",encoding="utf-8")
            (output/"meshes.json").write_text(json.dumps(results,indent=2)+"\n",encoding="utf-8")
            print(f"{label} {row['tet4']:,} Tet4, rim {row['rim_median_edge_mm']:.3f} mm",flush=True)
    assert sha256(step)==step_hash and sha256(mesher)==mesher_hash


if __name__=="__main__":
    main()
