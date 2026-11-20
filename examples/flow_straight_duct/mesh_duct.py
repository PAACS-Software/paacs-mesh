# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Generate three audited PAACS Prism6 wall-graded duct meshes."""

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from mesh_visual_data import read_archive, sha256

LEVELS = ((.60, 3, 1.15, 8), (.30, 4, .90, 12), (.15, 5, .65, 16))
GROWTH = 1.22


def summarize_archive(archive, receipt):
    mesh = read_archive(archive)
    points = mesh.points
    prisms = np.array([cell.ids for cell in mesh.cells if cell.kind == 7], dtype=np.int64)
    volume_tags = {cell.tag for cell in mesh.cells if cell.dim == 3}
    face_counts = Counter(cell.tag for cell in mesh.cells if cell.dim == 2)
    expected = {tag: receipt["face_counts"][side] for side, tag in receipt["face_tags"].items()}
    if len(prisms)!=receipt["prism6"] or volume_tags!={receipt["volume_tag"]} or face_counts!=expected:
        raise ValueError("Prism or original CAD tag counts disagree with archive")
    xyz = points[prisms]
    # A straight Prism6 is the measured in-plane triangle translated along z.
    shift = xyz[:,3:] - xyz[:,:3]
    z_low, z_high = xyz[:,0,2], xyz[:,3,2]
    if (np.any(shift[:,:,:2] != 0) or not np.all(shift[:,:,2] > 0) or
            not np.all(xyz[:,:3,2] == z_low[:,None]) or not np.all(xyz[:,3:,2] == z_high[:,None])):
        raise ValueError("Expected straight positive axial prisms")
    lower = xyz[:,:3]
    normal = np.cross(lower[:,1]-lower[:,0], lower[:,2]-lower[:,0])
    volume = .5*normal[:,2]*shift[:,0,2]
    if not np.all(volume > 0):
        raise ValueError("Nonpositive straight prism volume")
    np.testing.assert_allclose(volume.sum(),4000.,rtol=1e-12)
    stations = [np.unique(points[:,axis]) for axis in range(3)]
    widths = [np.diff(axis) for axis in stations]
    if (not np.array_equal(stations[0],stations[1]) or
            len(widths[0])!=receipt["cross_section_intervals_per_axis"] or
            len(widths[2])!=receipt["axial_layers"]):
        raise ValueError("Cross-section/axial grid differs from planned sweep")
    # Each axial interval must contain exactly the same distinct section triangles.
    _, xy_ids = np.unique(points[:,:2],axis=0,return_inverse=True)
    projected = np.sort(xy_ids[prisms[:,:3]],axis=1)
    reference = np.unique(projected[z_low==stations[2][0]],axis=0)
    for start,end in zip(stations[2],stations[2][1:]):
        on_layer = projected[(z_low==start)&(z_high==end)]
        if len(on_layer)!=len(reference) or not np.array_equal(np.unique(on_layer,axis=0),reference):
            raise ValueError("Prism layer changed section topology")
    first = [float(widths[0][0]),float(widths[0][-1]),float(widths[1][0]),float(widths[1][-1])]
    increments = np.concatenate([widths[0][1:]/widths[0][:-1],widths[1][1:]/widths[1][:-1]])
    maximum_growth = float(max(increments.max(),1/increments.min()))
    return dict(
        archive_sha256=sha256(archive), model_unit="mm", metres_per_unit=.001,
        dimensions_mm=np.ptp(points, axis=0).tolist(), axial_stations_mm=stations[2].tolist(),
        cross_section_stations_mm=stations[0].tolist(), cross_section_cells=2*(len(stations[0])-1)**2,
        first_height_mm=min(first), first_height_by_wall_mm=dict(zip(
            ("x_minus","x_plus","y_minus","y_plus"), first)),
        measured_near_wall_growth=float(widths[0][1]/widths[0][0]),
        maximum_adjacent_growth=maximum_growth,
        minimum_prism_volume_mm3=float(volume.min()), summed_prism_volume_mm3=float(volume.sum()),
    )


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step",type=Path,required=True)
    parser.add_argument("--mesher",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    output=args.output.resolve()
    output.mkdir(parents=True)
    step,mesher=args.step.resolve(),args.mesher.resolve()
    step_hash,executable_hash=sha256(step),sha256(mesher)
    env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    receipts=[]
    for level,(first,rows,core,layers) in enumerate(LEVELS,1):
        label=f"level-{level}"
        folder=output/label
        folder.mkdir()
        archive=folder/"mesh.paacsmesh"
        command=[str(mesher),str(step),str(archive),str(first),str(GROWTH),str(rows),str(core),str(layers)]
        process=subprocess.run(command,env=env,capture_output=True,text=True,errors="replace",timeout=90)
        (folder/"mesh.stdout.log").write_text(process.stdout,encoding="utf-8")
        (folder/"mesh.stderr.log").write_text(process.stderr,encoding="utf-8")
        if process.returncode:
            raise ValueError(f"PAACS refused {label}: {process.stderr.strip()}")
        report=json.loads(process.stdout.splitlines()[-1])
        receipt=dict(label=label,level=level,command=command,archive=str(archive),
                     first_height_target_mm=first,growth_target=GROWTH,core_size_cap_mm=core,
                     step_sha256=step_hash,mesher_sha256=executable_hash,
                     **report,**summarize_archive(archive,report))
        if not receipt["linear_volume_audit_valid"] or not np.allclose(
                list(receipt["first_height_by_wall_mm"].values()),first,atol=1e-12,rtol=0):
            raise ValueError("Archive first wall height or volume audit disagrees with request")
        receipts.append(receipt)
        (folder/"mesh-receipt.json").write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
        (output/"meshes.json").write_text(json.dumps(receipts,indent=2)+"\n",encoding="utf-8")
        print(f"{label}: {receipt['prism6']:,} Prism6, first height {first:g} mm",flush=True)
    assert sha256(step)==step_hash and sha256(mesher)==executable_hash


if __name__=="__main__":
    main()
