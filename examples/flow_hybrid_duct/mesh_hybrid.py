# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Generate three reproducible PAACS hybrid meshes of the same straight STEP duct."""

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

LEVELS = ((.60, 3, 6, 8), (.30, 4, 9, 12), (.15, 5, 12, 16))
GROWTH = 1.22


def inspect(archive, receipt):
    mesh = read_archive(archive)
    p = mesh.points
    np.testing.assert_allclose(p.min(axis=0),[-5,-5,0],rtol=0,atol=1e-12)
    np.testing.assert_allclose(p.max(axis=0),[5,5,40],rtol=0,atol=1e-12)
    volume = [cell for cell in mesh.cells if cell.dim == 3]
    boundary = [cell for cell in mesh.cells if cell.dim == 2]
    count = Counter(cell.kind for cell in volume)
    tags = receipt["face_tags"]
    if count != {5: receipt["tet4"], 7: receipt["prism6"]} or {cell.tag for cell in volume}!={receipt["volume_tag"]}:
        raise ValueError("Archive contains unexpected volume kinds or boundary identity")
    if {cell.tag for cell in boundary} != set(tags.values()):
        raise ValueError("Archive lost an original STEP face tag")

    def cap(name):
        volume_cells = [cell for cell in boundary if cell.tag==tags[name]]
        points = np.unique([node for cell in volume_cells for node in cell.ids])
        signature = Counter((cell.kind,tuple(sorted(p[i,:2].astype("<f8").tobytes() for i in cell.ids)))
                            for cell in volume_cells)
        return len(points), signature

    inlet, outlet = cap("inlet"), cap("outlet")
    if inlet!=outlet or not inlet[1] or set(kind for kind,_ in inlet[1])!={2,4}:
        raise ValueError("Hybrid pressure-periodic cap nodes or face patterns differ")
    prisms = np.array([cell.ids for cell in volume if cell.kind==7])
    corners = p[prisms]
    if not np.array_equal(corners[:,:3,2],corners[:,3:,2]):
        raise ValueError("A wall prism changed its axial station")
    wall_depth = receipt["wall_height_mm"]
    inner = 5-wall_depth
    r = np.max(np.abs(corners[:,:,:2]),axis=2)
    lower, upper = r[:,:3].max(axis=1), r[:,3:].max(axis=1)
    np.testing.assert_allclose(r[:,:3],np.broadcast_to(lower[:,None],r[:,:3].shape),rtol=0,atol=1e-10)
    np.testing.assert_allclose(r[:,3:],np.broadcast_to(upper[:,None],r[:,3:].shape),rtol=0,atol=1e-10)
    offsets = np.unique(np.round(np.r_[5-lower,5-upper],12))
    heights = np.diff(offsets)
    if (len(heights)!=receipt["prism_wall_layers"] or np.any(upper>=lower) or
            not np.allclose([offsets[0],offsets[-1]],[0,wall_depth],atol=1e-10,rtol=0)):
        raise ValueError("Wall prism radial extents do not match the schedule")
    return dict(
        archive_sha256=sha256(archive), model_unit="mm", metres_per_unit=.001,
        dimensions_mm=np.ptp(p,axis=0).tolist(), periodic_cap_points=inlet[0],
        periodic_cap_faces=sum(inlet[1].values()), cap_tri3_count=sum(n for (kind,_),n in inlet[1].items() if kind==2),
        cap_quad4_count=sum(n for (kind,_),n in inlet[1].items() if kind==4),
        wall_first_height_mm=float(heights[0]), wall_offsets_mm=offsets.tolist(), wall_layer_heights_mm=heights.tolist(),
        maximum_adjacent_wall_growth=float(np.max(heights[1:]/heights[:-1])),
        core_width_mm=2*inner, axial_stations_mm=np.unique(p[:,2]).tolist(),
    )


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step",type=Path,required=True)
    parser.add_argument("--mesher",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    output=args.output.resolve()
    output.mkdir(parents=True)
    step,exe=args.step.resolve(),args.mesher.resolve()
    step_hash,exe_hash=sha256(step),sha256(exe)
    env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    rows=[]
    for level,(first,wall,core,axial) in enumerate(LEVELS,1):
        label=f"level-{level}"
        folder=output/label
        folder.mkdir()
        archive=folder/"mesh.paacsmesh"
        command=[str(exe),str(step),str(archive),str(first),str(GROWTH),str(wall),str(core),str(axial)]
        proc=subprocess.run(command,env=env,capture_output=True,text=True,errors="replace",timeout=90)
        (folder/"mesh.stdout.log").write_text(proc.stdout,encoding="utf-8")
        (folder/"mesh.stderr.log").write_text(proc.stderr,encoding="utf-8")
        if proc.returncode:
            raise ValueError(f"PAACS refused {label}: {proc.stderr.strip()}")
        report=json.loads(proc.stdout.splitlines()[-1])
        row=dict(label=label,level=level,archive=str(archive),command=command,
                 step_sha256=step_hash,mesher_sha256=exe_hash,first_height_target_mm=first,growth_target=GROWTH,
                 **report,**inspect(archive,report))
        if not row["linear_volume_audit_valid"] or abs(row["wall_first_height_mm"]-first)>1e-10:
            raise ValueError("Hybrid wall height or volume audit disagrees with request")
        rows.append(row)
        (folder/"mesh-receipt.json").write_text(json.dumps(row,indent=2)+"\n",encoding="utf-8")
        (output/"meshes.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")
        print(f"{label}: {row['prism6']:,} Prism6 + {row['tet4']:,} Tet4, first height {first:g} mm",flush=True)
    assert sha256(step)==step_hash and sha256(exe)==exe_hash


if __name__=="__main__":
    main()
