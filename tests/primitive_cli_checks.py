#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""CLI diagnostics and failed-request preservation of existing output bytes."""
import json
import subprocess
import sys
import tempfile
from pathlib import Path
exe=str(Path(sys.argv[1]).resolve())
with tempfile.TemporaryDirectory() as temp:
    output=Path(temp)/"mesh.paacsmesh"
    count=0
    for name in ["cylinder","cone","frustum","tube","torus","elbow"]:
        run=subprocess.run([exe,name,".3",str(output)],capture_output=True,text=True,check=True)
        assert json.loads(run.stdout)["audit_passed"] and output.read_bytes()[:8]==b"PAACSM8\0"
        count+=1
    for name in ["cylinder","cone","frustum","tube","torus","elbow"]:
        for family in ["tet4","tet10"]:
            run=subprocess.run([exe,name,".3",str(output),"--elements",family],capture_output=True,text=True,check=True)
            r=json.loads(run.stdout)
            assert r[family+"s"]>0 and r["hexes"]==r["pyramids"]==0 and r["audit_passed"]
            if family=="tet10": assert r["quadratic_connectivity_verified"] and r["certified_tet10s"]==r["tet10s"]
            count+=1
    for shape in ["torus", "elbow"]:
        for family in ["tet4", "tet10"]:
            for strategy in ["auto", "uniform", "swept"]:
                run=subprocess.run([exe,shape,".3",str(output),"--elements",family,
                                    "--circular-strategy",strategy],capture_output=True,text=True,check=True)
                r=json.loads(run.stdout)
                assert r["uniform_circular"] == (strategy != "swept")
                assert not r["feature_analysis_enabled"] and r["feature_working_bytes"]==0
                count+=1
    for shape,extra in [("cylinder",[]),("elbow",[]),("torus",[])]:
        base=[exe,shape,".3",str(output),"--elements","tet4",*extra]
        subprocess.run(base,capture_output=True,text=True,check=True)
        original=output.read_bytes()
        run=subprocess.run([*base,"--detect-features","on"],capture_output=True,text=True,check=True)
        r=json.loads(run.stdout)
        assert output.read_bytes()==original and r["feature_analysis_enabled"]
        assert (r["feature_edges"]>0)==(shape!="torus")
        assert r["feature_working_bytes"]>0
        run=subprocess.run([*base,"--detect-features","on","--feature-min-length","10"],
                           capture_output=True,text=True,check=True)
        assert json.loads(run.stdout)["feature_edges"]==0 and output.read_bytes()==original
        count+=2
    run=subprocess.run([exe,"cone",".3",str(output),"--elements","tet4","--warn-quality",".99"],capture_output=True,text=True,check=True)
    r=json.loads(run.stdout)
    assert r["quality_warning"] and r["quality_examples"] and "quality warning" in run.stderr
    count+=1
    for family in ["tet4", "tet10"]:
        for strategy in ["auto", "uniform", "swept"]:
            run=subprocess.run([exe,"cylinder",".3",str(output),"--elements",family,
                                "--cylinder-strategy",strategy],capture_output=True,text=True,check=True)
            r=json.loads(run.stdout)
            assert r["uniform_cylinder"] == (strategy != "swept")
            if strategy != "swept":
                assert r["cylinder_geometry_limited"]
                assert 0 < r["cylinder_axial_spacing"] <= r["cylinder_radial_spacing"] <= .3
            count+=1
    for name,size,extra in [
        ("torus",".3",["--elements","tet4","--circular-strategy","unknown"]),
        ("torus",".3",["--circular-strategy","uniform"]),
        ("torus",".3",["--elements","tet4","--inner-radius",".5","--circular-strategy","uniform"]),
        ("cylinder",".3",["--elements","tet4","--circular-strategy","uniform"]),
        ("cylinder",".3",["--detect-features","maybe"]),
        ("cylinder",".3",["--detect-features","on","--feature-angle","nan"]),
        ("cylinder",".3",["--detect-features","on","--feature-angle-tolerance","-1"]),
        ("cylinder",".3",["--detect-features","on","--feature-min-length","-1"]),
        ("cylinder",".3",["--detect-features","on","--feature-max-edges","1"]),
        ("cylinder",".3",["--elements","tet4","--cylinder-strategy","unknown"]),
        ("cylinder",".3",["--cylinder-strategy","uniform"]),
        ("cone",".3",["--elements","tet4","--cylinder-strategy","uniform"]),
        ("cylinder",".3",["--elements","tet8"]),
        ("cone",".3",["--cone-strategy","unknown"]),
        ("cone",".3",["--quality-policy","ignore"]),
        ("cone",".3",["--elements","tet4","--warn-quality",".99","--quality-policy","reject"]),
        ("cylinder",".3",["--warn-quality","nan"]),
        ("cylinder",".3",["--cone-strategy","graded"]),
        ("unknown",".3",[]),("cone","nan",[]),("torus",".3junk",[]),
        ("torus",".3",["--major-radius","1"]),("tube",".3",["--inner-radius","2"]),
        ("cylinder",".3",["--threads","0"]),("cone",".3",["--max-nodes","3"]),
        ("torus",".3",["--max-working-bytes","1"]),("cone",".3",["--min-quality",".99"]),
        ("cylinder",".3",["--volume-error",".000001"]),("cone",".3",["--radius"]),
        ("elbow",".3",["--angle-degrees","360"]),("cone",".3",["--nonsense","3"])]:
        output.write_bytes(b"existing mesh must survive")
        run=subprocess.run([exe,name,size,str(output),*extra],capture_output=True,text=True)
        assert run.returncode and output.read_bytes()==b"existing mesh must survive",(name,size,extra,run)
        count+=1
print(f"{count} primitive CLI checks passed")
