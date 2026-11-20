#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Fresh-process scaling measurements; includes independent qualification.
Run with no concurrent build/test workload. Standard library only.
"""
import argparse
import json
import platform
import statistics
import subprocess
import sys
from pathlib import Path

if len(sys.argv)>1 and sys.argv[1]=="--child":
    result=subprocess.run(sys.argv[2:],capture_output=True,text=True,timeout=240)
    peak=None
    try:
        import resource
        peak=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        if platform.system()=="Darwin":peak/=1024
    except ImportError:
        pass
    print(json.dumps(dict(returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,peak_rss_kib=peak)))
    sys.exit(0)

p=argparse.ArgumentParser()
p.add_argument("executable",type=Path)
p.add_argument("--output",type=Path,default=Path("primitive-scaling"))
p.add_argument("--repeats",type=int,default=3)
p.add_argument("--threads",type=int,nargs="+",default=[1,4])
a=p.parse_args()
if a.repeats<1 or any(t<1 for t in a.threads):p.error("positive repeats/threads required")
a.output.mkdir(parents=True,exist_ok=True)
matrix={"cylinder":[.5,.25,.125,.0625],"cone":[.5,.25,.125,.0625],
        "torus":[.8,.4,.2,.15],"tube":[.5,.25,.125,.0625]}
rows=[]
for name,factors in matrix.items():
    for factor in factors:
        for repeat in range(a.repeats):
            for threads in a.threads if repeat%2==0 else list(reversed(a.threads)):
                cmd=[str(a.executable.resolve()),"--case",name,"--size-factor",str(factor),"--threads",str(threads),"--single"]
                measured=json.loads(subprocess.check_output([sys.executable,__file__,"--child",*cmd],text=True,timeout=245))
                if measured["returncode"]:raise RuntimeError(measured["stderr"])
                row=json.loads(measured["stdout"]);row.update(size_factor=factor,repeat=repeat,process_peak_rss_kib=measured["peak_rss_kib"])
                row["mesher_seconds"]=sum(row[k] for k in ["generation_seconds","validation_seconds","audit_seconds"])
                rows.append(row)
        print(f"passed {name}, size factor {factor}",flush=True)
summary=[]
for name,factors in matrix.items():
    for factor in factors:
        selected=[r for r in rows if r["case"]==name and r["size_factor"]==factor]
        if len({r["hash"] for r in selected})!=1:raise RuntimeError("mesh nondeterminism")
        for threads in a.threads:
            group=[r for r in selected if r["threads"]==threads]
            row=dict(group[0]);row.pop("repeat")
            for k in ["generation_seconds","validation_seconds","audit_seconds","mesher_seconds","wall_seconds"]:row[k]=statistics.median(r[k] for r in group)
            peaks=[r["process_peak_rss_kib"] for r in group if r["process_peak_rss_kib"] is not None]
            row["process_peak_rss_kib"]=max(peaks) if peaks else None
            summary.append(row)
(a.output/"raw.json").write_text(json.dumps({"platform":platform.platform(),"measurement":"Fresh processes; phase timers exclude independent test oracles and streaming hash; wall/RSS include them. RSS excludes driver. Alternating thread order, median times and maximum RSS. No GPU.","records":rows},indent=2)+"\n")
(a.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(f"{len(rows)} runs passed with stable mesh hashes")
