#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Isolated fresh-process timing/RSS for profile strategies and slab reuse."""
import argparse
import json
import resource
import statistics
import subprocess
import sys
from pathlib import Path

def main():
    if len(sys.argv)>1 and sys.argv[1]=='--one':
        result=subprocess.run(sys.argv[2:],capture_output=True,text=True)
        print(json.dumps(dict(returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,
                              peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)))
        return
    parser=argparse.ArgumentParser();parser.add_argument('executable',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args()
    configurations=[('start',family,strategy,local,1) for family in ['tet4','hex']
                    for strategy in ['sweep','volume'] for local in [.18,.12,.09]]
    configurations += [('fine',family,strategy,.12,1) for family in ['tet4','hex'] for strategy in ['sweep','volume']]
    configurations += [(mode,'tet10',strategy,.12,1) for mode in ['start','fine'] for strategy in ['sweep','volume']]
    configurations += [(mode,family,'volume',.12,0) for mode in ['start','fine'] for family in ['tet4','tet10','hex']]
    rows=[]
    for mode,family,strategy,local,reuse in configurations:
        for repeat in range(3):
            for threads in ([1,4] if repeat%2==0 else [4,1]):
                command=[str(args.executable.resolve()),mode,family,strategy,str(local),str(threads),str(reuse)]
                result=subprocess.run([sys.executable,__file__,'--one',*command],capture_output=True,text=True,check=True)
                receipt=json.loads(result.stdout)
                if receipt['returncode']:raise RuntimeError(receipt)
                row=json.loads(receipt['stdout']);row.update(repeat=repeat,peak_rss_kib=receipt['peak_rss_kib']);rows.append(row)
        print(mode,family,strategy,local,reuse,flush=True)
    summary=[]
    for configuration in configurations:
        mode,family,strategy,local,reuse=configuration
        group=[r for r in rows if (r['mode'],r['family'],r['strategy'],r['local'],r['reuse'])==configuration]
        assert len({r['hash'] for r in group})==1,'repeat/thread mesh mismatch'
        for threads in [1,4]:
            runs=[r for r in group if r['threads']==threads];row=runs[0].copy();row.pop('repeat')
            row['seconds']=statistics.median(r['seconds'] for r in runs)
            row['minimum_seconds']=min(r['seconds'] for r in runs);row['maximum_seconds']=max(r['seconds'] for r in runs)
            row['peak_rss_kib']=max(r['peak_rss_kib'] for r in runs);row['repeats']=len(runs);summary.append(row)
    comparisons=[]
    for mode in ['start','fine']:
        for family in ['tet4','tet10','hex']:
            for threads in [1,4]:
                group=[r for r in summary if r['mode']==mode and r['family']==family and r['threads']==threads and r['strategy']=='volume' and r['local']==.12]
                cached=next(r for r in group if r['reuse']);uncached=next(r for r in group if not r['reuse'])
                comparisons.append(dict(mode=mode,family=family,threads=threads,
                                        speedup=uncached['seconds']/cached['seconds'],identical_hash=cached['hash']==uncached['hash']))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(raw=rows,summary=summary,reuse_comparisons=comparisons,notes=[
        f'{len(rows)} fresh processes; 3 repeats/configuration/thread count. Run only after all builds, tests and renders finish.',
        'seconds measures the full mesh call, including planning, geometry checks and boundary audit; compilation, hashing and optional archive I/O are outside this interval.',
        'Linux peak RSS includes process startup and hashing; small jobs have a launcher measurement floor. This is not the conservative working-byte admission estimate.',
        'local refinement and uniform-fine baselines are different meshes; no equal-error solver or competitor speed claim.',
        'Reused slabs are affine images of validated connectivity with identical two-dimensional sections. All final world-space geometry and boundary checks still run.',
        'Volume sizing keeps a common outer polygon. Growth controls axial slab spacing, not all 3D edges or adjacent-volume sizes.',
        'Parallel generation/validation is deterministic at tested 1/4 threads; planning and Delaunay insertion remain serial. No GPU or many-core scaling claim.'
    ]),indent=2)+'\n')
if __name__=='__main__':main()
