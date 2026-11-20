#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Fresh-process Tet4/Tet10 measurements. Run with no competing build/tests."""
import argparse
import json
import resource
import statistics
import subprocess
import sys
from pathlib import Path

def main():
    if len(sys.argv)>1 and sys.argv[1]=='--one':
        run=subprocess.run(sys.argv[2:],capture_output=True,text=True)
        print(json.dumps(dict(returncode=run.returncode,stdout=run.stdout,stderr=run.stderr,
                              peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)))
        return
    p=argparse.ArgumentParser()
    p.add_argument('executable',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--repeats',type=int,default=3)
    a=p.parse_args();rows=[]
    configurations=[('cylinder','tet4',False),('cylinder','tet10',False),
                    ('cone','tet4',False),('cone','tet10',False),('cone','tet4',True)]
    for shape,family,swept in configurations:
        factors=[1,.5,.25,.125] if shape=='cylinder' else [1,.5,.25,1/6]
        for factor in factors:
            for repeat in range(a.repeats):
                for threads in ([1,4] if repeat%2==0 else [4,1]):
                    command=[str(a.executable.resolve()),'--case',shape,'--elements',family,
                             '--size-factor',str(factor),'--threads',str(threads),'--single']
                    if swept:command+=['--swept']
                    wrapper=subprocess.run([sys.executable,__file__,'--one',*command],capture_output=True,text=True,check=True)
                    receipt=json.loads(wrapper.stdout)
                    if receipt['returncode']:raise RuntimeError(receipt)
                    row=json.loads(receipt['stdout']);row.update(strategy='swept' if swept else 'automatic',
                         size_factor=factor,repeat=repeat,peak_rss_kib=receipt['peak_rss_kib'])
                    row['mesher_seconds']=sum(row[k] for k in ['generation_seconds','validation_seconds','audit_seconds'])
                    rows.append(row)
            print(shape,family,'swept' if swept else 'automatic',factor,flush=True)
    summary=[]
    keys=sorted({(r['case'],r['element_family'],r['strategy'],r['size_factor']) for r in rows})
    for case,family,strategy,factor in keys:
        allruns=[r for r in rows if (r['case'],r['element_family'],r['strategy'],r['size_factor'])==(case,family,strategy,factor)]
        assert len({r['hash'] for r in allruns})==1,'thread/repeat nondeterminism'
        for threads in [1,4]:
            run=[r for r in allruns if r['threads']==threads];row=run[0].copy()
            for key in ['mesher_seconds','generation_seconds','validation_seconds','audit_seconds','wall_seconds']:
                row[key]=statistics.median(r[key] for r in run)
            row['peak_rss_kib']=max(r['peak_rss_kib'] for r in run)
            row.pop('repeat');row['repeats']=len(run);summary.append(row)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(dict(raw=rows,summary=summary,notes=[
        'CPU mesher only; no GPU route. GCC16 Release, C++26.',
        'Times exclude archive serialization. Whole-process peak RSS includes independent geometric checks and hashing.',
        'Reported phase medians need not sum to median total. Tet10 validation includes elevation, edge sorting, certification and volume integration.',
        'Graded and swept cones use the same requested size and tolerances, with different cell counts; this is not an equal-cell-count comparison.',
        'Default resource ceilings are retained. All repeated and cross-thread mesh hashes match.'
    ]),indent=2)+'\n')
if __name__=='__main__':main()
