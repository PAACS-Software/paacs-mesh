#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Isolated processes; identical requests for swept and uniform cylinders."""
import argparse,json,resource,statistics,subprocess,sys
from pathlib import Path

def main():
    if len(sys.argv)>1 and sys.argv[1]=='--one':
        run=subprocess.run(sys.argv[2:],capture_output=True,text=True)
        print(json.dumps(dict(returncode=run.returncode,stdout=run.stdout,stderr=run.stderr,
                              peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)))
        return
    p=argparse.ArgumentParser();p.add_argument('executable',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--repeats',type=int,default=3);a=p.parse_args();rows=[]
    for factor in [1,.5,.25,.125]:
        for family in ['tet4','tet10']:
            for repeat in range(a.repeats):
                for strategy in (['swept','uniform'] if repeat%2==0 else ['uniform','swept']):
                    for threads in ([1,4] if repeat%2==0 else [4,1]):
                        cmd=[str(a.executable.resolve()),'--case','cylinder','--elements',family,
                             '--size-factor',str(factor),'--threads',str(threads),'--single']
                        if strategy=='swept':cmd+=['--cylinder-swept']
                        run=subprocess.run([sys.executable,__file__,'--one',*cmd],capture_output=True,text=True,check=True)
                        receipt=json.loads(run.stdout)
                        if receipt['returncode']:raise RuntimeError(receipt)
                        row=json.loads(receipt['stdout']);row.update(strategy=strategy,size=.3*factor,
                            repeat=repeat,peak_rss_kib=receipt['peak_rss_kib'])
                        row['mesher_seconds']=sum(row[k] for k in ['generation_seconds','validation_seconds','audit_seconds'])
                        rows.append(row)
            print(factor,family,flush=True)
    summary=[]
    for size,family,strategy in sorted({(r['size'],r['element_family'],r['strategy']) for r in rows}):
        group=[r for r in rows if (r['size'],r['element_family'],r['strategy'])==(size,family,strategy)]
        assert len({r['hash'] for r in group})==1,'non-deterministic mesh'
        for threads in [1,4]:
            runs=[r for r in group if r['threads']==threads];row=runs[0].copy();row.pop('repeat')
            for key in ['mesher_seconds','generation_seconds','validation_seconds','audit_seconds','wall_seconds']:
                row[key]=statistics.median(r[key] for r in runs)
            row['peak_rss_kib']=max(r['peak_rss_kib'] for r in runs);row['repeats']=len(runs);summary.append(row)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(dict(raw=rows,summary=summary,notes=[
        'Radius 1, height 2; same requested nominal size and geometry tolerances. Different strategies produce different counts.',
        'GCC16 Release C++26, strict FP. No competing builds/tests during measurements.',
        'Mesher time excludes independent test oracles/hash and file I/O. RSS includes test oracles/hash; fresh subprocess per sample.',
        'Tet10 uses straight-sided geometry. This run measures mesh generation.',
        'All repeat and 1/4-thread hashes agree. Phase medians may not sum to the median total.'
    ]),indent=2)+'\n')
if __name__=='__main__':main()
