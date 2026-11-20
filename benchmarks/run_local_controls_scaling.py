#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Fresh-process CPU/time/RSS measurements; run with no competing work."""
import argparse,json,resource,statistics,subprocess,sys
from pathlib import Path

def main():
    if len(sys.argv)>1 and sys.argv[1]=='--one':
        run=subprocess.run(sys.argv[2:],capture_output=True,text=True)
        print(json.dumps(dict(returncode=run.returncode,stdout=run.stdout,stderr=run.stderr,
                              peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)))
        return
    p=argparse.ArgumentParser();p.add_argument('executable',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    cases=[(kind,base,base/4,1.3) for base in [.5,.25,.125] for kind in ['uniform-fine','cube','sphere','face']]
    cases += [(kind,n,.01,1.2) for n in [8,16,32] for kind in ['structured','unstructured','mixed']]
    rows=[]
    for kind,size,local,growth in cases:
        for repeat in range(3):
            for threads in ([1,4] if repeat%2==0 else [4,1]):
                cmd=[str(a.executable.resolve()),kind,str(size),str(local),str(growth),str(threads)]
                receipt=json.loads(subprocess.run([sys.executable,__file__,'--one',*cmd],capture_output=True,text=True,check=True).stdout)
                if receipt['returncode']:
                    raise RuntimeError(receipt)
                row=json.loads(receipt['stdout']);row.update(repeat=repeat,parameter=size,peak_rss_kib=receipt['peak_rss_kib']);rows.append(row)
        print(kind,size,flush=True)
    summary=[]
    for kind,parameter in sorted({(r['case'],r['parameter']) for r in rows}):
        group=[r for r in rows if (r['case'],r['parameter'])==(kind,parameter)]
        assert len({r['hash'] for r in group})==1,'repeat/thread mesh mismatch'
        for threads in [1,4]:
            runs=[r for r in group if r['threads']==threads];result=runs[0].copy();result.pop('repeat')
            for key in list(result):
                if key.endswith('Seconds'):
                    result[key]=statistics.median(r[key] for r in runs)
            result['peak_rss_kib']=max(r['peak_rss_kib'] for r in runs);result['repeats']=len(runs);summary.append(result)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(dict(raw=rows,summary=summary,notes=[
        '126 fresh processes; 3 repeats per configuration/thread count. No competing builds or test/render jobs.',
        'wallSeconds measures the complete mesher call including validation/audit. Field compilation, source triangulation, mesh hashing and optional I/O are excluded.',
        'RSS includes process startup, source generation and hashing. Small processes have a measurement floor from the launcher.',
        'Uniform-fine uses the local target everywhere; adaptive controls refine the requested region and grade into the background. These are different meshes, not equal-error solver comparisons.',
        'Layer growth is physical thickness ratio; adaptive actualGrowth is cbrt(volume) across face-neighbor tetrahedra. These metrics are not interchangeable.',
        'All repeat and one/four-thread mesh hashes agree.'
    ]),indent=2)+'\n')

if __name__=='__main__':main()
