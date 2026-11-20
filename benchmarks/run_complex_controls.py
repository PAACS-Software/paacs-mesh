#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Fresh-process complex-geometry sizing benchmarks; mesh-call timing only."""
import argparse,json,subprocess,statistics,resource,hashlib
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('executable',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
rows=[];summaries=[]
for case in ['i_beam','multi_bore','star_duct','cooling_channels','notched_frame','thin_bracket','mapped_block']:
 for mode in ['start','fine']:
  for family in ['native','tet4']:
   for threads in [1,4]:
    run=[]
    for repeat in range(3):
     # wait4 measures this child rather than cumulative earlier-child RSS.
     import os,tempfile
     with tempfile.TemporaryFile() as out,tempfile.TemporaryFile() as err:
      child=subprocess.Popen([str(a.executable.resolve()),case,mode,family,str(threads)],stdout=out,stderr=err)
      _,status,usage=os.wait4(child.pid,0);child.returncode=os.waitstatus_to_exitcode(status)
      out.seek(0);err.seek(0);text=out.read().decode();error=err.read().decode()
     if child.returncode: raise RuntimeError((case,mode,family,threads,child.returncode,error))
     item=json.loads(text);item.update(repeat=repeat,peak_rss_kib=usage.ru_maxrss);run.append(item);rows.append(item)
    if len({x['hash'] for x in run})!=1:raise RuntimeError('repeat nondeterminism')
    summary={k:v for k,v in run[0].items() if k not in ['repeat','seconds','peak_rss_kib']}
    summary.update(median_seconds=statistics.median(x['seconds'] for x in run),minimum_seconds=min(x['seconds'] for x in run),maximum_seconds=max(x['seconds'] for x in run),peak_rss_kib=max(x['peak_rss_kib'] for x in run))
    summaries.append(summary)
    print(case,mode,family,threads,f"{summary['median_seconds']*1000:.3f} ms",flush=True)
for r in summaries:
 alt=next(x for x in summaries if all(x[k]==r[k] for k in ['case','mode','family']) and x['threads']!=r['threads'])
 if alt['hash']!=r['hash']:raise RuntimeError('thread nondeterminism')
a.output.write_text(json.dumps(dict(executable_sha256=hashlib.sha256(a.executable.read_bytes()).hexdigest(),repeats=3,sample_count=len(rows),summary_count=len(summaries),samples=rows,summaries=summaries,all_repeat_and_thread_hashes_equal=True),indent=2)+'\n')
