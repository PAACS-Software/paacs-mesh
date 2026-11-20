#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Isolated fresh-process spatial sizing benchmark; no external mesher claims."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess

p = argparse.ArgumentParser()
p.add_argument('executable', type=Path)
p.add_argument('output', type=Path)
p.add_argument('--repeats', type=int, default=3)
a = p.parse_args()
if a.repeats < 1:
    p.error('repeats must be positive')
rows, summaries = [], []
cases = ['point', 'segment', 'triangle', 'box', 'ball', 'minimum', 'maximum', 'smooth', 'fine', 'coarse']
for case in cases:
    for threads in [1, 4]:
        run = []
        for repeat in range(a.repeats):
            result = subprocess.run([str(a.executable.resolve()), case, str(threads)],
                                    capture_output=True, text=True, check=True)
            item = json.loads(result.stdout)
            item['repeat'] = repeat
            run.append(item)
            rows.append(item)
        if len({x['hash'] for x in run}) != 1:
            raise RuntimeError('repeat nondeterminism')
        summary = {k: v for k, v in run[0].items() if k not in ['repeat', 'seconds']}
        summary.update(median_seconds=statistics.median(x['seconds'] for x in run),
                       minimum_seconds=min(x['seconds'] for x in run),
                       maximum_seconds=max(x['seconds'] for x in run))
        summaries.append(summary)
        print(case, threads, f"{summary['median_seconds'] * 1000:.3f} ms", flush=True)
for case in cases:
    if len({x['hash'] for x in summaries if x['case'] == case}) != 1:
        raise RuntimeError('thread nondeterminism')
report = dict(executable_sha256=hashlib.sha256(a.executable.read_bytes()).hexdigest(),
              repeats=a.repeats, sample_count=len(rows), summary_count=len(summaries),
              timing_scope='meshAdaptiveBox call including required boundary audit; excludes field compilation, extra size audit, serialization and process startup',
              memory_scope='reported logical output and planner working payloads; not process RSS or allocator peak',
              baseline_scope='fine is constant 0.065: matched minimum only for the point fixture; coarse is constant 0.4',
              all_repeat_and_thread_hashes_equal=True, samples=rows, summaries=summaries)
a.output.parent.mkdir(parents=True, exist_ok=True)
a.output.write_text(json.dumps(report, indent=2) + '\n')
