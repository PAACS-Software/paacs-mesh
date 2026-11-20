#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Isolated-process geometry qualification and RSS measurements (no dependencies)."""
import argparse
import csv
import json
import platform
import statistics
import subprocess
import sys
from pathlib import Path

# A fresh wrapper process has exactly one measured child. Using cumulative
# RUSAGE_CHILDREN in the long-lived driver would report earlier cases' peaks.
if len(sys.argv)>1 and sys.argv[1]=="--measure-child":
    result=subprocess.run(sys.argv[2:],text=True,capture_output=True,timeout=180)
    peak=None
    try:
        import resource
        peak=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        if platform.system()=="Darwin":peak/=1024
        peak=int(peak)
    except ImportError:
        pass
    print(json.dumps({"returncode":result.returncode,"stdout":result.stdout,
                      "stderr":result.stderr,"peak_rss_kib":peak}))
    sys.exit(0)

parser = argparse.ArgumentParser()
parser.add_argument("executable", type=Path)
parser.add_argument("--output", type=Path, default=Path("geometry-results"))
parser.add_argument("--repeats", type=int, default=3)
parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
args = parser.parse_args()
if args.repeats < 1 or any(t < 1 for t in args.threads):
    parser.error("repeats and threads must be positive")
exe = args.executable.resolve()
names = subprocess.check_output([str(exe), "--list"], text=True).splitlines()
args.output.mkdir(parents=True, exist_ok=True)
rows = []
for case, factor in [(n, 1.0) for n in names] + [("block_unit", .2), ("block_unit", .08)]:
    for threads in args.threads:
        for repeat in range(args.repeats):
            command = [str(exe), "--case", case, "--single", "--threads", str(threads),
                       "--size-factor", str(factor)]
            wrapper = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                                      "--measure-child", *command], text=True,
                                     capture_output=True, check=True, timeout=185)
            measured=json.loads(wrapper.stdout)
            records = [json.loads(line) for line in measured["stdout"].splitlines()
                       if line.startswith("{")]
            if len(records) != 1:
                raise RuntimeError(f"invalid result for {case}: {measured['stderr']}")
            row = records[0]
            row.update(repeat=repeat, size_factor=factor, returncode=measured["returncode"])
            row["process_peak_rss_kib"] = measured["peak_rss_kib"]
            rows.append(row)
            if not row.get("passed") or measured["returncode"]:
                (args.output / "failed.json").write_text(json.dumps(rows, indent=2) + "\n")
                raise RuntimeError(f"qualification failed: {case}: {measured['stderr']}")
raw = {"platform": platform.platform(), "python": platform.python_version(),
       "measurement": "Each invocation is isolated. RSS includes geometry construction, meshing, independent validation, quality analysis, and archive hashing. Stage timers measure only their named mesher phases. No warmup; medians across fresh processes.",
       "records": rows}
(args.output / "raw.json").write_text(json.dumps(raw, indent=2) + "\n")
summary = []
for case, factor in dict.fromkeys((r["case"], r["size_factor"]) for r in rows):
    hashes = {r["hash"] for r in rows if r["case"] == case and r["size_factor"] == factor}
    if len(hashes) != 1:
        raise RuntimeError(f"nondeterministic mesh: {case} factor={factor}")
    for threads in args.threads:
        selected = [r for r in rows if r["case"] == case and r["size_factor"] == factor
                    and r["threads"] == threads]
        row = {k: selected[0][k] for k in ["case", "strategy", "threads", "nodes", "cells",
                 "volume_cells", "output_bytes", "minimum_quality", "mean_quality", "hash"]}
        row["size_factor"] = factor
        for key in ["analysis_seconds", "generation_seconds", "audit_seconds", "wall_seconds"]:
            row[key] = statistics.median(r[key] for r in selected)
        peaks=[r["process_peak_rss_kib"] for r in selected if r["process_peak_rss_kib"] is not None]
        row["process_peak_rss_kib"] = max(peaks) if peaks else None
        summary.append(row)
(args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
with (args.output / "summary.csv").open("w", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=list(summary[0]))
    writer.writeheader()
    writer.writerows(summary)
print(f"Passed {len(rows)} isolated runs; stable hashes across repeats/thread counts.")
