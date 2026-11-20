#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Compare the frozen rc.2 and production audits in isolated processes."""
import argparse
import json
import platform
import statistics
import subprocess
import sys
from pathlib import Path

if len(sys.argv) > 1 and sys.argv[1] == "--child":
    result = subprocess.run(sys.argv[2:], text=True, capture_output=True, timeout=120)
    peak = None
    try:
        import resource
        peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        if platform.system() == "Darwin":
            peak /= 1024
    except ImportError:
        pass
    print(json.dumps(dict(code=result.returncode, output=result.stdout, error=result.stderr,
                         peak_rss_kib=peak)))
    sys.exit(0)

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("executable", type=Path)
p.add_argument("--output", type=Path, default=Path("audit-results"))
p.add_argument("--repeats", type=int, default=5)
a = p.parse_args()
if a.repeats < 1:
    p.error("repeats must be positive")
a.output.mkdir(parents=True, exist_ok=True)
rows = []
for n in [25, 63, 100]:
    for shuffled in [False, True]:
        for repeat in range(a.repeats):
            # Alternate order to reduce systematic order/host-load bias.
            methods = ["reference", "bucket"] if repeat % 2 == 0 else ["bucket", "reference"]
            for method in methods:
                command = [str(a.executable.resolve()), str(n), method]
                if shuffled:
                    command += ["shuffle"]
                raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()),
                                               "--child", *command], text=True, timeout=125)
                result = json.loads(raw)
                if result["code"]:
                    raise RuntimeError(result)
                row = json.loads(result["output"])
                row.update(repeat=repeat, process_peak_rss_kib=result["peak_rss_kib"])
                rows.append(row)
summary = []
for n in [25, 63, 100]:
    for shuffled in [False, True]:
        selected = [r for r in rows if r["n"] == n and r["shuffled"] == shuffled]
        for key in ["cells", "face_uses", "interior_faces", "exterior_faces"]:
            if len({r[key] for r in selected}) != 1:
                raise RuntimeError("reference diagnostic mismatch")
        info = dict(n=n, shuffled=shuffled, cells=n**3)
        for method in ["reference", "bucket"]:
            samples = [r for r in selected if r["method"] == method]
            info[method + "_seconds"] = statistics.median(r["seconds"] for r in samples)
            info[method + "_temporary_bytes"] = samples[0]["temporary_bytes"]
            peaks = [r["process_peak_rss_kib"] for r in samples if r["process_peak_rss_kib"] is not None]
            info[method + "_peak_rss_kib"] = max(peaks) if peaks else None
        info["speedup"] = info["reference_seconds"] / info["bucket_seconds"]
        info["temporary_reduction_fraction"] = 1-info["bucket_temporary_bytes"]/info["reference_temporary_bytes"]
        summary.append(info)
(a.output / "raw.json").write_text(json.dumps({"platform": platform.platform(), "records": rows}, indent=2)+"\n")
(a.output / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
print(f"{len(rows)} isolated comparisons passed; medians and whole-process peak RSS recorded.")
