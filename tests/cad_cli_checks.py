#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Exercise publication, strict options and preservation of existing output."""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("executable", type=Path)
p.add_argument("--fixtures", type=Path, default=Path(__file__).resolve().parent / "fixtures" / "cad")
a = p.parse_args()
rows = []
with tempfile.TemporaryDirectory(prefix="paacs-cli-") as directory:
    output = Path(directory) / "held.paacsmesh"
    def run(name, fixture, options, succeeds):
        previous = output.read_bytes() if output.exists() else None
        result = subprocess.run([str(a.executable.resolve()), str(a.fixtures / fixture), ".3",
                                 str(output), *options], capture_output=True, text=True, timeout=120)
        if (result.returncode == 0) != succeeds:
            raise RuntimeError(f"{name}: code={result.returncode}: {result.stderr}")
        row = dict(case=name, passed=True, returncode=result.returncode)
        if succeeds:
            data = [json.loads(s) for s in result.stdout.splitlines() if s.startswith("{")]
            assert len(data) == 1 and output.read_bytes().startswith(b"PAACSM8\0")
            row["result"] = data[0]
        else:
            assert output.read_bytes() == previous, f"{name}: replaced existing output"
            row["diagnostic"] = result.stderr.strip()
        rows.append(row)
        return row.get("result")
    base = run("ordinary_holes", "two_holes.step", [], True)
    for ratio in ["4", ".25"]:
        r = run("graded_" + ratio, "two_holes.step", ["--layers", "12", "--layer-ratio", ratio], True)
        assert r["layers"] == 12 and r["layer_ratio"] == float(ratio)
    r = run("reverse_source", "two_holes.step", ["--layers", "12", "--source-face", str(base["target_face"])], True)
    assert r["source_face"] == base["target_face"]
    run("mapped_grade", "box.brep", ["--layers", "12", "--layer-ratio", "4"], True)
    run("perforated_plate", "perforated_plate.step", [], True)
    cases = [
        ("zero_layers", "box.step", ["--layers", "0"]),
        ("negative_layers", "box.step", ["--layers", "-2"]),
        ("trailing_layers", "box.step", ["--layers", "2junk"]),
        ("nan_ratio", "box.step", ["--layer-ratio", "nan"]),
        ("negative_ratio", "box.step", ["--layer-ratio", "-1"]),
        ("excess_growth", "box.step", ["--layers", "2", "--layer-ratio", "4"]),
        ("unknown_source", "box.step", ["--source-face", "999"]),
        ("tet_layers", "box.step", ["--strategy", "tet", "--layers", "12"]),
        ("schedule_budget", "box.step", ["--layers", "4294967295", "--max-nodes", "100"]),
        ("audit_budget", "box.step", ["--max-audit-bytes", "1"]),
        ("quality_floor", "box.step", ["--min-quality", "1"]),
        ("curved_refusal", "cylinder.step", []),
        ("assembly_refusal", "assembly.step", []),
        ("missing_file", "absent.step", []),
    ]
    for name, fixture, options in cases:
        run(name, fixture, options, False)
print(json.dumps({"passed": len(rows), "cases": rows}, indent=2))
