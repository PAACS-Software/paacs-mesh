#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Verify native CAD request boundaries and existing-archive preservation."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile


def invalid_number(token):
    raise ValueError(f"nonfinite report number: {token}")


def report_from(stdout):
    reports = []
    for line in stdout.splitlines():
        try:
            value = json.loads(line, parse_constant=invalid_number)
        except ValueError:
            continue
        if isinstance(value, dict) and value.get("schema") == "paacs-native-cad-report-v1":
            reports.append(value)
    if len(reports) != 1:
        raise RuntimeError(f"Expected one native report, got {len(reports)}: {stdout[-4000:]}")
    return reports[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--fixtures", type=Path,
                        default=Path(__file__).resolve().parent / "fixtures" / "cad")
    args = parser.parse_args()
    executable = args.executable.resolve()
    rows = []
    with tempfile.TemporaryDirectory(prefix="paacs-native-cad-cli-") as directory:
        output = Path(directory) / "held.paacsmesh"
        output.write_bytes(b"existing archive sentinel")

        def run(name, fixture, options=(), success=False, reason=None, size=".3"):
            previous = output.read_bytes()
            command = [str(executable), str(args.fixtures / fixture), size, str(output),
                       "--native-cad", "--report-json", *options]
            child = subprocess.run(command, capture_output=True, text=True, timeout=60)
            receipt = report_from(child.stdout)
            if (child.returncode == 0) != success:
                raise RuntimeError(f"{name}: code={child.returncode}, report={receipt}, stderr={child.stderr}")
            if success:
                if receipt.get("status") != "succeeded" or receipt.get("reason") != "none":
                    raise RuntimeError(f"{name}: wrong success report: {receipt}")
                if receipt.get("mesh_published") is not True or receipt.get("archive_published") is not True:
                    raise RuntimeError(f"{name}: missing publication: {receipt}")
                if not output.read_bytes().startswith(b"PAACSM8\0"):
                    raise RuntimeError(f"{name}: archive magic missing")
                if receipt.get("nodes", 0) <= 0 or receipt.get("cells", 0) <= 0:
                    raise RuntimeError(f"{name}: empty successful mesh")
                if receipt.get("millimetres_per_model_unit") != 1:
                    raise RuntimeError(f"{name}: STEP transfer units missing")
            else:
                if receipt.get("mesh_published") is not False or receipt.get("archive_published") is not False:
                    raise RuntimeError(f"{name}: refusal published mesh/archive: {receipt}")
                if output.read_bytes() != previous:
                    raise RuntimeError(f"{name}: replaced existing archive")
                normalized = str(receipt.get("reason", "")).replace("_", "-")
                if normalized != reason:
                    raise RuntimeError(f"{name}: expected {reason}, got {normalized}: {receipt}")
            rows.append({"case": name, "passed": True, "returncode": child.returncode,
                         "archive_preserved": not success, "result": receipt})
            return receipt

        cylinder = run("native_cylinder", "cylinder.step", success=True)
        if cylinder.get("curved_faces") != 1 or cylinder["bodies"][0].get("route") != "cylinder":
            raise RuntimeError(f"cylinder did not consume native curved CAD: {cylinder}")
        sized = run("curved_face_and_box", "cylinder.step",
                    ("--face-size", "1", ".38", "--size-box", "-.3", "-.3", ".7", ".3", ".3", "1.3", ".31",
                     "--min-size", ".34", "--growth", "1.55", "--max-work", "100000000"),
                    success=True, size=".7")
        sized_body = sized["bodies"][0]
        if (sized.get("native_spatial_sizing") is not True or sized.get("audit_valid") is not True
                or sized.get("face_size_controls") != 1 or sized.get("size_box_controls") != 1
                or sized_body.get("route") != "native_spatial_sizing"
                or sized_body.get("local_delaunay_refinement_nodes", 0) == 0):
            raise RuntimeError(f"curved controls did not refine and audit the native volume: {sized}")
        standard = run("planar_box", "box.step", success=True, size="3e-1")
        standard_bytes = output.read_bytes()
        advanced = run("advanced_recognized_box", "box.step", ("--recovery-mode", "advanced"), success=True)
        quality_disabled = run("advanced_box_with_regular_quality_disabled", "box.step",
                               ("--recovery-mode", "advanced", "--no-improve-regular-quality"), success=True)
        quality_enabled = run("standard_box_with_regular_quality_enabled", "box.step",
                              ("--improve-regular-quality",), success=True)
        if standard.get("improve_regular_quality") is not False or \
                advanced.get("improve_regular_quality") is not True or \
                quality_disabled.get("improve_regular_quality") is not False or \
                quality_enabled.get("improve_regular_quality") is not True or \
                output.read_bytes() != standard_bytes:
            raise RuntimeError("regular quality policies did not preserve the recognized box route")
        reset = run("explicit_standard_box", "box.step",
                    ("--recovery-mode", "advanced", "--recovery-mode", "standard"), success=True)
        if standard["recovery_mode"] != "standard" or advanced["recovery_mode"] != "advanced" or \
                reset["recovery_mode"] != "standard" or output.read_bytes() != standard_bytes or \
                advanced["bodies"][0]["constraint_recovery"]["diagnostic"] != "not_attempted":
            raise RuntimeError("recovery mode did not preserve the standard or recognized route")
        local_delaunay = run("advanced_box_with_local_delaunay", "box.step",
                             ("--recovery-mode", "advanced", "--live-material-recovery",
                              "--local-delaunay-refinement"), success=True)
        body = local_delaunay["bodies"][0]
        if (standard.get("local_delaunay_refinement") is not False
                or advanced.get("local_delaunay_refinement") is not False
                or local_delaunay.get("local_delaunay_refinement") is not True
                or body.get("local_delaunay_support_nodes") != 0
                or body.get("local_delaunay_refinement_nodes") != 0
                or body.get("local_delaunay_quality_edits") != 0
                or body.get("local_delaunay_used_global_fallback") is not False
                or output.read_bytes() != standard_bytes):
            raise RuntimeError("local Delaunay option did not preserve the recognized box route")
        assembly = run("separated_assembly", "assembly.step", success=True)
        if assembly.get("solids") != 2 or len(assembly.get("bodies", [])) != 2:
            raise RuntimeError(f"assembly bodies were lost: {assembly}")
        sphere = run("native_surface_sphere", "sphere.step", ("--warning-quality", ".25"), success=True)
        if sphere["bodies"][0].get("route") != "native_surface_star" or sphere["bodies"][0].get("native_volume_reason") != "none":
            raise RuntimeError(f"sphere fallback certificate missing: {sphere}")
        refusals = [
            ("disabled_sphere_fallback", "sphere.step", ("--no-native-fallback",), "unsupported-body"),
            ("malformed_brep", "malformed.brep", ("--brep",), "parse-failed"),
            ("unknown_face", "cylinder.step", ("--face-size", "999999", ".1"), "invalid-request"),
            ("unknown_spatial_face", "cylinder.step",
             ("--native-sizing", "--face-size", "999999", ".1"), "invalid-request"),
            ("inverted_size_box", "box.step",
             ("--size-box", ".4", ".2", ".2", ".3", ".8", ".8", ".1"), "invalid-request"),
            ("nonfinite_size_box", "box.step",
             ("--size-box", ".2", ".2", ".2", ".8", ".8", ".8", "nan"), "invalid-request"),
            ("spatial_sizing_assembly", "assembly.step",
             ("--size-box", ".2", ".2", ".2", ".8", ".8", ".8", ".1"), "unsupported-control"),
            ("strict_quality", "cylinder.step", ("--min-quality", "1"), "quality-requirement-not-met"),
            ("zero_work", "cylinder.step", ("--max-work", "0"), "budget-exceeded"),
            ("one_work", "cylinder.step", ("--max-work", "1"), "budget-exceeded"),
            ("working_storage", "cylinder.step", ("--max-working-bytes", "1"), "budget-exceeded"),
            ("audit_storage", "cylinder.step", ("--max-audit-bytes", "1"), "budget-exceeded"),
            ("node_budget", "cylinder.step", ("--max-nodes", "1"), "budget-exceeded"),
            ("cell_budget", "cylinder.step", ("--max-cells", "1"), "budget-exceeded"),
            ("output_budget", "cylinder.step", ("--max-bytes", "1"), "budget-exceeded"),
            ("fallback_cell_budget", "sphere.step", ("--max-cells", "1"), "budget-exceeded"),
            ("fallback_node_budget", "sphere.step", ("--max-nodes", "1"), "budget-exceeded"),
            ("fallback_working_storage", "sphere.step", ("--max-working-bytes", "1"), "budget-exceeded"),
            ("fallback_strict_quality", "sphere.step", ("--min-quality", "1"), "quality-requirement-not-met"),
            ("file_budget", "box.step", ("--max-file-bytes", "1"), "file-budget"),
            ("invalid_budget_number", "box.step", ("--max-nodes", "-1"), "invalid-request"),
            ("invalid_element_flag", "box.step", ("--elements", "tet10"), "invalid-request"),
            ("missing_element_flag", "box.step", ("--elements",), "invalid-request"),
            ("invalid_recovery_mode", "box.step", ("--recovery-mode", "fast"), "invalid-request"),
            ("missing_recovery_mode", "box.step", ("--recovery-mode",), "invalid-request"),
            ("nan_growth", "cylinder.step", ("--growth", "nan"), "invalid-request"),
            ("underflow_box_coordinate", "box.step",
             ("--size-box", "1e-999", ".2", ".2", ".8", ".8", ".8", ".1"), "invalid-request"),
            ("minimum_without_control", "cylinder.step", ("--min-size", ".1"), "unsupported-control"),
        ]
        for name, fixture, options, reason in refusals:
            run(name, fixture, options, reason=reason)
        for name, size in (("nan_size", "nan"), ("trailing_size", ".3junk"),
                           ("leading_space", " .3"), ("trailing_space", ".3 "),
                           ("leading_plus", "+.3"), ("hex_size", "0x1p-2"),
                           ("decimal_comma", "0,3"), ("overflow_size", "1e999"),
                           ("zero_size", "0")):
            run(name, "box.step", reason="invalid-request", size=size)
        # Use a private copy: the input itself must survive output alias errors.
        source = Path(directory) / 'source.step'
        source_bytes = (args.fixtures / 'cylinder.step').read_bytes()
        source.write_bytes(source_bytes)
        alias = Path(directory) / 'source-alias.step'
        os.link(source, alias)
        for name, destination in [('same_input_output', source), ('hardlink_input_output', alias)]:
            child = subprocess.run([str(executable), str(source), '.3', str(destination),
                                    '--native-cad', '--report-json'],
                                   capture_output=True, text=True, timeout=60)
            receipt = report_from(child.stdout)
            if child.returncode != 2 or receipt.get('reason') != 'invalid_request' or \
                    receipt.get('mesh_published') or receipt.get('archive_published') or \
                    source.read_bytes() != source_bytes or destination.read_bytes() != source_bytes:
                raise RuntimeError(f'{name}: source alias was not preserved: {receipt}')
            rows.append({'case': name, 'passed': True, 'returncode': child.returncode,
                         'source_preserved': True, 'result': receipt})
    print(json.dumps({"schema": "paacs-native-cad-cli-checks-v1", "passed": len(rows),
                      "cases": rows}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
