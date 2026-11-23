"""Reproduce selected-wall refinement with an installed PAACS Mesh wheel.

Run with the interpreter whose installed package is being measured, preferably
``python -I benchmarks/refinement_diagnostics.py --output build/quality-feedback/NAME``.
Each case is a separate bounded process. Native validity policy is unchanged.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import struct
import subprocess
import sys
import time


TET_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
WIDTHS = (2, 3, 3, 6, 4, 4, 10, 6, 8, 5)


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def stats(values):
    ordered = sorted(values)
    if not ordered:
        return {"count": 0}
    def q(fraction):
        position = (len(ordered) - 1) * fraction
        lower = int(position)
        upper = min(lower + 1, len(ordered) - 1)
        return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])
    mean = statistics.fmean(ordered)
    return {"count": len(ordered), "minimum": ordered[0], "p01": q(.01),
            "p05": q(.05), "median": q(.5), "p95": q(.95), "maximum": ordered[-1],
            "mean": mean, "coefficient_of_variation": statistics.pstdev(ordered) / mean if mean else None}


def measure(path):
    """Independently measure public archive data; native audits remain authoritative."""
    data = Path(path).read_bytes()
    if len(data) < 40 or len(data) > 512 * 1024 * 1024:
        raise ValueError("Unexpected archive byte count")
    magic, version, flags, nodes, cells, connectivity = struct.unpack_from("<8sIIQQQ", data)
    if (magic, version, flags) != (b"PAACSM8\0", 1, 0) or nodes > 5_000_000 or cells > 10_000_000:
        raise ValueError("Unexpected archive header")
    points, node_tags, at = [], [], 40
    for _ in range(nodes):
        x, y, z, dimension, tag = struct.unpack_from("<dddii", data, at)
        if not all(math.isfinite(value) for value in (x, y, z)):
            raise ValueError("Nonfinite coordinate")
        points.append((x, y, z))
        node_tags.append((dimension, tag))
        at += 32
    face_edges, face_cells = {}, {}
    qualities, volumes, sizes, ratios, worst = [], [], [], [], []
    observed_connectivity = 0
    by_radius = {"r_0_to_025": [], "r_025_to_05": [], "r_05_to_075": [], "r_075_to_1": []}
    for cell_id in range(cells):
        kind, dimension, tag = struct.unpack_from("<B3xii", data, at)
        at += 12
        if kind >= len(WIDTHS):
            raise ValueError("Unexpected cell type")
        width = WIDTHS[kind]
        indices = struct.unpack_from(f"<{width}q", data, at)
        at += 8 * width
        observed_connectivity += width
        if not all(0 <= index < nodes for index in indices):
            raise ValueError("Invalid connectivity")
        if kind == 2 and dimension == 2:
            face_cells[tag] = face_cells.get(tag, 0) + 1
            face_edges.setdefault(tag, set()).update(tuple(sorted((indices[i], indices[(i + 1) % 3]))) for i in range(3))
        if kind != 5:
            continue
        xyz = [points[index] for index in indices]
        lengths2 = [sum((xyz[a][axis] - xyz[b][axis]) ** 2 for axis in range(3)) for a, b in TET_EDGES]
        a, b, c = [tuple(xyz[index][axis] - xyz[0][axis] for axis in range(3)) for index in (1, 2, 3)]
        signed_volume = (a[0] * (b[1] * c[2] - b[2] * c[1]) +
                         a[1] * (b[2] * c[0] - b[0] * c[2]) +
                         a[2] * (b[0] * c[1] - b[1] * c[0])) / 6
        volume = abs(signed_volume)
        quality = 12 * (3 * volume) ** (2 / 3) / sum(lengths2) if sum(lengths2) else 0
        center = tuple(sum(point[axis] for point in xyz) / 4 for axis in range(3))
        radius = math.hypot(center[0], center[1])
        key = ("r_0_to_025" if radius < .25 else "r_025_to_05" if radius < .5
               else "r_05_to_075" if radius < .75 else "r_075_to_1")
        by_radius[key].append(quality)
        qualities.append(quality)
        volumes.append(volume)
        sizes.append((6 * math.sqrt(2) * volume) ** (1 / 3))
        ratios.append(math.sqrt(max(lengths2) / min(lengths2)) if min(lengths2) else 1e300)
        worst.append({"cell_id": cell_id, "body_tag": tag, "mean_ratio": quality,
                      "vertices": [{"node_id": index, "point": points[index], "classification": node_tags[index]}
                                   for index in indices],
                      "centroid": center, "centroid_radius": radius,
                      "wall_distance": abs(1 - radius), "cap_distance": min(abs(center[2]), abs(2 - center[2]))})
    if at != len(data) or observed_connectivity != connectivity:
        raise ValueError("Archive count mismatch")
    return {"measurement_policy": "Euclidean chord lengths, not geodesic sizes or strict sizing certificates; Tet4 mean ratio independently recomputed; radius bins describe this owned R1/H2 fixture only",
            "nodes": nodes, "cells": cells, "tet4": len(qualities),
            "node_classifications": dict(sorted(Counter(f"{dimension}:{tag}" for dimension, tag in node_tags).items())),
            "archive_sha256": hashlib.sha256(data).hexdigest(),
            "face_triangles": {str(tag): count for tag, count in sorted(face_cells.items())},
            "face_unique_edge_lengths_mm": {str(tag): stats([math.dist(points[a], points[b]) for a, b in edges]) for tag, edges in sorted(face_edges.items())},
            "face_longest_edges": {str(tag): sorted([
                {"nodes": [a, b], "points": [points[a], points[b]], "length_mm": math.dist(points[a], points[b]),
                 "axial_difference_mm": abs(points[a][2] - points[b][2])}
                for a, b in edges], key=lambda edge: edge["length_mm"], reverse=True)[:20]
                for tag, edges in sorted(face_edges.items())},
            "tet4_mean_ratio": stats(qualities),
            "quality_counts_below": {str(value): sum(quality < value for quality in qualities) for value in (.1, .2, .3, .5)},
            "tet4_volume_mm3": stats(volumes), "tet4_regular_volume_equivalent_edge_mm": stats(sizes),
            "tet4_longest_to_shortest_edge_ratio": stats(ratios),
            "quality_by_centroid_radius": {key: stats(values) for key, values in by_radius.items()},
            "worst_tetrahedra": sorted(worst, key=lambda row: row["mean_ratio"])[:20]}


def options(args):
    return dict(assembly="conforming", workers=1, validation_workers=1,
                min_size=None, growth=1.3, native_sizing=None,
                warning_quality=None, minimum_quality=None, advancing_interior_front=None,
                max_nodes=500_000, max_cells=3_000_000, max_work=100_000_000,
                max_output_bytes=256 * 1024 * 1024,
                max_working_bytes=1024 * 1024 * 1024, max_audit_bytes=512 * 1024 * 1024)


def worker(args):
    import paacs_mesh as pm
    directory = args.output.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    receipt = {"schema": "paacs-refinement-case-v1", "route": args.route,
               "wall_size_mm": args.wall_size, "global_size_mm": .5,
               "package_version": pm.__version__, "native_version": pm.native_version,
               "package_path": str(Path(pm.__file__).resolve()), "interpreter": sys.executable,
               "python_version": sys.version, "host": platform.platform(),
               "source": str(args.source.resolve()), "source_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
               "requested_options": options(args), "native_deadline_seconds": args.max_seconds,
               "retry_policy": {"max_attempts": 1, "allow_advanced": False, "finer_factors": []},
               "status": "incomplete"}
    try:
        imported = time.perf_counter()
        model = pm.read_step(args.source)
        receipt["import_wall_seconds"] = time.perf_counter() - imported
        selection = model.select_faces(surface="cylinder", radius=1, expect=1)
        save(directory / "source-inspection.json", {"info": model.info,
             "faces": model.inspect("faces", topology=True),
             "selection": {"tags": selection.tags, "query": selection.query, "guards": selection.guards}})
        mesh_started = time.perf_counter()
        if args.route == "direct":
            controls = [] if args.wall_size is None else [pm.FaceSize(tag, args.wall_size) for tag in selection.tags]
            mesh = model.mesh(.5, face_sizes=controls, recovery="standard", max_seconds=args.max_seconds, **options(args))
        else:
            recipe = pm.create_recipe(model, .5,
                selected_face_sizes=[] if args.wall_size is None else [(selection, args.wall_size)],
                max_attempts=1, allow_advanced=False, finer_factors=(), max_seconds=args.max_seconds,
                **options(args))
            pm.save_recipe(model, recipe, directory / "recipe.json")
            receipt["recipe_preparation_seconds"] = time.perf_counter() - mesh_started
            result = pm.run_saved_recipe(directory / "recipe.json", args.source)
            save(directory / "workflow-report.json", result.report)
            mesh = result.mesh
        receipt["call_wall_seconds"] = time.perf_counter() - mesh_started
        save(directory / "native-report.json", mesh.report)
        exported = time.perf_counter()
        mesh.write_archive(directory / "mesh.paacsmesh")
        mesh.write_vtk(directory / "mesh.vtk")
        receipt["export_wall_seconds"] = time.perf_counter() - exported
        analyzed = time.perf_counter()
        measurements = measure(directory / "mesh.paacsmesh")
        save(directory / "measurements.json", measurements)
        receipt["measurement_wall_seconds"] = time.perf_counter() - analyzed
        receipt["status"] = "succeeded"
        receipt["native_report"] = "native-report.json"
        receipt["measurements"] = "measurements.json"
    except Exception as error:
        receipt.update(status="failed", exception_type=type(error).__name__, exception=str(error),
                       reason=getattr(error, "reason", None), stage=getattr(error, "stage", None))
        if "mesh_started" in locals():
            receipt["call_wall_seconds"] = time.perf_counter() - mesh_started
        report = getattr(error, "report", None)
        if isinstance(report, dict):
            save(directory / "exception-report.json", report)
            if args.route == "direct":
                save(directory / "native-report.json", report)
            else:
                save(directory / "workflow-report.json", report)
                attempts = report.get("attempts", [])
                if attempts and isinstance(attempts[-1].get("native_report"), dict):
                    save(directory / "native-report.json", attempts[-1]["native_report"])
    finally:
        receipt["case_wall_seconds"] = time.perf_counter() - started
        save(directory / "case.json", receipt)
    return 0


def controller(args):
    directory = args.output.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    save(directory / "benchmark.json", {"schema": "paacs-refinement-benchmark-v1",
        "interpreter": sys.executable, "python_version": sys.version, "host": platform.platform(), "source": str(args.source.resolve()),
        "global_size_mm": .5, "wall_sizes_mm": [None, .5, .3, .2, .15, .1, .08],
        "routes": ["direct", "recipe"], "native_seconds_per_case": args.max_seconds,
        "hard_seconds_per_case": args.hard_seconds, "requested_options": options(args),
        "notes": ["Uniform global-only baseline has no FaceSize; wall .5 is an additional equal-size route control.",
                  "min_size=None means native effective minimum 0; unset face growth inherits 1.3.",
                  "No volume controls, quality relaxation, retry, smoothing override, or source mutation.",
                  "Saved recipe is source-bound, atomically saved then reopened through run_saved_recipe.",
                  "Subprocess end-to-end includes import, setup, meshing, export and measurements; call/native timings separate."]})
    rows = []
    for wall_size in [None, .5, .3, .2, .15, .1, .08]:
        for route in ("direct", "recipe"):
            label = "uniform" if wall_size is None else "wall-" + str(wall_size).replace(".", "p")
            name = label + "-" + route
            case_dir = directory / name
            command = [sys.executable, "-I", str(Path(__file__).resolve()), "--worker", "--source", str(args.source.resolve()),
                       "--output", str(case_dir), "--route", route, "--max-seconds", str(args.max_seconds)]
            if wall_size is not None:
                command += ["--wall-size", str(wall_size)]
            start = time.perf_counter()
            with (directory / (name + ".stdout.log")).open("wb") as stdout, (directory / (name + ".stderr.log")).open("wb") as stderr:
                try:
                    process = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=args.hard_seconds, check=False)
                    process_status = {"exit_code": process.returncode, "timeout": False}
                except subprocess.TimeoutExpired:
                    process_status = {"exit_code": None, "timeout": True}
            elapsed = time.perf_counter() - start
            row = {"case": name, "route": route, "wall_size_mm": wall_size,
                   "end_to_end_seconds": elapsed, **process_status}
            if (case_dir / "case.json").is_file():
                row.update(json.loads((case_dir / "case.json").read_text(encoding="utf-8")))
            else:
                row["status"] = "timeout" if process_status["timeout"] else "process_failed"
            if (case_dir / "native-report.json").is_file():
                report = json.loads((case_dir / "native-report.json").read_text(encoding="utf-8"))
                row.update(nodes=report.get("nodes"), tetrahedra=report.get("tetrahedra"),
                           reason=report.get("reason"), sizing=report.get("sizing"), native_seconds=report.get("seconds"),
                           quality=report.get("quality"), bodies=report.get("bodies"))
            rows.append(row)
            save(directory / "summary.json", rows)
            print(f"{name}: {row['status']} tets={row.get('tetrahedra')} seconds={elapsed:.3f}", flush=True)
    summarize(directory, rows)
    return 0 if all(row["status"] == "succeeded" for row in rows) else 1


def summarize(directory, rows=None):
    """Refresh independent measurements without rerunning or altering native results."""
    rows = rows if rows is not None else json.loads((directory / "summary.json").read_text(encoding="utf-8"))
    compact, pairs = [], []
    for row in rows:
        record = {key: row.get(key) for key in ("case", "status", "wall_size_mm", "nodes", "tetrahedra", "end_to_end_seconds", "call_wall_seconds", "reason")}
        case_dir = directory / row["case"]
        report_path = case_dir / "native-report.json"
        if report_path.is_file():
            report = json.loads(report_path.read_text(encoding="utf-8"))
            quality = report.get("quality", {})
            record.update(quality_minimum=quality.get("minimum"), quality_p05=quality.get("percentile_5"),
                          quality_mean=quality.get("mean"), native_meshing_seconds=report.get("seconds", {}).get("meshing_total"))
        archive_path = case_dir / "mesh.paacsmesh"
        if archive_path.is_file():
            measurements = measure(archive_path)
            save(case_dir / "measurements.json", measurements)
            wall = measurements["face_unique_edge_lengths_mm"].get("1", {})
            record.update(wall_edge_median_mm=wall.get("median"), wall_edge_maximum_mm=wall.get("maximum"),
                          wall_edge_cv=wall.get("coefficient_of_variation"),
                          wall_face_nodes=measurements["node_classifications"].get("2:1", 0),
                          elements_below_quality_01=measurements["quality_counts_below"]["0.1"])
        compact.append(record)
    for row in rows:
        if row["route"] != "direct":
            continue
        direct = directory / row["case"]
        recipe = directory / (row["case"].removesuffix("-direct") + "-recipe") / "mesh.paacsmesh"
        archive = direct / "mesh.paacsmesh"
        pairs.append({"wall_size_mm": row["wall_size_mm"],
                      "both_archives_available": archive.is_file() and recipe.is_file(),
                      "identical_archive_bytes": archive.read_bytes() == recipe.read_bytes() if archive.is_file() and recipe.is_file() else None})
    save(directory / "pair-comparison.json", pairs)
    save(directory / "compact-summary.json", compact)
    if compact:
        with (directory / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
            fields = list(dict.fromkeys(key for record in compact for key in record))
            writer = csv.DictWriter(stream, fields)
            writer.writeheader()
            writer.writerows(compact)
    return compact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1] / "tests/fixtures/cad/cylinder.step")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-seconds", type=float, default=90)
    parser.add_argument("--hard-seconds", type=float, default=150)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--route", choices=("direct", "recipe"), default="direct", help=argparse.SUPPRESS)
    parser.add_argument("--wall-size", type=float, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--summarize-existing", action="store_true", help="Refresh measurements from retained archives; do not rerun meshing")
    args = parser.parse_args()
    if (not args.source.is_file() or not math.isfinite(args.max_seconds) or not math.isfinite(args.hard_seconds)
            or args.max_seconds <= 0 or args.hard_seconds <= args.max_seconds):
        parser.error("source must exist and hard deadline must exceed positive native deadline")
    if args.summarize_existing:
        summarize(args.output.resolve())
        return 0
    return worker(args) if args.worker else controller(args)


if __name__ == "__main__":
    raise SystemExit(main())
