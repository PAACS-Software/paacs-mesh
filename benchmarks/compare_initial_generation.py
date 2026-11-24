"""Compare saved owned-fixture runs and render actual whole-Tet4 cutaways.

This script imports no native engine and performs no CAD queries. Numeric
tables retain hash-paired PAACS Python inspector results. Figure cells are
remeasured from their actual archives. An initial point cloud is never inferred
from a final archive. Output must be a new directory.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import struct
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
ARITY = (2, 3, 3, 6, 4, 4, 10, 6, 8, 5)
FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3))
OWNED = {"tests/fixtures/cad/thin_plate.step", "examples/elasticity_plate_hole/plate-hole.step",
         "tests/fixtures/cad/cylinder.step", "tests/fixtures/cad/sphere.step",
         "tests/fixtures/curved-reference/torus_rotated.step", "tests/fixtures/cad/through_hole.step",
         "tests/fixtures/cad/box.step", "tests/fixtures/cad/oblique_l.step"}
METRICS = ("mean_ratio", "scaled_jacobian", "minimum_dihedral_deg")
BG, INK, SUB = "#f4f7fa", "#17364b", "#526b7c"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_bytes(path, maximum):
    path = Path(path)
    if not path.is_file() or path.stat().st_size > maximum:
        raise ValueError(f"Missing or oversized input: {path.name}")
    value = path.read_bytes()
    if len(value) > maximum:
        raise ValueError("Input grew beyond its byte allowance")
    return value


def read_json(path, maximum=16 * 1024**2):
    raw = read_bytes(path, maximum)
    def reject(value):
        raise ValueError(f"Nonfinite JSON constant: {value}")
    return json.loads(raw.decode("utf-8-sig"), parse_constant=reject), digest(raw)


def relative(path):
    path = Path(path).resolve()
    return path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.name


def artifact(path, sha256):
    return {"path": relative(path), "sha256": sha256}


def check_deadline(started, args):
    if time.perf_counter() - started > args.max_seconds:
        raise ValueError("Comparison/figure deadline exceeded")


def signatures(path, args, started, *, include_tets=False):
    data = read_bytes(path, args.max_bytes)
    if len(data) < 40:
        raise ValueError("Truncated archive")
    magic, version, reserved, nodes, cells, connectivity = struct.unpack_from("<8sIIQQQ", data)
    if (magic, version, reserved) != (b"PAACSM8\0", 1, 0) or nodes > args.max_nodes or cells > args.max_cells:
        raise ValueError("Unsupported archive or count allowance exceeded")
    if 40 + 32 * nodes + 12 * cells + 8 * connectivity != len(data):
        raise ValueError("Archive length/count mismatch")
    modeled = 2 * len(data) + 512 * nodes + 512 * cells + 64 * 1024**2
    if modeled > args.max_working_bytes:
        raise ValueError("Archive analysis modeled-memory allowance exceeded")
    points, tags, encoded_nodes = [], [], []
    for i in range(nodes):
        if i % 4096 == 0:
            check_deadline(started, args)
        raw = data[40 + 32 * i:72 + 32 * i]
        x, y, z, dim, tag = struct.unpack("<3dii", raw)
        if not all(math.isfinite(v) for v in (x, y, z)) or not 0 <= dim <= 3 or tag < 0:
            raise ValueError("Invalid node coordinate or classification")
        points.append((x, y, z)); tags.append((dim, tag)); encoded_nodes.append(raw)
    at, observed_connectivity, tet_count = 40 + 32 * nodes, 0, 0
    boundary_nodes, material_nodes, body_counts, face_counts = set(), set(), Counter(), Counter()
    triangles, tets, tet_ids, oriented, unoriented = [], [], [], [], []
    ordered_triangles = hashlib.sha256()
    for cell_id in range(cells):
        if cell_id % 4096 == 0:
            check_deadline(started, args)
        kind, padding, dim, tag = struct.unpack_from("<B3sii", data, at); at += 12
        if kind >= len(ARITY) or padding != b"\0\0\0" or not 0 <= dim <= 3 or tag < 0:
            raise ValueError("Invalid cell metadata")
        ids = struct.unpack_from("<" + "q" * ARITY[kind], data, at); at += 8 * ARITY[kind]
        if any(i < 0 or i >= nodes for i in ids):
            raise ValueError("Invalid dense node index")
        observed_connectivity += len(ids)
        if kind == 2 and dim == 2:
            triangles.append((ids, tag)); boundary_nodes.update(ids); face_counts[tag] += 1
            vertices = [encoded_nodes[i][:24] for i in ids]
            prefix = struct.pack("<i", tag)
            oriented.append(prefix + min(b"".join(vertices[i:] + vertices[:i]) for i in range(3)))
            unoriented.append(prefix + b"".join(sorted(vertices)))
            ordered_triangles.update(struct.pack("<3qii", *ids, dim, tag))
        elif kind == 5 and dim == 3:
            tet_count += 1; body_counts[tag] += 1; material_nodes.update(ids)
            if include_tets:
                if tet_count > args.max_plot_tets:
                    raise ValueError("Selected figure exceeds Tet4 rendering allowance")
                tets.append(ids); tet_ids.append(cell_id)
        elif dim == 3:
            raise ValueError("Comparison supports Tet4 volume archives")
    if at != len(data) or observed_connectivity != connectivity:
        raise ValueError("Archive connectivity count mismatch")
    def cloud(indices):
        rows = [encoded_nodes[i] for i in indices]
        xyz = [row[:24] for row in rows]
        return {"nodes": len(rows), "ordered_xyz_sha256": digest(b"".join(xyz)),
                "multiset_xyz_sha256": digest(b"".join(sorted(xyz))),
                "ordered_xyz_and_classification_sha256": digest(b"".join(rows)),
                "multiset_xyz_and_classification_sha256": digest(b"".join(sorted(rows)))}
    inside = sorted(i for i in material_nodes if tags[i][0] == 3)
    receipt = {"archive": artifact(path, digest(data)), "nodes": nodes, "tetrahedra": tet_count,
        "body_tet4_counts": dict(sorted(body_counts.items())), "boundary_triangles": len(triangles),
        "face_tag_counts": dict(sorted(face_counts.items())),
        "final_all_points": cloud(range(nodes)), "final_boundary_nodes": cloud(sorted(boundary_nodes)),
        "final_classified_interior_material_nodes": cloud(inside),
        "surface": {"oriented_coordinates_and_face_tags_sha256": digest(b"".join(sorted(oriented))),
                    "unoriented_coordinates_and_face_tags_sha256": digest(b"".join(sorted(unoriented))),
                    "ordered_indexed_triangles_and_tags_sha256": ordered_triangles.hexdigest()},
        "initial_point_cloud": {"available": False,
            "reason": "Final mesh archive only; no stage point-coordinate snapshot supplied"},
        "signature_convention": "Exact encoded binary64 XYZ and native node classifications; triangle rotations preserve orientation; multiset hashes retain duplicates",
        "modeled_analysis_bytes": modeled,
        "scope": "Final archive identity/classifications; no CAD or intersection audit rerun"}
    return receipt, {"points": points, "tags": tags, "triangles": triangles, "tets": tets, "tet_ids": tet_ids}


def group_record(directory, case, cache):
    path = directory.parent / "results.json"
    if not path.is_file():
        return None, None
    if path not in cache:
        value, sha = read_json(path, 64 * 1024**2)
        if not isinstance(value, list):
            return None, None
        cache[path] = value, sha
    rows, sha = cache[path]
    matches = [r for r in rows if r.get("case_id") == case.get("case_id")
               and r.get("repetition", 1) == case.get("repetition", 1)
               and r.get("mesh_sha256") == case.get("mesh_sha256")]
    if len(matches) != 1:
        raise ValueError("Missing or ambiguous hash-paired process receipt")
    row = matches[0]
    for key in ("status", "source_sha256", "call_seconds"):
        if row.get(key) != case.get(key):
            raise ValueError("Process and case receipt disagree")
    return row, artifact(path, sha)


def load_case(directory, args, started, cache):
    directory = Path(directory).resolve()
    case, case_sha = read_json(directory / "case.json")
    spec = case.get("specification", {})
    fixture = spec.get("fixture")
    if fixture not in OWNED:
        raise ValueError("Only explicitly owned fixture paths are shareable in this comparison")
    actual_source = digest(read_bytes(ROOT / fixture, 32 * 1024**2))
    if case.get("source_sha256") != actual_source:
        raise ValueError("Owned source differs from the saved case hash")
    report_name = case.get("native_report_file")
    if report_name not in ("native-report.json", "refusal-report.json"):
        raise ValueError("Case has no recognized Native report pairing")
    report, report_sha = read_json(directory / report_name)
    if case.get("native_report_sha256") != report_sha or report.get("source_sha256") != actual_source:
        raise ValueError("Native report hash/source differs from case receipt")
    process, process_artifact = group_record(directory, case, cache)
    identity = case.get("development_identity", {})
    success = case.get("status") == "succeeded"
    attempts = case.get("external_attempts")
    if attempts is None and isinstance(case.get("retry_count"), int):
        attempts = case["retry_count"] + 1
    result = {"case_id": case.get("case_id"), "status": case.get("status"),
        "case": artifact(directory / "case.json", case_sha), "native_report": artifact(directory / report_name, report_sha),
        "process_receipt": process_artifact, "source": {"fixture": fixture, "sha256": actual_source},
        "binary": {"native_extension_sha256": identity.get("native_extension_sha256", case.get("native_extension_sha256")),
                   "source_snapshot_sha256": identity.get("source_snapshot_sha256"),
                   "source_commit": identity.get("source_commit"),
                   "source_worktree_clean": identity.get("source_worktree_clean"),
                   "scope": "Selected development binary and separate source snapshot; not installed-wheel or release qualification"},
        "requested": {"size_mm": case.get("requested_size_mm", spec.get("size")),
                      "options": case.get("options", case.get("limits")), "sizing": report.get("sizing"),
                      "fidelity": report.get("fidelity")},
        "timings": {"api_call_seconds": case.get("call_seconds"), "native_seconds": report.get("seconds"),
                    "worker_seconds": case.get("worker_seconds"),
                    "process_seconds": process.get("process_seconds") if process else None,
                    "external_attempts": attempts,
                    "external_retry_count": attempts - 1 if isinstance(attempts, int) else None,
                    "scope": "Single saved attempt; Native time includes internal retries. Nested retry receipts must not be added again."},
        "native": {key: report.get(key) for key in ("status", "reason", "resources", "workers", "validation_workers",
                   "audit", "quality", "mesh_scope", "interbody_nonoverlap_certified", "interferences")},
        "bodies": [{key: b.get(key) for key in ("solid_tag", "route", "reason", "cad_volume", "mesh_volume",
                    "relative_volume_error", "maximum_boundary_node_error", "maximum_boundary_sample_error",
                    "surface_chord_error_limit", "fidelity", "route_resolution", "source_shell", "initial_front",
                    "initial_topology_quality", "constraint_recovery", "volume_sizing", "surface_spacing")}
                   for b in report.get("bodies", [])]}
    if not success:
        result["archive"] = None
        return result
    if (report.get("status") != "succeeded" or report.get("mesh_published") is not True
            or any(report.get(key, {}).get("valid") is not True or report.get(key, {}).get("complete") is not True
                   for key in ("audit", "quality"))):
        raise ValueError("Saved success lacks its Native validity gate")
    signature, _ = signatures(directory / "mesh.paacsmesh", args, started)
    if case.get("mesh_sha256") != signature["archive"]["sha256"]:
        raise ValueError("Archive differs from saved case hash")
    for record in (case, report):
        for key, measured in (("nodes", signature["nodes"]), ("tetrahedra", signature["tetrahedra"])):
            if record.get(key) is not None and record[key] != measured:
                raise ValueError("Saved counts differ from actual archive")
    summary, summary_sha = read_json(directory / "quality-summary.json")
    if summary.get("mesh_sha256") != signature["archive"]["sha256"] or summary.get("node_count") != signature["nodes"]:
        raise ValueError("Inspector summary differs from the saved archive")
    if (summary.get("orientation_counts", {}).get("positive") != signature["tetrahedra"]
            or summary.get("orientation_counts", {}).get("inverted") != 0
            or summary.get("orientation_counts", {}).get("degenerate") != 0):
        raise ValueError("Saved inspector summary does not report every Tet4 positive")
    for metric in METRICS:
        values = summary.get("quality", {}).get(metric, {})
        if values.get("count") != signature["tetrahedra"] or any(
                not isinstance(values.get(key), (int, float)) or not math.isfinite(values[key])
                for key in ("minimum", "p05", "median")):
            raise ValueError("Incomplete retained quality distribution")
    result.update(archive=signature, quality=summary.get("quality"), achieved_sizing=summary.get("sizing"),
        thickness_capture=case.get("thickness_capture"),
        inspector={"artifact": artifact(directory / "quality-summary.json", summary_sha),
                   "metric_conventions": summary.get("metric_conventions"),
                   "scope": "Retained PAACS Python inspector measurements, hash-paired to archive; not a new Native CAD/intersection audit"})
    return result


def pair_summary(label, before, after):
    if before["source"] != after["source"]:
        raise ValueError("A before/after pair must use the identical owned source")
    result = {"label": label, "before": before, "after": after,
              "same_requested_size": before["requested"]["size_mm"] == after["requested"]["size_mm"],
              "interpretation": "Different counts, boundary sampling, achieved sizes or CAD error prevent a pure speed/quality equivalence claim."}
    if before.get("archive") and after.get("archive"):
        a, b = before["archive"], after["archive"]
        result["identities"] = {"same_final_ordered_points_and_classifications":
            a["final_all_points"]["ordered_xyz_and_classification_sha256"] == b["final_all_points"]["ordered_xyz_and_classification_sha256"],
            "same_final_point_multiset_and_classifications":
            a["final_all_points"]["multiset_xyz_and_classification_sha256"] == b["final_all_points"]["multiset_xyz_and_classification_sha256"],
            "same_exact_oriented_surface_and_face_tags": a["surface"]["oriented_coordinates_and_face_tags_sha256"] == b["surface"]["oriented_coordinates_and_face_tags_sha256"],
            "same_indexed_source_triangles_and_tags": a["surface"]["ordered_indexed_triangles_and_tags_sha256"] == b["surface"]["ordered_indexed_triangles_and_tags_sha256"],
            "same_complete_archive": a["archive"]["sha256"] == b["archive"]["sha256"],
            "initial_point_cloud_comparison_available": False}
    return result


def draw_figure(selections, records, output, args, started):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    spec = importlib.util.spec_from_file_location("initial_generation_quality", ROOT / "python/paacs_mesh/quality_metrics.py")
    quality = importlib.util.module_from_spec(spec); spec.loader.exec_module(quality)
    if len({records[d]["source"]["sha256"] for _, d in selections}) != 1:
        raise ValueError("One figure must use the same owned geometry for its cameras/cut")
    if any(records[d]["source"]["fixture"] != "tests/fixtures/cad/thin_plate.step" for _, d in selections):
        raise ValueError("The four-panel thin-part showcase supports the owned thin plate fixture")
    fig = plt.figure(figsize=(16, 11), facecolor=BG)
    fig.text(.035, .965, "PAACS  /  Initial Tet4 construction", fontsize=25, weight="bold", color=INK)
    fig.text(.035, .932, "Stored whole elements · common camera and cut · full-mesh red x-ray outliers", fontsize=12, color=SUB)
    panels = []
    bounds = None
    for index, (label, directory) in enumerate(selections):
        check_deadline(started, args)
        record = records[directory]
        if record.get("archive") is None:
            raise ValueError("Only Native-published saved meshes can be rendered")
        sig, geometry = signatures(directory / "mesh.paacsmesh", args, started, include_tets=True)
        if sig["archive"] != record["archive"]["archive"]:
            raise ValueError("Archive changed since comparison")
        points = np.asarray(geometry["points"], dtype=float)
        connection = np.asarray(geometry["tets"], dtype=np.int64)
        ids = np.asarray(geometry["tet_ids"], dtype=np.int64)
        corners = points[connection]
        if bounds is None:
            bounds = points.min(axis=0), points.max(axis=0)
        elif not np.array_equal(bounds, (points.min(axis=0), points.max(axis=0))):
            raise ValueError("Figure cases have different physical bounds")
        values = []
        for local, xyz in enumerate(corners):
            if local % 256 == 0:
                check_deadline(started, args)
            measured = quality.tet4_metrics(xyz)
            if measured["orientation"] != "positive" or measured["range_limited"]:
                raise ValueError("Figure encountered nonpositive/range-limited cell")
            values.append(measured["mean_ratio"])
        values = np.asarray(values)
        retained = record["quality"]["mean_ratio"]
        for key, actual in zip(("minimum", "p05", "median"), np.quantile(values, [0, .05, .5])):
            if not math.isclose(float(actual), retained[key], rel_tol=1e-10, abs_tol=1e-15):
                raise ValueError("Fresh figure metrics disagree with retained inspector summary")
        highlights = np.flatnonzero(values < args.highlight_below)
        if len(highlights) > args.max_highlights:
            raise ValueError("Full-mesh highlight allowance exceeded")
        lo, hi = bounds
        cut = float((lo[1] + hi[1]) / 2)
        kept = corners.mean(axis=1)[:, 1] >= cut
        exterior = {}
        for local in np.flatnonzero(kept):
            for face in FACES:
                vertices = tuple(int(connection[local, i]) for i in face); key = tuple(sorted(vertices))
                if key in exterior:
                    exterior[key][2] += 1
                else:
                    exterior[key] = [vertices, int(local), 1]
        if any(v[2] > 2 for v in exterior.values()):
            raise ValueError("Retained whole cells have nonmanifold incidence")
        exposed = sorted((k, v) for k, v in exterior.items() if v[2] == 1)
        if not exposed or len(exposed) > args.max_faces:
            raise ValueError("Whole-element exposed-face allowance exceeded")
        display = points.copy(); display[:, 2] = lo[2] + (points[:, 2] - lo[2]) * args.z_scale
        ax = fig.add_subplot(2, 2, index + 1, projection="3d", computed_zorder=False); ax.set_facecolor(BG)
        surface = [display[list(v)] for v, _ in geometry["triangles"]]
        ax.add_collection3d(Poly3DCollection(surface, facecolors="#8eb8c3", edgecolors="none", alpha=.07, zorder=0))
        polygons = [display[list(row[1][0])] for row in exposed]
        owners = np.asarray([row[1][1] for row in exposed])
        ax.add_collection3d(Poly3DCollection(polygons, facecolors=plt.get_cmap("viridis")(Normalize(0, .15)(values[owners])),
            edgecolors="#17364b80", linewidths=.24, zorder=1))
        red = [display[connection[local, list(face)]] for local in highlights for face in FACES]
        if red:
            ax.add_collection3d(Poly3DCollection(red, facecolors="#ea4d3a", edgecolors="#8c2118", linewidths=.8, alpha=.65, zorder=10))
            for local in highlights[:8]:
                center = display[connection[local]].mean(axis=0)
                ax.text(*center, f"#{ids[local]}", fontsize=7, color="#8c2118", zorder=20)
        spans = (hi - lo) * np.array([1, 1, args.z_scale]); padding = max(spans) * .015
        ax.set(xlim=(lo[0]-padding, hi[0]+padding), ylim=(lo[1]-padding, hi[1]+padding),
               zlim=(lo[2]-padding, lo[2]+spans[2]+padding))
        ax.set_box_aspect(spans); ax.view_init(elev=24, azim=-60); ax.set_proj_type("ortho"); ax.set_axis_off()
        ax.set_title(label, fontsize=14, color=INK, weight="bold", pad=-2)
        q = record["quality"]["mean_ratio"]
        times = record["timings"]
        native = (times["native_seconds"] or {}).get("meshing_total")
        def number(value):
            return "unavailable" if value is None else f"{value:.3g}"
        ax.text2D(.02, -.03, f"{len(ids):,} Tet4 · MR min {q['minimum']:.3g} / p05 {q['p05']:.3g} / median {q['median']:.3g}\n"
            f"Native {number(native)} s · API {number(times['api_call_seconds'])} s · process {number(times['process_seconds'])} s\n"
            f"{len(highlights)} red cells below {args.highlight_below:g}", transform=ax.transAxes, fontsize=9, color=SUB, linespacing=1.6)
        panels.append({"label": label, "case": record["case"], "archive": sig["archive"],
            "cut": {"axis": "y", "value": cut, "rule": "retain complete Tet4 whose centroid y >= cut"},
            "kept_cell_ids_sha256": digest(np.asarray(ids[kept], dtype="<i8").tobytes()),
            "exposed_face_owner_mapping_sha256": digest(np.asarray([(int(ids[v[1]]), *k) for k, v in exposed], dtype="<i8").tobytes()),
            "highlighted_cell_ids": [int(ids[i]) for i in highlights],
            "highlighted_mean_ratios": [float(values[i]) for i in highlights],
            "fresh_figure_metric_scope": "PAACS Python tet4_metrics from actual unscaled archive coordinates",
            "kept_tetrahedra": int(kept.sum()), "exposed_faces": len(exposed)})
    fig.subplots_adjust(left=.02, right=.90, bottom=.14, top=.88, wspace=.04, hspace=.34)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, .15), cmap="viridis"), cax=fig.add_axes([.925, .30, .014, .38]))
    bar.set_label("Mean ratio; scale 0–0.15 (ideal = 1)", fontsize=10, color=INK)
    remaining = records[selections[-1][1]]["quality"]["mean_ratio"]["minimum"]
    fig.text(.035, .100, f"Cleanup-on minimum MR {remaining:.3g} remains poor; thin anisotropy persists. Red marks numerical slivers < {args.highlight_below:g}, not all poor cells.", fontsize=10, color=SUB)
    fig.text(.035, .075, f"Every panel: actual Z extent {bounds[1][2]-bounds[0][2]:g} mm; Z enlarged ×{args.z_scale:g} for display only. Stored whole cells only, never clipped or rebuilt.", fontsize=10, color=SUB)
    fig.text(.035, .050, "Red x-ray cells include the full mesh, including hidden/removed-half cells. All measurements use unscaled coordinates.", fontsize=10, color=SUB)
    fig.text(.035, .025, "Saved Native checks; no CAD audit rerun. Zero red cells and improved shape do not establish solver adequacy. Timings are single attempts.", fontsize=10, color=SUB)
    figure = output / "initial-generation-thin.png"
    fig.savefig(figure, dpi=args.dpi, facecolor=BG); plt.close(fig)
    return {"schema": "paacs-initial-generation-figure-v1", "figure": {"path": figure.name, "sha256": digest(figure.read_bytes())},
        "camera": {"elevation_degrees": 24, "azimuth_degrees": -60, "projection": "orthographic"},
        "display_z_scale": args.z_scale, "measurements_use_actual_unscaled_coordinates": True,
        "highlight_below_mean_ratio": args.highlight_below, "panels": panels,
        "quality_source_sha256": digest((ROOT / "python/paacs_mesh/quality_metrics.py").read_bytes())}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair", nargs=3, action="append", metavar=("LABEL", "BEFORE_DIR", "AFTER_DIR"), required=True)
    parser.add_argument("--figure-case", nargs=2, action="append", default=[], metavar=("LABEL", "CASE_DIR"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-seconds", type=float, default=240)
    parser.add_argument("--max-bytes", type=int, default=128 * 1024**2)
    parser.add_argument("--max-working-bytes", type=int, default=1024**3)
    parser.add_argument("--max-nodes", type=int, default=100_000)
    parser.add_argument("--max-cells", type=int, default=500_000)
    parser.add_argument("--max-plot-tets", type=int, default=20_000)
    parser.add_argument("--max-faces", type=int, default=30_000)
    parser.add_argument("--max-highlights", type=int, default=2048)
    parser.add_argument("--max-report-bytes", type=int, default=2 * 1024**2)
    parser.add_argument("--highlight-below", type=float, default=1e-6)
    parser.add_argument("--z-scale", type=float, default=40)
    parser.add_argument("--dpi", type=int, default=160)
    args = parser.parse_args(argv)
    if not 1 <= len(args.pair) <= 12 or len(args.figure_case) not in (0, 4):
        parser.error("Require1–12 pairs and either zero or four figure cases")
    for name in ("max_seconds", "max_bytes", "max_working_bytes", "max_nodes", "max_cells", "max_plot_tets",
                 "max_faces", "max_highlights", "max_report_bytes", "z_scale", "dpi"):
        if not math.isfinite(getattr(args, name)) or getattr(args, name) <= 0:
            parser.error(f"{name} must be positive and finite")
    if not 0 < args.highlight_below < 1 or 16 * 11 * args.dpi**2 > 8_000_000:
        parser.error("Invalid highlight threshold or image pixel allowance exceeded")
    output = args.output.resolve()
    if output.exists() or not output.parent.is_dir():
        parser.error("Output must be a new directory under an existing parent")
    started, cache, records = time.perf_counter(), {}, {}
    paths = {Path(p).resolve() for _, a, b in args.pair for p in (a, b)}
    paths.update(Path(p).resolve() for _, p in args.figure_case)
    for path in sorted(paths):
        records[path] = load_case(path, args, started, cache)
    comparisons = [pair_summary(label, records[Path(a).resolve()], records[Path(b).resolve()]) for label, a, b in args.pair]
    result = {"schema": "paacs-initial-generation-comparison-v1", "pairs": comparisons,
        "scope": "Owned-fixture saved runs; final identities plus retained inspector quality, achieved sizing and Native checks. No CAD/solver audits rerun.",
        "timing_note": "Single saved attempts; not medians. All process costs use exact source/archive/case-paired process receipts.",
        "initial_point_cloud_note": "Final-point signatures do not certify the initial admitted point cloud; absent stage snapshots remain unavailable.",
        "script_sha256": digest(Path(__file__).read_bytes())}
    payload = json.dumps(result, indent=2, allow_nan=False).encode() + b"\n"
    if len(payload) > args.max_report_bytes:
        raise ValueError("Comparison report output allowance exceeded")
    with tempfile.TemporaryDirectory(prefix=".initial-generation-", dir=output.parent) as temporary:
        staging = Path(temporary) / "output"; staging.mkdir()
        (staging / "comparison.json").write_bytes(payload)
        columns = ("pair", "side", "case", "status", "tetrahedra", "api_seconds", "native_seconds", "process_seconds",
                   "mean_ratio_minimum", "mean_ratio_p05", "mean_ratio_median", "scaled_jacobian_minimum", "scaled_jacobian_p05",
                   "scaled_jacobian_median", "dihedral_minimum_degrees", "dihedral_p05_degrees", "dihedral_median_degrees")
        with (staging / "comparison.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns); writer.writeheader()
            for pair in comparisons:
                for side in ("before", "after"):
                    r = pair[side]; q = r.get("quality", {}); timing = r["timings"]
                    writer.writerow({"pair": pair["label"], "side": side, "case": r["case_id"], "status": r["status"],
                        "tetrahedra": (r.get("archive") or {}).get("tetrahedra"), "api_seconds": timing["api_call_seconds"],
                        "native_seconds": (timing["native_seconds"] or {}).get("meshing_total"), "process_seconds": timing["process_seconds"],
                        "mean_ratio_minimum": q.get("mean_ratio", {}).get("minimum"), "mean_ratio_p05": q.get("mean_ratio", {}).get("p05"),
                        "mean_ratio_median": q.get("mean_ratio", {}).get("median"),
                        "scaled_jacobian_minimum": q.get("scaled_jacobian", {}).get("minimum"),
                        "scaled_jacobian_p05": q.get("scaled_jacobian", {}).get("p05"),
                        "scaled_jacobian_median": q.get("scaled_jacobian", {}).get("median"),
                        "dihedral_minimum_degrees": q.get("minimum_dihedral_deg", {}).get("minimum"),
                        "dihedral_p05_degrees": q.get("minimum_dihedral_deg", {}).get("p05"),
                        "dihedral_median_degrees": q.get("minimum_dihedral_deg", {}).get("median")})
        if args.figure_case:
            provenance = draw_figure([(label, Path(p).resolve()) for label, p in args.figure_case], records, staging, args, started)
            encoded = json.dumps(provenance, indent=2, allow_nan=False).encode() + b"\n"
            if len(encoded) > args.max_report_bytes:
                raise ValueError("Figure provenance output allowance exceeded")
            (staging / "figure-provenance.json").write_bytes(encoded)
        check_deadline(started, args)
        staging.rename(output)
    print(json.dumps({"output": relative(output), "pairs": len(comparisons), "figures": int(bool(args.figure_case)),
                      "analysis_seconds": time.perf_counter() - started}))


if __name__ == "__main__":
    main()
