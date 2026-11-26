"""Inspect accepted qualification archives and render actual whole-element cuts.

Requires Matplotlib for figures. No native meshing or CAD module is imported.
Independent file checks do not certify global intersections or solver adequacy.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
from collections import defaultdict
from contextlib import closing
from fractions import Fraction
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import statistics
import sys
import time
import types

ROOT = Path(__file__).resolve().parents[2]
FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3))


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load(path):
    require(Path(path).stat().st_size <= 16*1024**2, "Receipt exceeds 16 MiB")
    return json.loads(Path(path).read_text(encoding="utf-8-sig"),
                      parse_constant=lambda value: require(False, "Nonfinite JSON: "+value))


class Limits:
    def __init__(self, args):
        self.args, self.started, self.work, self.retained, self.peak = args, time.perf_counter(), 0, 0, 0

    def tick(self, count=1):
        self.work += count
        require(self.work <= self.args.max_work, "Inspection work allowance exhausted")
        require(time.perf_counter()-self.started <= self.args.max_seconds, "Cooperative inspection deadline exhausted")

    def admit(self, extra):
        self.peak = max(self.peak, self.retained+extra)
        require(self.peak <= self.args.max_working_bytes, "Modeled working payload allowance exhausted")


def modules():
    package = types.ModuleType("_paacs_thin_hollow_file")
    package.__path__ = [str(ROOT/"python/paacs_mesh")]
    sys.modules[package.__name__] = package
    archive = importlib.import_module(package.__name__+".mesh_archive")
    metrics = importlib.import_module(package.__name__+".quality_metrics")

    class FileArchive(archive.ArchiveView):
        def _load_source(self, source):
            self.original_path = Path(source).resolve(strict=True)
            require(self.original_path.stat().st_size <= self.max_bytes//2, "Archive snapshot exceeds byte allowance")
            self.path = Path(self._temporary.name)/"mesh.paacsmesh"
            with self.original_path.open("rb") as incoming, self.path.open("xb") as outgoing:
                remaining = self.original_path.stat().st_size
                while remaining:
                    self._deadline()
                    block = incoming.read(min(1024**2, remaining))
                    require(block, "Archive shrank during snapshot")
                    outgoing.write(block)
                    remaining -= len(block)
                require(not incoming.read(1), "Archive grew during snapshot")
            self.protected_paths.append(self.original_path)
    return archive, metrics, FileArchive


def metadata(path):
    path = Path(path).resolve(strict=True)
    receipt = load(path)
    job, worker = receipt["job"], receipt["worker"]
    require(receipt.get("state") == "complete" and receipt.get("status") == "native_accepted"
            and receipt.get("worker_pid_verified") is True and not receipt.get("hard_timeout")
            and receipt.get("input_guards_unchanged") is True and worker.get("input_guards_unchanged") is True,
            "Controller is not a completed, guarded Native acceptance")
    report_path = Path(job["output"])/worker["native_report_file"]
    require(sha(report_path) == worker["native_report_sha256"], "Native report hash changed")
    report = load(report_path)
    require(report.get("status") == "succeeded" and report.get("mesh_published") is True
            and all(report.get(s, {}).get(k) is True for s in ("audit", "quality") for k in ("complete", "valid")),
            "Native publication/audit/quality gates did not pass")
    source, archive = Path(job["source"]).resolve(strict=True), Path(job["output"])/"mesh.paacsmesh"
    require(sha(source) == job["source_sha256"] == worker["source_sha256"] == report["source_sha256"], "Source binding changed")
    require(sha(archive) == worker["archive_sha256"] and worker["native_sha256"] == job["native_sha256"], "Archive/native identity changed")
    bodies, faces = worker["original_body_tags"], worker["original_face_tags"]
    require(len(bodies) == 1 and faces and all(isinstance(t, int) and not isinstance(t, bool) and t > 0 for t in bodies+faces), "One original body with positive face tags required")
    guards = {str(p): sha(p) for p in (path, report_path, source, archive)}
    manifests = [Path(p) for p in job["input_sha256"] if Path(p).name == "manifest.json"]
    require(len(manifests) == 1, "Unique owned fixture manifest unavailable")
    manifest = manifests[0]
    require(sha(manifest) == job["input_sha256"][str(manifest)], "Owned fixture manifest changed")
    require(sum(f == job["specification"] for f in load(manifest)["fixtures"]) == 1, "Job specification differs from owned manifest")
    guards[str(manifest)] = sha(manifest)
    return {"controller": str(path), "job": job, "worker": worker, "report": report, "archive": archive,
            "bodies": bodies, "faces": sorted(set(faces)), "guards": guards}


def winding(ids):
    return -1 if sum(ids[i] > ids[j] for i in range(3) for j in range(i+1, 3)) % 2 else 1


def determinant(points):
    a, b, c, d = [[Fraction(x) for x in p] for p in points]
    u, v, w = [[x-y for x, y in zip(p, a)] for p in (b, c, d)]
    return u[0]*(v[1]*w[2]-v[2]*w[1])-u[1]*(v[0]*w[2]-v[2]*w[0])+u[2]*(v[0]*w[1]-v[1]*w[0])


def distribution(values):
    values = sorted(values)
    return {"minimum": values[0], "p01": values[int((len(values)-1)*.01)], "median": statistics.median(values),
            "p95": values[int((len(values)-1)*.95)], "maximum": values[-1], "count": len(values)}


def inspect(meta, archive_module, metrics, FileArchive, limits):
    args, begun = limits.args, time.perf_counter()
    with closing(FileArchive(meta["archive"], max_nodes=args.max_nodes, max_cells=args.max_cells,
                             max_bytes=args.max_archive_bytes, max_work=args.max_work, max_seconds=args.max_seconds)) as source:
        require(source.mesh_sha256 == meta["worker"]["archive_sha256"], "Immutable archive snapshot differs from accepted receipt")
        payload = source.node_count*256+source.cell_count*4096+16*1024**2
        limits.admit(payload)
        limits.tick(source.node_count+source.cell_count+sum(n*archive_module.CELL_WIDTHS[k] for k, n in enumerate(source.type_counts)))
        nodes, classes = [], []
        for i in range(source.node_count):
            limits.tick(); point, cls = source.node(i)
            require(cls["tag"] > 0 and (cls["dimension"] != 2 or cls["tag"] in meta["faces"])
                    and (cls["dimension"] != 3 or cls["tag"] == meta["bodies"][0]), "Node CAD classification differs")
            nodes.append(point); classes.append(cls)
        cells, faces, boundary, exact_volume = [], defaultdict(list), {}, Fraction(0)
        for cell_id in range(source.cell_count):
            kind, ids, cls = source.cell(cell_id); limits.tick(1+len(ids))
            require(len(set(ids)) == len(ids), "Repeated element node index")
            if archive_module.CELL_NAMES[kind] == "tet4":
                require(cls == {"dimension": 3, "tag": meta["bodies"][0]}, "Tet4 original body differs")
                xyz = [nodes[i] for i in ids]; det = determinant(xyz)
                require(det > 0, "Exact binary-coordinate orientation is nonpositive")
                exact_volume += det/6; measured = metrics.tet4_metrics(xyz)
                require(measured["orientation"] == "positive" and not measured["range_limited"]
                        and all(isinstance(v, (int, float, str, bool)) and (not isinstance(v, float) or math.isfinite(v)) for v in measured.values()), "Unrepresentable element metric")
                cells.append({"cell_id": cell_id, "node_ids": list(ids), "body_tag": cls["tag"],
                              "centroid": [math.fsum(p[j]*.25 for p in xyz) for j in range(3)], "metrics": measured})
                for local in FACES:
                    oriented = tuple(ids[j] for j in local); faces[tuple(sorted(oriented))].append((cell_id, oriented))
            elif archive_module.CELL_NAMES[kind] == "tri3":
                require(cls["dimension"] == 2 and cls["tag"] in meta["faces"], "Tri3 original face differs")
                key = tuple(sorted(ids)); require(key not in boundary, "Duplicate boundary triangle")
                boundary[key] = (tuple(ids), cls["tag"])
            else:
                raise ValueError("Only single-body Tet4/Tri3 archives are supported")
        require(cells, "No tetrahedra")
        exposed = set()
        for key, owners in faces.items():
            limits.tick(); require(1 <= len(owners) <= 2, "Nonmanifold tetrahedron face")
            if len(owners) == 2:
                require(winding(owners[0][1]) != winding(owners[1][1]), "Internal face winding differs")
            else:
                exposed.add(key)
                require(key in boundary and winding(owners[0][1]) == winding(boundary[key][0]), "Exterior winding or connectivity differs")
        require(exposed == set(boundary) and sorted({tag for _, tag in boundary.values()}) == meta["faces"], "Exterior or original CAD face tags missing")
    native_quality = meta["report"]["quality"]
    minimum = min(c["metrics"]["mean_ratio"] for c in cells)
    native_minimum = native_quality.get("minimum_quality", native_quality.get("minimum"))
    require(isinstance(native_minimum, (int, float)) and math.isfinite(native_minimum), "Native minimum quality unavailable")
    agrees = abs(minimum-native_minimum) <= max(1e-12, abs(native_minimum)*1e-9)
    volume = float(exact_volume); native_volume = native_quality["total_volume"]
    require(abs(volume-native_volume) <= max(1e-8, abs(native_volume)*1e-10), "Archive and Native volumes differ")
    summary = {"case": meta["job"]["id"], "profile": meta["job"]["profile"], "source_sha256": meta["job"]["source_sha256"],
               "archive_sha256": meta["worker"]["archive_sha256"], "native_sha256": meta["job"]["native_sha256"], "input_sha256": meta["guards"],
               "nodes": len(nodes), "tet4": len(cells), "tri3": len(boundary), "body_tags": meta["bodies"], "face_tags": meta["faces"],
               "units": meta["report"].get("units", "unknown"), "fixture_specification": meta["job"]["specification"],
               "requested_h_mm": meta["job"]["global_h_mm"], "exact_positive_orientations": len(cells), "opposite_internal_faces": True,
               "exterior_matches_tri3": True, "exact_volume_mm3": volume, "native_quality_binding": {"native_minimum": native_minimum, "independent_minimum": minimum, "agrees": agrees},
               "quality": {k: distribution(c["metrics"][k] for c in cells) for k in ("mean_ratio", "minimum_dihedral_deg", "scaled_jacobian", "edge_median", "edge_max")},
               "below_mean_ratio_0_1": sum(c["metrics"]["mean_ratio"] < .1 for c in cells), "mesh_call_seconds": meta["worker"]["mesh_call_seconds"],
               "controller_seconds": load(meta["controller"])["process_seconds"], "file_inspection_seconds": time.perf_counter()-begun, "cells": cells}
    limits.retained += payload
    return {"summary": summary, "nodes": nodes, "faces": faces, "by_id": {c["cell_id"]: c for c in cells}}


def fixture_display(data, limits):
    """Validate the pinned rigid placement, then change display coordinates only."""
    spec = data["meta"]["job"]["specification"]
    matrix = spec.get("placement_matrix")
    require(isinstance(matrix, list) and len(matrix) == 4 and all(isinstance(row, list) and len(row) == 4 for row in matrix), "Fixture placement must be a4by4 matrix")
    require(all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for row in matrix for v in row), "Fixture placement must contain finite numbers")
    matrix = [[float(v) for v in row] for row in matrix]
    require(max(abs(matrix[3][j]-(1. if j == 3 else 0.)) for j in range(4)) <= 1e-12, "Fixture placement has a nonhomogeneous last row")
    rotation = [row[:3] for row in matrix[:3]]
    error = max(abs(math.fsum(rotation[k][i]*rotation[k][j] for k in range(3))-(1. if i == j else 0.)) for i in range(3) for j in range(3))
    a, b, c = rotation
    det = a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
    require(error <= 1e-10 and det > 0 and abs(det-1.) <= 3e-10, "Fixture placement is not a right-handed orthonormal rigid transform")
    translation = [matrix[i][3] for i in range(3)]
    def inverse(point):
        limits.tick()
        result = [math.fsum(rotation[k][j]*(point[k]-translation[k]) for k in range(3)) for j in range(3)]
        require(all(math.isfinite(v) for v in result), "Nonfinite fixture display coordinate")
        return result
    require(spec.get("profile") in ("round", "square", "rectangle", "rounded-rectangle"), "Unknown owned profile bounds")
    width, height, length = [spec.get(k) for k in ("outer_width_mm", "outer_height_mm", "length_mm")]
    require(all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) and v > 0 for v in (width, height, length)), "Invalid owned profile dimensions")
    nodes = [inverse(p) for p in data["nodes"]]
    centers = {cell["cell_id"]: inverse(cell["centroid"]) for cell in data["summary"]["cells"]}
    lo = [min(p[j] for p in nodes) for j in range(3)]; hi = [max(p[j] for p in nodes) for j in range(3)]
    expected = [[-width/2, -height/2, 0.], [width/2, height/2, length]]
    chord = data["meta"]["job"]["kwargs"].get("max_chord_error")
    if chord is None:
        chord = max((b.get("surface_chord_error_limit", 0.) for b in data["meta"]["report"].get("bodies", [])), default=0.)
    require(isinstance(chord, (int, float)) and math.isfinite(chord) and chord >= 0, "Finite display-bound chord allowance unavailable")
    # A circular tessellation need not contain cardinal extrema. The declared
    # CAD chord limit bounds that display-box difference; it does not weaken
    # the original orientation, topology or geometry admission checks.
    tolerance = max(1e-6, max(width, height, length)*1e-8, min(chord, .01*min(width, height, length)))
    require(all(abs(actual[j]-target[j]) <= tolerance for actual, target in zip((lo, hi), expected) for j in range(3)), "Transformed archive bounds disagree with the pinned fixture dimensions")
    inverse_matrix = [[rotation[k][j] for k in range(3)]+[-math.fsum(rotation[k][j]*translation[k] for k in range(3))] for j in range(3)]+[[0.,0.,0.,1.]]
    frame = {"fixture_to_world": matrix, "world_to_fixture": inverse_matrix,
             "relation": "p_fixture = transpose(R) * (p_archive_world - translation)",
             "orthonormal_max_error": error, "determinant": det, "expected_fixture_bounds_mm": expected,
             "measured_fixture_bounds_mm": [lo, hi], "bounds_tolerance_mm": tolerance,
             "metric_coordinates": "original archive world coordinates; not recomputed for display", "ids": "original archive cell/node/body IDs"}
    return nodes, centers, frame


def render(datasets, labels, output, limits):
    os.environ.setdefault("MPLCONFIGDIR", str(output/".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    args = limits.args; cmap, norm = plt.get_cmap("viridis"), Normalize(0, 1)
    extra = int(18*args.dpi)*int(5.2*len(datasets)*args.dpi)*16+64*1024**2
    limits.admit(extra); fig = plt.figure(figsize=(18, 5.2*len(datasets)), facecolor="#f3f6f9")
    manifests, frames, rendered = [], [], 0
    try:
        fig.suptitle("PAACS Mesh / owned hollow STEP extrusions", x=.035, ha="left", y=.985, fontsize=20, weight="bold")
        for row, data in enumerate(datasets):
            summary = data["summary"]
            extra += 128*(len(data["nodes"])+len(summary["cells"])); limits.admit(extra)
            nodes, centers, frame = fixture_display(data, limits); frames.append(frame)
            lo = [min(p[j] for p in nodes) for j in range(3)]; hi = [max(p[j] for p in nodes) for j in range(3)]; span = [hi[j]-lo[j] for j in range(3)]
            require(min(span) > 0, "Degenerate rendering bounds")
            cells = summary["cells"]; all_ids = {c["cell_id"] for c in cells}
            expected = frame["expected_fixture_bounds_mm"]
            cut_y, cut_z = [(expected[0][j]+expected[1][j])/2 for j in (1, 2)]
            selections = (all_ids, {c["cell_id"] for c in cells if centers[c["cell_id"]][1] <= cut_y},
                          {c["cell_id"] for c in cells if centers[c["cell_id"]][2] <= cut_z})
            for column, selected in enumerate(selections):
                polygons, colors, visible = [], [], []
                for key, owners in sorted(data["faces"].items()):
                    limits.tick(); chosen = [o for o in owners if o[0] in selected]
                    if len(chosen) != 1 or (column == 2 and len(owners) == 1):
                        continue
                    rendered += 1; require(rendered <= args.max_render_faces, "Total rendered-face allowance exhausted")
                    extra += 2048; limits.admit(extra)
                    cell_id, ids = chosen[0]; polygons.append([nodes[i] if column < 2 else nodes[i][:2] for i in ids])
                    colors.append(cmap(norm(data["by_id"][cell_id]["metrics"]["mean_ratio"])))
                    visible.append({"cell_id": cell_id, "node_ids": list(ids), "body_tag": summary["body_tags"][0], "original_exterior": len(owners) == 1})
                ax = fig.add_subplot(len(datasets), 3, row*3+column+1, projection="3d" if column < 2 else None, facecolor="#f3f6f9")
                if column < 2:
                    ax.add_collection3d(Poly3DCollection(polygons, facecolors=colors, edgecolors="#335361", linewidths=.15))
                    for j, setter in enumerate((ax.set_xlim, ax.set_ylim, ax.set_zlim)): setter(lo[j]-.025*span[j], hi[j]+.025*span[j])
                    ax.set_box_aspect(span); ax.view_init(elev=28, azim=35); ax.set_axis_off()
                else:
                    order = sorted(range(len(visible)), key=lambda i: (math.fsum(nodes[n][2] for n in visible[i]["node_ids"])/3, visible[i]["cell_id"]))
                    polygons, colors, visible = ([v[i] for i in order] for v in (polygons, colors, visible))
                    ax.add_collection(PolyCollection(polygons, facecolors=colors, edgecolors="#335361", linewidths=.2))
                    ax.set_xlim(lo[0]-.04*span[0], hi[0]+.04*span[0]); ax.set_ylim(lo[1]-.04*span[1], hi[1]+.04*span[1]); ax.set_aspect("equal")
                    ax.set_xlabel("fixture x / mm"); ax.set_ylabel("fixture y / mm")
                title = ("Full mesh / fixture frame", "Whole-Tet fixture-y half-cut", "Exposed fixture-z cut / end projection")[column]
                ax.set_title(title, fontsize=11)
                manifests.append({"dataset": row, "coordinate_frame": "fixture", "display_frame_index": row,
                                  "kind": ("full", "centroid_y_half", "centroid_z_half_projected_xy")[column],
                                  "cut_threshold_fixture_mm": None if column == 0 else (cut_y if column == 1 else cut_z),
                                  "selected_cell_ids": sorted(selected), "visible_faces": visible})
            q = summary["quality"]; text = f'{labels[row]} / {summary["tet4"]:,} Tet4 / call {summary["mesh_call_seconds"]:.2f}s / h {summary["requested_h_mm"]:g}mm; achieved median edge {q["edge_median"]["median"]:.3g}mm\nMR min/median {q["mean_ratio"]["minimum"]:.3g}/{q["mean_ratio"]["median"]:.3f}; minimum dihedral {q["minimum_dihedral_deg"]["minimum"]:.3g} degrees'
            fig.text(.035, .93-row*.80/len(datasets), text, fontsize=10, weight="bold")
        fig.subplots_adjust(left=.025, right=.98, top=.86, bottom=.13, hspace=.45, wspace=.16)
        bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=fig.add_axes([.31, .065, .38, .012]), orientation="horizontal")
        bar.set_label("Tet4 mean ratio / ideal = 1 / fixed scale for every mesh", fontsize=9)
        fig.text(.035, .02, "Validated fixture-frame display; world-coordinate metrics and archive IDs retained. Whole elements and projected end faces; no new CAD/intersection/solver proof.", fontsize=8)
        path = output/"hollow-meshes.png"; temporary = output/"hollow-meshes.tmp.png"
        fig.savefig(temporary, dpi=args.dpi, metadata={"Software": "PAACS file-only owned qualification figures"})
        require(temporary.stat().st_size < args.max_output_bytes, "Figure output allowance exhausted"); os.replace(temporary, path)
        return {"path": str(path), "sha256": sha(path), "bytes": path.stat().st_size, "panels": manifests, "display_frames": frames, "quality_scale": [0, 1]}
    finally:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controllers", nargs="+", type=Path, required=True); parser.add_argument("--labels", nargs="+")
    parser.add_argument("--output", type=Path, required=True); parser.add_argument("--max-seconds", type=float, default=240)
    for name, default in (("max-nodes",250000),("max-cells",500000),("max-work",50000000),("max-working-bytes",2*1024**3),
                          ("max-output-bytes",128*1024**2),("max-archive-bytes",128*1024**2),("max-render-faces",100000),("dpi",180)):
        parser.add_argument("--"+name, type=int, default=default)
    args = parser.parse_args(); require(1 <= len(args.controllers) <= 6, "Supply one to six accepted controllers")
    require(args.labels is None or len(args.labels) == len(args.controllers), "One label per controller required")
    require(args.labels is None or all(0 < len(label) <= 100 and "\n" not in label for label in args.labels), "Labels limited to 100 single-line characters")
    require(all(v > 0 and math.isfinite(v) for k, v in vars(args).items() if k.startswith("max_")) and 60 <= args.dpi <= 300, "Invalid resource bounds")
    output = args.output.resolve(); require(not output.exists(), "Choose a fresh output directory")
    limits = Limits(args); archive, metrics, FileArchive = modules()
    code = {str(Path(p).resolve()): sha(p) for p in (__file__, archive.__file__, metrics.__file__, ROOT/"python/paacs_mesh/_cad_common.py")}
    datasets = []
    for path in args.controllers:
        meta = metadata(path); limits.tick(0); data = inspect(meta, archive, metrics, FileArchive, limits); data["meta"] = meta; datasets.append(data)
    output.mkdir(parents=True, exist_ok=False)
    figure = render(datasets, args.labels or [d["summary"]["case"] for d in datasets], output, limits)
    require(not any(n.endswith("._native") for n in sys.modules), "Native module unexpectedly imported")
    require(all(sha(p) == digest for d in datasets for p, digest in d["meta"]["guards"].items()) and all(sha(p) == h for p, h in code.items()), "An input or implementation changed")
    agrees = all(d["summary"]["native_quality_binding"]["agrees"] for d in datasets)
    report = {"schema": "paacs-owned-hollow-file-inspection-v1", "native_execution": False, "native_quality_binding_agrees": agrees,
              "datasets": [d["summary"] for d in datasets], "figure": figure, "implementation_sha256": code, "metric_conventions": metrics.METRIC_CONVENTIONS,
              "elapsed_seconds": time.perf_counter()-limits.started, "work": limits.work, "maximum_modeled_payload_bytes": limits.peak,
              "scope": "Every Tet4 measured with exact orientation and boundary-incidence/tag checks. Native CAD admission inherited. No global-intersection, continuous CAD-error or solver adequacy proof."}
    temporary = output/"inspection.json.tmp"; count = figure["bytes"]+sum(p.stat().st_size for p in (output/".matplotlib").rglob("*") if p.is_file())
    with temporary.open("wb") as stream:
        for part in json.JSONEncoder(separators=(",", ":"), allow_nan=False).iterencode(report):
            encoded = part.encode(); count += len(encoded); require(count <= args.max_output_bytes, "Aggregate output allowance exhausted"); stream.write(encoded)
    os.replace(temporary, output/"inspection.json"); print(json.dumps({"figure": str(output/"hollow-meshes.png"), "inspection": str(output/"inspection.json"), "native_execution": False, "quality_binding_agrees": agrees}))
    return 0 if agrees else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(json.dumps({"status": "refused", "native_execution": False, "error": str(error)}), file=sys.stderr)
        raise SystemExit(2)
