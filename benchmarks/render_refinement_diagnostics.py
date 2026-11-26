"""Render retained refinement_diagnostics.py results without rerunning meshing.

Requires the development mesh-inspection API. Standard-library SVG and a
self-contained HTML/SVG gallery are always produced. CairoSVG is optional for
PNG output; it is not a PAACS Mesh runtime dependency. Native reports are copied
verbatim and labelled saved reports, without a new native audit.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
import base64
from html import escape
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import tempfile
import time


CASES = ("uniform-direct", "wall-0p2-direct", "wall-0p08-direct", "wall-0p1-direct")
MIB = 1024 * 1024
SAVED_BASIS = "hash_paired_saved_native_report_not_recertified"


def _pairs(values):
    result = {}
    for key, value in values:
        if key in result:
            raise ValueError(f"Duplicate JSON field: {key}")
        result[key] = value
    return result


def _finite_float(value):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Nonfinite JSON number")
    return result


def _constant(value):
    raise ValueError(f"Nonfinite JSON constant: {value}")


class Budget:
    def __init__(self, args):
        self.args = args
        self.started = time.perf_counter()
        self.input_bytes = self.output_bytes = 0

    def check(self):
        if time.perf_counter() - self.started > self.args.total_seconds:
            raise ValueError("Rendering exceeded its cooperative total-seconds allowance")

    def read(self, path, *, metadata=False):
        self.check()
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"Input must be a regular local file: {path.name}")
        size = path.stat().st_size
        if metadata and size > min(MIB, self.args.max_bytes // 64):
            raise ValueError("Metadata exceeds its bounded JSON parsing allowance")
        if size > self.args.max_input_bytes - self.input_bytes:
            raise ValueError("Selected inputs exceed max-input-bytes")
        with path.open("rb") as stream:
            data = stream.read(size + 1)
        if len(data) != size:
            raise ValueError("Input changed while it was read")
        self.input_bytes += size
        self.check()
        return data

    def json(self, path):
        result = json.loads(self.read(path, metadata=True), object_pairs_hook=_pairs,
                            parse_constant=_constant, parse_float=_finite_float)
        if not isinstance(result, dict):
            raise ValueError("Expected a JSON object")
        return result

    def write(self, path, data):
        self.check()
        if len(data) > self.args.max_output_bytes - self.output_bytes:
            raise ValueError("Generated files exceed max-output-bytes")
        with path.open("xb") as stream:
            stream.write(data)
        self.output_bytes += len(data)
        self.check()

    def copy_archive(self, source, destination):
        self.check()
        if source.is_symlink() or not source.is_file():
            raise ValueError("A retained archive must be a regular local file")
        size = source.stat().st_size
        if size > min(self.args.max_input_bytes - self.input_bytes,
                      self.args.max_output_bytes - self.output_bytes):
            raise ValueError("Archive exceeds the remaining input/output byte allowance")
        digest = hashlib.sha256()
        with source.open("rb") as incoming, destination.open("xb") as outgoing:
            remaining = size
            while remaining:
                self.check()
                data = incoming.read(min(MIB, remaining))
                if not data:
                    raise ValueError("Archive changed during copy")
                outgoing.write(data)
                digest.update(data)
                remaining -= len(data)
            if incoming.read(1):
                raise ValueError("Archive grew during copy")
        self.input_bytes += size
        self.output_bytes += size
        return digest.hexdigest()

    def generated(self, path):
        size = path.stat().st_size
        if size > self.args.max_output_bytes - self.output_bytes:
            raise ValueError("Generated image exceeds the remaining output allowance")
        self.output_bytes += size
        self.check()

    def encoded_json(self, value):
        return (json.dumps(value, indent=2, allow_nan=False, sort_keys=True) + "\n").encode("utf-8")


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(MIB), b""):
            digest.update(data)
    return digest.hexdigest()


def _padded(extent, fraction):
    # Normalize before subtracting or padding large finite placements.
    magnitude = max(abs(value) for point in extent for value in point) or 1.
    local = [[value / magnitude for value in point] for point in extent]
    largest = max(b - a for a, b in zip(*local)) or 1.
    result = [[(value + sign * fraction * largest) * magnitude for value in point]
              for point, sign in zip(local, (-1, 1))]
    if not all(math.isfinite(value) for point in result for value in point):
        raise ValueError("Padded camera frame exceeds the finite coordinate range; supply explicit bounds")
    return result


def _gallery(images, root, budget):
    svg_bytes = sum((root / row["path"]).stat().st_size for row in images)
    if 12 * svg_bytes + MIB > budget.args.max_bytes:
        raise ValueError("Gallery encoding exceeds max-bytes; increase the explicit allowance")
    encoded = [base64.b64encode((root / row["path"]).read_bytes()).decode("ascii") for row in images]
    width, height = 2464, 2150
    title = "PAACS Mesh · retained refinement diagnostics"
    subtitle = "Actual archived geometry and native IDs. Shared full-view frame; final panel is an explicit close-up."
    note = "Saved native reports, not independently re-audited. X-ray overlays include every selected Tet4; gold uses original whole-cell faces."
    metadata = {"schema": "paacs-refinement-gallery-v1", "native_reaudit_performed": False,
                "images": [{key: row[key] for key in ("path", "mesh_sha256", "source_sha256", "units", "highlight_count", "svg_sha256")} for row in images]}
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           f'<title>{escape(title)}</title><metadata>{escape(json.dumps(metadata, sort_keys=True))}</metadata>',
           '<rect width="100%" height="100%" fill="#f3f6fa"/>',
           f'<text x="24" y="53" font-family="sans-serif" font-size="32" font-weight="bold" fill="#17364b">{escape(title)}</text>',
           f'<text x="24" y="82" font-family="sans-serif" font-size="18" fill="#17364b">{escape(subtitle)}</text>',
           f'<text x="24" y="108" font-family="sans-serif" font-size="16" fill="#17364b">{escape(note)}</text>']
    html = ['<!doctype html><html><head><meta charset="utf-8"><title>' + escape(title) + '</title>',
            '<style>body{margin:0;padding:24px;background:#f3f6fa;color:#17364b;font-family:Arial,sans-serif}h1{margin:0 0 8px}p{margin:5px 0}.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px;margin-top:24px}article{background:white;border:1px solid #c4d1db}h2{font-size:18px;margin:12px 16px 4px}article p{font-size:14px;margin:0 16px 8px}img{display:block;width:100%;height:auto}footer{margin-top:18px}@media(max-width:800px){.grid{grid-template-columns:1fr}}</style></head><body>',
            f'<h1>{escape(title)}</h1><p>{escape(subtitle)}</p><p>{escape(note)}</p><section class="grid">']
    for index, (row, data) in enumerate(zip(images, encoded)):
        x, y = 24 + (index % 2) * 1216, 134 + (index // 2) * 980
        heading, caption = row["heading"], row["caption"]
        svg += [f'<rect x="{x}" y="{y}" width="1200" height="964" fill="white" stroke="#c4d1db"/>',
                f'<text x="{x + 16}" y="{y + 28}" font-family="sans-serif" font-size="20" font-weight="bold" fill="#17364b">{escape(heading)}</text>',
                f'<text x="{x + 16}" y="{y + 51}" font-family="sans-serif" font-size="14" fill="#17364b">{escape(caption)}</text>',
                f'<image x="{x + 1}" y="{y + 64}" width="1198" height="899" href="data:image/svg+xml;base64,{data}"/>']
        html.append(f'<article><h2>{escape(heading)}</h2><p>{escape(caption)}</p><img alt="{escape(heading)}" src="data:image/svg+xml;base64,{data}"></article>')
    footer = "Exact original coordinates, mesh/source hashes, units and per-cell measured values are embedded in each SVG. Shape measurements are not a simulation certificate."
    svg += [f'<text x="24" y="2110" font-family="sans-serif" font-size="14" fill="#17364b">{escape(footer)}</text>', '</svg>']
    html += [f'</section><footer>{escape(footer)}</footer></body></html>']
    budget.write(root / "gallery.svg", "\n".join(svg).encode("utf-8"))
    budget.write(root / "gallery.html", "\n".join(html).encode("utf-8"))


def render(args):
    from paacs_mesh.mesh_inspection import open_mesh
    from paacs_mesh.workflow import _rename_new_directory

    source = args.benchmark.resolve(strict=True)
    if not source.is_dir():
        raise ValueError("benchmark must be a result directory")
    output = args.output.absolute()
    if os.path.lexists(output):
        raise FileExistsError("Output directory must be new")
    parent = output.parent.resolve(strict=True)
    output = parent / output.name
    budget = Budget(args)
    cairo = None
    png_status = {"requested": args.png, "status": "disabled"}
    if args.png != "none":
        try:
            cairo = importlib.import_module("cairosvg")
            png_status.update(status="available", renderer="CairoSVG", version=cairo.__version__)
        except (ImportError, OSError) as error:
            if args.png == "cairosvg":
                raise ValueError("CairoSVG and its native Cairo dependency are required for --png cairosvg") from error
            png_status.update(status="unavailable", reason="optional_cairosvg_not_available")
    with tempfile.TemporaryDirectory(prefix=".paacs-refinement-", dir=parent) as temporary:
        staging = Path(temporary).resolve() / "result"
        if staging.parent.parent != parent:
            raise ValueError("Staging directory escaped its output parent")
        staging.mkdir()
        (staging / "bundles").mkdir()
        inputs, summaries = [], {}
        common, source_hash, units = None, None, None
        limits = dict(max_bytes=args.max_bytes, max_work=args.max_work, max_seconds=args.max_seconds,
                      max_cells=args.max_cells, max_nodes=args.max_nodes)
        for case in CASES:
            case_dir = (source / case).resolve(strict=True)
            if not case_dir.is_relative_to(source):
                raise ValueError("Case directory escaped its benchmark input")
            report_path = case_dir / "native-report.json"
            report_data = budget.read(report_path, metadata=True)
            report = json.loads(report_data, object_pairs_hook=_pairs,
                                parse_constant=_constant, parse_float=_finite_float)
            if (not isinstance(report, dict) or report.get("mesh_published") is not True or report.get("reason") != "none"
                    or any(not isinstance(report.get(name), dict) or report[name].get("valid") is not True
                           or report[name].get("complete") is not True for name in ("audit", "quality"))):
                raise ValueError(f"{case} lacks a complete reported-success native result")
            identity, length_units = report.get("source_sha256"), report.get("units")
            if not isinstance(identity, str) or len(identity) != 64 or any(c not in "0123456789abcdef" for c in identity):
                raise ValueError("Every case requires a known lowercase source SHA-256")
            if not isinstance(length_units, str) or not length_units:
                raise ValueError("Every case requires reported length units")
            if source_hash is None:
                source_hash, units = identity, length_units
            if (source_hash, units) != (identity, length_units):
                raise ValueError("Selected cases have different CAD sources or units")
            bundle = staging / "bundles" / case
            bundle.mkdir()
            archive_hash = budget.copy_archive(case_dir / "mesh.paacsmesh", bundle / "mesh.paacsmesh")
            measurement_path = case_dir / "measurements.json"
            if measurement_path.exists():
                expected = budget.json(measurement_path).get("archive_sha256")
                if expected is not None and expected != archive_hash:
                    raise ValueError(f"{case} archive differs from its recorded measurements")
            budget.write(bundle / "report.json", report_data)
            manifest = {"schema": "paacs-mesh-inspection-bundle-v1", "mesh_sha256": archive_hash,
                        "report_sha256": hashlib.sha256(report_data).hexdigest(), "source_sha256": identity,
                        "units": units, "validation_basis": SAVED_BASIS, "native_reaudit_performed": False}
            budget.write(bundle / "inspection.json", budget.encoded_json(manifest))
            with open_mesh(bundle, **limits) as inspected:
                summary = inspected.summary(spatial_bins=1)
                summaries[case] = summary
                box = summary["bounds"]
                common = box if common is None else [[min(common[0][i], box[0][i]) for i in range(3)],
                                                     [max(common[1][i], box[1][i]) for i in range(3)]]
                inputs.append({"case": case, **inspected.identity(), "report_sha256": manifest["report_sha256"],
                               "validation": inspected.validation, "native_reaudit_performed": False})
            budget.check()
        frame = [args.view_bounds[:3], args.view_bounds[3:]] if args.view_bounds else _padded(common, .05)
        images = []
        for index, case in enumerate(CASES):
            with open_mesh(staging / "bundles" / case, **limits) as inspected:
                options = dict(metric="mean_ratio", highlight_below=.1, view_bounds=frame)
                filename = ("overview.svg", "wall-0p2-cut.svg", "wall-0p08-cut.svg", "worst-cell-closeup.svg")[index]
                if index in (1, 2):
                    options.update(cut_axis="y", cut_value=0)
                if index == 3:
                    worst = inspected.query(sort_by="mean_ratio", limit=1, associations=False)["cells"]
                    if not worst:
                        raise ValueError("The fine-wall case has no supported Tet4 for a close-up")
                    options["highlight_cell_ids"] = [worst[0]["cell_id"]]
                    options["view_bounds"] = ([args.closeup_bounds[:3], args.closeup_bounds[3:]] if args.closeup_bounds
                                               else _padded(worst[0]["bounds"], .25))
                destination = staging / filename
                receipt = inspected.visualize(destination, max_output_bytes=min(32 * MIB, args.max_output_bytes - budget.output_bytes), **options)
                budget.generated(destination)
                receipt["path"] = filename
                receipt["input_case"] = case
                receipt["native_reaudit_performed"] = False
                receipt["validation_basis"] = SAVED_BASIS
                if index == 0:
                    heading = "Overview · global size " + str(summaries[case]["sizing"]["requested"]["global_size"]) + " " + units
                elif index == 3:
                    heading = "Worst Tet4 · original cell " + str(receipt["highlighted_cell_ids"][0])
                else:
                    targets = sorted({control["size"] for control in summaries[case]["sizing"]["requested"]["face_sizes"]})
                    heading = "Whole-cell cut · wall target " + ", ".join(str(value) for value in targets) + " " + units
                receipt["heading"] = heading
                receipt["caption"] = f'{summaries[case]["cell_counts"]["tet4"]:,} Tet4 · {receipt["highlight_count"]:,} measured highlights · exact original mesh'
                images.append(receipt)
            budget.check()
        _gallery(images, staging, budget)
        if cairo is not None:
            png_status["files"] = []
            for filename in [row["path"] for row in images] + ["gallery.svg"]:
                budget.check()
                pixels = 2464 * 2150 if filename == "gallery.svg" else 1200 * 900
                if 16 * pixels > args.max_bytes:
                    raise ValueError("Optional PNG raster workspace exceeds max-bytes")
                destination = staging / Path(filename).with_suffix(".png")
                with destination.open("xb") as stream:
                    class BoundedWriter:
                        def write(self, data):
                            budget.check()
                            if len(data) > args.max_output_bytes - budget.output_bytes:
                                raise ValueError("PNG exceeds the remaining output byte allowance")
                            count = stream.write(data)
                            budget.output_bytes += count
                            return count
                    cairo.svg2png(url=str(staging / filename), write_to=BoundedWriter())
                budget.check()
                png_status["files"].append({"file": destination.name, "sha256": _sha256(destination),
                                            "bytes": destination.stat().st_size,
                                            "source_svg_sha256": _sha256(staging / filename)})
            png_status["status"] = "rendered"
        result = {"schema": "paacs-refinement-render-v1", "native_reaudit_performed": False,
                  "source_sha256": source_hash, "units": units, "inputs": inputs, "images": images,
                  "gallery": {name: {"file": name, "sha256": _sha256(staging / name), "bytes": (staging / name).stat().st_size}
                              for name in ("gallery.svg", "gallery.html")}, "png": png_status,
                  "limits": {key: getattr(args, key) for key in ("max_input_bytes", "max_output_bytes", "max_bytes", "max_cells", "max_nodes", "max_work", "max_seconds", "total_seconds")},
                  "scope": "Native reports copied verbatim from retained attempts; inspection and rendering are new measurements, not native recertification",
                  "input_bytes_read": budget.input_bytes}
        budget.write(staging / "render-manifest.json", budget.encoded_json(result))
        budget.check()
        _rename_new_directory(staging, output)
    return {"output": str(output), "manifest": "render-manifest.json", "svg_images": len(images), "png": png_status}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("benchmark", type=Path, help="result directory from refinement_diagnostics.py")
    parser.add_argument("--output", type=Path, required=True, help="new output directory; never replaces an existing path")
    parser.add_argument("--png", choices=("auto", "none", "cairosvg"), default="auto")
    parser.add_argument("--view-bounds", type=float, nargs=6, help="fixed common x0 y0 z0 x1 y1 z1 world frame")
    parser.add_argument("--closeup-bounds", type=float, nargs=6, help="explicit worst-cell world frame")
    for name, default in (("max-input-bytes", 256 * MIB), ("max-output-bytes", 64 * MIB),
                          ("max-bytes", 256 * MIB), ("max-cells", 1_000_000),
                          ("max-nodes", 1_000_000), ("max-work", 20_000_000)):
        parser.add_argument("--" + name, type=int, default=default)
    parser.add_argument("--max-seconds", type=float, default=60, help="cooperative allowance per inspection operation")
    parser.add_argument("--total-seconds", type=float, default=300, help="cooperative allowance for the entire script")
    args = parser.parse_args(argv)
    try:
        for key in ("max_input_bytes", "max_output_bytes", "max_bytes", "max_cells", "max_nodes", "max_work", "max_seconds", "total_seconds"):
            value = getattr(args, key)
            if value <= 0 or not math.isfinite(value):
                raise ValueError(f"{key} must be positive and finite")
        print(json.dumps(render(args), allow_nan=False))
        return 0
    except (ValueError, TypeError, OverflowError, OSError, ImportError, KeyError, AttributeError) as error:
        parser.exit(2, f"Rendering refused: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
