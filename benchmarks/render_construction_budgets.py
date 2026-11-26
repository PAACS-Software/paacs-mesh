"""Render retained construction-budget runs without importing or meshing CAD.

Requires the development mesh-inspection API. SVG and self-contained HTML use
only the standard library; CairoSVG is optional. Example:
  python benchmarks/render_construction_budgets.py --case Legacy RUN_A --case "4x caps" RUN_B --output NEW_DIRECTORY --png none
Every input directory must contain case.json, native-report.json and
mesh.paacsmesh from construction_budget_trials.py. Reports remain saved
evidence, hash-paired with their mesh, without a new native audit.
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
import struct
import subprocess
import tempfile
import time

try:
    from .render_refinement_diagnostics import Budget, MIB, SAVED_BASIS, _constant, _finite_float, _padded, _pairs, _sha256
except ImportError:
    from render_refinement_diagnostics import Budget, MIB, SAVED_BASIS, _constant, _finite_float, _padded, _pairs, _sha256


def _hash(value, name):
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(f"{name} must be a recorded lowercase SHA-256")
    return value


def _number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


def _display(value):
    return "not recorded" if value is None else f"{value:.3g}"


def _global_budget(case):
    values = case.get("global_max_work"), case.get("global_max_working_bytes")
    if any(value is not None and (isinstance(value, bool) or not isinstance(value, int) or value <= 0)
           for value in values):
        raise ValueError("Recorded global allowances must be positive integers or unavailable")
    return values


def _native_report_pairing(case, report_bytes):
    fields = ("native_report_file", "native_report_sha256")
    if not any(name in case for name in fields):
        return {"basis": "legacy_unavailable", "verified": False}
    if case.get("native_report_file") != "native-report.json":
        raise ValueError("Saved native report filename must be native-report.json")
    expected = _hash(case.get("native_report_sha256"), "native_report_sha256")
    if hashlib.sha256(report_bytes).hexdigest() != expected:
        raise ValueError("Native report differs from the recorded file SHA-256")
    return {"basis": "saved_file_sha256", "verified": True, "sha256": expected}


def _resource_receipt(case, report):
    requested_work, requested_memory = _global_budget(case)
    resources = report.get("resources")
    requested_policy = case.get("options", {}).get("resource_policy")
    if resources is None:
        if requested_policy is not None:
            raise ValueError("Case resource policy lacks a saved native resource receipt")
        return {"basis": "legacy_requested_only", "policy": None,
                "requested_work": requested_work, "requested_working_bytes": requested_memory,
                "effective_work": None, "effective_working_bytes": None,
                "logical_work_cap_disabled": None}
    if not isinstance(resources, dict) or resources.get("policy") not in {"normal", "bounded"}:
        raise ValueError("Saved native resource policy is missing or invalid")
    policy = resources["policy"]
    if requested_policy is not None and requested_policy != policy:
        raise ValueError("Case and native report differ on resource policy")
    required = {"requested_max_work", "effective_max_work", "max_working_bytes", "logical_work_cap_disabled"}
    if not required.issubset(resources):
        raise ValueError("Saved native resource limits are incomplete")
    native_requested = resources["requested_max_work"]
    if native_requested is not None and (isinstance(native_requested, bool)
            or not isinstance(native_requested, int) or native_requested <= 0):
        raise ValueError("Saved native requested work limit is invalid")
    if native_requested != requested_work:
        raise ValueError("Case and native report differ on requested global work")
    effective_work, effective_memory = resources["effective_max_work"], resources["max_working_bytes"]
    if any(isinstance(value, bool) or not isinstance(value, int) or value <= 0
           for value in (effective_work, effective_memory)):
        raise ValueError("Saved native effective resource limits must be positive integers")
    disabled = resources["logical_work_cap_disabled"]
    if not isinstance(disabled, bool) or disabled != (policy == "normal" and native_requested is None):
        raise ValueError("Saved native work-cap state contradicts its policy and request")
    if requested_work is not None and effective_work != requested_work:
        raise ValueError("Saved native effective work differs from the requested global ceiling")
    if requested_work is None and policy == "bounded" and effective_work != 50_000_000:
        raise ValueError("Saved native bounded default work differs from its declared policy")
    if requested_work is None and policy == "normal" and effective_work < 2**60:
        raise ValueError("Saved native disabled work cap has an inconsistent practical ceiling")
    if requested_memory is not None and effective_memory != requested_memory:
        raise ValueError("Case and native report differ on the global working allowance")
    return {"basis": "saved_native_resources", "policy": policy,
            "requested_work": requested_work, "requested_working_bytes": requested_memory,
            "effective_work": effective_work, "effective_working_bytes": effective_memory,
            "logical_work_cap_disabled": disabled}


def _budget_identity(resources):
    if resources["basis"] == "legacy_requested_only":
        return (resources["basis"], resources["requested_work"], resources["requested_working_bytes"])
    return (resources["basis"], resources["policy"], resources["effective_work"],
            resources["effective_working_bytes"], resources["logical_work_cap_disabled"])


def _sizing_targets(case, report):
    if "wall_size_mm" not in case:
        raise ValueError("Saved case lacks its wall/global-only sizing declaration")
    wall_value = case.get("wall_size_mm")
    wall = None if wall_value is None else _number(wall_value, "wall_size_mm")
    global_size = _number(case.get("global_size_mm"), "global_size_mm")
    reported_size = _number(report.get("requested_size"), "requested_size")
    if (wall is not None and wall <= 0) or global_size <= 0 or reported_size != global_size:
        raise ValueError("Case and native report must agree on the positive global size")
    sizing = report.get("sizing", {})
    if not isinstance(sizing, dict) or not isinstance(sizing.get("face_sizes", []), list):
        raise ValueError("Saved native face-sizing controls are invalid")
    controls = sizing.get("face_sizes", [])
    if wall is None:
        if "face_sizes" not in sizing:
            raise ValueError("A global-only case needs its saved native face-sizing list")
        if controls:
            raise ValueError("A global-only case must not contain native face-sizing controls")
    elif not any(isinstance(control, dict) and control.get("size") == wall for control in controls):
        raise ValueError("The displayed wall target differs from the saved native face-sizing request")
    return wall, global_size


def _caption(stats):
    count = stats["tet4"]
    poor = stats["below_threshold"]
    resources = stats["effective_resources"]
    if resources["basis"] == "legacy_requested_only":
        work = f'effective unavailable (requested {_display(resources["requested_work"])})'
        memory = resources["requested_working_bytes"]
        memory = f'effective unavailable (requested {"unavailable" if memory is None else f"{memory / MIB:g} MiB"})'
    else:
        work = ("no practical work ceiling" if resources["logical_work_cap_disabled"]
                else _display(resources["effective_work"])) + f' ({resources["policy"]})'
        memory = f'{resources["effective_working_bytes"] / MIB:g} MiB modeled'
    return [
        f'{count:,} Tet4; {stats["nodes"]:,} nodes; {stats["native_seconds"]:.2f} s saved native time; one attempt',
        f'Mean ratio: min {_display(stats["mean_ratio"]["minimum"])}; p05 {_display(stats["mean_ratio"]["p05"])}; median {_display(stats["mean_ratio"]["median"])}',
        f'Below {stats["threshold"]:g}: {poor:,}/{count:,} ({100 * poor / count if count else 0:.2f}%); dihedral p05 {_display(stats["dihedral"]["p05"])} deg; SJ p05 {_display(stats["jacobian"]["p05"])}',
        f'Global work {work}; working allowance {memory}',
        f'Native SHA {stats["native_sha256"][:12]}; saved report file SHA {"verified" if stats.get("native_report_pairing", {}).get("verified") else "unavailable (legacy)"}',
    ]


def _gallery(rows, directory, budget, args, warnings):
    sizes = sum((directory / row["path"]).stat().st_size for row in rows)
    if 16 * sizes + 2 * MIB > args.max_bytes:
        raise ValueError("Self-contained gallery encoding exceeds max-bytes")
    columns = min(args.columns, len(rows))
    lines = (len(rows) + columns - 1) // columns
    width, height = 32 + columns * 1216, 210 + lines * 1080
    title = "PAACS Mesh / construction budget evidence"
    note = "Actual archived Tet4 meshes; same camera and whole-element cut. Saved native results, not independently re-audited."
    metadata = {"schema": "paacs-construction-gallery-v1", "native_reaudit_performed": False,
                "warnings": warnings, "images": [{key: row[key] for key in ("path", "mesh_sha256", "svg_sha256", "input_case")} for row in rows]}
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           f'<title>{escape(title)}</title><metadata>{escape(json.dumps(metadata, sort_keys=True))}</metadata>',
           '<rect width="100%" height="100%" fill="#f3f6fa"/>',
           f'<text x="24" y="48" font-family="sans-serif" font-size="30" font-weight="bold" fill="#17364b">{escape(title)}</text>',
           f'<text x="24" y="77" font-family="sans-serif" font-size="17" fill="#17364b">{escape(note)}</text>']
    html = ['<!doctype html><html><head><meta charset="utf-8"><title>' + escape(title) + '</title>',
            '<style>body{margin:0;padding:24px;background:#f3f6fa;color:#17364b;font-family:Arial,sans-serif}h1{margin:0 0 8px}p{margin:5px 0}.warning{font-weight:bold;color:#8b3424}.grid{display:grid;grid-template-columns:repeat(' + str(columns) + ',minmax(0,1fr));gap:16px;margin-top:20px}article{background:white;border:1px solid #c4d1db}h2{font-size:18px;margin:12px 14px 8px}article p{font-size:13px;margin:0 14px 4px}img{display:block;width:100%;height:auto}footer{margin-top:18px}@media(max-width:1000px){.grid{grid-template-columns:1fr}}</style></head><body>',
            f'<h1>{escape(title)}</h1><p>{escape(note)}</p>']
    for i, warning in enumerate(warnings):
        svg.append(f'<text x="24" y="{102 + 20 * i}" font-family="sans-serif" font-size="15" fill="#8b3424">{escape(warning)}</text>')
        html.append(f'<p class="warning">{escape(warning)}</p>')
    html.append('<section class="grid">')
    for index, row in enumerate(rows):
        x, y = 24 + (index % columns) * 1216, 150 + (index // columns) * 1080
        payload = base64.b64encode((directory / row["path"]).read_bytes()).decode("ascii")
        heading = row["heading"]
        captions = _caption(row["stats"])
        svg += [f'<rect x="{x}" y="{y}" width="1200" height="1064" fill="white" stroke="#c4d1db"/>',
                f'<text x="{x + 14}" y="{y + 29}" font-family="sans-serif" font-size="20" font-weight="bold" fill="#17364b">{escape(heading)}</text>']
        for line, caption in enumerate(captions):
            svg.append(f'<text x="{x + 14}" y="{y + 53 + 20 * line}" font-family="sans-serif" font-size="14" fill="#17364b">{escape(caption)}</text>')
        svg.append(f'<image x="{x + 1}" y="{y + 145}" width="1198" height="899" href="data:image/svg+xml;base64,{payload}"/>')
        html.append('<article><h2>' + escape(heading) + '</h2>' + ''.join('<p>' + escape(c) + '</p>' for c in captions)
                    + '<img alt="' + escape(heading) + '" src="data:image/svg+xml;base64,' + payload + '"></article>')
    footer = "Larger caps do not guarantee better quality. Auto/Thorough scheduling stays opt-in. Shape measurements do not establish simulation convergence."
    svg += [f'<text x="24" y="{height - 20}" font-family="sans-serif" font-size="14" fill="#17364b">{escape(footer)}</text>', '</svg>']
    html += ['</section><footer>' + escape(footer) + '</footer></body></html>']
    budget.write(directory / "gallery.svg", '\n'.join(svg).encode("utf-8"))
    budget.write(directory / "gallery.html", '\n'.join(html).encode("utf-8"))
    return width, height


def _png(directory, rows, gallery_size, budget, args):
    status = {"requested": args.png, "status": "disabled"}
    if args.png == "none":
        return status
    if args.png == "browser":
        return _browser_png(directory, rows, gallery_size, budget, args)
    try:
        cairo = importlib.import_module("cairosvg")
    except (ImportError, OSError) as error:
        if args.png == "cairosvg":
            raise ValueError("CairoSVG and native Cairo are required for --png cairosvg") from error
        return {"requested": args.png, "status": "unavailable", "reason": "optional_cairosvg_not_available"}
    status.update(status="rendered", renderer="CairoSVG", version=cairo.__version__, files=[])
    for filename in [row["path"] for row in rows] + ["gallery.svg"]:
        budget.check()
        pixels = math.prod(gallery_size) if filename == "gallery.svg" else 1200 * 900
        if 16 * pixels > args.max_bytes:
            raise ValueError("Optional PNG workspace exceeds max-bytes")
        destination = directory / Path(filename).with_suffix(".png")
        with destination.open("xb") as stream:
            class Writer:
                def write(self, data):
                    budget.check()
                    if len(data) > args.max_output_bytes - budget.output_bytes:
                        raise ValueError("PNG exceeds max-output-bytes")
                    count = stream.write(data)
                    budget.output_bytes += count
                    return count
            cairo.svg2png(url=str(directory / filename), write_to=Writer())
        budget.check()
        status["files"].append({"file": destination.name, "sha256": _sha256(destination),
                                "bytes": destination.stat().st_size})
    return status


def _browser_png(directory, rows, gallery_size, budget, args):
    """Use an explicitly supplied existing Chromium browser, without user data."""
    if args.browser_executable is None or not args.browser_executable.is_file():
        raise ValueError("--png browser requires an existing --browser-executable")
    browser = args.browser_executable.resolve(strict=True)
    status = {"requested": "browser", "status": "rendered", "renderer": str(browser),
              "profile": "fresh workspace-only temporary profile; no user profile", "files": []}
    # Browser helper processes can briefly retain Crashpad files after the
    # screenshot process exits. Keep this private scratch OUTSIDE deliverables;
    # a locked cleanup must not turn a complete image into partial publication.
    profile_parent = directory.parent.parent.resolve()
    with tempfile.TemporaryDirectory(prefix=".paacs-browser-render-", dir=profile_parent,
                                     ignore_cleanup_errors=True) as temporary:
        profile = Path(temporary).resolve()
        if profile.parent != profile_parent or not profile.name.startswith(".paacs-browser-render-"):
            raise ValueError("Browser profile escaped the staging directory")
        for filename in [row["path"] for row in rows] + ["gallery.svg"]:
            budget.check()
            # A lingering browser helper must not turn the next invocation
            # into an asynchronous request to an already running instance.
            item_profile = profile / Path(filename).stem
            width, height = gallery_size if filename == "gallery.svg" else (1200, 900)
            pixels = width * height
            if 32 * pixels + 8 * (directory / filename).stat().st_size > args.max_bytes:
                raise ValueError("Browser raster payload allowance exceeds max-bytes (not a process RSS cap)")
            if 5 * pixels + 65536 > args.max_output_bytes - budget.output_bytes:
                raise ValueError("Worst-case PNG output exceeds the remaining byte allowance")
            destination = directory / Path(filename).with_suffix(".png")
            command = [str(browser), "--headless=new", "--disable-gpu", "--hide-scrollbars", "--no-first-run",
                       "--no-default-browser-check", "--disable-extensions", "--disable-background-networking",
                       "--disable-component-update", "--disable-sync", "--disable-crash-reporter",
                       "--force-device-scale-factor=1", f"--user-data-dir={item_profile}",
                       f"--window-size={width},{height}", f"--screenshot={destination}",
                       (directory / filename).resolve().as_uri()]
            timeout = min(60, args.total_seconds - (time.perf_counter() - budget.started))
            if timeout <= 0:
                raise ValueError("No rendering deadline remains for the browser")
            try:
                result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                        timeout=timeout, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            except subprocess.TimeoutExpired as error:
                raise ValueError("Browser rasterization exceeded its bounded subprocess deadline") from error
            if result.returncode or not destination.is_file():
                raise ValueError(f"Browser rasterization failed with exit code {result.returncode}")
            with destination.open("rb") as stream:
                header = stream.read(24)
            if (len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n"
                    or struct.unpack(">II", header[16:24]) != (width, height)):
                raise ValueError("Browser did not produce the requested PNG dimensions")
            budget.generated(destination)
            status["files"].append({"file": destination.name, "sha256": _sha256(destination),
                                    "bytes": destination.stat().st_size, "width": width, "height": height,
                                    "source_svg_sha256": _sha256(directory / filename)})
    status["profile_cleanup"] = "deferred" if profile.exists() else "complete"
    if profile.exists():
        status["retained_private_profile"] = str(profile)
    return status


def render(args):
    if args.development_native is not None:
        try:
            from .construction_budget_trials import load_development
        except ImportError:
            from construction_budget_trials import load_development
        load_development(args)
    import paacs_mesh as pm
    from paacs_mesh.workflow import _rename_new_directory

    output = args.output.absolute()
    if os.path.lexists(output):
        raise FileExistsError("Output directory must be new")
    parent = output.parent.resolve(strict=True)
    output = parent / output.name
    budget = Budget(args)
    limits = dict(max_bytes=args.max_bytes, max_cells=args.max_cells, max_nodes=args.max_nodes,
                  max_work=args.max_work, max_seconds=args.max_seconds)
    with tempfile.TemporaryDirectory(prefix=".paacs-construction-render-", dir=parent) as temporary:
        staging = Path(temporary).resolve() / "result"
        if staging.parent.parent != parent:
            raise ValueError("Staging escaped its output parent")
        staging.mkdir()
        (staging / "bundles").mkdir()
        records, shared, hashes, globals_seen = [], None, set(), set()
        unknown_resources = False
        unknown_report_pairing = False
        common = None
        for index, (label, raw_directory) in enumerate(args.case):
            directory = Path(raw_directory).resolve(strict=True)
            if not directory.is_dir():
                raise ValueError("Every case input must be a result directory")
            case = budget.json(directory / "case.json")
            report_bytes = budget.read(directory / "native-report.json", metadata=True)
            report_pairing = _native_report_pairing(case, report_bytes)
            report = json.loads(report_bytes, object_pairs_hook=_pairs, parse_constant=_constant, parse_float=_finite_float)
            if (case.get("schema") != "paacs-construction-budget-case-v1" or case.get("status") != "succeeded"
                    or case.get("native_status") != "succeeded" or case.get("retry_count") != 0
                    or case.get("workers") != 1 or case.get("validation_workers") != 1):
                raise ValueError("Only saved successful single-worker, zero-retry construction runs are supported")
            if (not isinstance(report, dict) or report.get("mesh_published") is not True or report.get("reason") != "none"
                    or any(not isinstance(report.get(name), dict) or report[name].get("complete") is not True
                           or report[name].get("valid") is not True for name in ("audit", "quality"))):
                raise ValueError("Case lacks complete reported native audit and quality success")
            source = _hash(report.get("source_sha256"), "source_sha256")
            if case.get("source_sha256") != source or report.get("units") != "mm" or not report.get("kernel_version"):
                raise ValueError("Saved case/report source, millimetre units and kernel identity must agree")
            wall, global_size = _sizing_targets(case, report)
            native_seconds = _number(case.get("native_seconds"), "native_seconds")
            seconds = report.get("seconds", {})
            if not isinstance(seconds, dict) or seconds.get("meshing_total") != native_seconds:
                raise ValueError("Case and native report differ on the saved meshing time")
            options, construction = case.get("options"), report.get("construction")
            if not isinstance(options, dict) or not isinstance(construction, dict):
                raise ValueError("Saved construction settings are missing")
            for name, value in options.items():
                if name == "resource_policy":
                    continue  # Resource intent has its own native report block.
                recorded = "effort" if name == "construction_effort" else name
                if construction.get(recorded) != value:
                    raise ValueError("Case and native report differ on construction settings")
            resources = _resource_receipt(case, report)
            identity = (source, report["units"], report["kernel_version"], global_size)
            if shared is not None and identity != shared:
                raise ValueError("Cases differ in CAD source, units, kernel or requested global size")
            shared = identity
            native_hash = _hash(case.get("native_extension_sha256"), "native_extension_sha256")
            hashes.add(native_hash)
            globals_seen.add(_budget_identity(resources))
            unknown_resources |= resources["basis"] == "legacy_requested_only"
            unknown_report_pairing |= not report_pairing["verified"]
            case_id = f"case-{index + 1:02}"
            bundle = staging / "bundles" / case_id
            bundle.mkdir()
            archive_hash = budget.copy_archive(directory / "mesh.paacsmesh", bundle / "mesh.paacsmesh")
            if archive_hash != _hash(case.get("mesh_sha256"), "mesh_sha256"):
                raise ValueError("Archive differs from the recorded case mesh SHA-256")
            budget.write(bundle / "report.json", report_bytes)
            budget.write(bundle / "case.json", budget.encoded_json(case))
            manifest = {"schema": "paacs-mesh-inspection-bundle-v1", "mesh_sha256": archive_hash,
                        "report_sha256": hashlib.sha256(report_bytes).hexdigest(), "source_sha256": source,
                        "units": report["units"], "validation_basis": SAVED_BASIS, "native_reaudit_performed": False}
            budget.write(bundle / "inspection.json", budget.encoded_json(manifest))
            with pm.inspect_mesh(bundle, **limits) as inspected:
                summary = inspected.summary(spatial_bins=1)
                count = summary["cell_counts"]["tet4"]
                if case.get("nodes") != inspected.node_count or case.get("tetrahedra") != count:
                    raise ValueError("Recorded counts differ from the actual archive")
                for field, actual in (("nodes", inspected.node_count), ("tetrahedra", count)):
                    if field in report and (isinstance(report[field], bool)
                            or not isinstance(report[field], int) or report[field] != actual):
                        raise ValueError("Saved native report counts differ from the actual archive")
                poor = inspected.query(quality_below=args.highlight_below, limit=1, associations=False)["total"]
                extent = summary["bounds"]
                if extent is None:
                    raise ValueError("An empty archive cannot be compared")
                common = extent if common is None else [[min(common[0][a], extent[0][a]) for a in range(3)],
                                                        [max(common[1][a], extent[1][a]) for a in range(3)]]
                stats = {"tet4": count, "nodes": inspected.node_count, "wall_size_mm": wall,
                         "global_size_mm": global_size,
                         "native_seconds": native_seconds,
                         "mean_ratio": summary["quality"]["mean_ratio"],
                         "dihedral": summary["quality"]["minimum_dihedral_deg"],
                         "jacobian": summary["quality"]["scaled_jacobian"], "threshold": args.highlight_below,
                         "below_threshold": poor, "native_sha256": native_hash,
                         "global_budgets": _global_budget(case), "effective_resources": resources,
                         "native_report_pairing": report_pairing,
                         "profile": case.get("profile"),
                         "reported_peak_process_bytes_after_mesh": case.get("peak_process_bytes_after_mesh")}
                records.append({"input_case": case_id, "label": label, "input_directory": str(directory),
                                "heading": (f"Global size {global_size:g} mm" if wall is None else f"Wall {wall:g} mm") + f" / {label}", "stats": stats,
                                "mesh_sha256": archive_hash, "source_sha256": source,
                                "validation": inspected.validation, "native_reaudit_performed": False})
            budget.check()
        if len(hashes) > 1 and not args.allow_mixed_native:
            raise ValueError("Recorded native binaries differ; use --allow-mixed-native to show this explicitly")
        if len(globals_seen) > 1 and not args.allow_mixed_budgets:
            raise ValueError("Global work/memory allowances differ; use --allow-mixed-budgets to show this explicitly")
        warnings = []
        if len(hashes) > 1:
            warnings.append("Different recorded native binaries: retained-run evidence, not an isolated algorithm comparison.")
        if len(globals_seen) > 1:
            warnings.append("Resource policies or global limits differ; these are not equal-work performance comparisons.")
        if unknown_resources or unknown_report_pairing:
            unavailable = [name for name, missing in (("effective resource limits", unknown_resources),
                           ("saved report file hash", unknown_report_pairing)) if missing]
            warnings.append("Legacy evidence unavailable: " + " and ".join(unavailable) + "; corresponding claims remain unverified.")
        frame = [args.view_bounds[:3], args.view_bounds[3:]] if args.view_bounds else _padded(common, .05)
        axis = "xyz".index(args.cut_axis)
        cut = args.cut_value if args.cut_value is not None else common[0][axis] / 2 + common[1][axis] / 2
        images = []
        for record in records:
            with pm.inspect_mesh(staging / "bundles" / record["input_case"], **limits) as inspected:
                filename = record["input_case"] + ".svg"
                receipt = inspected.visualize(staging / filename, metric="mean_ratio", highlight_below=args.highlight_below,
                    cut_axis=args.cut_axis, cut_value=cut, view_bounds=frame, azimuth=args.azimuth, elevation=args.elevation,
                    width=1200, height=900, max_triangles=args.max_triangles, max_highlights=args.max_highlights,
                    max_output_bytes=min(32 * MIB, args.max_output_bytes - budget.output_bytes))
                budget.generated(staging / filename)
                images.append({**record, **receipt, "path": filename, "validation_basis": SAVED_BASIS,
                               "native_reaudit_performed": False})
            budget.check()
        gallery_size = _gallery(images, staging, budget, args, warnings)
        png = _png(staging, images, gallery_size, budget, args)
        result = {"schema": "paacs-construction-render-v1", "source_sha256": shared[0], "units": shared[1],
                  "kernel_version": shared[2], "recorded_native_sha256": sorted(hashes), "warnings": warnings,
                  "native_reaudit_performed": False, "validation_basis": SAVED_BASIS, "images": images,
                  "same_camera": {"view_bounds": frame, "azimuth": args.azimuth, "elevation": args.elevation},
                  "cut": {"axis": args.cut_axis, "value": cut, "rule": "retain whole Tet4 by centroid; expose original element faces"},
                  "timing_basis": "saved one-attempt native mesh times; rendering is not included", "png": png,
                  "gallery": {name: {"sha256": _sha256(staging / name), "bytes": (staging / name).stat().st_size}
                              for name in ("gallery.svg", "gallery.html")}}
        budget.write(staging / "render-manifest.json", budget.encoded_json(result))
        budget.check()
        _rename_new_directory(staging, output)
    return {"output": str(output), "manifest": "render-manifest.json", "images": len(images), "png": png}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", nargs=2, action="append", metavar=("LABEL", "DIRECTORY"), required=True)
    parser.add_argument("--output", type=Path, required=True, help="new directory; never replaces existing output")
    parser.add_argument("--columns", type=int, default=2)
    parser.add_argument("--allow-mixed-native", action="store_true", help="explicitly label differing recorded native SHA-256 values")
    parser.add_argument("--allow-mixed-budgets", action="store_true", help="explicitly label differing global work/memory allowances")
    parser.add_argument("--png", choices=("auto", "none", "cairosvg", "browser"), default="auto")
    parser.add_argument("--browser-executable", type=Path, help="existing Chromium/Edge executable; only used with --png browser")
    parser.add_argument("--view-bounds", type=float, nargs=6)
    parser.add_argument("--cut-axis", choices=("x", "y", "z"), default="y")
    parser.add_argument("--cut-value", type=float)
    parser.add_argument("--azimuth", type=float, default=35)
    parser.add_argument("--elevation", type=float, default=25)
    parser.add_argument("--highlight-below", type=float, default=.1)
    parser.add_argument("--development-native", type=Path)
    parser.add_argument("--occt-sdk", type=Path)
    for name, default in (("max-input-bytes", 256 * MIB), ("max-output-bytes", 128 * MIB),
                          ("max-bytes", 512 * MIB), ("max-cells", 1_000_000), ("max-nodes", 1_000_000),
                          ("max-work", 20_000_000), ("max-triangles", 100000), ("max-highlights", 10000)):
        parser.add_argument("--" + name, type=int, default=default)
    parser.add_argument("--max-seconds", type=float, default=60)
    parser.add_argument("--total-seconds", type=float, default=300)
    args = parser.parse_args(argv)
    try:
        if not 1 <= len(args.case) <= 8 or not 1 <= args.columns <= 4:
            raise ValueError("Provide 1..8 cases and 1..4 columns")
        for label, _ in args.case:
            if not label.strip() or len(label) > 80 or "\0" in label:
                raise ValueError("Each panel label must be nonempty and at most 80 characters")
        for name in ("max_input_bytes", "max_output_bytes", "max_bytes", "max_cells", "max_nodes", "max_work",
                     "max_triangles", "max_highlights", "max_seconds", "total_seconds"):
            value = getattr(args, name)
            if value <= 0 or not math.isfinite(value):
                raise ValueError(f"{name} must be positive and finite")
        for value in (args.azimuth, args.elevation, args.highlight_below, *(args.view_bounds or []),
                      *([] if args.cut_value is None else [args.cut_value])):
            if not math.isfinite(value):
                raise ValueError("Camera, cut and threshold values must be finite")
        print(json.dumps(render(args), allow_nan=False))
        return 0
    except (ValueError, TypeError, OSError, ImportError, KeyError, AttributeError, OverflowError) as error:
        parser.exit(2, f"Rendering refused: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
