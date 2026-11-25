"""Inspect CAD and mesh quality, compare attempts, or mesh STEP with PAACS Mesh."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
from contextlib import contextmanager
import ctypes
import json
import os
from pathlib import Path
import sys

from . import (
    BoxSize, CadImportError, ExportError, FaceSize, MeshingError, SelectionError,
    WorkflowError, InspectionError, compare_meshes, open_mesh, publish_result, read_step, run_saved_recipe,
)
from ._cad_common import MAX_JSON_BYTES, json_text, load_json
from .diagnostics import capabilities, diagnose, format_diagnostics
from .mesh_archive import CELL_NAMES, atomic_bytes
from .mesh_inspection import METRICS
from .cad import (_CONSTRUCTION_LIMITS, _DEFAULT_CONSTRUCTION_EFFORT, _DEFAULT_RESOURCE_POLICY,
                  _RESOURCE_POLICIES)


_COMMANDS = ("capabilities", "inspect", "diagnose", "run", "inspect-mesh", "visualize", "compare")


def _flush_output():
    """Flush Python and C runtimes before changing the process stdout handle."""
    sys.stdout.flush()
    sys.stderr.flush()
    # The extension uses the host C runtime, not Python's text streams. Windows
    # extensions may use UCRT even when ctypes' default runtime is msvcrt.
    for library in (("ucrtbase", "msvcrt") if os.name == "nt" else (None,)):
        try:
            runtime = ctypes.CDLL(library)
            flush = runtime.fflush
            flush.argtypes, flush.restype = [ctypes.c_void_p], ctypes.c_int
            flush(None)
        except (OSError, AttributeError):
            continue


@contextmanager
def _native_output_to_stderr():
    """Keep native transfer messages outside a command's JSON stdout stream."""
    _flush_output()
    saved = os.dup(1)
    try:
        os.dup2(2, 1)
        yield
    finally:
        try:
            _flush_output()
        finally:
            os.dup2(saved, 1)
            os.close(saved)


def _json(value, *, stream=None):
    # Redirected Windows streams can use a legacy code page. JSON escapes keep
    # the machine stream ASCII-compatible while preserving Unicode on decode.
    payload = json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False)
    if len(payload) > MAX_JSON_BYTES:
        payload = json.dumps(value, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
    if len(payload) > MAX_JSON_BYTES:
        raise InspectionError("output_budget", "Machine output exceeds the 8 MiB JSON budget")
    print(payload, file=stream)


def _construction_cap(value):
    try:
        result = int(value)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError("construction caps must be positive integers") from error
    if not 1 <= result <= 2 ** 63 - 1:
        raise argparse.ArgumentTypeError("construction caps must be positive integers below 2^63")
    return result


def _fidelity_number(value, name):
    from .cad import _fidelity_options
    try:
        return _fidelity_options(**{name: float(value)})[name]
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def _curvature_points(value):
    from .cad import _curvature_value
    try:
        return _curvature_value(int(value))
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def _resource_cap(value):
    try:
        result = int(value)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError("resource limits must be positive integers") from error
    if not 1 <= result <= 2 ** 63 - 1:
        raise argparse.ArgumentTypeError("resource limits must be positive integers below 2^63")
    return result


def _construction_arguments(parser):
    parser.add_argument("--resource-policy", choices=_RESOURCE_POLICIES, default=_DEFAULT_RESOURCE_POLICY,
                        help="normal or bounded resource admission (default: %(default)s); explicit limits still apply")
    parser.add_argument("--construction-effort", choices=("legacy", "auto", "thorough"), default=_DEFAULT_CONSTRUCTION_EFFORT,
                        help="bounded construction allowance policy (default: %(default)s); validity checks stay enabled")
    for name in _CONSTRUCTION_LIMITS:
        parser.add_argument("--" + name.replace("_", "-"), type=_construction_cap,
                            help="optional positive construction cap; existing global budgets still apply")


def _legacy(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog="Discovery commands: " + ", ".join(_COMMANDS) + ". "
               "Use 'paacs-mesh COMMAND --help' for their options.",
    )
    parser.add_argument("input", type=Path, help="STEP file; transferred to mm")
    parser.add_argument("output", type=Path, help=".vtk or .paacsmesh destination")
    parser.add_argument("--size", type=float, required=True, help="global target size in mm")
    parser.add_argument("--recovery", choices=("standard", "advanced"), default="standard")
    parser.add_argument("--assembly", choices=("preserve_parts", "conforming"), default="preserve_parts")
    parser.add_argument("--face-size", nargs=2, action="append", default=[],
                        metavar=("TAG", "SIZE"), help="refine a CAD face; repeat for more faces")
    parser.add_argument("--size-box", nargs=7, type=float, action="append", default=[],
                        metavar=("X0", "Y0", "Z0", "X1", "Y1", "Z1", "SIZE"),
                        help="refine inside a box in mm; repeat for more regions")
    parser.add_argument("--min-size", type=float, help="minimum requested field size in mm")
    parser.add_argument("--max-chord-error", type=lambda value: _fidelity_number(value, "max_chord_error"),
                        help="positive CAD chord-error limit in mm; omitted uses the route default")
    parser.add_argument("--max-relative-volume-error", type=lambda value: _fidelity_number(value, "max_relative_volume_error"),
                        help="CAD relative-volume limit strictly between zero and one")
    parser.add_argument("--curvature-points-per-turn", type=_curvature_points,
                        help="integer 0..64; zero disables optional curvature sampling")
    parser.add_argument("--initial-topology-quality", action=argparse.BooleanOptionalAction, default=None,
                        help="severe initial-topology cleanup; default depends on resource policy, does not force the front")
    parser.add_argument("--native-source-fallback", action=argparse.BooleanOptionalAction, default=False,
                        help="experimental bounded native source after supported chart/wire refusal; preserves Native validity checks")
    parser.add_argument("--boundary-first", action=argparse.BooleanOptionalAction, default=False,
                        help="experimental live Tet4 order: recover boundary before interior insertion; does not force a constructor")
    parser.add_argument("--trimmed-ray-witnesses", action=argparse.BooleanOptionalAction, default=False,
                        help="opt in to strict serial CAD ray validation, with ordinary classifier fallback")
    parser.add_argument("--growth", type=float, default=1.3,
                        help="size-field growth, at least 1 (default: 1.3)")
    parser.add_argument("--native-sizing", action=argparse.BooleanOptionalAction, default=None,
                        help="force spatial sizing or retain legacy cap sizing; automatic with controls")
    parser.add_argument("--advancing-interior-front", action=argparse.BooleanOptionalAction, default=None,
                        help="select advancing initial Tet4 placement or the legacy constructor; automatic by default")
    parser.add_argument("--max-work", type=_resource_cap, help="explicit native logical work ceiling; overrides the resource policy")
    parser.add_argument("--max-working-bytes", type=_resource_cap, help="modeled native working memory allowance in bytes")
    parser.add_argument("--max-nodes", type=_resource_cap, help="maximum output nodes")
    parser.add_argument("--max-cells", type=_resource_cap, help="maximum output cells")
    parser.add_argument("--max-seconds", type=float, help="cooperative native cancellation deadline")
    parser.add_argument("--warning-quality", type=float,
                        help="preferred mean-ratio quality in [0,1]; 0 disables warnings")
    parser.add_argument("--minimum-quality", type=float,
                        help="strict minimum mean-ratio quality in [0,1]")
    _construction_arguments(parser)
    args = parser.parse_args(argv)
    try:
        with _native_output_to_stderr():
            model = read_step(args.input)
            face_sizes = [FaceSize(tag=int(tag), size=float(size)) for tag, size in args.face_size]
            volume_controls = [BoxSize(lower=tuple(box[:3]), upper=tuple(box[3:6]), size=box[6])
                               for box in args.size_box]
            mesh = model.mesh(args.size, recovery=args.recovery, assembly=args.assembly,
                              face_sizes=face_sizes, volume_controls=volume_controls,
                              min_size=args.min_size, growth=args.growth,
                              native_sizing=args.native_sizing, max_work=args.max_work,
                              max_working_bytes=args.max_working_bytes, max_nodes=args.max_nodes,
                              max_cells=args.max_cells, resource_policy=args.resource_policy,
                              advancing_interior_front=args.advancing_interior_front,
                              max_seconds=args.max_seconds, warning_quality=args.warning_quality,
                              minimum_quality=args.minimum_quality,
                              max_chord_error=args.max_chord_error,
                              max_relative_volume_error=args.max_relative_volume_error,
                              curvature_points_per_turn=args.curvature_points_per_turn,
                              initial_topology_quality=args.initial_topology_quality,
                              native_source_fallback=args.native_source_fallback,
                              boundary_first=args.boundary_first,
                              trimmed_ray_witnesses=args.trimmed_ray_witnesses,
                              construction_effort=args.construction_effort,
                              **{name: getattr(args, name) for name in _CONSTRUCTION_LIMITS})
            mesh.write(args.output)
        _json(mesh.report)
        return 0
    except (CadImportError, MeshingError, ExportError) as error:
        print(json.dumps({"status": "failed", "stage": error.stage,
                          "reason": error.reason, "report": error.report},
                         allow_nan=False), file=sys.stderr)
        return 1
    except (ValueError, TypeError, OverflowError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 2


class _CommandParser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError(message)


def _command_parser(command):
    parser = _CommandParser(prog=f"paacs-mesh {command}")
    parser.add_argument("--json", action="store_true", help="write one JSON document to stdout")
    if command in ("inspect-mesh", "visualize", "compare"):
        return _mesh_parser(parser, command)
    if command == "capabilities":
        parser.description = "Discover supported operations and their limitations."
    elif command == "inspect":
        parser.description = "Inspect one bounded page of source-bound CAD entities."
        parser.add_argument("input", type=Path, help="STEP file; transferred to mm")
        parser.add_argument("--kind", choices=("faces", "bodies", "edges", "vertices"), default="faces")
        parser.add_argument("--offset", type=int, default=0)
        parser.add_argument("--limit", type=int, default=100, help="entity page size, 1 through 1000")
        parser.add_argument("--topology", action="store_true", help="include bounded topology uses")
    elif command == "diagnose":
        parser.description = "Summarize an existing native or workflow report without rerunning its audits."
        parser.add_argument("report", type=Path, help="native mesh or saved workflow JSON report")
    else:
        parser.description = "Run a saved source-checked recipe and atomically create a new result bundle."
        parser.add_argument("recipe", type=Path, help="saved PAACS recipe JSON")
        parser.add_argument("--source", type=Path, required=True, help="explicit local STEP input")
        parser.add_argument("--output", type=Path, required=True, help="new result directory; never overwritten")
    return parser


def _mesh_parser(parser, command):
    parser.add_argument("input", type=Path, help="native archive or hash-paired inspection bundle")
    parser.add_argument("--max-bytes", type=int, default=256 * 1024 * 1024,
                        help="owned storage/scratch allowance per snapshot or operation; compare retains two snapshots")
    parser.add_argument("--max-cells", type=int, default=1_000_000)
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    parser.add_argument("--max-work", type=int, default=20_000_000)
    parser.add_argument("--max-seconds", type=float, default=60., help="cooperative deadline per inspection operation")
    if command == "inspect-mesh":
        parser.description = "Measure and query stable native cell IDs without recertifying an archive."
        parser.add_argument("--cell", type=int)
        parser.add_argument("--offset", type=int, default=0)
        parser.add_argument("--limit", type=int, default=100)
        parser.add_argument("--metric", choices=METRICS, default="mean_ratio")
        parser.add_argument("--quality-below", type=float)
        parser.add_argument("--region", nargs=6, type=float, metavar=("X0", "Y0", "Z0", "X1", "Y1", "Z1"))
        parser.add_argument("--spatial-mode", choices=("centroid", "aabb_overlap"), default="centroid")
        parser.add_argument("--body", type=int, action="append")
        parser.add_argument("--cell-type", choices=("all",) + CELL_NAMES, default="tet4")
        parser.add_argument("--orientation", choices=("positive", "inverted", "degenerate", "unsupported"))
        parser.add_argument("--sort-by", choices=("cell_id",) + METRICS, default="cell_id")
        parser.add_argument("--descending", action="store_true")
        parser.add_argument("--no-associations", action="store_true")
        parser.add_argument("--source", type=Path, help="matching STEP for optional trimmed-face centroid distances")
        parser.add_argument("--cad-distance", type=float, help="inclusive maximum centroid-to-selected-face distance")
        parser.add_argument("--cad-face", type=int, action="append")
        parser.add_argument("--max-cad-work", type=int, default=100000)
        parser.add_argument("--max-cad-seconds", type=float, default=10.)
        parser.add_argument("--cad-all-faces", action="store_true", help="include every selected face measurement in each returned record")
        parser.add_argument("--summary", action="store_true")
        parser.add_argument("--spatial-bins", type=int, default=4)
        parser.add_argument("--vtk", type=Path, help="atomically export the full snapshot with quality fields and native IDs, independent of page filters")
    elif command == "visualize":
        parser.description = "Render actual element IDs as deterministic SVG with explicit interior x-ray highlights."
        parser.add_argument("--output", type=Path, required=True)
        parser.add_argument("--source", type=Path, help="optional matching STEP source for display-tessellation context")
        parser.add_argument("--color-by", choices=METRICS, default="mean_ratio")
        parser.add_argument("--highlight-below", type=float, default=.2)
        parser.add_argument("--cell", action="append", type=int, help="explicit native cell IDs to isolate")
        parser.add_argument("--cut-axis", choices=("x", "y", "z"))
        parser.add_argument("--cut-value", type=float)
        parser.add_argument("--view-bounds", nargs=6, type=float)
        parser.add_argument("--azimuth", type=float, default=35.)
        parser.add_argument("--elevation", type=float, default=25.)
        parser.add_argument("--width", type=int, default=1200)
        parser.add_argument("--height", type=int, default=900)
        parser.add_argument("--max-triangles", type=int, default=100000)
        parser.add_argument("--max-highlights", type=int, default=10000)
        parser.add_argument("--max-output-bytes", type=int, default=32 * 1024 * 1024)
    else:
        parser.description = "Compare two separately reported native-valid attempts on the same CAD source."
        parser.add_argument("second", type=Path, help="second hash-paired inspection bundle")
        parser.add_argument("--spatial-bins", type=int, default=4)
        parser.add_argument("--output", type=Path, help="optional atomic JSON comparison export")
    return parser


def _human_mesh(report):
    print(f"PAACS {report['schema']}")
    if "mesh_sha256" in report:
        print(f"Mesh SHA-256: {report['mesh_sha256']}")
        print(f"Units: {report.get('units') or 'not stored in bare archive'}")
    if "cells" in report:
        print(f"Showing {len(report['cells'])} of {report['total']} matching elements.")
        for cell in report["cells"]:
            print(f"Cell {cell['cell_id']} ({cell['cell_type']}), body {cell['body_tag']}, {cell['orientation']}")
            print("  Metrics: " + json.dumps(cell["metrics"], allow_nan=False))
            print(f"  Nodes: {cell['node_ids']}; centroid: {cell['centroid']}")
        if report["next_offset"] is not None:
            print(f"Next page: --offset {report['next_offset']}")
    else:
        for key, value in report.items():
            if key not in ("schema", "mesh_sha256", "units", "metric_conventions"):
                print(f"{key.replace('_', ' ')}: {json.dumps(value, indent=2, allow_nan=False)}")


def _mesh_command(command, args):
    limits = {name: getattr(args, name) for name in ("max_bytes", "max_cells", "max_nodes", "max_work", "max_seconds")}
    with open_mesh(args.input, **limits) as inspection:
        if command == "inspect-mesh":
            if args.summary:
                if (args.cell is not None or args.quality_below is not None or args.region is not None
                        or args.body is not None or args.orientation is not None or args.cell_type != "tet4"
                        or args.offset != 0 or args.limit != 100 or args.sort_by != "cell_id" or args.descending
                        or args.metric != "mean_ratio" or args.spatial_mode != "centroid" or args.no_associations
                        or args.source is not None or args.cad_distance is not None or args.cad_face is not None or args.cad_all_faces):
                    raise ValueError("--summary describes the entire snapshot; use filtered pages or its spatial distributions for localization")
                report = inspection.summary(spatial_bins=args.spatial_bins)
            else:
                with _native_output_to_stderr():
                    cad = read_step(args.source, max_file_bytes=min(args.max_bytes // 2, 128 * 1024 * 1024)) if args.source else None
                report = inspection.query(cell_id=args.cell, offset=args.offset, limit=args.limit,
                    metric=args.metric, quality_below=args.quality_below,
                    region=[args.region[:3], args.region[3:]] if args.region is not None else None,
                    spatial_mode=args.spatial_mode, body_tags=args.body,
                    cell_type=None if args.cell_type == "all" else args.cell_type,
                    orientation=args.orientation, sort_by=args.sort_by, descending=args.descending,
                    associations=not args.no_associations, cad=cad, cad_distance=args.cad_distance,
                    cad_face_tags=args.cad_face, max_cad_work=args.max_cad_work, max_cad_seconds=args.max_cad_seconds,
                    cad_include_faces=args.cad_all_faces)
            if args.vtk is not None:
                report["exports"] = {"vtk": str(inspection.write_vtk(args.vtk))}
        elif command == "visualize":
            with _native_output_to_stderr():
                cad = read_step(args.source, max_file_bytes=min(args.max_bytes // 2, 128 * 1024 * 1024)) if args.source else None
                report = inspection.visualize(args.output, metric=args.color_by, highlight_below=args.highlight_below,
                    highlight_cell_ids=args.cell, cad=cad, cut_axis=args.cut_axis, cut_value=args.cut_value,
                    view_bounds=[args.view_bounds[:3], args.view_bounds[3:]] if args.view_bounds else None,
                    azimuth=args.azimuth, elevation=args.elevation, width=args.width, height=args.height,
                    max_triangles=args.max_triangles, max_highlights=args.max_highlights,
                    max_output_bytes=args.max_output_bytes)
        else:
            with open_mesh(args.second, **limits) as second:
                report = compare_meshes(inspection, second, spatial_bins=args.spatial_bins,
                                       max_work=args.max_work, max_seconds=args.max_seconds)
                if args.output is not None:
                    protected = inspection._archive.protected_paths + second._archive.protected_paths
                    atomic_bytes(args.output, json_text(report).encode("utf-8"), protected=protected)
        if args.json:
            _json(report)
        else:
            _human_mesh(report)
    return 0


def _human_capabilities(report):
    print(f"PAACS Mesh capabilities ({report.get('package_version', 'unknown version')})")
    # Use the same registry as machine clients, with readable field labels and
    # indentation rather than a separate hand-maintained capability list.
    for key, value in report.items():
        print(f"{key.replace('_', ' ').capitalize()}:")
        rendered = json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)
        for line in rendered.splitlines():
            print(f"  {line}")


def _human_inspection(report):
    print(f"CAD {report['kind']} ({report['units']})")
    print(f"Source SHA-256: {report['source_sha256']}")
    print(f"Showing {len(report['entities'])} of {report['total']} entities at offset {report['offset']}.")
    singular = {"faces": "Face", "bodies": "Body", "edges": "Edge", "vertices": "Vertex"}[report["kind"]]
    for record in report["entities"]:
        geometry = record.get("surface", record.get("curve"))
        print(f"{singular} {record['tag']}" + (f" ({geometry})" if geometry else ""))
        for key, value in record.items():
            if key == "tag":
                continue
            label = key.replace("_", " ").capitalize()
            print(f"  {label}: {json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True)}")
    if report["next_offset"] is not None:
        print(f"Next page: --offset {report['next_offset']}")


def _command_error(error, *, machine, native):
    message = str(error)
    report = {
        "schema": "paacs-mesh-cli-error-v1", "status": "failed",
        "stage": getattr(error, "stage", "input"),
        "reason": getattr(error, "reason", "invalid_argument"),
        "message": message[:4096], "message_truncated": len(message) > 4096,
        "report": getattr(error, "report", {}),
    }
    if machine:
        try:
            _json(report)
        except InspectionError:
            report["report"] = {}
            report["report_omitted_due_to_output_budget"] = True
            _json(report)
    else:
        print(f"{report['stage']}: {report['reason']}: {report['message']}", file=sys.stderr)
    return 1 if native else 2


def _advanced(command, argv):
    machine = "--json" in argv
    try:
        args = _command_parser(command).parse_args(argv)
        if command in ("inspect-mesh", "visualize", "compare"):
            return _mesh_command(command, args)
        if command == "capabilities":
            report = capabilities()
            if args.json:
                _json(report)
            else:
                _human_capabilities(report)
            return 0
        if command == "inspect":
            with _native_output_to_stderr():
                report = read_step(args.input).inspect(args.kind, offset=args.offset,
                    limit=args.limit, topology=args.topology)
            if args.json:
                _json(report)
            else:
                _human_inspection(report)
            return 0
        if command == "diagnose":
            report = diagnose(load_json(args.report))
            if args.json:
                _json(report)
            else:
                print(format_diagnostics(report))
            return 0 if report["status"] in ("valid", "valid_with_warnings") else 1
        if os.path.lexists(args.output):
            raise WorkflowError("output_exists", "Publication never replaces an existing output", stage="publication")
        with _native_output_to_stderr():
            result = run_saved_recipe(args.recipe, args.source)
        diagnostics = diagnose(result.report)
        if diagnostics["status"] not in ("valid", "valid_with_warnings"):
            if args.json:
                _json(diagnostics)
            else:
                print(format_diagnostics(diagnostics), file=sys.stderr)
            return 1
        with _native_output_to_stderr():
            publication = publish_result(result, args.output)
        report = {"schema": "paacs-mesh-run-v1", "status": "succeeded",
                  "publication": publication, "diagnostics": diagnostics}
        if args.json:
            _json(report)
        else:
            print(format_diagnostics(diagnostics))
            print(f"Saved result bundle: {publication['path']}")
        return 0
    except (CadImportError, MeshingError, ExportError, SelectionError, WorkflowError, InspectionError) as error:
        return _command_error(error, machine=machine,
                              native=not isinstance(error, SelectionError) or error.reason != "invalid_argument")
    except (ValueError, TypeError, OverflowError, OSError) as error:
        return _command_error(error, machine=machine, native=False)


def main(argv=None):
    """Return 0 for success, 1 for a refused operation, and 2 for invalid input."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in _COMMANDS:
        return _advanced(argv[0], argv[1:])
    return _legacy(argv)


if __name__ == "__main__":
    raise SystemExit(main())
