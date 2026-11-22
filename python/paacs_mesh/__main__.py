"""Inspect STEP geometry, run saved recipes or mesh directly with PAACS Mesh."""

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
    WorkflowError, publish_result, read_step, run_saved_recipe,
)
from ._cad_common import load_json
from .diagnostics import capabilities, diagnose, format_diagnostics


_COMMANDS = ("capabilities", "inspect", "diagnose", "run")


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
    print(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), file=stream)


def _legacy(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog="Discovery commands: capabilities, inspect, diagnose, run. "
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
    parser.add_argument("--growth", type=float, default=1.3,
                        help="size-field growth, at least 1 (default: 1.3)")
    parser.add_argument("--native-sizing", action=argparse.BooleanOptionalAction, default=None,
                        help="force spatial sizing or retain legacy cap sizing; automatic with controls")
    parser.add_argument("--advancing-interior-front", action=argparse.BooleanOptionalAction, default=None,
                        help="select advancing initial Tet4 placement or the legacy constructor; automatic by default")
    parser.add_argument("--max-work", type=int, help="native meshing and validation work budget")
    parser.add_argument("--max-seconds", type=float, help="cooperative native cancellation deadline")
    parser.add_argument("--warning-quality", type=float,
                        help="preferred mean-ratio quality in [0,1]; 0 disables warnings")
    parser.add_argument("--minimum-quality", type=float,
                        help="strict minimum mean-ratio quality in [0,1]")
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
                              advancing_interior_front=args.advancing_interior_front,
                              max_seconds=args.max_seconds, warning_quality=args.warning_quality,
                              minimum_quality=args.minimum_quality)
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
    report = {
        "schema": "paacs-mesh-cli-error-v1", "status": "failed",
        "stage": getattr(error, "stage", "input"),
        "reason": getattr(error, "reason", "invalid_argument"),
        "message": str(error), "report": getattr(error, "report", {}),
    }
    if machine:
        _json(report)
    else:
        print(f"{report['stage']}: {report['reason']}: {report['message']}", file=sys.stderr)
    return 1 if native else 2


def _advanced(command, argv):
    machine = "--json" in argv
    try:
        args = _command_parser(command).parse_args(argv)
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
    except (CadImportError, MeshingError, ExportError, SelectionError, WorkflowError) as error:
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
