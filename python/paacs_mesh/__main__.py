"""Mesh a STEP solid using the same native API as Python callers."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
import json
from pathlib import Path
import sys

from . import BoxSize, CadImportError, ExportError, FaceSize, MeshingError, read_step


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
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
    parser.add_argument("--max-work", type=int, help="native meshing and validation work budget")
    parser.add_argument("--max-seconds", type=float, help="cooperative native cancellation deadline")
    parser.add_argument("--warning-quality", type=float,
                        help="preferred mean-ratio quality in [0,1]; 0 disables warnings")
    parser.add_argument("--minimum-quality", type=float,
                        help="strict minimum mean-ratio quality in [0,1]")
    args = parser.parse_args(argv)
    try:
        model = read_step(args.input)
        face_sizes = [FaceSize(tag=int(tag), size=float(size)) for tag, size in args.face_size]
        volume_controls = [BoxSize(lower=tuple(box[:3]), upper=tuple(box[3:6]), size=box[6])
                           for box in args.size_box]
        mesh = model.mesh(args.size, recovery=args.recovery, assembly=args.assembly,
                          face_sizes=face_sizes, volume_controls=volume_controls,
                          min_size=args.min_size, growth=args.growth,
                          native_sizing=args.native_sizing, max_work=args.max_work,
                          max_seconds=args.max_seconds, warning_quality=args.warning_quality,
                          minimum_quality=args.minimum_quality)
        mesh.write(args.output)
        print(json.dumps(mesh.report, indent=2, allow_nan=False))
        return 0
    except (CadImportError, MeshingError, ExportError) as error:
        print(json.dumps({"status": "failed", "stage": error.stage,
                          "reason": error.reason, "report": error.report},
                         allow_nan=False), file=sys.stderr)
        return 1
    except (ValueError, TypeError, OverflowError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
