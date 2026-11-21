"""Mesh a STEP solid using the same native API as Python callers."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
import json
from pathlib import Path
import sys

from . import CadImportError, ExportError, MeshingError, read_step


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="STEP file; transferred to mm")
    parser.add_argument("output", type=Path, help=".vtk or .paacsmesh destination")
    parser.add_argument("--size", type=float, required=True, help="global target size in mm")
    parser.add_argument("--recovery", choices=("standard", "advanced"), default="standard")
    parser.add_argument("--assembly", choices=("preserve_parts", "conforming"), default="preserve_parts")
    parser.add_argument("--max-seconds", type=float, help="cooperative native cancellation deadline")
    args = parser.parse_args(argv)
    try:
        model = read_step(args.input)
        mesh = model.mesh(args.size, recovery=args.recovery, assembly=args.assembly,
                          max_seconds=args.max_seconds)
        mesh.write(args.output)
        print(json.dumps(mesh.report, indent=2, allow_nan=False))
        return 0
    except (CadImportError, MeshingError, ExportError) as error:
        print(json.dumps({"status": "failed", "stage": error.stage,
                          "reason": error.reason, "report": error.report},
                         allow_nan=False), file=sys.stderr)
        return 1
    except (ValueError, TypeError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
