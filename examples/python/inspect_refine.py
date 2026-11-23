"""Inspect an owned STEP fixture, compare two explicit sizing recipes, and verify them.

Run from a source checkout with PAACS Mesh installed::

    python examples/python/inspect_refine.py --output out/inspect-refine

The default cylinder was created for this project. An alternative STEP file must
contain exactly one cylindrical face of radius 1 mm; ambiguous selections stop.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse
import json
from pathlib import Path

import paacs_mesh as pm


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def verify(model, result):
    """Keep shape warnings separate from mandatory native validity checks."""
    report = result.mesh.report
    if result.source_sha256 != model.source_sha256 or report["source_sha256"] != model.source_sha256:
        raise RuntimeError("The result no longer identifies the imported STEP bytes")
    if not report["mesh_published"] or report["reason"] != "none":
        raise RuntimeError("The native mesher did not publish a successful result")
    for name in ("audit", "quality"):
        if not report[name]["complete"] or not report[name]["valid"]:
            raise RuntimeError(f"The mandatory {name} check was incomplete or invalid")
    if sorted(result.mesh.body_tags) != sorted(model.body_tags):
        raise RuntimeError("Body tags changed")
    if sorted(result.mesh.face_tags) != sorted(model.face_tags):
        raise RuntimeError("Face tags changed")
    return {
        "nodes": result.mesh.node_count,
        "tetrahedra": result.mesh.tetrahedron_count,
        "quality": report["quality"],
        "seconds_including_retries": result.report["seconds"],
        "attempts": len(result.report["attempts"]),
        "cad_error_by_body": [
            {key: body[key] for key in (
                "solid_tag", "relative_volume_error", "maximum_boundary_node_error",
                "maximum_boundary_sample_error", "surface_chord_error_limit",
            )}
            for body in report["bodies"]
        ],
        "face_tags": result.mesh.face_tags,
        "body_tags": result.mesh.body_tags,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path,
                        default=Path(__file__).resolve().parents[2] / "tests/fixtures/cad/cylinder.step")
    parser.add_argument("--output", type=Path, default=Path("out/inspect-refine"),
                        help="new output directory; existing destinations are refused")
    args = parser.parse_args(argv)
    model = pm.read_step(args.source)

    # Inspect before choosing a tag. This is one bounded page, with continuation
    # recorded if a different input has more than 100 faces.
    inspection = model.inspect("faces", limit=100, topology=True)
    wall = model.select_faces(surface="cylinder", radius=1.0, expect=1)
    model.define_group("cylinder_wall", wall, group_id="inspect-refine-cylinder-wall")

    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output / "inspection.json", inspection)
    model.preview(wall).write_svg(args.output / "selected-wall.svg")
    comparison = {"source_sha256": model.source_sha256, "units": "mm", "runs": {}}

    # The second attempt is a deliberate comparison, not an automatic diagnosis
    # that smaller elements will improve shape quality or simulation accuracy.
    for name, controls in (("baseline", []), ("wall-refined", [(wall, .22)])):
        recipe = pm.create_recipe(
            model, .45, selected_face_sizes=controls,
            assembly="conforming", native_sizing=True, min_size=.18, growth=1.3,
            workers=1, validation_workers=1, max_work=100_000_000,
            max_attempts=1, max_seconds=90, warning_quality=.25,
        )
        path = args.output / f"{name}.recipe.json"
        pm.save_recipe(model, recipe, path)
        # Re-import from an explicit local path; the saved display name is never
        # treated as a path, and recipe identity/property guards are rechecked.
        result = pm.run_saved_recipe(path, args.source)
        diagnostics = pm.diagnose(result.mesh)
        write_json(args.output / f"{name}.diagnostics.json", diagnostics)
        print(name)
        print(pm.format_diagnostics(diagnostics))
        comparison["runs"][name] = verify(model, result)
        pm.publish_result(result, args.output / name)

    write_json(args.output / "comparison.json", comparison)
    print(f"Verified source identity, native validity and CAD tags; saved {args.output}")
    print("Compare measured quality and CAD errors before choosing a sizing recipe.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
