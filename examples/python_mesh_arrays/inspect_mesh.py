#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Inspect native Tet4 arrays and save coordinates, connectivity and CAD tags."""

import argparse
import json
from pathlib import Path

import numpy as np
import paacs_mesh as pm


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--size", type=float, default=0.3, help="target size in mm")
    parser.add_argument("--face-tag", type=int, help="CAD face tag to inspect")
    args = parser.parse_args()
    if args.output.suffix.lower() != ".npz":
        parser.error("output must use the .npz extension")
    if args.output.exists():
        parser.error("output already exists; choose a new destination")

    mesh = pm.read_step(args.input).mesh(args.size)
    arrays = mesh.to_numpy()
    tets = arrays["cells"]["tet4"]
    coordinates = arrays["points"][tets["connectivity"]]
    edges = coordinates[:, 1:] - coordinates[:, :1]
    volumes = np.einsum("ij,ij->i", edges[:, 0],
                        np.cross(edges[:, 1], edges[:, 2])) / 6
    summary = {
        "units": arrays["units"],
        "nodes": len(arrays["points"]),
        "tetrahedra": len(volumes),
        "minimum_signed_volume": float(volumes.min()),
        "total_volume": float(volumes.sum()),
        "native_total_volume": mesh.report["quality"]["total_volume"],
        "body_tags": mesh.body_tags,
        "mesh_scope": mesh.report["mesh_scope"],
        "interbody_nonoverlap_certified": mesh.report["interbody_nonoverlap_certified"],
    }
    if args.face_tag is not None:
        if args.face_tag not in mesh.face_tags:
            parser.error("face tag is absent from the mesh")
        triangles = arrays["cells"].get("tri3")
        summary["selected_face"] = args.face_tag
        summary["selected_triangles"] = int(np.count_nonzero(
            triangles["entity_tags"] == args.face_tag)) if triangles is not None else 0

    saved = {name: arrays[name] for name in (
        "points", "node_entity_dimensions", "node_entity_tags")}
    for kind, block in arrays["cells"].items():
        saved.update({f"{kind}_{name}": value for name, value in block.items()})
    summary["array_bytes"] = sum(value.nbytes for value in saved.values())
    saved["metadata"] = np.asarray(json.dumps(summary, allow_nan=False))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("xb") as stream:
        np.savez_compressed(stream, **saved)
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
