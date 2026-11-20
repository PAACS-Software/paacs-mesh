# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Write and reopen the C3D4 mesh include used by the plate study."""

from collections import defaultdict
from pathlib import Path
import re

import numpy as np

from scripts.mesh_visual_data import Cell, Mesh


# C3D4 element surface identifiers, expressed in outward winding for a
# positive determinant. The side IDs follow the Abaqus solid element table.
C3D4_SIDES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))


def exterior_sides(tetrahedra):
    sides = {}
    closed = set()
    for element, tet in enumerate(tetrahedra, 1):
        for side, shape in enumerate(C3D4_SIDES, 1):
            face = tuple(tet[i] for i in shape)
            key = tuple(sorted(face))
            if key in closed:
                raise ValueError("Nonmanifold tetrahedron face")
            if key in sides:
                del sides[key]
                closed.add(key)
            else:
                sides[key] = (element, side)
    return sides


def write_solver_input(mesh, path):
    tets = [cell for cell in mesh.cells if cell.dim == 3]
    surface = [cell for cell in mesh.cells if cell.dim == 2]
    if not tets or any(cell.kind != 5 for cell in tets) or any(cell.kind != 2 for cell in surface):
        raise ValueError("The plate solver include requires Tet4 with a Tri3 boundary")
    exterior = exterior_sides([cell.ids for cell in tets])
    bodies, surfaces = defaultdict(list), defaultdict(list)
    for i, cell in enumerate(tets, 1):
        bodies[cell.tag].append(i)
    for cell in surface:
        owner, side = exterior.pop(tuple(sorted(cell.ids)), (None, None))
        if owner is None:
            raise ValueError("A CAD-tagged triangle has no unique exterior Tet4 side")
        surfaces[cell.tag].append((owner, side))
    if exterior:
        raise ValueError("Some exterior tetrahedron faces lack CAD tags")
    with Path(path).open("x", encoding="ascii", newline="\n") as f:
        f.write("** PAACS mesh include. Coordinates in mm.\n*NODE\n")
        for i, point in enumerate(mesh.points, 1):
            f.write(f"{i}, {point[0]:.17g}, {point[1]:.17g}, {point[2]:.17g}\n")
        f.write("*ELEMENT, TYPE=C3D4\n")
        for i, cell in enumerate(tets, 1):
            f.write(", ".join(map(str, (i, *(node+1 for node in cell.ids)))) + "\n")
        for tag, elements in sorted(bodies.items()):
            f.write(f"*ELSET, ELSET=CAD_BODY_{tag}\n")
            for start in range(0, len(elements), 16):
                f.write(", ".join(map(str, elements[start:start+16])) + "\n")
        for tag, faces in sorted(surfaces.items()):
            f.write(f"*SURFACE, NAME=CAD_FACE_{tag}, TYPE=ELEMENT\n")
            for element, side in sorted(faces):
                f.write(f"{element}, S{side}\n")


def read_solver_input(path):
    """Read the plate mesh include subset and recover the original CAD groups."""
    nodes, elements, bodies, surfaces = {}, {}, {}, {}
    block, tag = None, None
    path = Path(path)
    if path.stat().st_size > 128 << 20:
        raise ValueError("Plate mesh include exceeds byte budget")
    for line in path.read_text(encoding="ascii").splitlines():
        line = line.strip().upper()
        if not line or line.startswith("**"):
            continue
        if line.startswith("*"):
            if line == "*NODE":
                block = "node"
            elif line == "*ELEMENT, TYPE=C3D4":
                block = "element"
            elif found := re.fullmatch(r"\*ELSET, ELSET=CAD_BODY_(\d+)", line):
                block, tag = "body", int(found[1])
                bodies.setdefault(tag, [])
            elif found := re.fullmatch(r"\*SURFACE, NAME=CAD_FACE_(\d+), TYPE=ELEMENT", line):
                block, tag = "surface", int(found[1])
                surfaces.setdefault(tag, [])
            else:
                raise ValueError(f"Unsupported plate mesh include keyword {line}")
            continue
        fields = line.split(",")
        if block == "node":
            index = int(fields[0])
            coordinates = tuple(map(float, fields[1:]))
            if len(coordinates) != 3 or index in nodes:
                raise ValueError("Invalid solver node record")
            nodes[index] = coordinates
        elif block == "element":
            index, *conn = map(int, fields)
            if len(conn) != 4 or len(set(conn)) != 4 or index in elements:
                raise ValueError("Invalid C3D4 record")
            elements[index] = conn
        elif block == "surface" and len(fields) == 2 and re.fullmatch(r"S[1-4]", fields[1].strip()):
            surfaces[tag].append((int(fields[0]), int(fields[1].strip()[1])))
        elif block == "body":
            bodies[tag].extend(map(int, fields))
        else:
            raise ValueError("Input record outside a solver keyword block")
    point_ids = {tag: i for i, tag in enumerate(sorted(nodes))}
    points = np.array([nodes[tag] for tag in sorted(nodes)], dtype=float)
    if not len(points) or not np.isfinite(points).all():
        raise ValueError("Missing or nonfinite solver coordinates")
    owner = {}
    for tag, group in bodies.items():
        for element in group:
            if element in owner:
                raise ValueError("A solver element belongs to multiple body groups")
            owner[element] = tag
    if not elements or owner.keys() != elements.keys():
        raise ValueError("Solver body groups do not cover all C3D4 elements")
    tets = {tag: tuple(point_ids[node] for node in conn) for tag, conn in sorted(elements.items())}
    cells = [Cell(5, 3, owner[tag], conn) for tag, conn in tets.items()]
    exported_faces = set()
    for tag, group in sorted(surfaces.items()):
        for element, side in group:
            face = tuple(tets[element][i] for i in C3D4_SIDES[side-1])
            key = tuple(sorted(face))
            if key in exported_faces:
                raise ValueError("Duplicate solver surface facet")
            exported_faces.add(key)
            cells.append(Cell(2, 2, tag, face))
    return Mesh(points, cells, Path(path))
