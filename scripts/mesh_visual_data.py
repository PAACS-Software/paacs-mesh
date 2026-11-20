# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Measured mesh data and plane intersections for mesh_visual_report.py.

The corner/edge/face orderings mirror include/paacs/topo/ElementTraits.hpp.
No tetrahedralization or generated replacement mesh is used for display.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import math
import struct

import numpy as np

# PAACS CellType, not ElementType: Prism6 precedes Hex8 in the archive.
ARITY = (2, 3, 3, 6, 4, 4, 10, 6, 8, 5)
NAMES = ('Line2', 'Line3', 'Tri3', 'Tri6', 'Quad4', 'Tet4', 'Tet10',
         'Prism6', 'Hex8', 'Pyramid5')
DIMENSION = (1, 1, 2, 2, 2, 3, 3, 3, 3, 3)
CORNERS = (2, 2, 3, 3, 4, 4, 4, 6, 8, 5)
TET_EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
EDGES = {
    5: TET_EDGES, 6: TET_EDGES,
    7: ((0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)),
    8: ((0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),
        (0,4),(1,5),(2,6),(3,7)),
    9: ((0,1),(1,2),(2,3),(3,0),(0,4),(1,4),(2,4),(3,4)),
}
TET_FACES = ((1,2,3),(0,3,2),(0,1,3),(0,2,1))
FACES = {
    5: TET_FACES, 6: TET_FACES,
    7: ((0,2,1),(3,4,5),(0,1,4,3),(1,2,5,4),(2,0,3,5)),
    8: ((0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),
        (2,3,7,6),(3,0,4,7)),
    9: ((0,3,2,1),(0,1,4),(1,2,4),(2,3,4),(3,0,4)),
}
# The fixed diagnostic decompositions in src/algo/MixedMeshQuality.cpp.
# These only measure quality; the displayed mesh is never replaced by them.
QUALITY_SUBTETS = {
    7: ((0,1,2,3),(1,4,2,3),(2,4,5,3)),
    8: ((0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)),
    9: ((0,1,2,4),(0,2,3,4)),
}
VTK_KIND = {3:0,21:1,5:2,22:3,9:4,10:5,24:6,13:7,12:8,14:9}


@dataclass
class Cell:
    kind: int
    dim: int
    tag: int
    ids: tuple


@dataclass
class Mesh:
    points: np.ndarray
    cells: list
    source: Path | None = None

    @property
    def bounds(self):
        return self.points.min(axis=0), self.points.max(axis=0)

    @property
    def scale(self):
        lo, hi = self.bounds
        return float(np.linalg.norm(hi-lo))


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_archive(path, max_bytes=512 << 20, max_nodes=2_000_000,
                 max_cells=4_000_000):
    path = Path(path)
    if path.stat().st_size > max_bytes:
        raise ValueError('archive exceeds visualization byte budget')
    data = path.read_bytes()
    if len(data) < 40 or data[:8] != b'PAACSM8\0':
        raise ValueError('invalid PAACS archive header')
    version, reserved, n, c, connectivity = struct.unpack_from('<IIQQQ', data, 8)
    if version != 1 or reserved:
        raise ValueError('unsupported archive version or nonzero reserved field')
    if n > max_nodes or c > max_cells:
        raise ValueError('archive exceeds visualization node/cell budget')
    if 40 + 32*n + 12*c + 8*connectivity != len(data):
        raise ValueError('archive count/byte mismatch')
    nodes = np.frombuffer(data, dtype=np.dtype([('xyz','<f8',(3,)),
                          ('dim','<i4'),('tag','<i4')]), count=n, offset=40)
    points = nodes['xyz'].copy()
    if not np.isfinite(points).all() or np.any(nodes['dim'] < 0) or \
            np.any(nodes['dim'] > 3) or np.any(nodes['tag'] < 0):
        raise ValueError('invalid node coordinates/classification')
    cells, at, used = [], 40 + 32*n, 0
    for _ in range(c):
        if at + 12 > len(data):
            raise ValueError('truncated cell header')
        kind, r0, r1, r2, dim, tag = struct.unpack_from('<BBBBii', data, at)
        if kind >= len(ARITY) or r0 or r1 or r2 or not 0 <= dim <= 3 or tag < 0:
            raise ValueError('invalid cell type/classification/reserved fields')
        width = ARITY[kind]
        if at + 12 + 8*width > len(data):
            raise ValueError('truncated cell connectivity')
        ids = struct.unpack_from('<'+'q'*width, data, at+12)
        if min(ids) < 0 or max(ids) >= n:
            raise ValueError('out-of-range node index')
        cells.append(Cell(kind, dim, tag, ids))
        at += 12 + 8*width
        used += width
    if at != len(data) or used != connectivity:
        raise ValueError('trailing data or inconsistent connectivity count')
    return Mesh(points, cells, path.resolve())


def read_vtk(path):
    """PAACS legacy ASCII UNSTRUCTURED_GRID; does not silently linearize types."""
    path = Path(path)
    lines = path.read_text(encoding='utf-8').splitlines()
    if len(lines) < 4 or lines[2].strip() != 'ASCII' or \
            lines[3].strip() != 'DATASET UNSTRUCTURED_GRID':
        raise ValueError('expected legacy ASCII VTK UNSTRUCTURED_GRID')
    tokens, at = ' '.join(lines[4:]).split(), 0
    if tokens[at] != 'POINTS':
        raise ValueError('missing VTK POINTS')
    n = int(tokens[at+1]); at += 3
    points = np.asarray(tokens[at:at+3*n], dtype=float).reshape(n, 3); at += 3*n
    if tokens[at] != 'CELLS':
        raise ValueError('missing VTK CELLS')
    c, expected = int(tokens[at+1]), int(tokens[at+2]); at += 3
    ids, used = [], 0
    for _ in range(c):
        width = int(tokens[at]); at += 1
        ids.append(tuple(map(int, tokens[at:at+width])))
        at += width; used += width+1
    if used != expected or tokens[at] != 'CELL_TYPES' or int(tokens[at+1]) != c:
        raise ValueError('VTK connectivity/type count mismatch')
    at += 2
    kinds = [VTK_KIND[int(x)] for x in tokens[at:at+c]]
    cells = [Cell(kind, DIMENSION[kind], 0, conn) for kind, conn in zip(kinds, ids)]
    if not np.isfinite(points).all() or any(len(x.ids) != ARITY[x.kind] or
            min(x.ids) < 0 or max(x.ids) >= n for x in cells):
        raise ValueError('invalid VTK mesh')
    return Mesh(points, cells, path.resolve())


def read_mesh(path):
    return read_vtk(path) if Path(path).suffix.lower() == '.vtk' else read_archive(path)


def read_geometry(path):
    path = Path(path)
    if path.suffix.lower() != '.json':
        return read_mesh(path)
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    points = np.asarray(obj['points'], dtype=float)
    tags = obj.get('face_tags', [0]*len(obj['triangles']))
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all() or \
            len(tags) != len(obj['triangles']):
        raise ValueError('invalid tessellation coordinates/tags')
    cells = [Cell(2, 2, int(tag), tuple(map(int, ids)))
             for ids, tag in zip(obj['triangles'], tags)]
    if any(len(x.ids) != 3 or min(x.ids) < 0 or max(x.ids) >= len(points) for x in cells):
        raise ValueError('invalid tessellation triangle')
    return Mesh(points, cells, path.resolve())


def volume_data(mesh):
    """Corner polyhedra only, after planar/convex and Tet10-midpoint checks.

    Mixed volume is the divergence sum over the retained trait face loops.
    Quality is independently calculated mean ratio for Tet4/Tet10 corners,
    and the native fixed SubTetMeanRatio diagnostic for linear mixed cells.
    """
    volumes, quality, volume_cells, exterior, counts = [], [], [], {}, {}
    tolerance = max(mesh.scale*2e-10, np.finfo(float).tiny)
    tet_indices = [i for i,c in enumerate(mesh.cells) if c.kind in (5,6)]
    tet_measurements = {}
    if tet_indices:
        p = mesh.points[np.asarray([mesh.cells[i].ids[:4] for i in tet_indices])]
        v = np.linalg.det(p[:,1:]-p[:,0,None])/6
        if not np.isfinite(v).all() or np.any(v <= 0):
            raise ValueError('nonpositive tetrahedral corner volume')
        edge_sum = sum(np.sum((p[:,a]-p[:,b])**2, axis=1) for a,b in TET_EDGES)
        q = 12*(3*v)**(2/3)/edge_sum
        tet_measurements = {ci:(float(vol),float(quality))
                            for ci,vol,quality in zip(tet_indices,v,q)}
        high = [i for i in tet_indices if mesh.cells[i].kind == 6]
        if high:
            hp = mesh.points[np.asarray([mesh.cells[i].ids for i in high])]
            for e,(a,b) in enumerate(TET_EDGES):
                if np.max(np.linalg.norm(hp[:,4+e]-(hp[:,a]+hp[:,b])*.5, axis=1)) > tolerance:
                    raise ValueError('curved Tet10 requires a curved-cell renderer; corner section refused')
    for ci, cell in enumerate(mesh.cells):
        counts[NAMES[cell.kind]] = counts.get(NAMES[cell.kind], 0) + 1
        if cell.kind not in FACES:
            continue
        if cell.kind in (5,6):
            volume, measured_quality = tet_measurements[ci]
            volumes.append(volume); quality.append(measured_quality); volume_cells.append(ci)
            for face in TET_FACES:
                ids = tuple(cell.ids[j] for j in face)
                key = tuple(sorted(ids))
                if key not in exterior:
                    exterior[key] = [ids, 1]
                else:
                    exterior[key][1] += 1
            continue
        p = mesh.points[list(cell.ids[:CORNERS[cell.kind]])]
        if cell.kind == 6:
            for e, (a,b) in enumerate(TET_EDGES):
                if np.linalg.norm(mesh.points[cell.ids[4+e]]-(p[a]+p[b])*.5) > tolerance:
                    raise ValueError('curved Tet10 requires a curved-cell renderer; corner section refused')
        cell_scale = float(np.linalg.norm(np.ptp(p, axis=0)))
        local_tolerance = max(cell_scale*2e-10, np.finfo(float).tiny)
        contributions = []
        origin = p[0]
        for face in FACES[cell.kind]:
            poly = p[list(face)]
            normal = np.cross(poly[1]-poly[0], poly[2]-poly[0])
            length = np.linalg.norm(normal)
            if not length > 0:
                raise ValueError(f'degenerate face in cell {ci}')
            distance = (p-poly[0]) @ (normal/length)
            if np.max(distance) > local_tolerance or \
                    np.max(np.abs(distance[list(face)])) > local_tolerance:
                raise ValueError(f'nonconvex/nonplanar/inverted corner cell {ci}; exact section refused')
            for k in range(1, len(face)-1):
                a,b,d = poly[0]-origin, poly[k]-origin, poly[k+1]-origin
                contributions.append(float(np.dot(a, np.cross(b,d)))/6)
            ids = tuple(cell.ids[j] for j in face)
            key = tuple(sorted(ids))
            if key not in exterior:
                exterior[key] = [ids, 1]
            else:
                exterior[key][1] += 1
        volume = math.fsum(contributions)
        if not volume > 0 or not math.isfinite(volume):
            raise ValueError(f'nonpositive corner volume in cell {ci}')
        volumes.append(volume); volume_cells.append(ci)
        subtet_quality, subtet_volume = [], []
        for tet in QUALITY_SUBTETS[cell.kind]:
            t = p[list(tet)]
            tv = float(np.linalg.det(t[1:]-t[0]))/6
            if tv <= 0:
                raise ValueError(f'nonpositive native diagnostic subtet in cell {ci}')
            edge_sum = math.fsum(float(np.sum((t[a]-t[b])**2)) for a,b in TET_EDGES)
            subtet_quality.append(12*(3*tv)**(2/3)/edge_sum)
            subtet_volume.append(tv)
        if abs(math.fsum(subtet_volume)-volume) > volume*1e-9:
            raise ValueError('native quality decomposition disagrees with corner polyhedron volume')
        quality.append(min(subtet_quality))
    if any(count > 2 for _,count in exterior.values()):
        raise ValueError('nonmanifold volume-face incidence')
    boundary = [ids for ids,count in exterior.values() if count == 1]
    v, q = np.asarray(volumes), np.asarray(quality)
    if not len(v):
        raise ValueError('mesh contains no supported volume cells')
    sizes = np.cbrt(v)
    tetq = q[[mesh.cells[i].kind in (5,6) for i in volume_cells]]
    metrics = {
        'nodes':len(mesh.points), 'cells':len(mesh.cells), 'volume_cells':len(v),
        'cell_types':counts, 'positive_corner_volumes':True,
        'convex_planar_corner_cells':True, 'exterior_volume_faces':len(boundary),
        'volume_sum':math.fsum(v.tolist()), 'volume_cv':float(v.std()/v.mean()),
        'volume_min':float(v.min()), 'volume_max':float(v.max()),
        'volume_max_min':float(v.max()/v.min()),
        'equivalent_size_definition':'cube root of positive cell volume; not target edge size',
        'equivalent_size_cv':float(sizes.std()/sizes.mean()),
        'volume_quantiles':np.quantile(v,[0,.1,.5,.9,1]).tolist(),
        'equivalent_size_quantiles':np.quantile(sizes,[0,.1,.5,.9,1]).tolist(),
        'tet_mean_ratio_samples':len(tetq),
        'tet_mean_ratio_min':float(tetq.min()) if len(tetq) else None,
        'tet_mean_ratio_quantiles':np.quantile(tetq,[0,.1,.5,.9,1]).tolist() if len(tetq) else [],
        'quality_metric':'Tet corner mean ratio; linear mixed cells use native fixed minimum SubTetMeanRatio diagnostic',
        'subtet_mean_ratio_min':float(q.min()),
        'subtet_mean_ratio_quantiles':np.quantile(q,[0,.1,.5,.9,1]).tolist(),
        'high_order_scope':'Tet10 sections/quality use corners only; straight midsides verified; faceted geometry',
        'audit_scope':'Independent corner geometry and face multiplicity checks, not replacement for native full audit',
    }
    return v, q, volume_cells, boundary, metrics


def section(mesh, volume_cells, volumes, axis=1, value=None):
    lo, hi = mesh.bounds
    if value is None:
        value = (lo[axis]+hi[axis])*.5
    axes = [i for i in range(3) if i != axis]
    polygons, cell_numbers, seen = [], [], set()
    intersections = {}
    duplicate_coplanar, intersected = 0, 0
    for vi, ci in enumerate(volume_cells):
        cell = mesh.cells[ci]
        p = mesh.points[list(cell.ids[:CORNERS[cell.kind]])]
        distances = p[:,axis]-value
        if p[:,axis].min() > value or p[:,axis].max() < value:
            continue
        vertices = {}
        for a,b in EDGES[cell.kind]:
            da, db = distances[a], distances[b]
            if da == 0:
                vertices[('v',cell.ids[a])] = p[a].copy()
            if (da < 0 < db) or (db < 0 < da):
                key = ('e',*sorted((cell.ids[a],cell.ids[b])))
                if key not in intersections:
                    intersections[key] = p[a]+(p[b]-p[a])*(-da)/(db-da)
                    intersections[key][axis] = value
                vertices[key] = intersections[key]
            if db == 0:
                vertices[('v',cell.ids[b])] = p[b].copy()
        if len(vertices) < 3:
            continue
        poly = np.asarray(list(vertices.values()))
        projected = poly[:,axes]
        center = projected.mean(axis=0)
        order = np.argsort(np.arctan2(projected[:,1]-center[1], projected[:,0]-center[0]))
        poly = poly[order]
        area = polygon_area(poly[:,axes])
        if area == 0:
            continue
        intersected += 1
        # A section coinciding with a shared whole face belongs to two cells.
        # Draw that one physical polygon once and record the deduplication.
        key = tuple(sorted(vertices))
        if key in seen:
            duplicate_coplanar += 1
            continue
        seen.add(key); polygons.append(poly); cell_numbers.append(vi)
    areas = [polygon_area(poly[:,axes]) for poly in polygons]
    metrics = {'axis':'xyz'[axis], 'value':float(value), 'polygon_count':len(polygons),
               'intersected_cells_before_shared_face_dedup':intersected,
               'coplanar_shared_polygons_removed':duplicate_coplanar,
               'section_area':math.fsum(areas),
               'maximum_plane_residual':max((float(np.max(np.abs(p[:,axis]-value)))
                                              for p in polygons), default=0),
               'method':'Stored-coordinate convex cell-edge/plane intersections; shared mesh identities deduplicated without plane snapping'}
    return polygons, np.asarray(cell_numbers, dtype=int), metrics


def polygon_area(poly):
    # Translate to avoid cancellation for small cells far from the origin.
    p = poly-poly[0]
    return .5*abs(float(np.dot(p[:,0], np.roll(p[:,1],-1))-
                           np.dot(p[:,1], np.roll(p[:,0],-1))))


def clip_polygon(poly, axis, value, keep=1):
    """Clip an existing exterior face, never replace it by fabricated tets."""
    out = []
    for i,p in enumerate(poly):
        q = poly[(i+1)%len(poly)]
        dp,dq = keep*(p[axis]-value), keep*(q[axis]-value)
        if dp >= 0:
            out.append(p)
        if (dp >= 0) != (dq >= 0):
            out.append(p+(q-p)*(-dp)/(dq-dp))
    return np.asarray(out)


def self_test():
    points = np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
    mesh = Mesh(points, [Cell(5,3,1,(0,1,2,3))])
    v, _, ids, _, _ = volume_data(mesh)
    assert abs(v[0]-1/6) < 1e-15
    polygons, _, metric = section(mesh, ids, v, axis=2, value=.5)
    assert len(polygons) == 1 and len(polygons[0]) == 3
    assert abs(metric['section_area']-.125) < 1e-15
    # A nearby positive thin tetrahedron is not intersected merely because
    # its vertices fall inside a display epsilon around the section plane.
    plane = .13
    thin = np.array([[plane+1e-13,0,0],[plane+2e-13,0,0],
                     [plane+1e-13,1,0],[plane+1e-13,0,1]])
    tiny = Mesh(thin,[Cell(5,3,1,(0,1,2,3))])
    cuts,_,_ = section(tiny,[0],np.ones(1),axis=0,value=plane)
    assert len(cuts) == 0
    thin[0,0] = plane-1e-13
    tiny = Mesh(thin,[Cell(5,3,1,(0,1,2,3))])
    cuts,_,metric = section(tiny,[0],np.ones(1),axis=0,value=plane)
    assert len(cuts) == 1 and len(cuts[0]) == 3
    assert abs(metric['section_area']-.125) < 1e-15
    # This plane cuts four actual edges: the result is a quadrilateral, not
    # an arbitrary triangle from a centroid-selected tetrahedron.
    rotated = points.copy(); rotated[:,2] = [0,0,1,1]
    mesh = Mesh(rotated, [Cell(5,3,1,(0,1,2,3))])
    v, _, ids, _, _ = volume_data(mesh)
    polygons, _, _ = section(mesh, ids, v, axis=2, value=.5)
    assert len(polygons[0]) == 4
    cube = np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0],
                     [0,0,1],[1,0,1],[1,1,1],[0,1,1]])
    mesh = Mesh(cube, [Cell(8,3,1,tuple(range(8)))])
    v, q, ids, boundary, _ = volume_data(mesh)
    assert abs(v[0]-1) < 1e-15 and len(boundary) == 6
    assert abs(q[0]-0.7559526299369238) < 1e-14
    _, _, metric = section(mesh, ids, v, axis=1, value=.5)
    assert abs(metric['section_area']-1) < 1e-15
    joined = Mesh(np.vstack((points,[[0,0,-1]])),
                  [Cell(5,3,1,(0,1,2,3)),Cell(5,3,2,(0,2,1,4))])
    v, _, ids, _, _ = volume_data(joined)
    polygons, _, metric = section(joined,ids,v,axis=2,value=0)
    assert len(polygons) == 1 and metric['coplanar_shared_polygons_removed'] == 1
    assert abs(metric['section_area']-.5) < 1e-15
    midsides = np.array([(points[a]+points[b])*.5 for a,b in TET_EDGES])
    curved = np.vstack((points, midsides)); curved[4,2] += .01
    try:
        volume_data(Mesh(curved, [Cell(6,3,1,tuple(range(10)))]))
    except ValueError as error:
        assert 'curved Tet10' in str(error)
    else:
        raise AssertionError('curved Tet10 must not be described as an exact corner section')
    return {'passed':True, 'checks':['positive tet volume','analytic triangular section',
            'near-plane thin tetrahedron has no phantom section','thin crossing section retained',
            'four-edge quadrilateral section','planar Hex8 volume/section',
            'native-defined Hex8 SubTetMeanRatio','shared-face section deduplication',
            'curved Tet10 refused']}
