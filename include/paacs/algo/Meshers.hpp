// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// The 1D -> 2D -> 3D meshing pipeline. Each stage consumes the geometry
// abstraction + sizing spec and appends classified nodes/cells to the Mesh.
//
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/sizing/SizeField.hpp"

#include <set>

namespace paacs::algo {

struct TetHandoffMetrics;

// 1D: discretize every model edge into Line2 segments sized by spec.h().
// Node spacing adapts to the sizing field sampled along the curve.
void meshEdges(const geo::Model& model, const sizing::MeshSpec& spec,
               mesh::Mesh& out);

// 2D: mesh every face with Tri3 via Delaunay in the parametric domain,
// reusing the 1D boundary nodes so the surface mesh is conforming.
void meshFaces(const geo::Model& model, const sizing::MeshSpec& spec,
               mesh::Mesh& out);

// 3D: tetrahedralize every solid (Tet4), reusing surface nodes and inserting
// interior points on a sizing-driven lattice. Prototype note: interior
// Delaunay does not yet do constrained boundary recovery on concave solids.
void meshVolumes(const geo::Model& model, const sizing::MeshSpec& spec,
                 mesh::Mesh& out);

// CK-3 constrained-v2 volume entry. The accepted legacy meshVolumes() path
// remains unchanged; this additive route validates a canonical Tet candidate,
// imports the retained subset into private authoritative topology, releases
// transient candidate storage, then appends the verified cells.
void meshVolumesAuthoritative(const geo::Model& model,
                              const sizing::MeshSpec& spec,
                              mesh::Mesh& out,
                              TetHandoffMetrics& metrics);

// Convenience: full pipeline.
// Face tags whose geometry is curved (nonzero curvature bound sampled over
// the parameter domain). Peeling stopgaps are restricted to tets touching
// these faces: chordal hull caps exist only on curved boundaries, while on
// planar solids "low quality all-boundary" tets can be the legitimate
// domain (thin plates) and must never be removed.
inline std::set<int> curvedFaceTags(const geo::Model& model) {
    std::set<int> tags;
    for (const geo::GFace* f : model.faces()) {
        Vec2 lo, hi;
        f->paramBounds(lo, hi);
        const Vec2 samples[5] = {
            {(lo.u + hi.u) / 2, (lo.v + hi.v) / 2},
            {lo.u + 0.25 * (hi.u - lo.u), lo.v + 0.25 * (hi.v - lo.v)},
            {lo.u + 0.75 * (hi.u - lo.u), lo.v + 0.25 * (hi.v - lo.v)},
            {lo.u + 0.25 * (hi.u - lo.u), lo.v + 0.75 * (hi.v - lo.v)},
            {lo.u + 0.75 * (hi.u - lo.u), lo.v + 0.75 * (hi.v - lo.v)}};
        for (const Vec2& uv : samples)
            if (f->curvatureAt(uv) > 0) { tags.insert(f->tag()); break; }
    }
    return tags;
}

inline void meshModel(const geo::Model& model, const sizing::MeshSpec& spec,
                      mesh::Mesh& out) {
    core::setThreadCount(spec.threads);
    // Scale-aware node coincidence tolerance (PM-NUM-002): absolute 1e-9
    // collapsed models smaller than ~1e-9 to a single node. Derive from the
    // model bounding box instead. (Full Tolerances policy is milestone M3.)
    {
        Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
        for (const geo::GSolid* s : model.solids()) {
            Vec3 a, b;
            s->bbox(a, b);
            lo.x = std::min(lo.x, a.x); hi.x = std::max(hi.x, b.x);
            lo.y = std::min(lo.y, a.y); hi.y = std::max(hi.y, b.y);
            lo.z = std::min(lo.z, a.z); hi.z = std::max(hi.z, b.z);
        }
        const double diag = (hi - lo).norm();
        if (diag > 0 && diag < 1e300 && out.dedupTolerance() < 0)
            out.setDedupTolerance(1e-9 * diag);
    }
    meshEdges(model, spec, out);
    meshFaces(model, spec, out);
    meshVolumes(model, spec, out);
}

// Order elevation: Tet4 -> Tet10, Tri3 -> Tri6, Line2 -> Line3, deduplicating
// shared midside nodes. Midside nodes on classified boundary entities are
// snapped back to the true geometry (validity-aware curving hook lives here).
struct ElevateResult {
    Index midsides = 0;            // midside nodes created
    Index snapped = 0;             // midsides projected onto CAD geometry
    Index straightened = 0;        // midsides blended back for validity
    Index invalidAfter = 0;        // tet10 still invalid (should be 0)
};
// Tet4 -> Tet10 (and Tri3 -> Tri6, Line2 -> Line3). Midsides classified via
// the owning-cell fallback rule and snapped to CAD geometry; elements whose
// sampled Jacobian goes non-positive from snapping are straightened back
// just enough to restore validity (geometric fidelity yields to validity -
// proper high-order untangling is the roadmap upgrade).
ElevateResult elevateToQuadratic(const geo::Model& model, mesh::Mesh& m);

} // namespace paacs::algo

namespace paacs::io {
// Legacy VTK writer (unstructured grid) - loads in ParaView for inspection.
void writeVtk(const mesh::Mesh& m, const std::string& path,
              bool volumeOnly = false);
// MSH 4.1 ASCII with entity blocks from mesh classification.
void writeMsh(const mesh::Mesh& m, const std::string& path);
} // namespace paacs::io
