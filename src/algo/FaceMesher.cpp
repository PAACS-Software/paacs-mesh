// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Meshers.hpp"
#include <cstdint>
#include "paacs/algo/Delaunay2D.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace paacs::algo {

using mesh::CellType;
using mesh::Classification;

// Surface meshing strategy: work in *scaled* parametric coordinates
// (u*su, v*sv) with (su, sv) = GFace::metricScale(), an approximate isometry
// for planes and cylinders. Delaunay shape quality and spacing decisions in
// that metric then transfer to 3D. (General NURBS need the full first
// fundamental form - roadmap item for the OCCT backend.)
void meshFaces(const geo::Model& model, const sizing::MeshSpec& spec,
               mesh::Mesh& out) {
    for (const geo::GFace* f : model.faces()) {
        const Vec2 sc = f->metricScale();

        // 1. Boundary nodes: mesh nodes classified on this face's boundary
        //    edges and their end vertices.
        std::unordered_set<int> loopEdges, loopVerts;
        for (int et : f->boundaryEdgeTags()) {
            loopEdges.insert(et);
            if (const geo::GEdge* e = model.edge(et)) {
                loopVerts.insert(e->startVertexTag());
                loopVerts.insert(e->endVertexTag());
            }
        }

        std::vector<Index> meshIds;   // mesh node id per local point
        std::vector<Vec2>  uv;        // raw parametric coordinates
        std::vector<Vec2>  w;         // scaled (metric) coordinates
        auto pushPoint = [&](Index id, const Vec2& q) {
            meshIds.push_back(id);
            uv.push_back(q);
            w.push_back({q.u * sc.u, q.v * sc.v});
        };

        for (Index i = 0; i < out.nodeCount(); ++i) {
            const auto& c = out.nodeClass(i);
            const bool onLoop =
                (c.dim == 1 && loopEdges.count(c.tag)) ||
                (c.dim == 0 && loopVerts.count(c.tag));
            if (!onLoop) continue;
            pushPoint(i, f->paramAt(out.node(i)));
        }
        if (meshIds.size() < 3) continue;

        // 2. Interior lattice: staggered rows in metric space, curvature-aware
        //    spacing, rejected inside trimmed regions or within 0.6*h of any
        //    existing point.
        Vec2 lo, hi;
        f->paramBounds(lo, hi);
        const Vec2 mlo{lo.u * sc.u, lo.v * sc.v};
        const Vec2 mhi{hi.u * sc.u, hi.v * sc.v};

        auto hAtParam = [&](const Vec2& q) {
            return spec.hGeom(f->pointAt(q), f->curvatureAt(q));
        };
        const double hRef = hAtParam({(lo.u + hi.u) / 2, (lo.v + hi.v) / 2});

        int row = 0;
        long jrow = 0;
        for (double mv = mlo.v + hRef; mv < mhi.v - 0.4 * hRef;
             mv += hRef * 0.866) {
            const double shift = (row++ % 2) ? hRef * 0.5 : 0.0;
            long jcol = 0;
            for (double mu = mlo.u + hRef * 0.7 + shift;
                 mu < mhi.u - 0.4 * hRef; mu += hRef) {
                // Structured lattices on curved surfaces create near-coplanar
                // point patterns => Delaunay slivers. Deterministic hash
                // jitter breaks the regularity (planar faces stay exact).
                double ju = 0, jv = 0;
                const double kap = f->curvatureAt({mu / sc.u, mv / sc.v});
                if (kap > 1e-12) {
                    const auto hsh = [](long a, long b) {
                        std::uint64_t x = (std::uint64_t)a * 2654435761ULL ^
                                          (std::uint64_t)b * 40503ULL;
                        x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33;
                        return (double)(x & 0xffff) / 65535.0 - 0.5;
                    };
                    ju = 0.35 * hRef * hsh(jrow, jcol);
                    jv = 0.35 * hRef * hsh(jcol + 7919, jrow);
                }
                ++jcol;
                const Vec2 q{(mu + ju) / sc.u, (mv + jv) / sc.v};
                if (!f->containsParam(q)) continue;
                const double hq = hAtParam(q);
                const Vec2 mq{mu + ju, mv + jv};
                bool tooClose = false;
                for (size_t k = 0; k < w.size() && !tooClose; ++k)
                    tooClose = (w[k] - mq).norm() < 0.6 * hq;
                if (tooClose) continue;
                pushPoint(out.addNodeUnique(f->pointAt(q),
                                            Classification{2, f->tag()}),
                          q);
            }
            ++jrow;
        }

        // 3. Delaunay in the metric plane (near-isometric => good 3D shape).
        //    Prototype: convex parametric loops; concave loops need
        //    constrained edge recovery (next milestone).
        auto tris = delaunay2D(w);

        // 4. Filter triangles whose centroid leaves the trimmed region.
        for (const auto& t : tris) {
            const Vec2 cen{
                (uv[size_t(t.v[0])].u + uv[size_t(t.v[1])].u + uv[size_t(t.v[2])].u) / 3,
                (uv[size_t(t.v[0])].v + uv[size_t(t.v[1])].v + uv[size_t(t.v[2])].v) / 3};
            if (cen.u < lo.u || cen.u > hi.u || cen.v < lo.v || cen.v > hi.v)
                continue;
            if (!f->containsParam(cen)) continue;
            out.addCell(CellType::Tri3,
                        {meshIds[size_t(t.v[0])], meshIds[size_t(t.v[1])],
                         meshIds[size_t(t.v[2])]},
                        Classification{2, f->tag()});
        }
    }
}

} // namespace paacs::algo
