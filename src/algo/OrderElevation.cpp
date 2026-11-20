// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Meshers.hpp"
#include "paacs/mesh/Quality.hpp"

#include <map>
#include <utility>

namespace paacs::algo {

using mesh::CellType;
using mesh::Classification;
using mesh::Mesh;

namespace {

// Midside classification. If both endpoints live on the same entity the
// midside belongs to it. Otherwise the ONLY safe answer is the entity of the
// cell that owns this edge: two endpoints on different CAD entities say
// nothing about where the edge between them lives (e.g. a face-interior mesh
// edge can join nodes classified on two different boundary curves - its
// midside belongs to the FACE, not to either curve). Choosing "the higher-dim
// endpoint" here used to snap such midsides onto the wrong entity, moving
// them by up to ~0.16 on a unit box. Because cells are processed in
// dimension order (Line3 before Tri6 before Tet10) and midsides are cached
// per node pair, true curve midsides are created by the Line3 pass and are
// already in the cache when faces and solids ask for them.
Classification midClass(const Classification& a, const Classification& b,
                        const Classification& owningCell) {
    if (a.dim == b.dim && a.tag == b.tag) return a;
    return owningCell;
}

} // namespace

ElevateResult elevateToQuadratic(const geo::Model& model, Mesh& m) {
    ElevateResult res;
    std::map<std::pair<Index, Index>, Index> midOf;

    auto midNode = [&](Index a, Index b,
                       const Classification& owningCell) -> Index {
        auto key = a < b ? std::make_pair(a, b) : std::make_pair(b, a);
        auto it = midOf.find(key);
        if (it != midOf.end()) return it->second;

        Classification cls = midClass(m.nodeClass(a), m.nodeClass(b),
                                      owningCell);
        Vec3 p = (m.node(a) + m.node(b)) * 0.5;

        // Geometry snapping: midside nodes classified on a curved edge/face
        // are pulled onto the true geometry. This is the hook where
        // validity-aware curving (Bezier Jacobian bounds) will attach.
        if (cls.dim == 2) {
            if (const geo::GFace* f = model.face(cls.tag))
                p = f->pointAt(f->paramAt(p));
        } else if (cls.dim == 1) {
            if (const geo::GEdge* e = model.edge(cls.tag))
                p = e->closestPoint(p);
        }
        Index id = m.addNode(p, cls);
        midOf[key] = id;
        ++res.midsides;
        if ((p - (m.node(a) + m.node(b)) * 0.5).norm() > 1e-14) ++res.snapped;
        return id;
    };

    const Index nCells = m.cellCount();
    Mesh elevated;

    // Copy all original nodes 1:1 (ids preserved).
    for (Index i = 0; i < m.nodeCount(); ++i)
        elevated.addNode(m.node(i), m.nodeClass(i));

    // Midside creation must happen in `m`'s id space; simplest correct
    // prototype approach: create midside nodes in m, then rebuild cells.
    struct NewCell { CellType t; std::vector<Index> conn; Classification c; };
    std::vector<NewCell> newCells;
    newCells.reserve(static_cast<size_t>(nCells));

    for (Index ci = 0; ci < nCells; ++ci) {
        const auto& c = m.cell(ci);
        const Index* v = m.cellNodes(ci);
        switch (c.type) {
        case CellType::Line2:
            newCells.push_back({CellType::Line3,
                {v[0], v[1], midNode(v[0], v[1], c.onEntity)}, c.onEntity});
            break;
        case CellType::Tri3:
            newCells.push_back({CellType::Tri6,
                {v[0], v[1], v[2],
                 midNode(v[0], v[1], c.onEntity),
                 midNode(v[1], v[2], c.onEntity),
                 midNode(v[2], v[0], c.onEntity)},
                c.onEntity});
            break;
        case CellType::Tet4:
            // MSH/VTK tet10 ordering: 4 corners then edges
            // 01, 12, 20, 03, 13, 23.
            newCells.push_back({CellType::Tet10,
                {v[0], v[1], v[2], v[3],
                 midNode(v[0], v[1], c.onEntity),
                 midNode(v[1], v[2], c.onEntity),
                 midNode(v[2], v[0], c.onEntity),
                 midNode(v[0], v[3], c.onEntity),
                 midNode(v[1], v[3], c.onEntity),
                 midNode(v[2], v[3], c.onEntity)},
                c.onEntity});
            break;
        default:
            newCells.push_back({c.type,
                std::vector<Index>(v, v + mesh::nodesPerCell(c.type)),
                c.onEntity});
        }
    }

    // Copy midside nodes appended to m into `elevated`, then swap in cells.
    for (Index i = elevated.nodeCount(); i < m.nodeCount(); ++i)
        elevated.addNode(m.node(i), m.nodeClass(i));
    for (const auto& nc : newCells)
        elevated.addCell(nc.t, nc.conn, nc.c);

    m = std::move(elevated);

    // Validity pass: snapping midsides onto strongly curved geometry can
    // locally invert the isoparametric map even though every corner tet is
    // positive (measured: 6 inverted tet10 on a coarse cylinder). Blend the
    // offending midsides toward their linear midpoints - halving the
    // curvature each round monotonically approaches the straight (valid)
    // element - until all sampled Jacobians are positive. Deterministic:
    // fixed cell order, fixed blend factor.
    for (int round = 0; round < 8; ++round) {
        Index fixedThisRound = 0;
        for (Index ci = 0; ci < m.cellCount(); ++ci) {
            if (m.cell(ci).type != CellType::Tet10) continue;
            const Index* v = m.cellNodes(ci);
            Vec3 p[10];
            for (int k = 0; k < 10; ++k) p[k] = m.node(v[k]);
            if (mesh::tet10MinJacobianRatio(p) > 1e-6) continue;
            static const int ed[6][2] = {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
            for (int e = 0; e < 6; ++e) {
                const Vec3 lin =
                    (m.node(v[ed[e][0]]) + m.node(v[ed[e][1]])) * 0.5;
                const Vec3 cur = m.node(v[4 + e]);
                if ((cur - lin).norm2() < 1e-28) continue;
                m.setNode(v[4 + e], cur + (lin - cur) * 0.5);
                ++res.straightened;
            }
            ++fixedThisRound;
        }
        if (!fixedThisRound) break;
    }
    res.invalidAfter = mesh::countInvalidTet10(m, 0.0);
    return res;
}

} // namespace paacs::algo
