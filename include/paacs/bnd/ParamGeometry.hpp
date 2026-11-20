// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// ============================================================================
// BND-1 STAGE A2 - source-embedding geometric preflight (validation only).
// Design authority: docs/PMBND001_DESIGN.md sec 12 (A2 slice) + sec 5a/5b.
//
// A2 proves the SOURCE Tri3 embedding is geometrically sane in the
// AUTHORITATIVE parameter coordinates before any candidate generation:
//   - nondegeneracy      (no zero-area / collinear triangle)
//   - consistent orientation (all triangles same sign)
//   - proper crossings   (no two edges cross transversally)
//   - T-junctions        (no vertex lies in an edge's open interior)
// All tests are EXACT: predicates::orient2d in stored coordinates, never
// re-derived paramAt (sec 5b), never workspace Vec3 (sec 5a). This is a
// separate concern from the combinatorial preflight, so it lives in its
// own header; SegmentRecovery.hpp stays combinatorial. NO MUTATION PATH,
// NO CANDIDATE GENERATION - both remain NO-GO until A2 closes.
//
// A2 REQUIRES a workspace that already passed the A1 combinatorial
// preflight and manifold audit (validated Tri3, single-chart cells,
// total authoritative parameter map). It reads only const views + the
// ConstraintSet, and binds its result to the complete OwnershipStamp.
// ============================================================================
#include <algorithm>
#include <array>
#include <optional>
#include <vector>

#include "paacs/bnd/SegmentRecovery.hpp"
#include "paacs/core/Types.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/topo/Workspace.hpp"

namespace paacs::bnd {

enum class EmbeddingStatus {
    Ok,
    StaleStamp,             // ConstraintSet stamp no longer matches ws
    MissingParameter,       // an active node lacks an authoritative coord
    DegenerateTriangle,     // zero signed area (collinear) in param space
    InconsistentOrientation,// triangles do not all share one sign
    EdgeCrossing,           // two edges cross transversally
    TJunction,              // a vertex lies in the open interior of an edge
    TriangleOverlap,        // a triangle interior overlaps another's
};

struct EmbeddingIssue {
    EmbeddingStatus status{};
    // Up to two cells and/or two nodes describing the offending config;
    // optional so "not applicable" is never a zero sentinel.
    std::optional<topo::CellHandle> cellA, cellB;
    std::optional<topo::NodeHandle> nodeA, nodeB;
};

namespace detail {

// An undirected edge in the source mesh, carried with its two endpoints'
// authoritative coordinates so every test is exact and re-derivation
// free. Endpoints ordered by full-handle operator< for canonical keys.
struct ParamEdge {
    topo::NodeHandle a, b;
    Vec2 pa, pb;
    ChartId chart = 0;
};

// Exact "does point q lie on the CLOSED segment [pa,pb]" in stored
// coordinates: collinear (orient2d == Zero) AND within the exact
// bounding box. Endpoints included.
inline bool onClosedSegment(const Vec2& pa, const Vec2& pb, const Vec2& q) {
    if (predicates::orient2d(pa, pb, q) != predicates::Sign::Zero)
        return false;
    const double lo = std::min(pa.u, pb.u), hi = std::max(pa.u, pb.u);
    const double lov = std::min(pa.v, pb.v), hiv = std::max(pa.v, pb.v);
    return q.u >= lo && q.u <= hi && q.v >= lov && q.v <= hiv;
}

inline bool sameEndpoint(const topo::NodeHandle& x, const topo::NodeHandle& y) {
    return x == y;   // full-handle identity
}

// Exact "is q STRICTLY inside triangle (a,b,c)" in stored coordinates:
// all three orient2d signs equal and none Zero. Orientation-agnostic
// (works for CW or CCW triangles) by requiring the three signs to agree.
inline bool strictlyInside(const Vec2& a, const Vec2& b, const Vec2& c,
                           const Vec2& q) {
    const auto s1 = predicates::orient2d(a, b, q);
    const auto s2 = predicates::orient2d(b, c, q);
    const auto s3 = predicates::orient2d(c, a, q);
    if (s1 == predicates::Sign::Zero || s2 == predicates::Sign::Zero ||
        s3 == predicates::Sign::Zero)
        return false;                 // on an edge/vertex: not STRICT
    return s1 == s2 && s2 == s3;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// The A2 preflight. Returns every issue found (deterministically ordered),
// or empty for a clean embedding. Conservative and total: it examines all
// active Tri3 cells and the unique undirected edge set derived from them.
// ---------------------------------------------------------------------------
inline std::vector<EmbeddingIssue> preflightSourceEmbedding(
    const topo::TopologyWorkspace& ws, const ConstraintSet& cs) {
    std::vector<EmbeddingIssue> out;
    if (!cs.stamp.matches(ws)) {
        out.push_back({EmbeddingStatus::StaleStamp, {}, {}, {}, {}});
        return out;
    }
    const auto cells = ws.cells().activeBySequence();

    // Gather per-cell authoritative coordinates; a missing coord is fatal
    // to the geometric tests (A1 totality should already guarantee it, but
    // A2 must not read an absent coordinate).
    struct CellGeom {
        topo::CellHandle h;
        std::array<topo::NodeHandle, 3> n;
        std::array<Vec2, 3> p;
        ChartId chart = 0;              // A1 guarantees a cell is single-chart
    };
    std::vector<CellGeom> geom;
    geom.reserve(cells.size());
    for (const auto& ch : cells) {
        const auto cd = ws.cells().get(ch);
        CellGeom g;
        g.h = ch;
        bool ok = true;
        for (int i = 0; i < 3; ++i) {
            g.n[std::size_t(i)] = cd.nodes[std::size_t(i)];
            const ParamPoint* pp = cs.paramOf(cd.nodes[std::size_t(i)]);
            if (!pp) {
                out.push_back({EmbeddingStatus::MissingParameter, ch, {},
                               cd.nodes[std::size_t(i)], {}});
                ok = false;
                break;
            }
            g.p[std::size_t(i)] = pp->p;
            if (i == 0) g.chart = pp->chart;   // A1: all three agree
        }
        if (ok) geom.push_back(g);
    }
    if (!out.empty()) return out;   // cannot run geometry without coords

    // 1. Nondegeneracy + 2. consistent orientation, PER CHART (sec 5b:
    //    cross-chart coordinate comparison is meaningless).
    std::unordered_map<ChartId, predicates::Sign> ref;
    for (const auto& g : geom) {
        const auto s = predicates::orient2d(g.p[0], g.p[1], g.p[2]);
        if (s == predicates::Sign::Zero) {
            out.push_back({EmbeddingStatus::DegenerateTriangle, g.h, {},
                           {}, {}});
            continue;
        }
        const auto it = ref.find(g.chart);
        if (it == ref.end())
            ref.emplace(g.chart, s);
        else if (s != it->second)
            out.push_back({EmbeddingStatus::InconsistentOrientation, g.h,
                           {}, {}, {}});
    }

    // Build the unique undirected edge set (canonical full-handle key).
    struct EdgeKey {
        topo::NodeHandle lo, hi;
        bool operator==(const EdgeKey& o) const {
            return lo == o.lo && hi == o.hi;
        }
    };
    struct EdgeKeyHash {
        std::size_t operator()(const EdgeKey& k) const {
            const std::hash<topo::NodeHandle> hh;
            return hh(k.lo) * 1000003u + hh(k.hi);
        }
    };
    std::unordered_map<EdgeKey, detail::ParamEdge, EdgeKeyHash> edges;
    for (const auto& g : geom) {
        for (int i = 0; i < 3; ++i) {
            topo::NodeHandle x = g.n[std::size_t(i)],
                             y = g.n[std::size_t((i + 1) % 3)];
            Vec2 px = g.p[std::size_t(i)], py = g.p[std::size_t((i + 1) % 3)];
            if (y < x) {
                std::swap(x, y);
                std::swap(px, py);
            }
            edges.insert(
                {{x, y}, detail::ParamEdge{x, y, px, py, g.chart}});
        }
    }
    std::vector<detail::ParamEdge> E;
    E.reserve(edges.size());
    for (const auto& [k, e] : edges) E.push_back(e);
    // deterministic order: by canonical endpoints
    std::sort(E.begin(), E.end(),
              [](const detail::ParamEdge& x, const detail::ParamEdge& y) {
                  if (!(x.a == y.a)) return x.a < y.a;
                  return x.b < y.b;
              });

    // Collect the unique vertex set with coordinates (for T-junctions).
    struct Vtx {
        topo::NodeHandle h;
        Vec2 p;
        ChartId chart;
    };
    std::vector<Vtx> V;
    {
        std::unordered_set<topo::NodeHandle> seen;
        for (const auto& g : geom)
            for (int i = 0; i < 3; ++i)
                if (seen.insert(g.n[std::size_t(i)]).second)
                    V.push_back({g.n[std::size_t(i)], g.p[std::size_t(i)],
                                 g.chart});
        std::sort(V.begin(), V.end(),
                  [](const Vtx& x, const Vtx& y) { return x.h < y.h; });
    }

    // 3. Proper crossings: every unordered pair of edges that do NOT share
    //    an endpoint. Transversal crossing <=> each segment strictly
    //    separates the other's endpoints (strict orientation signs).
    for (std::size_t i = 0; i < E.size(); ++i) {
        for (std::size_t j = i + 1; j < E.size(); ++j) {
            const auto& e1 = E[i];
            const auto& e2 = E[j];
            if (e1.chart != e2.chart) continue;   // sec 5b: chart-local
            if (detail::sameEndpoint(e1.a, e2.a) ||
                detail::sameEndpoint(e1.a, e2.b) ||
                detail::sameEndpoint(e1.b, e2.a) ||
                detail::sameEndpoint(e1.b, e2.b))
                continue;   // shared vertex: not a proper crossing
            const auto d1 = predicates::orient2d(e1.pa, e1.pb, e2.pa);
            const auto d2 = predicates::orient2d(e1.pa, e1.pb, e2.pb);
            const auto d3 = predicates::orient2d(e2.pa, e2.pb, e1.pa);
            const auto d4 = predicates::orient2d(e2.pa, e2.pb, e1.pb);
            const bool straddle1 = d1 != d2 && d1 != predicates::Sign::Zero &&
                                   d2 != predicates::Sign::Zero;
            const bool straddle2 = d3 != d4 && d3 != predicates::Sign::Zero &&
                                   d4 != predicates::Sign::Zero;
            if (straddle1 && straddle2)
                out.push_back({EmbeddingStatus::EdgeCrossing, {}, {},
                               e1.a, e2.a});
        }
    }

    // 4. T-junctions: a vertex strictly interior to some edge it is not an
    //    endpoint of (on the closed segment but not equal to either end).
    for (const auto& v : V) {
        for (const auto& e : E) {
            if (v.chart != e.chart) continue;     // sec 5b: chart-local
            if (detail::sameEndpoint(v.h, e.a) ||
                detail::sameEndpoint(v.h, e.b))
                continue;
            if (detail::onClosedSegment(e.pa, e.pb, v.p))
                out.push_back({EmbeddingStatus::TJunction, {}, {}, v.h, e.a});
        }
    }

    // 5. Triangle-interior overlap: for each unordered pair of triangles
    //    that share NO vertex, reject if any vertex of one lies strictly
    //    inside the other (a nested or overlapping triangle - invalid
    //    embedding the edge/vertex tests miss). Deterministic attribution:
    //    the cell whose vertex is inside is cellA, the container cellB.
    for (std::size_t i = 0; i < geom.size(); ++i) {
        for (std::size_t j = 0; j < geom.size(); ++j) {
            if (i == j) continue;
            const auto& gi = geom[i];   // vertex donor
            const auto& gj = geom[j];   // container
            if (gi.chart != gj.chart) continue;   // sec 5b: chart-local
            for (int a = 0; a < 3; ++a) {
                // Skip only a donor vertex SHARED with the container - a
                // shared vertex is legally on the container's boundary; a
                // NON-shared donor vertex strictly inside is the defect
                // (this catches same-side nesting on a shared edge, which
                // skipping the whole pair would miss).
                bool shared = false;
                for (int b = 0; b < 3; ++b)
                    if (gi.n[std::size_t(a)] == gj.n[std::size_t(b)]) {
                        shared = true;
                        break;
                    }
                if (shared) continue;
                if (detail::strictlyInside(gj.p[0], gj.p[1], gj.p[2],
                                           gi.p[std::size_t(a)])) {
                    out.push_back({EmbeddingStatus::TriangleOverlap, gi.h,
                                   gj.h, gi.n[std::size_t(a)], {}});
                    break;              // one witness per (donor,container)
                }
            }
        }
    }

    // Deterministic total order over the whole issue vector.
    std::stable_sort(
        out.begin(), out.end(),
        [](const EmbeddingIssue& a, const EmbeddingIssue& b) {
            if (a.status != b.status)
                return static_cast<int>(a.status) < static_cast<int>(b.status);
            auto cmpCell = [](const std::optional<topo::CellHandle>& x,
                              const std::optional<topo::CellHandle>& y) {
                if (x.has_value() != y.has_value())
                    return int(x.has_value()) - int(y.has_value());
                if (x && !(*x == *y)) return *x < *y ? -1 : 1;
                return 0;
            };
            auto cmpNode = [](const std::optional<topo::NodeHandle>& x,
                              const std::optional<topo::NodeHandle>& y) {
                if (x.has_value() != y.has_value())
                    return int(x.has_value()) - int(y.has_value());
                if (x && !(*x == *y)) return *x < *y ? -1 : 1;
                return 0;
            };
            if (int c = cmpCell(a.cellA, b.cellA)) return c < 0;
            if (int c = cmpCell(a.cellB, b.cellB)) return c < 0;
            if (int c = cmpNode(a.nodeA, b.nodeA)) return c < 0;
            if (int c = cmpNode(a.nodeB, b.nodeB)) return c < 0;
            return false;
        });
    return out;
}

}  // namespace paacs::bnd
