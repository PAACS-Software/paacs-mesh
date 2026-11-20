// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// ============================================================================
// BND-1 - constrained SEGMENT recovery (2D, Tri3, existing endpoints).
// Design authority: docs/PMBND001_DESIGN.md (amended section 12 scope).
//
// STAGE A1 ONLY (foundation slice): value types, ownership stamping,
// whole-set preflight, the BND-owned VALIDATED temporary Tri3 adjacency
// index, and segment-chain validation/discovery. NO MUTATION PATH EXISTS
// IN THIS FILE. The source-embedding geometric preflight (nondegeneracy,
// orientation, crossings, T-junctions) is Stage A2; candidate generation
// and the cavity commit are later stages, mutating exclusively through
// the frozen M4-B.8 cavity primitive.
//
// IDENTITY RULE (Stage-A review): FULL handles everywhere - every map is
// keyed by the complete NodeHandle/CellHandle (slot, generation, pool),
// via the existing std::hash<Handle> specialization; reducing a handle
// to slot+generation would let a foreign handle with matching fields
// alias local data. The adjacency index binds to the COMPLETE ownership
// stamp, not the revision alone.
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "paacs/core/Types.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/topo/Workspace.hpp"

namespace paacs::bnd {

using ConstraintVertexId = std::uint64_t;
using ConstraintEntityId = std::uint64_t;
using ChartId = std::uint64_t;

// Authoritative stored parameter coordinates (sec 5b): one entry per node,
// never re-derived; the owning chart disambiguates seams/periodic faces.
struct ParamPoint {
    Vec2 p{};
    ChartId chart = 0;
};

struct ConstraintSegment {
    ConstraintEntityId id = 0;
    ConstraintVertexId a = 0, b = 0;
    ChartId chart = 0;
    mesh::Classification cls{};        // applied to recovered edges:
                                       // dim == 1 and tag >= 0 (sec 12)
};

// Concrete ownership stamp (sec 12): the ACTUAL public identities - all
// four pool ids plus the topology revision.
struct OwnershipStamp {
    std::uint64_t nodePool = 0, cellPool = 0, edgePool = 0, facePool = 0;
    std::uint64_t revision = 0;
    static OwnershipStamp of(const topo::TopologyWorkspace& ws) {
        return {ws.nodes().poolId(), ws.cells().poolId(),
                ws.edges().poolId(), ws.faces().poolId(),
                ws.topologyRevision()};
    }
    bool matches(const topo::TopologyWorkspace& ws) const {
        const auto now = of(ws);
        return nodePool == now.nodePool && cellPool == now.cellPool &&
               edgePool == now.edgePool && facePool == now.facePool &&
               revision == now.revision;
    }
};

// BND-2A additive (design §5.2): the persistent AUTHORITATIVE symbolic
// parameter of a CREATED Steiner node. A flat POD (no variant) so it can live
// in the BND-1 ConstraintSet without a forward-dependency on the BND-2 symbolic
// header. `affine == false` => an existing/stored node at `storedUv`; otherwise
// t = num/2^level on the canonical undirected endpoints (lo,hi) by
// ConstraintVertexId. BND-1 never reads this map, so its behaviour is unchanged.
struct AuthoritativeParamRecord {
    bool affine = false;
    Vec2 storedUv{};
    ConstraintVertexId lo = 0, hi = 0;
    std::uint64_t num = 0;
    std::uint32_t level = 0;
    ChartId chart = 0;
};

struct ConstraintSet {
    std::unordered_map<ConstraintVertexId, topo::NodeHandle> vertices;
    std::unordered_map<topo::NodeHandle, ParamPoint> params;  // FULL key
    std::vector<ConstraintSegment> segments;
    OwnershipStamp stamp;
    // BND-2A additive persistent authoritative-parameter store for created
    // Steiner nodes (design §5.2); populated atomically by the recovery route.
    std::unordered_map<topo::NodeHandle, AuthoritativeParamRecord> authParams;
    const ParamPoint* paramOf(const topo::NodeHandle& h) const {
        const auto it = params.find(h);
        return it == params.end() ? nullptr : &it->second;
    }
    const AuthoritativeParamRecord* authParamOf(const topo::NodeHandle& h) const {
        const auto it = authParams.find(h);
        return it == authParams.end() ? nullptr : &it->second;
    }
};

// ---------------------------------------------------------------------------
// Shared manifold-Tri3 audit: used by BOTH the whole-set preflight and
// the adjacency-index builder. useCount <= 2 alone is NOT manifoldness -
// a pinched vertex has manifold edges with a disconnected node link (the
// same invariant M4-B.8 Stage A had to add).
// ---------------------------------------------------------------------------
enum class Tri3AuditResult {
    Ok,
    NonTri3Cell,
    BadEdgeCount,          // cellEdges() span size != 3
    OverusedEdge,          // more than two owners
    DuplicateTriangle,     // two live cells, same 3 nodes (any order)
    PinchedVertex,         // node link is not one simple path/cycle
};

inline Tri3AuditResult auditTri3Manifold(
    const topo::TopologyWorkspace& ws) {
    std::unordered_map<std::int32_t, int> owners;   // edge slot -> count
    std::unordered_map<topo::NodeHandle,
                       std::vector<std::array<topo::NodeHandle, 2>>>
        link;                                       // node -> opposite edges
    // Orientation-insensitive canonical triangle key: the three complete
    // handles sorted. Two coincident Tri3 cells (identical OR reversed
    // connectivity) form a multigraph the degree/connectivity audit would
    // otherwise misread as a valid cycle - they must be rejected FIRST,
    // and specifically, not incidentally (review).
    using TriKey = std::array<topo::NodeHandle, 3>;
    struct TriKeyHash {
        std::size_t operator()(const TriKey& k) const {
            std::size_t s = 0;
            const std::hash<topo::NodeHandle> hh;
            for (const auto& v : k)
                s = s * 1000003u + hh(v);   // full-handle hash, all 3 fields
            return s;
        }
    };
    std::unordered_set<TriKey, TriKeyHash> triKeys;   // full-handle equality
    for (const auto& ch : ws.cells().activeBySequence()) {
        const auto cd = ws.cells().get(ch);
        if (cd.type != topo::ElementType::Tri3)
            return Tri3AuditResult::NonTri3Cell;
        const auto ce = ws.edges().cellEdges(ws.cells(), ch);
        if (ce.size() != 3) return Tri3AuditResult::BadEdgeCount;
        TriKey key{cd.nodes[0], cd.nodes[1], cd.nodes[2]};
        std::sort(key.begin(), key.end());   // orientation-insensitive,
                                             // full-handle operator<
        if (!triKeys.insert(key).second)
            return Tri3AuditResult::DuplicateTriangle;
        for (const auto& eh : ce)
            if (++owners[eh.slot] > 2) return Tri3AuditResult::OverusedEdge;
        for (int i = 0; i < 3; ++i)
            link[cd.nodes[std::size_t(i)]].push_back(
                {cd.nodes[std::size_t((i + 1) % 3)],
                 cd.nodes[std::size_t((i + 2) % 3)]});
    }
    // Node-link audit: the opposite edges of a node's incident triangles
    // must form ONE simple path or cycle.
    for (const auto& [n, segs] : link) {
        std::unordered_map<topo::NodeHandle, int> deg;
        for (const auto& s : segs) {
            if (++deg[s[0]] > 2 || ++deg[s[1]] > 2)
                return Tri3AuditResult::PinchedVertex;
        }
        int odd = 0;
        for (const auto& [v, d] : deg)
            if (d == 1) ++odd;
        if (!(odd == 0 || odd == 2)) return Tri3AuditResult::PinchedVertex;
        // connectivity: walk from any segment; all segments reachable
        std::unordered_map<topo::NodeHandle,
                           std::vector<topo::NodeHandle>>
            adj;
        for (const auto& s : segs) {
            adj[s[0]].push_back(s[1]);
            adj[s[1]].push_back(s[0]);
        }
        std::unordered_set<topo::NodeHandle> seen;
        std::vector<topo::NodeHandle> stack{segs[0][0]};
        while (!stack.empty()) {
            const auto v = stack.back();
            stack.pop_back();
            if (!seen.insert(v).second) continue;
            for (const auto& w : adj[v]) stack.push_back(w);
        }
        if (seen.size() != adj.size()) return Tri3AuditResult::PinchedVertex;
    }
    return Tri3AuditResult::Ok;
}

// ---------------------------------------------------------------------------
// Whole-set preflight (sec 5 / sec 12): Rejected = zero revision change is
// structural here - Stage A1 has no mutation path at all.
// ---------------------------------------------------------------------------
enum class InvalidInputReason {
    StaleStamp,
    DuplicateEntityId,
    UnknownEndpoint,
    StaleHandle,
    ForeignHandle,
    NonfiniteCoordinate,
    MissingParameter,       // parameter-map TOTALITY over active nodes
    ChartMismatch,
    InvalidClassification,  // requires dim == 1 AND tag >= 0
    DegenerateSegment,
    NonTri3Cell,
    NonmanifoldIncidence,   // bad edge count or > 2 owners
    DuplicateTriangle,      // coincident live Tri3 cells
    PinchedVertex,
};

struct PreflightIssue {
    InvalidInputReason reason{};
    std::optional<ConstraintEntityId> entity;   // no zero sentinels
    std::optional<topo::NodeHandle> node;
};

inline std::vector<PreflightIssue> preflightConstraintSet(
    const topo::TopologyWorkspace& ws, const ConstraintSet& cs) {
    std::vector<PreflightIssue> out;
    if (!cs.stamp.matches(ws)) {
        out.push_back({InvalidInputReason::StaleStamp, {}, {}});
        return out;                  // nothing else is trustworthy
    }
    // vertex association handles: owned and current
    for (const auto& [vid, h] : cs.vertices) {
        if (h.pool != ws.nodes().poolId())
            out.push_back({InvalidInputReason::ForeignHandle, {}, h});
        else if (!ws.nodes().isCurrent(h))
            out.push_back({InvalidInputReason::StaleHandle, {}, h});
    }
    // parameter-map ENTRY audit (sec 12 / review: extra, stale, foreign,
    // nonfinite entries are all faults, not just gaps)
    for (const auto& [h, pp] : cs.params) {
        if (h.pool != ws.nodes().poolId())
            out.push_back({InvalidInputReason::ForeignHandle, {}, h});
        else if (!ws.nodes().isCurrent(h))
            out.push_back({InvalidInputReason::StaleHandle, {}, h});
        if (!std::isfinite(pp.p.u) || !std::isfinite(pp.p.v))
            out.push_back({InvalidInputReason::NonfiniteCoordinate, {}, h});
    }
    // TOTALITY over every ACTIVE node (a candidate cavity may use any)
    for (const auto& nh : ws.nodes().activeBySequence())
        if (!cs.paramOf(nh))
            out.push_back({InvalidInputReason::MissingParameter, {}, nh});
    // chart consistency ACROSS candidate Tri3 cells (Stage-A1
    // conservative rule, sec 12: every active cell single-chart)
    for (const auto& ch : ws.cells().activeBySequence()) {
        const auto cd = ws.cells().get(ch);
        const ParamPoint* p0 = cs.paramOf(cd.nodes[0]);
        const ParamPoint* p1 = cs.paramOf(cd.nodes[1]);
        const ParamPoint* p2 = cs.paramOf(cd.nodes[2]);
        if (p0 && p1 && p2 &&
            !(p0->chart == p1->chart && p1->chart == p2->chart)) {
            out.push_back({InvalidInputReason::ChartMismatch, {}, cd.nodes[0]});
        }
    }
    // per-segment checks
    std::unordered_set<ConstraintEntityId> seen;
    for (const auto& s : cs.segments) {
        if (!seen.insert(s.id).second)
            out.push_back({InvalidInputReason::DuplicateEntityId, s.id, {}});
        const auto ia = cs.vertices.find(s.a);
        const auto ib = cs.vertices.find(s.b);
        if (ia == cs.vertices.end() || ib == cs.vertices.end()) {
            out.push_back({InvalidInputReason::UnknownEndpoint, s.id, {}});
            continue;
        }
        const ParamPoint* pa = cs.paramOf(ia->second);
        const ParamPoint* pb = cs.paramOf(ib->second);
        if (!pa || !pb) continue;    // already reported as MissingParameter
        if (pa->chart != s.chart || pb->chart != s.chart)
            out.push_back({InvalidInputReason::ChartMismatch, s.id, {}});
        if (s.a == s.b || (pa->p.u == pb->p.u && pa->p.v == pb->p.v))
            out.push_back({InvalidInputReason::DegenerateSegment, s.id, {}});
        if (s.cls.dim != 1 || s.cls.tag < 0)
            out.push_back({InvalidInputReason::InvalidClassification, s.id,
                           {}});
    }
    // manifold Tri3 audit (shared with the index builder)
    switch (auditTri3Manifold(ws)) {
        case Tri3AuditResult::Ok:
            break;
        case Tri3AuditResult::NonTri3Cell:
            out.push_back({InvalidInputReason::NonTri3Cell, {}, {}});
            break;
        case Tri3AuditResult::BadEdgeCount:
        case Tri3AuditResult::OverusedEdge:
            out.push_back({InvalidInputReason::NonmanifoldIncidence, {}, {}});
            break;
        case Tri3AuditResult::DuplicateTriangle:
            out.push_back({InvalidInputReason::DuplicateTriangle, {}, {}});
            break;
        case Tri3AuditResult::PinchedVertex:
            out.push_back({InvalidInputReason::PinchedVertex, {}, {}});
            break;
    }
    // Deterministic diagnostics (review): vertices/params iterate in
    // unordered_map order, so canonicalize by a documented TOTAL key -
    // reason, then entity presence/value, then full node identity -
    // before returning. Two logically identical sets now yield
    // bit-identical issue vectors regardless of insertion order.
    std::stable_sort(
        out.begin(), out.end(),
        [](const PreflightIssue& a, const PreflightIssue& b) {
            if (a.reason != b.reason)
                return static_cast<int>(a.reason) <
                       static_cast<int>(b.reason);
            const bool ae = a.entity.has_value(), be = b.entity.has_value();
            if (ae != be) return ae < be;
            if (ae && *a.entity != *b.entity) return *a.entity < *b.entity;
            const bool an = a.node.has_value(), bn = b.node.has_value();
            if (an != bn) return an < bn;
            if (an && !(*a.node == *b.node))   // full-handle total order
                return *a.node < *b.node;
            return false;
        });
    return out;
}

// ---------------------------------------------------------------------------
// BND-owned VALIDATED temporary Tri3 adjacency index (sec 12). The face-
// based workspace adjacency is empty for Tri3; this index replaces it,
// binds to the COMPLETE ownership stamp, and REFUSES construction on any
// non-Tri3 cell, wrong edge count, three-owner edge, or pinched vertex.
// Rebuilt before each entity's candidate search (each successful prefix
// commit moves the revision, which invalidates the stamp).
// ---------------------------------------------------------------------------
struct TriAdjacencyIndex {
    OwnershipStamp stamp;            // COMPLETE binding, not revision alone
    // FULL CellHandle key -> the neighbor across each of its 3 edges,
    // ALIGNED with cellEdges() order; a default handle marks boundary.
    std::unordered_map<topo::CellHandle, std::array<topo::CellHandle, 3>>
        neighbors;
    // FULL NodeHandle key -> incident cells (the star)
    std::unordered_map<topo::NodeHandle, std::vector<topo::CellHandle>>
        star;
    bool current(const topo::TopologyWorkspace& ws) const {
        return stamp.matches(ws);
    }
};

inline std::optional<TriAdjacencyIndex> buildTriAdjacency(
    const topo::TopologyWorkspace& ws) {
    if (auditTri3Manifold(ws) != Tri3AuditResult::Ok) return std::nullopt;
    TriAdjacencyIndex idx;
    idx.stamp = OwnershipStamp::of(ws);
    struct Owner {
        topo::CellHandle c{};
        int local = -1;
    };
    std::unordered_map<std::int32_t, std::array<Owner, 2>> owners;
    for (const auto& ch : ws.cells().activeBySequence()) {
        idx.neighbors[ch] = {topo::CellHandle{}, topo::CellHandle{},
                             topo::CellHandle{}};
        const auto ce = ws.edges().cellEdges(ws.cells(), ch);
        int li = 0;
        for (const auto& eh : ce) {
            auto& o = owners[eh.slot];
            if (o[0].local < 0)
                o[0] = {ch, li};
            else
                o[1] = {ch, li};     // audit guarantees no third owner
            ++li;
        }
        const auto cd = ws.cells().get(ch);
        for (int i = 0; i < 3; ++i)
            idx.star[cd.nodes[std::size_t(i)]].push_back(ch);
    }
    for (const auto& [slot, o] : owners) {
        if (o[1].local < 0) continue;                 // boundary edge
        idx.neighbors[o[0].c][std::size_t(o[0].local)] = o[1].c;
        idx.neighbors[o[1].c][std::size_t(o[1].local)] = o[0].c;
    }
    return idx;
}

// ---------------------------------------------------------------------------
// Segment-chain validation (sec 1): exact tests in STORED constraint-
// domain coordinates - orient2d == Zero collinearity, exact bounds, and
// a strictly increasing deterministic dominant-axis parameter.
// ---------------------------------------------------------------------------
enum class ChainStatus {
    Valid,
    Empty,
    EndpointMismatch,
    DuplicateNode,
    MissingParameter,
    ChartMismatch,
    OffSegment,
    NonIncreasingParameter,
    MissingEdge,
    EdgeNotRequired,
    ClassificationMismatch,
};

struct ChainCheck {
    ChainStatus status = ChainStatus::Valid;
    std::size_t at = 0;              // failing chain index (node or edge)
};

namespace detail {
// Deterministic dominant axis of the segment direction: u unless
// |dv| > |du| (ties keep u - documented canonical rule, sec 7).
inline int dominantAxis(const Vec2& a, const Vec2& b) {
    const double du = std::fabs(b.u - a.u), dv = std::fabs(b.v - a.v);
    return dv > du ? 1 : 0;
}
inline double axisCoord(const Vec2& p, int axis) {
    return axis == 0 ? p.u : p.v;
}
}  // namespace detail

inline ChainCheck validateSegmentChain(
    const topo::TopologyWorkspace& ws, const ConstraintSet& cs,
    const ConstraintSegment& seg,
    const std::vector<topo::NodeHandle>& chain) {
    if (chain.size() < 2) return {ChainStatus::Empty, 0};
    const auto ia = cs.vertices.find(seg.a);
    const auto ib = cs.vertices.find(seg.b);
    if (ia == cs.vertices.end() || ib == cs.vertices.end() ||
        !(chain.front() == ia->second) || !(chain.back() == ib->second))
        return {ChainStatus::EndpointMismatch, 0};
    const ParamPoint* pa = cs.paramOf(ia->second);
    const ParamPoint* pb = cs.paramOf(ib->second);
    if (!pa || !pb) return {ChainStatus::MissingParameter, 0};
    const int axis = detail::dominantAxis(pa->p, pb->p);
    const double ta = detail::axisCoord(pa->p, axis);
    const double tb = detail::axisCoord(pb->p, axis);
    const double dir = tb > ta ? 1.0 : -1.0;   // nondegenerate (preflight)
    std::unordered_set<topo::NodeHandle> seen;
    double prev = dir * ta;
    for (std::size_t i = 0; i < chain.size(); ++i) {
        const auto& nh = chain[i];
        if (!seen.insert(nh).second)
            return {ChainStatus::DuplicateNode, i};
        const ParamPoint* pp = cs.paramOf(nh);
        if (!pp) return {ChainStatus::MissingParameter, i};
        if (pp->chart != seg.chart) return {ChainStatus::ChartMismatch, i};
        // exact collinearity + exact bounds in stored coordinates
        if (predicates::orient2d(pa->p, pb->p, pp->p) !=
            predicates::Sign::Zero)
            return {ChainStatus::OffSegment, i};
        const double lo = std::min(pa->p.u, pb->p.u),
                     hi = std::max(pa->p.u, pb->p.u);
        const double lov = std::min(pa->p.v, pb->p.v),
                     hiv = std::max(pa->p.v, pb->p.v);
        if (pp->p.u < lo || pp->p.u > hi || pp->p.v < lov || pp->p.v > hiv)
            return {ChainStatus::OffSegment, i};
        // strictly increasing dominant-axis parameter
        const double t = dir * detail::axisCoord(pp->p, axis);
        if (i == 0) {
            prev = t;
        } else {
            if (!(t > prev)) return {ChainStatus::NonIncreasingParameter, i};
            prev = t;
        }
    }
    for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
        const auto eh =
            ws.edges().find(ws.nodes(), chain[i], chain[i + 1]);
        if (!eh) return {ChainStatus::MissingEdge, i};
        const auto d = ws.edges().get(*eh);
        if (!d.requiredConstraint) return {ChainStatus::EdgeNotRequired, i};
        if (d.cls.dim != seg.cls.dim || d.cls.tag != seg.cls.tag)
            return {ChainStatus::ClassificationMismatch, i};
    }
    return {ChainStatus::Valid, 0};
}

// Chain DISCOVERY: walk from map(a) along required, exactly-on-segment
// edges with strictly increasing dominant-axis parameter, requiring a
// UNIQUE next node at every step. Returns the chain INCLUDING both
// endpoints, or nullopt. Refuses a stale or foreign index.
inline std::optional<std::vector<topo::NodeHandle>> discoverChain(
    const topo::TopologyWorkspace& ws, const ConstraintSet& cs,
    const TriAdjacencyIndex& idx, const ConstraintSegment& seg) {
    if (!idx.current(ws)) return std::nullopt;
    const auto ia = cs.vertices.find(seg.a);
    const auto ib = cs.vertices.find(seg.b);
    if (ia == cs.vertices.end() || ib == cs.vertices.end())
        return std::nullopt;
    const ParamPoint* pa = cs.paramOf(ia->second);
    const ParamPoint* pb = cs.paramOf(ib->second);
    if (!pa || !pb) return std::nullopt;
    const int axis = detail::dominantAxis(pa->p, pb->p);
    const double dir =
        detail::axisCoord(pb->p, axis) > detail::axisCoord(pa->p, axis)
            ? 1.0
            : -1.0;
    std::vector<topo::NodeHandle> chain{ia->second};
    std::unordered_set<topo::NodeHandle> visited{ia->second};
    const std::size_t bound = ws.nodes().activeBySequence().size() + 1;
    while (!(chain.back() == ib->second) && chain.size() <= bound) {
        const auto& cur = chain.back();
        const ParamPoint* pc = cs.paramOf(cur);
        if (!pc) return std::nullopt;
        const double tc = dir * detail::axisCoord(pc->p, axis);
        std::optional<topo::NodeHandle> next;
        const auto sit = idx.star.find(cur);
        if (sit == idx.star.end()) return std::nullopt;
        std::unordered_set<topo::NodeHandle> tried;
        for (const auto& ch : sit->second) {
            const auto cd = ws.cells().get(ch);
            for (int i = 0; i < 3; ++i) {
                const auto& cand = cd.nodes[std::size_t(i)];
                if (cand == cur) continue;
                if (!tried.insert(cand).second) continue;
                if (visited.count(cand)) continue;
                const ParamPoint* pn = cs.paramOf(cand);
                if (!pn || pn->chart != seg.chart) continue;
                if (predicates::orient2d(pa->p, pb->p, pn->p) !=
                    predicates::Sign::Zero)
                    continue;
                const double tn = dir * detail::axisCoord(pn->p, axis);
                if (!(tn > tc)) continue;
                if (tn > dir * detail::axisCoord(pb->p, axis)) continue;
                const auto eh = ws.edges().find(ws.nodes(), cur, cand);
                if (!eh) continue;
                if (!ws.edges().get(*eh).requiredConstraint) continue;
                if (next) return std::nullopt;   // ambiguous: not a chain
                next = cand;
            }
        }
        if (!next) return std::nullopt;
        visited.insert(*next);
        chain.push_back(*next);
    }
    if (!(chain.back() == ib->second)) return std::nullopt;
    return chain;
}

}  // namespace paacs::bnd
