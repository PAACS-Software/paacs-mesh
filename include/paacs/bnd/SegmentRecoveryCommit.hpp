// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// ============================================================================
// BND-1 candidate-generation slice, increment 3 - THE MUTATION PATH.
// Design authority: docs/PMBND001_DESIGN.md sec 13 rev 3 (ACCEPTED).
//
// Implements: 13.3-0 normative routing, 13.3/13.3a mark-only cavity with
// full-star disk closure and its own disk audit, 13.3b classification
// policy (aligned with the frozen ensure-only legality), 13.4 the
// RecoveryCandidate contract, 13.4a the six-gate candidate preview
// validator, 13.5 single-classification rule, 13.6 failure tiers,
// 13.7 prefix stopping rule, 13.8/13.8a the TRANSVERSE strip walk with
// its canonical initial cell, and the 13.9 GEOMETRIC ear rank over VALID
// ears only (an ear is eligible only when its apex is convex, no other
// polygon vertex lies inside the ear or ON the new diagonal, and the
// diagonal crosses no nonadjacent polygon edge).
//
// Outcome discipline for the transverse path: the walk's structured
// degeneracies (interior-node hit, collinear overlap, tangential touch,
// boundary reached, walk guard, crossed required edge) are Refusals with
// zero mutation. Once BND-1's own disk audit has ACCEPTED the strip, a
// failure to triangulate it is NOT a policy refusal but a candidate-
// generation failure: Failed(InvalidCandidateGeometry), tier 1 of 13.6
// (nothing prepared -> FULL observable digest unchanged). There is no
// TransverseNotImplemented outcome; that placeholder enum was removed.
//
// ALL mutation flows through the FROZEN M4-B.8 cavity transaction. This
// header performs no pool writes of its own and makes ZERO frozen edits.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <unordered_set>
#include <array>
#include <map>
#include <optional>
#include <vector>

#include "paacs/bnd/SegmentRecoverySession.hpp"
#include "paacs/topo/Cavity.hpp"

namespace paacs::bnd {

// ---------------------------------------------------------------------------
// 13.5 outcomes (sec 5). "Unattempted" is reporting metadata, not a sixth.
// ---------------------------------------------------------------------------
enum class Outcome { AlreadyPresent, Recovered, Refused, Failed, Rejected };

enum class RefusedReason {
    None,
    PartialChain,             // 13.3-0 case 3
    ClosureNotADisk,          // 13.3a
    MixedCellClassification,  // 13.5
    ConflictingEdgeMark,      // 13.3b
    InteriorNodeHit,          // 13.8
    CrossedRequiredEdge,      // 13.8: retiring a required edge
    UnknownChart,
    CollinearOverlapWithEdge, // 13.8
    TangentialTouch,          // 13.8
    BoundaryReached,          // 13.8
    NoUniqueInitialCell,      // 13.8a
    WalkGuard,                // 13.8
};

enum class FailedReason {
    None,
    StampChainBroken,           // 13.2 four-way identity            (tier 1)
    InvalidCandidateGeometry,   // 13.4a validator rejected          (tier 1)
    CavityRejected,             // frozen stage/validate logic_error (tier 1)
    CandidateAllocFailed,       // alloc failure BEFORE prepare():
                                //   routing, candidate build, or the
                                //   cavity stage/validate zone       (tier 1)
    CavityThrew,                // alloc failure INSIDE prepare()/commit()
                                //   - the ONLY prepare-zone reason   (tier 2)
    IndexRebuildFailed,         // 13.2a - owned by the NEXT entity   (tier 1)
};

struct EntityResult {
    ConstraintEntityId entity = 0;
    Outcome outcome{};
    RefusedReason refused = RefusedReason::None;
    FailedReason failed = FailedReason::None;
};

struct RecoveryReport {
    std::vector<EntityResult> processed;          // ascending entity id
    std::vector<ConstraintEntityId> unattempted;  // 13.7
    bool rejected = false;                        // whole-set preflight
    SessionStatus rejectReason{};                 // preserved (13.7)
};

// 13.4 - the candidate contract, built entirely off-workspace.
struct RecoveryCandidate {
    OwnershipStamp stamp;
    ConstraintEntityId entity = 0;
    ChartId chart = 0;
    std::vector<topo::CellHandle> removed;
    std::vector<topo::TriCavityCell> added;
    mesh::Classification cellClass{};
    std::vector<topo::CavityEdgeRequirement> requirements;
    std::vector<topo::NodeHandle> resultChain;    // EVIDENCE only
    std::vector<topo::NodeHandle> auditBoundary;  // EVIDENCE only
};

inline int detail_chain_axis(const Vec2& a, const Vec2& b) {
    const double du = std::fabs(b.u - a.u), dv = std::fabs(b.v - a.v);
    return dv > du ? 1 : 0;
}
inline double chain_axis_coord(const Vec2& p, int axis) {
    return axis == 0 ? p.u : p.v;
}

namespace detail {

inline const Vec2& pOf(const ConstraintSet& cs, const topo::NodeHandle& h) {
    static const Vec2 kZero{};
    const ParamPoint* pp = cs.paramOf(h);
    return pp ? pp->p : kZero;
}

// 13.3a - the canonical closure: full star of every INTERIOR chain node,
// plus the cells incident to the first and last chain edges.
inline std::vector<topo::CellHandle> canonicalClosure(
    const topo::TopologyWorkspace& ws, const TriAdjacencyIndex& idx,
    const std::vector<topo::NodeHandle>& chain) {
    std::vector<topo::CellHandle> out;
    std::unordered_set<topo::CellHandle> seen;
    auto addCell = [&](const topo::CellHandle& c) {
        if (seen.insert(c).second) out.push_back(c);
    };
    for (std::size_t i = 1; i + 1 < chain.size(); ++i) {
        const auto it = idx.star.find(chain[i]);
        if (it == idx.star.end()) continue;
        for (const auto& c : it->second) addCell(c);
    }
    auto addEdgeCells = [&](const topo::NodeHandle& x,
                            const topo::NodeHandle& y) {
        const auto it = idx.star.find(x);
        if (it == idx.star.end()) return;
        for (const auto& c : it->second) {
            const auto cd = ws.cells().get(c);
            for (int k = 0; k < 3; ++k)
                if (cd.nodes[std::size_t(k)] == y) {
                    addCell(c);
                    break;
                }
        }
    };
    if (chain.size() >= 2) {
        addEdgeCells(chain[0], chain[1]);
        addEdgeCells(chain[chain.size() - 2], chain.back());
    }
    // canonical order: by full CellHandle, so the candidate is
    // reproducible independent of star iteration order
    std::sort(out.begin(), out.end());
    return out;
}

// 13.3a - BND-1's OWN disk audit of the closure, before any handoff:
// edge-connected, one simple boundary cycle, single node link, chi == 1.
inline bool closureIsDisk(const topo::TopologyWorkspace& ws,
                          const std::vector<topo::CellHandle>& cells,
                          std::vector<topo::NodeHandle>* boundaryOut) {
    if (boundaryOut) boundaryOut->clear();
    if (cells.empty()) return false;
    std::unordered_set<topo::CellHandle> inSet(cells.begin(), cells.end());
    std::unordered_map<std::int32_t, std::vector<topo::CellHandle>> byEdge;
    for (const auto& c : cells)
        for (const auto& eh : ws.edges().cellEdges(ws.cells(), c))
            byEdge[eh.slot].push_back(c);
    // 1. edge-connectivity over shared edges
    std::unordered_set<topo::CellHandle> seen;
    std::vector<topo::CellHandle> stack{cells.front()};
    while (!stack.empty()) {
        const auto c = stack.back();
        stack.pop_back();
        if (!seen.insert(c).second) continue;
        for (const auto& eh : ws.edges().cellEdges(ws.cells(), c))
            for (const auto& n : byEdge[eh.slot])
                if (inSet.count(n) && !seen.count(n)) stack.push_back(n);
    }
    if (seen.size() != cells.size()) return false;
    // 2. Euler: V - E + F == 1
    std::unordered_set<topo::NodeHandle> V;
    std::unordered_set<std::int32_t> E;
    for (const auto& c : cells) {
        const auto cd = ws.cells().get(c);
        for (int k = 0; k < 3; ++k) V.insert(cd.nodes[std::size_t(k)]);
        for (const auto& eh : ws.edges().cellEdges(ws.cells(), c))
            E.insert(eh.slot);
    }
    if (static_cast<long long>(V.size()) - static_cast<long long>(E.size()) +
            static_cast<long long>(cells.size()) !=
        1)
        return false;
    // 3. NODE-LINK audit: for every node, the opposite edges of its
    // incident triangles IN THE CLOSURE must form ONE path or cycle
    // (degree <= 2 everywhere AND a single connected component).
    {
        std::unordered_map<topo::NodeHandle,
                           std::vector<std::array<topo::NodeHandle, 2>>>
            link;
        for (const auto& c : cells) {
            const auto cd = ws.cells().get(c);
            for (int i = 0; i < 3; ++i)
                link[cd.nodes[std::size_t(i)]].push_back(
                    {cd.nodes[std::size_t((i + 1) % 3)],
                     cd.nodes[std::size_t((i + 2) % 3)]});
        }
        for (const auto& [n, segs] : link) {
            std::unordered_map<topo::NodeHandle,
                               std::vector<topo::NodeHandle>>
                adj;
            std::unordered_map<topo::NodeHandle, int> deg;
            for (const auto& s : segs) {
                if (++deg[s[0]] > 2 || ++deg[s[1]] > 2) return false;
                adj[s[0]].push_back(s[1]);
                adj[s[1]].push_back(s[0]);
            }
            if (adj.empty()) return false;
            std::unordered_set<topo::NodeHandle> ls;
            std::vector<topo::NodeHandle> st{segs[0][0]};
            while (!st.empty()) {
                const auto v = st.back();
                st.pop_back();
                if (!ls.insert(v).second) continue;
                for (const auto& w : adj[v]) st.push_back(w);
            }
            if (ls.size() != adj.size()) return false;   // pinched
        }
    }
    // 4. BOUNDARY: collect boundary edges (used once inside the set) and
    // TRAVERSE them, proving they form exactly ONE cycle.
    std::vector<std::array<topo::NodeHandle, 2>> bedges;
    for (const auto& [slot, owners] : byEdge) {
        int inCount = 0;
        for (const auto& c : owners)
            if (inSet.count(c)) ++inCount;
        if (inCount != 1) continue;
        for (const auto& c : owners) {
            if (!inSet.count(c)) continue;
            const auto ce = ws.edges().cellEdges(ws.cells(), c);
            const auto cd = ws.cells().get(c);
            for (int k = 0; k < 3; ++k)
                if (ce[std::size_t(k)].slot == slot)
                    bedges.push_back({cd.nodes[std::size_t(k)],
                                      cd.nodes[std::size_t((k + 1) % 3)]});
            break;
        }
    }
    if (bedges.empty()) return false;
    std::unordered_map<topo::NodeHandle, std::vector<topo::NodeHandle>> badj;
    for (const auto& e : bedges) {
        badj[e[0]].push_back(e[1]);
        badj[e[1]].push_back(e[0]);
    }
    for (const auto& [n, nb] : badj)
        if (nb.size() != 2) return false;         // every boundary node deg 2
    // single cycle: walk it and require ALL boundary nodes visited
    std::vector<topo::NodeHandle> ring;
    {
        topo::NodeHandle start = badj.begin()->first;
        for (const auto& [n, nb] : badj)
            if (n < start) start = n;             // canonical entry point
        topo::NodeHandle prev{}, cur = start;
        std::unordered_set<topo::NodeHandle> visited;
        while (visited.insert(cur).second) {
            ring.push_back(cur);
            const auto& nb = badj[cur];
            const topo::NodeHandle nxt =
                (nb[0] == prev) ? nb[1] : nb[0];
            prev = cur;
            cur = nxt;
        }
        if (!(cur == start)) return false;        // not a closed cycle
        if (ring.size() != badj.size()) return false;   // more than one
    }
    if (boundaryOut) *boundaryOut = ring;
    return true;
}

// Exact: is q collinear with segment a..b and inside its closed extent?
inline bool pointOnSegment(const Vec2& a, const Vec2& b, const Vec2& q) {
    if (predicates::orient2d(a, b, q) != predicates::Sign::Zero) return false;
    const double lo = std::min(a.u, b.u), hi = std::max(a.u, b.u);
    const double lov = std::min(a.v, b.v), hiv = std::max(a.v, b.v);
    return q.u >= lo && q.u <= hi && q.v >= lov && q.v <= hiv;
}

// Exact: do open segments a..b and c..d cross PROPERLY (interiors meet,
// no endpoint incidence, no collinear overlap)?
inline bool segmentsProperlyCross(const Vec2& a, const Vec2& b, const Vec2& c,
                                  const Vec2& d) {
    const auto d1 = predicates::orient2d(a, b, c);
    const auto d2 = predicates::orient2d(a, b, d);
    const auto d3 = predicates::orient2d(c, d, a);
    const auto d4 = predicates::orient2d(c, d, b);
    const bool s1 = d1 != d2 && d1 != predicates::Sign::Zero &&
                    d2 != predicates::Sign::Zero;
    const bool s2 = d3 != d4 && d3 != predicates::Sign::Zero &&
                    d4 != predicates::Sign::Zero;
    return s1 && s2;
}

// 13.9 - the FULL validity condition an ear at ring index `i` must satisfy
// BEFORE it competes for the geometric rank. Ranking a convex apex without
// these checks let a concave polygon emit a triangulation the six-gate
// candidate validator then rejected, while a different (equally reachable)
// triangulation of the same polygon passed - the review reproduction. A
// valid ear requires ALL of:
//   (a) the apex is convex in the accepted winding;
//   (b) no OTHER polygon vertex lies strictly inside the ear triangle;
//   (c) no OTHER polygon vertex lies ON the new diagonal prev..next;
//   (d) the diagonal crosses NO nonadjacent polygon edge.
inline bool earValid(const ConstraintSet& cs,
                     const std::vector<topo::NodeHandle>& ring, std::size_t i,
                     predicates::Sign want) {
    const std::size_t n = ring.size();
    const auto& prev = ring[(i + n - 1) % n];
    const auto& cur = ring[i];
    const auto& next = ring[(i + 1) % n];
    const Vec2 &pp = pOf(cs, prev), &pc = pOf(cs, cur), &pn = pOf(cs, next);
    if (predicates::orient2d(pp, pc, pn) != want) return false;   // (a)
    for (std::size_t j = 0; j < n; ++j) {
        if (ring[j] == prev || ring[j] == cur || ring[j] == next) continue;
        const Vec2& pq = pOf(cs, ring[j]);
        const auto s1 = predicates::orient2d(pp, pc, pq);
        const auto s2 = predicates::orient2d(pc, pn, pq);
        const auto s3 = predicates::orient2d(pn, pp, pq);
        if (s1 != predicates::Sign::Zero && s1 == s2 && s2 == s3)
            return false;                                          // (b)
        if (pointOnSegment(pp, pn, pq)) return false;              // (c)
    }
    for (std::size_t e = 0; e < n; ++e) {                          // (d)
        const auto& r = ring[e];
        const auto& s = ring[(e + 1) % n];
        // edges incident to the diagonal's endpoints are adjacent by
        // construction and share a point, never a proper crossing.
        if (r == prev || r == next || s == prev || s == next) continue;
        if (segmentsProperlyCross(pp, pn, pOf(cs, r), pOf(cs, s)))
            return false;
    }
    return true;
}

// 13.9 - GEOMETRIC ear clipping. Polygon given as an ordered ring of
// handles; returns triangles as index triples into the ring, or nullopt
// if no valid ear exists (degenerate / non-simple input). Only ears that
// pass earValid() compete, and among them the APEX with the
// lexicographically smallest (u, v) is clipped first (never handle order).
inline std::optional<std::vector<std::array<topo::NodeHandle, 3>>> earClip(
    const ConstraintSet& cs, std::vector<topo::NodeHandle> ring,
    predicates::Sign want) {
    std::vector<std::array<topo::NodeHandle, 3>> out;
    while (ring.size() > 3) {
        std::optional<std::size_t> best;
        Vec2 bestApex{};
        for (std::size_t i = 0; i < ring.size(); ++i) {
            if (!earValid(cs, ring, i, want)) continue;
            const Vec2& pc = pOf(cs, ring[i]);
            // 13.9: rank by APEX coordinate, never by handle/slot order
            if (!best || pc.u < bestApex.u ||
                (pc.u == bestApex.u && pc.v < bestApex.v)) {
                best = i;
                bestApex = pc;
            }
        }
        if (!best) return std::nullopt;
        const std::size_t i = *best;
        out.push_back({ring[(i + ring.size() - 1) % ring.size()], ring[i],
                       ring[(i + 1) % ring.size()]});
        ring.erase(ring.begin() + static_cast<std::ptrdiff_t>(i));
    }
    if (ring.size() == 3) out.push_back({ring[0], ring[1], ring[2]});
    return out;
}

// Exact winding of a simple polygon: the lexicographically smallest
// (u, v) vertex is always STRICTLY CONVEX on the hull, so the turn taken
// there IS the polygon's winding sign - computed exactly by orient2d, with
// no signed-area epsilon. Zero means the extreme vertex is collinear with
// both neighbours (a degenerate/non-simple ring): reported as no winding.
inline predicates::Sign polygonWinding(
    const ConstraintSet& cs, const std::vector<topo::NodeHandle>& ring) {
    const std::size_t n = ring.size();
    std::size_t lo = 0;
    for (std::size_t i = 1; i < n; ++i) {
        const Vec2& p = pOf(cs, ring[i]);
        const Vec2& b = pOf(cs, ring[lo]);
        if (p.u < b.u || (p.u == b.u && p.v < b.v)) lo = i;
    }
    return predicates::orient2d(pOf(cs, ring[(lo + n - 1) % n]),
                                pOf(cs, ring[lo]),
                                pOf(cs, ring[(lo + 1) % n]));
}

// Clip a polygon and normalize every emitted triangle to the accepted
// chart orientation. CRITICAL: convexity ("is this an ear?") is defined
// relative to the polygon's OWN winding, never the chart's `want`. The two
// sides of a strip wind oppositely by construction, so one side winds
// against `want`; testing ears with `want` there would mark REFLEX
// vertices as convex and clip them, producing an overlapping triangulation
// the six-gate validator then rejects (the review reproduction). We take
// the polygon's true winding, clip valid ears in that winding, then swap
// each triangle to `want` for the chart.
inline std::optional<std::vector<std::array<topo::NodeHandle, 3>>>
earClipOriented(const ConstraintSet& cs,
                const std::vector<topo::NodeHandle>& ring,
                predicates::Sign want) {
    if (ring.size() < 3) return std::nullopt;
    const predicates::Sign winding = polygonWinding(cs, ring);
    if (winding == predicates::Sign::Zero) return std::nullopt;
    auto tris = earClip(cs, ring, winding);
    if (!tris) return std::nullopt;
    for (auto& t : *tris) {
        const auto s = predicates::orient2d(pOf(cs, t[0]), pOf(cs, t[1]),
                                            pOf(cs, t[2]));
        if (s == predicates::Sign::Zero) return std::nullopt;
        if (s != want) std::swap(t[1], t[2]);
    }
    return tris;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// 13.8/13.8a - the TRANSVERSE STRIP WALK and its canonical initial cell.
// All tests are exact orient2d in the authoritative stored coordinates.
// ---------------------------------------------------------------------------
namespace detail {

struct WalkResult {
    bool ok = false;
    RefusedReason why = RefusedReason::None;
    std::vector<topo::CellHandle> strip;   // crossed triangles, in order
    std::vector<topo::EdgeHandle> crossed; // the edges crossed, in order
    std::vector<topo::NodeHandle> left;    // boundary polygon, one side
    std::vector<topo::NodeHandle> right;   // boundary polygon, other side
};

// 13.8a: the unique incident triangle whose OPEN WEDGE at `a` strictly
// contains direction a->b. Zero, several, or any zero orientation is a
// structured refusal.
inline bool initialCell(const topo::TopologyWorkspace& ws,
                        const ConstraintSet& cs, const TriAdjacencyIndex& idx,
                        const topo::NodeHandle& a, const topo::NodeHandle& b,
                        topo::CellHandle& out, RefusedReason& why) {
    const auto sit = idx.star.find(a);
    if (sit == idx.star.end()) {
        why = RefusedReason::NoUniqueInitialCell;
        return false;
    }
    const Vec2& pa = pOf(cs, a);
    const Vec2& pb = pOf(cs, b);
    int found = 0;
    for (const auto& ch : sit->second) {
        const auto cd = ws.cells().get(ch);
        topo::NodeHandle p{}, q{};
        int k = 0;
        for (int i = 0; i < 3; ++i)
            if (!(cd.nodes[std::size_t(i)] == a)) {
                (k++ == 0 ? p : q) = cd.nodes[std::size_t(i)];
            }
        if (k != 2) continue;
        const Vec2& pp = pOf(cs, p);
        const Vec2& pq = pOf(cs, q);
        // direction a->b strictly inside the wedge (a,p),(a,q)
        const auto s1 = predicates::orient2d(pa, pp, pb);
        const auto s2 = predicates::orient2d(pa, pb, pq);
        const auto s3 = predicates::orient2d(pa, pp, pq);
        if (s1 == predicates::Sign::Zero || s2 == predicates::Sign::Zero) {
            // b's direction is collinear with an incident edge: this is
            // the vertex-hit / collinear case, detected at step one.
            why = RefusedReason::CollinearOverlapWithEdge;
            return false;
        }
        if (s3 == predicates::Sign::Zero) continue;   // degenerate fan cell
        if (s1 == s3 && s2 == s3) {
            out = ch;
            ++found;
        }
    }
    if (found == 1) return true;
    why = RefusedReason::NoUniqueInitialCell;
    return false;
}

// The walk itself (13.8). Collects the crossed strip and the two boundary
// chains; refuses on every degeneracy the design enumerates.
inline WalkResult transverseWalk(const topo::TopologyWorkspace& ws,
                                 const ConstraintSet& cs,
                                 const TriAdjacencyIndex& idx,
                                 const topo::NodeHandle& a,
                                 const topo::NodeHandle& b) {
    WalkResult r;
    topo::CellHandle cur{};
    if (!initialCell(ws, cs, idx, a, b, cur, r.why)) return r;
    const Vec2& pa = pOf(cs, a);
    const Vec2& pb = pOf(cs, b);
    const std::size_t guard = ws.cells().activeBySequence().size() + 1;
    std::vector<topo::NodeHandle> leftChain{a}, rightChain{a};
    topo::NodeHandle fromU{}, fromV{};   // the edge just crossed
    bool haveEdge = false;
    while (r.strip.size() <= guard) {
        r.strip.push_back(cur);
        const auto cd = ws.cells().get(cur);
        // reached b?
        bool hasB = false;
        for (int i = 0; i < 3; ++i)
            if (cd.nodes[std::size_t(i)] == b) hasB = true;
        if (hasB) {
            // close the two chains at b
            leftChain.push_back(b);
            rightChain.push_back(b);
            r.left = leftChain;
            r.right = rightChain;
            r.ok = true;
            return r;
        }
        // find the unique OTHER edge properly crossed by the open segment
        const auto ce = ws.edges().cellEdges(ws.cells(), cur);
        int crossedLocal = -1;
        topo::NodeHandle cu{}, cv{};
        int crossings = 0;
        for (int i = 0; i < 3; ++i) {
            const auto& x = cd.nodes[std::size_t(i)];
            const auto& y = cd.nodes[std::size_t((i + 1) % 3)];
            if (haveEdge && ((x == fromU && y == fromV) ||
                             (x == fromV && y == fromU)))
                continue;                       // the edge we entered by
            const Vec2& px = pOf(cs, x);
            const Vec2& py = pOf(cs, y);
            const auto d1 = predicates::orient2d(pa, pb, px);
            const auto d2 = predicates::orient2d(pa, pb, py);
            const auto d3 = predicates::orient2d(px, py, pa);
            const auto d4 = predicates::orient2d(px, py, pb);
            if (d1 == predicates::Sign::Zero ||
                d2 == predicates::Sign::Zero) {
                // the segment passes exactly through a mesh vertex
                if (!(x == a || y == a || x == b || y == b)) {
                    r.why = RefusedReason::InteriorNodeHit;
                    return r;
                }
                continue;
            }
            const bool s1 = d1 != d2;
            const bool s2 = d3 != d4 && d3 != predicates::Sign::Zero &&
                            d4 != predicates::Sign::Zero;
            if (s1 && s2) {
                crossedLocal = i;
                cu = x;
                cv = y;
                ++crossings;
            }
        }
        if (crossings == 0) {
            r.why = RefusedReason::TangentialTouch;
            return r;
        }
        if (crossings > 1) {
            r.why = RefusedReason::WalkGuard;
            return r;
        }
        // record the two crossed-edge endpoints on their respective sides
        {
            const auto sideU = predicates::orient2d(pa, pb, pOf(cs, cu));
            if (sideU == predicates::Sign::Positive) {
                if (leftChain.empty() || !(leftChain.back() == cu))
                    leftChain.push_back(cu);
                if (rightChain.empty() || !(rightChain.back() == cv))
                    rightChain.push_back(cv);
            } else {
                if (rightChain.empty() || !(rightChain.back() == cu))
                    rightChain.push_back(cu);
                if (leftChain.empty() || !(leftChain.back() == cv))
                    leftChain.push_back(cv);
            }
        }
        // 13.8: recovery RETIRES every crossed interior edge, so a
        // crossed edge that is already REQUIRED must be refused here as a
        // BND policy decision - never handed to the frozen cavity.
        {
            const auto eh = ws.edges().find(ws.nodes(), cu, cv);
            if (eh) {
                if (ws.edges().get(*eh).requiredConstraint) {
                    r.why = RefusedReason::CrossedRequiredEdge;
                    return r;
                }
                r.crossed.push_back(*eh);
            }
        }
        // step across
        const auto nb = idx.neighbors.find(cur);
        if (nb == idx.neighbors.end()) {
            r.why = RefusedReason::BoundaryReached;
            return r;
        }
        const auto next = nb->second[std::size_t(crossedLocal)];
        if (!ws.cells().isCurrent(next)) {
            r.why = RefusedReason::BoundaryReached;   // crossed a boundary
            return r;
        }
        fromU = cu;
        fromV = cv;
        haveEdge = true;
        cur = next;
        (void)ce;
    }
    r.why = RefusedReason::WalkGuard;
    return r;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// 13.4a - the six-gate candidate preview validator, supplied to the FROZEN
// Phase-1 validate() hook. Plain CavityPreview values only; no workspace,
// no pools. Captures the session snapshot + accepted chart orientation.
// ---------------------------------------------------------------------------
inline std::function<bool(const topo::CavityPreview&)> makeCandidateValidator(
    const ConstraintSet& cs, ChartId chart, predicates::Sign want) {
    return [&cs, chart, want](const topo::CavityPreview& pv) -> bool {
        // 1. complete SAME-CHART parameter coverage
        for (const auto& tri : pv.addedConnectivity)
            for (const auto& n : tri) {
                const ParamPoint* pp = cs.paramOf(n);
                if (!pp || pp->chart != chart) return false;
            }
        // 2/3. nondegeneracy + accepted chart orientation
        for (const auto& tri : pv.addedConnectivity) {
            const auto s = predicates::orient2d(detail::pOf(cs, tri[0]),
                                                detail::pOf(cs, tri[1]),
                                                detail::pOf(cs, tri[2]));
            if (s == predicates::Sign::Zero) return false;
            if (s != want) return false;
        }
        // 4/5/6. crossings, T-junctions, overlap among the ADDED cells
        struct Seg { topo::NodeHandle a, b; };
        std::vector<Seg> E;
        for (const auto& tri : pv.addedConnectivity)
            for (int i = 0; i < 3; ++i) {
                topo::NodeHandle x = tri[std::size_t(i)],
                                 y = tri[std::size_t((i + 1) % 3)];
                if (y < x) std::swap(x, y);
                bool dup = false;
                for (const auto& e : E)
                    if (e.a == x && e.b == y) { dup = true; break; }
                if (!dup) E.push_back({x, y});
            }
        for (std::size_t i = 0; i < E.size(); ++i)
            for (std::size_t j = i + 1; j < E.size(); ++j) {
                if (E[i].a == E[j].a || E[i].a == E[j].b ||
                    E[i].b == E[j].a || E[i].b == E[j].b)
                    continue;
                const Vec2 &a1 = detail::pOf(cs, E[i].a),
                           &b1 = detail::pOf(cs, E[i].b),
                           &a2 = detail::pOf(cs, E[j].a),
                           &b2 = detail::pOf(cs, E[j].b);
                const auto d1 = predicates::orient2d(a1, b1, a2);
                const auto d2 = predicates::orient2d(a1, b1, b2);
                const auto d3 = predicates::orient2d(a2, b2, a1);
                const auto d4 = predicates::orient2d(a2, b2, b1);
                const bool s1 = d1 != d2 && d1 != predicates::Sign::Zero &&
                                d2 != predicates::Sign::Zero;
                const bool s2 = d3 != d4 && d3 != predicates::Sign::Zero &&
                                d4 != predicates::Sign::Zero;
                if (s1 && s2) return false;               // 4. crossing
            }
        std::vector<topo::NodeHandle> V;
        for (const auto& tri : pv.addedConnectivity)
            for (const auto& n : tri) {
                bool dup = false;
                for (const auto& v : V)
                    if (v == n) { dup = true; break; }
                if (!dup) V.push_back(n);
            }
        for (const auto& v : V)
            for (const auto& e : E) {
                if (v == e.a || v == e.b) continue;
                const Vec2 &pa = detail::pOf(cs, e.a),
                           &pb = detail::pOf(cs, e.b), &pv2 = detail::pOf(cs, v);
                if (predicates::orient2d(pa, pb, pv2) != predicates::Sign::Zero)
                    continue;
                const double lo = std::min(pa.u, pb.u), hi = std::max(pa.u, pb.u);
                const double lov = std::min(pa.v, pb.v),
                             hiv = std::max(pa.v, pb.v);
                if (pv2.u >= lo && pv2.u <= hi && pv2.v >= lov && pv2.v <= hiv)
                    return false;                          // 5. T-junction
            }
        for (std::size_t i = 0; i < pv.addedConnectivity.size(); ++i)
            for (std::size_t j = 0; j < pv.addedConnectivity.size(); ++j) {
                if (i == j) continue;
                const auto& gi = pv.addedConnectivity[i];
                const auto& gj = pv.addedConnectivity[j];
                for (int a = 0; a < 3; ++a) {
                    bool shared = false;
                    for (int b = 0; b < 3; ++b)
                        if (gi[std::size_t(a)] == gj[std::size_t(b)]) {
                            shared = true;
                            break;
                        }
                    if (shared) continue;
                    const auto s1 = predicates::orient2d(
                        detail::pOf(cs, gj[0]), detail::pOf(cs, gj[1]),
                        detail::pOf(cs, gi[std::size_t(a)]));
                    const auto s2 = predicates::orient2d(
                        detail::pOf(cs, gj[1]), detail::pOf(cs, gj[2]),
                        detail::pOf(cs, gi[std::size_t(a)]));
                    const auto s3 = predicates::orient2d(
                        detail::pOf(cs, gj[2]), detail::pOf(cs, gj[0]),
                        detail::pOf(cs, gi[std::size_t(a)]));
                    if (s1 != predicates::Sign::Zero && s1 == s2 && s2 == s3)
                        return false;                      // 6. overlap
                }
            }
        return true;
    };
}


// ---------------------------------------------------------------------------
// 13.3-0 routing + 13.3/13.3a/13.3b the mark-only path.
// ---------------------------------------------------------------------------
enum class Routing { AlreadyPresent, MarkOnly, PartialChain, Transverse };

// Classify an entity per the normative routing (13.3-0).
inline Routing routeEntity(const topo::TopologyWorkspace& ws,
                           const ConstraintSet& cs,
                           const TriAdjacencyIndex& idx,
                           const ConstraintSegment& seg,
                           std::vector<topo::NodeHandle>& chainOut) {
    // A complete GEOMETRIC chain: nodes exist along the segment and each
    // consecutive pair is a current edge (marks not required here).
    const auto ia = cs.vertices.find(seg.a);
    const auto ib = cs.vertices.find(seg.b);
    if (ia == cs.vertices.end() || ib == cs.vertices.end())
        return Routing::Transverse;
    const ParamPoint* pa = cs.paramOf(ia->second);
    const ParamPoint* pb = cs.paramOf(ib->second);
    if (!pa || !pb) return Routing::Transverse;
    const int axis = detail_chain_axis(pa->p, pb->p);
    const double dir =
        chain_axis_coord(pb->p, axis) > chain_axis_coord(pa->p, axis) ? 1.0
                                                                      : -1.0;
    // Walk geometric neighbours strictly increasing along the segment.
    std::vector<topo::NodeHandle> chain{ia->second};
    std::unordered_set<topo::NodeHandle> visited{ia->second};
    const std::size_t bound = ws.nodes().activeBySequence().size() + 1;
    bool progressed = true;
    while (!(chain.back() == ib->second) && chain.size() <= bound &&
           progressed) {
        progressed = false;
        const auto& cur = chain.back();
        const ParamPoint* pc = cs.paramOf(cur);
        if (!pc) break;
        const double tc = dir * chain_axis_coord(pc->p, axis);
        std::optional<topo::NodeHandle> next;
        const auto sit = idx.star.find(cur);
        if (sit == idx.star.end()) break;
        for (const auto& ch : sit->second) {
            const auto cd = ws.cells().get(ch);
            for (int i = 0; i < 3; ++i) {
                const auto& cand = cd.nodes[std::size_t(i)];
                if (cand == cur || visited.count(cand)) continue;
                const ParamPoint* pn = cs.paramOf(cand);
                if (!pn || pn->chart != seg.chart) continue;
                if (predicates::orient2d(pa->p, pb->p, pn->p) !=
                    predicates::Sign::Zero)
                    continue;
                const double tn = dir * chain_axis_coord(pn->p, axis);
                if (!(tn > tc)) continue;
                if (tn > dir * chain_axis_coord(pb->p, axis)) continue;
                if (!ws.edges().find(ws.nodes(), cur, cand)) continue;
                if (next && !(*next == cand)) return Routing::PartialChain;
                next = cand;
            }
        }
        if (next) {
            visited.insert(*next);
            chain.push_back(*next);
            progressed = true;
        }
    }
    if (!(chain.back() == ib->second)) {
        // no complete chain: if ANY on-segment edge was found the chain is
        // PARTIAL (13.3-0 case 3); otherwise transverse recovery.
        return chain.size() > 1 ? Routing::PartialChain : Routing::Transverse;
    }
    chainOut = chain;
    // marks: AlreadyPresent only when EVERY chain edge is required AND its
    // classification already matches (sec 5).
    bool allMarked = true;
    for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
        const auto eh = ws.edges().find(ws.nodes(), chain[i], chain[i + 1]);
        if (!eh) return Routing::PartialChain;
        const auto d = ws.edges().get(*eh);
        if (!d.requiredConstraint || d.cls.dim != seg.cls.dim ||
            d.cls.tag != seg.cls.tag) {
            allMarked = false;
            break;
        }
    }
    return allMarked ? Routing::AlreadyPresent : Routing::MarkOnly;
}

// 13.3b - frozen-aligned per-edge legality, resolved as a BND policy.
inline bool markLegal(const topo::TopologyWorkspace& ws,
                      const topo::EdgeHandle& eh,
                      const mesh::Classification& want) {
    const auto d = ws.edges().get(eh);
    // EXACT default only ({-1,-1}), matching the frozen primitive: a
    // partially-negative classification like {-1,-2} is NOT default.
    const mesh::Classification kDefault{};
    const bool isDefault =
        d.cls.dim == kDefault.dim && d.cls.tag == kDefault.tag;
    if (d.requiredConstraint)
        return d.cls.dim == want.dim && d.cls.tag == want.tag;
    return isDefault || (d.cls.dim == want.dim && d.cls.tag == want.tag);
}

// Build the mark-only candidate (13.3/13.3a/13.3b/13.5).
inline std::optional<RecoveryCandidate> buildMarkOnlyCandidate(
    const topo::TopologyWorkspace& ws, const RecoverySession& sess,
    const ConstraintSegment& seg, const std::vector<topo::NodeHandle>& chain,
    RefusedReason& why) {
    const auto& cs = sess.input();
    const auto closure = detail::canonicalClosure(ws, sess.index(), chain);
    std::vector<topo::NodeHandle> ring;
    if (closure.empty() || !detail::closureIsDisk(ws, closure, &ring)) {
        why = RefusedReason::ClosureNotADisk;
        return std::nullopt;
    }
    // 13.5: one identical classification across the removed strip
    mesh::Classification cls = ws.cells().get(closure.front()).cls;
    for (const auto& c : closure) {
        const auto cd = ws.cells().get(c);
        if (cd.cls.dim != cls.dim || cd.cls.tag != cls.tag) {
            why = RefusedReason::MixedCellClassification;
            return std::nullopt;
        }
    }
    // 13.3b: every chain edge must be legally markable
    for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
        const auto eh = ws.edges().find(ws.nodes(), chain[i], chain[i + 1]);
        if (!eh || !markLegal(ws, *eh, seg.cls)) {
            why = RefusedReason::ConflictingEdgeMark;
            return std::nullopt;
        }
    }
    RecoveryCandidate cand;
    cand.stamp = sess.stamp();
    cand.entity = seg.id;
    cand.chart = seg.chart;
    cand.removed = closure;
    cand.cellClass = cls;
    for (const auto& c : closure) {             // EXACT reproduction
        const auto cd = ws.cells().get(c);
        cand.added.push_back(topo::TriCavityCell{
            {cd.nodes[0], cd.nodes[1], cd.nodes[2]}, cd.cls});
    }
    for (std::size_t i = 0; i + 1 < chain.size(); ++i)
        cand.requirements.push_back(
            topo::CavityEdgeRequirement{chain[i], chain[i + 1], seg.cls});
    cand.resultChain = chain;
    cand.auditBoundary = ring;        // EVIDENCE (13.4), now populated
    (void)cs;
    return cand;
}

// 13.8/13.9 - build the TRANSVERSE candidate: walk the strip, ear-clip
// both boundary polygons in the parameter domain, and assemble the
// RecoveryCandidate. Off-workspace throughout; no mutation here.
// `fail` carries the ONE candidate-generation FAILURE this builder can
// raise: a strip that passed BND-1's own disk audit but could not be
// triangulated is Failed(InvalidCandidateGeometry), never a policy
// Refusal. Everything upstream of the disk audit is a structured Refusal
// reported through `why`. `fail` stays None on every Refusal.
inline std::optional<RecoveryCandidate> buildTransverseCandidate(
    const topo::TopologyWorkspace& ws, const RecoverySession& sess,
    const ConstraintSegment& seg, const topo::NodeHandle& a,
    const topo::NodeHandle& b, predicates::Sign want, RefusedReason& why,
    FailedReason& fail) {
    const auto& cs = sess.input();
    const auto w = detail::transverseWalk(ws, cs, sess.index(), a, b);
    if (!w.ok) {
        why = w.why;
        return std::nullopt;
    }
    // 13.5: one identical classification across the removed strip
    mesh::Classification cls = ws.cells().get(w.strip.front()).cls;
    for (const auto& c : w.strip) {
        const auto cd = ws.cells().get(c);
        if (cd.cls.dim != cls.dim || cd.cls.tag != cls.tag) {
            why = RefusedReason::MixedCellClassification;
            return std::nullopt;
        }
    }
    // 13.3a: the strip must be a disk before handoff (BND-1 audits it)
    std::vector<topo::NodeHandle> ring;
    if (!detail::closureIsDisk(ws, w.strip, &ring)) {
        why = RefusedReason::ClosureNotADisk;
        return std::nullopt;
    }
    // 13.9: ear-clip each side. Each polygon is the boundary chain plus
    // the new segment edge (a..b), so the two fans share exactly that edge.
    auto clipSide = [&](std::vector<topo::NodeHandle> chain)
        -> std::optional<std::vector<std::array<topo::NodeHandle, 3>>> {
        if (chain.size() < 3) return std::nullopt;
        return detail::earClipOriented(cs, chain, want);
    };
    auto leftTris = clipSide(w.left);
    auto rightTris = clipSide(w.right);
    if (!leftTris || !rightTris) {
        // The strip already passed BND-1's disk audit above, so this is
        // NOT a policy refusal: it is a candidate-generation failure of a
        // strip we accepted. Report it as Failed(InvalidCandidateGeometry)
        // (13.6 tier 1 - nothing is prepared, so the caller can assert
        // FULL observable-digest equality), never Refused(ClosureNotADisk).
        fail = FailedReason::InvalidCandidateGeometry;
        return std::nullopt;
    }
    RecoveryCandidate cand;
    cand.stamp = sess.stamp();
    cand.entity = seg.id;
    cand.chart = seg.chart;
    cand.removed = w.strip;
    std::sort(cand.removed.begin(), cand.removed.end());
    cand.cellClass = cls;
    for (const auto& tri : *leftTris)
        cand.added.push_back(topo::TriCavityCell{{tri[0], tri[1], tri[2]},
                                                 cls});
    for (const auto& tri : *rightTris)
        cand.added.push_back(topo::TriCavityCell{{tri[0], tri[1], tri[2]},
                                                 cls});
    cand.requirements.push_back(
        topo::CavityEdgeRequirement{a, b, seg.cls});
    cand.resultChain = {a, b};
    cand.auditBoundary = ring;
    return cand;
}

// ---------------------------------------------------------------------------
// The per-entity driver: candidate -> frozen cavity transaction.
// 13.6 failure tiers are inherited from the frozen primitive.
// ---------------------------------------------------------------------------
inline EntityResult commitCandidate(topo::TopologyWorkspace& ws,
                                    RecoverySession& sess,
                                    const RecoveryCandidate& cand,
                                    predicates::Sign want) {
    EntityResult r;
    r.entity = cand.entity;
    if (!sess.handoffAllowed(ws, cand.stamp)) {
        r.outcome = Outcome::Failed;
        r.failed = FailedReason::StampChainBroken;
        return r;
    }
    try {
        topo::TriCavitySpec spec;
        spec.removed = cand.removed;
        spec.added = cand.added;
        spec.requirements = cand.requirements;
        auto tx = ws.beginCavityTransaction();
        // PHASE SPLIT (13.6): stage() + validate() run BEFORE prepare(), so
        // an allocation failure here prepared NOTHING (tier 1). We catch it
        // separately as CandidateAllocFailed so the ONLY reason that ever
        // escapes the prepare/commit zone is CavityThrew (tier 2). This is
        // what lets the probe attribute tiers by reason without ambiguity.
        try {
            (void)tx.stage(spec);
            // 13.4a: the six-gate validator runs BEFORE prepare()
            if (!tx.validate(makeCandidateValidator(sess.input(), cand.chart,
                                                    want))) {
                tx.abandon();
                r.outcome = Outcome::Failed;
                r.failed = FailedReason::InvalidCandidateGeometry;
                return r;
            }
        } catch (const std::bad_alloc&) {
            // stage/validate zone: nothing prepared (tier 1)
            r.outcome = Outcome::Failed;
            r.failed = FailedReason::CandidateAllocFailed;
            return r;
        }
        (void)tx.commit();   // prepare() + commit(): the tier-2 zone
    } catch (const std::logic_error&) {
        // frozen stage/validate rejection: nothing prepared (tier 1)
        r.outcome = Outcome::Failed;
        r.failed = FailedReason::CavityRejected;
        return r;
    } catch (const std::bad_alloc&) {
        // prepare()/commit() allocation failure: the ONLY tier-2 reason
        r.outcome = Outcome::Failed;
        r.failed = FailedReason::CavityThrew;
        return r;
    } catch (...) {
        r.outcome = Outcome::Failed;
        r.failed = FailedReason::CavityThrew;
        return r;
    }
    r.outcome = Outcome::Recovered;
    return r;
}

// COMPLETE whole-set rejection (13.7): EVERY input segment gets a Rejected
// outcome row in `processed` (ascending id). `unattempted` stays EMPTY -
// under 13.7 an entity is EITHER a processed outcome OR an unattempted id,
// never both; a whole-set rejection itemises every entity as Rejected, so
// none remain unattempted. No entity is dropped, none is double-counted.
// Callers invoke this after reserving the report, so the push_backs are
// allocation-free; the wrapper's backstop invokes it inside a nested guard
// for the pathological case where even that reservation was exhausted
// (there it may allocate, hence the guard).
inline void fillWholeSetRejection(RecoveryReport& rep, const ConstraintSet& cs,
                                  SessionStatus why) {
    rep.rejected = true;
    rep.rejectReason = why;
    rep.processed.clear();
    rep.unattempted.clear();          // 13.7: all entities are itemised below
    for (const auto& s : cs.segments)
        rep.processed.push_back(EntityResult{s.id, Outcome::Rejected,
                                             RefusedReason::None,
                                             FailedReason::None});
    std::sort(rep.processed.begin(), rep.processed.end(),
              [](const auto& a, const auto& b) { return a.entity < b.entity; });
}

// ---------------------------------------------------------------------------
// 13.7 - the whole-set entry point with the prefix stopping rule.
// ---------------------------------------------------------------------------
inline void recoverSegmentsImpl(topo::TopologyWorkspace& ws,
                                const ConstraintSet& cs,
                                RecoveryReport& rep) {
    // Reserve the WHOLE report up front, BEFORE any other allocation, so
    // that recording every entity's outcome and the unattempted range is
    // ALLOCATION-FREE thereafter. Every failure handler below can then emit
    // a COMPLETE report (each input entity appears as an outcome or in the
    // unattempted list) without risking a further allocation. The +1 covers
    // the 13.2a rebuild-attributed next-entity row. (If even this first
    // reservation is exhausted, the exception propagates to the wrapper's
    // backstop, which records a structured whole-set rejection - the sole
    // case in which the report cannot be itemised, because nothing at all
    // could be allocated.)
    rep.processed.reserve(cs.segments.size() + 1);
    rep.unattempted.reserve(cs.segments.size());

    // COMPLETE whole-set rejection, alloc-free thanks to the reservation.
    auto rejectWholeSet = [&](SessionStatus why) {
        fillWholeSetRejection(rep, cs, why);
    };

    SessionStatus st{};
    std::optional<RecoverySession> sessOpt;
    try {
        sessOpt = RecoverySession::open(ws, cs, st);
    } catch (const std::bad_alloc&) {
        // Allocation failure during the whole-set preflight: nothing was
        // attempted and nothing changed (13.6 tier 1).
        st = SessionStatus::IndexBuildFailed;
        sessOpt.reset();
    }
    if (!sessOpt) {                       // whole-set Rejected, zero change
        rejectWholeSet(st);               // COMPLETE, structured reason
        return;
    }
    // NO session move: `RecoverySession sess = std::move(*sessOpt)` was the
    // cross-platform leak site (noexcept on libstdc++, allocates on MSVC).
    // Operate through the session already stored in the optional instead.
    RecoverySession& sess = *sessOpt;
    std::vector<ConstraintSegment> segs;
    std::size_t i = 0;
    try {
        segs = sess.input().segments;
        std::sort(segs.begin(), segs.end(),
                  [](const auto& a, const auto& b) { return a.id < b.id; });
    } catch (const std::bad_alloc&) {
        // Could not even enumerate the segment set: a COMPLETE whole-set
        // rejection (every entity accounted for), not a bare flag.
        rejectWholeSet(SessionStatus::IndexBuildFailed);
        return;
    }
    try {
    for (; i < segs.size(); ++i) {
        const auto& seg = segs[i];
        const auto want = sess.orientationOf(seg.chart);
        EntityResult r;
        r.entity = seg.id;
        if (!want) {
            r.outcome = Outcome::Refused;
            r.refused = RefusedReason::UnknownChart;
            rep.processed.push_back(r);
            break;
        }
        std::vector<topo::NodeHandle> chain;
        Routing route;
        try {
            route = routeEntity(ws, sess.input(), sess.index(), seg, chain);
        } catch (const std::bad_alloc&) {
            r.outcome = Outcome::Failed;      // routing: before prepare()
            r.failed = FailedReason::CandidateAllocFailed;   // 13.6 tier 1
            rep.processed.push_back(r);
            break;
        }
        if (route == Routing::AlreadyPresent) {
            r.outcome = Outcome::AlreadyPresent;
            rep.processed.push_back(r);
            continue;                                  // 13.7: continue
        }
        if (route == Routing::PartialChain) {
            r.outcome = Outcome::Refused;
            r.refused = RefusedReason::PartialChain;
            rep.processed.push_back(r);
            break;                                     // 13.7: stop
        }
        RefusedReason why = RefusedReason::None;
        FailedReason fail = FailedReason::None;   // 13.6: candidate-gen fail
        std::optional<RecoveryCandidate> cand;
        try {
            if (route == Routing::Transverse) {
                const auto ia = sess.input().vertices.find(seg.a);
                const auto ib = sess.input().vertices.find(seg.b);
                if (ia == sess.input().vertices.end() ||
                    ib == sess.input().vertices.end()) {
                    r.outcome = Outcome::Refused;
                    r.refused = RefusedReason::PartialChain;
                    rep.processed.push_back(r);
                    break;
                }
                cand = buildTransverseCandidate(ws, sess, seg, ia->second,
                                                ib->second, *want, why, fail);
            } else {
                cand = buildMarkOnlyCandidate(ws, sess, seg, chain, why);
            }
        } catch (const std::bad_alloc&) {
            r.outcome = Outcome::Failed;      // candidate build: before prepare()
            r.failed = FailedReason::CandidateAllocFailed;   // 13.6 tier 1
            rep.processed.push_back(r);
            break;
        }
        if (!cand) {
            if (fail != FailedReason::None) {
                // 13.6 tier 1: an accepted strip we could not triangulate.
                r.outcome = Outcome::Failed;
                r.failed = fail;
            } else {
                r.outcome = Outcome::Refused;
                r.refused = why;
            }
            rep.processed.push_back(r);
            break;                                     // 13.7: stop
        }
        r = commitCandidate(ws, sess, *cand, *want);
        rep.processed.push_back(r);
        if (r.outcome != Outcome::Recovered) break;    // 13.7: stop
        // 13.2a: rebuild for the NEXT entity only; a failure is owned by
        // that next entity, never by the one that just committed.
        bool advanced = true;
        if (i + 1 < segs.size()) {
            try {
                advanced = sess.advanceAfterCommit(ws);
            } catch (const std::bad_alloc&) {
                advanced = false;         // session poisons itself below
            }
        }
        if (i + 1 < segs.size() && !advanced) {
            EntityResult nxt;
            nxt.entity = segs[i + 1].id;
            nxt.outcome = Outcome::Failed;
            nxt.failed = FailedReason::IndexRebuildFailed;
            rep.processed.push_back(nxt);
            ++i;
            break;
        }
    }
    } catch (const std::bad_alloc&) {
        // Residual escape from the loop body. commitCandidate() catches all
        // of its own exceptions, so nothing from the prepare/commit zone can
        // reach here; a residual escape is therefore pre-prepare (tier 1).
        // The report is already reserved, so recording it is alloc-free and
        // the trailing unattempted fill below keeps the report COMPLETE.
        if (i < segs.size()) {
            EntityResult r;
            r.entity = segs[i].id;
            r.outcome = Outcome::Failed;
            r.failed = FailedReason::CandidateAllocFailed;   // 13.6 tier 1
            rep.processed.push_back(r);
        }
    }
    for (std::size_t j = i + 1; j < segs.size(); ++j)
        rep.unattempted.push_back(segs[j].id);
}

// 13.7 - public entry point. recoverSegmentsImpl reserves the whole report
// up front and no longer moves the session out of its optional (the two
// sites that could allocate-and-throw outside a phase guard, one of which
// leaked on MSVC), so every allocation failure it can encounter after that
// reservation is turned into a COMPLETE structured report. The only failure
// the impl cannot pre-reserve for is the pathological exhaustion of that
// FIRST tiny reservation itself. This wrapper is the backstop for exactly
// that residual case: with no prefix recorded it re-attempts a COMPLETE
// itemised whole-set rejection (allocation may have recovered by now), and
// only if THAT too is refused does it fall back to a bare structured
// rejection flag. Either way the five-outcome interface never leaks and, in
// every case where a report can be allocated at all, every input entity is
// accounted for (proven by requireComplete() in probe_bnd).
inline RecoveryReport recoverSegments(topo::TopologyWorkspace& ws,
                                      const ConstraintSet& cs) {
    RecoveryReport rep;
    try {
        recoverSegmentsImpl(ws, cs, rep);
    } catch (const std::bad_alloc&) {
        if (rep.processed.empty()) {
            try {
                fillWholeSetRejection(rep, cs, SessionStatus::IndexBuildFailed);
            } catch (...) {
                // truly nothing can be allocated: bare structured rejection
                rep = RecoveryReport{};
                rep.rejected = true;
                rep.rejectReason = SessionStatus::IndexBuildFailed;
            }
        }
    }
    return rep;
}

}  // namespace paacs::bnd
// BND-1 increment 3 (transverse recovery: valid-ear rank + Failed geometry).
