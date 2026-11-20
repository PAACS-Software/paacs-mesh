// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// ============================================================================
// BND-1 candidate-generation slice, increment 1 (VALIDATION ONLY).
// Design authority: docs/PMBND001_DESIGN.md sec 13 rev 3 (ACCEPTED).
//
// This header implements sec 13.1 (whole-set PSLG preflight), 13.1a
// (injective vertex mapping), and 13.2/13.2a (the BND-owned
// RecoverySession stamp chain + orientationByChart). It contains NO
// mutation path: no cavity transaction, no walk, no ear clipping. Those
// arrive in increment 2 under the same accepted scope.
//
// Everything is chart-local and exact (sec 5a/5b): comparisons use
// predicates::orient2d over the AUTHORITATIVE stored ConstraintSet
// coordinates, never re-derived, never workspace Vec3. Full-handle
// identity throughout (the A1 discipline).
// ============================================================================
#include <algorithm>
#include <map>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "paacs/bnd/ParamGeometry.hpp"
#include "paacs/bnd/SegmentRecovery.hpp"

namespace paacs::bnd {

// ---------------------------------------------------------------------------
// 13.1 / 13.1a - whole-set PSLG preflight over the ConstraintSet ITSELF.
// (A1 validates each segment against the MESH; this validates the
// constraint set against itself.) Rejected => zero revision change,
// which is structural here: this code never mutates.
// ---------------------------------------------------------------------------
enum class PslgStatus {
    Ok,
    NonInjectiveVertexMapping,  // 13.1a: two ids -> one handle
    DuplicateSegment,
    CoincidentVertices,         // distinct handles, identical coords
    ConstraintCrossing,         // proper transversal, no shared endpoint
    ConstraintTJunction,        // vertex strictly interior to a segment
    CollinearOverlap,           // collinear, sharing more than a point
};

struct PslgIssue {
    PslgStatus status{};
    std::optional<ConstraintEntityId> entityA, entityB;
    std::optional<topo::NodeHandle> node;
};

namespace detail {

struct SegGeom {
    ConstraintEntityId id;
    ChartId chart;
    topo::NodeHandle na, nb;
    Vec2 pa, pb;
};

// Exact: do collinear segments overlap in more than a single point?
// Precondition: all four endpoints are collinear (checked by caller).
inline bool collinearOverlapMoreThanPoint(const SegGeom& s, const SegGeom& t) {
    // Project onto the dominant axis of s and compare closed intervals.
    const int axis = detail::dominantAxis(s.pa, s.pb);
    auto c = [&](const Vec2& p) { return detail::axisCoord(p, axis); };
    const double s0 = std::min(c(s.pa), c(s.pb)), s1 = std::max(c(s.pa), c(s.pb));
    const double t0 = std::min(c(t.pa), c(t.pb)), t1 = std::max(c(t.pa), c(t.pb));
    const double lo = std::max(s0, t0), hi = std::min(s1, t1);
    return lo < hi;                 // strict: touching at one point is fine
}

}  // namespace detail

inline std::vector<PslgIssue> preflightPslg(const ConstraintSet& cs) {
    std::vector<PslgIssue> out;

    // 13.1a: the ConstraintVertexId -> NodeHandle mapping must be INJECTIVE.
    {
        std::unordered_map<topo::NodeHandle, ConstraintVertexId> seenHandle;
        std::vector<std::pair<ConstraintVertexId, topo::NodeHandle>> sorted(
            cs.vertices.begin(), cs.vertices.end());
        std::sort(sorted.begin(), sorted.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        for (const auto& [vid, h] : sorted)
            if (!seenHandle.emplace(h, vid).second)
                out.push_back({PslgStatus::NonInjectiveVertexMapping, vid,
                               seenHandle[h], h});
    }

    // Resolve segments to exact geometry; skip any whose endpoints or
    // parameters are missing (A1 owns those diagnostics).
    std::vector<detail::SegGeom> S;
    S.reserve(cs.segments.size());
    for (const auto& seg : cs.segments) {
        const auto ia = cs.vertices.find(seg.a);
        const auto ib = cs.vertices.find(seg.b);
        if (ia == cs.vertices.end() || ib == cs.vertices.end()) continue;
        const ParamPoint* pa = cs.paramOf(ia->second);
        const ParamPoint* pb = cs.paramOf(ib->second);
        if (!pa || !pb) continue;
        S.push_back({seg.id, seg.chart, ia->second, ib->second, pa->p, pb->p});
    }
    std::sort(S.begin(), S.end(), [](const auto& a, const auto& b) {
        return a.id < b.id;                    // deterministic iteration
    });

    // Duplicate segments: same canonical endpoint pair (full handles).
    {
        std::map<std::pair<topo::NodeHandle, topo::NodeHandle>,
                 ConstraintEntityId>
            seen;
        for (const auto& s : S) {
            auto key = s.na < s.nb ? std::make_pair(s.na, s.nb)
                                   : std::make_pair(s.nb, s.na);
            const auto it = seen.find(key);
            if (it != seen.end())
                out.push_back({PslgStatus::DuplicateSegment, s.id, it->second,
                               {}});
            else
                seen.emplace(key, s.id);
        }
    }

    // Coincident constraint VERTICES: distinct handles, identical stored
    // coordinates within one chart.
    {
        std::vector<std::pair<topo::NodeHandle, ParamPoint>> V(
            cs.params.begin(), cs.params.end());
        std::sort(V.begin(), V.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        std::unordered_set<topo::NodeHandle> constrained;
        for (const auto& [vid, h] : cs.vertices) constrained.insert(h);
        for (std::size_t i = 0; i < V.size(); ++i)
            for (std::size_t j = i + 1; j < V.size(); ++j) {
                if (!constrained.count(V[i].first) ||
                    !constrained.count(V[j].first))
                    continue;                   // only CONSTRAINT vertices
                if (V[i].second.chart != V[j].second.chart) continue;
                if (V[i].second.p.u == V[j].second.p.u &&
                    V[i].second.p.v == V[j].second.p.v)
                    out.push_back({PslgStatus::CoincidentVertices, {}, {},
                                   V[j].first});
            }
    }

    // Pairwise segment tests, chart-local, exact.
    for (std::size_t i = 0; i < S.size(); ++i) {
        for (std::size_t j = i + 1; j < S.size(); ++j) {
            const auto& s = S[i];
            const auto& t = S[j];
            if (s.chart != t.chart) continue;   // sec 5b
            const bool shares = s.na == t.na || s.na == t.nb ||
                                s.nb == t.na || s.nb == t.nb;
            const auto d1 = predicates::orient2d(s.pa, s.pb, t.pa);
            const auto d2 = predicates::orient2d(s.pa, s.pb, t.pb);
            const auto d3 = predicates::orient2d(t.pa, t.pb, s.pa);
            const auto d4 = predicates::orient2d(t.pa, t.pb, s.pb);
            // collinear overlap (all four orientations zero)
            if (d1 == predicates::Sign::Zero && d2 == predicates::Sign::Zero &&
                d3 == predicates::Sign::Zero && d4 == predicates::Sign::Zero) {
                if (detail::collinearOverlapMoreThanPoint(s, t))
                    out.push_back({PslgStatus::CollinearOverlap, s.id, t.id,
                                   {}});
                continue;
            }
            if (!shares) {
                // proper transversal crossing
                const bool str1 = d1 != d2 && d1 != predicates::Sign::Zero &&
                                  d2 != predicates::Sign::Zero;
                const bool str2 = d3 != d4 && d3 != predicates::Sign::Zero &&
                                  d4 != predicates::Sign::Zero;
                if (str1 && str2)
                    out.push_back({PslgStatus::ConstraintCrossing, s.id, t.id,
                                   {}});
            }
            // unsplit T-junction: an endpoint of one strictly interior to
            // the other (excluding shared endpoints).
            auto tj = [&](const detail::SegGeom& host,
                          const detail::SegGeom& guest,
                          const topo::NodeHandle& gh, const Vec2& gp) {
                if (gh == host.na || gh == host.nb) return;
                if (predicates::orient2d(host.pa, host.pb, gp) !=
                    predicates::Sign::Zero)
                    return;
                const int axis = detail::dominantAxis(host.pa, host.pb);
                const double c0 = detail::axisCoord(host.pa, axis),
                             c1 = detail::axisCoord(host.pb, axis),
                             cg = detail::axisCoord(gp, axis);
                if (cg > std::min(c0, c1) && cg < std::max(c0, c1))
                    out.push_back({PslgStatus::ConstraintTJunction, host.id,
                                   guest.id, gh});
            };
            tj(s, t, t.na, t.pa);
            tj(s, t, t.nb, t.pb);
            tj(t, s, s.na, s.pa);
            tj(t, s, s.nb, s.pb);
        }
    }

    // Deterministic total order (sec 12 discipline).
    std::stable_sort(out.begin(), out.end(),
                     [](const PslgIssue& a, const PslgIssue& b) {
                         if (a.status != b.status)
                             return static_cast<int>(a.status) <
                                    static_cast<int>(b.status);
                         if (a.entityA != b.entityA) return a.entityA < b.entityA;
                         if (a.entityB != b.entityB) return a.entityB < b.entityB;
                         const bool an = a.node.has_value(),
                                    bn = b.node.has_value();
                         if (an != bn) return an < bn;
                         if (an && !(*a.node == *b.node)) return *a.node < *b.node;
                         return false;
                     });
    return out;
}

// ---------------------------------------------------------------------------
// 13.2 / 13.2a - the BND-owned RecoverySession.
// ---------------------------------------------------------------------------
enum class SessionStatus {
    Ok,
    InputStampMismatch,     // caller's ConstraintSet stamp is not current
    A1Rejected,
    A2Rejected,
    PslgRejected,
    IndexBuildFailed,
};

#if defined(PAACS_TOPO_TESTING)
namespace detail {
inline bool& recoverySessionIndexBuildActive() noexcept {
    static thread_local bool active = false;
    return active;
}
struct RecoverySessionIndexBuildMarker {
    RecoverySessionIndexBuildMarker() noexcept {
        recoverySessionIndexBuildActive() = true;
    }
    ~RecoverySessionIndexBuildMarker() {
        recoverySessionIndexBuildActive() = false;
    }
};
}  // namespace detail
#endif

class RecoverySession {
  public:
    // Opens a session: validates the caller's stamp ONCE, runs A1, A2 and
    // the PSLG preflight, derives orientationByChart, builds the validated
    // index, and takes an IMMUTABLE SESSION-OWNED SNAPSHOT of the input
    // (13.2). The caller's ConstraintSet is never modified and may be
    // destroyed afterwards: all later candidate generation uses the
    // snapshot, never the caller's object.
    static std::optional<RecoverySession> open(
        const topo::TopologyWorkspace& ws, const ConstraintSet& cs,
        SessionStatus& status) {
        if (!cs.stamp.matches(ws)) {
            status = SessionStatus::InputStampMismatch;
            return std::nullopt;
        }
        if (!preflightConstraintSet(ws, cs).empty()) {
            status = SessionStatus::A1Rejected;
            return std::nullopt;
        }
        if (!preflightSourceEmbedding(ws, cs).empty()) {
            status = SessionStatus::A2Rejected;
            return std::nullopt;
        }
        if (!preflightPslg(cs).empty()) {
            status = SessionStatus::PslgRejected;
            return std::nullopt;
        }
#if defined(PAACS_TOPO_TESTING)
        detail::RecoverySessionIndexBuildMarker indexBuildMarker;
#endif
        auto idx = buildTriAdjacency(ws);
        if (!idx) {
            status = SessionStatus::IndexBuildFailed;
            return std::nullopt;
        }
        RecoverySession s;
        s.input_ = cs;                    // IMMUTABLE snapshot (13.2)
        s.callerStamp_ = cs.stamp;        // preserved UNCHANGED
        s.stamp_ = OwnershipStamp::of(ws);   // evolving, separate
        s.index_ = std::move(*idx);
        s.orientationByChart_ = deriveOrientation(ws, s.input_);
        s.valid_ = true;
        status = SessionStatus::Ok;
        return s;
    }

    // The session-owned authoritative input. All candidate generation
    // MUST read this, never the caller's (possibly destroyed) object.
    const ConstraintSet& input() const noexcept { return input_; }
    // The caller's ORIGINAL stamp, preserved unchanged for the record.
    const OwnershipStamp& callerStamp() const noexcept { return callerStamp_; }
    // The EVOLVING session stamp.
    const OwnershipStamp& stamp() const noexcept { return stamp_; }
    const TriAdjacencyIndex& index() const noexcept { return index_; }
    bool valid() const noexcept { return valid_; }
    std::optional<predicates::Sign> orientationOf(ChartId c) const {
        const auto it = orientationByChart_.find(c);
        return it == orientationByChart_.end()
                   ? std::nullopt
                   : std::optional<predicates::Sign>{it->second};
    }
    std::size_t chartCount() const noexcept {
        return orientationByChart_.size();
    }

    // 13.2: the four-way identity required before EVERY handoff. A
    // POISONED session (see advanceAfterCommit) can never allow a handoff.
    bool handoffAllowed(const topo::TopologyWorkspace& ws,
                        const OwnershipStamp& candidate) const {
        if (!valid_) return false;
        const auto now = OwnershipStamp::of(ws);
        auto eq = [](const OwnershipStamp& x, const OwnershipStamp& y) {
            return x.nodePool == y.nodePool && x.cellPool == y.cellPool &&
                   x.edgePool == y.edgePool && x.facePool == y.facePool &&
                   x.revision == y.revision;
        };
        return eq(candidate, index_.stamp) && eq(index_.stamp, stamp_) &&
               eq(stamp_, now);
    }

    // 13.2a: advance AFTER a commit. SECURITY CONTRACT (review): this is
    // NOT a generic resynchronizer. Adopting a changed workspace REQUIRES
    // re-establishing the validation the session was opened on - A1 and A2
    // against the session's OWN snapshot - because an arbitrary external
    // mutation may have invalidated either. (The PSLG result cannot change:
    // it is a property of the immutable snapshot alone.) If revalidation or
    // the index rebuild fails, the session is POISONED: it returns false,
    // valid_ becomes false permanently, and handoffAllowed() can never be
    // true again. The caller attributes that Failed to the NEXT entity
    // (13.2a), never to the one that just committed.
    // NOTE (review): this method ALLOCATES (snapshot copy, A1/A2, maps).
    // Any allocation failure POISONS the session and reports false rather
    // than escaping the five-outcome interface; the caller attributes the
    // Failed to the NEXT entity (13.2a).
    bool advanceAfterCommit(const topo::TopologyWorkspace& ws) noexcept {
        try {
            return advanceImpl(ws);
        } catch (...) {
            valid_ = false;
            return false;
        }
    }

  private:
    bool advanceImpl(const topo::TopologyWorkspace& ws) {
        if (!valid_) return false;
        // Validate the snapshot against the CURRENT workspace. The
        // snapshot keeps the caller's original stamp by contract, so the
        // validation view carries the current one.
        ConstraintSet view = input_;
        view.stamp = OwnershipStamp::of(ws);
        if (!preflightConstraintSet(ws, view).empty()) {
            valid_ = false;
            return false;
        }
        if (!preflightSourceEmbedding(ws, view).empty()) {
            valid_ = false;
            return false;
        }
        auto idx = buildTriAdjacency(ws);
        if (!idx) {
            valid_ = false;
            return false;
        }
        index_ = std::move(*idx);
        stamp_ = OwnershipStamp::of(ws);
        orientationByChart_ = deriveOrientation(ws, view);
        return true;
    }

  public:
    RecoverySession(const RecoverySession&) = default;
    RecoverySession(RecoverySession&&) = default;
    RecoverySession& operator=(RecoverySession&&) = default;

  private:
    RecoverySession() = default;

    // 13.4a: A2 returns no orientation artifact, so the SESSION derives
    // and retains the accepted per-chart sign.
    static std::unordered_map<ChartId, predicates::Sign> deriveOrientation(
        const topo::TopologyWorkspace& ws, const ConstraintSet& cs) {
        std::unordered_map<ChartId, predicates::Sign> m;
        for (const auto& ch : ws.cells().activeBySequence()) {
            const auto cd = ws.cells().get(ch);
            const ParamPoint* p0 = cs.paramOf(cd.nodes[0]);
            const ParamPoint* p1 = cs.paramOf(cd.nodes[1]);
            const ParamPoint* p2 = cs.paramOf(cd.nodes[2]);
            if (!p0 || !p1 || !p2) continue;
            m.emplace(p0->chart,
                      predicates::orient2d(p0->p, p1->p, p2->p));
        }
        return m;
    }

    ConstraintSet input_;                 // immutable session-owned snapshot
    OwnershipStamp callerStamp_;          // caller's original, unchanged
    OwnershipStamp stamp_;                // evolving session stamp
    TriAdjacencyIndex index_;
    std::unordered_map<ChartId, predicates::Sign> orientationByChart_;
    bool valid_ = false;
};

}  // namespace paacs::bnd
