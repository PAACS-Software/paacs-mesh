// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// BND-1 allocation/rollback probe (standalone; NOT built under sanitizers).
// Design authority: docs/PMBND001_DESIGN.md sec 13 rev 3.
//
//   fresh-workspace-per-k sweep: EVERY allocation of a whole mark-only
//   recovery (session open, PSLG, routing, closure+disk audit, candidate,
//   cavity stage/validate/prepare/commit) fails in turn. After each
//   failure the outcome must be one of the five (never an escaped
//   exception), the logical state must be intact, byte invariants must
//   hold, and a disarmed retry on the SAME workspace must recover.
#include <cstdio>
#include <cstdlib>
#include <map>
#include <new>
#include <set>

#include "paacs/bnd/SegmentRecoveryCommit.hpp"
#include "workspace_digest.hpp"

static long gFailCountdown = 0;
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t n) {
    if (gFailCountdown > 0 && --gFailCountdown == 0) throw std::bad_alloc();
    if (void* p = std::malloc(n ? n : 1)) return p;
    throw std::bad_alloc();
}
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

using namespace paacs;
static int gBad = 0;
#define REQUIRE(c)                                                        \
    do {                                                                  \
        if (!(c)) {                                                       \
            std::printf("FAIL %d: %s\n", __LINE__, #c);                   \
            ++gBad;                                                       \
        }                                                                 \
    } while (0)

struct Fix {
    std::unique_ptr<topo::TopologyWorkspace> ws;
    std::vector<topo::NodeHandle> n;
    std::vector<topo::CellHandle> c;
    bnd::ConstraintSet cs;
};

// COMPLETENESS (13.7): every input entity appears EXACTLY once - either as
// a processed outcome OR in the unattempted list, never both, never twice,
// never dropped, never invented. This checks MULTIPLICITIES and
// DISJOINTNESS (a plain set-union would hide the double-reporting the
// reviewer found). Checked on every recoverSegments() call below.
static void requireComplete(const bnd::RecoveryReport& rep,
                            const bnd::ConstraintSet& cs) {
    std::map<bnd::ConstraintEntityId, int> pc, uc;
    for (const auto& r : rep.processed) pc[r.entity]++;
    for (const auto& u : rep.unattempted) uc[u]++;
    for (const auto& kv : pc) REQUIRE(kv.second == 1);   // no dup in processed
    for (const auto& kv : uc) REQUIRE(kv.second == 1);   // no dup in unattempted
    for (const auto& kv : uc) REQUIRE(pc.find(kv.first) == pc.end());  // disjoint
    std::set<bnd::ConstraintEntityId> want, got;
    for (const auto& s : cs.segments) want.insert(s.id);
    for (const auto& kv : pc) got.insert(kv.first);
    for (const auto& kv : uc) got.insert(kv.first);
    REQUIRE(got == want);          // exactly the input id set (none invented)
    // each entity accounted for EXACTLY once across the two lists
    REQUIRE(rep.processed.size() + rep.unattempted.size() == cs.segments.size());
}

// The 2x2 fan: interior chain node with a full star, so the mark-only
// closure is a genuine disk exercising every code path.
static Fix build() {
    Fix f;
    topo::NodePool np;
    const double C[5][2] = {{0, 0}, {2, 0}, {4, 0}, {2, 2}, {2, -2}};
    for (const auto& q : C) f.n.push_back(np.add({q[0], q[1], 0}, {2, 1}));
    topo::CellPool cp;
    const int T[4][3] = {{0, 1, 3}, {1, 2, 3}, {1, 0, 4}, {2, 1, 4}};
    for (const auto& t : T)
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{f.n[std::size_t(t[0])],
                                          f.n[std::size_t(t[1])],
                                          f.n[std::size_t(t[2])]},
            {2, 5}, np));
    topo::EdgePool ep;
    for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
    f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    f.cs.stamp = bnd::OwnershipStamp::of(*f.ws);
    for (std::size_t i = 0; i < 5; ++i)
        f.cs.params[f.n[i]] = bnd::ParamPoint{{C[i][0], C[i][1]}, 7};
    f.cs.vertices[1] = f.n[0];
    f.cs.vertices[2] = f.n[2];
    f.cs.segments.push_back({10, 1, 2, 7, {1, 9}});
    return f;
}

// Two disjoint quads and TWO entities, BOTH recovered by the TRANSVERSE
// route (each segment crosses its quad's diagonal). recoverSegments()
// performs a rebuild BETWEEN them, so the 13.2a attribution is reached
// with a TRANSVERSE prefix already committed - not the mark-only route.
static Fix build2() {
    Fix f;
    topo::NodePool np;
    const double C[8][2] = {{0, 0},  {2, 0},  {2, 2},  {0, 2},
                            {30, 0}, {32, 0}, {32, 2}, {30, 2}};
    for (const auto& q : C) f.n.push_back(np.add({q[0], q[1], 0}, {2, 1}));
    topo::CellPool cp;
    const int T[4][3] = {{0, 1, 2}, {0, 2, 3},     // quad A, diagonal n0-n2
                         {4, 5, 6}, {4, 6, 7}};     // quad B, diagonal n4-n6
    for (const auto& t : T)
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{f.n[std::size_t(t[0])],
                                          f.n[std::size_t(t[1])],
                                          f.n[std::size_t(t[2])]},
            {2, 5}, np));
    topo::EdgePool ep;
    for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
    f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    f.cs.stamp = bnd::OwnershipStamp::of(*f.ws);
    for (std::size_t i = 0; i < 8; ++i)
        f.cs.params[f.n[i]] = bnd::ParamPoint{{C[i][0], C[i][1]}, 7};
    f.cs.vertices[1] = f.n[1];    // quad A: segment n1 -> n3 crosses n0-n2
    f.cs.vertices[2] = f.n[3];
    f.cs.segments.push_back({10, 1, 2, 7, {1, 9}});
    f.cs.vertices[3] = f.n[5];    // quad B: segment n5 -> n7 crosses n4-n6
    f.cs.vertices[4] = f.n[7];
    f.cs.segments.push_back({20, 3, 4, 7, {1, 9}});
    return f;
}

// A quad whose diagonal the constrained segment CROSSES: exercises the
// initial-cell wedge test, the strip walk, boundary extraction, BOTH
// ear-clipping passes, candidate construction, validation and commit.
static Fix buildTransverse() {
    Fix f;
    topo::NodePool np;
    const double C[4][2] = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
    for (const auto& q : C) f.n.push_back(np.add({q[0], q[1], 0}, {2, 1}));
    topo::CellPool cp;
    const int T[2][3] = {{0, 1, 2}, {0, 2, 3}};
    for (const auto& t : T)
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{f.n[std::size_t(t[0])],
                                          f.n[std::size_t(t[1])],
                                          f.n[std::size_t(t[2])]},
            {2, 5}, np));
    topo::EdgePool ep;
    for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
    f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    f.cs.stamp = bnd::OwnershipStamp::of(*f.ws);
    for (std::size_t i = 0; i < 4; ++i)
        f.cs.params[f.n[i]] = bnd::ParamPoint{{C[i][0], C[i][1]}, 7};
    f.cs.vertices[1] = f.n[1];
    f.cs.vertices[2] = f.n[3];
    f.cs.segments.push_back({10, 1, 2, 7, {1, 9}});
    return f;
}

int main() {
    // Baseline: an undisturbed recovery succeeds.
    {
        Fix f = build();
        f.ws->rebuildAdjacency();
        const auto rep = bnd::recoverSegments(*f.ws, f.cs);
        REQUIRE(!rep.rejected);
        REQUIRE(rep.processed.size() == 1);
        REQUIRE(rep.processed[0].outcome == bnd::Outcome::Recovered);
        requireComplete(rep, f.cs);
    }
    long k = 1, faults = 0, escaped = 0;
    for (;; ++k) {
        Fix f = build();
        f.ws->rebuildAdjacency();
        const auto L0 = paacs_tests::observableLogicalDigest(*f.ws);
        bool faulted = false, recovered = false;
        bnd::RecoveryReport rep;
        gFailCountdown = k;
        try {
            rep = bnd::recoverSegments(*f.ws, f.cs);
            gFailCountdown = 0;
            recovered = !rep.processed.empty() &&
                        rep.processed[0].outcome == bnd::Outcome::Recovered;
        } catch (const std::bad_alloc&) {
            gFailCountdown = 0;
            faulted = true;
            ++escaped;      // an ESCAPE: the five-outcome interface leaked
        }
        gFailCountdown = 0;
        if (!faulted && !recovered) ++faults;   // reported, not escaped
        // COMPLETENESS: even an allocation-failed run reports every entity
        if (!faulted) requireComplete(rep, f.cs);
        // byte invariants hold either way
        {
            const auto m = measureWorkspaceMemory(*f.ws);
            REQUIRE(m.consistent());
            REQUIRE(f.ws->nodes().stats().bytesConsistent());
            REQUIRE(f.ws->cells().stats().bytesConsistent());
            REQUIRE(f.ws->edges().stats().bytesConsistent());
            REQUIRE(f.ws->faces().stats().bytesConsistent());
        }
        if (!recovered) {
            // nothing committed -> the LOGICAL state is intact
            REQUIRE(paacs_tests::observableLogicalDigest(*f.ws) == L0);
            // disarmed retry on the SAME workspace succeeds
            auto cs2 = f.cs;
            cs2.stamp = bnd::OwnershipStamp::of(*f.ws);
            const auto again = bnd::recoverSegments(*f.ws, cs2);
            REQUIRE(!again.processed.empty());
            REQUIRE(again.processed[0].outcome == bnd::Outcome::Recovered);
            f.ws->rebuildAdjacency();
        } else {
            break;          // k exceeded the allocation count
        }
    }
    // ---- two-entity sweep: reaches the BETWEEN-ENTITY index rebuild ----
    long k2 = 1, both = 0, partial = 0, escaped2 = 0, rebuildFailed = 0;
    for (;; ++k2) {
        Fix f = build2();
        f.ws->rebuildAdjacency();
        bnd::RecoveryReport rep;
        bool faulted = false;
        gFailCountdown = k2;
        try {
            rep = bnd::recoverSegments(*f.ws, f.cs);
        } catch (const std::bad_alloc&) {
            ++escaped2;
            faulted = true;
        }
        gFailCountdown = 0;
        if (faulted) continue;
        requireComplete(rep, f.cs);   // both entities always accounted for
        std::size_t rec = 0;
        bool sawRebuild = false;
        for (const auto& r : rep.processed) {
            if (r.outcome == bnd::Outcome::Recovered) ++rec;
            if (r.failed == bnd::FailedReason::IndexRebuildFailed) {
                sawRebuild = true;
                ++rebuildFailed;
                // 13.2a EXACT result vector: entity 10 stays Recovered,
                // entity 20 carries the rebuild Failed, nothing else.
                REQUIRE(rep.processed.size() == 2);
                REQUIRE(rep.processed[0].entity == 10);
                REQUIRE(rep.processed[0].outcome == bnd::Outcome::Recovered);
                REQUIRE(rep.processed[1].entity == 20);
                REQUIRE(rep.processed[1].outcome == bnd::Outcome::Failed);
                REQUIRE(rep.unattempted.empty());   // 20 WAS attempted
            }
        }
        if (rec == 2) {
            ++both;
            break;                        // k exceeded the allocations
        }
        if (rec == 1) {
            ++partial;
            // Entity 10's TRANSVERSE recovery committed: the crossing
            // segment n1-n3 now EXISTS as a required, classified edge, and
            // the crossed diagonal n0-n2 is gone. A later failure must
            // preserve that transverse prefix EXACTLY.
            f.ws->rebuildAdjacency();
            const auto a = f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[3]);
            REQUIRE(bool(a));
            REQUIRE(f.ws->edges().get(*a).requiredConstraint);
            REQUIRE(f.ws->edges().get(*a).cls.dim == 1);
            REQUIRE(f.ws->edges().get(*a).cls.tag == 9);
            REQUIRE(!f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[2]));
            // entity 20's quad is untouched: its crossing segment n5-n7 is
            // NOT yet an edge and its diagonal n4-n6 still stands.
            REQUIRE(!f.ws->edges().find(f.ws->nodes(), f.n[5], f.n[7]));
            REQUIRE(bool(f.ws->edges().find(f.ws->nodes(), f.n[4], f.n[6])));
            if (sawRebuild) {
                // the rebuild-attributed case must ALSO show the prefix mark
                REQUIRE(f.ws->edges().get(*a).requiredConstraint);
            }
        }
    }
    // ---- TRANSVERSE sweep: every allocation of the crossing route ----
    // Each fault is classified by tier (13.6): tier 1 (candidate / stage /
    // validator rejection - NOTHING prepared) must leave the FULL
    // observable digest unchanged; tier 2 (preparation / allocation
    // failure) may leave permitted capacity residue, so only the LOGICAL
    // digest is guaranteed equal. We check the STRONGER invariant for the
    // tier-1 cases rather than the logical digest generically, and require
    // BOTH tiers to be exercised by the seam battery.
    long k3 = 1, tFaults = 0, tEscaped = 0, tier1 = 0, tier2 = 0;
    for (;; ++k3) {
        Fix f = buildTransverse();
        f.ws->rebuildAdjacency();
        const auto L0 = paacs_tests::observableLogicalDigest(*f.ws);
        const auto D0 = paacs_tests::observableDigest(*f.ws);   // FULL
        bnd::RecoveryReport rep;
        bool faulted = false;
        gFailCountdown = k3;
        try {
            rep = bnd::recoverSegments(*f.ws, f.cs);
        } catch (const std::bad_alloc&) {
            ++tEscaped;
            faulted = true;
        }
        gFailCountdown = 0;
        if (faulted) continue;
        requireComplete(rep, f.cs);   // entity accounted for on every fault
        const bool rec = !rep.processed.empty() &&
                         rep.processed[0].outcome == bnd::Outcome::Recovered;
        {
            const auto m = measureWorkspaceMemory(*f.ws);
            REQUIRE(m.consistent());
            REQUIRE(f.ws->nodes().stats().bytesConsistent());
            REQUIRE(f.ws->cells().stats().bytesConsistent());
            REQUIRE(f.ws->edges().stats().bytesConsistent());
            REQUIRE(f.ws->faces().stats().bytesConsistent());
        }
        if (!rec) {
            ++tFaults;
            // PHASE-attribute the fault from the REPORTED reason, not from a
            // digest coincidence (a tier-2 failure can also leave no
            // residue, so full-digest equality does not by itself identify
            // the phase). After the reason phase-split, CavityThrew is the
            // ONLY reason raised from the prepare()/commit() zone (13.6
            // tier 2). EVERY other outcome prepared NOTHING (tier 1):
            // Rejected (whole-set), Refused, StampChainBroken,
            // InvalidCandidateGeometry, CavityRejected, and
            // CandidateAllocFailed (routing / candidate build / the
            // stage+validate zone, all BEFORE prepare()). We then assert the
            // digest invariant appropriate to THAT phase.
            const bool preparePhase =
                !rep.processed.empty() &&
                rep.processed[0].outcome == bnd::Outcome::Failed &&
                rep.processed[0].failed == bnd::FailedReason::CavityThrew;
            if (preparePhase) {
                ++tier2;   // tier 2: logical digest intact, residue permitted
                REQUIRE(paacs_tests::observableLogicalDigest(*f.ws) == L0);
            } else {
                ++tier1;   // tier 1: FULL observable digest unchanged
                REQUIRE(paacs_tests::observableDigest(*f.ws) == D0);
            }
            auto cs2 = f.cs;
            cs2.stamp = bnd::OwnershipStamp::of(*f.ws);
            const auto again = bnd::recoverSegments(*f.ws, cs2);
            REQUIRE(!again.processed.empty());
            REQUIRE(again.processed[0].outcome == bnd::Outcome::Recovered);
            f.ws->rebuildAdjacency();
            // the recovered segment really exists and is marked
            const auto e =
                f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[3]);
            REQUIRE(bool(e));
            REQUIRE(f.ws->edges().get(*e).requiredConstraint);
        } else {
            break;
        }
    }
    std::printf("transverse tiers: %ld tier1 (full-digest), %ld tier2 "
                "(logical-digest)\n",
                tier1, tier2);
    REQUIRE(tier1 >= 1);   // the early candidate-side faults are tier 1
    REQUIRE(tier2 >= 1);   // the prepare-side faults are tier 2
    std::printf("bnd sweep: %ld reported failures, %ld escaped exceptions\n",
                faults, escaped);
    std::printf("transverse: %ld reported failures, %ld escapes\n",
                tFaults, tEscaped);
    REQUIRE(tEscaped == 0);
    REQUIRE(tFaults >= 10);
    std::printf("two-entity: %ld partial-prefix, %ld rebuild-attributed, "
                "%ld escapes\n",
                partial, rebuildFailed, escaped2);
    REQUIRE(escaped2 == 0);
    REQUIRE(partial >= 1);   // a later-entity failure preserved the prefix
    REQUIRE(faults >= 10);
    REQUIRE(escaped == 0);   // the five-outcome interface never leaks
    if (gBad) {
        std::printf("BND PROBE FAILED (%d)\n", gBad);
        return 1;
    }
    std::printf("BND PROBE PASSED\n");
    return 0;
}
