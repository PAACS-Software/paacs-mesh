// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// BND-2A / ASCR probes (standalone; NOT built under sanitizers):
//   1. no_allocation_in_phase3 - global operator new counts every
//      allocation while detail::phase3Active() is set; the COMPOSITE
//      Phase-3 (commitNodeBatch + commitCavityPrepared) must count ZERO
//      for BOTH a single-Steiner fan AND a multi-Steiner batch, proving
//      the atomic Steiner commit is allocation-free on ANY stdlib.
//   2. prepare-side allocation sweep - every allocation of begin/stage/
//      prepare fails in turn on a FRESH workspace per injected index k;
//      after each failure the LOGICAL state is restored, byte invariants
//      hold, and a disarmed retry commits.
//   3. commit-zone allocation sweep - faults inside commit()'s throwing
//      zone (BEFORE Phase 3); the transaction rolls back cleanly.
// Fixtures: a quad disk with an EXHAUSTED-generation removed cell (mixing
// fresh materialization + in-commit reuse + incidence-erase), retriangulated
// EITHER into a single central fan (1 node, 4 cells) OR into a two-Steiner
// disk with a staged-staged edge (2 nodes, 6 cells), so the batch path is
// genuinely exercised - not only the single-center fan.
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <new>

#include "paacs/topo/SteinerCavity.hpp"
#include "paacs/bnd/SteinerRecovery.hpp"
#include "workspace_digest.hpp"

static long gPhase3Allocs = 0;
static long gFailCountdown = 0;
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t n) {
    if (paacs::topo::detail::phase3Active()) ++gPhase3Allocs;
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

using namespace paacs::topo;
namespace bnd = paacs::bnd;
using paacs::Vec2;
using paacs::Vec3;
static int gBad = 0;
#define REQUIRE(c)                                                        \
    do {                                                                  \
        if (!(c)) {                                                       \
            std::printf("FAIL %d: %s\n", __LINE__, #c);                   \
            ++gBad;                                                       \
        }                                                                 \
    } while (0)

struct Fix {
    std::unique_ptr<TopologyWorkspace> ws;
    std::vector<NodeHandle> n;
    std::vector<CellHandle> c;
    WorkspaceStamp stamp;
};
// Quad disk; the FIRST removed triangle is generation-exhausted.
static Fix buildQuad() {
    Fix f;
    NodePool np;
    f.n.push_back(np.add({0, 0, 0}, {2, 1}));
    f.n.push_back(np.add({2, 0, 0}, {2, 1}));
    f.n.push_back(np.add({2, 2, 0}, {2, 1}));
    f.n.push_back(np.add({0, 2, 0}, {2, 1}));
    CellPool cp;
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[1], f.n[2]}, {2, 5}, np));
    f.c[0] = paacs::topo::detail::TransactionAccessTestProbe::
        exhaustLiveCell(cp, f.c[0]);
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[2], f.n[3]}, {2, 5}, np));
    EdgePool ep;
    for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
    f.ws = std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), FacePool{});
    f.stamp = WorkspaceStamp::of(*f.ws);
    return f;
}
static SteinerCavitySpec fanSpec(const Fix& f) {
    SteinerCavitySpec s;
    s.stamp = f.stamp;
    s.removed = f.c;
    s.steinerNodes.push_back({0, {1, 1, 0}, {2, 1}});
    auto tri = [&](std::uint32_t a, std::uint32_t b) {
        SteinerCavityCell c;
        c.nodes = {SteinerNodeRef::Existing(f.n[a]),
                   SteinerNodeRef::Existing(f.n[b]), SteinerNodeRef::Staged(0)};
        c.cls = {2, 5};
        return c;
    };
    s.added = {tri(0, 1), tri(1, 2), tri(2, 3), tri(3, 0)};
    s.requirements.push_back(
        {SteinerNodeRef::Existing(f.n[0]), SteinerNodeRef::Staged(0), {1, 9}});
    return s;
}
static SteinerCavitySpec multiSpec(const Fix& f) {
    SteinerCavitySpec s;
    s.stamp = f.stamp;
    s.removed = f.c;
    s.steinerNodes.push_back({0, {0.7, 1, 0}, {2, 1}});
    s.steinerNodes.push_back({1, {1.3, 1, 0}, {2, 1}});
    auto A = SteinerNodeRef::Staged(0), B = SteinerNodeRef::Staged(1);
    auto E = [&](std::uint32_t i) { return SteinerNodeRef::Existing(f.n[i]); };
    auto cell = [&](SteinerNodeRef a, SteinerNodeRef b, SteinerNodeRef c) {
        SteinerCavityCell x; x.nodes = {a, b, c}; x.cls = {2, 5}; return x;
    };
    s.added = {cell(E(0), E(1), A), cell(E(1), B, A), cell(E(1), E(2), B),
               cell(E(2), E(3), B), cell(E(3), A, B), cell(E(3), E(0), A)};
    s.requirements.push_back({A, B, {1, 9}});
    return s;
}

using SpecFn = std::function<SteinerCavitySpec(const Fix&)>;

static void runSweeps(const char* label, const SpecFn& mk,
                      std::size_t expectCells, std::size_t expectNodes,
                      long& faultsOut, long& faults2Out) {
    // 1. no allocation in composite Phase 3
    {
        Fix f = buildQuad();
        f.ws->rebuildAdjacency();
        const auto spec = mk(f);
        auto tx = f.ws->beginSteinerCavityTransaction();
        (void)tx.stage(spec);
        tx.prepare();
        gPhase3Allocs = 0;
        const auto out = tx.commit();
        REQUIRE(gPhase3Allocs == 0);
        REQUIRE(out.createdNodes.size() == expectNodes);
        REQUIRE(out.resultCells.size() == expectCells);
        f.ws->rebuildAdjacency();
    }
    // 2. prepare-side allocation sweep
    long k = 1, faults = 0;
    for (;; ++k) {
        Fix f = buildQuad();
        f.ws->rebuildAdjacency();
        const auto spec = mk(f);
        const auto L0 = paacs_tests::observableLogicalDigest(*f.ws);
        bool faulted = false;
        {
            gFailCountdown = k;
            try {
                auto tx = f.ws->beginSteinerCavityTransaction();
                (void)tx.stage(spec);
                tx.prepare();
                gFailCountdown = 0;
                tx.abandon();
            } catch (const std::bad_alloc&) {
                faulted = true;
            }
            gFailCountdown = 0;
        }
        REQUIRE(paacs_tests::observableLogicalDigest(*f.ws) == L0);
        {
            const auto m = measureWorkspaceMemory(*f.ws);
            REQUIRE(m.consistent());
            REQUIRE(f.ws->nodes().stats().bytesConsistent());
            REQUIRE(f.ws->cells().stats().bytesConsistent());
            REQUIRE(f.ws->edges().stats().bytesConsistent());
        }
        {
            auto tx = f.ws->beginSteinerCavityTransaction();
            (void)tx.stage(mk(f));
            const auto out = tx.commit();
            REQUIRE(out.resultCells.size() == expectCells);
            REQUIRE(out.createdNodes.size() == expectNodes);
            f.ws->rebuildAdjacency();
            REQUIRE(f.ws->adjacencyReady());
        }
        if (faulted) ++faults; else break;
    }
    // 3. commit-zone allocation sweep
    long k2 = 1, faults2 = 0;
    for (;; ++k2) {
        Fix f = buildQuad();
        f.ws->rebuildAdjacency();
        const auto spec = mk(f);
        const auto L0 = paacs_tests::observableLogicalDigest(*f.ws);
        bool faulted = false;
        {
            auto tx = f.ws->beginSteinerCavityTransaction();
            (void)tx.stage(spec);
            tx.prepare();
            gFailCountdown = k2;
            try {
                const auto out = tx.commit();
                gFailCountdown = 0;
                REQUIRE(out.resultCells.size() == expectCells);
            } catch (const std::bad_alloc&) {
                faulted = true;
                gFailCountdown = 0;
                tx.abandon();
            }
        }
        if (faulted) {
            ++faults2;
            REQUIRE(paacs_tests::observableLogicalDigest(*f.ws) == L0);
            const auto m = measureWorkspaceMemory(*f.ws);
            REQUIRE(m.consistent());
            REQUIRE(f.ws->nodes().stats().bytesConsistent());
            REQUIRE(f.ws->cells().stats().bytesConsistent());
            REQUIRE(f.ws->edges().stats().bytesConsistent());
            {
                auto tx = f.ws->beginSteinerCavityTransaction();
                (void)tx.stage(mk(f));
                (void)tx.commit();
                f.ws->rebuildAdjacency();
            }
        } else {
            break;
        }
    }
    std::printf("%s: %ld prepare-side, %ld commit-zone faults\n", label,
                faults, faults2);
    faultsOut = faults;
    faults2Out = faults2;
}

// ---- public BND single-Steiner ROUTE sweep -------------------------------
struct MockProv : bnd::GeometryProvider {
    bnd::GeometryProviderStamp st{42, 1};
    bnd::GeometryProviderStamp stamp() const override { return st; }
    Vec3 pointAt(bnd::ChartId, const Vec2& uv) const override {
        return {uv.u, uv.v, 0};
    }
};
struct RouteFix {
    std::unique_ptr<TopologyWorkspace> ws;
    std::vector<NodeHandle> n;
    bnd::ConstraintSet cs;
    bnd::ConstraintSegment seg;
};
static RouteFix buildRouteFix() {
    RouteFix f;
    NodePool np;
    f.n.push_back(np.add({0, 0, 0}, {2, 1}));
    f.n.push_back(np.add({2, 0, 0}, {2, 1}));
    f.n.push_back(np.add({2, 2, 0}, {2, 1}));
    f.n.push_back(np.add({0, 2, 0}, {2, 1}));
    CellPool cp;
    cp.add(ElementType::Tri3,
           std::vector<NodeHandle>{f.n[0], f.n[1], f.n[2]}, {2, 5}, np);
    cp.add(ElementType::Tri3,
           std::vector<NodeHandle>{f.n[0], f.n[2], f.n[3]}, {2, 5}, np);
    EdgePool ep;
    for (const auto& ch : cp.activeBySequence()) ep.insertCellEdges(np, cp, ch);
    f.ws = std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), FacePool{});
    const Vec2 uv[4] = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
    for (int i = 0; i < 4; ++i) {
        f.cs.vertices[100 + std::uint64_t(i)] = f.n[std::size_t(i)];
        f.cs.params[f.n[std::size_t(i)]] = bnd::ParamPoint{uv[i], 7};
    }
    f.seg.id = 500;
    f.seg.a = 101; f.seg.b = 103; f.seg.chart = 7; f.seg.cls = {1, 9};
    f.cs.segments.push_back(f.seg);
    f.cs.stamp = bnd::OwnershipStamp::of(*f.ws);
    return f;
}
static bool isTier2(bnd::SteinerFail f) {
    return f == bnd::SteinerFail::StagedNodeAllocFailed ||
           f == bnd::SteinerFail::CavityThrew;
}
// Multi-Steiner fan (K=3): segment A-B crosses P1-Q, P2-Q, P3-Q -> chain of 5
// (A, s1, s2, s3, B). Exercises the multi-cell strip + multi-staged commit.
static RouteFix buildMultiRouteFix() {
    RouteFix f;
    NodePool np;
    NodeHandle A = np.add({0, 0, 0}, {2, 1});
    NodeHandle P1 = np.add({1, 1, 0}, {2, 1});
    NodeHandle P2 = np.add({2, 1, 0}, {2, 1});
    NodeHandle P3 = np.add({3, 1, 0}, {2, 1});
    NodeHandle B = np.add({4, 0, 0}, {2, 1});
    NodeHandle Q = np.add({2, -1, 0}, {2, 1});
    CellPool cp;
    auto tri = [&](NodeHandle a, NodeHandle b, NodeHandle c) {
        cp.add(ElementType::Tri3, std::vector<NodeHandle>{a, b, c}, {2, 5}, np);
    };
    tri(A, P1, Q); tri(P1, P2, Q); tri(P2, P3, Q); tri(P3, B, Q);
    EdgePool ep;
    for (const auto& ch : cp.activeBySequence()) ep.insertCellEdges(np, cp, ch);
    f.ws = std::make_unique<TopologyWorkspace>(std::move(np), std::move(cp),
                                               std::move(ep), FacePool{});
    auto V = [&](std::uint64_t id, NodeHandle h, double x, double y) {
        f.cs.vertices[id] = h; f.cs.params[h] = bnd::ParamPoint{{x, y}, 7};
    };
    V(100, A, 0, 0); V(101, B, 4, 0); V(200, P1, 1, 1); V(201, P2, 2, 1);
    V(202, P3, 3, 1); V(300, Q, 2, -1);
    f.seg.id = 500; f.seg.a = 100; f.seg.b = 101; f.seg.chart = 7; f.seg.cls = {1, 9};
    f.cs.segments.push_back(f.seg);
    f.cs.stamp = bnd::OwnershipStamp::of(*f.ws);
    return f;
}

// Existing-node reuse fan: interior node X (param only, not a constraint vertex)
// on the segment is reused -> chain of 5 (A, s1, X, s2, B), one atomic commit.
static RouteFix buildReuseRouteFix() {
    RouteFix f;
    NodePool np;
    NodeHandle A = np.add({0, 0, 0}, {2, 1});
    NodeHandle U1 = np.add({1, 1, 0}, {2, 1});
    NodeHandle L1 = np.add({1, -1, 0}, {2, 1});
    // A reused interior chain node must already carry the prescribed
    // dimension-1 classification; recovery never silently reclassifies it.
    NodeHandle X = np.add({2, 0, 0}, {1, 9});
    NodeHandle U2 = np.add({3, 1, 0}, {2, 1});
    NodeHandle L2 = np.add({3, -1, 0}, {2, 1});
    NodeHandle B = np.add({4, 0, 0}, {2, 1});
    CellPool cp;
    auto tri = [&](NodeHandle a, NodeHandle b, NodeHandle c) {
        cp.add(ElementType::Tri3, std::vector<NodeHandle>{a, b, c}, {2, 5}, np);
    };
    tri(A, U1, L1); tri(U1, X, L1); tri(U1, U2, X);
    tri(X, U2, L2); tri(X, L2, L1); tri(U2, B, L2);
    EdgePool ep;
    for (const auto& ch : cp.activeBySequence()) ep.insertCellEdges(np, cp, ch);
    f.ws = std::make_unique<TopologyWorkspace>(std::move(np), std::move(cp),
                                               std::move(ep), FacePool{});
    auto V = [&](std::uint64_t id, NodeHandle h, double x, double y) {
        f.cs.vertices[id] = h; f.cs.params[h] = bnd::ParamPoint{{x, y}, 7};
    };
    V(100, A, 0, 0); V(101, B, 4, 0); V(200, U1, 1, 1); V(201, L1, 1, -1);
    V(202, U2, 3, 1); V(203, L2, 3, -1);
    f.cs.params[X] = bnd::ParamPoint{{2, 0}, 7};   // reused: param only
    f.seg.id = 500; f.seg.a = 100; f.seg.b = 101; f.seg.chart = 7; f.seg.cls = {1, 9};
    f.cs.segments.push_back(f.seg);
    f.cs.stamp = bnd::OwnershipStamp::of(*f.ws);
    return f;
}

static long routeSweep(const char* label, RouteFix (*mk)(),
                       std::size_t chainLen, const bnd::SteinerPolicy& pol) {
    long k = 1, faults = 0, tier1 = 0, tier2 = 0;
    for (;; ++k) {
        RouteFix f = mk();
        f.ws->rebuildAdjacency();
        // Put the live parameter map immediately below its rehash threshold.
        // Inserting the staged-node entry into this LIVE map therefore changes
        // bucket state. Erasing the entry on failure cannot restore that state,
        // while a genuine off-to-the-side shadow leaves it untouched.
        f.cs.params.max_load_factor(f.cs.params.load_factor() + 0.01F);
        MockProv prov;
        const auto full0 = paacs_tests::observableDigest(*f.ws);
        const auto log0 = paacs_tests::observableLogicalDigest(*f.ws);
        const std::size_t nodes0 = f.ws->nodes().activeBySequence().size();
        const std::size_t params0 = f.cs.params.size();
        const std::size_t auth0 = f.cs.authParams.size();
        const std::size_t paramBuckets0 = f.cs.params.bucket_count();
        const std::size_t authBuckets0 = f.cs.authParams.bucket_count();
        gFailCountdown = k;
        // The public route is LEAK-PROOF: it must NEVER throw - it maps every
        // allocation failure to a structured outcome.
        bnd::SteinerResult res =
            bnd::recoverSteinerSegment(*f.ws, f.cs, f.seg, prov, prov.st, pol);
        gFailCountdown = 0;
        if (res.status != bnd::SteinerStatus::Recovered) {
            ++faults;
            // exact status/reason: a swept allocation fault is always Failed
            // with an allocation-class reason.
            REQUIRE(res.status == bnd::SteinerStatus::Failed);
            REQUIRE(res.fail == bnd::SteinerFail::CandidateAllocFailed ||
                    res.fail == bnd::SteinerFail::StagedNodeAllocFailed ||
                    res.fail == bnd::SteinerFail::CavityThrew);
            // TIER attribution: tier 1 -> FULL digest equality; tier 2 ->
            // LOGICAL digest equality (capacity residue permitted).
            if (isTier2(res.fail)) {
                ++tier2;
                REQUIRE(paacs_tests::observableLogicalDigest(*f.ws) == log0);
            } else {
                ++tier1;
                REQUIRE(paacs_tests::observableDigest(*f.ws) == full0);
            }
            // NO staged node became current, and NO metadata was published.
            REQUIRE(f.ws->nodes().activeBySequence().size() == nodes0);
            REQUIRE(f.cs.params.size() == params0);
            REQUIRE(f.cs.authParams.size() == auth0);
            // Metadata is BND state, not an ASCR pool: a failed operation may
            // not leave rehash/capacity residue from inserting into the LIVE
            // maps. A prepared shadow published by no-throw swap preserves
            // these values exactly.
            REQUIRE(f.cs.params.bucket_count() == paramBuckets0);
            REQUIRE(f.cs.authParams.bucket_count() == authBuckets0);
            REQUIRE(measureWorkspaceMemory(*f.ws).consistent());
            REQUIRE(f.ws->nodes().stats().bytesConsistent());
            REQUIRE(f.ws->cells().stats().bytesConsistent());
            REQUIRE(f.ws->edges().stats().bytesConsistent());
            // disarmed retry on the SAME workspace recovers AND publishes.
            auto r2 =
                bnd::recoverSteinerSegment(*f.ws, f.cs, f.seg, prov, prov.st, pol);
            if (r2.status != bnd::SteinerStatus::Recovered) {
                std::printf(
                    "retry[%s] failed at allocation %ld: status=%d reject=%d "
                    "refuse=%d fail=%d chain=%zu\n",
                    label, k, static_cast<int>(r2.status),
                    static_cast<int>(r2.reject), static_cast<int>(r2.refuse),
                    static_cast<int>(r2.fail), r2.recoveredChain.size());
            }
            REQUIRE(r2.status == bnd::SteinerStatus::Recovered);
            REQUIRE(r2.recoveredChain.size() == chainLen);
            if (r2.recoveredChain.size() > 1)
                REQUIRE(f.cs.paramOf(r2.recoveredChain[1]) != nullptr);
            REQUIRE(f.cs.stamp.matches(*f.ws));
            f.ws->rebuildAdjacency();
            REQUIRE(f.ws->adjacencyReady());
        } else {
            // completed without a fault at k: sweep exhausted the route
            REQUIRE(res.recoveredChain.size() == chainLen);
            REQUIRE(f.cs.paramOf(res.recoveredChain[1]) != nullptr);
            REQUIRE(f.cs.stamp.matches(*f.ws));
            f.ws->rebuildAdjacency();
            REQUIRE(f.ws->adjacencyReady());
            break;
        }
    }
    std::printf("route[%s]: %ld faults (tier1 %ld / tier2 %ld), 0 escaped\n",
                label, faults, tier1, tier2);
    return faults;
}

// The set-level public wrapper must be just as leak-proof as its per-entity
// helper. In particular its report reservations and segment-vector copy occur
// before the first entity, so faults there must become a COMPLETE structured
// report rather than escape the five-outcome interface.
static long setRouteSweep() {
    long k = 1, faults = 0, escapes = 0;
    for (;; ++k) {
        RouteFix f = buildRouteFix();
        f.ws->rebuildAdjacency();
        MockProv prov;
        const auto full0 = paacs_tests::observableDigest(*f.ws);
        const auto log0 = paacs_tests::observableLogicalDigest(*f.ws);
        bnd::SteinerSegmentsReport rep;
        bool escaped = false;
        gFailCountdown = k;
        try {
            rep = bnd::recoverSteinerSegments(
                *f.ws, f.cs, prov, prov.st);
        } catch (...) {
            escaped = true;
        }
        gFailCountdown = 0;
        if (escaped) {
            ++escapes;
            ++faults;
            REQUIRE(paacs_tests::observableDigest(*f.ws) == full0);
        } else {
            REQUIRE(rep.processed.size() + rep.unattempted.size() == 1);
            REQUIRE(rep.processed.size() == 1);
            if (!rep.processed.empty() &&
                rep.processed[0].result.status !=
                    bnd::SteinerStatus::Recovered) {
                ++faults;
                const auto reason = rep.processed[0].result.fail;
                REQUIRE(rep.processed[0].result.status ==
                        bnd::SteinerStatus::Failed);
                if (isTier2(reason))
                    REQUIRE(paacs_tests::observableLogicalDigest(*f.ws) ==
                            log0);
                else
                    REQUIRE(paacs_tests::observableDigest(*f.ws) == full0);
            } else {
                break;
            }
        }
        // Every failed/escaped injection must leave the same fixture reusable.
        const auto retry =
            bnd::recoverSteinerSegments(*f.ws, f.cs, prov, prov.st);
        REQUIRE(retry.processed.size() == 1);
        REQUIRE(retry.unattempted.empty());
        if (!retry.processed.empty())
            REQUIRE(retry.processed[0].result.status ==
                    bnd::SteinerStatus::Recovered);
    }
    REQUIRE(escapes == 0);
    std::printf("set-route: %ld faults / %ld escaped exceptions\n",
                faults, escapes);
    return faults;
}

// A genuine two-entity orchestration sweep. Entity 400 is already present and
// already marked, so it performs NO commit. Entity 500 is the transverse
// single-Steiner route. An allocation failure between them must never be
// attributed as 13.2a IndexRebuildFailed: there was no committed prefix to
// rebuild after. This distinguishes real post-commit handoff from an
// unconditional per-entity session reopen.
static long setRouteTwoEntitySweep() {
    long k = 1, faults = 0, escapes = 0, falseRebuilds = 0;
    for (;; ++k) {
        RouteFix f = buildRouteFix();
        f.ws->rebuildAdjacency();
        bnd::ConstraintSegment present;
        present.id = 400;
        present.a = 100;
        present.b = 101;
        present.chart = 7;
        present.cls = {1, 9};
        const auto eh =
            f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[1]);
        REQUIRE(eh.has_value());
        if (eh) {
            f.ws->testEdgesMutable().setRequiredConstraint(*eh, true);
            f.ws->testEdgesMutable().setClassification(*eh, present.cls);
        }
        f.cs.segments = {present, f.seg};
        f.cs.stamp = bnd::OwnershipStamp::of(*f.ws);
        MockProv prov;
        const auto full0 = paacs_tests::observableDigest(*f.ws);
        bnd::SteinerSegmentsReport rep;
        bool escaped = false;
        gFailCountdown = k;
        try {
            rep = bnd::recoverSteinerSegments(
                *f.ws, f.cs, prov, prov.st);
        } catch (...) {
            escaped = true;
        }
        gFailCountdown = 0;
        if (escaped) {
            ++escapes;
            ++faults;
            REQUIRE(paacs_tests::observableDigest(*f.ws) == full0);
        } else {
            REQUIRE(rep.processed.size() + rep.unattempted.size() == 2);
            if (rep.processed.size() >= 2 &&
                rep.processed[0].result.status ==
                    bnd::SteinerStatus::AlreadyPresent &&
                rep.processed[1].result.status ==
                    bnd::SteinerStatus::Failed &&
                rep.processed[1].result.fail ==
                    bnd::SteinerFail::IndexRebuildFailed) {
                ++falseRebuilds;
            }
            const bool recovered =
                rep.processed.size() == 2 &&
                rep.processed[0].result.status ==
                    bnd::SteinerStatus::AlreadyPresent &&
                rep.processed[1].result.status ==
                    bnd::SteinerStatus::Recovered;
            if (recovered) break;
            ++faults;
        }
        // A fresh retry is sufficient here: faults can occur during the second
        // entity after it commits, so the faulted workspace need not remain a
        // pre-operation fixture in every injection case.
    }
    REQUIRE(escapes == 0);
    REQUIRE(falseRebuilds == 0);
    std::printf(
        "set-route-2: %ld faults / %ld escapes / %ld false rebuilds\n",
        faults, escapes, falseRebuilds);
    return faults;
}
// Committed-prefix 13.2a sweep: TWO disjoint quads, entity 500 genuinely
// Recovered (real commit), entity 501 pending. When the post-commit session
// rebuild for entity 501 faults, entity 501 MUST be Failed(IndexRebuildFailed)
// with the committed prefix (500 Recovered) preserved and a complete report; a
// disarmed retry then recovers both.
static long setRouteCommittedPrefixSweep() {
    long k = 1, faults = 0, escapes = 0, idxAttrib = 0;
    for (;; ++k) {
        // two disjoint quads
        NodePool np;
        std::vector<NodeHandle> n;
        auto addn = [&](double x, double y) {
            n.push_back(np.add({x, y, 0}, {2, 1}));
        };
        addn(0, 0); addn(2, 0); addn(2, 2); addn(0, 2);
        addn(10, 0); addn(12, 0); addn(12, 2); addn(10, 2);
        CellPool cp;
        auto tri = [&](int a, int b, int c) {
            cp.add(ElementType::Tri3,
                   std::vector<NodeHandle>{n[std::size_t(a)], n[std::size_t(b)],
                                           n[std::size_t(c)]},
                   {2, 5}, np);
        };
        tri(0, 1, 2); tri(0, 2, 3); tri(4, 5, 6); tri(4, 6, 7);
        EdgePool ep;
        for (const auto& ch : cp.activeBySequence())
            ep.insertCellEdges(np, cp, ch);
        auto ws = std::make_unique<TopologyWorkspace>(
            std::move(np), std::move(cp), std::move(ep), FacePool{});
        ws->rebuildAdjacency();
        bnd::ConstraintSet cs;
        const Vec2 uv[8] = {{0, 0}, {2, 0}, {2, 2}, {0, 2},
                            {10, 0}, {12, 0}, {12, 2}, {10, 2}};
        for (int i = 0; i < 8; ++i) {
            cs.vertices[100 + std::uint64_t(i)] = n[std::size_t(i)];
            cs.params[n[std::size_t(i)]] = bnd::ParamPoint{uv[i], 7};
        }
        cs.segments = {{500, 101, 103, 7, {1, 9}}, {501, 105, 107, 7, {1, 9}}};
        cs.stamp = bnd::OwnershipStamp::of(*ws);
        MockProv prov;
        const auto L0 = paacs_tests::observableLogicalDigest(*ws);
        bnd::SteinerSegmentsReport rep;
        bool escaped = false;
        gFailCountdown = k;
        try {
            rep = bnd::recoverSteinerSegments(*ws, cs, prov, prov.st);
        } catch (...) {
            escaped = true;
        }
        gFailCountdown = 0;
        if (escaped) { ++escapes; ++faults; continue; }
        // completeness + disjointness ALWAYS
        REQUIRE(rep.processed.size() + rep.unattempted.size() == 2);
        {
            std::map<bnd::ConstraintEntityId, int> seen;
            for (const auto& e : rep.processed) ++seen[e.entity];
            for (const auto id : rep.unattempted) ++seen[id];
            REQUIRE(seen.size() == 2);
            for (const auto& [id, c2] : seen) { (void)id; REQUIRE(c2 == 1); }
        }
        const bool bothRecovered =
            rep.processed.size() == 2 &&
            rep.processed[0].result.status == bnd::SteinerStatus::Recovered &&
            rep.processed[1].result.status == bnd::SteinerStatus::Recovered;
        // the committed-prefix index-rebuild attribution: 500 Recovered, 501
        // Failed(IndexRebuildFailed).
        if (rep.processed.size() == 2 &&
            rep.processed[0].result.status == bnd::SteinerStatus::Recovered &&
            rep.processed[1].result.status == bnd::SteinerStatus::Failed &&
            rep.processed[1].result.fail ==
                bnd::SteinerFail::IndexRebuildFailed) {
            ++idxAttrib;
            REQUIRE(rep.processed[0].entity == 500);
            REQUIRE(rep.processed[1].entity == 501);
            // committed prefix preserved: entity 500's new node is live.
            REQUIRE(ws->nodes().activeBySequence().size() == 9);
            // disarmed retry recovers the pending entity.
            auto retry = bnd::recoverSteinerSegments(*ws, cs, prov, prov.st);
            REQUIRE(!retry.processed.empty());
            REQUIRE(retry.processed.back().result.status ==
                    bnd::SteinerStatus::Recovered);
        }
        if (bothRecovered) break;
        (void)L0;
        ++faults;
    }
    REQUIRE(escapes == 0);
    REQUIRE(idxAttrib >= 1);   // the attribution path was genuinely reached
    std::printf(
        "set-route-prefix: %ld faults / %ld escapes / %ld index-attrib\n",
        faults, escapes, idxAttrib);
    return faults;
}

enum class PrefixFlavor { Multi, Reuse };
struct PrefixSetFix {
    std::unique_ptr<TopologyWorkspace> ws;
    bnd::ConstraintSet cs;
    std::size_t nodesAfterPrefix = 0;
    std::size_t prefixChainSize = 0;
    std::size_t prefixCreated = 0;
};

// Two genuine set fixtures. Entity 500 is either a K=3 multi-Steiner strip or
// a reuse strip (A,s1,X,s2,B); entity 501 is a disjoint single-Steiner quad.
// Thus every post-500 session rebuild is preceded by a real multi/reuse commit.
static PrefixSetFix buildPrefixSetFix(PrefixFlavor flavor) {
    PrefixSetFix f;
    NodePool np;
    CellPool cp;
    bnd::ConstraintSet cs;
    auto tri = [&](NodeHandle a, NodeHandle b, NodeHandle c) {
        cp.add(ElementType::Tri3, std::vector<NodeHandle>{a, b, c},
               {2, 5}, np);
    };
    auto vertex = [&](std::uint64_t id, NodeHandle h, double x, double y) {
        cs.vertices[id] = h;
        cs.params[h] = bnd::ParamPoint{{x, y}, 7};
    };

    if (flavor == PrefixFlavor::Multi) {
        const NodeHandle A = np.add({0, 0, 0}, {2, 1});
        const NodeHandle P1 = np.add({1, 1, 0}, {2, 1});
        const NodeHandle P2 = np.add({2, 1, 0}, {2, 1});
        const NodeHandle P3 = np.add({3, 1, 0}, {2, 1});
        const NodeHandle B = np.add({4, 0, 0}, {2, 1});
        const NodeHandle Q = np.add({2, -1, 0}, {2, 1});
        tri(A, P1, Q);
        tri(P1, P2, Q);
        tri(P2, P3, Q);
        tri(P3, B, Q);
        vertex(100, A, 0, 0);
        vertex(101, B, 4, 0);
        vertex(200, P1, 1, 1);
        vertex(201, P2, 2, 1);
        vertex(202, P3, 3, 1);
        vertex(203, Q, 2, -1);
        f.nodesAfterPrefix = 13;  // 10 initial + three staged
        f.prefixChainSize = 5;
        f.prefixCreated = 3;
    } else {
        const NodeHandle A = np.add({0, 0, 0}, {2, 1});
        const NodeHandle U1 = np.add({1, 1, 0}, {2, 1});
        const NodeHandle L1 = np.add({1, -1, 0}, {2, 1});
        const NodeHandle X = np.add({2, 0, 0}, {1, 9});
        const NodeHandle U2 = np.add({3, 1, 0}, {2, 1});
        const NodeHandle L2 = np.add({3, -1, 0}, {2, 1});
        const NodeHandle B = np.add({4, 0, 0}, {2, 1});
        tri(A, U1, L1);
        tri(U1, X, L1);
        tri(U1, U2, X);
        tri(X, U2, L2);
        tri(X, L2, L1);
        tri(U2, B, L2);
        vertex(100, A, 0, 0);
        vertex(101, B, 4, 0);
        vertex(200, U1, 1, 1);
        vertex(201, L1, 1, -1);
        vertex(202, U2, 3, 1);
        vertex(203, L2, 3, -1);
        cs.params[X] = bnd::ParamPoint{{2, 0}, 7};
        f.nodesAfterPrefix = 13;  // 11 initial + two staged; X is reused
        f.prefixChainSize = 5;
        f.prefixCreated = 2;
    }

    // Disjoint entity 501: a square triangulated on 110-112, while the
    // prescribed 111-113 diagonal requires one Steiner point.
    const NodeHandle C0 = np.add({10, 0, 0}, {2, 1});
    const NodeHandle C1 = np.add({12, 0, 0}, {2, 1});
    const NodeHandle C2 = np.add({12, 2, 0}, {2, 1});
    const NodeHandle C3 = np.add({10, 2, 0}, {2, 1});
    // Match the clockwise source-embedding orientation of the prefix fixture.
    tri(C0, C2, C1);
    tri(C0, C3, C2);
    vertex(110, C0, 10, 0);
    vertex(111, C1, 12, 0);
    vertex(112, C2, 12, 2);
    vertex(113, C3, 10, 2);

    EdgePool ep;
    for (const auto& ch : cp.activeBySequence())
        ep.insertCellEdges(np, cp, ch);
    f.ws = std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), FacePool{});
    f.ws->rebuildAdjacency();
    cs.segments = {{500, 100, 101, 7, {1, 9}},
                   {501, 111, 113, 7, {1, 9}}};
    cs.stamp = bnd::OwnershipStamp::of(*f.ws);
    f.cs = std::move(cs);
    return f;
}

static void requireCompleteDisjoint(const bnd::SteinerSegmentsReport& rep,
                                    std::size_t count) {
    REQUIRE(rep.processed.size() + rep.unattempted.size() == count);
    std::map<bnd::ConstraintEntityId, int> seen;
    for (const auto& e : rep.processed) ++seen[e.entity];
    for (const auto id : rep.unattempted) ++seen[id];
    REQUIRE(seen.size() == count);
    for (const auto& [id, multiplicity] : seen) {
        (void)id;
        REQUIRE(multiplicity == 1);
    }
}

static long setRouteGeneralPrefixSweep(const char* label,
                                       PrefixFlavor flavor) {
    long k = 1, faults = 0, escapes = 0, indexAttributions = 0;
    bnd::SteinerPolicy policy;
    policy.splitBudget = 8;
    policy.minFeatureLevel = 52;
    {
        PrefixSetFix clean = buildPrefixSetFix(flavor);
        MockProv provider;
        const auto baseline = bnd::recoverSteinerSegments(
            *clean.ws, clean.cs, provider, provider.st, policy);
        requireCompleteDisjoint(baseline, 2);
        const bool cleanSuccess =
            baseline.processed.size() == 2 &&
            baseline.processed[0].result.status ==
                bnd::SteinerStatus::Recovered &&
            baseline.processed[1].result.status ==
                bnd::SteinerStatus::Recovered;
        REQUIRE(cleanSuccess);
        if (!cleanSuccess) {
            if (!baseline.processed.empty())
                std::printf(
                    "set-route-prefix[%s] clean failure: processed=%zu "
                    "first status=%d reject=%d refuse=%d fail=%d\n",
                    label, baseline.processed.size(),
                    static_cast<int>(baseline.processed[0].result.status),
                    static_cast<int>(baseline.processed[0].result.reject),
                    static_cast<int>(baseline.processed[0].result.refuse),
                    static_cast<int>(baseline.processed[0].result.fail));
            return 0;
        }
    }
    for (;; ++k) {
        PrefixSetFix f = buildPrefixSetFix(flavor);
        MockProv prov;
        bnd::SteinerSegmentsReport rep;
        bool escaped = false;
        gFailCountdown = k;
        try {
            rep = bnd::recoverSteinerSegments(
                *f.ws, f.cs, prov, prov.st, policy);
        } catch (...) {
            escaped = true;
        }
        gFailCountdown = 0;
        if (escaped) {
            ++escapes;
            ++faults;
            continue;
        }

        requireCompleteDisjoint(rep, 2);
        const bool bothRecovered =
            rep.processed.size() == 2 &&
            rep.processed[0].entity == 500 &&
            rep.processed[0].result.status ==
                bnd::SteinerStatus::Recovered &&
            rep.processed[1].entity == 501 &&
            rep.processed[1].result.status ==
                bnd::SteinerStatus::Recovered;

        if (rep.processed.size() == 2 &&
            rep.processed[0].entity == 500 &&
            rep.processed[0].result.status ==
                bnd::SteinerStatus::Recovered &&
            rep.processed[1].entity == 501 &&
            rep.processed[1].result.status == bnd::SteinerStatus::Failed &&
            rep.processed[1].result.fail ==
                bnd::SteinerFail::IndexRebuildFailed) {
            ++indexAttributions;
            const auto prefixChain = rep.processed[0].result.recoveredChain;
            REQUIRE(prefixChain.size() == f.prefixChainSize);
            REQUIRE(rep.processed[0].result.createdParams.size() ==
                    f.prefixCreated);
            REQUIRE(f.ws->nodes().activeBySequence().size() ==
                    f.nodesAfterPrefix);
            for (const auto& h : prefixChain)
                REQUIRE(f.ws->nodes().isCurrent(h));
            for (std::size_t i = 0; i + 1 < prefixChain.size(); ++i) {
                const auto edge = f.ws->edges().find(
                    f.ws->nodes(), prefixChain[i], prefixChain[i + 1]);
                REQUIRE(edge.has_value());
                if (edge) {
                    const auto ed = f.ws->edges().get(*edge);
                    REQUIRE(ed.requiredConstraint);
                    REQUIRE(ed.cls.dim == 1 && ed.cls.tag == 9);
                }
            }

            // Disarmed SAME-workspace retry: the committed prefix routes as
            // AlreadyPresent with the identical chain, then 501 recovers.
            const auto retry = bnd::recoverSteinerSegments(
                *f.ws, f.cs, prov, prov.st, policy);
            requireCompleteDisjoint(retry, 2);
            REQUIRE(retry.processed.size() == 2);
            if (retry.processed.size() == 2) {
                REQUIRE(retry.processed[0].entity == 500);
                REQUIRE(retry.processed[0].result.status ==
                        bnd::SteinerStatus::AlreadyPresent);
                REQUIRE(retry.processed[0].result.recoveredChain ==
                        prefixChain);
                REQUIRE(retry.processed[1].entity == 501);
                REQUIRE(retry.processed[1].result.status ==
                        bnd::SteinerStatus::Recovered);
            }
        }
        if (bothRecovered) break;
        ++faults;
    }
    REQUIRE(escapes == 0);
    REQUIRE(indexAttributions >= 1);
    std::printf(
        "set-route-prefix[%s]: %ld faults / %ld escapes / %ld "
        "index-attrib\n",
        label, faults, escapes, indexAttributions);
    return faults;
}

// Provider-call failure attribution: pointAt bad_alloc -> CandidateAllocFailed;
// other throw -> GeometryProviderError; domain -> ProviderDomainRejected; stamp
// throws are contained identically. Full digest unchanged in every case.
static void routeProviderFailures() {
    struct P : bnd::GeometryProvider {
        bnd::GeometryProviderStamp st{42, 1};
        int mode = 0;   // 1 bad_alloc, 2 other, 3 domain, 4 stamp bad_alloc,
                        // 5 stamp other
        bnd::GeometryProviderStamp stamp() const override {
            if (mode == 4) throw std::bad_alloc();
            if (mode == 5) throw std::runtime_error("stamp");
            return st;
        }
        Vec3 pointAt(bnd::ChartId, const Vec2&) const override {
            if (mode == 1) throw std::bad_alloc();
            if (mode == 2) throw std::runtime_error("point");
            if (mode == 3) throw bnd::GeometryProviderDomainError("dom");
            return {1, 1, 0};
        }
    };
    struct Case { int mode; bnd::SteinerStatus st; bnd::SteinerFail f;
                  bnd::SteinerReject rj; };
    const Case cases[] = {
        {1, bnd::SteinerStatus::Failed, bnd::SteinerFail::CandidateAllocFailed,
         bnd::SteinerReject::None},
        {2, bnd::SteinerStatus::Failed, bnd::SteinerFail::GeometryProviderError,
         bnd::SteinerReject::None},
        {3, bnd::SteinerStatus::Rejected, bnd::SteinerFail::None,
         bnd::SteinerReject::ProviderDomainRejected},
        {4, bnd::SteinerStatus::Failed, bnd::SteinerFail::CandidateAllocFailed,
         bnd::SteinerReject::None},
        {5, bnd::SteinerStatus::Failed, bnd::SteinerFail::GeometryProviderError,
         bnd::SteinerReject::None},
    };
    for (const auto& cse : cases) {
        RouteFix f = buildRouteFix();
        f.ws->rebuildAdjacency();
        P p;
        p.mode = cse.mode;
        const auto full0 = paacs_tests::observableDigest(*f.ws);
        const std::size_t params0 = f.cs.params.size();
        bnd::SteinerResult res;
        bool escaped = false;
        try {
            res = bnd::recoverSteinerSegment(*f.ws, f.cs, f.seg, p, p.st);
        } catch (...) {
            escaped = true;
        }
        REQUIRE(!escaped);
        REQUIRE(res.status == cse.st);
        REQUIRE(res.fail == cse.f);
        REQUIRE(res.reject == cse.rj);
        REQUIRE(paacs_tests::observableDigest(*f.ws) == full0);
        REQUIRE(f.cs.params.size() == params0);   // no metadata published
    }
    std::printf("route provider-failure taxonomy: 5/5, 0 escaped\n");
}

int main(int argc, char**) {
    // Keep qualification progress visible even if an injected allocation
    // exposes a hard process failure on a platform runtime.
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    // Focused developer mode avoids replaying the older sweeps while
    // diagnosing the two new committed-prefix fixtures. CTest invokes the
    // no-argument path below and therefore still runs the complete battery.
    if (argc > 1) {
        const long multi =
            setRouteGeneralPrefixSweep("multi", PrefixFlavor::Multi);
        const long reuse =
            setRouteGeneralPrefixSweep("reuse", PrefixFlavor::Reuse);
        REQUIRE(multi >= 10);
        REQUIRE(reuse >= 10);
        if (gBad)
            std::printf("BND2 PREFIX PROBE FAILED (%d)\n", gBad);
        else
            std::printf("BND2 PREFIX PROBE PASSED\n");
        return gBad ? 1 : 0;
    }
    long fa = 0, fb = 0, ga = 0, gb = 0;
    runSweeps("single-fan", fanSpec, 4, 1, fa, fb);
    runSweeps("multi-steiner", multiSpec, 6, 2, ga, gb);
    bnd::SteinerPolicy singlePol;                        // default: budget 1
    bnd::SteinerPolicy multiPol; multiPol.splitBudget = 8;
    const long rfaults = routeSweep("single", buildRouteFix, 3, singlePol);
    const long rmfaults = routeSweep("multi", buildMultiRouteFix, 5, multiPol);
    const long rrfaults = routeSweep("reuse", buildReuseRouteFix, 5, multiPol);
    const long sfaults = setRouteSweep();
    const long s2faults = setRouteTwoEntitySweep();
    const long spfaults = setRouteCommittedPrefixSweep();
    const long smpfaults =
        setRouteGeneralPrefixSweep("multi", PrefixFlavor::Multi);
    const long srpfaults =
        setRouteGeneralPrefixSweep("reuse", PrefixFlavor::Reuse);
    routeProviderFailures();
    REQUIRE(fa >= 10);
    REQUIRE(fb >= 1);
    REQUIRE(ga >= 10);
    REQUIRE(gb >= 1);
    REQUIRE(rfaults >= 10);
    REQUIRE(rmfaults >= 10);   // multi-Steiner route allocation sweep
    REQUIRE(rrfaults >= 10);   // existing-node reuse route allocation sweep
    REQUIRE(sfaults >= 10);
    REQUIRE(s2faults >= 10);
    REQUIRE(spfaults >= 10);
    REQUIRE(smpfaults >= 10);
    REQUIRE(srpfaults >= 10);
    if (gBad) {
        std::printf("BND2 PROBE FAILED (%d)\n", gBad);
        return 1;
    }
    std::printf("BND2 PROBE PASSED\n");
    return 0;
}
