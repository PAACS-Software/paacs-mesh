// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// M4-B.8 Stage C probes (standalone; NOT built under sanitizers):
//   1. no_allocation_in_commit_prepared - global operator new counts
//      every allocation while detail::phase3Active() is set; the cavity
//      commit's Phase-3 tail must count ZERO.
//   2. allocation_sweep - every allocation of the WHOLE cavity
//      operation (begin/stage/plan/prepare/commit throwing zone) fails
//      in turn on a FRESH workspace per injected index k; after each
//      failure the logical state is restored, byte invariants hold, and
//      a disarmed retry on the SAME workspace commits successfully.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>

#include "paacs/topo/Cavity.hpp"
#include "workspace_digest.hpp"

static long gPhase3Allocs = 0;
static long gFailCountdown = 0;   // probe-armed bad_alloc injection
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
};
// WIDEST fixture (Stage-C review): the EXHAUSTED 2x1 strip - four
// cells, TWO created + TWO retired edges, a retained interior edge,
// one generation-exhausted slot (fresh materialization + fresh
// incidence entry) mixed with in-commit reuse, and a requirement on a
// created edge.
static Fix buildStrip() {
    Fix f;
    NodePool np;
    f.n.push_back(np.add({0, 0, 0}, {2, 1}));
    f.n.push_back(np.add({1, 0, 0}, {2, 1}));
    f.n.push_back(np.add({2, 0, 0}, {2, 1}));
    f.n.push_back(np.add({0, 1, 0}, {2, 1}));
    f.n.push_back(np.add({1, 1, 0}, {2, 1}));
    f.n.push_back(np.add({2, 1, 0}, {2, 1}));
    CellPool cp;
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[1], f.n[4]}, {2, 5}, np));
    f.c[0] = paacs::topo::detail::TransactionAccessTestProbe::
        exhaustLiveCell(cp, f.c[0]);
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[4], f.n[3]}, {2, 5}, np));
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[1], f.n[2], f.n[5]}, {2, 5}, np));
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[1], f.n[5], f.n[4]}, {2, 5}, np));
    EdgePool ep;
    for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
    ep.setRequiredConstraint(*ep.find(np, f.n[0], f.n[1]), true);
    ep.setClassification(*ep.find(np, f.n[0], f.n[1]), {1, 9});
    f.ws = std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), FacePool{});
    return f;
}
static TriCavitySpec stripSpec(const Fix& f) {
    TriCavitySpec s;
    s.removed = f.c;
    s.added = {TriCavityCell{{f.n[0], f.n[1], f.n[3]}, {2, 5}},
               TriCavityCell{{f.n[1], f.n[4], f.n[3]}, {2, 5}},
               TriCavityCell{{f.n[1], f.n[2], f.n[4]}, {2, 5}},
               TriCavityCell{{f.n[2], f.n[5], f.n[4]}, {2, 5}}};
    s.requirements = {CavityEdgeRequirement{f.n[1], f.n[3], {1, 7}}};
    return s;
}
int main() {
    // ---- 1. no_allocation_in_commit_prepared -------------------------
    {
        Fix f = buildStrip();                 // WIDEST fixture
        f.ws->rebuildAdjacency();
        const auto rev0 = f.ws->topologyRevision();
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(stripSpec(f));
        tx.prepare();
        gPhase3Allocs = 0;
        const auto out = tx.commit();
        REQUIRE(gPhase3Allocs == 0);
        REQUIRE(out.resultCells.size() == 4);
        REQUIRE(out.requiredEdges.size() == 1);
        REQUIRE(f.ws->topologyRevision() == rev0 + 1);
        f.ws->rebuildAdjacency();
        std::printf("phase3 allocations during commitCavityPrepared: %ld\n",
                    gPhase3Allocs);
    }
    // ---- 2. allocation_sweep (fresh workspace per k) ------------------
    long k = 1;
    long faults = 0;
    for (;; ++k) {
        Fix f = buildStrip();                 // widest fixture
        f.ws->rebuildAdjacency();
        const auto L0 = paacs_tests::observableLogicalDigest(*f.ws);
        bool faulted = false;
        {
            gFailCountdown = k;   // arm ONLY around the operation
            try {
                auto tx = f.ws->beginCavityTransaction();
                (void)tx.stage(stripSpec(f));
                (void)tx.computePlan();
                tx.prepare();
                gFailCountdown = 0;          // completed without fault:
                tx.abandon();                // sweep of the PREPARE side
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
            REQUIRE(f.ws->faces().stats().bytesConsistent());
        }
        // disarmed retry on the SAME workspace: full commit succeeds
        {
            auto tx = f.ws->beginCavityTransaction();
            (void)tx.stage(stripSpec(f));
            const auto out = tx.commit();
            REQUIRE(out.resultCells.size() == 4);
            REQUIRE(f.ws->state() ==
                    WorkspaceState::CommittedAdjacencyStale);
            f.ws->rebuildAdjacency();
            REQUIRE(f.ws->adjacencyReady());
        }
        if (faulted)
            ++faults;
        else
            break;                            // k exceeded prepare-side
    }
    // ---- 2b. commit-side sweep: fault inside commit()'s throwing zone
    long k2 = 1;
    long faults2 = 0;
    for (;; ++k2) {
        Fix f = buildStrip();                 // widest fixture
        f.ws->rebuildAdjacency();
        const auto L0 = paacs_tests::observableLogicalDigest(*f.ws);
        bool faulted = false;
        {
            auto tx = f.ws->beginCavityTransaction();
            (void)tx.stage(stripSpec(f));
            tx.prepare();                     // disarmed preparation
            gFailCountdown = k2;              // arm the COMMIT zone only
            try {
                const auto out = tx.commit();
                gFailCountdown = 0;
                REQUIRE(out.resultCells.size() == 4);
            } catch (const std::bad_alloc&) {
                faulted = true;
                gFailCountdown = 0;
                tx.abandon();                 // placeholders cleaned
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
            REQUIRE(f.ws->faces().stats().bytesConsistent());
            {   // clean retry to full commit on the SAME workspace
                auto tx = f.ws->beginCavityTransaction();
                (void)tx.stage(stripSpec(f));
                (void)tx.commit();
                f.ws->rebuildAdjacency();
            }
        } else {
            break;                            // commit completed: Phase 3
        }                                     // itself never allocates
    }
    std::printf("sweep: %ld prepare-side faults, %ld commit-zone faults\n",
                faults, faults2);
    REQUIRE(faults >= 10);
    REQUIRE(faults2 >= 1);
    if (gBad) {
        std::printf("CAVITY PROBE FAILED (%d)\n", gBad);
        return 1;
    }
    std::printf("CAVITY PROBE PASSED\n");
    return 0;
}
