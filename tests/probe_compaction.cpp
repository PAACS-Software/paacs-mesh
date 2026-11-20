// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Standalone COMPACTION allocation-failure sweep (NOT part of the test
// binary).
//
// Build (canonical target: `cmake --build . --target paacs_probe_compaction`;
// standalone equivalent, WITH the standard warning flags and no line
// continuations - a trailing backslash inside a // comment triggers
// -Wcomment):
//   g++ -std=c++20 -Wall -Wextra -Wpedantic -DPAACS_TOPO_TESTING -g -O1
//       -Iinclude src/topo/ElementTraits.cpp tests/probe_compaction.cpp
//       -o probe_compaction
//   ./probe_compaction
//
// C2 review blocker: the named failure proof only injected throws at
// seven stage BOUNDARIES. This probe sweeps EVERY allocation that
// compactWorkspace performs — map insertion, partially built edge/face
// maps, incidence reconstruction, workspace/adjacency construction and
// V2 hash temporaries alike — by failing the k-th global operator new
// for every k in 1..N, verifying full source-state preservation after
// each bad_alloc, and finishing with a clean successful compaction.
//
// (Run WITHOUT sanitizers — ASan replaces global operator new.)
#include <cstdio>
#include <cstdlib>
#include <new>
#include <vector>

#include "paacs/topo/Compaction.hpp"
#include "workspace_digest.hpp"

// --- global allocation countdown ------------------------------------------
static bool gArmed = false;
static long gAllocs = 0;
static long gFailAt = 0;     // 0 = count only

// GCC inlines our operator new (which obtains storage from malloc) and
// then sees operator delete calling free on a pointer it has attributed
// to "operator new", so it reports -Wmismatched-new-delete. The pairing
// is in fact correct: these replacements ARE the global allocator, and
// malloc/free are each other's counterpart. The diagnostic is a GCC
// artifact of replacing global operator new, unavoidable for an
// allocation-counting probe, so it is suppressed NARROWLY around these
// definitions only - never for library or test code.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t n) {
    if (gArmed) {
        ++gAllocs;
        if (gFailAt > 0 && gAllocs == gFailAt) throw std::bad_alloc();
    }
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

// FULL observable-state digest — the SHARED implementation from
// tests/workspace_digest.hpp, identical to the one the test binary
// uses (C2 review: the earlier ad-hoc digest here omitted the five pool
// bindings and the V1 hash, so "full source digest" overclaimed).
// Computed with the counter DISARMED, since it allocates.
static std::vector<std::uint64_t> digest(const TopologyWorkspace& ws) {
    return paacs_tests::observableDigest(ws);
}

// Churned source WITH a retired slot, a freed slot and non-default
// metadata: the widest fixture, so the sweep covers the most code.
static std::unique_ptr<TopologyWorkspace> buildSource() {
    NodePool np;
    const NodeHandle doomed = np.add({7, 7, 7}, {0, 8});
    const NodeHandle freed = np.add({8, 8, 8}, {0, 9});
    std::vector<NodeHandle> n;
    for (int i = 0; i < 5; ++i)
        n.push_back(np.add({double(i), double(i % 2), double(i % 3)},
                           {3, 1}));
    np.testRetire(doomed);
    np.remove(freed);
    CellPool cp;
    auto t1 = cp.add(ElementType::Tet4,
                     std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                     {3, 5}, np);
    auto t2 = cp.add(ElementType::Tet4,
                     std::vector<NodeHandle>{n[0], n[2], n[1], n[4]},
                     {3, 6}, np);
    EdgePool ep;
    ep.insertCellEdges(np, cp, t1);
    ep.insertCellEdges(np, cp, t2);
    FacePool fp;
    fp.insertCellFaces(np, cp, ep, t1);
    fp.insertCellFaces(np, cp, ep, t2);
    np.setClassification(n[0], {2, 21});
    auto e01 = *ep.find(np, n[0], n[1]);
    ep.setClassification(e01, {1, 31});
    ep.setRequiredConstraint(e01, true);
    auto shared = *fp.find(np, std::vector<NodeHandle>{n[0], n[1], n[2]});
    fp.setClassification(shared, {2, 41});
    return std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), std::move(fp));
}

int main() {
    auto src = buildSource();
    const auto d0 = digest(*src);

    // 1. Count every allocation a successful compaction makes.
    gAllocs = 0;
    gFailAt = 0;
    gArmed = true;
    { (void)compactWorkspace(*src); }
    gArmed = false;
    const long total = gAllocs;
    std::printf("compaction allocations to sweep: %ld\n", total);
    REQUIRE(total > 0);
    REQUIRE(digest(*src) == d0);

    // 2. Fail the k-th allocation for EVERY k, and prove the source is
    //    observably unchanged after each failure. RAII destruction of
    //    the temporaries is the entire exception guarantee.
    long threw = 0, completed = 0;
    for (long k = 1; k <= total; ++k) {
        gAllocs = 0;
        gFailAt = k;
        gArmed = true;
        bool caught = false;
        try {
            (void)compactWorkspace(*src);
            ++completed;
        } catch (const std::bad_alloc&) {
            caught = true;
        } catch (...) {
            gArmed = false;
            std::printf("FAIL: k=%ld threw a non-bad_alloc exception\n", k);
            ++gBad;
            continue;
        }
        gArmed = false;
        if (caught) ++threw;
        if (digest(*src) != d0) {
            std::printf("FAIL: source mutated after failure at k=%ld\n", k);
            ++gBad;
            break;
        }
    }
    std::printf("swept %ld allocations: %ld threw bad_alloc, %ld completed\n",
                total, threw, completed);
    REQUIRE(threw + completed == total);
    REQUIRE(threw > 0);
    // Every injection point must actually be reachable as a failure:
    // an allocation that never throws would mean the counter missed it.
    REQUIRE(completed == 0);
    REQUIRE(digest(*src) == d0);

    // 3. A clean compaction still succeeds afterwards.
    gArmed = false;
    gFailAt = 0;
    auto res = compactWorkspace(*src);
    REQUIRE(res.workspace != nullptr);
    REQUIRE(res.workspace->semanticTopologyHash() ==
            src->semanticTopologyHash());
    REQUIRE(res.nodeMap.size() == src->nodes().stats().live);
    REQUIRE(res.cellMap.size() == src->cells().stats().live);
    REQUIRE(res.edgeMap.size() == src->edges().stats().live);
    REQUIRE(res.faceMap.size() == src->faces().stats().live);
    REQUIRE(digest(*src) == d0);

    std::printf(gBad ? "%d FAILURE(S)\n" : "COMPACTION PROBE PASSED\n",
                gBad);
    return gBad ? 1 : 0;
}
