// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Standalone Phase-3 allocation probe (NOT part of the test binary).
//
// Build (canonical target: `cmake --build . --target paacs_probe_transaction`;
// standalone equivalent, WITH the standard warning flags and no line
// continuations - a trailing backslash inside a // comment triggers
// -Wcomment):
//   g++ -std=c++20 -Wall -Wextra -Wpedantic -DPAACS_TOPO_TESTING -g -O1
//       -Iinclude src/topo/ElementTraits.cpp tests/probe_transaction.cpp
//       -o probe_transaction
//   ./probe_transaction
//
// Proves [topo.transaction.no_allocation_during_add_commit] and
// [topo.transaction.no_allocation_during_remove_commit] with a global
// operator-new counter armed by detail::phase3Active() (a testing-only
// flag set around Phase 3 by TopologyTransaction::commit). Coverage per
// the reviewer's list: fresh cell/edge/face slots; reused slots in all
// three pools; mixed reused+fresh; and a REPLACE commit (stage 4: net
// deltas with retained + created + retiring entities). (Run WITHOUT
// sanitizers - ASan
// replaces global operator new.)
#include <cstdio>
#include <cstdlib>
#include <new>

#include "paacs/topo/Workspace.hpp"

static long gPhase3Allocs = 0;
// See probe_compaction.cpp: -Wmismatched-new-delete is a GCC artifact of
// replacing global operator new with malloc/free. Suppressed narrowly
// around these definitions only.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t n) {
    if (paacs::topo::detail::phase3Active()) ++gPhase3Allocs;
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
#define REQUIRE(c) do { if (!(c)) { std::printf("FAIL %d: %s\n", __LINE__, #c); ++gBad; } } while (0)

int main() {
    NodePool np;
    std::vector<NodeHandle> n;
    for (int i = 0; i < 5; ++i)
        n.push_back(np.add({double(i), double(i % 2), double(i % 3)}, {3,1}));
    CellPool cp;
    auto c1 = cp.add(ElementType::Tet4,
                     std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                     {3, 1}, np);
    EdgePool ep; ep.insertCellEdges(np, cp, c1);
    FacePool fp; fp.insertCellFaces(np, cp, ep, c1);
    TopologyWorkspace ws(std::move(np), std::move(cp), std::move(ep),
                         std::move(fp));

    // 1. ADD with FRESH cell/edge/face slots (never freed before).
    {
        gPhase3Allocs = 0;
        auto tx = ws.beginTransaction();
        (void)tx.addCell(ElementType::Tet4,
                         std::vector<NodeHandle>{n[0], n[2], n[1], n[4]},
                         {3, 2});
        auto res = tx.commit();
        std::printf("fresh add: phase3 allocs = %ld\n", gPhase3Allocs);
        REQUIRE(gPhase3Allocs == 0);
        ws.rebuildAdjacency();

        // 2. REMOVE commit: zero allocations.
        TopologyPolicy pol; pol.permitOrphanNodes = true;
        gPhase3Allocs = 0;
        auto tx2 = ws.beginTransaction(pol);
        tx2.removeCell(res.cell);
        (void)tx2.commit();
        std::printf("remove: phase3 allocs = %ld\n", gPhase3Allocs);
        REQUIRE(gPhase3Allocs == 0);
        ws.rebuildAdjacency();
    }
    // 3. ADD with REUSED slots in all three pools (the remove above freed
    //    one cell slot, three edge slots, three face slots) plus fresh
    //    ones needed beyond the free lists = MIXED reused+fresh.
    {
        gPhase3Allocs = 0;
        auto tx = ws.beginTransaction();
        (void)tx.addCell(ElementType::Tet4,
                         std::vector<NodeHandle>{n[0], n[2], n[1], n[4]},
                         {3, 3});
        auto res = tx.commit();
        std::printf("mixed reuse add: phase3 allocs = %ld\n", gPhase3Allocs);
        REQUIRE(gPhase3Allocs == 0);
        REQUIRE(res.topologyChanged);
        ws.rebuildAdjacency();
        REQUIRE(ws.adjacency().cellCount() == 2);
        REQUIRE(ws.adjacency().boundaryFaces().size() == 6);
    }
    // 4. REPLACE commit (apex swap: 3 retained edges + 3 creates, 1
    //    retained face + 3 creates, 3 retirements each): zero allocations.
    {
        NodePool np2; std::vector<NodeHandle> m;
        for (int i = 0; i < 5; ++i)
            m.push_back(np2.add({double(i), 1, double(i % 2)}, {3, 1}));
        CellPool cp2;
        auto c = cp2.add(ElementType::Tet4,
                         std::vector<NodeHandle>{m[0], m[1], m[2], m[3]},
                         {3, 1}, np2);
        EdgePool ep2; ep2.insertCellEdges(np2, cp2, c);
        FacePool fp2; fp2.insertCellFaces(np2, cp2, ep2, c);
        TopologyWorkspace ws2(std::move(np2), std::move(cp2),
                              std::move(ep2), std::move(fp2));
        TopologyPolicy pol; pol.permitOrphanNodes = true;
        gPhase3Allocs = 0;
        auto tx = ws2.beginTransaction(pol);
        tx.replaceCellConnectivity(
            c, std::vector<NodeHandle>{m[0], m[1], m[2], m[4]});
        auto res = tx.commit();
        std::printf("replace: phase3 allocs = %ld\n", gPhase3Allocs);
        REQUIRE(gPhase3Allocs == 0);
        REQUIRE(res.cell == c);              // handle preserved
        ws2.rebuildAdjacency();
        REQUIRE(ws2.adjacency().boundaryFaces().size() == 4);
    }
    // 5. METADATA commit (classification field write): zero allocations.
    {
        NodePool np3; std::vector<NodeHandle> m;
        for (int i = 0; i < 4; ++i)
            m.push_back(np3.add({double(i), 2, 0}, {3, 1}));
        CellPool cp3;
        auto c = cp3.add(ElementType::Tet4,
                         std::vector<NodeHandle>{m[0], m[1], m[2], m[3]},
                         {3, 1}, np3);
        EdgePool ep3; ep3.insertCellEdges(np3, cp3, c);
        FacePool fp3; fp3.insertCellFaces(np3, cp3, ep3, c);
        TopologyWorkspace ws3(std::move(np3), std::move(cp3),
                              std::move(ep3), std::move(fp3));
        gPhase3Allocs = 0;
        auto tx = ws3.beginTransaction();
        tx.setCellClassification(c, paacs::mesh::Classification{3, 77});
        (void)tx.commit();
        std::printf("metadata: phase3 allocs = %ld\n", gPhase3Allocs);
        REQUIRE(gPhase3Allocs == 0);
    }
    // 6+7. NODE lifecycle commits: zero allocations for add and remove.
    {
        NodePool np4; std::vector<NodeHandle> m;
        for (int i = 0; i < 2; ++i)
            m.push_back(np4.add({double(i), 3, 0}, {1, 1}));
        CellPool cp4;
        auto c = cp4.add(ElementType::Segment2,
                         std::vector<NodeHandle>{m[0], m[1]}, {1, 1}, np4);
        EdgePool ep4; ep4.insertCellEdges(np4, cp4, c);
        FacePool fp4;
        TopologyWorkspace ws4(std::move(np4), std::move(cp4),
                              std::move(ep4), std::move(fp4));
        TopologyPolicy pol; pol.permitOrphanNodes = true;
        gPhase3Allocs = 0;
        auto tx = ws4.beginTransaction(pol);
        tx.addNode({5, 5, 5}, {0, 9});
        auto res = tx.commit();
        std::printf("addNode: phase3 allocs = %ld\n", gPhase3Allocs);
        REQUIRE(gPhase3Allocs == 0);
        ws4.rebuildAdjacency();
        gPhase3Allocs = 0;
        auto tx2 = ws4.beginTransaction();
        tx2.removeNode(res.node);
        (void)tx2.commit();
        std::printf("removeNode: phase3 allocs = %ld\n", gPhase3Allocs);
        REQUIRE(gPhase3Allocs == 0);
    }
    std::printf(gBad ? "%d FAILURE(S)\n" : "TRANSACTION PROBE PASSED\n",
                gBad);
    return gBad ? 1 : 0;
}
