// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// PAACS Mesh unit tests - M4-B.5 cell ownership / adjacency. Separate TU
// (reviewer: keep sanitizer test builds fast by splitting topology tests).
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <memory>
#include <span>
#include <thread>
#include <vector>

#include "paacs/topo/Workspace.hpp"

#include "test_harness.hpp"

using namespace paacs;

static std::vector<std::uint64_t> adjacencyTokens(
    const topo::CellAdjacency& adj, const topo::CellPool& cells) {
    std::vector<std::uint64_t> out;
    auto handle = [&](const auto& h) {
        out.push_back(static_cast<std::uint64_t>(h.slot));
        out.push_back(static_cast<std::uint64_t>(h.gen));
        out.push_back(h.pool);
    };
    out.push_back(adj.boundNodePool());
    out.push_back(adj.boundCellPool());
    out.push_back(adj.boundEdgePool());
    out.push_back(adj.boundFacePool());
    out.push_back(static_cast<std::uint64_t>(adj.cellCount()));
    for (const auto& cell : cells.activeBySequence()) {
        handle(cell);
        out.push_back(adj.hasCell(cells, cell) ? 1U : 0U);
        const auto faces = adj.cellFaces(cells, cell);
        out.push_back(static_cast<std::uint64_t>(faces.size()));
        for (const auto& face : faces) handle(face);
        const auto neighbors = adj.neighbors(cells, cell);
        out.push_back(static_cast<std::uint64_t>(neighbors.size()));
        for (const auto& neighbor : neighbors) {
            handle(neighbor.cell);
            handle(neighbor.via);
        }
    }
    out.push_back(static_cast<std::uint64_t>(adj.boundaryFaces().size()));
    for (const auto& boundary : adj.boundaryFaces()) {
        handle(boundary.face);
        handle(boundary.owner);
    }
    out.push_back(static_cast<std::uint64_t>(adj.interiorFaces().size()));
    for (const auto& interior : adj.interiorFaces()) {
        handle(interior.face);
        handle(interior.a);
        handle(interior.b);
    }
    out.push_back(static_cast<std::uint64_t>(adj.nonManifoldFaces().size()));
    for (const auto& nonManifold : adj.nonManifoldFaces()) {
        handle(nonManifold.face);
        out.push_back(static_cast<std::uint64_t>(nonManifold.cells.size()));
        for (const auto& cell : nonManifold.cells) handle(cell);
    }
    return out;
}

// ---------------------------------------------------------------------------
static void testCellAdjacency() {
    TEST("topo: M4-B.5 cell adjacency (relations, ownership, invalidation)");
    using topo::NodePool; using topo::CellPool; using topo::EdgePool;
    using topo::FacePool; using topo::CellAdjacency; using topo::ElementType;
    using topo::NodeHandle; using topo::CellHandle; using topo::FaceClass;
    using topo::TopologyWorkspace;

    auto prep = [](NodePool& np, std::vector<NodeHandle>& n, int nn) {
        for (int i=0;i<nn;++i)
            n.push_back(np.add({double(i),double(i%3),double(i%2)},{3,1}));
    };
    auto addCell = [](NodePool& np, CellPool& cp, EdgePool& ep, FacePool& fp,
                      auto type, std::initializer_list<NodeHandle> hs) {
        const std::vector<NodeHandle> v(hs);
        auto c = cp.add(type, v, {3,1}, np);
        ep.insertCellEdges(np,cp,c);
        fp.insertCellFaces(np,cp,ep,c);
        return c;
    };

    // (A) [claim: shared face => ONE neighbor relation; boundary faces have
    // one owner] Two Tet4 sharing one face.
    {
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,5);
        auto cA = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[1],n[2],n[3]});
        auto cB = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[2],n[1],n[4]});
        auto adj = CellAdjacency::build(np,cp,ep,fp);
        CHECK(adj.cellCount()==2);
        auto nA = adj.neighbors(cp,cA); auto nB = adj.neighbors(cp,cB);
        CHECK(nA.size()==1 && nA[0].cell==cB);       // exactly one relation
        CHECK(nB.size()==1 && nB[0].cell==cA);
        CHECK(nA[0].via==nB[0].via);                 // through the same face
        CHECK(adj.boundaryFaces().size()==6);
        for (const auto& bf : adj.boundaryFaces()) { // one owner each
            CHECK(cp.isCurrent(bf.owner));
            CHECK(fp.classify(bf.face)==FaceClass::Boundary);
        }
        CHECK(adj.interiorFaces().size()==1);
        CHECK(adj.nonManifoldFaces().empty());
        CHECK(adj.checkValid(np,cp,ep,fp));
        CHECK(topo::cellAdjacencyReady(np,cp,ep,fp,adj));
    }

    // (B) [claim: interior faces have two owners with OPPOSITE orientation]
    {
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,5);
        auto cA = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[1],n[2],n[3]});
        auto cB = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[2],n[1],n[4]});
        (void)cA; (void)cB;
        auto adj = CellAdjacency::build(np,cp,ep,fp);
        CHECK(adj.interiorFaces().size()==1);
        const auto& itf = adj.interiorFaces()[0];
        CHECK(cp.isCurrent(itf.a) && cp.isCurrent(itf.b));
        const auto fd = fp.get(itf.face);            // frozen face layer
        CHECK(fd.incident.size()==2);
        CHECK(fd.incident[0].sense + fd.incident[1].sense == 0);  // opposite
    }

    // (C) [claim: nonmanifold faces explicitly represented, NO fabricated
    // pairwise neighbor relations] Three tets on one face.
    {
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,6);
        auto c1 = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[1],n[2],n[3]});
        auto c2 = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[1],n[2],n[4]});
        auto c3 = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[1],n[2],n[5]});
        auto adj = CellAdjacency::build(np,cp,ep,fp);
        CHECK(adj.nonManifoldFaces().size()==1);
        CHECK(adj.nonManifoldFaces()[0].cells.size()==3);
        CHECK(adj.neighbors(cp,c1).empty());         // no fabricated pairs
        CHECK(adj.neighbors(cp,c2).empty());
        CHECK(adj.neighbors(cp,c3).empty());
        CHECK(adj.interiorFaces().empty());
        CHECK(adj.checkValid(np,cp,ep,fp));          // nonmanifold is VALID
    }

    // (D) [claim: removing a cell before removing adjacency is DETECTED]
    {
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,5);
        auto cA = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[1],n[2],n[3]});
        auto cB = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[2],n[1],n[4]});
        (void)cA;
        auto adj = CellAdjacency::build(np,cp,ep,fp);
        fp.removeCellFaces(cp,cB); ep.removeCellEdges(cp,cB); cp.remove(cB);
        CHECK(!adj.checkValid(np,cp,ep,fp));
        const auto rep = topo::auditCellAdjacency(np,cp,ep,fp,adj);
        bool ghost=false;
        for (const auto& d : rep.items)
            ghost |= d.code==diag::Code::GhostAdjacency;
        CHECK(ghost);                                // specific diagnostic
        CHECK(adj.neighbors(cp,cB).empty());         // stale handle: nothing
        CHECK(!adj.hasCell(cp,cB));
    }

    // (C2) [claim: the approved capacity snapshot reports every exact
    // retained category, saturates per-field, and is an immutable read]
    {
        {
            CellAdjacency empty;
            const auto snapshot = empty.capacitySnapshot();
            CHECK(snapshot.complete);
            CHECK(snapshot.entryCount == 0);
            CHECK(snapshot.entryFaceCapacity == 0);
            CHECK(snapshot.entryNeighborCapacity == 0);
            CHECK(snapshot.boundaryCapacity == 0);
            CHECK(snapshot.interiorCapacity == 0);
            CHECK(snapshot.nonManifoldCapacity == 0);
            CHECK(snapshot.nonManifoldCellCapacity == 0);
        }
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,6);
        (void)addCell(np,cp,ep,fp,ElementType::Tet4,
                      {n[0],n[1],n[2],n[3]});
        (void)addCell(np,cp,ep,fp,ElementType::Tet4,
                      {n[0],n[1],n[2],n[4]});
        (void)addCell(np,cp,ep,fp,ElementType::Tet4,
                      {n[0],n[1],n[2],n[5]});
        auto adj = CellAdjacency::build(np,cp,ep,fp);
        const auto logical = adjacencyTokens(adj,cp);
        const auto snapshot = adj.capacitySnapshot();
        CHECK(snapshot.complete);
        CHECK(snapshot.entryBucketCount > 0);
        CHECK(snapshot.entryCount == 3);
        CHECK(snapshot.entryFaceCapacity >= 12);
        CHECK(snapshot.entryNeighborCapacity == 0);
        CHECK(snapshot.boundaryCapacity >= adj.boundaryFaces().size());
        CHECK(snapshot.interiorCapacity >= adj.interiorFaces().size());
        CHECK(snapshot.nonManifoldCapacity >=
              adj.nonManifoldFaces().size());
        CHECK(snapshot.nonManifoldCellCapacity >= 3);

        auto sameSnapshot = [](const auto& a, const auto& b) {
            return a.entryBucketCount == b.entryBucketCount &&
                   a.entryCount == b.entryCount &&
                   a.entryFaceCapacity == b.entryFaceCapacity &&
                   a.entryNeighborCapacity == b.entryNeighborCapacity &&
                   a.boundaryCapacity == b.boundaryCapacity &&
                   a.interiorCapacity == b.interiorCapacity &&
                   a.nonManifoldCapacity == b.nonManifoldCapacity &&
                   a.nonManifoldCellCapacity == b.nonManifoldCellCapacity &&
                   a.complete == b.complete;
        };
        CHECK(sameSnapshot(snapshot,adj.capacitySnapshot()));
        std::atomic<bool> concurrentReads{true};
        std::array<std::thread,4> readers;
        for (auto& reader : readers)
            reader = std::thread([&] {
                for (int i=0;i<1000;++i)
                    if (!sameSnapshot(snapshot,adj.capacitySnapshot()))
                        concurrentReads=false;
            });
        for (auto& reader : readers) reader.join();
        CHECK(concurrentReads.load());
        CHECK(adjacencyTokens(adj,cp) == logical);
        CHECK(adj.checkValid(np,cp,ep,fp));

        std::size_t total=(std::numeric_limits<std::size_t>::max)()-2;
        CHECK(!CellAdjacency::testCapacityCheckedAdd(total,3));
        CHECK(total==(std::numeric_limits<std::size_t>::max)());
        topo::faultCountdown()=1;
        CHECK(sameSnapshot(snapshot,adj.capacitySnapshot()));
        CHECK(topo::faultCountdown()==1);
        topo::faultCountdown()=0;

        NodePool npI; CellPool cpI; EdgePool epI; FacePool fpI;
        std::vector<NodeHandle> ni; prep(npI,ni,5);
        (void)addCell(npI,cpI,epI,fpI,ElementType::Tet4,
                      {ni[0],ni[1],ni[2],ni[3]});
        (void)addCell(npI,cpI,epI,fpI,ElementType::Tet4,
                      {ni[0],ni[2],ni[1],ni[4]});
        TopologyWorkspace workspace(std::move(npI),std::move(cpI),
                                    std::move(epI),std::move(fpI));
        const auto semanticBefore=workspace.semanticTopologyHash();
        const auto strictBefore=workspace.strictWorkspaceHash();
        const auto interior=workspace.adjacency().capacitySnapshot();
        CHECK(interior.complete);
        CHECK(interior.entryCount==2);
        CHECK(interior.entryNeighborCapacity>=2);
        CHECK(interior.interiorCapacity>=1);
        CHECK(interior.boundaryCapacity>=6);
        CHECK(interior.nonManifoldCapacity==0);
        CHECK(interior.nonManifoldCellCapacity==0);
        CHECK(workspace.semanticTopologyHash()==semanticBefore);
        CHECK(workspace.strictWorkspaceHash()==strictBefore);
    }

    // (E) [claim: reused cell slots cannot create ghost adjacency]
    {
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,6);
        auto cA = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[1],n[2],n[3]});
        auto cB = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[2],n[1],n[4]});
        (void)cA;
        auto adj = CellAdjacency::build(np,cp,ep,fp);
        fp.removeCellFaces(cp,cB); ep.removeCellEdges(cp,cB); cp.remove(cB);
        // Re-add a cell; the pool reuses the freed slot with a new gen.
        auto cNew = cp.add(ElementType::Tet4,{n[0],n[2],n[1],n[5]},{3,1},np);
        CHECK(cNew.slot == cB.slot);                 // slot actually reused
        CHECK(cNew.gen != cB.gen);
        CHECK(!adj.hasCell(cp,cNew));                // no ghost adjacency
        CHECK(adj.neighbors(cp,cNew).empty());
        CHECK(adj.cellFaces(cp,cNew).empty());
        CHECK(!adj.checkValid(np,cp,ep,fp));         // snapshot stale anyway
    }

    // (F) [claim: stale face handles inside the snapshot are rejected]
    // Remove ONLY one cell's faces (cell stays current): recorded face
    // handles for its exclusive faces become stale.
    {
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,5);
        auto cA = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[1],n[2],n[3]});
        auto cB = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[2],n[1],n[4]});
        (void)cA;
        auto adj = CellAdjacency::build(np,cp,ep,fp);
        fp.removeCellFaces(cp,cB);                   // faces gone, cell stays
        CHECK(!adj.checkValid(np,cp,ep,fp));
        const auto rep = topo::auditCellAdjacency(np,cp,ep,fp,adj);
        CHECK(!rep.acceptable(diag::Severity::Warning));
    }

    // (G) [claim: foreign pools are rejected]
    {
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,5);
        auto cA = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[1],n[2],n[3]});
        auto adj = CellAdjacency::build(np,cp,ep,fp);
        CellPool other;                              // foreign cell pool
        CHECK(!adj.hasCell(other,cA));
        CHECK(adj.neighbors(other,cA).empty());
        NodePool np2; CellPool cp2; EdgePool ep2; FacePool fp2;
        const auto rep = topo::auditCellAdjacency(np2,cp2,ep2,fp2,adj);
        bool bind=false;
        for (const auto& d : rep.items)
            bind |= d.code==diag::Code::AdjacencyPoolBindingMismatch;
        CHECK(bind);
        CHECK(!adj.checkValid(np2,cp2,ep2,fp2));
    }

    // (H) [claim: gate enforced; allocation failure leaves state unchanged]
    {
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,5);
        // Gate: no faces inserted yet for this cell -> build must throw.
        auto c0 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        ep.insertCellEdges(np,cp,c0);
        bool gateThrew=false;
        {
            FacePool empty;
            // faceTopologyReady on an EMPTY face pool is trivially true,
            // but the COMPLETENESS gate (review blocker) must reject: a
            // current 3D cell (c0) has no face incidence.
            bool incThrew=false;
            try { (void)CellAdjacency::build(np,cp,ep,empty); }
            catch (const std::invalid_argument&) { incThrew=true; }
            CHECK(incThrew);
        }
        fp.insertCellFaces(np,cp,ep,c0);
        auto cB = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[2],n[1],n[4]});
        (void)cB;
        // Ghost the face layer (remove a cell without face cleanup in a
        // scratch copy is impossible; instead corrupt via seam) -> build
        // must throw on unready faces.
        fp.testSetHead(fp.testMaterialized());
        try { (void)CellAdjacency::build(np,cp,ep,fp); }
        catch (const std::invalid_argument&) { gateThrew=true; }
        CHECK(gateThrew);
    }
    {
        // Allocation failure: rebuild under fault injection leaves the
        // original adjacency unchanged and still valid.
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,5);
        auto cA = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[1],n[2],n[3]});
        auto cB = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[2],n[1],n[4]});
        auto adj = CellAdjacency::build(np,cp,ep,fp);
        const auto original = adjacencyTokens(adj, cp);
        CHECK(noexcept(adj.swap(adj)));
        for (int stage=1; stage<=10; ++stage) {
            topo::faultCountdown()=stage;
            bool threw=false;
            try { adj.rebuild(np,cp,ep,fp); }
            catch (const std::bad_alloc&) { threw=true; }
            topo::faultCountdown()=0;
            CHECK(adj.checkValid(np,cp,ep,fp));      // unchanged or rebuilt
            CHECK(adj.neighbors(cp,cA).size()==1);
            CHECK(adj.neighbors(cp,cB).size()==1);
            CHECK(adjacencyTokens(adj, cp) == original);
            if (threw) {
                // The same object remains reusable after every failed
                // scratch construction.
                adj.rebuild(np,cp,ep,fp);
                CHECK(adjacencyTokens(adj, cp) == original);
            }
        }
    }

    // (I) [claim: corrupt face incidence does not crash the adjacency
    // audit] All face corruption seams -> audit REPORTS, never crashes.
    {
        auto scenario = [&](auto&& corrupt) {
            NodePool np; CellPool cp; EdgePool ep; FacePool fp;
            std::vector<NodeHandle> n;
            for (int i=0;i<5;++i) n.push_back(np.add({double(i),0,0},{3,1}));
            auto cA = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
            ep.insertCellEdges(np,cp,cA);
            fp.insertCellFaces(np,cp,ep,cA);
            auto adj = CellAdjacency::build(np,cp,ep,fp);
            corrupt(fp);
            const auto rep = topo::auditCellAdjacency(np,cp,ep,fp,adj);
            CHECK(!rep.acceptable(diag::Severity::Warning));
            CHECK(!adj.checkValid(np,cp,ep,fp));
        };
        scenario([](topo::FacePool& f){ f.testSetHead(f.testMaterialized()); });
        scenario([](topo::FacePool& f){ f.testSetTail(-7); });
        scenario([](topo::FacePool& f){ f.testSetNext(f.testHead(), 99999); });
        scenario([](topo::FacePool& f){
            auto v=f.activeBySequence(); f.testDropTableEntry(v[0]); });
        scenario([](topo::FacePool& f){
            auto v=f.activeBySequence(); f.testFlipSense(v[0], 0); });
    }

    // (J) [claim: adjacency-record corruption is caught semantically]
    {
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,5);
        auto cA = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[1],n[2],n[3]});
        auto cB = addCell(np,cp,ep,fp,ElementType::Tet4,{n[0],n[2],n[1],n[4]});
        {   // phantom neighbor
            auto adj = CellAdjacency::build(np,cp,ep,fp);
            adj.testAddPhantomNeighbor(cA.slot,
                {cB, adj.neighbors(cp,cA)[0].via});
            CHECK(!adj.checkValid(np,cp,ep,fp));
            const auto rep = topo::auditCellAdjacency(np,cp,ep,fp,adj);
            bool mm=false;
            for (const auto& d : rep.items)
                mm |= d.code==diag::Code::AdjacencyMismatch;
            CHECK(mm);
        }
        {   // wrong recorded generation
            auto adj = CellAdjacency::build(np,cp,ep,fp);
            adj.testSetEntryGen(cA.slot, cA.gen + 7);
            CHECK(!adj.checkValid(np,cp,ep,fp));
        }
        {   // dropped entry
            auto adj = CellAdjacency::build(np,cp,ep,fp);
            adj.testDropEntry(cB.slot);
            CHECK(!adj.checkValid(np,cp,ep,fp));
        }
    }

    // (K) [claim: mixed elements adjoin through quad faces too] Two Hex8
    // sharing a quad face.
    {
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,12);
        auto h1 = addCell(np,cp,ep,fp,ElementType::Hex8,
            {n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]});
        // Second hex stacked on top, sharing the {4,5,6,7} face reversed.
        auto h2 = addCell(np,cp,ep,fp,ElementType::Hex8,
            {n[4],n[5],n[6],n[7],n[8],n[9],n[10],n[11]});
        auto adj = CellAdjacency::build(np,cp,ep,fp);
        auto n1 = adj.neighbors(cp,h1); auto n2 = adj.neighbors(cp,h2);
        CHECK(n1.size()==1 && n1[0].cell==h2);
        CHECK(n2.size()==1 && n2[0].cell==h1);
        CHECK(adj.interiorFaces().size()==1);
        CHECK(adj.boundaryFaces().size()==10);       // 6+6-2 shared
        CHECK(adj.checkValid(np,cp,ep,fp));
    }

    // (L) [claim: readiness = internal consistency + COMPLETE incidence
    // coverage; review blocker] Empty and partial face incidence reject;
    // completing it succeeds; 2D no-face cells are a documented skip.
    {
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,5);
        auto c1 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        auto c2 = cp.add(ElementType::Tet4,{n[0],n[2],n[1],n[4]},{3,1},np);
        ep.insertCellEdges(np,cp,c1); ep.insertCellEdges(np,cp,c2);
        // Empty FacePool: faceTopologyReady trivially true, build rejects.
        CHECK(topo::faceTopologyReady(np,cp,ep,fp));
        bool emptyThrew=false;
        try { (void)CellAdjacency::build(np,cp,ep,fp); }
        catch (const std::invalid_argument&) { emptyThrew=true; }
        CHECK(emptyThrew);
        // Partial: faces only for c1 -> build rejects; audit emits the
        // specific AdjacencyIncompleteIncidence code.
        fp.insertCellFaces(np,cp,ep,c1);
        bool partialThrew=false;
        try { (void)CellAdjacency::build(np,cp,ep,fp); }
        catch (const std::invalid_argument&) { partialThrew=true; }
        CHECK(partialThrew);
        {
            CellAdjacency none;
            const auto rep = topo::auditCellAdjacency(np,cp,ep,fp,none);
            bool inc=false;
            for (const auto& d : rep.items)
                inc |= d.code==diag::Code::AdjacencyIncompleteIncidence;
            CHECK(inc);
        }
        // Completing incidence makes build succeed.
        fp.insertCellFaces(np,cp,ep,c2);
        auto adj = CellAdjacency::build(np,cp,ep,fp);
        CHECK(adj.cellCount()==2);
        CHECK(adj.neighbors(cp,c1).size()==1);
        CHECK(adj.checkValid(np,cp,ep,fp));
        CHECK(topo::cellAdjacencyReady(np,cp,ep,fp,adj));
    }
    {
        // 2D documented skip: a Tri3 (faceCount 0) never gates adjacency
        // and carries no entry.
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,7);
        auto tri = cp.add(ElementType::Tri3,{n[0],n[1],n[2]},{2,1},np);
        auto tet = cp.add(ElementType::Tet4,{n[3],n[4],n[5],n[6]},{3,1},np);
        ep.insertCellEdges(np,cp,tri); ep.insertCellEdges(np,cp,tet);
        fp.insertCellFaces(np,cp,ep,tet);
        auto adj = CellAdjacency::build(np,cp,ep,fp);
        CHECK(adj.cellCount()==1);
        CHECK(!adj.hasCell(cp,tri));
        CHECK(adj.hasCell(cp,tet));
        CHECK(adj.checkValid(np,cp,ep,fp));
    }

    // (M) [claim: move construction/assignment qualification] Moved-to
    // adjacency answers identically; moved-from is safely not-ready;
    // rebuild after move works.
    {
        NodePool np; CellPool cp; EdgePool ep; FacePool fp;
        std::vector<NodeHandle> n; prep(np,n,5);
        auto cA = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        auto cB = cp.add(ElementType::Tet4,{n[0],n[2],n[1],n[4]},{3,1},np);
        ep.insertCellEdges(np,cp,cA); ep.insertCellEdges(np,cp,cB);
        fp.insertCellFaces(np,cp,ep,cA); fp.insertCellFaces(np,cp,ep,cB);
        CellAdjacency src = CellAdjacency::build(np,cp,ep,fp);
        CellAdjacency dst(std::move(src));
        CHECK(dst.checkValid(np,cp,ep,fp));            // moved-to valid
        CHECK(dst.hasCell(cp,cA));
        CHECK(dst.neighbors(cp,cA).size()==1);
        CHECK(dst.neighbors(cp,cA)[0].cell==cB);
        CHECK(dst.boundaryFaces().size()==6);
        CHECK(topo::cellAdjacencyReady(np,cp,ep,fp,dst));
        // Moved-from: safely not-ready, answers nothing, never crashes.
        CHECK(!src.checkValid(np,cp,ep,fp));
        CHECK(!src.hasCell(cp,cA));
        CHECK(src.neighbors(cp,cA).empty());
        CHECK(!topo::cellAdjacencyReady(np,cp,ep,fp,src));
        // Rebuild after move works, both directions.
        src.rebuild(np,cp,ep,fp);
        CHECK(src.checkValid(np,cp,ep,fp));
        dst.rebuild(np,cp,ep,fp);
        CHECK(dst.checkValid(np,cp,ep,fp));
    }

    // (N) [BND-2A Increment-5 reviewed amendment] A failed workspace
    // adjacency rebuild preserves the stale state, previous snapshot,
    // topology revision and old-snapshot attribution. Every injected
    // direct-build stage is independently retried on the same workspace.
    {
        int failures = 0;
        for (int stage = 1; stage <= 14; ++stage) {
            NodePool np;
            std::vector<NodeHandle> n;
            prep(np, n, 5);
            CellPool cp;
            const auto c1 = cp.add(
                ElementType::Tet4, {n[0], n[1], n[2], n[3]}, {3, 1}, np);
            EdgePool ep;
            ep.insertCellEdges(np, cp, c1);
            FacePool fp;
            fp.insertCellFaces(np, cp, ep, c1);
            auto ws = std::make_unique<topo::TopologyWorkspace>(
                std::move(np), std::move(cp), std::move(ep), std::move(fp));

            const auto oldAttribution = ws->adjacencyBuiltFromRevision();
            auto tx = ws->beginTransaction();
            const std::array<NodeHandle, 4> addedNodes{
                n[0], n[2], n[1], n[4]};
            (void)tx.addCell(
                ElementType::Tet4, addedNodes, {3, 2});
            const auto committed = tx.commit();
            CHECK(committed.topologyChanged);
            CHECK(ws->state() ==
                  topo::WorkspaceState::CommittedAdjacencyStale);
            const auto revision = ws->topologyRevision();
            const auto semantic = ws->semanticTopologyHash();
            const auto strict = ws->strictWorkspaceHash();
            const auto oldSnapshot =
                adjacencyTokens(ws->debugStaleAdjacency(), ws->cells());

            topo::faultCountdown() = stage;
            bool threw = false;
            try {
                ws->rebuildAdjacency();
            } catch (const std::bad_alloc&) {
                threw = true;
                ++failures;
            }
            topo::faultCountdown() = 0;
            if (threw) {
                CHECK(ws->state() ==
                      topo::WorkspaceState::CommittedAdjacencyStale);
                CHECK(ws->topologyRevision() == revision);
                CHECK(ws->adjacencyBuiltFromRevision() == oldAttribution);
                CHECK(adjacencyTokens(ws->debugStaleAdjacency(),
                                      ws->cells()) == oldSnapshot);
                CHECK(ws->semanticTopologyHash() == semantic);
                CHECK(ws->strictWorkspaceHash() == strict);
                ws->rebuildAdjacency();
            }
            CHECK(ws->state() == topo::WorkspaceState::Ready);
            CHECK(ws->topologyRevision() == revision);
            CHECK(ws->adjacencyBuiltFromRevision() == revision);
            CHECK(ws->adjacencyReady());
            CHECK(ws->adjacency().checkValid(
                ws->nodes(), ws->cells(), ws->edges(), ws->faces()));
        }
        CHECK(failures >= 10);
    }
}

// Registration hook consumed by test_main.cpp.
namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* adjacencyTests(std::size_t& count) {
    static const Entry kAdj[] = {
        {"topo.adjacency",            testCellAdjacency},
    };
    count = sizeof(kAdj) / sizeof(kAdj[0]);
    return kAdj;
}
} // namespace paacs_tests
