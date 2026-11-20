// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// PAACS Mesh unit tests - M4-B.6A stages 1-2: prepared-mutation
// foundations and the workspace transaction state machine. Separate TU
// (keeps sanitizer builds fast).
#include <cstdio>
#include <bit>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <vector>

#include "paacs/topo/Workspace.hpp"
#include "paacs/topo/Cavity.hpp"
#include "paacs/topo/Compaction.hpp"
#include "workspace_digest.hpp"

#include "test_harness.hpp"

using namespace paacs;

// Probe + digest come from the shared testing header.

// Non-friend detection concept for the inaccessibility proof
// [topo.transaction.prepared_access_not_public].
template <class A>
concept PreparedPubliclyCallable = requires(A& a) {
    paacs::topo::detail::TransactionAccess::allocPrepared(a);
};

// ---------------------------------------------------------------------------
static void testTransactionFoundations() {
    TEST("topo: M4-B.6A stage 1 - capacity + prepared-access foundations");
    using topo::HandleAllocator; using topo::CellTag; using topo::CellHandle;
    using Probe = topo::detail::TransactionAccessTestProbe;

    // [topo.transaction.prepared_access_not_public] Compile-time proof:
    // from this NON-friend context the private primitives are
    // inaccessible (access failure makes the requires-expression
    // unsatisfied). The testing probe (a sanctioned friend) can call them.
    {
        // Access failure must land in the concept's immediate context so
        // it evaluates FALSE (a non-dependent requires-expression would
        // make the private access a hard error instead).
        constexpr bool publiclyCallable =
            PreparedPubliclyCallable<HandleAllocator<CellTag>>;
        static_assert(!publiclyCallable,
            "prepared primitives must not be callable outside "
            "TopologyTransaction/test probe");
        CHECK(!publiclyCallable);
    }

    // [topo.transaction.slot_capacity_prepared] reserveMaterializedCapacity
    // changes NO accounting, sequences, or existing handles - only
    // capacity - and the reserved capacity is honored by fresh prepared
    // allocs (the testing precondition aborts on any unprepared push).
    {
        HandleAllocator<CellTag> a;
        std::vector<CellHandle> pre;
        for (int i = 0; i < 5; ++i) pre.push_back(a.alloc());
        a.free(pre[1]);                     // one reusable slot
        const auto mat = a.materializedSlots();
        const auto liv = a.liveCount();
        const auto fre = a.reusableFreeSlots();
        const auto ret = a.retiredSlots();
        const auto tot = a.totalCreated();
        a.reserveMaterializedCapacity(mat + 64);
        CHECK(a.materializedSlots() == mat);         // NOT materialized
        CHECK(a.liveCount() == liv);
        CHECK(a.reusableFreeSlots() == fre);
        CHECK(a.retiredSlots() == ret);
        CHECK(a.totalCreated() == tot);              // sequences untouched
        for (const auto& h : pre)
            if (h.slot != pre[1].slot) CHECK(a.isCurrent(h));
        // Fresh prepared allocs into the reserved capacity: no abort.
        a.reserveFreeCapacity(65);
        std::vector<CellHandle> fresh;
        for (int i = 0; i < 64; ++i) fresh.push_back(Probe::alloc(a));
        CHECK(a.liveCount() == liv + 64);
        for (const auto& h : fresh) CHECK(a.isCurrent(h));
        // Prepared frees into reserved free-list capacity: no abort.
        for (const auto& h : fresh) Probe::free(a, h);
        CHECK(a.liveCount() == liv);
    }

    // [topo.transaction.reused_handle_commit_noalloc /
    //  fresh_handle_commit_noalloc - in-suite form] Prepared alloc pops
    // the free list first (reuse), then uses reserved fresh capacity;
    // counters and sequence semantics match the public allocator.
    {
        HandleAllocator<CellTag> a;
        auto h0 = a.alloc();
        const auto seq0 = a.seqOf(h0);
        a.free(h0);                          // slot 0 reusable
        a.reserveMaterializedCapacity(a.materializedSlots() + 2);
        a.reserveFreeCapacity(4);
        auto r = Probe::alloc(a);            // must REUSE slot 0
        CHECK(r.slot == h0.slot);
        CHECK(r.gen == h0.gen + 1);          // generation advanced by free
        CHECK(a.seqOf(r) > seq0);            // fresh sequence
        auto f1 = Probe::alloc(a);           // fresh slot 1
        auto f2 = Probe::alloc(a);           // fresh slot 2
        CHECK(f1.slot == 1 && f2.slot == 2);
        CHECK(a.materializedSlots() == 3);
        CHECK(a.isCurrent(r) && a.isCurrent(f1) && a.isCurrent(f2));
    }

    // [topo.transaction.slot_exhaustion_before_commit - stage-1 form]
    // Beyond-INT32_MAX capacity requests throw BEFORE any state change.
    {
        HandleAllocator<CellTag> a;
        auto h = a.alloc();
        bool threw = false;
        try {
            a.reserveMaterializedCapacity(
                static_cast<std::size_t>(
                    std::numeric_limits<std::int32_t>::max()) + 1);
        } catch (const std::length_error&) { threw = true; }
        CHECK(threw);
        CHECK(a.materializedSlots() == 1);
        CHECK(a.isCurrent(h));
        CHECK(a.totalCreated() == 1);
    }

    // [topo.transaction.sequence_exhaustion_before_commit - stage-1 form]
    // The Phase-2 precheck path: at forced sequence exhaustion the
    // ordinary allocator THROWS before mutation, and the condition is
    // observable to the planner (totalCreated == UINT64_MAX). The
    // prepared primitive additionally carries an aborting testing
    // precondition for the same condition (unreachable in-suite by
    // design; asserted by code inspection + the abort seam).
    {
        HandleAllocator<CellTag> a;
        auto h = a.alloc();
        a.testForceSequenceMax();
        CHECK(a.totalCreated() ==
              std::numeric_limits<std::uint64_t>::max());   // plannable
        bool threw = false;
        try { (void)a.alloc(); } catch (const std::length_error&) { threw = true; }
        if (!threw) { try { (void)a.alloc(); } catch (...) { threw = true; } }
        CHECK(threw);
        CHECK(a.isCurrent(h));               // no mutation on failure
        CHECK(a.liveCount() == 1);
    }

    // [topo.transaction.prepared_commit_primitives_noexcept] The
    // static_asserts inside the probe prove noexcept at compile time;
    // exercise both at runtime for behavior.
    {
        HandleAllocator<CellTag> a;
        a.reserveMaterializedCapacity(1);
        a.reserveFreeCapacity(1);
        auto h = Probe::alloc(a);
        CHECK(a.isCurrent(h));
        Probe::free(a, h);
        CHECK(!a.isCurrent(h));
        CHECK(a.reusableFreeSlots() == 1);
    }
}

// ---------------------------------------------------------------------------
static void testTransactionStateMachine() {
    TEST("topo: M4-B.6A stage 2 - workspace transaction state machine");
    using topo::NodePool; using topo::CellPool; using topo::EdgePool;
    using topo::FacePool; using topo::TopologyWorkspace;
    using topo::WorkspaceState; using topo::ElementType;
    using topo::NodeHandle;

    auto makeWorkspace = []() {
        NodePool np; std::vector<NodeHandle> n;
        for (int i = 0; i < 5; ++i)
            n.push_back(np.add({double(i), double(i % 2), 0}, {3, 1}));
        CellPool cp;
        auto cA = cp.add(ElementType::Tet4,
                         std::vector<NodeHandle>{n[0],n[1],n[2],n[3]},
                         {3,1}, np);
        auto cB = cp.add(ElementType::Tet4,
                         std::vector<NodeHandle>{n[0],n[2],n[1],n[4]},
                         {3,1}, np);
        EdgePool ep;
        ep.insertCellEdges(np, cp, cA); ep.insertCellEdges(np, cp, cB);
        FacePool fp;
        fp.insertCellFaces(np, cp, ep, cA);
        fp.insertCellFaces(np, cp, ep, cB);
        return TopologyWorkspace(std::move(np), std::move(cp),
                                 std::move(ep), std::move(fp));
    };

    // Adoption: Ready, revision 1, adjacency current and queryable.
    {
        auto ws = makeWorkspace();
        CHECK(ws.state() == WorkspaceState::Ready);
        CHECK(ws.topologyRevision() == 1);
        CHECK(ws.adjacencyReady());
        CHECK(ws.adjacency().cellCount() == 2);
        CHECK(ws.nodes().stats().live == 5);
        (void)ws.semanticTopologyHash();
        (void)ws.strictWorkspaceHash();
    }

    // [topo.transaction.single_active] exactly one active transaction;
    // no reuse after commit/abandon.
    {
        auto ws = makeWorkspace();
        auto tx = ws.beginTransaction();
        CHECK(ws.state() == WorkspaceState::TransactionStaging);
        bool second = false;
        try { auto tx2 = ws.beginTransaction(); }
        catch (const std::logic_error&) { second = true; }
        CHECK(second);
        auto res = tx.commit();
        CHECK(!res.topologyChanged);
        CHECK(res.revision == 1);
        bool reuse = false;
        try { tx.abandon(); } catch (const std::logic_error&) { reuse = true; }
        CHECK(reuse);                        // no reuse after commit
        auto tx3 = ws.beginTransaction();    // new transaction now fine
        tx3.abandon();
    }

    // [topo.transaction.query_rejected_while_active] visibility contract:
    // every workspace topology read rejects during an active transaction.
    {
        auto ws = makeWorkspace();
        auto tx = ws.beginTransaction();
        int rejected = 0;
        auto expectReject = [&](auto&& fn) {
            try { fn(); } catch (const std::logic_error&) { ++rejected; }
        };
        expectReject([&]{ (void)ws.nodes(); });
        expectReject([&]{ (void)ws.cells(); });
        expectReject([&]{ (void)ws.edges(); });
        expectReject([&]{ (void)ws.faces(); });
        expectReject([&]{ (void)ws.adjacency(); });
        expectReject([&]{ (void)ws.debugStaleAdjacency(); });
        expectReject([&]{ (void)ws.semanticTopologyHash(); });
        expectReject([&]{ (void)ws.strictWorkspaceHash(); });
        expectReject([&]{ ws.rebuildAdjacency(); });
        CHECK(rejected == 9);
        CHECK(!ws.adjacencyReady());         // reports unavailable
        CHECK(tx.active());                  // preview API allowed
        (void)tx.policy();
        tx.abandon();
    }

    // [topo.transaction.abandon_restores_visibility] abandon() and the
    // DESTRUCTOR both restore full visibility unchanged.
    {
        auto ws = makeWorkspace();
        const auto sem = ws.semanticTopologyHash();
        const auto str = ws.strictWorkspaceHash();
        { auto tx = ws.beginTransaction();
          CHECK(ws.state() == WorkspaceState::TransactionStaging);
        }                                    // destructor abandons
        CHECK(ws.state() == WorkspaceState::Ready);
        CHECK(ws.adjacencyReady());
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.strictWorkspaceHash() == str);
        auto tx = ws.beginTransaction();
        tx.abandon();                        // explicit abandon
        CHECK(ws.state() == WorkspaceState::Ready);
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.strictWorkspaceHash() == str);
        CHECK(ws.adjacency().cellCount() == 2);
    }

    // [topo.transaction.no_production_sequence_consumption_on_failure]
    // begin + abandon (and empty commit) consume NO production sequences
    // in ANY pool and change neither hash.
    {
        auto ws = makeWorkspace();
        const auto sem = ws.semanticTopologyHash();
        const auto str = ws.strictWorkspaceHash();
        const auto tn = ws.nodes().stats().totalCreated;
        const auto tc = ws.cells().stats().totalCreated;
        const auto te = ws.edges().stats().totalCreated;
        const auto tf = ws.faces().stats().totalCreated;
        for (int i = 0; i < 3; ++i) {
            auto tx = ws.beginTransaction();
            if (i % 2) tx.abandon(); else (void)tx.commit();
        }
        CHECK(ws.nodes().stats().totalCreated == tn);
        CHECK(ws.cells().stats().totalCreated == tc);
        CHECK(ws.edges().stats().totalCreated == te);
        CHECK(ws.faces().stats().totalCreated == tf);
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.strictWorkspaceHash() == str);
        CHECK(ws.topologyRevision() == 1);
    }

    // Hash split sanity: semantic and strict differ (strict includes
    // slot/gen), and both are deterministic across recomputation.
    {
        auto ws = makeWorkspace();
        CHECK(ws.semanticTopologyHash() != ws.strictWorkspaceHash());
        CHECK(ws.semanticTopologyHash() == ws.semanticTopologyHash());
        CHECK(ws.strictWorkspaceHash() == ws.strictWorkspaceHash());
    }

    // [topo.transaction.workspace_rejects_missing_edge_incidence /
    //  workspace_rejects_partial_edge_incidence /
    //  workspace_accepts_complete_2d_incidence] Elevated rule: the
    // constructor invariant covers EDGE incidence for cells with edges
    // (1D/2D included), not only the 3D face rule (fourth recurrence of
    // readiness-vs-coverage, caught by the previous hash test's own
    // invalid state).
    {
        auto buildPools = [](int withEdges, NodePool& np, CellPool& cp,
                             EdgePool& ep) {
            std::vector<NodeHandle> n;
            for (int i = 0; i < 4; ++i)
                n.push_back(np.add({double(i), 0, 0}, {1, 1}));
            auto c1 = cp.add(ElementType::Segment2,
                             std::vector<NodeHandle>{n[0], n[1]}, {1, 1}, np);
            auto c2 = cp.add(ElementType::Segment2,
                             std::vector<NodeHandle>{n[2], n[3]}, {1, 1}, np);
            if (withEdges >= 1) ep.insertCellEdges(np, cp, c1);
            if (withEdges >= 2) ep.insertCellEdges(np, cp, c2);
        };
        auto adoptionRejected = [&](int withEdges) {
            NodePool np; CellPool cp; EdgePool ep; FacePool fp;
            buildPools(withEdges, np, cp, ep);
            try {
                TopologyWorkspace ws(std::move(np), std::move(cp),
                                     std::move(ep), std::move(fp));
                return false;
            } catch (const std::invalid_argument&) { return true; }
        };
        CHECK(adoptionRejected(0));          // missing all edge incidence
        CHECK(adoptionRejected(1));          // partial edge incidence
        {                                    // complete 2D incidence
            NodePool np; CellPool cp; EdgePool ep; FacePool fp;
            buildPools(2, np, cp, ep);
            TopologyWorkspace ws(std::move(np), std::move(cp),
                                 std::move(ep), std::move(fp));
            CHECK(ws.adjacencyReady());      // empty adjacency, ready
            CHECK(topo::workspaceIncidenceComplete(
                ws.nodes(), ws.cells(), ws.edges(), ws.faces()));
            (void)ws.semanticTopologyHash();
        }
    }

    // [topo.transaction.hash_detects_edge_incidence_change] Sensitivity
    // proof on a COMPLETE, adopted workspace (reviewer-preferred seam
    // route): swap the two cells' complete incidence lists via the
    // testing seam - nodes, cells, edge records, sequences and use
    // counts all constant, only the cell->edge ASSIGNMENT changes.
    // Both hashes must change; swapping back must restore them exactly.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i)
            n.push_back(np.add({double(i), 0, 0}, {1, 1}));
        CellPool cp;
        auto c1 = cp.add(ElementType::Segment2,
                         std::vector<NodeHandle>{n[0], n[1]}, {1, 1}, np);
        auto c2 = cp.add(ElementType::Segment2,
                         std::vector<NodeHandle>{n[2], n[3]}, {1, 1}, np);
        EdgePool ep;
        ep.insertCellEdges(np, cp, c1);
        ep.insertCellEdges(np, cp, c2);      // COMPLETE incidence
        FacePool fp;
        TopologyWorkspace ws(std::move(np), std::move(cp),
                             std::move(ep), std::move(fp));
        // HASH-INSTRUMENT CORRUPTION TESTING (reviewer wording): the
        // seam deliberately creates an invalid semantic assignment; the
        // UNGATED instrument seam proves sensitivity, while the PUBLIC
        // hash accessors are gated on workspaceTopologyReady and must
        // REJECT the same corrupted state.
        const auto sem0 = ws.testHashUngated(false);
        const auto str0 = ws.testHashUngated(true);
        CHECK(ws.semanticTopologyHash() == sem0);   // gated == instrument
        ws.testEdgesMutable().testSwapCellIncidence(c1.slot, c2.slot);
        CHECK(ws.testHashUngated(false) != sem0);   // assignment visible
        CHECK(ws.testHashUngated(true) != str0);
        bool gatedRejected = false;
        try { (void)ws.semanticTopologyHash(); }
        catch (const std::logic_error&) { gatedRejected = true; }
        CHECK(gatedRejected);                       // public path rejects
        ws.testEdgesMutable().testSwapCellIncidence(c1.slot, c2.slot);
        CHECK(ws.testHashUngated(false) == sem0);   // restored exactly
        CHECK(ws.testHashUngated(true) == str0);
        CHECK(ws.strictWorkspaceHash() == str0);    // public path healthy
    }

    // [topo.transaction.strict_hash_detects_sequence_consumption] State C
    // adds and removes an extra node before adoption: live content is
    // IDENTICAL to state D, but a production sequence was consumed.
    // Semantic hashes must match; strict hashes must differ.
    {
        auto makeWs = [](bool consumeSequence) {
            NodePool np;
            std::vector<NodeHandle> n;
            for (int i = 0; i < 4; ++i)
                n.push_back(np.add({double(i), 0, 0}, {3, 1}));
            if (consumeSequence) {
                auto extra = np.add({9, 9, 9}, {3, 1});
                np.remove(extra);            // gone, but sequence consumed
            }
            CellPool cp;
            auto c = cp.add(ElementType::Tet4,
                            std::vector<NodeHandle>{n[0],n[1],n[2],n[3]},
                            {3,1}, np);
            EdgePool ep; ep.insertCellEdges(np, cp, c);
            FacePool fp; fp.insertCellFaces(np, cp, ep, c);
            return TopologyWorkspace(std::move(np), std::move(cp),
                                     std::move(ep), std::move(fp));
        };
        auto wsC = makeWs(true);
        auto wsD = makeWs(false);
        CHECK(wsC.semanticTopologyHash() == wsD.semanticTopologyHash());
        CHECK(wsC.strictWorkspaceHash() != wsD.strictWorkspaceHash());
    }
}

// ---------------------------------------------------------------------------
static void testTransactionCellOps() {
    TEST("topo: M4-B.6A stage 3 - atomic add/remove cell");
    using topo::NodePool; using topo::CellPool; using topo::EdgePool;
    using topo::FacePool; using topo::TopologyWorkspace;
    using topo::WorkspaceState; using topo::ElementType;
    using topo::NodeHandle; using topo::CellHandle;

    struct Fixture {
        std::vector<NodeHandle> n;
        CellHandle c1;
        std::unique_ptr<TopologyWorkspace> ws;
        Fixture() {
            NodePool np;
            for (int i = 0; i < 5; ++i)
                n.push_back(np.add({double(i), double(i % 2), double(i % 3)},
                                   {3, 1}));
            CellPool cp;
            c1 = cp.add(ElementType::Tet4,
                        std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                        {3, 1}, np);
            EdgePool ep; ep.insertCellEdges(np, cp, c1);
            FacePool fp; fp.insertCellFaces(np, cp, ep, c1);
            ws = std::make_unique<TopologyWorkspace>(
                std::move(np), std::move(cp), std::move(ep), std::move(fp));
        }
    };
    const std::vector<NodeHandle> kNoNodes;

    // [topo.transaction.add_cell_atomic + boundary_to_interior +
    //  edge_face_incidence_recomputed + adjacency_revision_invalidated +
    //  postcommit_workspace_ready] Add a second tet sharing face
    //  {n0,n1,n2} with opposite winding.
    {
        Fixture f; auto& ws = *f.ws;
        CHECK(ws.adjacency().boundaryFaces().size() == 4);
        const auto rev0 = ws.topologyRevision();
        auto tx = ws.beginTransaction();
        (void)tx.addCell(ElementType::Tet4,
                         std::vector<NodeHandle>{f.n[0], f.n[2], f.n[1],
                                                 f.n[4]},
                         {3, 2});
        auto res = tx.commit();
        CHECK(res.topologyChanged);
        CHECK(res.revision == rev0 + 1);
        CHECK(ws.topologyRevision() == rev0 + 1);
        CHECK(ws.state() == WorkspaceState::CommittedAdjacencyStale);
        const CellHandle c2 = res.cell;
        CHECK(ws.cells().isCurrent(c2));
        CHECK(ws.cells().get(c2).cls.tag == 2);
        // [adjacency_revision_invalidated]
        CHECK(!ws.adjacencyReady());
        CHECK(ws.adjacencyBuiltFromRevision() == rev0);
        bool staleRejected = false;
        try { (void)ws.adjacency(); }
        catch (const std::logic_error&) { staleRejected = true; }
        CHECK(staleRejected);
        (void)ws.debugStaleAdjacency();      // explicit accessor still legal
        // [postcommit_workspace_ready]
        CHECK(topo::workspaceTopologyReady(ws.nodes(), ws.cells(),
                                           ws.edges(), ws.faces()));
        // [edge_face_incidence_recomputed]
        CHECK(ws.edges().hasCellEdges(ws.cells(), c2));
        CHECK(ws.edges().cellEdges(ws.cells(), c2).size() == 6);
        CHECK(ws.faces().hasCellFaces(ws.cells(), c2));
        CHECK(ws.faces().cellFaces(ws.cells(), c2).size() == 4);
        CHECK(ws.edges().liveCount() == 9);  // 6 + 6 - 3 shared
        CHECK(ws.faces().liveCount() == 7);  // 4 + 4 - 1 shared
        // [boundary_to_interior] rebuild and verify shared face now has
        // two owners with opposite senses; boundary count 6.
        ws.rebuildAdjacency();
        CHECK(ws.state() == WorkspaceState::Ready);
        CHECK(ws.adjacencyReady());
        CHECK(ws.adjacency().boundaryFaces().size() == 6);
        CHECK(ws.adjacency().neighbors(ws.cells(), f.c1).size() == 1);
        CHECK(ws.adjacency().neighbors(ws.cells(), f.c1)[0].cell == c2);
        auto shared = ws.faces().find(
            ws.nodes(), std::vector<NodeHandle>{f.n[0], f.n[1], f.n[2]});
        CHECK(shared.has_value());
        const auto fd = ws.faces().get(*shared);
        CHECK(fd.incident.size() == 2);
        CHECK(fd.incident[0].sense * fd.incident[1].sense == -1);
        // Both hashes changed deterministically and are recomputable.
        CHECK(ws.semanticTopologyHash() == ws.semanticTopologyHash());

        // [topo.transaction.remove_cell_atomic + interior_to_boundary]
        // (removing c2 legitimately orphans n4 - it exists only for c2 -
        // so this transaction runs under permitOrphanNodes).
        const auto semBefore = ws.semanticTopologyHash();
        topo::TopologyPolicy permitOrphans;
        permitOrphans.permitOrphanNodes = true;
        auto tx2 = ws.beginTransaction(permitOrphans);
        tx2.removeCell(c2);
        auto res2 = tx2.commit();
        CHECK(res2.topologyChanged);
        CHECK(ws.topologyRevision() == rev0 + 2);
        CHECK(!ws.cells().isCurrent(c2));
        CHECK(ws.edges().liveCount() == 6);  // c2-only edges retired
        CHECK(ws.faces().liveCount() == 4);  // c2-only faces retired
        CHECK(!ws.edges().hasCellEdges(ws.cells(), c2));
        CHECK(!ws.faces().hasCellFaces(ws.cells(), c2));
        CHECK(ws.semanticTopologyHash() != semBefore);
        ws.rebuildAdjacency();
        CHECK(ws.adjacency().boundaryFaces().size() == 4);   // back to 1 tet
        const auto fd2 = ws.faces().get(*shared);
        CHECK(fd2.incident.size() == 1);     // interior -> boundary
        CHECK(topo::workspaceTopologyReady(ws.nodes(), ws.cells(),
                                           ws.edges(), ws.faces()));
        // Survivors intact.
        CHECK(ws.cells().isCurrent(f.c1));
        for (const auto& nh : f.n) CHECK(ws.nodes().isCurrent(nh));
    }

    // [topo.transaction.add_cell_failure_each_phase2_stage +
    //  failed_add_hashes_unchanged + failed_operations_consume_no_sequences
    //  + surviving_handles_unchanged] Fault-inject every Phase-2 stage.
    {
        Fixture f; auto& ws = *f.ws;
        const auto sem = ws.semanticTopologyHash();
        const auto str = ws.strictWorkspaceHash();
        const auto tn = ws.nodes().stats().totalCreated;
        const auto tc = ws.cells().stats().totalCreated;
        int failures = 0;
        for (int stage = 1; stage <= 8; ++stage) {
            auto tx = ws.beginTransaction();
            (void)tx.addCell(ElementType::Tet4,
                             std::vector<NodeHandle>{f.n[0], f.n[2], f.n[1],
                                                     f.n[4]},
                             {3, 2});
            topo::faultCountdown() = stage;
            bool threw = false;
            try { (void)tx.commit(); }
            catch (const std::bad_alloc&) { threw = true; }
            topo::faultCountdown() = 0;
            if (threw) {
                ++failures;
                tx.abandon();
            }
            CHECK(ws.state() == WorkspaceState::Ready ||
                  ws.state() == WorkspaceState::CommittedAdjacencyStale);
            if (threw) {
                CHECK(ws.semanticTopologyHash() == sem);
                CHECK(ws.strictWorkspaceHash() == str);
            }
            if (!threw) {                    // countdown outlived Phase 2:
                topo::TopologyPolicy pol;    // undo the committed add
                pol.permitOrphanNodes = true;
                const auto added = ws.cells().activeBySequence().back();
                auto tx3 = ws.beginTransaction(pol);
                tx3.removeCell(added);
                (void)tx3.commit();
                ws.rebuildAdjacency();
            }
        }
        CHECK(failures >= 5);                // every real stage was hit
        CHECK(ws.nodes().stats().totalCreated == tn);
        CHECK(ws.cells().stats().totalCreated >= tc);   // successes consumed
        for (const auto& nh : f.n) CHECK(ws.nodes().isCurrent(nh));
        CHECK(ws.cells().isCurrent(f.c1));
        ws.rebuildAdjacency();
        CHECK(ws.adjacency().cellCount() == 1);
    }

    // Pure-failure sequence accounting: with a stage that ALWAYS fails
    // early (stage 1), sequences and hashes never move.
    {
        Fixture f; auto& ws = *f.ws;
        const auto sem = ws.semanticTopologyHash();
        const auto str = ws.strictWorkspaceHash();
        const auto tc = ws.cells().stats().totalCreated;
        const auto te = ws.edges().stats().totalCreated;
        const auto tf = ws.faces().stats().totalCreated;
        for (int i = 0; i < 3; ++i) {
            auto tx = ws.beginTransaction();
            (void)tx.addCell(ElementType::Tet4,
                             std::vector<NodeHandle>{f.n[0], f.n[2], f.n[1],
                                                     f.n[4]},
                             {3, 2});
            topo::faultCountdown() = 1;
            bool threw = false;
            try { tx.prepare(); } catch (const std::bad_alloc&) { threw = true; }
            topo::faultCountdown() = 0;
            CHECK(threw);
            tx.abandon();
        }
        CHECK(ws.cells().stats().totalCreated == tc);
        CHECK(ws.edges().stats().totalCreated == te);
        CHECK(ws.faces().stats().totalCreated == tf);
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.strictWorkspaceHash() == str);
    }

    // [topo.transaction.remove_cell_failure_each_phase2_stage +
    //  failed_remove_hashes_unchanged]
    {
        Fixture f; auto& ws = *f.ws;
        const auto sem = ws.semanticTopologyHash();
        const auto str = ws.strictWorkspaceHash();
        for (int stage = 1; stage <= 2; ++stage) {
            topo::TopologyPolicy pol;
            pol.permitOrphanNodes = true;
            auto tx2 = ws.beginTransaction(pol);
            tx2.removeCell(f.c1);
            topo::faultCountdown() = stage;
            bool threw = false;
            try { tx2.prepare(); } catch (const std::bad_alloc&) { threw = true; }
            topo::faultCountdown() = 0;
            CHECK(threw);
            tx2.abandon();
            CHECK(ws.cells().isCurrent(f.c1));
            CHECK(ws.semanticTopologyHash() == sem);
            CHECK(ws.strictWorkspaceHash() == str);
        }
        // And the orphan policy itself: default policy rejects in Phase 1.
        auto tx = ws.beginTransaction();
        bool orphanRejected = false;
        try { tx.removeCell(f.c1); }
        catch (const std::invalid_argument&) { orphanRejected = true; }
        CHECK(orphanRejected);
        tx.abandon();
    }

    // [topo.transaction.abandoned_add_has_no_effect +
    //  abandoned_remove_has_no_effect +
    //  queries_reject_during_prepared_state + second_commit_rejected +
    //  commit_after_abandon_rejected]
    {
        Fixture f; auto& ws = *f.ws;
        const auto sem = ws.semanticTopologyHash();
        const auto str = ws.strictWorkspaceHash();
        const auto edgeTable = ws.edges().stats().live;
        {
            auto tx = ws.beginTransaction();
            (void)tx.addCell(ElementType::Tet4,
                             std::vector<NodeHandle>{f.n[0], f.n[2], f.n[1],
                                                     f.n[4]},
                             {3, 2});
            tx.prepare();                    // placeholders LIVE
            CHECK(ws.state() == WorkspaceState::TransactionPrepared);
            int rejected = 0;
            auto expectReject = [&](auto&& fn) {
                try { fn(); } catch (const std::logic_error&) { ++rejected; }
            };
            expectReject([&]{ (void)ws.nodes(); });
            expectReject([&]{ (void)ws.edges(); });
            expectReject([&]{ (void)ws.semanticTopologyHash(); });
            expectReject([&]{ (void)ws.adjacency(); });
            CHECK(rejected == 4);
            tx.abandon();                    // placeholders erased
        }
        CHECK(ws.state() == WorkspaceState::Ready);
        CHECK(ws.edges().stats().live == edgeTable);
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.strictWorkspaceHash() == str);
        {
            topo::TopologyPolicy pol;
            pol.permitOrphanNodes = true;
            auto tx = ws.beginTransaction(pol);
            tx.removeCell(f.c1);
            tx.prepare();
            tx.abandon();                    // abandoned remove: no effect
        }
        CHECK(ws.cells().isCurrent(f.c1));
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.strictWorkspaceHash() == str);
        {
            auto tx = ws.beginTransaction();
            (void)tx.commit();               // empty commit consumes tx
            bool second = false;
            try { (void)tx.commit(); } catch (const std::logic_error&) { second = true; }
            CHECK(second);                   // second_commit_rejected
        }
        {
            auto tx = ws.beginTransaction();
            tx.abandon();
            bool after = false;
            try { (void)tx.commit(); } catch (const std::logic_error&) { after = true; }
            CHECK(after);                    // commit_after_abandon_rejected
        }
        // One operation per transaction (v1 restriction).
        {
            auto tx = ws.beginTransaction();
            (void)tx.addCell(ElementType::Tet4,
                             std::vector<NodeHandle>{f.n[0], f.n[2], f.n[1],
                                                     f.n[4]},
                             {3, 2});
            bool secondOp = false;
            try {
                (void)tx.addCell(ElementType::Tet4,
                                 std::vector<NodeHandle>{f.n[0], f.n[2],
                                                         f.n[1], f.n[4]},
                                 {3, 2});
            } catch (const std::logic_error&) { secondOp = true; }
            CHECK(secondOp);
            tx.abandon();
        }
    }
    // [topo.transaction.successful_first_add_binds_pools +
    //  failed_add_preserves_pool_bindings +
    //  abandoned_add_preserves_pool_bindings +
    //  binding_change_not_hidden_by_hash_evidence]
    // Pool ownership binding is BEHAVIORAL state (review blocker): a
    // failed or abandoned transaction must preserve it, and the hash
    // pair cannot see it - so it is tested DIRECTLY via the probe.
    {
        using Probe = topo::detail::TransactionAccessTestProbe;
        auto makeEmptyWs = [](std::vector<NodeHandle>& n) {
            NodePool np;
            for (int i = 0; i < 2; ++i)
                n.push_back(np.add({double(i), 0, 0}, {1, 1}));
            CellPool cp; EdgePool ep; FacePool fp;   // all UNBOUND
            return std::make_unique<TopologyWorkspace>(
                std::move(np), std::move(cp), std::move(ep), std::move(fp));
        };
        auto allUnbound = [](const auto& b) {
            return b.cellNode == 0 && b.edgeNode == 0 && b.edgeCell == 0 &&
                   b.faceNode == 0 && b.faceCell == 0;
        };
        // Failed Phase-2 preparation preserves (absence of) bindings.
        {
            std::vector<NodeHandle> n;
            auto ws = makeEmptyWs(n);
            for (int stage = 1; stage <= 6; ++stage) {
                auto tx = ws->beginTransaction();
                (void)tx.addCell(ElementType::Segment2,
                                 std::vector<NodeHandle>{n[0], n[1]},
                                 {1, 1});
                topo::faultCountdown() = stage;
                bool threw = false;
                try { tx.prepare(); }
                catch (const std::bad_alloc&) { threw = true; }
                topo::faultCountdown() = 0;
                if (!threw) { tx.abandon(); continue; }
                tx.abandon();
                CHECK(allUnbound(Probe::bindings(ws->cells(), ws->edges(),
                                                 ws->faces())));
            }
        }
        // Abandon after SUCCESSFUL prepare also preserves bindings
        // (binding now happens only in Phase 3).
        {
            std::vector<NodeHandle> n;
            auto ws = makeEmptyWs(n);
            auto tx = ws->beginTransaction();
            (void)tx.addCell(ElementType::Segment2,
                             std::vector<NodeHandle>{n[0], n[1]}, {1, 1});
            tx.prepare();
            tx.abandon();
            CHECK(allUnbound(Probe::bindings(ws->cells(), ws->edges(),
                                             ws->faces())));
            CHECK(topo::workspaceTopologyReady(ws->nodes(), ws->cells(),
                                               ws->edges(), ws->faces()));
        }
        // Successful FIRST add binds all pools; and the hash evidence
        // cannot see the binding change (semantic hash identical across
        // empty-unbound vs emptied-but-bound states), which is exactly
        // why bindings get their own direct proof.
        {
            std::vector<NodeHandle> n;
            auto ws = makeEmptyWs(n);
            const auto sem0 = ws->semanticTopologyHash();
            CHECK(allUnbound(Probe::bindings(ws->cells(), ws->edges(),
                                             ws->faces())));
            auto tx = ws->beginTransaction();
            (void)tx.addCell(ElementType::Segment2,
                             std::vector<NodeHandle>{n[0], n[1]}, {1, 1});
            auto res = tx.commit();
            const auto b = Probe::bindings(ws->cells(), ws->edges(),
                                           ws->faces());
            CHECK(b.cellNode != 0);          // successful_first_add_binds
            CHECK(b.edgeNode != 0 && b.edgeCell != 0);
            CHECK(b.faceNode != 0 && b.faceCell != 0);
            ws->rebuildAdjacency();
            // Remove it again: live content returns to empty...
            topo::TopologyPolicy pol; pol.permitOrphanNodes = true;
            auto tx2 = ws->beginTransaction(pol);
            tx2.removeCell(res.cell);
            (void)tx2.commit();
            ws->rebuildAdjacency();
            // ...semantic hash equals the pre-binding state, while the
            // bindings verifiably CHANGED: hash evidence cannot stand in
            // for direct behavioral-state checks.
            CHECK(ws->semanticTopologyHash() == sem0);
            const auto b2 = Probe::bindings(ws->cells(), ws->edges(),
                                            ws->faces());
            CHECK(b2.cellNode != 0 && b2.edgeNode != 0 && b2.faceNode != 0);
        }
    }
    // [topo.transaction.reused_slots_require_no_materialized_growth]
    // (review blocker: Phase 2 reserved a fresh slot for every logical
    // entity, ignoring free-list reuse). A pure-reuse transaction must
    // not grow allocator storage OR materialized counts across
    // prepare/abandon, and its commit must actually reuse the slots.
    {
        using Probe = topo::detail::TransactionAccessTestProbe;
        std::vector<NodeHandle> n;
        NodePool np;
        for (int i = 0; i < 2; ++i)
            n.push_back(np.add({double(i), 0, 0}, {1, 1}));
        CellPool cp; EdgePool ep; FacePool fp;
        TopologyWorkspace ws(std::move(np), std::move(cp), std::move(ep),
                             std::move(fp));
        topo::TopologyPolicy orphanOk; orphanOk.permitOrphanNodes = true;
        auto tx0 = ws.beginTransaction();
        (void)tx0.addCell(ElementType::Segment2,
                          std::vector<NodeHandle>{n[0], n[1]}, {1, 1});
        auto res0 = tx0.commit();
        ws.rebuildAdjacency();
        auto tr = ws.beginTransaction(orphanOk);
        tr.removeCell(res0.cell);
        (void)tr.commit();
        ws.rebuildAdjacency();
        // One reusable cell slot and one reusable edge slot now exist.
        const auto cellCap = Probe::slotCap(ws.cells());
        const auto edgeCap = Probe::slotCap(ws.edges());
        const auto cellMat = Probe::materialized(ws.cells());
        const auto edgeMat = Probe::materialized(ws.edges());
        {
            auto tx = ws.beginTransaction();
            (void)tx.addCell(ElementType::Segment2,
                             std::vector<NodeHandle>{n[0], n[1]}, {1, 1});
            tx.prepare();
            tx.abandon();
        }
        CHECK(Probe::slotCap(ws.cells()) == cellCap);   // NO growth
        CHECK(Probe::slotCap(ws.edges()) == edgeCap);
        CHECK(Probe::materialized(ws.cells()) == cellMat);
        CHECK(Probe::materialized(ws.edges()) == edgeMat);
        // Commit reuses the freed slots.
        auto tx = ws.beginTransaction();
        (void)tx.addCell(ElementType::Segment2,
                         std::vector<NodeHandle>{n[0], n[1]}, {1, 1});
        auto res = tx.commit();
        CHECK(res.cell.slot == res0.cell.slot);
        CHECK(res.cell.gen == res0.cell.gen + 1);
        CHECK(Probe::materialized(ws.cells()) == cellMat);
    }

    // [topo.transaction.mixed_reused_fresh_capacity_exact +
    //  reuse_allowed_at_materialized_slot_limit +
    //  fresh_slot_exhaustion_rejected_before_prepare]
    // Unit-level EXACT accounting through the probe: the plan states the
    // exact reusable/fresh split, the reservation target equals
    // materialized + fresh, and pure-reuse plans have fresh == 0 (so the
    // reservation target equals materializedSlots(), which early-returns
    // BEFORE the INT32_MAX check - a reusable slot permits a transaction
    // even at the materialized-slot limit; the fresh > 0 boundary throw
    // is the stage-1 slot_exhaustion proof).
    {
        using Probe = topo::detail::TransactionAccessTestProbe;
        NodePool np; std::vector<NodeHandle> n;
        for (int i = 0; i < 9; ++i)
            n.push_back(np.add({double(i), double(i % 2), double(i % 3)},
                               {3, 1}));
        CellPool cp;
        auto t1 = cp.add(ElementType::Tet4,
                         std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                         {3, 1}, np);
        auto t2 = cp.add(ElementType::Tet4,
                         std::vector<NodeHandle>{n[0], n[2], n[1], n[4]},
                         {3, 1}, np);
        (void)t1;
        EdgePool ep; ep.insertCellEdges(np, cp, t1);
        ep.insertCellEdges(np, cp, t2);
        FacePool fp; fp.insertCellFaces(np, cp, ep, t1);
        fp.insertCellFaces(np, cp, ep, t2);
        // Free tet2's unique entities: 1 cell, 3 edges, 3 faces.
        fp.removeCellFaces(cp, t2);
        ep.removeCellEdges(cp, t2);
        cp.remove(t2);
        const auto edgeMat = Probe::materialized(ep);   // 9
        const auto faceMat = Probe::materialized(fp);   // 7
        CHECK(edgeMat == 9 && faceMat == 7);
        // MIXED: a tet on four brand-new nodes needs 6 new edge records
        // (3 reusable slots -> fresh 3) and 4 new faces (3 reusable ->
        // fresh 1); the cell reuses its slot (fresh 0).
        {
            auto plan = Probe::planAdd(np, cp, ep, fp, ElementType::Tet4,
                                       std::vector<NodeHandle>{n[5], n[6],
                                                               n[7], n[8]},
                                       {3, 2});
            Probe::prepareAdd(np, cp, ep, fp, plan);
            CHECK(plan.cellsAcct.requested == 1);
            CHECK(plan.cellsAcct.reusable == 1);
            CHECK(plan.cellsAcct.fresh == 0);
            CHECK(plan.edgesAcct.requested == 6);
            CHECK(plan.edgesAcct.reusable == 3);
            CHECK(plan.edgesAcct.fresh == 3);
            CHECK(plan.edgesAcct.slotCapacity == edgeMat + 3);   // EXACT
            CHECK(plan.facesAcct.requested == 4);
            CHECK(plan.facesAcct.reusable == 3);
            CHECK(plan.facesAcct.fresh == 1);
            CHECK(plan.facesAcct.slotCapacity == faceMat + 1);
            auto ch = Probe::commitAdd(np, cp, ep, fp, plan);
            CHECK(cp.isCurrent(ch));
            // Exactly three edge handles reused old slots, three fresh.
            std::size_t reused = 0, fresh = 0;
            for (const auto& eh : ep.cellEdges(cp, ch))
                (static_cast<std::size_t>(eh.slot) < edgeMat ? reused
                                                             : fresh)++;
            CHECK(reused == 3 && fresh == 3);
            CHECK(Probe::materialized(ep) == edgeMat + 3);
            CHECK(Probe::materialized(fp) == faceMat + 1);
            CHECK(topo::edgeTopologyReady(np, cp, ep));
            CHECK(topo::faceTopologyReady(np, cp, ep, fp));
            // Undo for the next block.
            fp.removeCellFaces(cp, ch);
            ep.removeCellEdges(cp, ch);
            cp.remove(ch);
        }
        // PURE REUSE at the (conceptual) limit: re-adding tet2 needs
        // 3 edge and 3 face creations, all satisfiable from free lists.
        {
            auto plan = Probe::planAdd(np, cp, ep, fp, ElementType::Tet4,
                                       std::vector<NodeHandle>{n[0], n[2],
                                                               n[1], n[4]},
                                       {3, 1});
            Probe::prepareAdd(np, cp, ep, fp, plan);
            CHECK(plan.cellsAcct.fresh == 0);
            CHECK(plan.edgesAcct.fresh == 0);
            CHECK(plan.facesAcct.fresh == 0);
            CHECK(plan.cellsAcct.slotCapacity ==
                  Probe::materialized(cp));   // target == materialized:
            CHECK(plan.edgesAcct.slotCapacity ==
                  Probe::materialized(ep));   // early-returns before the
                                              // INT32_MAX check, so reuse
                                              // is allowed AT the limit
            Probe::abandonAdd(ep, fp, plan);
        }
        // reserveMaterializedCapacity(materialized) is a no-throw no-op
        // (the fresh == 0 reservation form).
        {
            topo::HandleAllocator<topo::CellTag> a;
            (void)a.alloc();
            a.reserveMaterializedCapacity(a.materializedSlots());
            CHECK(a.materializedSlots() == 1);
        }
        // CONDITIONAL sequence checks: with edge AND face sequences
        // forced exhausted, an op that allocates no edge/face handles
        // still plans (Tri3 whose 3 edges all exist); an op needing new
        // edges is rejected with the structured error.
        {
            Probe::forceEdgeSeqMax(ep);
            Probe::forceFaceSeqMax(fp);
            auto plan = Probe::planAdd(np, cp, ep, fp, ElementType::Tri3,
                                       std::vector<NodeHandle>{n[0], n[1],
                                                               n[2]},
                                       {2, 1});
            CHECK(plan.edgesAcct.requested == 0);   // not yet accounted
            bool rejected = false;
            try {
                (void)Probe::planAdd(np, cp, ep, fp, ElementType::Tri3,
                                     std::vector<NodeHandle>{n[0], n[1],
                                                             n[5]},
                                     {2, 1});
            } catch (const std::length_error&) { rejected = true; }
            CHECK(rejected);                 // new edges -> rejected
        }
    }
    (void)kNoNodes;
}

// ---------------------------------------------------------------------------
static void testTransactionReplaceAndPolicy() {
    TEST("topo: M4-B.6A stage 4 - connectivity replacement + policy + validator");
    using topo::NodePool; using topo::CellPool; using topo::EdgePool;
    using topo::FacePool; using topo::TopologyWorkspace;
    using topo::WorkspaceState; using topo::ElementType;
    using topo::NodeHandle; using topo::TransactionPreview;

    auto makeTetWs = [](std::vector<NodeHandle>& n, int extraNodes) {
        NodePool np;
        for (int i = 0; i < 4 + extraNodes; ++i)
            n.push_back(np.add({double(i), double(i % 2), double(i % 3)},
                               {3, 1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4,
                        std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                        {3, 1}, np);
        EdgePool ep; ep.insertCellEdges(np, cp, c);
        FacePool fp; fp.insertCellFaces(np, cp, ep, c);
        return std::pair(std::make_unique<TopologyWorkspace>(
                             std::move(np), std::move(cp), std::move(ep),
                             std::move(fp)),
                         c);
    };

    // [replace_connectivity_atomic + replace_preserves_cell_handle +
    //  replace_preserves_cell_sequence + replace_net_edge_delta_exact +
    //  replace_net_face_delta_exact + replace_postcommit_workspace_ready +
    //  validator_receives_preview_only + validator_not_called_during_commit]
    // Apex swap n3 -> n4: base edges/face retained; apex entities replaced.
    {
        std::vector<NodeHandle> n;
        auto [wsp, c] = makeTetWs(n, 1);
        auto& ws = *wsp;
        topo::TopologyPolicy pol; pol.permitOrphanNodes = true;   // n3 drops
        const auto seq0 = ws.cells().seqOf(c);
        const auto cellsCreated = ws.cells().stats().totalCreated;
        const auto rev0 = ws.topologyRevision();
        auto tx = ws.beginTransaction(pol);
        tx.replaceCellConnectivity(
            c, std::vector<NodeHandle>{n[0], n[1], n[2], n[4]});
        int validatorCalls = 0;
        const bool ok = tx.validate([&](const TransactionPreview& v) {
            ++validatorCalls;
            // Preview VALUES only - exact NET deltas.
            CHECK(v.kind == TransactionPreview::Kind::Replace);
            CHECK(v.cell == c);
            CHECK(v.retainedEdges == 3);         // base 01,02,12
            CHECK(v.newEdgeRecords == 3);        // 04,14,24
            CHECK(v.reusedEdgeIncidence == 0);
            CHECK(v.removedEdgeIncidence == 3);  // 03,13,23
            CHECK(v.retiringEdges == 3);         // sole user
            CHECK(v.retainedFaces == 1);         // base {0,1,2}
            CHECK(v.newFaceRecords == 3);
            CHECK(v.removedFaceIncidence == 3);
            CHECK(v.retiringFaces == 3);
            return true;
        });
        CHECK(ok);
        CHECK(validatorCalls == 1);
        auto res = tx.commit();
        CHECK(validatorCalls == 1);              // NOT called during commit
        CHECK(res.topologyChanged);
        // Handle + sequence PRESERVED: same pool/slot/gen/creation seq.
        CHECK(res.cell == c);
        CHECK(ws.cells().isCurrent(c));
        CHECK(ws.cells().seqOf(c) == seq0);
        CHECK(ws.cells().stats().totalCreated == cellsCreated);
        const auto snap = ws.cells().get(c);
        CHECK(snap.nodes[3] == n[4]);
        CHECK(ws.topologyRevision() == rev0 + 1);
        CHECK(ws.state() == WorkspaceState::CommittedAdjacencyStale);
        // Net materialized growth exactly the fresh creates.
        CHECK(ws.edges().liveCount() == 6);
        CHECK(ws.faces().liveCount() == 4);
        CHECK(topo::workspaceTopologyReady(ws.nodes(), ws.cells(),
                                           ws.edges(), ws.faces()));
        ws.rebuildAdjacency();
        CHECK(ws.adjacency().boundaryFaces().size() == 4);
        // Retained base edge still current; apex edges retired.
        CHECK(ws.edges().find(ws.nodes(), n[0], n[1]).has_value());
        CHECK(!ws.edges().find(ws.nodes(), n[0], n[3]).has_value());
        CHECK(ws.edges().find(ws.nodes(), n[0], n[4]).has_value());
    }

    // [replace_updates_orientation] Pure permutation (n0,n1,n2,n3) ->
    // (n1,n0,n2,n3): every corner SET retained (0 creates, 0 removes),
    // and the retained faces' senses are recomputed and updated in place.
    {
        std::vector<NodeHandle> n;
        auto [wsp, c] = makeTetWs(n, 0);
        auto& ws = *wsp;
        auto base = ws.faces().find(
            ws.nodes(), std::vector<NodeHandle>{n[0], n[1], n[2]});
        CHECK(base.has_value());
        CHECK(ws.faces().get(*base).incident[0].sense == +1);
        auto tx = ws.beginTransaction();
        tx.replaceCellConnectivity(
            c, std::vector<NodeHandle>{n[1], n[0], n[2], n[3]});
        const auto v = tx.preview();
        CHECK(v.retainedEdges == 6 && v.newEdgeRecords == 0);
        CHECK(v.retainedFaces == 4 && v.newFaceRecords == 0);
        CHECK(v.removedEdgeIncidence == 0 && v.removedFaceIncidence == 0);
        (void)tx.commit();
        CHECK(ws.faces().get(*base).incident[0].sense == -1);   // flipped
        CHECK(topo::workspaceTopologyReady(ws.nodes(), ws.cells(),
                                           ws.edges(), ws.faces()));
        ws.rebuildAdjacency();
        CHECK(ws.adjacency().checkValid(ws.nodes(), ws.cells(), ws.edges(),
                                        ws.faces()));
    }

    // [nonmanifold_reject_policy + nonmanifold_allow_policy] Third cell
    // on a shared face: Reject throws in Phase 1; AllowAndDiagnose
    // commits REPRESENTABLE nonmanifold topology, the composite gate
    // stays true, and the face audit explicitly identifies it.
    {
        std::vector<NodeHandle> n;
        auto [wsp, c1] = makeTetWs(n, 2);
        auto& ws = *wsp;
        (void)c1;
        auto tx0 = ws.beginTransaction();
        (void)tx0.addCell(ElementType::Tet4,
                          std::vector<NodeHandle>{n[0], n[2], n[1], n[4]},
                          {3, 1});
        (void)tx0.commit();
        ws.rebuildAdjacency();
        // Reject policy (default): Phase-1 structured rejection.
        {
            auto tx = ws.beginTransaction();
            bool rejected = false;
            try {
                (void)tx.addCell(ElementType::Tet4,
                                 std::vector<NodeHandle>{n[0], n[1], n[2],
                                                         n[5]},
                                 {3, 1});
            } catch (const std::invalid_argument&) { rejected = true; }
            CHECK(rejected);
            tx.abandon();
        }
        // AllowAndDiagnose: representable, ready, and DIAGNOSED.
        {
            topo::TopologyPolicy allow;
            allow.nonManifold = topo::NonManifoldPolicy::AllowAndDiagnose;
            auto tx = ws.beginTransaction(allow);
            (void)tx.addCell(ElementType::Tet4,
                             std::vector<NodeHandle>{n[0], n[1], n[2],
                                                     n[5]},
                             {3, 1});
            auto res = tx.commit();
            CHECK(res.topologyChanged);
            CHECK(topo::workspaceTopologyReady(ws.nodes(), ws.cells(),
                                               ws.edges(), ws.faces()));
            auto shared = ws.faces().find(
                ws.nodes(), std::vector<NodeHandle>{n[0], n[1], n[2]});
            CHECK(shared.has_value());
            CHECK(ws.faces().get(*shared).incident.size() == 3);
            const auto rep = topo::auditFaceTopology(ws.nodes(), ws.cells(),
                                                     ws.edges(), ws.faces());
            bool diagnosed = false;
            for (const auto& d : rep.items)
                diagnosed |= d.code == diag::Code::NonManifoldFace;
            CHECK(diagnosed);                // representable != silently ok
            CHECK(rep.acceptable(diag::Severity::Warning));   // Info only
            ws.rebuildAdjacency();           // adjacency handles it too
            CHECK(ws.adjacency().nonManifoldFaces().size() == 1);
            CHECK(ws.adjacency().checkValid(ws.nodes(), ws.cells(),
                                            ws.edges(), ws.faces()));
        }
    }

    // [interior_orientation_reject_policy + stage-3 product-form fix
    //  regression] Same-winding second cell rejected; and after the
    //  face's CREATOR is removed (survivor sense -1), a +1 newcomer is
    //  correctly ACCEPTED (the old check assumed the survivor was +1).
    {
        std::vector<NodeHandle> n;
        auto [wsp, cA] = makeTetWs(n, 1);
        auto& ws = *wsp;
        {
            auto tx = ws.beginTransaction();
            bool rejected = false;
            try {   // SAME winding as cA across the shared face
                (void)tx.addCell(ElementType::Tet4,
                                 std::vector<NodeHandle>{n[0], n[1], n[2],
                                                         n[4]},
                                 {3, 1});
            } catch (const std::invalid_argument&) { rejected = true; }
            CHECK(rejected);
            tx.abandon();
        }
        // Regression: add opposite tet B (sense -1), remove creator A,
        // then add C with A's winding (sense +1 vs survivor -1: product
        // -1 -> LEGAL).
        auto txB = ws.beginTransaction();
        (void)txB.addCell(ElementType::Tet4,
                          std::vector<NodeHandle>{n[0], n[2], n[1], n[4]},
                          {3, 1});
        (void)txB.commit();
        ws.rebuildAdjacency();
        topo::TopologyPolicy pol; pol.permitOrphanNodes = true;
        auto txR = ws.beginTransaction(pol);
        txR.removeCell(cA);
        (void)txR.commit();
        ws.rebuildAdjacency();
        auto shared = ws.faces().find(
            ws.nodes(), std::vector<NodeHandle>{n[0], n[1], n[2]});
        CHECK(shared.has_value());
        CHECK(ws.faces().get(*shared).incident.size() == 1);
        CHECK(ws.faces().get(*shared).incident[0].sense == -1);   // survivor
        auto txC = ws.beginTransaction();
        (void)txC.addCell(ElementType::Tet4,
                          std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                          {3, 1});
        auto resC = txC.commit();            // old check would have thrown
        CHECK(resC.topologyChanged);
        const auto fd = ws.faces().get(*shared);
        CHECK(fd.incident.size() == 2);
        CHECK(fd.incident[0].sense * fd.incident[1].sense == -1);
        CHECK(topo::workspaceTopologyReady(ws.nodes(), ws.cells(),
                                           ws.edges(), ws.faces()));
    }

    // [validator_rejects_before_prepare + validator_cannot_query_workspace]
    {
        std::vector<NodeHandle> n;
        auto [wsp, c] = makeTetWs(n, 1);
        auto& ws = *wsp;
        topo::TopologyPolicy pol; pol.permitOrphanNodes = true;
        const auto sem = ws.semanticTopologyHash();
        const auto str = ws.strictWorkspaceHash();
        auto tx = ws.beginTransaction(pol);
        tx.replaceCellConnectivity(
            c, std::vector<NodeHandle>{n[0], n[1], n[2], n[4]});
        int wsQueriesRejected = 0;
        const bool ok = tx.validate([&](const TransactionPreview&) {
            try { (void)ws.nodes(); }        // preview only: workspace
            catch (const std::logic_error&) { ++wsQueriesRejected; }
            try { (void)ws.semanticTopologyHash(); }
            catch (const std::logic_error&) { ++wsQueriesRejected; }
            return false;                    // REJECT
        });
        CHECK(!ok);
        CHECK(wsQueriesRejected == 2);       // cannot query workspace
        CHECK(ws.state() == WorkspaceState::TransactionStaging);   // no prep
        bool commitRefused = false;
        try { (void)tx.commit(); }
        catch (const std::logic_error&) { commitRefused = true; }
        CHECK(commitRefused);
        tx.abandon();
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.strictWorkspaceHash() == str);
        CHECK(ws.cells().get(c).nodes[3] == n[3]);   // untouched
    }

    // [replace_failure_each_phase2_stage + failed_replace_hashes_unchanged
    //  + failed_replace_bindings_unchanged +
    //  failed_replace_consumes_no_sequences +
    //  abandoned_replace_has_no_effect +
    //  surviving_unaffected_handles_unchanged]
    {
        using Probe = topo::detail::TransactionAccessTestProbe;
        std::vector<NodeHandle> n;
        auto [wsp, c] = makeTetWs(n, 1);
        auto& ws = *wsp;
        topo::TopologyPolicy pol; pol.permitOrphanNodes = true;
        const auto sem = ws.semanticTopologyHash();
        const auto str = ws.strictWorkspaceHash();
        const auto tc = ws.cells().stats().totalCreated;
        const auto te = ws.edges().stats().totalCreated;
        const auto tf = ws.faces().stats().totalCreated;
        const auto b0 = Probe::bindings(ws.cells(), ws.edges(), ws.faces());
        auto e01 = ws.edges().find(ws.nodes(), n[0], n[1]);
        for (int stage = 1; stage <= 3; ++stage) {
            auto tx = ws.beginTransaction(pol);
            tx.replaceCellConnectivity(
                c, std::vector<NodeHandle>{n[0], n[1], n[2], n[4]});
            topo::faultCountdown() = stage;
            bool threw = false;
            try { tx.prepare(); } catch (const std::bad_alloc&) { threw = true; }
            topo::faultCountdown() = 0;
            CHECK(threw);
            tx.abandon();
            CHECK(ws.semanticTopologyHash() == sem);
            CHECK(ws.strictWorkspaceHash() == str);
        }
        // Abandon after successful prepare: no effect either.
        {
            auto tx = ws.beginTransaction(pol);
            tx.replaceCellConnectivity(
                c, std::vector<NodeHandle>{n[0], n[1], n[2], n[4]});
            tx.prepare();
            tx.abandon();
        }
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.strictWorkspaceHash() == str);
        CHECK(ws.cells().stats().totalCreated == tc);
        CHECK(ws.edges().stats().totalCreated == te);
        CHECK(ws.faces().stats().totalCreated == tf);
        const auto b1 = Probe::bindings(ws.cells(), ws.edges(), ws.faces());
        CHECK(b1.cellNode == b0.cellNode && b1.edgeNode == b0.edgeNode &&
              b1.edgeCell == b0.edgeCell && b1.faceNode == b0.faceNode &&
              b1.faceCell == b0.faceCell);
        CHECK(ws.cells().isCurrent(c));
        CHECK(ws.cells().get(c).nodes[3] == n[3]);
        CHECK(e01.has_value() && ws.edges().isCurrent(*e01));
        CHECK(topo::workspaceTopologyReady(ws.nodes(), ws.cells(),
                                           ws.edges(), ws.faces()));
    }
}

// ---------------------------------------------------------------------------
static void testTransactionMetadata() {
    TEST("topo: M4-B.6B - metadata transactions (existing classification IDs)");
    using topo::NodePool; using topo::CellPool; using topo::EdgePool;
    using topo::FacePool; using topo::TopologyWorkspace;
    using topo::WorkspaceState; using topo::ElementType;
    using topo::NodeHandle; using topo::TransactionPreview;
    using mesh::Classification;

    auto makeWs = [](std::vector<NodeHandle>& n) {
        NodePool np;
        for (int i = 0; i < 4; ++i)
            n.push_back(np.add({double(i), double(i % 2), 0}, {3, 1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4,
                        std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                        {3, 1}, np);
        EdgePool ep; ep.insertCellEdges(np, cp, c);
        FacePool fp; fp.insertCellFaces(np, cp, ep, c);
        return std::pair(std::make_unique<TopologyWorkspace>(
                             std::move(np), std::move(cp), std::move(ep),
                             std::move(fp)),
                         c);
    };

    // [node/cell/edge/face_classification_atomic +
    //  metadata_preserves_handles/connectivity/incidence +
    //  metadata_changes_hashes_deterministically +
    //  metadata_revision_invalidates_adjacency + edge constraint]
    {
        std::vector<NodeHandle> n;
        auto [wsp, c] = makeWs(n);
        auto& ws = *wsp;
        auto e01 = ws.edges().find(ws.nodes(), n[0], n[1]);
        auto base = ws.faces().find(
            ws.nodes(), std::vector<NodeHandle>{n[0], n[1], n[2]});
        CHECK(e01.has_value() && base.has_value());
        const auto sem0 = ws.semanticTopologyHash();
        const auto str0 = ws.strictWorkspaceHash();
        // Representative UNAFFECTED entities, captured before any update.
        const auto otherNodeCls = ws.nodes().get(n[1]).cls;
        const auto otherEdge = *ws.edges().find(ws.nodes(), n[1], n[2]);
        const auto otherEdgeCls = ws.edges().get(otherEdge).cls;
        const auto otherEdgeConstraint =
            ws.edges().get(otherEdge).requiredConstraint;
        const auto otherFace = *ws.faces().find(
            ws.nodes(), std::vector<NodeHandle>{n[0], n[1], n[3]});
        const auto otherFaceCls = ws.faces().get(otherFace).cls;
        const auto conn0 = ws.cells().get(c);
        const auto edges0 = std::vector<topo::EdgeHandle>(
            ws.edges().cellEdges(ws.cells(), c).begin(),
            ws.edges().cellEdges(ws.cells(), c).end());

        auto setAndVerify = [&](auto stage, auto getVal, auto expect) {
            const auto rev = ws.topologyRevision();
            auto tx = ws.beginTransaction();
            stage(tx);
            auto res = tx.commit();
            CHECK(res.topologyChanged);
            CHECK(ws.topologyRevision() == rev + 1);
            // [metadata_revision_invalidates_adjacency]
            CHECK(ws.state() == WorkspaceState::CommittedAdjacencyStale);
            CHECK(!ws.adjacencyReady());
            ws.rebuildAdjacency();
            CHECK(ws.adjacencyReady());
            CHECK(getVal() == expect);
        };
        setAndVerify([&](auto& tx) {
            tx.setNodeClassification(n[0], Classification{2, 7}); },
            [&] { return ws.nodes().get(n[0]).cls.tag; }, 7);
        setAndVerify([&](auto& tx) {
            tx.setCellClassification(c, Classification{3, 9}); },
            [&] { return ws.cells().get(c).cls.tag; }, 9);
        setAndVerify([&](auto& tx) {
            tx.setEdgeClassification(*e01, Classification{1, 5}); },
            [&] { return ws.edges().get(*e01).cls.tag; }, 5);
        setAndVerify([&](auto& tx) {
            tx.setFaceClassification(*base, Classification{2, 4}); },
            [&] { return ws.faces().get(*base).cls.tag; }, 4);
        setAndVerify([&](auto& tx) {
            tx.setEdgeRequiredConstraint(*e01, true); },
            [&] { return ws.edges().get(*e01).requiredConstraint ? 1 : 0; },
            1);
        // NO implicit propagation: only the named entities changed.
        // (Review correction: the original form here contained an
        // unconditional "|| true" - an assertion structurally incapable
        // of failing. Unaffected values are now captured BEFORE the
        // operations and compared AFTER.)
        CHECK(ws.nodes().get(n[1]).cls.dim == otherNodeCls.dim);
        CHECK(ws.nodes().get(n[1]).cls.tag == otherNodeCls.tag);
        CHECK(ws.edges().get(otherEdge).cls.dim == otherEdgeCls.dim);
        CHECK(ws.edges().get(otherEdge).cls.tag == otherEdgeCls.tag);
        CHECK(ws.edges().get(otherEdge).requiredConstraint ==
              otherEdgeConstraint);
        CHECK(ws.faces().get(otherFace).cls.dim == otherFaceCls.dim);
        CHECK(ws.faces().get(otherFace).cls.tag == otherFaceCls.tag);
        // [metadata_preserves_handles/connectivity/incidence]
        CHECK(ws.cells().isCurrent(c));
        CHECK(ws.nodes().isCurrent(n[0]) && ws.edges().isCurrent(*e01));
        CHECK(ws.faces().isCurrent(*base));
        const auto conn1 = ws.cells().get(c);
        for (std::uint8_t i = 0; i < conn0.count; ++i)
            CHECK(conn1.nodes[i] == conn0.nodes[i]);
        const auto ce = ws.edges().cellEdges(ws.cells(), c);
        CHECK(ce.size() == edges0.size());
        for (std::size_t i = 0; i < edges0.size(); ++i)
            CHECK(ce[i] == edges0[i]);
        CHECK(ws.faces().cellFaces(ws.cells(), c).size() == 4);
        // [metadata_changes_hashes_deterministically] hashes moved, and
        // restoring every value restores BOTH hashes exactly.
        CHECK(ws.semanticTopologyHash() != sem0);
        CHECK(ws.strictWorkspaceHash() != str0);
        std::vector<std::function<void(topo::TopologyTransaction&)>> undos;
        undos.push_back([&](topo::TopologyTransaction& tx) {
            tx.setNodeClassification(n[0], Classification{3, 1}); });
        undos.push_back([&](topo::TopologyTransaction& tx) {
            tx.setCellClassification(c, Classification{3, 1}); });
        undos.push_back([&](topo::TopologyTransaction& tx) {
            tx.setEdgeClassification(*e01, Classification{}); });
        undos.push_back([&](topo::TopologyTransaction& tx) {
            tx.setFaceClassification(*base, Classification{}); });
        undos.push_back([&](topo::TopologyTransaction& tx) {
            tx.setEdgeRequiredConstraint(*e01, false); });
        for (auto& undo : undos) {
            auto tx = ws.beginTransaction();
            undo(tx);
            (void)tx.commit();
            ws.rebuildAdjacency();
        }
        CHECK(ws.semanticTopologyHash() == sem0);
        CHECK(ws.strictWorkspaceHash() == str0);
    }

    // [no-op rule] new value == current value: topologyChanged=false,
    // revision unchanged, adjacency REMAINS ready (never went stale).
    {
        std::vector<NodeHandle> n;
        auto [wsp, c] = makeWs(n);
        auto& ws = *wsp;
        const auto rev = ws.topologyRevision();
        const auto sem = ws.semanticTopologyHash();
        auto tx = ws.beginTransaction();
        tx.setCellClassification(c, Classification{3, 1});   // current value
        CHECK(tx.preview().noOp);
        auto res = tx.commit();
        CHECK(!res.topologyChanged);
        CHECK(ws.topologyRevision() == rev);
        CHECK(ws.state() == WorkspaceState::Ready);
        CHECK(ws.adjacencyReady());          // never invalidated
        CHECK(ws.semanticTopologyHash() == sem);
    }

    // [metadata_validator_rejects + metadata_abandon_has_no_effect +
    //  metadata_failure_has_no_effect + preview values]
    {
        std::vector<NodeHandle> n;
        auto [wsp, c] = makeWs(n);
        auto& ws = *wsp;
        const auto sem = ws.semanticTopologyHash();
        const auto str = ws.strictWorkspaceHash();
        {
            auto tx = ws.beginTransaction();
            tx.setCellClassification(c, Classification{3, 42});
            const bool ok = tx.validate([&](const TransactionPreview& v) {
                CHECK(v.kind == TransactionPreview::Kind::Metadata);
                CHECK(v.metaTarget == TransactionPreview::MetaTarget::Cell);
                CHECK(v.cell == c);
                CHECK(v.oldClassification.tag == 1);
                CHECK(v.newClassification.tag == 42);
                CHECK(!v.noOp);
                return false;                // REJECT
            });
            CHECK(!ok);
            bool refused = false;
            try { (void)tx.commit(); }
            catch (const std::logic_error&) { refused = true; }
            CHECK(refused);
            tx.abandon();
        }
        {                                    // abandon (incl. after prepare)
            auto tx = ws.beginTransaction();
            tx.setCellClassification(c, Classification{3, 42});
            tx.prepare();
            tx.abandon();
        }
        {                                    // Phase-2 failure
            auto tx = ws.beginTransaction();
            tx.setCellClassification(c, Classification{3, 42});
            topo::faultCountdown() = 1;
            bool threw = false;
            try { tx.prepare(); } catch (const std::bad_alloc&) { threw = true; }
            topo::faultCountdown() = 0;
            CHECK(threw);
            tx.abandon();
        }
        CHECK(ws.cells().get(c).cls.tag == 1);   // untouched throughout
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.strictWorkspaceHash() == str);
        CHECK(ws.state() == WorkspaceState::Ready);
    }

    // [topo.transaction.metadata_rejects_stale_and_foreign_targets]
    // Phase 3 uses unchecked slot access, so Phase-1 target validation
    // is a critical safety precondition: every metadata family must
    // reject stale AND foreign targets, leaving revision, bindings,
    // hashes and adjacency readiness untouched, workspace readable
    // after abandon.
    {
        using Probe = topo::detail::TransactionAccessTestProbe;
        std::vector<NodeHandle> n;
        auto [wsp, c] = makeWs(n);
        auto& ws = *wsp;
        auto e01 = *ws.edges().find(ws.nodes(), n[0], n[1]);
        auto base = *ws.faces().find(
            ws.nodes(), std::vector<NodeHandle>{n[0], n[1], n[2]});
        // Stale handles: alter generation. Foreign handles: alter pool id.
        NodeHandle staleN = n[0]; staleN.gen += 1;
        NodeHandle foreignN = n[0]; foreignN.pool += 1;
        topo::CellHandle staleC = c; staleC.gen += 1;
        topo::CellHandle foreignC = c; foreignC.pool += 1;
        topo::EdgeHandle staleE = e01; staleE.gen += 1;
        topo::EdgeHandle foreignE = e01; foreignE.pool += 1;
        topo::FaceHandle staleF = base; staleF.gen += 1;
        topo::FaceHandle foreignF = base; foreignF.pool += 1;
        const auto rev = ws.topologyRevision();
        const auto sem = ws.semanticTopologyHash();
        const auto str = ws.strictWorkspaceHash();
        const auto b0 = Probe::bindings(ws.cells(), ws.edges(), ws.faces());
        const std::vector<
            std::function<void(topo::TopologyTransaction&)>> attempts = {
            [&](auto& tx) { tx.setNodeClassification(staleN, {2, 7}); },
            [&](auto& tx) { tx.setNodeClassification(foreignN, {2, 7}); },
            [&](auto& tx) { tx.setCellClassification(staleC, {3, 7}); },
            [&](auto& tx) { tx.setCellClassification(foreignC, {3, 7}); },
            [&](auto& tx) { tx.setEdgeClassification(staleE, {1, 7}); },
            [&](auto& tx) { tx.setEdgeClassification(foreignE, {1, 7}); },
            [&](auto& tx) { tx.setFaceClassification(staleF, {2, 7}); },
            [&](auto& tx) { tx.setFaceClassification(foreignF, {2, 7}); },
            [&](auto& tx) { tx.setEdgeRequiredConstraint(staleE, true); },
            [&](auto& tx) { tx.setEdgeRequiredConstraint(foreignE, true); },
        };
        int rejected = 0;
        for (const auto& attempt : attempts) {
            auto tx = ws.beginTransaction();
            try { attempt(tx); }
            catch (const std::invalid_argument&) { ++rejected; }
            tx.abandon();
            CHECK(ws.state() == WorkspaceState::Ready);   // readable
            CHECK(ws.adjacencyReady());
        }
        CHECK(rejected == 10);               // every family, both modes
        CHECK(ws.topologyRevision() == rev);
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.strictWorkspaceHash() == str);
        const auto b1 = Probe::bindings(ws.cells(), ws.edges(), ws.faces());
        CHECK(b1.cellNode == b0.cellNode && b1.edgeNode == b0.edgeNode &&
              b1.edgeCell == b0.edgeCell && b1.faceNode == b0.faceNode &&
              b1.faceCell == b0.faceCell);
        CHECK(ws.nodes().get(n[0]).cls.tag == 1);   // nothing written
        CHECK(ws.edges().get(e01).requiredConstraint == false);
    }
}

// ---------------------------------------------------------------------------
static void testTransactionNodeLifecycle() {
    TEST("topo: M4-B.6C - node lifecycle (create, then remove)");
    using topo::NodePool; using topo::CellPool; using topo::EdgePool;
    using topo::FacePool; using topo::TopologyWorkspace;
    using topo::WorkspaceState; using topo::ElementType;
    using topo::NodeHandle; using topo::TransactionPreview;
    using Probe = topo::detail::TransactionAccessTestProbe;

    auto makeWs = [](std::vector<NodeHandle>& n) {
        NodePool np;
        for (int i = 0; i < 4; ++i)
            n.push_back(np.add({double(i), double(i % 2), 0}, {3, 1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4,
                        std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                        {3, 1}, np);
        EdgePool ep; ep.insertCellEdges(np, cp, c);
        FacePool fp; fp.insertCellFaces(np, cp, ep, c);
        return std::pair(std::make_unique<TopologyWorkspace>(
                             std::move(np), std::move(cp), std::move(ep),
                             std::move(fp)),
                         c);
    };
    topo::TopologyPolicy orphanOk; orphanOk.permitOrphanNodes = true;

    // [add_node_atomic + policy definition + add_node_fresh_slot_exact +
    //  add_node_validator_rejects + add_node_abandon_has_no_effect +
    //  add_node_failure_preserves_state]
    {
        std::vector<NodeHandle> n;
        auto [wsp, c] = makeWs(n);
        (void)c;
        auto& ws = *wsp;
        // POLICY (defined clearly): default permitOrphanNodes == false
        // REJECTS creation - a new node is necessarily isolated in a
        // one-operation transaction.
        {
            auto tx = ws.beginTransaction();
            bool rejected = false;
            try { tx.addNode({9, 9, 9}, {0, 3}); }
            catch (const std::invalid_argument&) { rejected = true; }
            CHECK(rejected);
            tx.abandon();
        }
        const auto sem = ws.semanticTopologyHash();
        const auto str = ws.strictWorkspaceHash();
        const auto tn = ws.nodes().stats().totalCreated;
        const auto nodeCap = Probe::slotCap(ws.nodes());
        const auto nodeMat = Probe::materialized(ws.nodes());
        // Validator rejection: preview VALUES; commit refuses.
        {
            auto tx = ws.beginTransaction(orphanOk);
            tx.addNode({9, 9, 9}, {0, 3});
            const bool ok = tx.validate([&](const TransactionPreview& v) {
                CHECK(v.kind == TransactionPreview::Kind::AddNode);
                CHECK(v.nodeCoord.x == 9);
                CHECK(v.newClassification.tag == 3);
                return false;
            });
            CHECK(!ok);
            bool refused = false;
            try { (void)tx.commit(); }
            catch (const std::logic_error&) { refused = true; }
            CHECK(refused);
            tx.abandon();
        }
        // Failure + abandon preserve state.
        {
            auto tx = ws.beginTransaction(orphanOk);
            tx.addNode({9, 9, 9}, {0, 3});
            topo::faultCountdown() = 1;
            bool threw = false;
            try { tx.prepare(); } catch (const std::bad_alloc&) { threw = true; }
            topo::faultCountdown() = 0;
            CHECK(threw);
            tx.abandon();
        }
        {
            auto tx = ws.beginTransaction(orphanOk);
            tx.addNode({9, 9, 9}, {0, 3});
            tx.prepare();
            tx.abandon();
        }
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.strictWorkspaceHash() == str);
        CHECK(ws.nodes().stats().totalCreated == tn);
        // Atomic creation: fresh slot EXACT (materialized +1), revision
        // increments, adjacency invalidated, no implicit entities.
        const auto rev = ws.topologyRevision();
        const auto cellsLive = ws.cells().stats().live;
        const auto edgesLive = ws.edges().liveCount();
        const auto facesLive = ws.faces().liveCount();
        auto tx = ws.beginTransaction(orphanOk);
        tx.addNode({9, 9, 9}, {0, 3});
        auto res = tx.commit();
        CHECK(res.topologyChanged);
        CHECK(ws.topologyRevision() == rev + 1);
        CHECK(ws.state() == WorkspaceState::CommittedAdjacencyStale);
        CHECK(!ws.adjacencyReady());
        CHECK(ws.nodes().isCurrent(res.node));
        CHECK(ws.nodes().get(res.node).cls.tag == 3);
        CHECK(ws.nodes().get(res.node).coord.x == 9);
        CHECK(Probe::materialized(ws.nodes()) == nodeMat + 1);   // fresh
        CHECK(Probe::slotCap(ws.nodes()) >= nodeCap);
        CHECK(ws.cells().stats().live == cellsLive);   // no implicit
        CHECK(ws.edges().liveCount() == edgesLive);    // creations
        CHECK(ws.faces().liveCount() == facesLive);
        // Composite gate: an orphan node does not violate readiness.
        CHECK(topo::workspaceTopologyReady(ws.nodes(), ws.cells(),
                                           ws.edges(), ws.faces()));
        ws.rebuildAdjacency();
        CHECK(ws.adjacencyReady());

        // [remove_unreferenced_node_atomic +
        //  remove_node_reuses_slot_generation + add_node_reuses_slot]
        const auto slot0 = res.node.slot;
        const auto gen0 = res.node.gen;
        auto tx2 = ws.beginTransaction();
        tx2.removeNode(res.node);            // orphan: removable
        auto res2 = tx2.commit();
        CHECK(res2.topologyChanged);
        CHECK(!ws.nodes().isCurrent(res.node));
        ws.rebuildAdjacency();
        auto tx3 = ws.beginTransaction(orphanOk);
        tx3.addNode({8, 8, 8}, {0, 4});
        auto res3 = tx3.commit();
        CHECK(res3.node.slot == slot0);      // slot reused,
        CHECK(res3.node.gen == gen0 + 1);    // generation advanced
        CHECK(Probe::materialized(ws.nodes()) == nodeMat + 1);   // no growth
        ws.rebuildAdjacency();
    }

    // [remove_referenced_node_rejected (+ per-layer defense-in-depth) +
    //  remove_node_failure_preserves_state]
    {
        std::vector<NodeHandle> n;
        auto [wsp, c] = makeWs(n);
        (void)c;
        auto& ws = *wsp;
        const auto sem = ws.semanticTopologyHash();
        // Cell-referenced node: rejected (audit-valid states always hit
        // this check first, because valid edges/faces only reference
        // cell-connectivity nodes).
        {
            auto tx = ws.beginTransaction();
            bool rejected = false;
            try { tx.removeNode(n[0]); }
            catch (const std::invalid_argument&) { rejected = true; }
            CHECK(rejected);
            tx.abandon();
        }
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.nodes().isCurrent(n[0]));
        // Failure path on a legitimate removal.
        auto txa = ws.beginTransaction(orphanOk);
        txa.addNode({7, 7, 7}, {0, 5});
        auto orphan = txa.commit().node;
        ws.rebuildAdjacency();
        {
            auto tx = ws.beginTransaction();
            tx.removeNode(orphan);
            topo::faultCountdown() = 1;
            bool threw = false;
            try { tx.prepare(); } catch (const std::bad_alloc&) { threw = true; }
            topo::faultCountdown() = 0;
            CHECK(threw);
            tx.abandon();
            CHECK(ws.nodes().isCurrent(orphan));
        }
    }

    // [remove_node_with_edge_reference_rejected /
    //  remove_node_with_face_reference_rejected] Defense-in-depth at
    // unit level via the probe: deliberately PARTIAL pools (cell removed
    // through the frozen public ops without cascading) leave live edges/
    // faces referencing nodes no cell references - the exact
    // corrupt-adjacent states these branches guard against.
    {
        // Edge branch: cell gone, edges kept.
        NodePool np; std::vector<NodeHandle> n;
        for (int i = 0; i < 2; ++i)
            n.push_back(np.add({double(i), 0, 0}, {1, 1}));
        CellPool cp;
        auto c = cp.add(ElementType::Segment2,
                        std::vector<NodeHandle>{n[0], n[1]}, {1, 1}, np);
        EdgePool ep; ep.insertCellEdges(np, cp, c);
        ep.removeCellEdges(cp, c);           // drop incidence first...
        ep.insertCellEdges(np, cp, c);       // ...recreate cleanly
        FacePool fp;
        cp.remove(c);                        // cell gone; edge LIVE
        bool edgeRejected = false;
        try { (void)Probe::planRemoveNode(np, cp, ep, fp, n[0]); }
        catch (const std::invalid_argument&) { edgeRejected = true; }
        CHECK(edgeRejected);
    }
    {
        // Face branch: cell and edges gone, faces kept.
        NodePool np; std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i)
            n.push_back(np.add({double(i), double(i % 2), 0}, {3, 1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4,
                        std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                        {3, 1}, np);
        EdgePool ep; ep.insertCellEdges(np, cp, c);
        FacePool fp; fp.insertCellFaces(np, cp, ep, c);
        ep.removeCellEdges(cp, c);           // edges gone
        cp.remove(c);                        // cell gone; faces LIVE
        bool faceRejected = false;
        try { (void)Probe::planRemoveNode(np, cp, ep, fp, n[0]); }
        catch (const std::invalid_argument&) { faceRejected = true; }
        CHECK(faceRejected);
    }

    // [topo.transaction.remove_node_rejects_stale_and_foreign]
    // Phase 3 uses prepared unchecked access, so the Phase-1 target
    // validation is a critical precondition (closure-review
    // qualification): stale and foreign node handles must be rejected,
    // with revision, adjacency readiness, BOTH hashes and pool bindings
    // unchanged, and the workspace readable after abandon.
    {
        std::vector<NodeHandle> n;
        auto [wsp, c] = makeWs(n);
        (void)c;
        auto& ws = *wsp;
        ws.rebuildAdjacency();
        auto txa = ws.beginTransaction(orphanOk);
        txa.addNode({6, 6, 6}, {0, 2});
        auto orphan = txa.commit().node;     // a legitimately removable
        ws.rebuildAdjacency();               // node to derive stale/foreign
        NodeHandle stale = orphan; stale.gen += 1;
        NodeHandle foreign = orphan; foreign.pool += 1;
        const auto rev = ws.topologyRevision();
        const auto sem = ws.semanticTopologyHash();
        const auto str = ws.strictWorkspaceHash();
        const auto b0 = Probe::bindings(ws.cells(), ws.edges(), ws.faces());
        for (const NodeHandle& bad : {stale, foreign}) {
            auto tx = ws.beginTransaction();
            bool rejected = false;
            try { tx.removeNode(bad); }
            catch (const std::invalid_argument&) { rejected = true; }
            CHECK(rejected);
            tx.abandon();
            CHECK(ws.state() == topo::WorkspaceState::Ready);   // readable
            CHECK(ws.adjacencyReady());
            CHECK(ws.topologyRevision() == rev);
            CHECK(ws.semanticTopologyHash() == sem);
            CHECK(ws.strictWorkspaceHash() == str);
            const auto b1 = Probe::bindings(ws.cells(), ws.edges(),
                                            ws.faces());
            CHECK(b1.cellNode == b0.cellNode && b1.edgeNode == b0.edgeNode &&
                  b1.edgeCell == b0.edgeCell && b1.faceNode == b0.faceNode &&
                  b1.faceCell == b0.faceCell);
        }
        CHECK(ws.nodes().isCurrent(orphan)); // untouched throughout
    }

    // [finite-coordinate rule] Nonfinite coordinates are a Phase-1
    // structured rejection (recorded decision).
    {
        std::vector<NodeHandle> n;
        auto [wsp, c] = makeWs(n);
        (void)c;
        auto& ws = *wsp;
        const auto sem = ws.semanticTopologyHash();
        const double nan = std::nan("");
        const double inf = std::numeric_limits<double>::infinity();
        int rejected = 0;
        for (const paacs::Vec3& bad :
             {paacs::Vec3{nan, 0, 0}, paacs::Vec3{0, inf, 0},
              paacs::Vec3{0, 0, -inf}}) {
            auto tx = ws.beginTransaction(orphanOk);
            try { tx.addNode(bad, {0, 1}); }
            catch (const std::invalid_argument&) { ++rejected; }
            tx.abandon();
        }
        CHECK(rejected == 3);
        CHECK(ws.semanticTopologyHash() == sem);
        CHECK(ws.state() == topo::WorkspaceState::Ready);
    }
}

// ---------------------------------------------------------------------------
static void testHashV2() {
    TEST("topo: M4-B.7 stage H - PAACS_WORKSPACE_HASH_V2 (rank-based)");
    using topo::NodePool; using topo::CellPool; using topo::EdgePool;
    using topo::FacePool; using topo::TopologyWorkspace;
    using topo::ElementType; using topo::NodeHandle;

    // [topo.hash.v2_rank_invariance] Identical logical content AND
    // identical live order, different sequence histories (gaps from
    // removed churn in ALL FOUR pools): V2 equal, V1 and strict differ.
    // SCOPE: invariance to historical sequence gaps preserving live
    // order - NOT permutation independence.
    {
        auto buildFresh = [] {
            NodePool np; std::vector<NodeHandle> n;
            for (int i = 0; i < 4; ++i)
                np.add({double(i), double(i % 2), 0}, {3, 1});
            for (const NodeHandle& nh : np.activeBySequence())
                n.push_back(nh);
            CellPool cp;
            auto c = cp.add(ElementType::Tet4,
                            std::vector<NodeHandle>{n[0], n[1], n[2],
                                                    n[3]},
                            {3, 1}, np);
            EdgePool ep; ep.insertCellEdges(np, cp, c);
            FacePool fp; fp.insertCellFaces(np, cp, ep, c);
            return std::make_unique<TopologyWorkspace>(
                std::move(np), std::move(cp), std::move(ep),
                std::move(fp));
        };
        auto buildChurned = [] {
            NodePool np;
            std::vector<NodeHandle> j;
            for (int i = 0; i < 4; ++i)          // junk first: seqs 0..3
                j.push_back(np.add({9. + i, 9, 9}, {0, 9}));
            std::vector<NodeHandle> n;
            for (int i = 0; i < 4; ++i)          // real: seqs 4..7
                n.push_back(np.add({double(i), double(i % 2), 0},
                                   {3, 1}));
            CellPool cp;
            auto junk = cp.add(ElementType::Tet4,
                               std::vector<NodeHandle>{j[0], j[1], j[2],
                                                       j[3]},
                               {0, 9}, np);
            EdgePool ep; ep.insertCellEdges(np, cp, junk);
            FacePool fp; fp.insertCellFaces(np, cp, ep, junk);
            fp.removeCellFaces(cp, junk);        // churn every pool
            ep.removeCellEdges(cp, junk);
            cp.remove(junk);
            for (const NodeHandle& jh : j) np.remove(jh);
            auto c = cp.add(ElementType::Tet4,
                            std::vector<NodeHandle>{n[0], n[1], n[2],
                                                    n[3]},
                            {3, 1}, np);
            ep.insertCellEdges(np, cp, c);
            fp.insertCellFaces(np, cp, ep, c);
            return std::make_unique<TopologyWorkspace>(
                std::move(np), std::move(cp), std::move(ep),
                std::move(fp));
        };
        auto A = buildFresh();
        auto B = buildChurned();
        CHECK(A->semanticTopologyHashV2() == B->semanticTopologyHashV2());
        CHECK(A->semanticTopologyHashV1() != B->semanticTopologyHashV1());
        CHECK(A->strictWorkspaceHash() != B->strictWorkspaceHash());
    }

    // [topo.hash.v2_sense_normalization] Same logical mesh, but the
    // shared face is CREATED by the opposite cell (frozen public op
    // order), flipping the stored orientation and all raw senses
    // together: V2 equal (normalized to the lowest-cell-rank
    // incidence), V1 differs.
    {
        auto build = [](bool sharedFaceFromSecondCell) {
            NodePool np; std::vector<NodeHandle> n;
            for (int i = 0; i < 5; ++i)
                n.push_back(np.add({double(i), double(i % 2),
                                    double(i % 3)}, {3, 1}));
            CellPool cp;
            auto t1 = cp.add(ElementType::Tet4,
                             std::vector<NodeHandle>{n[0], n[1], n[2],
                                                     n[3]},
                             {3, 1}, np);
            auto t2 = cp.add(ElementType::Tet4,
                             std::vector<NodeHandle>{n[0], n[2], n[1],
                                                     n[4]},
                             {3, 1}, np);
            EdgePool ep;
            ep.insertCellEdges(np, cp, t1);
            ep.insertCellEdges(np, cp, t2);
            FacePool fp;
            if (sharedFaceFromSecondCell) {
                fp.insertCellFaces(np, cp, ep, t2);
                fp.insertCellFaces(np, cp, ep, t1);
            } else {
                fp.insertCellFaces(np, cp, ep, t1);
                fp.insertCellFaces(np, cp, ep, t2);
            }
            return std::make_unique<TopologyWorkspace>(
                std::move(np), std::move(cp), std::move(ep),
                std::move(fp));
        };
        auto A = build(false);
        auto B = build(true);
        // The stored orientation really flipped:
        std::vector<NodeHandle> an, bn;
        for (const NodeHandle& nh : A->nodes().activeBySequence())
            an.push_back(nh);
        for (const NodeHandle& nh : B->nodes().activeBySequence())
            bn.push_back(nh);
        auto fa = *A->faces().find(
            A->nodes(), std::vector<NodeHandle>{an[0], an[1], an[2]});
        auto fb = *B->faces().find(
            B->nodes(), std::vector<NodeHandle>{bn[0], bn[1], bn[2]});
        const auto da = A->faces().get(fa);
        const auto db = B->faces().get(fb);
        auto senseOfCellSlot = [](const auto& d, std::int32_t slot) {
            for (const auto& in : d.incident)
                if (in.cellSlot == slot) return int(in.sense);
            return 0;
        };
        CHECK(senseOfCellSlot(da, 0) == +1);     // A: t1 created it
        CHECK(senseOfCellSlot(db, 0) == -1);     // B: flipped storage
        CHECK(A->semanticTopologyHashV2() == B->semanticTopologyHashV2());
        CHECK(A->semanticTopologyHashV1() != B->semanticTopologyHashV1());
    }
}

// ---------------------------------------------------------------------------
static void testHashV2StreamContract() {
    TEST("topo: M4-B.7 stage H - V2 structured token-stream contract");
    using topo::NodePool; using topo::CellPool; using topo::EdgePool;
    using topo::FacePool; using topo::TopologyWorkspace;
    using topo::ElementType; using topo::NodeHandle;
    // The schema constants are HARDCODED here independently of the
    // implementation - this test IS the stream contract, so it must not
    // read the constants from the code under test.
    constexpr std::uint64_t SEC_NODES = 0x5EC0000000000001ULL;
    constexpr std::uint64_t SEC_CELLS = 0x5EC0000000000002ULL;
    constexpr std::uint64_t SEC_EDGES = 0x5EC0000000000003ULL;
    constexpr std::uint64_t SEC_CEINC = 0x5EC0000000000004ULL;
    constexpr std::uint64_t SEC_FACES = 0x5EC0000000000005ULL;
    constexpr std::uint64_t REC_NODE = 0x4E0DE00000000001ULL;
    constexpr std::uint64_t REC_CELL = 0x4E0DE00000000002ULL;
    constexpr std::uint64_t REC_EDGE = 0x4E0DE00000000003ULL;
    constexpr std::uint64_t REC_CELLEDGES = 0x4E0DE00000000004ULL;
    constexpr std::uint64_t REC_FACE = 0x4E0DE00000000005ULL;
    constexpr std::uint64_t NO_INC = 0xFFFFFFFFULL;
    const char* MARK = "PAACS_WORKSPACE_HASH_V2";

    // Small MIXED-element workspace: Tet4 + Tri3 sharing an edge.
    auto build = [] {
        NodePool np; std::vector<NodeHandle> n;
        for (int i = 0; i < 5; ++i)
            n.push_back(np.add({double(i), double(i % 2), double(i % 3)},
                               {3, 1}));
        CellPool cp;
        auto t1 = cp.add(ElementType::Tet4,
                         std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                         {3, 1}, np);
        auto t2 = cp.add(ElementType::Tri3,
                         std::vector<NodeHandle>{n[0], n[1], n[4]},
                         {2, 2}, np);
        EdgePool ep;
        ep.insertCellEdges(np, cp, t1);
        ep.insertCellEdges(np, cp, t2);
        FacePool fp;
        fp.insertCellFaces(np, cp, ep, t1);   // Tri3 has no faces:
        (void)t2;                             // insertCellFaces rejects 2D
        return std::make_unique<TopologyWorkspace>(
            std::move(np), std::move(cp), std::move(ep), std::move(fp));
    };
    auto ws = build();
    const auto tok = ws->testV2Tokens();

    // Independent structural parser over the emitted stream.
    std::size_t i = 0;
    auto next = [&] { CHECK(i < tok.size()); return tok[i++]; };
    for (const char* p = MARK; *p; ++p)                 // version marker
        CHECK(next() == static_cast<std::uint64_t>(*p));
    CHECK(next() == SEC_NODES);                         // section order
    const auto nNodes = next();
    CHECK(nNodes == 5);                                 // count == live
    for (std::uint64_t k = 0; k < nNodes; ++k) {
        CHECK(next() == REC_NODE);
        i += 5;                                         // x y z dim tag
        CHECK(i <= tok.size());
    }
    CHECK(next() == SEC_CELLS);
    const auto nCells = next();
    CHECK(nCells == 2);
    std::vector<std::uint64_t> connCounts;
    for (std::uint64_t k = 0; k < nCells; ++k) {
        CHECK(next() == REC_CELL);
        (void)next();                                   // element_type
        const auto cc = next();                         // EXPLICIT length
        connCounts.push_back(cc);
        for (std::uint64_t j = 0; j < cc; ++j) {
            const auto rank = next();
            CHECK(rank < nNodes);                       // ranks in range
        }
        i += 2;                                         // dim tag
    }
    CHECK(connCounts[0] == 4 && connCounts[1] == 3);    // Tet4 + Tri3
    CHECK(next() == SEC_EDGES);
    const auto nEdges = next();
    CHECK(nEdges == 8);                                 // 6 tet + 2 new tri
    std::uint64_t prevLo = 0, prevHi = 0;
    for (std::uint64_t k = 0; k < nEdges; ++k) {
        CHECK(next() == REC_EDGE);
        const auto lo = next(); const auto hi = next();
        CHECK(lo < hi);                                 // canonical
        if (k > 0)                                      // sorted order
            CHECK(prevLo < lo || (prevLo == lo && prevHi < hi));
        prevLo = lo; prevHi = hi;
        i += 3;                                         // dim tag constraint
    }
    CHECK(next() == SEC_CEINC);
    CHECK(next() == nCells);
    for (std::uint64_t k = 0; k < nCells; ++k) {
        CHECK(next() == REC_CELLEDGES);
        (void)next();                                   // cell_rank
        const auto ec = next();
        CHECK(ec != NO_INC);                            // both complete
        CHECK(ec == (k == 0 ? 6u : 3u));                // tet 6, tri 3
        i += 2 * ec;                                    // rank pairs
        CHECK(i <= tok.size());
    }
    CHECK(next() == SEC_FACES);
    const auto nFaces = next();
    CHECK(nFaces == 4);                                 // tet only
    for (std::uint64_t k = 0; k < nFaces; ++k) {
        CHECK(next() == REC_FACE);
        const auto corners = next();
        CHECK(corners == 3);                            // tri faces
        std::uint64_t prev = 0;
        for (std::uint64_t j = 0; j < corners; ++j) {
            const auto r = next();
            if (j > 0) CHECK(prev < r);                 // sorted key
            prev = r;
        }
        i += 2;                                         // dim tag
        const auto inc = next();
        CHECK(inc == 1);                                // boundary faces
        (void)next();                                   // cell_rank
        const auto sense = next();
        CHECK(sense == 1);                              // normalized +1
    }
    CHECK(i == tok.size());                             // no trailing data

    // Determinism across repeated builds.
    auto ws2 = build();
    CHECK(ws2->testV2Tokens() == tok);

    // Frozen golden for this mixed workspace: token count + stream FNV.
    std::uint64_t g = 1469598103934665603ULL;
    for (const auto v : tok)
        for (int b = 0; b < 8; ++b) {
            g ^= (v >> (8 * b)) & 0xFF;
            g *= 1099511628211ULL;
        }
    std::printf("      [v2 stream golden] tokens=%zu fnv=0x%016llx\n",
                tok.size(),
                static_cast<unsigned long long>(g));
    CHECK(tok.size() == 192);            // FROZEN GOLDEN
    CHECK(g == 0x9002078b362f95deULL);   // FROZEN GOLDEN
}

// ---------------------------------------------------------------------------
// --- M4-B.7 C1 shared fixtures ---------------------------------------------
namespace c1 {
using topo::NodePool; using topo::CellPool; using topo::EdgePool;
using topo::FacePool; using topo::TopologyWorkspace;
using topo::ElementType; using topo::NodeHandle; using topo::CellHandle;
using Probe = topo::detail::TransactionAccessTestProbe;

// Churned source: gaps in ALL FOUR pools, non-default metadata
// (node/edge/constraint/face), two tets sharing a face.
inline std::unique_ptr<TopologyWorkspace> buildChurned(
    std::vector<NodeHandle>* liveOut = nullptr) {
    NodePool np;
    std::vector<NodeHandle> j;
    for (int i = 0; i < 4; ++i) j.push_back(np.add({9. + i, 9, 9}, {0, 9}));
    std::vector<NodeHandle> n;
    for (int i = 0; i < 5; ++i)
        n.push_back(np.add({double(i), double(i % 2), double(i % 3)},
                           {3, 1}));
    CellPool cp;
    auto junk = cp.add(ElementType::Tet4,
                       std::vector<NodeHandle>{j[0], j[1], j[2], j[3]},
                       {0, 9}, np);
    EdgePool ep; ep.insertCellEdges(np, cp, junk);
    FacePool fp; fp.insertCellFaces(np, cp, ep, junk);
    fp.removeCellFaces(cp, junk);
    ep.removeCellEdges(cp, junk);
    cp.remove(junk);
    for (const NodeHandle& jh : j) np.remove(jh);
    auto t1 = cp.add(ElementType::Tet4,
                     std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                     {3, 5}, np);
    auto t2 = cp.add(ElementType::Tet4,
                     std::vector<NodeHandle>{n[0], n[2], n[1], n[4]},
                     {3, 6}, np);
    ep.insertCellEdges(np, cp, t1);
    ep.insertCellEdges(np, cp, t2);
    fp.insertCellFaces(np, cp, ep, t1);
    fp.insertCellFaces(np, cp, ep, t2);
    np.setClassification(n[0], {2, 21});
    auto e01 = *ep.find(np, n[0], n[1]);
    ep.setClassification(e01, {1, 31});
    ep.setRequiredConstraint(e01, true);
    auto shared = *fp.find(np, std::vector<NodeHandle>{n[0], n[1], n[2]});
    fp.setClassification(shared, {2, 41});
    if (liveOut) *liveOut = n;
    return std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), std::move(fp));
}
// Churn INCLUDING a retired slot (reviewer: churn must include retired
// slots, not just free ones). NodePool::testRetire is the frozen public
// testing seam; HandleAllocator::free implements retirement identically
// for every pool, so this exercises the shared retirement path.
inline std::unique_ptr<TopologyWorkspace> buildChurnedWithRetired(
    std::vector<NodeHandle>* liveOut, NodeHandle* retiredOut,
    NodeHandle* freedOut) {
    NodePool np;
    // A node that will be RETIRED (generation exhausted).
    const NodeHandle doomed = np.add({7, 7, 7}, {0, 8});
    // A node that will be FREED normally (slot reusable).
    const NodeHandle freed = np.add({8, 8, 8}, {0, 9});
    std::vector<NodeHandle> n;
    for (int i = 0; i < 5; ++i)
        n.push_back(np.add({double(i), double(i % 2), double(i % 3)},
                           {3, 1}));
    np.testRetire(doomed);               // retired in place
    np.remove(freed);                    // ordinary free slot
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
    if (liveOut) *liveOut = n;
    if (retiredOut) {
        *retiredOut = doomed;
        retiredOut->gen = std::numeric_limits<std::uint32_t>::max();
    }
    if (freedOut) *freedOut = freed;
    return std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), std::move(fp));
}
// Three tets sharing one face -> a 3-incidence NONMANIFOLD face.
inline std::unique_ptr<TopologyWorkspace> buildNonmanifold() {
    NodePool np; std::vector<NodeHandle> n;
    for (int i = 0; i < 6; ++i)
        n.push_back(np.add({double(i), double(i % 2), double(i % 3)},
                           {3, 1}));
    CellPool cp;
    auto a = cp.add(ElementType::Tet4,
                    std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                    {3, 1}, np);
    auto b = cp.add(ElementType::Tet4,
                    std::vector<NodeHandle>{n[0], n[2], n[1], n[4]},
                    {3, 1}, np);
    auto c = cp.add(ElementType::Tet4,
                    std::vector<NodeHandle>{n[0], n[1], n[2], n[5]},
                    {3, 1}, np);
    EdgePool ep;
    for (auto h : {a, b, c}) ep.insertCellEdges(np, cp, h);
    FacePool fp;
    for (auto h : {a, b, c}) fp.insertCellFaces(np, cp, ep, h);
    return std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), std::move(fp));
}

// Compare a face's incidence senses across the identity domains
// RELATIVE TO AN ANCHOR (the incidence whose source cell has the lowest
// creation sequence). Works for any incidence count, including
// nonmanifold faces.
template <class Res>
inline void checkSenseRelations(const TopologyWorkspace& src,
                                const TopologyWorkspace& dst,
                                const Res& res, const topo::FaceHandle& sfh,
                                const topo::FaceHandle& dfh) {
    const auto s = src.faces().get(sfh);
    const auto d = dst.faces().get(dfh);
    CHECK(s.incident.size() == d.incident.size());
    // Build (sourceCellSeq -> sense) and (mappedCellIdent -> sense).
    std::vector<std::pair<std::uint64_t, int>> ss;   // seq, sense
    std::vector<std::pair<std::uint64_t, int>> dd;
    for (const auto& si : s.incident) {
        CellHandle sc{};
        bool found = false;
        for (const CellHandle& ch : src.cells().activeBySequence())
            if (ch.slot == si.cellSlot && ch.gen == si.cellGen) {
                sc = ch; found = true; break;
            }
        CHECK(found);
        const auto dc = *res.cellMap.find(sc);
        ss.emplace_back(src.cells().seqOf(sc), int(si.sense));
        bool matched = false;
        for (const auto& di : d.incident)
            if (di.cellSlot == dc.slot && di.cellGen == dc.gen) {
                dd.emplace_back(src.cells().seqOf(sc), int(di.sense));
                matched = true;
                break;
            }
        CHECK(matched);                 // same incident cell set
    }
    std::sort(ss.begin(), ss.end());
    std::sort(dd.begin(), dd.end());
    CHECK(ss.size() == dd.size());
    if (ss.empty()) return;
    const int sAnchor = ss.front().second;
    const int dAnchor = dd.front().second;
    CHECK(sAnchor != 0 && dAnchor != 0);
    for (std::size_t i = 0; i < ss.size(); ++i) {
        CHECK(ss[i].first == dd[i].first);              // same cell
        CHECK(ss[i].second * sAnchor == dd[i].second * dAnchor);
    }
}

inline std::unique_ptr<TopologyWorkspace> buildCanonical() {
    NodePool np; std::vector<NodeHandle> n;
    for (int i = 0; i < 4; ++i)
        n.push_back(np.add({double(i), double(i % 2), 0}, {3, 1}));
    CellPool cp;
    auto c = cp.add(ElementType::Tet4,
                    std::vector<NodeHandle>{n[0], n[1], n[2], n[3]},
                    {3, 1}, np);
    EdgePool ep; ep.insertCellEdges(np, cp, c);
    FacePool fp; fp.insertCellFaces(np, cp, ep, c);
    return std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), std::move(fp));
}

// FULL observable-state digest per M4B7_DESIGN.md (C1 review: the
// previous snapshot omitted adjacency built-from revision, per-pool
// free/materialized/retired counts, memory reports, iteration orders,
// currentness, classifications and incidence signatures).
// observableDigest lives in the shared header so the standalone
// probes report the IDENTICAL digest (C2 review).
using paacs_tests::observableDigest;
using paacs_tests::observableLogicalDigest;
}  // namespace c1

// ---------------------------------------------------------------------------
static void testCompactionSemanticHashUnchanged() {
    TEST("topo.compaction.semantic_hash_unchanged");
    auto src = c1::buildChurned();
    const auto sem = src->semanticTopologyHash();
    auto res = topo::compactWorkspace(*src);
    CHECK(res.workspace->semanticTopologyHash() == sem);
    // and the full V2 token streams agree, not merely the fold
    CHECK(res.workspace->testV2Tokens() == src->testV2Tokens());
}

static void testCompactionNodeMapTotal() {
    TEST("topo.compaction.node_map_total");
    auto src = c1::buildChurned();
    auto res = topo::compactWorkspace(*src);
    CHECK(res.nodeMap.size() == src->nodes().stats().live);
    std::size_t seen = 0;
    for (const topo::NodeHandle& h : src->nodes().activeBySequence()) {
        CHECK(res.nodeMap.find(h).has_value());
        ++seen;
    }
    CHECK(seen == res.nodeMap.size());   // exactly once each
}

static void testCompactionCellMapTotal() {
    TEST("topo.compaction.cell_map_total");
    auto src = c1::buildChurned();
    auto res = topo::compactWorkspace(*src);
    CHECK(res.cellMap.size() == src->cells().stats().live);
    std::size_t seen = 0;
    for (const topo::CellHandle& h : src->cells().activeBySequence()) {
        CHECK(res.cellMap.find(h).has_value());
        ++seen;
    }
    CHECK(seen == res.cellMap.size());
}

static void testCompactionUnmappedHandleClasses() {
    TEST("topo.compaction.unmapped_handle_classes");
    // FREE-slot handles (a genuinely freed slot, per C1 review - the
    // earlier claim only built a modified-generation handle), stale
    // generations and foreign pools find NO mapping. The named
    // free/retired/stale proofs of the C2 list follow there; retirement
    // needs generation exhaustion and is deliberately left to C2.
    std::vector<topo::NodeHandle> live;
    auto src = c1::buildChurned(&live);
    src->rebuildAdjacency();
    topo::TopologyPolicy orphanOk; orphanOk.permitOrphanNodes = true;
    auto txa = src->beginTransaction(orphanOk);
    txa.addNode({4, 4, 4}, {0, 7});
    const auto tmp = txa.commit().node;      // live...
    src->rebuildAdjacency();
    auto txr = src->beginTransaction();
    txr.removeNode(tmp);                     // ...now a FREED slot
    (void)txr.commit();
    src->rebuildAdjacency();
    CHECK(!src->nodes().isCurrent(tmp));
    CHECK(src->nodes().stats().reusableFree >= 1);
    auto res = topo::compactWorkspace(*src);
    CHECK(!res.nodeMap.find(tmp).has_value());      // FREED slot
    topo::NodeHandle stale = live[0]; stale.gen += 1;
    CHECK(!res.nodeMap.find(stale).has_value());    // stale generation
    topo::NodeHandle foreign = live[0]; foreign.pool += 1;
    CHECK(!res.nodeMap.find(foreign).has_value());  // foreign pool
    CHECK(res.nodeMap.find(live[0]).has_value());   // control: live maps
}

static void testCompactionSourceOldHandlesRemainCurrent() {
    TEST("topo.compaction.source_old_handles_remain_current");
    auto src = c1::buildChurned();
    auto res = topo::compactWorkspace(*src);
    for (const topo::NodeHandle& h : src->nodes().activeBySequence())
        CHECK(src->nodes().isCurrent(h));
    for (const topo::CellHandle& h : src->cells().activeBySequence())
        CHECK(src->cells().isCurrent(h));
    for (const topo::EdgeHandle& h : src->edges().activeBySequence())
        CHECK(src->edges().isCurrent(h));
    for (const topo::FaceHandle& h : src->faces().activeBySequence())
        CHECK(src->faces().isCurrent(h));
    (void)res;
}

static void testCompactionOldHandlesForeign() {
    TEST("topo.compaction.old_handles_foreign_to_compacted_workspace");
    auto src = c1::buildChurned();
    auto res = topo::compactWorkspace(*src);
    auto& dst = *res.workspace;
    // Pool identity is the mechanism: the destination pools carry new
    // pool ids, so source handles are foreign to them.
    CHECK(dst.nodes().poolId() != src->nodes().poolId());
    CHECK(dst.cells().poolId() != src->cells().poolId());
    for (const topo::NodeHandle& h : src->nodes().activeBySequence())
        CHECK(!dst.nodes().isCurrent(h));
    for (const topo::CellHandle& h : src->cells().activeBySequence())
        CHECK(!dst.cells().isCurrent(h));
}

static void testCompactionMappedHandlesCurrent() {
    TEST("topo.compaction.mapped_handles_current");
    auto src = c1::buildChurned();
    auto res = topo::compactWorkspace(*src);
    auto& dst = *res.workspace;
    for (const topo::NodeHandle& h : src->nodes().activeBySequence())
        CHECK(dst.nodes().isCurrent(*res.nodeMap.find(h)));
    for (const topo::CellHandle& h : src->cells().activeBySequence())
        CHECK(dst.cells().isCurrent(*res.cellMap.find(h)));
}

static void testCompactionConnectivityRemapped() {
    TEST("topo.compaction.connectivity_remapped");
    auto src = c1::buildChurned();
    auto res = topo::compactWorkspace(*src);
    auto& dst = *res.workspace;
    for (const topo::CellHandle& h : src->cells().activeBySequence()) {
        const auto so = src->cells().get(h);
        const auto sd = dst.cells().get(*res.cellMap.find(h));
        CHECK(so.type == sd.type);
        CHECK(so.count == sd.count);
        for (std::uint8_t i = 0; i < so.count; ++i)
            CHECK(sd.nodes[i] == *res.nodeMap.find(so.nodes[i]));
        CHECK(so.cls.dim == sd.cls.dim && so.cls.tag == sd.cls.tag);
    }
}

static void testCompactionMetadataReconstructed() {
    TEST("topo.compaction.metadata_reconstructed");
    std::vector<topo::NodeHandle> n;
    auto src = c1::buildChurned(&n);
    auto res = topo::compactWorkspace(*src);
    auto& dst = *res.workspace;
    const auto a = *res.nodeMap.find(n[0]);
    const auto b = *res.nodeMap.find(n[1]);
    const auto c = *res.nodeMap.find(n[2]);
    CHECK(dst.nodes().get(a).cls.tag == 21);
    const auto de = *dst.edges().find(dst.nodes(), a, b);
    CHECK(dst.edges().get(de).cls.tag == 31);
    CHECK(dst.edges().get(de).requiredConstraint == true);
    const auto df = *dst.faces().find(
        dst.nodes(), std::vector<topo::NodeHandle>{a, b, c});
    CHECK(dst.faces().get(df).cls.tag == 41);
}

static void testCompactionWorkspaceReady() {
    TEST("topo.compaction.workspace_ready");
    auto src = c1::buildChurned();
    auto res = topo::compactWorkspace(*src);
    auto& dst = *res.workspace;
    CHECK(topo::workspaceTopologyReady(dst.nodes(), dst.cells(),
                                       dst.edges(), dst.faces()));
    CHECK(dst.state() == topo::WorkspaceState::Ready);
    CHECK(topo::workspaceIncidenceComplete(dst.nodes(), dst.cells(),
                                           dst.edges(), dst.faces()));
}

static void testCompactionAdjacencyRebuilt() {
    TEST("topo.compaction.adjacency_rebuilt");
    auto src = c1::buildChurned();
    auto res = topo::compactWorkspace(*src);
    auto& dst = *res.workspace;
    // The workspace constructor builds adjacency from the new revision
    // (C1 review: no redundant second build in compactWorkspace).
    CHECK(dst.adjacencyReady());
    CHECK(dst.adjacencyBuiltFromRevision() == dst.topologyRevision());
    CHECK(dst.adjacency().checkValid(dst.nodes(), dst.cells(),
                                     dst.edges(), dst.faces()));
    CHECK(dst.adjacency().boundaryFaces().size() ==
          src->adjacency().boundaryFaces().size());
}

static void testCompactionStrictHashDeterministicRebuild() {
    TEST("topo.compaction.strict_hash_matches_deterministic_rebuild");
    // Churned source: strict CHANGES.
    {
        auto src = c1::buildChurned();
        auto res = topo::compactWorkspace(*src);
        CHECK(res.workspace->semanticTopologyHash() ==
              src->semanticTopologyHash());
        CHECK(res.workspace->strictWorkspaceHash() !=
              src->strictWorkspaceHash());
    }
    // Already-canonical source: strict EQUALITY is legitimate; pool
    // ownership - not hash inequality - proves the new domain.
    {
        auto src = c1::buildCanonical();
        auto res = topo::compactWorkspace(*src);
        CHECK(res.workspace->semanticTopologyHash() ==
              src->semanticTopologyHash());
        CHECK(res.workspace->strictWorkspaceHash() ==
              src->strictWorkspaceHash());
        for (const topo::CellHandle& h : src->cells().activeBySequence()) {
            CHECK(!res.workspace->cells().isCurrent(h));
            CHECK(res.workspace->cells().isCurrent(*res.cellMap.find(h)));
        }
    }
    // The contract is equality with the deterministic rebuild: compact
    // twice, strict hashes agree in BOTH cases.
    {
        auto src = c1::buildChurned();
        auto r1 = topo::compactWorkspace(*src);
        auto r2 = topo::compactWorkspace(*src);
        CHECK(r1.workspace->strictWorkspaceHash() ==
              r2.workspace->strictWorkspaceHash());
    }
}

static void testCompactionDeterministicConstruction() {
    TEST("topo.compaction.deterministic_construction");
    auto src = c1::buildChurned();
    auto r1 = topo::compactWorkspace(*src);
    auto r2 = topo::compactWorkspace(*src);
    CHECK(r1.workspace->testV2Tokens() == r2.workspace->testV2Tokens());
    CHECK(r1.workspace->strictWorkspaceHash() ==
          r2.workspace->strictWorkspaceHash());
    // Structural determinism holds MODULO runtime ownership: each
    // compaction mints fresh pool ids by design, so the two results are
    // deliberately different OWNERSHIP domains while being identical in
    // semantic content and strict physical layout.
    CHECK(c1::observableDigest(*r1.workspace, false) ==
          c1::observableDigest(*r2.workspace, false));
    CHECK(r1.workspace->nodes().poolId() != r2.workspace->nodes().poolId());
    // Maps agree slot-for-slot and generation-for-generation.
    for (const topo::NodeHandle& h : src->nodes().activeBySequence()) {
        const auto a = *r1.nodeMap.find(h);
        const auto b = *r2.nodeMap.find(h);
        CHECK(a.slot == b.slot && a.gen == b.gen);
    }
    for (const topo::CellHandle& h : src->cells().activeBySequence()) {
        const auto a = *r1.cellMap.find(h);
        const auto b = *r2.cellMap.find(h);
        CHECK(a.slot == b.slot && a.gen == b.gen);
    }
}

static void testCompactionMemoryReportClasses() {
    TEST("topo.compaction.memory_report_classes");
    // C1 review blocker: the byte classes must PARTITION the total and
    // the estimates must be non-zero where table/incidence memory
    // genuinely exists (the old subtraction made them identically 0).
    auto src = c1::buildChurned();
    auto res = topo::compactWorkspace(*src);
    for (const topo::MemoryReport& m : {res.before, res.after}) {
        CHECK(m.consistent());
        CHECK(m.totalBytes ==
              m.recordBytes + m.allocatorBytes + m.estimatedOtherBytes);
        CHECK(m.recordBytes > 0);
        CHECK(m.allocatorBytes > 0);
        CHECK(m.estimatedOtherBytes > 0);   // edge/face tables+incidence
        CHECK(m.totalBytes ==
              m.nodeBytes + m.cellBytes + m.edgeBytes + m.faceBytes);
    }
    // Per-pool: the same partition holds, and the pools with hash tables
    // and incidence report non-zero estimates while the plain pools
    // report none.
    for (const auto& st : {src->nodes().stats(), src->cells().stats(),
                           src->edges().stats(), src->faces().stats()}) {
        CHECK(st.bytesConsistent());
        CHECK(st.payloadBytes == st.recordBytes + st.estimatedOtherBytes);
    }
    CHECK(src->nodes().stats().estimatedOtherBytes == 0);
    CHECK(src->cells().stats().estimatedOtherBytes == 0);
    CHECK(src->edges().stats().estimatedOtherBytes > 0);
    CHECK(src->faces().stats().estimatedOtherBytes > 0);
    // Churned source retains more than the rebuild.
    CHECK(res.after.totalBytes <= res.before.totalBytes);
}

static void testCompactionSourceSuccessObservablyUnchanged() {
    TEST("topo.compaction.source_success_observably_unchanged");
    auto src = c1::buildChurned();
    src->rebuildAdjacency();
    const auto d0 = c1::observableDigest(*src);
    auto res = topo::compactWorkspace(*src);
    CHECK(res.workspace != nullptr);
    CHECK(c1::observableDigest(*src) == d0);
}

static void testCompactionSourceFailureObservablyUnchanged() {
    TEST("topo.compaction.source_failure_observably_unchanged");
    auto src = c1::buildChurned();
    src->rebuildAdjacency();
    const auto d0 = c1::observableDigest(*src);
    int injected = 0;
    for (int stage = 1; stage <= 7; ++stage) {
        topo::detail::compactionFaultCountdown() = stage;
        bool threw = false;
        try { (void)topo::compactWorkspace(*src); }
        catch (const std::bad_alloc&) { threw = true; }
        topo::detail::compactionFaultCountdown() = 0;
        if (threw) ++injected;
        CHECK(c1::observableDigest(*src) == d0);   // FULL digest
    }
    CHECK(injected == 7);
    auto res = topo::compactWorkspace(*src);       // still works after
    CHECK(res.workspace->semanticTopologyHash() ==
          src->semanticTopologyHash());
    CHECK(c1::observableDigest(*src) == d0);
}

// --- M4-B.7 C2 -------------------------------------------------------------
static void testCompactionEdgeMapTotal() {
    TEST("topo.compaction.edge_map_total");
    auto src = c1::buildChurned();
    auto res = topo::compactWorkspace(*src);
    CHECK(res.edgeMap.size() == src->edges().stats().live);
    auto& dst = *res.workspace;
    std::size_t seen = 0;
    for (const topo::EdgeHandle& h : src->edges().activeBySequence()) {
        const auto m = res.edgeMap.find(h);
        CHECK(m.has_value());
        CHECK(dst.edges().isCurrent(*m));      // mapped -> current
        CHECK(!dst.edges().isCurrent(h));      // old -> foreign
        CHECK(src->edges().isCurrent(h));      // source untouched
        // Mapped through CANONICAL KEYS, never iteration order:
        const auto so = src->edges().get(h);
        const auto sd = dst.edges().get(*m);
        CHECK(sd.a == *res.nodeMap.find(so.a));
        CHECK(sd.b == *res.nodeMap.find(so.b));
        ++seen;
    }
    CHECK(seen == res.edgeMap.size());
}

static void testCompactionFaceMapTotal() {
    TEST("topo.compaction.face_map_total");
    auto src = c1::buildChurned();
    auto res = topo::compactWorkspace(*src);
    CHECK(res.faceMap.size() == src->faces().stats().live);
    auto& dst = *res.workspace;
    std::size_t seen = 0;
    for (const topo::FaceHandle& h : src->faces().activeBySequence()) {
        const auto m = res.faceMap.find(h);
        CHECK(m.has_value());
        CHECK(dst.faces().isCurrent(*m));
        CHECK(!dst.faces().isCurrent(h));
        CHECK(src->faces().isCurrent(h));
        const auto so = src->faces().get(h);
        const auto sd = dst.faces().get(*m);
        CHECK(so.cornerCount == sd.cornerCount);
        // The destination corners are exactly the MAPPED source corners
        // as a SET (stored order may differ when a different surviving
        // cell recreates the face). C2 review: the previous form sorted
        // two vectors and then compared only their sizes - vacuous,
        // since cornerCount equality already guaranteed it.
        std::vector<topo::NodeHandle> want, got;
        for (std::uint8_t i = 0; i < so.cornerCount; ++i) {
            want.push_back(*res.nodeMap.find(so.corners[i]));
            got.push_back(sd.corners[i]);
        }
        auto byIdent = [](const topo::NodeHandle& x,
                          const topo::NodeHandle& y) {
            return std::tie(x.slot, x.gen) < std::tie(y.slot, y.gen);
        };
        std::sort(want.begin(), want.end(), byIdent);
        std::sort(got.begin(), got.end(), byIdent);
        CHECK(want == got);              // handle-for-handle
        ++seen;
    }
    CHECK(seen == res.faceMap.size());
}

static void testCompactionFreeHandlesUnmapped() {
    TEST("topo.compaction.free_handles_unmapped");
    std::vector<topo::NodeHandle> live;
    topo::NodeHandle retired{}, freed{};
    auto src = c1::buildChurnedWithRetired(&live, &retired, &freed);
    CHECK(src->nodes().stats().reusableFree >= 1);
    auto res = topo::compactWorkspace(*src);
    CHECK(!res.nodeMap.find(freed).has_value());   // FREE slot: no map
    CHECK(res.nodeMap.find(live[0]).has_value());  // control
}

static void testCompactionRetiredHandlesUnmapped() {
    TEST("topo.compaction.retired_handles_unmapped");
    std::vector<topo::NodeHandle> live;
    topo::NodeHandle retired{}, freed{};
    auto src = c1::buildChurnedWithRetired(&live, &retired, &freed);
    CHECK(src->nodes().stats().retired == 1);      // genuinely retired
    auto res = topo::compactWorkspace(*src);
    CHECK(!res.nodeMap.find(retired).has_value()); // RETIRED: no map
    // Retired slots are not carried into the new domain at all.
    CHECK(res.workspace->nodes().stats().retired == 0);
    CHECK(res.workspace->nodes().stats().live == live.size());
    CHECK(res.workspace->nodes().stats().materializedSlots ==
          live.size());                            // dense
}

static void testCompactionStaleGenerationsUnmapped() {
    TEST("topo.compaction.stale_generations_unmapped");
    std::vector<topo::NodeHandle> live;
    auto src = c1::buildChurned(&live);
    auto res = topo::compactWorkspace(*src);
    for (const topo::NodeHandle& h : src->nodes().activeBySequence()) {
        topo::NodeHandle stale = h;
        stale.gen += 1;                            // surviving slot,
        CHECK(!res.nodeMap.find(stale).has_value());  // stale generation
        topo::NodeHandle foreign = h;
        foreign.pool += 1;
        CHECK(!res.nodeMap.find(foreign).has_value());
        CHECK(res.nodeMap.find(h).has_value());    // control
    }
}

static void testCompactionIncidenceRemapped() {
    TEST("topo.compaction.incidence_remapped");
    auto src = c1::buildChurned();
    auto res = topo::compactWorkspace(*src);
    auto& dst = *res.workspace;
    // Cell-edge incidence maps elementwise through the edge map.
    for (const topo::CellHandle& ch : src->cells().activeBySequence()) {
        const auto dch = *res.cellMap.find(ch);
        const auto se = src->edges().cellEdges(src->cells(), ch);
        const auto de = dst.edges().cellEdges(dst.cells(), dch);
        CHECK(se.size() == de.size());
        for (std::size_t i = 0; i < se.size(); ++i)
            CHECK(de[i] == *res.edgeMap.find(se[i]));   // traits order
        const auto sf = src->faces().cellFaces(src->cells(), ch);
        const auto df = dst.faces().cellFaces(dst.cells(), dch);
        CHECK(sf.size() == df.size());
        for (std::size_t i = 0; i < sf.size(); ++i)
            CHECK(df[i] == *res.faceMap.find(sf[i]));
    }
    // Face incidence maps through the cell map, with SEMANTIC
    // orientation relations preserved (raw stored senses may flip
    // globally per face - identity-domain data, normalized by V2).
    for (const topo::FaceHandle& fh : src->faces().activeBySequence()) {
        const auto dfh = *res.faceMap.find(fh);
        const auto s = src->faces().get(fh);
        const auto d = dst.faces().get(dfh);
        CHECK(s.incident.size() == d.incident.size());
        // Map source incidence (cellSlot -> sense) into destination
        // cell handles, then compare products pairwise.
        for (const auto& si : s.incident) {
            // find the source cell handle for this incidence
            topo::CellHandle sc{};
            bool found = false;
            for (const topo::CellHandle& c :
                 src->cells().activeBySequence())
                if (c.slot == si.cellSlot && c.gen == si.cellGen) {
                    sc = c; found = true; break;
                }
            CHECK(found);
            const auto dc = *res.cellMap.find(sc);
            bool matched = false;
            for (const auto& di : d.incident)
                if (di.cellSlot == dc.slot && di.cellGen == dc.gen)
                    matched = true;
            CHECK(matched);                 // same incident cell set
        }
        // Every sense compared RELATIVE TO AN ANCHOR incidence (the
        // lowest source cell sequence), which is exactly the V2
        // normalization: a global per-face flip is identity-domain
        // data, the relations are the semantics. Covers 2- AND
        // 3-incidence (nonmanifold) faces alike.
        c1::checkSenseRelations(*src, dst, res, fh, dfh);
    }
    // The same check on a NONMANIFOLD fixture (three cells on one face).
    {
        auto nm = c1::buildNonmanifold();
        auto r = topo::compactWorkspace(*nm);
        std::size_t threes = 0;
        for (const topo::FaceHandle& fh : nm->faces().activeBySequence()) {
            const auto dfh = *r.faceMap.find(fh);
            if (nm->faces().get(fh).incident.size() >= 3) ++threes;
            c1::checkSenseRelations(*nm, *r.workspace, r, fh, dfh);
        }
        CHECK(threes == 1);                 // the fixture really is NM
        CHECK(r.workspace->semanticTopologyHash() ==
              nm->semanticTopologyHash());
    }
}

static void testCompactionMetadataPreserved() {
    TEST("topo.compaction.metadata_preserved");
    std::vector<topo::NodeHandle> n;
    auto src = c1::buildChurned(&n);
    auto res = topo::compactWorkspace(*src);
    auto& dst = *res.workspace;
    // Every live entity's classification survives, checked THROUGH the
    // maps rather than by re-deriving keys.
    for (const topo::NodeHandle& h : src->nodes().activeBySequence()) {
        const auto s = src->nodes().get(h);
        const auto d = dst.nodes().get(*res.nodeMap.find(h));
        CHECK(s.cls.dim == d.cls.dim && s.cls.tag == d.cls.tag);
        CHECK(s.coord.x == d.coord.x && s.coord.y == d.coord.y &&
              s.coord.z == d.coord.z);
    }
    for (const topo::CellHandle& h : src->cells().activeBySequence()) {
        const auto s = src->cells().get(h);
        const auto d = dst.cells().get(*res.cellMap.find(h));
        CHECK(s.cls.dim == d.cls.dim && s.cls.tag == d.cls.tag);
    }
    for (const topo::EdgeHandle& h : src->edges().activeBySequence()) {
        const auto s = src->edges().get(h);
        const auto d = dst.edges().get(*res.edgeMap.find(h));
        CHECK(s.cls.dim == d.cls.dim && s.cls.tag == d.cls.tag);
        CHECK(s.requiredConstraint == d.requiredConstraint);
        CHECK(s.useCount == d.useCount);
    }
    for (const topo::FaceHandle& h : src->faces().activeBySequence()) {
        const auto s = src->faces().get(h);
        const auto d = dst.faces().get(*res.faceMap.find(h));
        CHECK(s.cls.dim == d.cls.dim && s.cls.tag == d.cls.tag);
        CHECK(s.shape == d.shape && s.cornerCount == d.cornerCount);
    }
}

static void testCompactionDeterministicAfterChurn() {
    TEST("topo.compaction.deterministic_after_churn");
    std::vector<topo::NodeHandle> live;
    topo::NodeHandle retired{}, freed{};
    auto src = c1::buildChurnedWithRetired(&live, &retired, &freed);
    CHECK(src->nodes().stats().retired == 1);        // RETIRED churn
    CHECK(src->nodes().stats().reusableFree >= 1);   // free churn
    auto r1 = topo::compactWorkspace(*src);
    auto r2 = topo::compactWorkspace(*src);
    // Deterministic modulo runtime ownership (fresh pool ids by design).
    CHECK(c1::observableDigest(*r1.workspace, false) ==
          c1::observableDigest(*r2.workspace, false));
    CHECK(r1.workspace->testV2Tokens() == r2.workspace->testV2Tokens());
    CHECK(r1.workspace->strictWorkspaceHash() ==
          r2.workspace->strictWorkspaceHash());
    // Semantics preserved across retired-slot churn; the destination is
    // dense and carries no retirement.
    CHECK(r1.workspace->semanticTopologyHash() ==
          src->semanticTopologyHash());
    CHECK(r1.workspace->nodes().stats().retired == 0);
    CHECK(r1.workspace->nodes().stats().reusableFree == 0);
    // ALL FOUR maps agree slot/generation-wise (C2 review: the earlier
    // form claimed this but checked only nodeMap).
    auto sameIdent = [](const auto& a, const auto& b) {
        return a.slot == b.slot && a.gen == b.gen;
    };
    for (const topo::NodeHandle& h : src->nodes().activeBySequence())
        CHECK(sameIdent(*r1.nodeMap.find(h), *r2.nodeMap.find(h)));
    for (const topo::CellHandle& h : src->cells().activeBySequence())
        CHECK(sameIdent(*r1.cellMap.find(h), *r2.cellMap.find(h)));
    for (const topo::EdgeHandle& h : src->edges().activeBySequence())
        CHECK(sameIdent(*r1.edgeMap.find(h), *r2.edgeMap.find(h)));
    for (const topo::FaceHandle& h : src->faces().activeBySequence())
        CHECK(sameIdent(*r1.faceMap.find(h), *r2.faceMap.find(h)));
    CHECK(r1.edgeMap.size() == src->edges().stats().live);
    CHECK(r1.faceMap.size() == src->faces().stats().live);
}

static void testCompactionMemoryReduced() {
    TEST("topo.compaction.memory_reduced");
    std::vector<topo::NodeHandle> live;
    topo::NodeHandle retired{}, freed{};
    auto src = c1::buildChurnedWithRetired(&live, &retired, &freed);
    auto res = topo::compactWorkspace(*src);
    CHECK(res.after.totalBytes < res.before.totalBytes);
    CHECK(res.savedBytes == res.before.totalBytes - res.after.totalBytes);
    CHECK(res.beneficial);                       // default policy
    CHECK(res.after.consistent() && res.before.consistent());
    // Thresholds can refuse a real but negligible gain.
    topo::CompactionOptions strict;
    strict.minSavedBytes = res.savedBytes + 1;
    auto res2 = topo::compactWorkspace(*src, strict);
    CHECK(res2.savedBytes == res.savedBytes);    // same measurement...
    CHECK(!res2.beneficial);                     // ...refused by policy
    CHECK(res2.workspace != nullptr);            // still returned
    topo::CompactionOptions frac;
    frac.minSavedFraction = 1.0;                 // impossible
    CHECK(!topo::compactWorkspace(*src, frac).beneficial);
}

static void testCompactionNoBenefitReported() {
    TEST("topo.compaction.no_benefit_reported");
    // An already-canonical workspace has nothing to reclaim: the result
    // reports beneficial=false and savedBytes=0 - NOT an error - and
    // still returns a valid compacted workspace.
    auto src = c1::buildCanonical();
    auto res = topo::compactWorkspace(*src);
    CHECK(!res.beneficial);
    CHECK(res.savedBytes == 0);
    CHECK(res.after.totalBytes >= res.before.totalBytes);
    CHECK(res.workspace != nullptr);
    CHECK(res.workspace->semanticTopologyHash() ==
          src->semanticTopologyHash());
    CHECK(res.workspace->state() == topo::WorkspaceState::Ready);
    CHECK(res.nodeMap.size() == src->nodes().stats().live);
    CHECK(res.edgeMap.size() == src->edges().stats().live);
}

static void testCompactionOptionDomain() {
    TEST("topo.compaction.option_domain_validated");
    std::vector<topo::NodeHandle> live;
    topo::NodeHandle retired{}, freed{};
    auto src = c1::buildChurnedWithRetired(&live, &retired, &freed);
    const double nan = std::nan("");
    const double inf = std::numeric_limits<double>::infinity();
    // Rejected BEFORE construction, with invalid_argument.
    for (const double bad : {-0.001, -1.0, nan, inf, -inf, 1.0001, 2.0}) {
        topo::CompactionOptions o;
        o.minSavedFraction = bad;
        bool rejected = false;
        try { (void)topo::compactWorkspace(*src, o); }
        catch (const std::invalid_argument&) { rejected = true; }
        CHECK(rejected);
    }
    // Boundary values are ACCEPTED.
    for (const double ok : {0.0, 1.0}) {
        topo::CompactionOptions o;
        o.minSavedFraction = ok;
        CHECK(topo::compactWorkspace(*src, o).workspace != nullptr);
    }
    // Exact threshold accepted; just above it rejected.
    auto base = topo::compactWorkspace(*src);
    CHECK(base.savedBytes > 0);
    const double frac = static_cast<double>(base.savedBytes) /
                        static_cast<double>(base.before.totalBytes);
    {
        topo::CompactionOptions exact;
        exact.minSavedFraction = frac;
        CHECK(topo::compactWorkspace(*src, exact).beneficial);
    }
    {
        topo::CompactionOptions above;
        above.minSavedFraction = std::nextafter(frac, 2.0);
        CHECK(!topo::compactWorkspace(*src, above).beneficial);
    }
    {
        topo::CompactionOptions exactBytes;
        exactBytes.minSavedBytes = base.savedBytes;
        CHECK(topo::compactWorkspace(*src, exactBytes).beneficial);
        exactBytes.minSavedBytes = base.savedBytes + 1;
        CHECK(!topo::compactWorkspace(*src, exactBytes).beneficial);
    }
}

// --- M4-B.8 STAGE A fixtures ------------------------------------------------
namespace cav {
using topo::NodePool; using topo::CellPool; using topo::EdgePool;
using topo::FacePool; using topo::TopologyWorkspace;
using topo::ElementType; using topo::NodeHandle; using topo::CellHandle;
using topo::TriCavityCell; using topo::TriCavitySpec;
using topo::CavityEdgeRequirement;

struct Fix {
    std::unique_ptr<TopologyWorkspace> ws;
    std::vector<NodeHandle> n;
    std::vector<CellHandle> c;
};
// Quad n0..n3 with diagonal n0n2: triangles A(n0,n1,n2), C(n0,n2,n3),
// consistently oriented. extraNodes appends isolated-but-used nodes for
// outside triangles per test.
inline Fix buildQuad(int extraNodes = 0,
                     bool outsideTriOnN1N3 = false,
                     bool markDiagonalRequired = false,
                     paacs::mesh::Classification edge01Cls = {},
                     bool edge01Required = false,
                     bool thirdTriOnDiagonal = false) {
    Fix f;
    NodePool np;
    f.n.push_back(np.add({0, 0, 0}, {2, 1}));
    f.n.push_back(np.add({1, 0, 0}, {2, 1}));
    f.n.push_back(np.add({1, 1, 0}, {2, 1}));
    f.n.push_back(np.add({0, 1, 0}, {2, 1}));
    for (int i = 0; i < extraNodes; ++i)
        f.n.push_back(np.add({2. + i, 0.5, 0}, {2, 1}));
    CellPool cp;
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[1], f.n[2]}, {2, 5}, np));
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[2], f.n[3]}, {2, 5}, np));
    EdgePool ep;
    ep.insertCellEdges(np, cp, f.c[0]);
    ep.insertCellEdges(np, cp, f.c[1]);
    if (outsideTriOnN1N3) {
        f.c.push_back(cp.add(ElementType::Tri3,
            std::vector<NodeHandle>{f.n[1], f.n[3], f.n[4]}, {2, 6}, np));
        ep.insertCellEdges(np, cp, f.c.back());
    }
    if (thirdTriOnDiagonal) {
        f.c.push_back(cp.add(ElementType::Tri3,
            std::vector<NodeHandle>{f.n[0], f.n[2], f.n[4]}, {2, 6}, np));
        ep.insertCellEdges(np, cp, f.c.back());
    }
    if (markDiagonalRequired)
        ep.setRequiredConstraint(*ep.find(np, f.n[0], f.n[2]), true);
    if (edge01Cls.dim != -1 || edge01Cls.tag != -1)
        ep.setClassification(*ep.find(np, f.n[0], f.n[1]), edge01Cls);
    if (edge01Required)
        ep.setRequiredConstraint(*ep.find(np, f.n[0], f.n[1]), true);
    f.ws = std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), FacePool{});
    return f;
}
// Diagonal-swap spec: R = {A, C}, M = {(n0,n1,n3), (n1,n2,n3)}.
inline TriCavitySpec swapSpec(const Fix& f) {
    TriCavitySpec s;
    s.removed = {f.c[0], f.c[1]};
    s.added = {TriCavityCell{{f.n[0], f.n[1], f.n[3]}, {2, 5}},
               TriCavityCell{{f.n[1], f.n[2], f.n[3]}, {2, 5}}};
    return s;
}
inline bool rejWith(const std::function<void()>& fn, const char* needle) {
    try { fn(); }
    catch (const std::invalid_argument& e) {
        return std::string(e.what()).find(needle) != std::string::npos;
    }
    return false;
}
}  // namespace cav

// ---------------------------------------------------------------------------
static void testCavityStateMachine() {
    TEST("topo.cavity.state_machine (begin/single-active/abandon/destroy)");
    auto f = cav::buildQuad();
    CHECK(f.ws->state() == topo::WorkspaceState::Ready);
    {
        auto tx = f.ws->beginCavityTransaction();
        CHECK(f.ws->state() == topo::WorkspaceState::TransactionStaging);
        bool threw = false;                     // second cavity txn
        try { (void)f.ws->beginCavityTransaction(); }
        catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
        threw = false;                          // ordinary txn also blocked
        try { (void)f.ws->beginTransaction(); }
        catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
        tx.abandon();
        CHECK(f.ws->state() == topo::WorkspaceState::Ready);
        threw = false;                          // use after abandon
        try { tx.stage(cav::swapSpec(f)); }
        catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
    }
    {                                           // destruction restores
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::swapSpec(f));
        CHECK(f.ws->state() == topo::WorkspaceState::TransactionStaging);
    }
    CHECK(f.ws->state() == topo::WorkspaceState::Ready);
    auto tx2 = f.ws->beginTransaction();        // ordinary txn works again
    tx2.abandon();
}

static void testCavityStageValidPreview() {
    TEST("topo.cavity.stage_valid_preview (diagonal swap plan)");
    auto f = cav::buildQuad();
    auto tx = f.ws->beginCavityTransaction();
    const auto& p = tx.stage(cav::swapSpec(f));
    CHECK(p.cellsRemoved == 2 && p.cellsAdded == 2);
    CHECK(p.edgesToRetain == 4);                // the quad boundary
    CHECK(p.edgesToRetire == 1);                // diagonal n0n2
    CHECK(p.edgesToCreate == 1);                // diagonal n1n3
    CHECK(p.interiorNodes.empty());
    CHECK(p.boundaryCycle.size() == 4);
    CHECK(p.boundaryCycle[0] == f.n[0]);        // lowest rank first
    CHECK(p.boundaryCycle[1] == f.n[1]);        // lexicographically
    CHECK(p.boundaryCycle[2] == f.n[2]);        // smaller direction
    CHECK(p.boundaryCycle[3] == f.n[3]);
    CHECK(!p.semanticNoop);
    for (const auto& row : p.edges) {           // net-delta rows
        if (row.a == f.n[0] && row.b == f.n[2])
            CHECK(row.initialUse == 2 && row.finalUse == 0);
        if (row.a == f.n[1] && row.b == f.n[3])
            CHECK(row.initialUse == 0 && row.finalUse == 2);
    }
    tx.abandon();
}

static void testCavityCanonicalizationIgnoresCallerOrder() {
    TEST("topo.cavity.canonicalization_ignores_caller_order");
    auto f = cav::buildQuad();
    auto spec1 = cav::swapSpec(f);
    auto spec2 = spec1;                          // shuffled + rotated
    std::swap(spec2.removed[0], spec2.removed[1]);
    std::swap(spec2.added[0], spec2.added[1]);
    spec2.added[0].nodes = {f.n[2], f.n[3], f.n[1]};   // rotation of
    spec2.added[1].nodes = {f.n[1], f.n[3], f.n[0]};   // the same tris
    auto tx1 = f.ws->beginCavityTransaction();
    const auto p1 = tx1.stage(spec1);
    tx1.abandon();
    auto tx2 = f.ws->beginCavityTransaction();
    const auto p2 = tx2.stage(spec2);
    tx2.abandon();
    CHECK(p1.boundaryCycle == p2.boundaryCycle);
    CHECK(p1.removedCells == p2.removedCells);
    CHECK(p1.addedConnectivity == p2.addedConnectivity);
    CHECK(p1.edges.size() == p2.edges.size());
    for (std::size_t i = 0; i < p1.edges.size(); ++i) {
        CHECK(p1.edges[i].a == p2.edges[i].a);
        CHECK(p1.edges[i].b == p2.edges[i].b);
        CHECK(p1.edges[i].finalUse == p2.edges[i].finalUse);
    }
}

static void testCavityValidatorRejection() {
    TEST("topo.cavity.validator_rejection (bool pattern; plain values)");
    auto f = cav::buildQuad();
    auto tx = f.ws->beginCavityTransaction();
    (void)tx.stage(cav::swapSpec(f));
    // The validator sees ONLY the plain-value preview; returning false
    // is ordinary control flow and marks the transaction rejected.
    CHECK(!tx.rejected());
    const bool ok = tx.validate([](const topo::CavityPreview& p) {
        return p.edgesToCreate == 99;            // pretend BND-1 said no
    });
    CHECK(!ok);
    CHECK(tx.rejected());
    // FROZEN semantics, directly (closure sync): a later validate(true)
    // returns TRUE (this call's result) while rejection stays STICKY.
    CHECK(tx.validate([](const topo::CavityPreview&) { return true; }));
    CHECK(tx.rejected());                    // still rejected
    tx.abandon();
    auto tx2 = f.ws->beginCavityTransaction();
    (void)tx2.stage(cav::swapSpec(f));
    CHECK(tx2.validate([](const topo::CavityPreview&) { return true; }));
    CHECK(!tx2.rejected());
    tx2.abandon();
}

static void testCavityPhase1DigestUnchanged() {
    TEST("topo.cavity.phase1_full_digest_unchanged (after abandonment)");
    auto f = cav::buildQuad();                   // plain valid quad
    f.ws->rebuildAdjacency();
    const auto d0 = c1::observableDigest(*f.ws);
    {                                            // staged + rejected + gone
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::swapSpec(f));
        (void)tx.validate([](const topo::CavityPreview&) { return false; });
    }
    CHECK(c1::observableDigest(*f.ws) == d0);
    {                                            // invalid staging attempt
        auto tx = f.ws->beginCavityTransaction();
        topo::TriCavitySpec bad;                 // empty removed
        bool threw = false;
        try { (void)tx.stage(bad); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }
    CHECK(c1::observableDigest(*f.ws) == d0);
}

static void testCavityDiskInvariantEnforced() {
    TEST("topo.cavity.disk_invariant_enforced (disconnected R rejected)");
    // Two triangles sharing nothing: not edge-connected.
    topo::NodePool np;
    std::vector<topo::NodeHandle> n;
    for (int i = 0; i < 6; ++i)
        n.push_back(np.add({double(i), double(i % 2), 0}, {2, 1}));
    topo::CellPool cp;
    auto a = cp.add(topo::ElementType::Tri3,
        std::vector<topo::NodeHandle>{n[0], n[1], n[2]}, {2, 5}, np);
    auto b = cp.add(topo::ElementType::Tri3,
        std::vector<topo::NodeHandle>{n[3], n[4], n[5]}, {2, 5}, np);
    topo::EdgePool ep;
    ep.insertCellEdges(np, cp, a);
    ep.insertCellEdges(np, cp, b);
    topo::TopologyWorkspace ws(std::move(np), std::move(cp), std::move(ep),
                               topo::FacePool{});
    auto tx = ws.beginCavityTransaction();
    topo::TriCavitySpec s;
    s.removed = {a, b};
    s.added = {topo::TriCavityCell{{n[0], n[1], n[2]}, {2, 5}},
               topo::TriCavityCell{{n[3], n[4], n[5]}, {2, 5}}};
    CHECK(cav::rejWith([&] { (void)tx.stage(s); }, "edge-connected"));
}

static void testCavityAnnulusRejected() {
    TEST("topo.cavity.boundary_cycle_is_simple (annulus: two loops)");
    // Triangulated square ring: 8 outer/inner nodes, 8 triangles,
    // TWO boundary cycles, Euler characteristic 0.
    topo::NodePool np;
    std::vector<topo::NodeHandle> o, in;
    o.push_back(np.add({0, 0, 0}, {2, 1}));
    o.push_back(np.add({4, 0, 0}, {2, 1}));
    o.push_back(np.add({4, 4, 0}, {2, 1}));
    o.push_back(np.add({0, 4, 0}, {2, 1}));
    in.push_back(np.add({1, 1, 0}, {2, 1}));
    in.push_back(np.add({3, 1, 0}, {2, 1}));
    in.push_back(np.add({3, 3, 0}, {2, 1}));
    in.push_back(np.add({1, 3, 0}, {2, 1}));
    topo::CellPool cp;
    std::vector<topo::CellHandle> cells;
    auto tri = [&](topo::NodeHandle x, topo::NodeHandle y,
                   topo::NodeHandle z) {
        cells.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{x, y, z}, {2, 5}, np));
    };
    for (int i = 0; i < 4; ++i) {
        const int j = (i + 1) % 4;
        tri(o[i], o[j], in[j]);
        tri(o[i], in[j], in[i]);
    }
    topo::EdgePool ep;
    for (const auto& c : cells) ep.insertCellEdges(np, cp, c);
    topo::TopologyWorkspace ws(std::move(np), std::move(cp), std::move(ep),
                               topo::FacePool{});
    auto tx = ws.beginCavityTransaction();
    topo::TriCavitySpec s;
    s.removed = cells;
    for (int i = 0; i < 4; ++i) {
        const int j = (i + 1) % 4;
        s.added.push_back(topo::TriCavityCell{{o[i], o[j], in[j]}, {2, 5}});
        s.added.push_back(topo::TriCavityCell{{o[i], in[j], in[i]}, {2, 5}});
    }
    // Review blocker 4: this fixture has TWO loops, so the simple-cycle
    // rule fires before Euler ever runs - pin THAT diagnostic only; the
    // chi rule has its own fixture (euler_characteristic_enforced).
    CHECK(cav::rejWith([&] { (void)tx.stage(s); }, "one simple cycle"));
}

static void testCavitySourceNonmanifoldEdgeRejected() {
    TEST("topo.cavity.source_nonmanifold_cavity_edge_rejected");
    auto f = cav::buildQuad(1, false, false, {}, false,
                            /*thirdTriOnDiagonal=*/true);
    CHECK(f.ws != nullptr);                      // ctor accepted it
    auto tx = f.ws->beginCavityTransaction();
    // Any cavity whose closure touches the 3-way diagonal must reject
    // BEFORE any classification is even attempted.
    CHECK(cav::rejWith([&] { (void)tx.stage(cav::swapSpec(f)); },
                       "initialUse > 2"));
}

static void testCavityInteriorEdgeOutsideUseRejected() {
    TEST("topo.cavity.candidate_interior_edge_with_outside_use_rejected");
    auto f = cav::buildQuad(1, /*outsideTriOnN1N3=*/true);
    // A surviving outside triangle already uses n1n3; M would make it
    // interior -> finalUse would be 3 -> nonmanifold. nodes(M) is a
    // subset of nodes(R), so ONLY the outsideUse rule can catch this.
    auto tx = f.ws->beginCavityTransaction();
    CHECK(cav::rejWith([&] { (void)tx.stage(cav::swapSpec(f)); },
                       "outsideUse != 0"));
}

static void testCavityOrphanNodeRejected() {
    TEST("topo.cavity.orphan_node_rejected (unconditional)");
    // Fan of four triangles around center n4; M drops the center.
    topo::NodePool np;
    std::vector<topo::NodeHandle> n;
    n.push_back(np.add({0, 0, 0}, {2, 1}));
    n.push_back(np.add({1, 0, 0}, {2, 1}));
    n.push_back(np.add({1, 1, 0}, {2, 1}));
    n.push_back(np.add({0, 1, 0}, {2, 1}));
    n.push_back(np.add({.5, .5, 0}, {2, 1}));
    topo::CellPool cp;
    std::vector<topo::CellHandle> cells;
    for (int i = 0; i < 4; ++i)
        cells.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{n[i], n[(i + 1) % 4], n[4]},
            {2, 5}, np));
    topo::EdgePool ep;
    for (const auto& c : cells) ep.insertCellEdges(np, cp, c);
    topo::TopologyWorkspace ws(std::move(np), std::move(cp), std::move(ep),
                               topo::FacePool{});
    auto tx = ws.beginCavityTransaction();
    topo::TriCavitySpec s;
    s.removed = cells;
    s.added = {topo::TriCavityCell{{n[0], n[1], n[2]}, {2, 5}},
               topo::TriCavityCell{{n[0], n[2], n[3]}, {2, 5}}};
    CHECK(cav::rejWith([&] { (void)tx.stage(s); }, "orphan"));
}

static void testCavityOutsideNodePulledInRejected() {
    TEST("topo.cavity.outside_node_pulled_in_rejected");
    auto f = cav::buildQuad(1);
    auto tx = f.ws->beginCavityTransaction();
    auto s = cav::swapSpec(f);
    s.added[0].nodes = {f.n[0], f.n[1], f.n[4]};   // n4 not in nodes(R)
    CHECK(cav::rejWith([&] { (void)tx.stage(s); }, "outside the"));
}

static void testCavityDuplicateOrStaleRemovalRejected() {
    TEST("topo.cavity.duplicate_or_stale_removal_rejected");
    auto f = cav::buildQuad();
    {
        auto tx = f.ws->beginCavityTransaction();
        auto s = cav::swapSpec(f);
        s.removed = {f.c[0], f.c[0]};
        CHECK(cav::rejWith([&] { (void)tx.stage(s); }, "duplicates"));
    }
    {
        auto tx = f.ws->beginCavityTransaction();
        auto s = cav::swapSpec(f);
        topo::CellHandle stale = f.c[0];
        stale.gen += 1;
        s.removed = {stale, f.c[1]};
        CHECK(cav::rejWith([&] { (void)tx.stage(s); }, "not current"));
    }
}

static void testCavityBoundaryViolationRejected() {
    TEST("topo.cavity.boundary_violation_rejected");
    auto f = cav::buildQuad();
    auto tx = f.ws->beginCavityTransaction();
    topo::TriCavitySpec s;                       // M covers half of R
    s.removed = {f.c[0], f.c[1]};
    s.added = {topo::TriCavityCell{{f.n[0], f.n[1], f.n[2]}, {2, 5}}};
    CHECK(cav::rejWith([&] { (void)tx.stage(s); },
                       "boundary(M) != boundary(R)"));
}

static void testCavityInsideSenseMismatchRejected() {
    TEST("topo.cavity.inside_sense_or_orientation_rejected");
    auto f = cav::buildQuad();
    // Flip one candidate triangle: either the interior opposite-sense
    // rule or the boundary inside-sense rule must fire.
    auto tx = f.ws->beginCavityTransaction();
    auto s = cav::swapSpec(f);
    std::swap(s.added[0].nodes[0], s.added[0].nodes[1]);   // flip
    CHECK(cav::rejWith([&] { (void)tx.stage(s); }, "sense"));
}

static void testCavityRequirementRules() {
    TEST("topo.cavity.requirement rules (ensure-only contract)");
    using RQ = topo::CavityEdgeRequirement;
    {   // valid: require the NEW diagonal; preview reflects it
        auto f = cav::buildQuad();
        auto tx = f.ws->beginCavityTransaction();
        auto s = cav::swapSpec(f);
        s.requirements = {RQ{f.n[1], f.n[3], {1, 7}}};
        const auto& p = tx.stage(s);
        bool found = false;
        for (const auto& r : p.edges)
            if (r.a == f.n[1] && r.b == f.n[3]) {
                found = true;
                CHECK(r.willBeRequired);
                CHECK(r.finalCls.dim == 1 && r.finalCls.tag == 7);
            }
        CHECK(found);
        tx.abandon();
    }
    {   // absent-from-final: requiring the RETIRING diagonal
        auto f = cav::buildQuad();
        auto tx = f.ws->beginCavityTransaction();
        auto s = cav::swapSpec(f);
        s.requirements = {RQ{f.n[0], f.n[2], {1, 7}}};
        CHECK(cav::rejWith([&] { (void)tx.stage(s); }, "absent"));
    }
    {   // duplicate contradictory
        auto f = cav::buildQuad();
        auto tx = f.ws->beginCavityTransaction();
        auto s = cav::swapSpec(f);
        s.requirements = {RQ{f.n[1], f.n[3], {1, 7}},
                          RQ{f.n[3], f.n[1], {1, 8}}};
        CHECK(cav::rejWith([&] { (void)tx.stage(s); }, "contradictory"));
    }
    {   // conflicting non-default classification on an existing edge
        auto f = cav::buildQuad(0, false, false, {1, 9});
        auto tx = f.ws->beginCavityTransaction();
        auto s = cav::swapSpec(f);
        s.requirements = {RQ{f.n[0], f.n[1], {1, 8}}};
        CHECK(cav::rejWith([&] { (void)tx.stage(s); }, "conflicting"));
    }
    {   // equal classification on that edge: accepted (idempotent)
        auto f = cav::buildQuad(0, false, false, {1, 9});
        auto tx = f.ws->beginCavityTransaction();
        auto s = cav::swapSpec(f);
        s.requirements = {RQ{f.n[0], f.n[1], {1, 9}}};
        (void)tx.stage(s);
        tx.abandon();
    }
    {   // already REQUIRED with different cls: rejected
        auto f = cav::buildQuad(0, false, false, {1, 9}, true);
        auto tx = f.ws->beginCavityTransaction();
        auto s = cav::swapSpec(f);
        s.requirements = {RQ{f.n[0], f.n[1], {1, 8}}};
        CHECK(cav::rejWith([&] { (void)tx.stage(s); },
                           "already constrained"));
    }
}

static void testCavityRequiredEdgeMayNotBeRetired() {
    TEST("topo.cavity.required_edge_may_not_be_retired");
    auto f = cav::buildQuad(0, false, /*markDiagonalRequired=*/true);
    auto tx = f.ws->beginCavityTransaction();
    // The swap retires n0n2, which is marked required -> reject.
    CHECK(cav::rejWith([&] { (void)tx.stage(cav::swapSpec(f)); },
                       "required edge would be retired"));
}

static void testCavitySemanticNoopDetected() {
    TEST("topo.cavity.semantic_noop_detected");
    auto f = cav::buildQuad();
    {   // identical oriented triangles + cls: no-op
        auto tx = f.ws->beginCavityTransaction();
        topo::TriCavitySpec s;
        s.removed = {f.c[0], f.c[1]};
        s.added = {topo::TriCavityCell{{f.n[1], f.n[2], f.n[0]}, {2, 5}},
                   topo::TriCavityCell{{f.n[2], f.n[3], f.n[0]}, {2, 5}}};
        CHECK(tx.stage(s).semanticNoop);
        tx.abandon();
    }
    {   // different classification: NOT a no-op
        auto tx = f.ws->beginCavityTransaction();
        topo::TriCavitySpec s;
        s.removed = {f.c[0], f.c[1]};
        s.added = {topo::TriCavityCell{{f.n[0], f.n[1], f.n[2]}, {2, 8}},
                   topo::TriCavityCell{{f.n[0], f.n[2], f.n[3]}, {2, 5}}};
        CHECK(!tx.stage(s).semanticNoop);
        tx.abandon();
    }
    {   // unsatisfied requirement: NOT a no-op
        auto tx = f.ws->beginCavityTransaction();
        topo::TriCavitySpec s;
        s.removed = {f.c[0], f.c[1]};
        s.added = {topo::TriCavityCell{{f.n[0], f.n[1], f.n[2]}, {2, 5}},
                   topo::TriCavityCell{{f.n[0], f.n[2], f.n[3]}, {2, 5}}};
        s.requirements = {topo::CavityEdgeRequirement{f.n[0], f.n[2],
                                                      {1, 7}}};
        CHECK(!tx.stage(s).semanticNoop);
        tx.abandon();
    }
}

// --- Stage-A repair round: closed/pinched/torus adversaries ----------------
namespace cav {
// Closed oriented octahedron: manifold everywhere, EMPTY boundary.
inline Fix buildOctahedron() {
    Fix f;
    NodePool np;
    for (int i = 0; i < 6; ++i)
        f.n.push_back(np.add({double(i), double(i % 2), double(i % 3)},
                             {2, 1}));
    static const int F[8][3] = {{0,1,2},{0,2,3},{0,3,4},{0,4,1},
                                {5,2,1},{5,3,2},{5,4,3},{5,1,4}};
    CellPool cp;
    for (const auto& fc : F)
        f.c.push_back(cp.add(ElementType::Tri3,
            std::vector<NodeHandle>{f.n[fc[0]], f.n[fc[1]], f.n[fc[2]]},
            {2, 5}, np));
    EdgePool ep;
    for (const auto& c : f.c) ep.insertCellEdges(np, cp, c);
    f.ws = std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), FacePool{});
    return f;
}
// Pinched punctured icosahedron: vertex 3 identified into vertex 0
// (disjoint neighbourhoods, so every edge stays use<=2), one face away
// from the pinch removed (nonempty 3-cycle boundary). The pinched
// vertex's link is TWO disjoint 5-cycles. chi = 0, but the node-link
// audit must fire FIRST - Euler alone cannot see a pinch that costs
// exactly the vertex a puncture refunds.
inline Fix buildPinchedIcosahedron() {
    Fix f;
    NodePool np;
    std::vector<NodeHandle> v;
    for (int i = 0; i < 12; ++i)
        v.push_back(np.add({double(i), double(i % 3), double(i % 5)},
                           {2, 1}));
    static const int F[20][3] = {
        {0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},
        {1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},
        {3,9,4},{3,4,2},{3,2,6},{3,6,8},{3,8,9},
        {4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1}};
    CellPool cp;
    for (int i = 0; i < 20; ++i) {
        if (i == 5) continue;                 // puncture {1,5,9}
        auto id = [&](int k) { return k == 3 ? v[0] : v[k]; };  // pinch
        f.c.push_back(cp.add(ElementType::Tri3,
            std::vector<NodeHandle>{id(F[i][0]), id(F[i][1]), id(F[i][2])},
            {2, 5}, np));
    }
    f.n = v;
    EdgePool ep;
    for (const auto& c : f.c) ep.insertCellEdges(np, cp, c);
    f.ws = std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), FacePool{});
    return f;
}
// Punctured 3x3 grid torus: manifold node links, ONE simple 3-cycle
// boundary, Euler characteristic -1.
inline Fix buildPuncturedTorus() {
    Fix f;
    NodePool np;
    std::vector<NodeHandle> v(9);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            v[std::size_t(3 * i + j)] =
                np.add({double(i), double(j), 0}, {2, 1});
    f.n = v;
    CellPool cp;
    bool skippedOne = false;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            const auto a = v[std::size_t(3 * i + j)];
            const auto b = v[std::size_t(3 * ((i + 1) % 3) + j)];
            const auto cN = v[std::size_t(3 * ((i + 1) % 3) + (j + 1) % 3)];
            const auto d = v[std::size_t(3 * i + (j + 1) % 3)];
            if (!skippedOne) { skippedOne = true; }   // puncture T1(0,0)
            else
                f.c.push_back(cp.add(ElementType::Tri3,
                    std::vector<NodeHandle>{a, b, cN}, {2, 5}, np));
            f.c.push_back(cp.add(ElementType::Tri3,
                std::vector<NodeHandle>{a, cN, d}, {2, 5}, np));
        }
    EdgePool ep;
    for (const auto& c : f.c) ep.insertCellEdges(np, cp, c);
    f.ws = std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), FacePool{});
    return f;
}
inline topo::TriCavitySpec readdAll(const Fix& f) {
    topo::TriCavitySpec s;
    s.removed = f.c;
    for (const auto& h : f.c) {
        const auto sn = f.ws->cells().get(h);
        s.added.push_back(topo::TriCavityCell{
            {sn.nodes[0], sn.nodes[1], sn.nodes[2]}, sn.cls});
    }
    return s;
}
}  // namespace cav

static void testCavityBoundaryCycleNonempty() {
    TEST("topo.cavity.boundary_cycle_nonempty (closed octahedron rejected)");
    auto f = cav::buildOctahedron();
    const auto spec = cav::readdAll(f);   // BEFORE begin: the frozen
                                          // visibility contract makes
                                          // ws accessors throw while a
                                          // transaction is active
    auto tx = f.ws->beginCavityTransaction();
    CHECK(cav::rejWith([&] { (void)tx.stage(spec); },
                       "boundary is empty"));
}

static void testCavityNodeLinkManifoldEnforced() {
    TEST("topo.cavity.node_link_manifold_enforced (pinched vertex)");
    auto f = cav::buildPinchedIcosahedron();
    CHECK(f.ws != nullptr);                     // ctor accepted the pinch
    const auto spec = cav::readdAll(f);       // BEFORE begin (visibility)
    auto tx = f.ws->beginCavityTransaction();
    // chi == 0 here, but the LINK audit must fire, not Euler: the pinch
    // has two disjoint 5-cycles at the identified vertex.
    CHECK(cav::rejWith([&] { (void)tx.stage(spec); }, "link"));
}

static void testCavityEulerCharacteristicEnforced() {
    TEST("topo.cavity.euler_characteristic_enforced (punctured torus)");
    auto f = cav::buildPuncturedTorus();
    const auto spec = cav::readdAll(f);       // BEFORE begin (visibility)
    auto tx = f.ws->beginCavityTransaction();
    // Connectivity, node links, nonempty boundary and the single simple
    // 3-cycle ALL pass on this fixture; ONLY chi == -1 can reject it.
    CHECK(cav::rejWith([&] { (void)tx.stage(spec); }, "Euler"));
}

static void testCavitySurvivingDuplicateRejected() {
    TEST("topo.cavity.surviving_duplicate_rejected (both orientations)");
    {   // REVERSED-orientation survivor: (a,c,b) vs candidate (a,b,c)
        topo::NodePool np;
        std::vector<topo::NodeHandle> n;
        for (int i = 0; i < 3; ++i)
            n.push_back(np.add({double(i), double(i % 2), 0}, {2, 1}));
        topo::CellPool cp;
        auto A = cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{n[0], n[1], n[2]}, {2, 5}, np);
        auto B = cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{n[0], n[2], n[1]}, {2, 5}, np);
        (void)B;
        topo::EdgePool ep;
        ep.insertCellEdges(np, cp, A);
        ep.insertCellEdges(np, cp, B);
        topo::TopologyWorkspace ws(std::move(np), std::move(cp),
                                   std::move(ep), topo::FacePool{});
        auto tx = ws.beginCavityTransaction();
        topo::TriCavitySpec s;
        s.removed = {A};
        s.added = {topo::TriCavityCell{{n[0], n[1], n[2]}, {2, 5}}};
        CHECK(cav::rejWith([&] { (void)tx.stage(s); },
                           "duplicates a surviving cell"));
    }
    {   // SAME-orientation survivor
        topo::NodePool np;
        std::vector<topo::NodeHandle> n;
        for (int i = 0; i < 3; ++i)
            n.push_back(np.add({double(i), double(i % 2), 0}, {2, 1}));
        topo::CellPool cp;
        auto A = cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{n[0], n[1], n[2]}, {2, 5}, np);
        auto B = cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{n[0], n[1], n[2]}, {2, 5}, np);
        (void)B;
        topo::EdgePool ep;
        ep.insertCellEdges(np, cp, A);
        ep.insertCellEdges(np, cp, B);
        topo::TopologyWorkspace ws(std::move(np), std::move(cp),
                                   std::move(ep), topo::FacePool{});
        auto tx = ws.beginCavityTransaction();
        topo::TriCavitySpec s;
        s.removed = {A};
        s.added = {topo::TriCavityCell{{n[0], n[1], n[2]}, {2, 5}}};
        // Either the duplicate rule or the no-op/sense structure may
        // reject; the DUPLICATE rule is the one that must exist, so pin
        // its message.
        CHECK(cav::rejWith([&] { (void)tx.stage(s); },
                           "duplicates a surviving cell"));
    }
}

static void testCavityCompositeGateEnforced() {
    TEST("topo.cavity.composite_gate_enforced (corrupted table refused)");
    auto f = cav::buildQuad();
    const auto diag = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[2]);
    f.ws->testEdgesMutable().testDropTableEntry(diag);
    auto tx = f.ws->beginCavityTransaction();
    // Pinned EXCLUSIVELY to the gate diagnostic (closure sync): if the
    // composite gate were deleted, the later unconditional guard would
    // otherwise keep this green. That guard stays as defense-in-depth
    // but is NOT claimed behaviourally proven by this fixture.
    bool gateFired = false;
    try { (void)tx.stage(cav::swapSpec(f)); }
    catch (const std::logic_error& e) {
        gateFired = std::string(e.what()).find("topology-ready") !=
                    std::string::npos;
    }
    CHECK(gateFired);
}

static void testCavityStateMachineExtended() {
    TEST("topo.cavity.state_machine_extended (move/mutex/stale-restore)");
    {   // move construction transfers the active transaction
        auto f = cav::buildQuad();
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::swapSpec(f));
        auto tx2 = std::move(tx);
        CHECK(f.ws->state() == topo::WorkspaceState::TransactionStaging);
        bool threw = false;                  // source is finished
        try { (void)tx.preview(); }
        catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
        CHECK(tx2.preview().cellsRemoved == 2);   // target is live
        tx2.abandon();
        CHECK(f.ws->state() == topo::WorkspaceState::Ready);
    }
    {   // ordinary transaction blocks a cavity transaction
        auto f = cav::buildQuad();
        auto tx = f.ws->beginTransaction();
        bool threw = false;
        try { (void)f.ws->beginCavityTransaction(); }
        catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
        tx.abandon();
    }
    {   // restoration to CommittedAdjacencyStale, not blindly Ready
        auto f = cav::buildQuad(1);
        f.ws->rebuildAdjacency();
        auto t1 = f.ws->beginTransaction();
        topo::TopologyPolicy orphanOk;
        orphanOk.permitOrphanNodes = true;
        t1.abandon();
        auto t2 = f.ws->beginTransaction(orphanOk);
        t2.addNode({9, 9, 9}, {0, 7});
        (void)t2.commit();
        CHECK(f.ws->state() ==
              topo::WorkspaceState::CommittedAdjacencyStale);
        {
            auto cx = f.ws->beginCavityTransaction();
            CHECK(f.ws->state() ==
                  topo::WorkspaceState::TransactionStaging);
        }
        CHECK(f.ws->state() ==
              topo::WorkspaceState::CommittedAdjacencyStale);
    }
}

static void testCavityRequirementOutsideClosureRejected() {
    TEST("topo.cavity.requirement_outside_closure_rejected");
    auto f = cav::buildQuad(1, true);       // outside tri (n1,n3,n4)
    auto tx = f.ws->beginCavityTransaction();
    topo::TriCavitySpec s;                  // cavity touches ONLY c0
    s.removed = {f.c[0]};
    s.added = {topo::TriCavityCell{{f.n[0], f.n[1], f.n[2]}, {2, 8}}};
    // wait: that duplicates nothing and keeps boundary; but requirement
    // names edge (n3,n4) - live, valid, and entirely OUTSIDE the closure.
    s.requirements = {topo::CavityEdgeRequirement{f.n[3], f.n[4], {1, 7}}};
    CHECK(cav::rejWith([&] { (void)tx.stage(s); },
                       "outside the cavity closure"));
}

// --- M4-B.8 STAGE B: the Phase-2 plan --------------------------------------
static void testCavityPlanSlotReuseExact() {
    TEST("topo.cavity.plan_slot_reuse_exact (balanced by Euler)");
    auto f = cav::buildQuad();
    auto tx = f.ws->beginCavityTransaction();
    (void)tx.stage(cav::swapSpec(f));
    const auto& p = tx.computePlan();
    // Cells: 2 removed, 2 added, no exhaustion -> ALL demand satisfied
    // by in-commit reuse; ZERO fresh slots, ZERO free-list consumption.
    CHECK(p.cells.demand == 2 && p.cells.plannedRetirements == 2);
    CHECK(p.cells.exhaustedRetirements == 0);
    CHECK(p.cells.reusableThisCommit == 2);
    CHECK(p.cells.inCommitReuse == 2);
    CHECK(p.cells.phase2Claims == 0);
    CHECK(p.cells.preExistingConsumed == 0);
    CHECK(p.cells.freshSlots == 0);
    CHECK(p.cells.balanced());
    // Edges: retire the old diagonal, create the new one - also 1:1.
    CHECK(p.edges.demand == 1 && p.edges.plannedRetirements == 1);
    CHECK(p.edges.inCommitReuse == 1);
    CHECK(p.edges.freshSlots == 0);
    CHECK(p.edges.balanced());
    CHECK(p.cells.sequenceDemand == 2 && p.edges.sequenceDemand == 1);
    CHECK(p.cells.remainingSequences > 2);   // fresh pools: vast space
    CHECK(p.cells.predictedGenerations.size() == 2);
    CHECK(p.cells.predictedGenerations[0] == 1);  // in-commit reuse: +1
    CHECK(p.cells.predictedGenerations[1] == 1);
    tx.abandon();
}

static void testCavityPlanPredictedOrder() {
    TEST("topo.cavity.plan_predicted_allocation_order (LIFO of frees)");
    auto f = cav::buildQuad();
    const auto c0slot = f.c[0].slot, c1slot = f.c[1].slot;
    const auto diag = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[2]);
    auto tx = f.ws->beginCavityTransaction();
    (void)tx.stage(cav::swapSpec(f));
    const auto& p = tx.computePlan();
    // Canonical removal order pushes [c0, c1]; LIFO pops predict
    // [c1, c0]. No Phase-2 claims precede them.
    CHECK(p.cells.predictedAllocationOrder.size() == 2);
    CHECK(p.cells.predictedAllocationOrder[0] == c1slot);
    CHECK(p.cells.predictedAllocationOrder[1] == c0slot);
    CHECK(p.edges.predictedAllocationOrder.size() == 1);
    CHECK(p.edges.predictedAllocationOrder[0] == diag.slot);
    tx.abandon();
}

static void testCavityPlanConstAndIdempotent() {
    TEST("topo.cavity.plan_is_const_and_idempotent");
    auto f = cav::buildQuad();
    f.ws->rebuildAdjacency();
    const auto d0 = c1::observableDigest(*f.ws);
    {
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::swapSpec(f));
        const auto& p1 = tx.computePlan();
        const auto& p2 = tx.computePlan();          // cached
        CHECK(&p1 == &p2);
        CHECK(p1.cells.predictedAllocationOrder ==
              p2.cells.predictedAllocationOrder);
    }
    CHECK(c1::observableDigest(*f.ws) == d0);       // const reads only
}

static void testCavityPlanFreeListInterplay() {
    TEST("topo.cavity.plan_free_list_interplay (pre-existing free slot)");
    // Free a cell slot FIRST through the ordinary frozen transaction,
    // then plan a cavity: the balanced demand still prefers in-commit
    // reuse, leaving the pre-existing free slot untouched.
    auto f = cav::buildQuad(1, true);      // + outside tri (n1,n3,n4)
    f.ws->rebuildAdjacency();
    {
        topo::TopologyPolicy orphanOk;     // n4 becomes orphaned
        orphanOk.permitOrphanNodes = true;
        auto tx = f.ws->beginTransaction(orphanOk);
        tx.removeCell(f.c[2]);             // frees the outside tri slot
        (void)tx.commit();
    }
    f.ws->rebuildAdjacency();
    const auto spec = [&] {                 // swap now legal: n1n3 edge
        return cav::swapSpec(f);            // no longer used outside
    }();
    auto tx = f.ws->beginCavityTransaction();
    (void)tx.stage(spec);
    const auto& p = tx.computePlan();
    CHECK(p.cells.freeListAvailable >= 1);  // the freed outside slot
    CHECK(p.cells.inCommitReuse == 2);      // reuse still preferred
    CHECK(p.cells.preExistingConsumed == 0);
    CHECK(p.cells.freshSlots == 0);
    CHECK(p.cells.finalFreeListSize == p.cells.freeListAvailable);
    tx.abandon();
}

static void testCavityPlanGuards() {
    TEST("topo.cavity.plan_guards (staged/rejected/noop)");
    auto f = cav::buildQuad();
    {   // nothing staged
        auto tx = f.ws->beginCavityTransaction();
        bool threw = false;
        try { (void)tx.computePlan(); }
        catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
        tx.abandon();
    }
    {   // rejected by validator
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::swapSpec(f));
        (void)tx.validate([](const topo::CavityPreview&) { return false; });
        bool threw = false;
        try { (void)tx.computePlan(); }
        catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
        tx.abandon();
    }
    {   // semantic no-op: an all-zero plan
        auto tx = f.ws->beginCavityTransaction();
        topo::TriCavitySpec s;
        s.removed = {f.c[0], f.c[1]};
        s.added = {topo::TriCavityCell{{f.n[1], f.n[2], f.n[0]}, {2, 5}},
                   topo::TriCavityCell{{f.n[2], f.n[3], f.n[0]}, {2, 5}}};
        CHECK(tx.stage(s).semanticNoop);
        const auto& p = tx.computePlan();
        CHECK(p.cells.demand == 0 && p.cells.plannedRetirements == 0);
        CHECK(p.edges.demand == 0 && p.edges.plannedRetirements == 0);
        CHECK(p.cells.predictedAllocationOrder.empty());
        tx.abandon();
    }
}

// --- Stage-B preparation round ---------------------------------------------
namespace cav {
// buildQuad with cell c0's slot made GENERATION-EXHAUSTED while live
// (disclosed testing seam), optionally with one pre-existing free slot.
inline Fix buildQuadExhausted(bool withFreeSlot) {
    // Pools are built OUTSIDE the workspace so the disclosed testing
    // seam can exhaust c0's slot generation BEFORE construction (there
    // is no frozen mutable-cells accessor, and none is added).
    Fix f;
    NodePool np;
    f.n.push_back(np.add({0, 0, 0}, {2, 1}));
    f.n.push_back(np.add({1, 0, 0}, {2, 1}));
    f.n.push_back(np.add({1, 1, 0}, {2, 1}));
    f.n.push_back(np.add({0, 1, 0}, {2, 1}));
    if (withFreeSlot) f.n.push_back(np.add({2, .5, 0}, {2, 1}));
    CellPool cp;
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[1], f.n[2]}, {2, 5}, np));
    // Exhaust BEFORE inserting incidence: the incidence entry records
    // the cell's generation, so it must see the exhausted handle.
    f.c[0] = paacs::topo::detail::TransactionAccessTestProbe::
        exhaustLiveCell(cp, f.c[0]);
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[2], f.n[3]}, {2, 5}, np));
    EdgePool ep;
    ep.insertCellEdges(np, cp, f.c[0]);
    ep.insertCellEdges(np, cp, f.c[1]);
    if (withFreeSlot) {
        f.c.push_back(cp.add(ElementType::Tri3,
            std::vector<NodeHandle>{f.n[1], f.n[3], f.n[4]}, {2, 6}, np));
        ep.insertCellEdges(np, cp, f.c.back());
    }
    f.ws = std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), FacePool{});
    if (withFreeSlot) {
        topo::TopologyPolicy orphanOk;
        orphanOk.permitOrphanNodes = true;
        auto tx = f.ws->beginTransaction(orphanOk);
        tx.removeCell(f.c[2]);          // frees the outside-tri slot
        (void)tx.commit();
    }
    return f;
}
}  // namespace cav

static void testCavityPlanExhaustedNoFreeSlot() {
    TEST("topo.cavity.plan_exhausted_retirement_no_free_slot");
    auto f = cav::buildQuadExhausted(false);
    auto tx = f.ws->beginCavityTransaction();
    (void)tx.stage(cav::swapSpec(f));
    const auto& p = tx.computePlan();
    // c0's slot retires on free -> one claim, satisfied FRESH (no free
    // list): freshSlots == exhaustedRetirements here, and the corrected
    // invariant preExistingConsumed + freshSlots == exhausted holds.
    CHECK(p.cells.exhaustedRetirements == 1);
    CHECK(p.cells.reusableThisCommit == 1);
    CHECK(p.cells.inCommitReuse == 1);
    CHECK(p.cells.phase2Claims == 1);
    CHECK(p.cells.freeListAvailable == 0);
    CHECK(p.cells.preExistingConsumed == 0);
    CHECK(p.cells.freshSlots == 1);
    CHECK(p.cells.balanced());
    // order: the fresh claim first, then the in-commit pop of c1's slot
    CHECK(p.cells.predictedAllocationOrder.size() == 2);
    CHECK(p.cells.predictedAllocationOrder[1] == f.c[1].slot);
    CHECK(p.cells.predictedGenerations[0] == 0);      // fresh slot
    tx.abandon();
}

static void testCavityPlanExhaustedWithFreeSlot() {
    TEST("topo.cavity.plan_exhausted_retirement_with_free_slot");
    auto f = cav::buildQuadExhausted(true);
    auto tx = f.ws->beginCavityTransaction();
    (void)tx.stage(cav::swapSpec(f));
    const auto& p = tx.computePlan();
    // Same exhaustion, but a pre-existing free slot satisfies the claim:
    // freshSlots == 0 != exhaustedRetirements == 1 - the reviewer's
    // corrected invariant is the one that holds.
    CHECK(p.cells.exhaustedRetirements == 1);
    CHECK(p.cells.phase2Claims == 1);
    CHECK(p.cells.freeListAvailable == 1);
    CHECK(p.cells.preExistingConsumed == 1);
    CHECK(p.cells.freshSlots == 0);
    CHECK(p.cells.preExistingConsumed + p.cells.freshSlots ==
          p.cells.exhaustedRetirements);
    CHECK(p.cells.balanced());
    CHECK(p.cells.predictedAllocationOrder[0] == f.c[2].slot);  // reused
    tx.abandon();
}

static void testCavitySequenceDemandValidation() {
    TEST("topo.cavity.sequence_demand_validated_before_prepare");
    {   // remaining == demand (0 == 0 at the exhaustion boundary): a
        // classification-only replacement creates NO edges, so a maxed
        // edge sequence space still satisfies the batched demand.
        auto f = cav::buildQuad();
        paacs::topo::detail::TransactionAccessTestProbe::forceEdgeSeqMax(
            f.ws->testEdgesMutable());
        auto tx = f.ws->beginCavityTransaction();
        topo::TriCavitySpec s;
        s.removed = {f.c[0], f.c[1]};
        s.added = {topo::TriCavityCell{{f.n[0], f.n[1], f.n[2]}, {2, 9}},
                   topo::TriCavityCell{{f.n[0], f.n[2], f.n[3]}, {2, 9}}};
        (void)tx.stage(s);
        const auto& p = tx.computePlan();
        CHECK(p.edges.sequenceDemand == 0);
        CHECK(p.edges.remainingSequences == 0);
        tx.prepare();                       // 0 <= 0: succeeds
        CHECK(tx.prepared());
        tx.abandon();
    }
    {   // one short: the swap needs one edge sequence, none remain -
        // prepare throws BEFORE any preparation (full digest equal).
        auto f = cav::buildQuad();
        paacs::topo::detail::TransactionAccessTestProbe::forceEdgeSeqMax(
            f.ws->testEdgesMutable());
        f.ws->rebuildAdjacency();
        const auto d0 = c1::observableDigest(*f.ws);
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::swapSpec(f));
        bool threw = false;
        try { tx.prepare(); }
        catch (const std::logic_error& e) {
            threw = std::string(e.what()).find("sequence") !=
                    std::string::npos;
        }
        CHECK(threw);
        CHECK(!tx.prepared());
        tx.abandon();
        CHECK(c1::observableDigest(*f.ws) == d0);   // NOTHING happened
    }
}

static void testCavityPreparedThenAbandonedLogicalState() {
    TEST("topo.cavity.prepared_then_abandoned_logical_state_restored");
    auto f = cav::buildQuad();
    f.ws->rebuildAdjacency();
    const auto L0 = c1::observableLogicalDigest(*f.ws);
    const auto full0 = c1::observableDigest(*f.ws);
    {
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::swapSpec(f));
        tx.prepare();
        CHECK(tx.prepared());
        tx.prepare();                       // idempotent
        CHECK(tx.prepared());
    }                                        // destruction abandons
    CHECK(f.ws->state() == topo::WorkspaceState::Ready);
    CHECK(c1::observableLogicalDigest(*f.ws) == L0);
    (void)full0;   // the FULL digest may legally differ (capacity residue)
    // clean retry: a fresh transaction reaches Prepared and abandons
    auto tx2 = f.ws->beginCavityTransaction();
    (void)tx2.stage(cav::swapSpec(f));
    tx2.prepare();
    CHECK(tx2.prepared());
    tx2.abandon();
    CHECK(c1::observableLogicalDigest(*f.ws) == L0);
}

static void testCavityPlanAfterMove() {
    TEST("topo.cavity.plan_and_prepare_after_move");
    auto f = cav::buildQuad();
    auto tx = f.ws->beginCavityTransaction();
    (void)tx.stage(cav::swapSpec(f));
    const auto order = tx.computePlan().cells.predictedAllocationOrder;
    auto tx2 = std::move(tx);
    // the captured allocator-state plan transferred with the transaction
    CHECK(tx2.computePlan().cells.predictedAllocationOrder == order);
    tx2.prepare();
    CHECK(tx2.prepared());
    tx2.abandon();
    CHECK(f.ws->state() == topo::WorkspaceState::Ready);
}

static void testCavityPrepareGuards() {
    TEST("topo.cavity.prepare_guards (unstaged/rejected)");
    auto f = cav::buildQuad();
    {
        auto tx = f.ws->beginCavityTransaction();
        bool threw = false;
        try { tx.prepare(); }
        catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
        tx.abandon();
    }
    {
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::swapSpec(f));
        (void)tx.validate([](const topo::CavityPreview&) { return false; });
        bool threw = false;
        try { tx.prepare(); }
        catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
        tx.abandon();
    }
}

// --- Stage-B repair round: state machine, failure sweep, boundaries --------
namespace cav {
// 2x1 strip: two quads (n0,n1,n4,n3), (n1,n2,n5,n4); four triangles.
// reDiagonal re-splits BOTH quads: creates 2 edges, retires 2, keeps
// interior n1n4. exhaustC0 makes cell A's slot generation-exhausted so
// preparation also exercises a fresh slot + fresh incidence entry.
inline Fix buildStrip(bool exhaustC0) {
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
    if (exhaustC0)
        f.c[0] = paacs::topo::detail::TransactionAccessTestProbe::
            exhaustLiveCell(cp, f.c[0]);
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[4], f.n[3]}, {2, 5}, np));
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[1], f.n[2], f.n[5]}, {2, 5}, np));
    f.c.push_back(cp.add(ElementType::Tri3,
        std::vector<NodeHandle>{f.n[1], f.n[5], f.n[4]}, {2, 5}, np));
    EdgePool ep;
    for (const auto& c : f.c) ep.insertCellEdges(np, cp, c);
    f.ws = std::make_unique<TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), FacePool{});
    return f;
}
inline TriCavitySpec reDiagonalSpec(const Fix& f) {
    TriCavitySpec s;
    s.removed = f.c;
    s.added = {TriCavityCell{{f.n[0], f.n[1], f.n[3]}, {2, 5}},
               TriCavityCell{{f.n[1], f.n[4], f.n[3]}, {2, 5}},
               TriCavityCell{{f.n[1], f.n[2], f.n[4]}, {2, 5}},
               TriCavityCell{{f.n[2], f.n[5], f.n[4]}, {2, 5}}};
    return s;
}
}  // namespace cav

static void testCavityPreparedStateMachine() {
    TEST("topo.cavity.prepared_state_machine (exact state; guards)");
    {   // from Ready
        auto f = cav::buildQuad();
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::swapSpec(f));
        tx.prepare();
        CHECK(f.ws->state() == topo::WorkspaceState::TransactionPrepared);
        bool threw = false;                 // reader queries rejected
        try { (void)f.ws->cells(); }
        catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
        threw = false;                      // validators are Phase 1 only
        try {
            (void)tx.validate(
                [](const topo::CavityPreview&) { return true; });
        } catch (const std::logic_error& e) {
            threw = std::string(e.what()).find("Phase 1") !=
                    std::string::npos;
        }
        CHECK(threw);
        tx.abandon();
        CHECK(f.ws->state() == topo::WorkspaceState::Ready);
    }
    {   // from CommittedAdjacencyStale: restored EXACTLY, not to Ready
        auto f = cav::buildQuad(1);
        f.ws->rebuildAdjacency();
        topo::TopologyPolicy orphanOk;
        orphanOk.permitOrphanNodes = true;
        auto t1 = f.ws->beginTransaction(orphanOk);
        t1.addNode({9, 9, 9}, {0, 7});
        (void)t1.commit();
        CHECK(f.ws->state() ==
              topo::WorkspaceState::CommittedAdjacencyStale);
        {
            auto cx = f.ws->beginCavityTransaction();
            (void)cx.stage(cav::swapSpec(f));
            cx.prepare();
            CHECK(f.ws->state() ==
                  topo::WorkspaceState::TransactionPrepared);
        }
        CHECK(f.ws->state() ==
              topo::WorkspaceState::CommittedAdjacencyStale);
    }
}

static void testCavityPhase2FailureSweep() {
    TEST("topo.cavity.phase2_failure_sweep (every internal seam)");
    // The exhausted strip exercises EVERY placement kind: two new-edge
    // table nodes, one fresh incidence entry, record growth, and both
    // allocator capacity classes. Fresh fixture per injected seam k;
    // after each failure: abandon, logical digest equal, byte partition
    // consistent, then retry on the SAME workspace to Prepared.
    int k = 1;
    int failures = 0;
    for (;; ++k) {
        auto f = cav::buildStrip(true);
        f.ws->rebuildAdjacency();
        const auto L0 = c1::observableLogicalDigest(*f.ws);
        const auto spec = cav::reDiagonalSpec(f);
        bool faulted = false;
        {
            auto tx = f.ws->beginCavityTransaction();
            (void)tx.stage(spec);
            paacs::topo::faultCountdown() = k;
            try {
                tx.prepare();
            } catch (const std::bad_alloc&) {
                faulted = true;
            }
            paacs::topo::faultCountdown() = 0;
            if (faulted) {
                CHECK(!tx.prepared());
                tx.abandon();
            } else {
                CHECK(tx.prepared());
                tx.abandon();
            }
        }
        CHECK(c1::observableLogicalDigest(*f.ws) == L0);
        {   // REAL byte-accounting invariants (review: totalBytes is
            // DEFINED as the pool sum, so asserting that identity was
            // vacuous and could never catch broken accounting)
            const auto m = paacs::topo::measureWorkspaceMemory(*f.ws);
            CHECK(m.consistent());
            CHECK(f.ws->nodes().stats().bytesConsistent());
            CHECK(f.ws->cells().stats().bytesConsistent());
            CHECK(f.ws->edges().stats().bytesConsistent());
            CHECK(f.ws->faces().stats().bytesConsistent());
        }
        if (faulted) {
            ++failures;
            // retry on that SAME workspace: reach Prepared, abandon
            auto tx = f.ws->beginCavityTransaction();
            (void)tx.stage(spec);
            tx.prepare();
            CHECK(tx.prepared());
            CHECK(f.ws->state() ==
                  topo::WorkspaceState::TransactionPrepared);
            tx.abandon();
            CHECK(c1::observableLogicalDigest(*f.ws) == L0);
        } else {
            break;                       // k exceeded the seam count
        }
    }
    CHECK(failures >= 5);                // entry, capacity, reserve,
                                         // 2 table, 2+ incidence seams
}

static void testCavitySequenceNonzeroBoundary() {
    TEST("topo.cavity.sequence_nonzero_boundary (demand>0, off-by-one)");
    using Probe = paacs::topo::detail::TransactionAccessTestProbe;
    {   // remaining == demand > 0: accepted (strip: edge demand 2)
        auto f = cav::buildStrip(false);
        Probe::setRemainingEdgeSequences(f.ws->testEdgesMutable(), 2);
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::reDiagonalSpec(f));
        const auto& p = tx.computePlan();
        CHECK(p.edges.sequenceDemand == 2);
        CHECK(p.edges.remainingSequences == 2);
        tx.prepare();
        CHECK(tx.prepared());
        tx.abandon();
    }
    {   // remaining == demand - 1: rejected BEFORE any preparation
        auto f = cav::buildStrip(false);
        Probe::setRemainingEdgeSequences(f.ws->testEdgesMutable(), 1);
        f.ws->rebuildAdjacency();
        const auto d0 = c1::observableDigest(*f.ws);
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::reDiagonalSpec(f));
        bool threw = false;
        try { tx.prepare(); }
        catch (const std::logic_error& e) {
            threw = std::string(e.what()).find("sequence") !=
                    std::string::npos;
        }
        CHECK(threw);
        CHECK(!tx.prepared());
        tx.abandon();
        CHECK(c1::observableDigest(*f.ws) == d0);   // FULL digest equal
    }
}

static void testCavityPreparedMoveTransfersPlaceholders() {
    TEST("topo.cavity.prepared_move_transfers_placeholders");
    // The exhausted strip places BOTH placeholder kinds: two new-edge
    // table nodes and one fresh incidence entry.
    auto f = cav::buildStrip(true);
    f.ws->rebuildAdjacency();
    const auto L0 = c1::observableLogicalDigest(*f.ws);
    {
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::reDiagonalSpec(f));
        tx.prepare();
        CHECK(f.ws->state() == topo::WorkspaceState::TransactionPrepared);
        auto tx2 = std::move(tx);                // move WHILE prepared
        CHECK(f.ws->state() == topo::WorkspaceState::TransactionPrepared);
        CHECK(tx2.prepared());                   // destination prepared
        bool threw = false;                      // source is finished
        try { (void)tx.prepared(), (void)tx.preview(); }
        catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
        tx2.abandon();                           // through the DESTINATION
    }
    CHECK(f.ws->state() == topo::WorkspaceState::Ready);
    CHECK(c1::observableLogicalDigest(*f.ws) == L0);
    // No placeholder survives: the two would-be diagonals resolve to
    // nothing in the edge table after abandonment.
    CHECK(!f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[3]));
    CHECK(!f.ws->edges().find(f.ws->nodes(), f.n[2], f.n[4]));
}

// --- Stage C (partial): semantic no-op commit ------------------------------
static void testCavitySemanticNoopCommit() {
    TEST("topo.cavity.semantic_noop_commit_preserves_state");
    // TWO pre-required cavity edges (01 via the fixture flag, 12 marked
    // below); the spec requests ONLY 01. An implementation that filled
    // requiredEdges from willBeRequired instead of the REQUESTED set
    // would include edge 12 and fail the exclusion assertion.
    auto f = cav::buildQuad(0, false, false, {1, 9}, /*edge01Req=*/true);
    f.ws->rebuildAdjacency();
    const auto e01 = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[1]);
    const auto rev0 = f.ws->topologyRevision();
    topo::TriCavitySpec s;
    s.removed = {f.c[0], f.c[1]};
    s.added = {topo::TriCavityCell{{f.n[1], f.n[2], f.n[0]}, {2, 5}},
               topo::TriCavityCell{{f.n[2], f.n[3], f.n[0]}, {2, 5}}};
    s.requirements = {topo::CavityEdgeRequirement{f.n[0], f.n[1], {1, 9}}};
    // second pre-required edge, NOT requested (must be excluded)
    const auto e12 = *f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[2]);
    f.ws->testEdgesMutable().setRequiredConstraint(e12, true);
    const auto full0 = c1::observableDigest(*f.ws);
    auto tx = f.ws->beginCavityTransaction();
    const auto& p = tx.stage(s);
    CHECK(p.semanticNoop);
    const auto out = tx.commit();
    CHECK(f.ws->state() == topo::WorkspaceState::Ready);   // EXACT prior
    CHECK(out.revision == rev0);
    CHECK(f.ws->topologyRevision() == rev0);                        // no churn
    CHECK(out.resultCells.size() == 2);                     // EXISTING
    CHECK(out.resultCells[0] == f.c[0]);                    // handles,
    CHECK(out.resultCells[1] == f.c[1]);                    // canonical
    CHECK(f.ws->cells().isCurrent(f.c[0]));
    CHECK(f.ws->cells().isCurrent(f.c[1]));
    CHECK(out.requiredEdges.size() == 1);   // EXACTLY the requested set:
    CHECK(out.requiredEdges[0] == e01);     // e12 is required-but-not-
    for (const auto& h : out.requiredEdges) //   -requested -> EXCLUDED
        CHECK(!(h == e12));
    CHECK(c1::observableDigest(*f.ws) == full0);            // FULL digest:
                                                            // zero churn
    bool threw = false;                                     // finished
    try { (void)tx.preview(); }
    catch (const std::logic_error&) { threw = true; }
    CHECK(threw);
}

static void testCavitySemanticNoopCommitFromStale() {
    TEST("topo.cavity.semantic_noop_commit_restores_stale_state");
    auto f = cav::buildQuad(1);
    f.ws->rebuildAdjacency();
    topo::TopologyPolicy orphanOk;
    orphanOk.permitOrphanNodes = true;
    auto t1 = f.ws->beginTransaction(orphanOk);
    t1.addNode({9, 9, 9}, {0, 7});
    (void)t1.commit();
    CHECK(f.ws->state() == topo::WorkspaceState::CommittedAdjacencyStale);
    topo::TriCavitySpec s;
    s.removed = {f.c[0], f.c[1]};
    s.added = {topo::TriCavityCell{{f.n[0], f.n[1], f.n[2]}, {2, 5}},
               topo::TriCavityCell{{f.n[0], f.n[2], f.n[3]}, {2, 5}}};
    auto tx = f.ws->beginCavityTransaction();
    CHECK(tx.stage(s).semanticNoop);
    (void)tx.commit();
    // "Preserve prior workspace state": STALE, not blindly Ready.
    CHECK(f.ws->state() == topo::WorkspaceState::CommittedAdjacencyStale);
}

static void testCavityPreparedNoopCommit() {
    TEST("topo.cavity.semantic_noop_commit_after_prepare");
    // prepare() -> TransactionPrepared -> commit(): the no-op path must
    // consume the Stage-B prepared state (an all-zero plan here, but
    // the PATH through prepared cleanup is what this exercises) and
    // restore the exact prior state and FULL digest.
    auto f = cav::buildQuad();
    f.ws->rebuildAdjacency();
    const auto rev0 = f.ws->topologyRevision();
    const auto full0 = c1::observableDigest(*f.ws);
    topo::TriCavitySpec s;
    s.removed = {f.c[0], f.c[1]};
    s.added = {topo::TriCavityCell{{f.n[1], f.n[2], f.n[0]}, {2, 5}},
               topo::TriCavityCell{{f.n[2], f.n[3], f.n[0]}, {2, 5}}};
    auto tx = f.ws->beginCavityTransaction();
    CHECK(tx.stage(s).semanticNoop);
    tx.prepare();
    CHECK(f.ws->state() == topo::WorkspaceState::TransactionPrepared);
    const auto out = tx.commit();
    CHECK(f.ws->state() == topo::WorkspaceState::Ready);
    CHECK(out.revision == rev0);
    CHECK(f.ws->topologyRevision() == rev0);
    CHECK(c1::observableDigest(*f.ws) == full0);
}

static void testCavitySubstantiveCommitSwap() {
    TEST("topo.cavity.substantive_commit (swap: net delta, handles, rev)");
    // Edge 01 pre-marked required {1,9}: outer-boundary constraint
    // metadata must survive with the SAME handle.
    auto f = cav::buildQuad(0, false, false, {1, 9}, true);
    f.ws->rebuildAdjacency();
    const auto rev0 = f.ws->topologyRevision();
    const auto full0 = c1::observableDigest(*f.ws);
    const auto e01 = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[1]);
    const auto e12 = *f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[2]);
    const auto e23 = *f.ws->edges().find(f.ws->nodes(), f.n[2], f.n[3]);
    const auto e03 = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[3]);
    const auto eDiag = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[2]);
    auto tx = f.ws->beginCavityTransaction();
    auto s = cav::swapSpec(f);
    s.requirements = {topo::CavityEdgeRequirement{f.n[1], f.n[3], {1, 7}}};
    (void)tx.stage(s);
    const auto planCells =
        tx.computePlan().cells.predictedAllocationOrder;
    const auto planCellGens = tx.computePlan().cells.predictedGenerations;
    const auto planEdges = tx.computePlan().edges.predictedAllocationOrder;
    const auto out = tx.commit();
    // revision exactly once; adjacency stale
    CHECK(f.ws->state() == topo::WorkspaceState::CommittedAdjacencyStale);
    CHECK(out.revision == rev0 + 1);
    CHECK(f.ws->topologyRevision() == rev0 + 1);
    CHECK(f.ws->adjacencyBuiltFromRevision() == rev0); // stale, valid
    CHECK(!f.ws->adjacencyReady());
    // predicted order AND generations agreed (also asserted in-helper)
    CHECK(out.resultCells.size() == 2);
    for (std::size_t j = 0; j < 2; ++j) {
        CHECK(out.resultCells[j].slot == planCells[j]);
        CHECK(out.resultCells[j].gen == planCellGens[j]);
    }
    // removed handles stale; result cells current with swapped conn
    CHECK(!f.ws->cells().isCurrent(f.c[0]));
    CHECK(!f.ws->cells().isCurrent(f.c[1]));
    CHECK(f.ws->cells().isCurrent(out.resultCells[0]));
    CHECK(f.ws->cells().isCurrent(out.resultCells[1]));
    // old diagonal retired; new diagonal created on the predicted slot
    CHECK(!f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[2]));
    const auto eNew = f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[3]);
    CHECK(bool(eNew));
    CHECK(eNew->slot == planEdges[0]);
    CHECK(!(eNew->slot == eDiag.slot) || eNew->gen != eDiag.gen);
    // boundary edge HANDLES retained; net use counts exact
    for (const auto& e : {e01, e12, e23, e03}) {
        CHECK(f.ws->edges().isCurrent(e));
        CHECK(f.ws->edges().get(e).useCount == 1);
    }
    CHECK(f.ws->edges().get(*eNew).useCount == 2);
    // outer-boundary constraint metadata SURVIVED on the same handle
    CHECK(f.ws->edges().get(e01).requiredConstraint);
    CHECK(f.ws->edges().get(e01).cls.dim == 1);
    CHECK(f.ws->edges().get(e01).cls.tag == 9);
    // requirement applied IN THE SAME COMMIT
    CHECK(out.requiredEdges.size() == 1);
    CHECK(out.requiredEdges[0] == *eNew);
    CHECK(f.ws->edges().get(*eNew).requiredConstraint);
    CHECK(f.ws->edges().get(*eNew).cls.dim == 1);
    CHECK(f.ws->edges().get(*eNew).cls.tag == 7);
    // substantive commit CHANGES the digest; topology valid while
    // adjacency stale: rebuild succeeds and audits pass afterwards
    CHECK(!(c1::observableDigest(*f.ws) == full0));
    f.ws->rebuildAdjacency();
    CHECK(f.ws->state() == topo::WorkspaceState::Ready);
}

static void testCavitySubstantiveCommitOutsideSnapshot() {
    TEST("topo.cavity.substantive_commit_outside_exactly_unchanged");
    // Outside triangle (n2,n3,n4) shares boundary edge n2n3 with the
    // cavity; its handle, connectivity, cls, other edges and incidence
    // must be EXACTLY unchanged, and n2n3 must keep its handle.
    cav::Fix f;
    {
        topo::NodePool np;
        for (int i = 0; i < 5; ++i)
            f.n.push_back(np.add({double(i % 3), double(i / 3 + i % 2), 0},
                                 {2, 1}));
        topo::CellPool cp;
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{f.n[0], f.n[1], f.n[2]},
            {2, 5}, np));
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{f.n[0], f.n[2], f.n[3]},
            {2, 5}, np));
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{f.n[3], f.n[2], f.n[4]},
            {2, 6}, np));
        topo::EdgePool ep;
        for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
        f.ws = std::make_unique<topo::TopologyWorkspace>(
            std::move(np), std::move(cp), std::move(ep),
            topo::FacePool{});
    }
    const auto outside = f.c[2];
    const auto e23 = *f.ws->edges().find(f.ws->nodes(), f.n[2], f.n[3]);
    const auto e34 = *f.ws->edges().find(f.ws->nodes(), f.n[3], f.n[4]);
    const auto e24 = *f.ws->edges().find(f.ws->nodes(), f.n[2], f.n[4]);
    const auto snapConn = f.ws->cells().get(outside);
    auto tx = f.ws->beginCavityTransaction();
    topo::TriCavitySpec s;                       // swap inside the quad
    s.removed = {f.c[0], f.c[1]};
    s.added = {topo::TriCavityCell{{f.n[0], f.n[1], f.n[3]}, {2, 5}},
               topo::TriCavityCell{{f.n[1], f.n[2], f.n[3]}, {2, 5}}};
    (void)tx.stage(s);
    (void)tx.commit();
    // outside cell: same handle, same payload
    CHECK(f.ws->cells().isCurrent(outside));
    const auto after = f.ws->cells().get(outside);
    CHECK(after.type == snapConn.type);
    for (int i = 0; i < 3; ++i)
        CHECK(after.nodes[std::size_t(i)] ==
              snapConn.nodes[std::size_t(i)]);
    CHECK(after.cls.dim == snapConn.cls.dim &&
          after.cls.tag == snapConn.cls.tag);
    // shared boundary edge: SAME handle, use back to 2
    CHECK(f.ws->edges().isCurrent(e23));
    CHECK(f.ws->edges().get(e23).useCount == 2);
    // untouched outside edges: same handles, unchanged counts
    CHECK(f.ws->edges().isCurrent(e34));
    CHECK(f.ws->edges().isCurrent(e24));
    CHECK(f.ws->edges().get(e34).useCount == 1);
    CHECK(f.ws->edges().get(e24).useCount == 1);
    f.ws->rebuildAdjacency();                    // audits pass
}

static void testCavitySubstantiveCommitExhausted() {
    TEST("topo.cavity.substantive_commit_exhausted_fixtures");
    {   // NO pre-existing free slot: the exhausted claim materializes
        // FRESH (gen 0); the exhausted slot is never reused.
        auto f = cav::buildQuadExhausted(false);
        const auto exhaustedSlot = f.c[0].slot;
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::swapSpec(f));
        const auto plan = tx.computePlan().cells;
        const auto out = tx.commit();
        CHECK(out.resultCells[0].slot == plan.predictedAllocationOrder[0]);
        CHECK(out.resultCells[0].gen == 0);            // FRESH slot
        for (const auto& ch : out.resultCells)
            CHECK(!(ch.slot == exhaustedSlot));        // retired forever
        f.ws->rebuildAdjacency();
    }
    {   // WITH a pre-existing free slot: it satisfies the claim; zero
        // fresh materialization.
        auto f = cav::buildQuadExhausted(true);
        const auto freedSlot = f.c[2].slot;
        auto tx = f.ws->beginCavityTransaction();
        (void)tx.stage(cav::swapSpec(f));
        const auto matBefore = tx.computePlan().cells.freshSlots;
        CHECK(matBefore == 0);
        const auto out = tx.commit();
        CHECK(out.resultCells[0].slot == freedSlot);   // reused free slot
        f.ws->rebuildAdjacency();
    }
}

static void testCavitySubstantiveThenNoopAndRetry() {
    TEST("topo.cavity.substantive_then_noop (V2 changes; noop after)");
    auto f = cav::buildQuad();
    f.ws->rebuildAdjacency();
    auto tx = f.ws->beginCavityTransaction();
    (void)tx.stage(cav::swapSpec(f));
    const auto out = tx.commit();
    f.ws->rebuildAdjacency();
    const auto full1 = c1::observableDigest(*f.ws);
    const auto rev1 = f.ws->topologyRevision();
    // Re-adding the SAME swapped triangles is now a semantic no-op.
    topo::TriCavitySpec s2;
    s2.removed = {out.resultCells[0], out.resultCells[1]};
    s2.added = {topo::TriCavityCell{{f.n[0], f.n[1], f.n[3]}, {2, 5}},
                topo::TriCavityCell{{f.n[1], f.n[2], f.n[3]}, {2, 5}}};
    auto tx2 = f.ws->beginCavityTransaction();
    CHECK(tx2.stage(s2).semanticNoop);
    const auto out2 = tx2.commit();
    CHECK(out2.revision == rev1);
    CHECK(c1::observableDigest(*f.ws) == full1);
    CHECK(f.ws->state() == topo::WorkspaceState::Ready);
}

// --- Stage-C evidence round: wider fixtures, V2 direct, exact snapshot ----
static void testCavitySubstantiveCommitStrip() {
    TEST("topo.cavity.substantive_commit_strip (4 cells, 2 create/retire)");
    // The exhausted strip commits with EVERYTHING at once: two table
    // placeholders, two retirements, a retained interior edge (n1n4),
    // one fresh slot + fresh incidence (exhausted c0), in-commit reuse,
    // and multiple incidence rewrites.
    auto f = cav::buildStrip(true);
    f.ws->rebuildAdjacency();
    const auto rev0 = f.ws->topologyRevision();
    const auto e14 = *f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[4]);
    auto tx = f.ws->beginCavityTransaction();
    (void)tx.stage(cav::reDiagonalSpec(f));
    const auto plan = tx.computePlan();
    const auto out = tx.commit();
    CHECK(out.resultCells.size() == 4);
    CHECK(f.ws->topologyRevision() == rev0 + 1);
    for (std::size_t j = 0; j < 4; ++j) {
        CHECK(out.resultCells[j].slot ==
              plan.cells.predictedAllocationOrder[j]);
        CHECK(out.resultCells[j].gen ==
              plan.cells.predictedGenerations[j]);
        CHECK(f.ws->cells().isCurrent(out.resultCells[j]));
    }
    // both old diagonals retired; both new ones exist with use 2
    CHECK(!f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[4]));
    CHECK(!f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[5]));
    const auto e13 = f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[3]);
    const auto e24 = f.ws->edges().find(f.ws->nodes(), f.n[2], f.n[4]);
    CHECK(bool(e13) && bool(e24));
    CHECK(f.ws->edges().get(*e13).useCount == 2);
    CHECK(f.ws->edges().get(*e24).useCount == 2);
    // the retained interior edge kept its HANDLE with use 2
    CHECK(f.ws->edges().isCurrent(e14));
    CHECK(f.ws->edges().get(e14).useCount == 2);
    f.ws->rebuildAdjacency();
    CHECK(f.ws->adjacencyReady());
}

static void testCavitySubstantiveV2Direct() {
    TEST("topo.cavity.substantive_replacement_changes_v2 (direct)");
    auto f = cav::buildQuad();
    f.ws->rebuildAdjacency();
    const auto v2Before = f.ws->semanticTopologyHashV2();
    const auto tokensBefore = f.ws->testV2Tokens();
    auto tx = f.ws->beginCavityTransaction();
    (void)tx.stage(cav::swapSpec(f));
    (void)tx.commit();
    // Direct semantic capture: revision/state changes cannot mask an
    // accidentally-unchanged V2 (the earlier full-digest assertion was
    // vacuous for this specific claim).
    CHECK(f.ws->semanticTopologyHashV2() != v2Before);
    CHECK(!(f.ws->testV2Tokens() == tokensBefore));
    f.ws->rebuildAdjacency();
}

static void testCavitySemanticDeterminismModuloOwnership() {
    TEST("topo.cavity.semantic_determinism_modulo_ownership");
    // FRESH quad vs an allocator-CHURNED quad (an extra node+triangle
    // committed and removed, shuffling slots, generations, sequences
    // and revisions) -> commit the SAME swap on both -> the semantic V2
    // stream and hash agree, while handles/strict history/revisions are
    // explicitly NOT required to (and do not) agree.
    auto fresh = cav::buildQuad();
    fresh.ws->rebuildAdjacency();
    auto churn = cav::buildQuad(1, true);   // + outside tri on n1n3,n4
    churn.ws->rebuildAdjacency();
    {
        topo::TopologyPolicy orphanOk;
        orphanOk.permitOrphanNodes = true;
        auto txc = churn.ws->beginTransaction(orphanOk);
        txc.removeCell(churn.c[2]);         // back to the plain quad,
        (void)txc.commit();                 // with different history
        churn.ws->rebuildAdjacency();
        auto txn = churn.ws->beginTransaction(orphanOk);
        txn.removeNode(churn.n[4]);
        (void)txn.commit();
        churn.ws->rebuildAdjacency();
    }
    CHECK(churn.ws->semanticTopologyHashV2() ==
          fresh.ws->semanticTopologyHashV2());      // equal pre-state
    auto t1 = fresh.ws->beginCavityTransaction();
    (void)t1.stage(cav::swapSpec(fresh));
    const auto o1 = t1.commit();
    auto t2 = churn.ws->beginCavityTransaction();
    (void)t2.stage(cav::swapSpec(churn));
    const auto o2 = t2.commit();
    CHECK(fresh.ws->semanticTopologyHashV2() ==
          churn.ws->semanticTopologyHashV2());
    CHECK(fresh.ws->testV2Tokens() == churn.ws->testV2Tokens());
    // modulo ownership and allocator history: these MAY differ, and do
    CHECK(o1.revision != o2.revision);
    CHECK(fresh.ws->strictWorkspaceHash() !=
          churn.ws->strictWorkspaceHash());
}

static void testCavityOutsideSnapshotExact() {
    TEST("topo.cavity.outside_snapshot_exact (ordered incidence+payload)");
    // Field-for-field: the outside cell's ORDERED cellEdges() relation
    // and every one of its edges' complete payloads (endpoints, cls,
    // required flag, use count) must be byte-identical after the commit
    // except the shared boundary edge's use count, which returns to 2.
    cav::Fix f;
    {
        topo::NodePool np;
        for (int i = 0; i < 5; ++i)
            f.n.push_back(np.add({double(i % 3), double(i / 3 + i % 2), 0},
                                 {2, 1}));
        topo::CellPool cp;
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{f.n[0], f.n[1], f.n[2]},
            {2, 5}, np));
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{f.n[0], f.n[2], f.n[3]},
            {2, 5}, np));
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{f.n[3], f.n[2], f.n[4]},
            {2, 6}, np));
        topo::EdgePool ep;
        for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
        ep.setRequiredConstraint(*ep.find(np, f.n[3], f.n[4]), true);
        ep.setClassification(*ep.find(np, f.n[3], f.n[4]), {1, 4});
        f.ws = std::make_unique<topo::TopologyWorkspace>(
            std::move(np), std::move(cp), std::move(ep),
            topo::FacePool{});
    }
    const auto outside = f.c[2];
    struct EdgeSnap {
        topo::EdgeHandle h;
        topo::NodeHandle a, b;
        paacs::mesh::Classification cls;
        bool req;
        std::uint32_t use;
    };
    std::vector<EdgeSnap> snap;
    {
        const auto ce = f.ws->edges().cellEdges(f.ws->cells(), outside);
        for (const auto& eh : ce) {
            const auto d = f.ws->edges().get(eh);
            snap.push_back(
                {eh, d.a, d.b, d.cls, d.requiredConstraint, d.useCount});
        }
    }
    CHECK(snap.size() == 3);
    auto tx = f.ws->beginCavityTransaction();
    topo::TriCavitySpec s;
    s.removed = {f.c[0], f.c[1]};
    s.added = {topo::TriCavityCell{{f.n[0], f.n[1], f.n[3]}, {2, 5}},
               topo::TriCavityCell{{f.n[1], f.n[2], f.n[3]}, {2, 5}}};
    (void)tx.stage(s);
    (void)tx.commit();
    const auto ce2 = f.ws->edges().cellEdges(f.ws->cells(), outside);
    CHECK(ce2.size() == 3);
    const auto e23 = *f.ws->edges().find(f.ws->nodes(), f.n[2], f.n[3]);
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(ce2[i] == snap[i].h);              // ORDERED incidence
        const auto d = f.ws->edges().get(ce2[i]);
        CHECK(d.a == snap[i].a);
        CHECK(d.b == snap[i].b);
        CHECK(d.cls.dim == snap[i].cls.dim &&
              d.cls.tag == snap[i].cls.tag);
        CHECK(d.requiredConstraint == snap[i].req);
        if (ce2[i] == e23)
            CHECK(d.useCount == 2);              // authorized: back to 2
        else
            CHECK(d.useCount == snap[i].use);
    }
    f.ws->rebuildAdjacency();
}

// Registration hook consumed by test_main.cpp.
namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* transactionTests(std::size_t& count) {
    static const Entry kTx[] = {
        {"topo.transaction.foundations",  testTransactionFoundations},
        {"topo.transaction.state",        testTransactionStateMachine},
        {"topo.transaction.cellops",      testTransactionCellOps},
        {"topo.transaction.replace",      testTransactionReplaceAndPolicy},
        {"topo.transaction.metadata",     testTransactionMetadata},
        {"topo.transaction.nodes",        testTransactionNodeLifecycle},
        {"topo.hash.v2",                  testHashV2},
        {"topo.hash.v2stream",            testHashV2StreamContract},
        {"topo.compaction.semantic_hash_unchanged",
         testCompactionSemanticHashUnchanged},
        {"topo.compaction.node_map_total", testCompactionNodeMapTotal},
        {"topo.compaction.cell_map_total", testCompactionCellMapTotal},
        {"topo.compaction.unmapped_handle_classes",
         testCompactionUnmappedHandleClasses},
        {"topo.compaction.source_old_handles_remain_current",
         testCompactionSourceOldHandlesRemainCurrent},
        {"topo.compaction.old_handles_foreign_to_compacted_workspace",
         testCompactionOldHandlesForeign},
        {"topo.compaction.mapped_handles_current",
         testCompactionMappedHandlesCurrent},
        {"topo.compaction.connectivity_remapped",
         testCompactionConnectivityRemapped},
        {"topo.compaction.metadata_reconstructed",
         testCompactionMetadataReconstructed},
        {"topo.compaction.workspace_ready", testCompactionWorkspaceReady},
        {"topo.compaction.adjacency_rebuilt",
         testCompactionAdjacencyRebuilt},
        {"topo.compaction.strict_hash_matches_deterministic_rebuild",
         testCompactionStrictHashDeterministicRebuild},
        {"topo.compaction.deterministic_construction",
         testCompactionDeterministicConstruction},
        {"topo.compaction.memory_report_classes",
         testCompactionMemoryReportClasses},
        {"topo.compaction.source_success_observably_unchanged",
         testCompactionSourceSuccessObservablyUnchanged},
        {"topo.compaction.source_failure_observably_unchanged",
         testCompactionSourceFailureObservablyUnchanged},
        {"topo.compaction.edge_map_total", testCompactionEdgeMapTotal},
        {"topo.compaction.face_map_total", testCompactionFaceMapTotal},
        {"topo.compaction.free_handles_unmapped",
         testCompactionFreeHandlesUnmapped},
        {"topo.compaction.retired_handles_unmapped",
         testCompactionRetiredHandlesUnmapped},
        {"topo.compaction.stale_generations_unmapped",
         testCompactionStaleGenerationsUnmapped},
        {"topo.compaction.incidence_remapped",
         testCompactionIncidenceRemapped},
        {"topo.compaction.metadata_preserved",
         testCompactionMetadataPreserved},
        {"topo.compaction.deterministic_after_churn",
         testCompactionDeterministicAfterChurn},
        {"topo.compaction.memory_reduced", testCompactionMemoryReduced},
        {"topo.compaction.no_benefit_reported",
         testCompactionNoBenefitReported},
        {"topo.compaction.option_domain_validated",
         testCompactionOptionDomain},
        {"topo.cavity.state_machine", testCavityStateMachine},
        {"topo.cavity.stage_valid_preview", testCavityStageValidPreview},
        {"topo.cavity.canonicalization_ignores_caller_order",
         testCavityCanonicalizationIgnoresCallerOrder},
        {"topo.cavity.validator_rejection", testCavityValidatorRejection},
        {"topo.cavity.phase1_full_digest_unchanged",
         testCavityPhase1DigestUnchanged},
        {"topo.cavity.disk_invariant_enforced",
         testCavityDiskInvariantEnforced},
        {"topo.cavity.boundary_cycle_is_simple",
         testCavityAnnulusRejected},
        {"topo.cavity.source_nonmanifold_cavity_edge_rejected",
         testCavitySourceNonmanifoldEdgeRejected},
        {"topo.cavity.candidate_interior_edge_with_outside_use_rejected",
         testCavityInteriorEdgeOutsideUseRejected},
        {"topo.cavity.orphan_node_rejected", testCavityOrphanNodeRejected},
        {"topo.cavity.outside_node_pulled_in_rejected",
         testCavityOutsideNodePulledInRejected},
        {"topo.cavity.duplicate_or_stale_removal_rejected",
         testCavityDuplicateOrStaleRemovalRejected},
        {"topo.cavity.boundary_violation_rejected",
         testCavityBoundaryViolationRejected},
        {"topo.cavity.inside_sense_or_orientation_rejected",
         testCavityInsideSenseMismatchRejected},
        {"topo.cavity.requirement_rules", testCavityRequirementRules},
        {"topo.cavity.required_edge_may_not_be_retired",
         testCavityRequiredEdgeMayNotBeRetired},
        {"topo.cavity.semantic_noop_detected",
         testCavitySemanticNoopDetected},
        {"topo.cavity.boundary_cycle_nonempty",
         testCavityBoundaryCycleNonempty},
        {"topo.cavity.node_link_manifold_enforced",
         testCavityNodeLinkManifoldEnforced},
        {"topo.cavity.euler_characteristic_enforced",
         testCavityEulerCharacteristicEnforced},
        {"topo.cavity.surviving_duplicate_rejected",
         testCavitySurvivingDuplicateRejected},
        {"topo.cavity.composite_gate_enforced",
         testCavityCompositeGateEnforced},
        {"topo.cavity.state_machine_extended",
         testCavityStateMachineExtended},
        {"topo.cavity.requirement_outside_closure_rejected",
         testCavityRequirementOutsideClosureRejected},
        {"topo.cavity.plan_slot_reuse_exact", testCavityPlanSlotReuseExact},
        {"topo.cavity.plan_predicted_allocation_order",
         testCavityPlanPredictedOrder},
        {"topo.cavity.plan_is_const_and_idempotent",
         testCavityPlanConstAndIdempotent},
        {"topo.cavity.plan_free_list_interplay",
         testCavityPlanFreeListInterplay},
        {"topo.cavity.plan_guards", testCavityPlanGuards},
        {"topo.cavity.plan_exhausted_retirement_no_free_slot",
         testCavityPlanExhaustedNoFreeSlot},
        {"topo.cavity.plan_exhausted_retirement_with_free_slot",
         testCavityPlanExhaustedWithFreeSlot},
        {"topo.cavity.sequence_demand_validated_before_prepare",
         testCavitySequenceDemandValidation},
        {"topo.cavity.prepared_then_abandoned_logical_state_restored",
         testCavityPreparedThenAbandonedLogicalState},
        {"topo.cavity.plan_and_prepare_after_move", testCavityPlanAfterMove},
        {"topo.cavity.prepare_guards", testCavityPrepareGuards},
        {"topo.cavity.prepared_state_machine",
         testCavityPreparedStateMachine},
        {"topo.cavity.phase2_failure_sweep", testCavityPhase2FailureSweep},
        {"topo.cavity.sequence_nonzero_boundary",
         testCavitySequenceNonzeroBoundary},
        {"topo.cavity.prepared_move_transfers_placeholders",
         testCavityPreparedMoveTransfersPlaceholders},
        {"topo.cavity.semantic_noop_commit_preserves_state",
         testCavitySemanticNoopCommit},
        {"topo.cavity.semantic_noop_commit_restores_stale_state",
         testCavitySemanticNoopCommitFromStale},
        {"topo.cavity.semantic_noop_commit_after_prepare",
         testCavityPreparedNoopCommit},
        {"topo.cavity.substantive_commit_swap",
         testCavitySubstantiveCommitSwap},
        {"topo.cavity.substantive_commit_outside_exactly_unchanged",
         testCavitySubstantiveCommitOutsideSnapshot},
        {"topo.cavity.substantive_commit_exhausted_fixtures",
         testCavitySubstantiveCommitExhausted},
        {"topo.cavity.substantive_then_noop",
         testCavitySubstantiveThenNoopAndRetry},
        {"topo.cavity.substantive_commit_strip",
         testCavitySubstantiveCommitStrip},
        {"topo.cavity.substantive_replacement_changes_v2",
         testCavitySubstantiveV2Direct},
        {"topo.cavity.semantic_determinism_modulo_ownership",
         testCavitySemanticDeterminismModuloOwnership},
        {"topo.cavity.outside_snapshot_exact",
         testCavityOutsideSnapshotExact},
    };
    count = sizeof(kTx) / sizeof(kTx[0]);
    return kTx;
}
} // namespace paacs_tests
