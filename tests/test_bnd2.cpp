// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// BND-2A / ASCR (Atomic Steiner Cavity Replacement) proofs. Design authority:
// docs/PMBND002_DESIGN.md rev 3. Registered in the main test binary via
// bnd2Tests().
//
// PART 1 - additive rev-3 combinatorial helpers (paacs::topo::cavaudit).
// PART 2 - SteinerCavityTransaction end-to-end: single- AND multi-Steiner
// atomic replacement, staged-staged edges + requirements, container-order
// independence, ownership-stamp validation, exact node-sequence boundaries,
// mixed reusable/fresh and retired node slots, full outside-cavity snapshot,
// the prepare-poison contract, and the adversarial Phase-1 rejections.
#include <cstdio>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <utility>
#include <memory>
#include <new>
#include <stdexcept>
#include <vector>

#include "paacs/topo/SteinerCavity.hpp"
#include "paacs/bnd/SteinerRecovery.hpp"
#include "workspace_digest.hpp"

namespace topo = paacs::topo;
namespace cav = paacs::topo::cavaudit;
namespace mesh = paacs::mesh;
namespace bnd = paacs::bnd;

extern int gFailures;
static const char* gCurrent = "";
#define TEST(name)                          \
    do {                                    \
        gCurrent = name;                    \
        std::printf("[test] %s\n", name);   \
    } while (0)
#define CHECK(cond)                                                        \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::printf("  FAIL %s:%d  %s  (in %s)\n", __FILE__, __LINE__, \
                        #cond, gCurrent);                                  \
            ++gFailures;                                                   \
        }                                                                  \
    } while (0)

namespace {

using cav::RankTri;
using cav::Rank;
using topo::NodeHandle;
using topo::CellHandle;
using topo::EdgeHandle;
using topo::SteinerNodeRef;

struct OutsideTriangleSnapshot {
    CellHandle cell;
    topo::ElementType type;
    mesh::Classification cellClass;
    std::array<NodeHandle, 3> nodes;
    std::array<topo::NodePool::NodeData, 3> nodeData;
    std::array<EdgeHandle, 3> edges;
    std::array<topo::EdgePool::EdgeData, 3> edgeData;
};

OutsideTriangleSnapshot snapshotOutsideTriangle(
    const topo::TopologyWorkspace& ws, CellHandle cell) {
    OutsideTriangleSnapshot out{};
    out.cell = cell;
    const auto cs = ws.cells().get(cell);
    out.type = cs.type;
    out.cellClass = cs.cls;
    for (std::size_t i = 0; i < 3; ++i) {
        out.nodes[i] = cs.nodes[i];
        out.nodeData[i] = ws.nodes().get(out.nodes[i]);
        const auto j = (i + 1) % 3;
        const auto edge = ws.edges().find(ws.nodes(), cs.nodes[i], cs.nodes[j]);
        CHECK(edge.has_value());
        if (edge) {
            out.edges[i] = *edge;
            out.edgeData[i] = ws.edges().get(*edge);
        }
    }
    return out;
}

bool outsideTriangleMatches(const topo::TopologyWorkspace& ws,
                            const OutsideTriangleSnapshot& before) {
    if (!ws.cells().isCurrent(before.cell)) return false;
    const auto cs = ws.cells().get(before.cell);
    if (cs.type != before.type ||
        cs.cls.dim != before.cellClass.dim ||
        cs.cls.tag != before.cellClass.tag)
        return false;
    for (std::size_t i = 0; i < 3; ++i) {
        if (!(cs.nodes[i] == before.nodes[i]) ||
            !ws.nodes().isCurrent(before.nodes[i]))
            return false;
        const auto nd = ws.nodes().get(before.nodes[i]);
        const auto& oldNd = before.nodeData[i];
        if (nd.coord.x != oldNd.coord.x || nd.coord.y != oldNd.coord.y ||
            nd.coord.z != oldNd.coord.z || nd.cls.dim != oldNd.cls.dim ||
            nd.cls.tag != oldNd.cls.tag)
            return false;
        const auto j = (i + 1) % 3;
        const auto edge =
            ws.edges().find(ws.nodes(), before.nodes[i], before.nodes[j]);
        if (!edge || !(*edge == before.edges[i])) return false;
        const auto ed = ws.edges().get(*edge);
        const auto& oldEd = before.edgeData[i];
        if (!(ed.a == oldEd.a) || !(ed.b == oldEd.b) ||
            ed.useCount != oldEd.useCount ||
            ed.requiredConstraint != oldEd.requiredConstraint ||
            ed.cls.dim != oldEd.cls.dim || ed.cls.tag != oldEd.cls.tag)
            return false;
    }
    return true;
}

// ===================================================================== PART 1
void testCavAuditDiskAndBoundary() {
    TEST("bnd2.cavaudit_disk_boundary (single tri, disk, staged fan, |M|!=|R|)");
    {
        std::vector<RankTri> t{{1, 2, 3}};
        CHECK(cav::auditCavityDiskCombinatorial(t));
        CHECK(cav::auditCavityNodeLink(t));
        std::vector<Rank> cyc;
        CHECK(cav::boundaryCycle(t, cyc));
        CHECK(cyc.size() == 3);
        CHECK(cyc[0] == 1);
    }
    std::vector<RankTri> R{{1, 2, 3}, {1, 3, 4}};
    CHECK(cav::auditCavityDiskCombinatorial(R));
    CHECK(cav::auditCavityNodeLink(R));
    const auto bR = cav::boundaryEdgeSet(R);
    CHECK(bR.size() == 4);
    {
        std::vector<Rank> cyc;
        CHECK(cav::boundaryCycle(R, cyc));
        const std::vector<Rank> expect{1, 2, 3, 4};
        CHECK(cyc == expect);
    }
    std::vector<RankTri> M{{1, 2, 100}, {2, 3, 100}, {3, 4, 100}, {4, 1, 100}};
    CHECK(cav::auditCavityDiskCombinatorial(M));
    CHECK(cav::auditCavityNodeLink(M));
    CHECK(cav::boundaryEdgeSet(M) == bR);
    CHECK(M.size() != R.size());
    {
        std::vector<RankTri> Mshuf{{100, 3, 2}, {4, 1, 100}, {1, 2, 100},
                                   {3, 4, 100}};
        CHECK(cav::auditCavityDiskCombinatorial(Mshuf));
        CHECK(cav::boundaryEdgeSet(Mshuf) == bR);
        std::vector<Rank> c1, c2;
        CHECK(cav::boundaryCycle(M, c1));
        CHECK(cav::boundaryCycle(Mshuf, c2));
        CHECK(c1 == c2);
    }
}

void testCavAuditRejections() {
    TEST("bnd2.cavaudit_rejections (self-edge, dup, disconnected, hole, pinch)");
    {
        std::vector<RankTri> bad{{1, 1, 2}};
        CHECK(!cav::auditCavityDiskCombinatorial(bad));
        CHECK(!cav::auditCavityNodeLink(bad));
        std::vector<Rank> cyc;
        CHECK(!cav::boundaryCycle(bad, cyc));
    }
    {
        std::vector<RankTri> empty;
        CHECK(!cav::auditCavityDiskCombinatorial(empty));
        CHECK(!cav::auditCavityNodeLink(empty));
        std::vector<Rank> cyc;
        CHECK(!cav::boundaryCycle(empty, cyc));
    }
    CHECK(!cav::auditCavityDiskCombinatorial(
        std::vector<RankTri>{{1, 2, 3}, {1, 2, 3}}));
    CHECK(!cav::auditCavityDiskCombinatorial(
        std::vector<RankTri>{{1, 2, 3}, {1, 3, 2}}));
    CHECK(!cav::auditCavityNodeLink(
        std::vector<RankTri>{{1, 2, 3}, {3, 2, 1}}));
    CHECK(!cav::auditCavityDiskCombinatorial(
        std::vector<RankTri>{{1, 2, 3}, {4, 5, 6}}));
    CHECK(!cav::auditCavityDiskCombinatorial(
        std::vector<RankTri>{{1, 2, 3}, {1, 2, 4}, {1, 2, 5}}));
    {
        std::vector<RankTri> pinch{{1, 2, 3}, {1, 4, 5}};
        CHECK(!cav::auditCavityDiskCombinatorial(pinch));
        CHECK(!cav::auditCavityNodeLink(pinch));
    }
    {
        std::vector<RankTri> tet{{1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}};
        std::vector<Rank> cyc;
        CHECK(!cav::boundaryCycle(tet, cyc));
        CHECK(!cav::auditCavityDiskCombinatorial(tet));
    }
    {
        std::vector<RankTri> ann{{1, 2, 5}, {2, 6, 5}, {2, 3, 6}, {3, 7, 6},
                                 {3, 4, 7}, {4, 8, 7}, {4, 1, 8}, {1, 5, 8}};
        std::vector<Rank> cyc;
        CHECK(!cav::boundaryCycle(ann, cyc));
        CHECK(!cav::auditCavityDiskCombinatorial(ann));
    }
}

// ===================================================================== PART 2
using Probe = paacs::topo::detail::TransactionAccessTestProbe;

struct Fix {
    std::unique_ptr<topo::TopologyWorkspace> ws;
    std::vector<NodeHandle> n;
    std::vector<CellHandle> c;
    topo::WorkspaceStamp stamp;   // captured at build (workspace readable)
};

// Unit-square disk: 4 corner nodes, 2 Tri3 across the diagonal n0-n2.
Fix buildQuad() {
    Fix f;
    topo::NodePool np;
    f.n.push_back(np.add({0, 0, 0}, {2, 1}));
    f.n.push_back(np.add({2, 0, 0}, {2, 1}));
    f.n.push_back(np.add({2, 2, 0}, {2, 1}));
    f.n.push_back(np.add({0, 2, 0}, {2, 1}));
    topo::CellPool cp;
    f.c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[1], f.n[2]}, {2, 5}, np));
    f.c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[2], f.n[3]}, {2, 5}, np));
    topo::EdgePool ep;
    for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
    f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    f.stamp = topo::WorkspaceStamp::of(*f.ws);   // readable (Ready)
    return f;
}

topo::SteinerCavitySpec fanSpec(const Fix& f, bool req = false) {
    topo::SteinerCavitySpec s;
    s.stamp = f.stamp;
    s.removed = f.c;
    s.steinerNodes.push_back({0, {1, 1, 0}, {2, 1}});   // localId 0
    auto tri = [&](std::uint32_t a, std::uint32_t b) {
        topo::SteinerCavityCell c;
        c.nodes = {SteinerNodeRef::Existing(f.n[a]),
                   SteinerNodeRef::Existing(f.n[b]),
                   SteinerNodeRef::Staged(0)};
        c.cls = {2, 5};
        return c;
    };
    s.added = {tri(0, 1), tri(1, 2), tri(2, 3), tri(3, 0)};
    if (req)
        s.requirements.push_back(
            {SteinerNodeRef::Existing(f.n[0]), SteinerNodeRef::Staged(0),
             {1, 9}});
    return s;
}

// Two interior Steiner nodes joined by a staged-staged edge. The localId ->
// coordinate mapping and every reference are FIXED (localId 0 = left (0.7,1),
// localId 1 = right (1.3,1)); `reversed` only permutes the steinerNodes VECTOR
// ORDER. Canonicalization by localId must therefore make both spellings
// produce a byte-identical committed topology (the determinism proof).
topo::SteinerCavitySpec multiSpec(const Fix& f, bool reversed,
                                  bool reqStagedStaged) {
    topo::SteinerCavitySpec s;
    s.stamp = f.stamp;
    s.removed = f.c;
    if (!reversed) {
        s.steinerNodes.push_back({0, {0.7, 1, 0}, {2, 1}});
        s.steinerNodes.push_back({1, {1.3, 1, 0}, {2, 1}});
    } else {   // same id->coord, reversed vector order
        s.steinerNodes.push_back({1, {1.3, 1, 0}, {2, 1}});
        s.steinerNodes.push_back({0, {0.7, 1, 0}, {2, 1}});
    }
    auto A = SteinerNodeRef::Staged(0);   // left
    auto B = SteinerNodeRef::Staged(1);   // right
    auto E = [&](std::uint32_t i) { return SteinerNodeRef::Existing(f.n[i]); };
    auto cell = [&](SteinerNodeRef a, SteinerNodeRef b, SteinerNodeRef c) {
        topo::SteinerCavityCell x;
        x.nodes = {a, b, c};
        x.cls = {2, 5};
        return x;
    };
    // {n0,n1,A},{n1,B,A},{n1,n2,B},{n2,n3,B},{n3,A,B},{n3,n0,A}
    s.added = {cell(E(0), E(1), A), cell(E(1), B, A), cell(E(1), E(2), B),
               cell(E(2), E(3), B), cell(E(3), A, B), cell(E(3), E(0), A)};
    if (reqStagedStaged)
        s.requirements.push_back({A, B, {1, 9}});
    return s;
}

void testAscrCommit() {
    TEST("bnd2.ascr_commit (quad -> central Steiner fan, atomic)");
    auto f = buildQuad();
    f.ws->rebuildAdjacency();
    const auto rev0 = f.ws->topologyRevision();
    auto tx = f.ws->beginSteinerCavityTransaction();
    const auto& pv = tx.stage(fanSpec(f));
    CHECK(pv.steinerCount == 1);
    CHECK(pv.cellsRemoved == 2);
    CHECK(pv.cellsAdded == 4);
    CHECK(pv.edgesToCreate == 4);
    CHECK(pv.edgesToRetain == 4);
    CHECK(pv.edgesToRetire == 1);
    CHECK(pv.boundaryRanks.size() == 4);
    const auto out = tx.commit();
    CHECK(out.createdNodes.size() == 1);
    CHECK(out.createdNodes[0].first == 0);
    CHECK(out.resultCells.size() == 4);
    CHECK(out.requiredEdges.empty());
    CHECK(out.revision == rev0 + 1);
    f.ws->rebuildAdjacency();
    CHECK(f.ws->adjacencyReady());
    const NodeHandle s = out.nodeOf(0);
    CHECK(f.ws->nodes().activeBySequence().size() == 5);
    const auto sd = f.ws->nodes().get(s);
    CHECK(sd.coord.x == 1 && sd.coord.y == 1);
    CHECK(f.ws->cells().activeBySequence().size() == 4);
    int spokeHits = 0;
    for (const auto& ch : f.ws->cells().activeBySequence()) {
        const auto sn = f.ws->cells().get(ch);
        for (int i = 0; i < 3; ++i)
            if (sn.nodes[std::size_t(i)] == s) ++spokeHits;
    }
    CHECK(spokeHits == 4);
    CHECK(!f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[2]).has_value());
    for (int i = 0; i < 4; ++i)
        CHECK(f.ws->edges().find(f.ws->nodes(), f.n[std::size_t(i)], s)
                  .has_value());
}

void testAscrMultiSteiner() {
    TEST("bnd2.ascr_multi_steiner (2 staged nodes + staged-staged edge)");
    auto f = buildQuad();
    f.ws->rebuildAdjacency();
    auto tx = f.ws->beginSteinerCavityTransaction();
    const auto& pv = tx.stage(multiSpec(f, false, /*reqSS=*/true));
    CHECK(pv.steinerCount == 2);
    CHECK(pv.cellsAdded == 6);
    CHECK(pv.cellsRemoved == 2);
    CHECK(pv.edgesToCreate == 7);   // 6 spokes + the S0-S1 interior edge
    CHECK(pv.edgesToRetire == 1);   // old diagonal
    const auto out = tx.commit();
    CHECK(out.createdNodes.size() == 2);
    const NodeHandle a = out.nodeOf(0), b = out.nodeOf(1);
    CHECK(!(a == b));
    f.ws->rebuildAdjacency();
    CHECK(f.ws->adjacencyReady());
    CHECK(f.ws->nodes().activeBySequence().size() == 6);
    CHECK(f.ws->cells().activeBySequence().size() == 6);
    // the staged-staged edge exists AND is a required constraint {1,9}
    const auto eh = f.ws->edges().find(f.ws->nodes(), a, b);
    CHECK(eh.has_value());
    if (eh) {
        const auto d = f.ws->edges().get(*eh);
        CHECK(d.requiredConstraint);
        CHECK(d.cls.dim == 1 && d.cls.tag == 9);
    }
    CHECK(out.requiredEdges.size() == 1);
}

void testAscrOrderIndependence() {
    TEST("bnd2.ascr_multi_steiner_order (all input vectors reversed -> same hash)");
    auto a = buildQuad();
    auto b = buildQuad();
    a.ws->rebuildAdjacency();
    b.ws->rebuildAdjacency();
    {  // A declares [localId 0, localId 1]
        auto tx = a.ws->beginSteinerCavityTransaction();
        (void)tx.stage(multiSpec(a, false, true));
        (void)tx.commit();
        a.ws->rebuildAdjacency();
    }
    {  // B declares [localId 1, localId 0] and reverses the other
       // order-insensitive input vectors, with the SAME structural content.
        auto tx = b.ws->beginSteinerCavityTransaction();
        auto spec = multiSpec(b, true, true);
        std::reverse(spec.removed.begin(), spec.removed.end());
        std::reverse(spec.added.begin(), spec.added.end());
        std::reverse(spec.requirements.begin(), spec.requirements.end());
        (void)tx.stage(spec);
        (void)tx.commit();
        b.ws->rebuildAdjacency();
    }
    // Canonicalization by localId makes the committed topology identical.
    CHECK(a.ws->semanticTopologyHash() == b.ws->semanticTopologyHash());
}

void testAscrMappingAndMixedSlots() {
    TEST("bnd2.ascr_mapping_mixed_slots (predicted==committed, reuse+fresh)");
    // Create exactly ONE free node slot. A 2-Steiner batch must therefore
    // consume that reusable slot and materialize one fresh slot.
    Fix f;
    topo::NodePool np;
    f.n.push_back(np.add({0, 0, 0}, {2, 1}));
    f.n.push_back(np.add({2, 0, 0}, {2, 1}));
    f.n.push_back(np.add({2, 2, 0}, {2, 1}));
    f.n.push_back(np.add({0, 2, 0}, {2, 1}));
    NodeHandle scratchA = np.add({9, 9, 0}, {2, 1});
    np.remove(scratchA);   // free slot
    topo::CellPool cp;
    f.c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[1], f.n[2]}, {2, 5}, np));
    f.c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<NodeHandle>{f.n[0], f.n[2], f.n[3]}, {2, 5}, np));
    topo::EdgePool ep;
    for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
    f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    f.stamp = topo::WorkspaceStamp::of(*f.ws);
    f.ws->rebuildAdjacency();
    auto spec = multiSpec(f, false, false);
    auto tx = f.ws->beginSteinerCavityTransaction();
    (void)tx.stage(spec);
    const auto out = tx.commit();
    CHECK(out.createdNodes.size() == 2);
    CHECK(out.nodeOf(0).slot == scratchA.slot);  // reusable slot consumed
    CHECK(out.nodeOf(0).gen == scratchA.gen + 1);
    CHECK(out.nodeOf(1).slot != scratchA.slot);  // genuinely fresh slot
    f.ws->rebuildAdjacency();
    CHECK(f.ws->adjacencyReady());
    // predicted==committed is asserted internally (commitNodeBatch aborts on
    // mismatch under PAACS_TOPO_TESTING); here we confirm the mapping is
    // usable and the two staged nodes materialized at their coords.
    const auto da = f.ws->nodes().get(out.nodeOf(0));
    const auto db = f.ws->nodes().get(out.nodeOf(1));
    CHECK(da.coord.x == 0.7 && db.coord.x == 1.3);
    CHECK(f.ws->cells().activeBySequence().size() == 6);
}

void testAscrRetiredSlotHistory() {
    TEST("bnd2.ascr_retired_slot_history (retired node slot in history)");
    topo::NodePool np;
    std::vector<NodeHandle> n;
    for (int i = 0; i < 4; ++i) n.push_back(np.add({0, 0, 0}, {2, 1}));
    // give the pool a genuinely RETIRED slot (generation exhausted)
    NodeHandle victim = np.add({5, 5, 0}, {2, 1});
    np.testRetire(victim);
    // real square corners
    np.setCoord(n[0], {0, 0, 0});
    np.setCoord(n[1], {2, 0, 0});
    np.setCoord(n[2], {2, 2, 0});
    np.setCoord(n[3], {0, 2, 0});
    topo::CellPool cp;
    std::vector<CellHandle> c;
    c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<NodeHandle>{n[0], n[1], n[2]}, {2, 5}, np));
    c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<NodeHandle>{n[0], n[2], n[3]}, {2, 5}, np));
    topo::EdgePool ep;
    for (const auto& ch : c) ep.insertCellEdges(np, cp, ch);
    Fix f;
    f.n = n;
    f.c = c;
    f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    f.stamp = topo::WorkspaceStamp::of(*f.ws);
    f.ws->rebuildAdjacency();
    auto spec = fanSpec(f);
    auto tx = f.ws->beginSteinerCavityTransaction();
    (void)tx.stage(spec);
    const auto out = tx.commit();
    CHECK(out.createdNodes.size() == 1);
    CHECK(out.nodeOf(0).slot != victim.slot);     // retired means never reused
    CHECK(!f.ws->nodes().isCurrent(victim));
    f.ws->rebuildAdjacency();
    CHECK(f.ws->adjacencyReady());   // retired slot is never reused
}

void testAscrPreparedMove() {
    TEST("bnd2.ascr_prepared_move (prepared ownership and placeholders transfer)");
    auto f = buildQuad();
    f.ws->rebuildAdjacency();
    const auto rev0 = f.ws->topologyRevision();
    auto source = f.ws->beginSteinerCavityTransaction();
    (void)source.stage(multiSpec(f, false, true));
    source.prepare();
    CHECK(source.prepared());
    auto target = std::move(source);
    CHECK(target.prepared());
    const auto out = target.commit();
    CHECK(out.createdNodes.size() == 2);
    CHECK(out.resultCells.size() == 6);
    CHECK(out.requiredEdges.size() == 1);
    CHECK(out.revision == rev0 + 1);
    f.ws->rebuildAdjacency();
    CHECK(f.ws->adjacencyReady());
}

void testAscrSeqBoundaries() {
    TEST("bnd2.ascr_seq_boundaries (remaining == demand / demand-1 / 0, N=2)");
    // demand == 2 staged nodes. Build the pool, set the node allocator's
    // remaining creation sequences to an exact value BEFORE moving it in.
    auto make = [](std::uint64_t remaining) {
        Fix f;
        topo::NodePool np;
        f.n.push_back(np.add({0, 0, 0}, {2, 1}));
        f.n.push_back(np.add({2, 0, 0}, {2, 1}));
        f.n.push_back(np.add({2, 2, 0}, {2, 1}));
        f.n.push_back(np.add({0, 2, 0}, {2, 1}));
        Probe::setRemainingNodeSequences(np, remaining);
        topo::CellPool cp;
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<NodeHandle>{f.n[0], f.n[1], f.n[2]}, {2, 5}, np));
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<NodeHandle>{f.n[0], f.n[2], f.n[3]}, {2, 5}, np));
        topo::EdgePool ep;
        for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
        f.ws = std::make_unique<topo::TopologyWorkspace>(
            std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
        f.stamp = topo::WorkspaceStamp::of(*f.ws);
        f.ws->rebuildAdjacency();
        return f;
    };
    // exactly demand (2): prepare + commit succeed
    {
        auto f = make(2);
        auto tx = f.ws->beginSteinerCavityTransaction();
        (void)tx.stage(multiSpec(f, false, false));
        bool ok = true;
        try { (void)tx.commit(); } catch (...) { ok = false; }
        CHECK(ok);
        f.ws->rebuildAdjacency();
        CHECK(f.ws->adjacencyReady());
    }
    // demand-1 (1): prepare rejects, workspace unchanged, tx poisoned
    {
        auto f = make(1);
        const auto d0 = paacs_tests::observableLogicalDigest(*f.ws);
        auto tx = f.ws->beginSteinerCavityTransaction();
        (void)tx.stage(multiSpec(f, false, false));
        bool threw = false;
        try { tx.prepare(); } catch (const std::exception&) { threw = true; }
        CHECK(threw);
        bool poisoned = false;
        try { (void)tx.commit(); } catch (const std::logic_error&) {
            poisoned = true;
        }
        CHECK(poisoned);
        tx.abandon();
        CHECK(paacs_tests::observableLogicalDigest(*f.ws) == d0);
    }
    // zero remaining: prepare rejects
    {
        auto f = make(0);
        auto tx = f.ws->beginSteinerCavityTransaction();
        (void)tx.stage(multiSpec(f, false, false));
        bool threw = false;
        try { tx.prepare(); } catch (const std::exception&) { threw = true; }
        CHECK(threw);
    }
}

void testAscrStampStale() {
    TEST("bnd2.ascr_stale_candidate_stamp (pre-mutation stamp refused)");
    auto f = buildQuad();
    f.ws->rebuildAdjacency();
    // capture a candidate stamp NOW
    topo::SteinerCavitySpec spec = fanSpec(f);   // stamp captured here
    const auto stamp0 = spec.stamp;
    // now perform an unrelated committed mutation (bumps the revision)
    {
        topo::TopologyPolicy pol;
        pol.permitOrphanNodes = true;
        auto tx = f.ws->beginTransaction(pol);
        tx.addNode({9, 9, 9}, {0, 7});
        (void)tx.commit();
        f.ws->rebuildAdjacency();
    }
    CHECK(!(stamp0 == topo::WorkspaceStamp::of(*f.ws)));   // revision moved
    // handles in `spec` are still current, but the stamp is stale -> refused
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    {
        auto tx = f.ws->beginSteinerCavityTransaction();
        bool threw = false;
        try { (void)tx.stage(spec); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
        tx.abandon();
    }
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testAscrStampForeign() {
    TEST("bnd2.ascr_foreign_stamp (stamp from another workspace refused)");
    auto f = buildQuad();
    auto g = buildQuad();
    f.ws->rebuildAdjacency();
    g.ws->rebuildAdjacency();
    topo::SteinerCavitySpec spec = fanSpec(f);
    spec.stamp = topo::WorkspaceStamp::of(*g.ws);   // foreign pool ids
    auto tx = f.ws->beginSteinerCavityTransaction();
    bool threw = false;
    try { (void)tx.stage(spec); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
}

void testAscrStaleForeignHandles() {
    TEST("bnd2.ascr_stale_and_foreign_handles (removed cell)");
    // stale removed handle: remove-then-superseded is hard; instead use a
    // FOREIGN removed cell (from another workspace) and a stale one via a
    // committed prior Steiner op that consumed the cells.
    {   // foreign removed cell
        auto f = buildQuad();
        auto g = buildQuad();
        f.ws->rebuildAdjacency();
        g.ws->rebuildAdjacency();
        topo::SteinerCavitySpec spec = fanSpec(f);
        spec.removed = g.c;   // cells from a different workspace
        auto tx = f.ws->beginSteinerCavityTransaction();
        bool threw = false;
        try { (void)tx.stage(spec); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }
    {   // stale removed handle: commit a fan, then reuse the old cell handles
        auto f = buildQuad();
        f.ws->rebuildAdjacency();
        const auto oldCells = f.c;
        {
            auto tx = f.ws->beginSteinerCavityTransaction();
            (void)tx.stage(fanSpec(f));
            (void)tx.commit();
            f.ws->rebuildAdjacency();
        }
        topo::SteinerCavitySpec spec;
        spec.stamp = topo::WorkspaceStamp::of(*f.ws);
        spec.removed = oldCells;                     // now stale
        spec.steinerNodes.push_back({0, {1, 1, 0}, {2, 1}});
        spec.added = {};                             // irrelevant; removed fails
        auto tx = f.ws->beginSteinerCavityTransaction();
        bool threw = false;
        try { (void)tx.stage(spec); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }
}

void testAscrOutsideSnapshot() {
    TEST("bnd2.ascr_outside_snapshot (full invariance of outside entities)");
    // 2x1 strip: cells c0,c1 = LEFT quad (removed), c2,c3 = RIGHT (outside).
    topo::NodePool np;
    std::vector<NodeHandle> n;
    n.push_back(np.add({0, 0, 0}, {2, 1}));   // 0
    n.push_back(np.add({1, 0, 0}, {2, 1}));   // 1
    n.push_back(np.add({2, 0, 0}, {2, 1}));   // 2
    n.push_back(np.add({0, 1, 0}, {2, 1}));   // 3
    n.push_back(np.add({1, 1, 0}, {2, 1}));   // 4
    n.push_back(np.add({2, 1, 0}, {2, 1}));   // 5
    topo::CellPool cp;
    std::vector<CellHandle> c;
    c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<NodeHandle>{n[0], n[1], n[4]}, {2, 5}, np));   // c0 left
    c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<NodeHandle>{n[0], n[4], n[3]}, {2, 5}, np));   // c1 left
    c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<NodeHandle>{n[1], n[2], n[5]}, {2, 6}, np));   // c2 right
    c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<NodeHandle>{n[1], n[5], n[4]}, {2, 6}, np));   // c3 right
    topo::EdgePool ep;
    for (const auto& ch : c) ep.insertCellEdges(np, cp, ch);
    Fix f;
    f.n = n;
    f.c = c;
    f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    f.ws->rebuildAdjacency();

    // Full snapshot of the OUTSIDE cells c2,c3 and every node/edge they
    // touch, including the shared cavity/outside boundary n1-n4.
    struct CellSnap { CellHandle h; topo::ElementType t; mesh::Classification cls;
                      std::array<NodeHandle, 3> nn; };
    auto snapCell = [&](CellHandle h) {
        const auto sn = f.ws->cells().get(h);
        return CellSnap{h, sn.type, sn.cls, {sn.nodes[0], sn.nodes[1], sn.nodes[2]}};
    };
    const CellSnap c2a = snapCell(f.c[2]), c3a = snapCell(f.c[3]);
    auto nodeSame = [&](NodeHandle h, const topo::NodePool::NodeData& d0) {
        const auto d = f.ws->nodes().get(h);
        return d.coord.x == d0.coord.x && d.coord.y == d0.coord.y &&
               d.coord.z == d0.coord.z && d.cls.dim == d0.cls.dim &&
               d.cls.tag == d0.cls.tag;
    };
    const auto n1a = f.ws->nodes().get(n[1]);
    const auto n2a = f.ws->nodes().get(n[2]);
    const auto n4a = f.ws->nodes().get(n[4]);
    const auto n5a = f.ws->nodes().get(n[5]);
    struct EdgeSnap {
        topo::EdgeHandle h;
        topo::EdgePool::EdgeData data;
    };
    auto edgeSnap = [&](NodeHandle a, NodeHandle b) {
        const auto eh = f.ws->edges().find(f.ws->nodes(), a, b);
        CHECK(eh.has_value());
        return EdgeSnap{eh ? *eh : topo::EdgeHandle{},
                        eh ? f.ws->edges().get(*eh)
                           : topo::EdgePool::EdgeData{}};
    };
    const auto e_12 = edgeSnap(n[1], n[2]);
    const auto e_25 = edgeSnap(n[2], n[5]);
    const auto e_15 = edgeSnap(n[1], n[5]);   // right interior diagonal
    const auto e_45 = edgeSnap(n[4], n[5]);
    const auto e_14 = edgeSnap(n[1], n[4]);   // shared cavity/outside edge
    const auto rev0 = f.ws->topologyRevision();

    // Steiner-fan the LEFT quad (n0,n1,n4,n3) around a center node.
    topo::SteinerCavitySpec s;
    s.stamp = topo::WorkspaceStamp::of(*f.ws);
    s.removed = {f.c[0], f.c[1]};
    s.steinerNodes.push_back({0, {0.5, 0.5, 0}, {2, 1}});
    auto E = [&](int i) { return SteinerNodeRef::Existing(n[std::size_t(i)]); };
    auto S = SteinerNodeRef::Staged(0);
    auto cell = [&](SteinerNodeRef a, SteinerNodeRef b) {
        topo::SteinerCavityCell x; x.nodes = {a, b, S}; x.cls = {2, 5};
        return x;
    };
    s.added = {cell(E(0), E(1)), cell(E(1), E(4)), cell(E(4), E(3)),
               cell(E(3), E(0))};
    {
        auto tx = f.ws->beginSteinerCavityTransaction();
        (void)tx.stage(s);
        const auto out = tx.commit();
        CHECK(out.resultCells.size() == 4);
        f.ws->rebuildAdjacency();
    }
    CHECK(f.ws->topologyRevision() == rev0 + 1);
    // outside cells unchanged (handle, type, cls, connectivity)
    auto cellSame = [&](const CellSnap& a) {
        if (!f.ws->cells().isCurrent(a.h)) return false;
        const auto sn = f.ws->cells().get(a.h);
        if (sn.type != a.t || sn.cls.dim != a.cls.dim || sn.cls.tag != a.cls.tag)
            return false;
        for (int i = 0; i < 3; ++i)
            if (!(sn.nodes[std::size_t(i)] == a.nn[std::size_t(i)])) return false;
        return true;
    };
    CHECK(cellSame(c2a));
    CHECK(cellSame(c3a));
    CHECK(nodeSame(n[1], n1a));
    CHECK(nodeSame(n[2], n2a));
    CHECK(nodeSame(n[4], n4a));
    CHECK(nodeSame(n[5], n5a));
    auto edgeSame = [&](NodeHandle a, NodeHandle b, const EdgeSnap& s0) {
        const auto eh = f.ws->edges().find(f.ws->nodes(), a, b);
        if (!eh) return false;
        if (!(*eh == s0.h)) return false;
        const auto d = f.ws->edges().get(*eh);
        return d.a == s0.data.a && d.b == s0.data.b &&
               d.useCount == s0.data.useCount &&
               d.requiredConstraint == s0.data.requiredConstraint &&
               d.cls.dim == s0.data.cls.dim &&
               d.cls.tag == s0.data.cls.tag;
    };
    CHECK(edgeSame(n[1], n[2], e_12));
    CHECK(edgeSame(n[2], n[5], e_25));
    CHECK(edgeSame(n[1], n[5], e_15));
    CHECK(edgeSame(n[4], n[5], e_45));
    CHECK(edgeSame(n[1], n[4], e_14));
}

void testAscrDeterminism() {
    TEST("bnd2.ascr_determinism (identical fixtures -> identical hash)");
    auto a = buildQuad();
    auto b = buildQuad();
    a.ws->rebuildAdjacency();
    b.ws->rebuildAdjacency();
    CHECK(a.ws->semanticTopologyHash() == b.ws->semanticTopologyHash());
    { auto tx = a.ws->beginSteinerCavityTransaction();
      (void)tx.stage(fanSpec(a, true)); (void)tx.commit();
      a.ws->rebuildAdjacency(); }
    { auto tx = b.ws->beginSteinerCavityTransaction();
      (void)tx.stage(fanSpec(b, true)); (void)tx.commit();
      b.ws->rebuildAdjacency(); }
    CHECK(a.ws->semanticTopologyHash() == b.ws->semanticTopologyHash());
}

void testAscrAbandonInvariance() {
    TEST("bnd2.ascr_abandon_invariance (stage+prepare then abandon -> no-op)");
    auto f = buildQuad();
    f.ws->rebuildAdjacency();
    const auto rev0 = f.ws->topologyRevision();
    const auto d0 = paacs_tests::observableLogicalDigest(*f.ws);
    {
        auto tx = f.ws->beginSteinerCavityTransaction();
        (void)tx.stage(fanSpec(f, true));
        tx.prepare();
        tx.abandon();
    }
    CHECK(f.ws->state() == topo::WorkspaceState::Ready);
    CHECK(f.ws->topologyRevision() == rev0);
    CHECK(paacs_tests::observableLogicalDigest(*f.ws) == d0);
    {
        auto tx = f.ws->beginSteinerCavityTransaction();
        (void)tx.stage(fanSpec(f));
        const auto out = tx.commit();
        CHECK(out.resultCells.size() == 4);
        f.ws->rebuildAdjacency();
        CHECK(f.ws->adjacencyReady());
    }
}

void testAscrValidatorReject() {
    TEST("bnd2.ascr_validator_reject (false validator -> commit refused)");
    auto f = buildQuad();
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    {
        auto tx = f.ws->beginSteinerCavityTransaction();
        (void)tx.stage(fanSpec(f));
        const bool ok = tx.validate(
            [](const topo::SteinerCavityPreview&) { return false; });
        CHECK(!ok);
        CHECK(tx.rejected());
        bool threw = false;
        try { (void)tx.commit(); } catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
    }
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

template <class Fn>
bool throwsInvalid(Fn&& fn) {
    try { fn(); }
    catch (const std::invalid_argument&) { return true; }
    catch (...) { return false; }
    return false;
}

void testSteinerSymbolicValidation() {
    TEST("bnd2.steiner_symbolic_validation (malformed tokens reject before shifts)");
    using bnd::AffineParam;
    using bnd::makeAffine;

    CHECK(throwsInvalid([&] { (void)makeAffine(7, 7, 1, 1); }));
    CHECK(throwsInvalid([&] { (void)makeAffine(1, 2, 1, 0); }));
    CHECK(throwsInvalid([&] { (void)makeAffine(1, 2, 0, 1); }));
    CHECK(throwsInvalid([&] { (void)makeAffine(1, 2, 2, 1); }));
    CHECK(throwsInvalid([&] {
        (void)makeAffine(1, 2, 1, bnd::kMaxSteinerLevel + 1);
    }));
    CHECK(throwsInvalid([&] {
        (void)bnd::compareDyadic(1, 0, 1, 1);
    }));

    // Validation must happen before narrowing level to uint8_t. Without
    // that check, level 257 aliases level 1 and enters the exact kernel as
    // a different, apparently valid symbolic point.
    const bnd::AuthoritativeParam malformed =
        AffineParam{1, 2, 1, 257};
    auto resolve = [](bnd::ConstraintVertexId id) {
        return id == 1 ? paacs::Vec2{0, 0} : paacs::Vec2{2, 0};
    };
    CHECK(throwsInvalid([&] {
        (void)bnd::toAffinePoint2(malformed, resolve);
    }));
    // Mapping-only conversion is still a public symbolic consumer: it must
    // enforce the same validation before shifting or calling the provider.
    const bnd::AuthoritativeParam overKernelLimit =
        AffineParam{1, 2, 1, bnd::kMaxSteinerLevel + 1};
    CHECK(throwsInvalid([&] {
        (void)bnd::toMappingVec2(overKernelLimit, resolve);
    }));
    const bnd::AuthoritativeParam unreduced = AffineParam{1, 2, 2, 2};
    CHECK(throwsInvalid([&] {
        (void)bnd::toMappingVec2(unreduced, resolve);
    }));
    const bnd::AuthoritativeParam noncanonical = AffineParam{2, 1, 1, 1};
    CHECK(throwsInvalid([&] {
        (void)bnd::toAffinePoint2(noncanonical, resolve);
    }));
}

void testSteinerLoweringUnknownId() {
    TEST("bnd2.steiner_lowering_guards (unknown refs and stale stamps reject)");
    bnd::SteinerRecoveryCandidate c;
    c.stamp = {1, 2, 3, 4, 5};
    c.staged.push_back(
        bnd::SteinerNode{10, bnd::StoredParam{{0.5, 0.5}}, 1,
                         {0.5, 0.5, 0}, {1, 3}});
    const topo::WorkspaceStamp same{1, 2, 3, 4, 5};
    bnd::TriNodeRefTriangle t;
    t.nodes = {bnd::TriNodeRef::Staged(11),  // not declared
               bnd::TriNodeRef::Staged(10),
               bnd::TriNodeRef::Staged(10)};
    c.added.push_back(t);
    CHECK(!bnd::lowerCandidate(c, same).has_value());

    c.added.clear();
    c.requirements.push_back(
        {bnd::TriNodeRef::Staged(10), bnd::TriNodeRef::Staged(99), {1, 3}});
    CHECK(!bnd::lowerCandidate(c, same).has_value());

    // A caller-supplied fresh topology stamp must not launder a candidate
    // built against an older session/workspace revision.
    c.requirements.clear();
    c.resultChain = {bnd::TriNodeRef::Staged(99)};
    CHECK(!bnd::lowerCandidate(c, same).has_value());

    c.resultChain.clear();
    const topo::WorkspaceStamp fresh{1, 2, 3, 4, 6};
    CHECK(!bnd::lowerCandidate(c, fresh).has_value());
}

void testSteinerSymbolicPositive() {
    TEST("bnd2.steiner_symbolic_positive (exact model; levels 1..52, boundaries)");
    using bnd::makeAffine;
    using bnd::AuthoritativeParam;
    using bnd::StoredParam;
    using Sign = paacs::predicates::Sign;
    std::map<bnd::ConstraintVertexId, paacs::Vec2> P{{10, {0, 0}},
                                                     {20, {4, 0}}};
    auto res = [&](bnd::ConstraintVertexId id) { return P.at(id); };
    const paacs::Vec2 loSeg{0, 0}, hiSeg{4, 0};   // segment on v == 0
    // reversed endpoint identity: {10,20}@1/4 == {20,10}@3/4
    auto a1 = makeAffine(10, 20, 1, 2), a2 = makeAffine(20, 10, 3, 2);
    CHECK(a1.lo == a2.lo && a1.hi == a2.hi && a1.num == a2.num &&
          a1.level == a2.level);
    // reduction 2/4 -> 1/2
    auto rr = makeAffine(10, 20, 2, 2);
    CHECK(rr.num == 1 && rr.level == 1);
    // full supported level range 1..52 with BOUNDARY numerators (1, 2^L-1):
    // exact collinearity (on the line) and exact monotone order.
    for (std::uint32_t L = 1; L <= bnd::kMaxSteinerLevel; ++L) {
        const std::uint64_t top = (std::uint64_t(1) << L) - 1;
        AuthoritativeParam pLo = makeAffine(10, 20, 1, L);
        AuthoritativeParam pHi = makeAffine(10, 20, top, L);
        CHECK(bnd::collinearWithSegment(loSeg, hiSeg, pLo, res) == Sign::Zero);
        CHECK(bnd::collinearWithSegment(loSeg, hiSeg, pHi, res) == Sign::Zero);
        if (L >= 2)   // at L==1 both reduce to 1/2 (top==1)
            CHECK(bnd::orderOnSegment(loSeg, hiSeg, pLo, pHi, res) ==
                  Sign::Negative);
    }
    // mixed Stored/Affine ordering: stored lo (t=0) before affine 1/2
    AuthoritativeParam sLo = StoredParam{{0, 0}}, aMid = makeAffine(10, 20, 1, 1);
    CHECK(bnd::orderOnSegment(loSeg, hiSeg, sLo, aMid, res) == Sign::Negative);
    // off-line stored point: nonzero sign
    AuthoritativeParam off = StoredParam{{2, 1}};
    CHECK(bnd::collinearWithSegment(loSeg, hiSeg, off, res) != Sign::Zero);
}

void testSteinerLoweringPositive() {
    TEST("bnd2.steiner_lowering_positive (reordered decls, staged resultChain)");
    bnd::SteinerRecoveryCandidate c;
    const topo::WorkspaceStamp ws{7, 8, 9, 10, 11};
    c.stamp = {7, 8, 9, 10, 11};   // OwnershipStamp matches ws
    // declare staged nodes OUT OF ORDER (id 20 before id 10)
    c.staged.push_back(bnd::SteinerNode{20, bnd::StoredParam{{0.7, 0.5}}, 1,
                                        {0.7, 0.5, 0}, {1, 3}});
    c.staged.push_back(bnd::SteinerNode{10, bnd::StoredParam{{0.3, 0.5}}, 1,
                                        {0.3, 0.5, 0}, {1, 3}});
    bnd::TriNodeRefTriangle t;
    t.nodes = {bnd::TriNodeRef::Staged(10), bnd::TriNodeRef::Staged(20),
               bnd::TriNodeRef::Existing(topo::NodeHandle{})};
    t.cls = {2, 5};
    c.added.push_back(t);
    // §1.1: resultChain carries a STAGED reference (cannot be a NodeHandle
    // before commit)
    c.resultChain = {bnd::TriNodeRef::Existing(topo::NodeHandle{}),
                     bnd::TriNodeRef::Staged(10),
                     bnd::TriNodeRef::Existing(topo::NodeHandle{})};
    CHECK(c.resultChain[1].staged && c.resultChain[1].stagedId == 10);

    auto low = bnd::lowerCandidate(c, ws);
    CHECK(low.has_value());
    if (low) {
        // ascending localId assignment independent of declaration order
        CHECK(low->idToLocal.size() == 2);
        CHECK(low->idToLocal[0].first == 10 && low->idToLocal[0].second == 0);
        CHECK(low->idToLocal[1].first == 20 && low->idToLocal[1].second == 1);
        CHECK(low->spec.steinerNodes.size() == 2);
        CHECK(low->spec.steinerNodes[0].coord.x == 0.3);   // id10 -> local0
        CHECK(low->spec.steinerNodes[1].coord.x == 0.7);   // id20 -> local1
        CHECK(low->spec.added.size() == 1);
        CHECK(low->spec.added[0].nodes[0].isStaged &&
              low->spec.added[0].nodes[0].localId == 0);
        CHECK(low->spec.added[0].nodes[1].isStaged &&
              low->spec.added[0].nodes[1].localId == 1);
    }
}

// ---- single-Steiner recovery ROUTE proofs -------------------------------
struct MockProvider : bnd::GeometryProvider {
    bnd::GeometryProviderStamp st{42, 1};
    enum Mode {
        Ok,
        Nonfinite,
        Domain,
        BadAlloc,
        Other,
        StampBadAlloc,
        StampOther
    } mode = Ok;
    bnd::GeometryProviderStamp stamp() const override {
        if (mode == StampBadAlloc) throw std::bad_alloc();
        if (mode == StampOther) throw std::runtime_error("provider stamp boom");
        return st;
    }
    paacs::Vec3 pointAt(bnd::ChartId, const paacs::Vec2& uv) const override {
        switch (mode) {
            case Nonfinite: return {std::nan(""), 0, 0};
            case Domain: throw bnd::GeometryProviderDomainError("out of domain");
            case BadAlloc: throw std::bad_alloc();
            case Other: throw std::runtime_error("provider boom");
            default: return {uv.u, uv.v, 0};
        }
    }
};

// Quad {n0,n1,n2,n3} split by diagonal n0-n2; the segment n1-n3 (the OTHER
// diagonal) is NOT present and needs one interior Steiner. vids: n0..n3 = 100..103.
struct QuadCS {
    Fix f;
    bnd::ConstraintSet cs;
    bnd::ConstraintSegment seg;
};
QuadCS buildQuadCS(bool reversed = false) {
    QuadCS q;
    q.f = buildQuad();   // n0(0,0) n1(2,0) n2(2,2) n3(0,2); cells across n0-n2
    for (int i = 0; i < 4; ++i)
        q.cs.vertices[100 + std::uint64_t(i)] = q.f.n[std::size_t(i)];
    const paacs::Vec2 uv[4] = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
    for (int i = 0; i < 4; ++i)
        q.cs.params[q.f.n[std::size_t(i)]] = bnd::ParamPoint{uv[i], 7};
    q.seg.id = 500;
    q.seg.a = reversed ? 103 : 101;   // n1 <-> n3
    q.seg.b = reversed ? 101 : 103;
    q.seg.chart = 7;
    q.seg.cls = {1, 9};
    q.cs.segments.push_back(q.seg);
    q.cs.stamp = bnd::OwnershipStamp::of(*q.f.ws);
    return q;
}

// A provider that throws on the Nth pointAt (to force a per-entity Failed in
// the middle of a set).
struct CountingProvider : bnd::GeometryProvider {
    bnd::GeometryProviderStamp st{42, 1};
    mutable int calls = 0;
    int throwOn = 0;   // 0 = never
    bnd::GeometryProviderStamp stamp() const override { return st; }
    paacs::Vec3 pointAt(bnd::ChartId, const paacs::Vec2& uv) const override {
        if (++calls == throwOn) throw std::runtime_error("nth pointAt");
        return {uv.u, uv.v, 0};
    }
};

void testSteinerSetLevel() {
    TEST("bnd2.steiner_set_level (complete report XOR + prefix stop)");
    // two DISJOINT quads far apart; segments 500 (quad1 diag), 501 (quad2
    // diag), 502 (quad1 boundary edge, already present).
    topo::NodePool np;
    std::vector<NodeHandle> n;
    auto add = [&](double x, double y) { n.push_back(np.add({x, y, 0}, {2, 1})); };
    add(0, 0); add(2, 0); add(2, 2); add(0, 2);            // quad1 n0..n3
    add(10, 0); add(12, 0); add(12, 2); add(10, 2);        // quad2 n4..n7
    topo::CellPool cp;
    auto tri = [&](int a, int b, int c) {
        cp.add(topo::ElementType::Tri3,
               std::vector<NodeHandle>{n[std::size_t(a)], n[std::size_t(b)],
                                       n[std::size_t(c)]},
               {2, 5}, np);
    };
    tri(0, 1, 2); tri(0, 2, 3); tri(4, 5, 6); tri(4, 6, 7);
    topo::EdgePool ep;
    for (const auto& ch : cp.activeBySequence()) ep.insertCellEdges(np, cp, ch);
    Fix f;
    f.n = n;
    f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    f.ws->rebuildAdjacency();
    bnd::ConstraintSet cs;
    const paacs::Vec2 uv[8] = {{0, 0}, {2, 0}, {2, 2}, {0, 2},
                               {10, 0}, {12, 0}, {12, 2}, {10, 2}};
    for (int i = 0; i < 8; ++i) {
        cs.vertices[100 + std::uint64_t(i)] = f.n[std::size_t(i)];
        cs.params[f.n[std::size_t(i)]] = bnd::ParamPoint{uv[i], 7};
    }
    cs.segments = {{500, 101, 103, 7, {1, 9}},   // quad1 diagonal (Steiner)
                   {501, 105, 107, 7, {1, 9}},   // quad2 diagonal (Steiner)
                   {502, 100, 101, 7, {1, 9}}};  // quad1 present boundary edge
    cs.stamp = bnd::OwnershipStamp::of(*f.ws);

    CountingProvider prov;
    prov.throwOn = 2;   // second Steiner mapping throws -> 501 Failed
    auto rep = bnd::recoverSteinerSegments(*f.ws, cs, prov, prov.st);
    // 500 Recovered (continue), 501 Failed (stop), 502 unattempted.
    CHECK(rep.processed.size() == 2);
    CHECK(rep.unattempted.size() == 1);
    if (rep.processed.size() == 2) {
        CHECK(rep.processed[0].entity == 500);
        CHECK(rep.processed[0].result.status == bnd::SteinerStatus::Recovered);
        CHECK(rep.processed[1].entity == 501);
        CHECK(rep.processed[1].result.status == bnd::SteinerStatus::Failed);
        CHECK(rep.processed[1].result.fail ==
              bnd::SteinerFail::GeometryProviderError);
    }
    if (rep.unattempted.size() == 1) CHECK(rep.unattempted[0] == 502);
    // completeness: each entity appears exactly once across the two lists.
    std::map<bnd::ConstraintEntityId, int> seen;
    for (const auto& e : rep.processed) ++seen[e.entity];
    for (const auto id : rep.unattempted) ++seen[id];
    CHECK(seen.size() == 3);
    for (const auto& [id, cnt] : seen) { (void)id; CHECK(cnt == 1); }
}

void testSteinerSetWholePreflight() {
    TEST("bnd2.steiner_set_whole_preflight (all rejected; provider not reached)");
    auto q = buildQuadCS();
    q.f.ws->rebuildAdjacency();
    q.cs.stamp = bnd::OwnershipStamp::of(*q.f.ws);
    // This prescribed diagonal properly crosses q.seg. The WHOLE set must be
    // rejected before any provider/per-entity work: every entity is itemized
    // Rejected and none is merely unattempted.
    q.cs.segments.push_back(
        bnd::ConstraintSegment{501, 100, 102, 7, {1, 10}});
    const auto d0 = paacs_tests::observableDigest(*q.f.ws);
    MockProvider prov;
    prov.mode = MockProvider::StampOther;   // must never win over preflight
    const auto rep =
        bnd::recoverSteinerSegments(*q.f.ws, q.cs, prov, prov.st);
    CHECK(rep.processed.size() == 2);
    CHECK(rep.unattempted.empty());
    if (!rep.processed.empty()) {
        CHECK(rep.processed[0].result.status ==
              bnd::SteinerStatus::Rejected);
        CHECK(rep.processed[0].result.reject ==
              bnd::SteinerReject::PslgRejected);
    }
    if (rep.processed.size() == 2) {
        CHECK(rep.processed[0].entity == 500);
        CHECK(rep.processed[1].entity == 501);
        for (const auto& e : rep.processed) {
            CHECK(e.result.status == bnd::SteinerStatus::Rejected);
            CHECK(e.result.reject == bnd::SteinerReject::PslgRejected);
        }
    }
    CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);
}

void testSteinerGates() {
    TEST("bnd2.steiner_gates (one discriminating adversary per reachable gate)");
    using AP = paacs::predicates::AffinePoint2;
    using bnd::SteinerGate;
    auto sp = [](paacs::Vec2 v) { return paacs::predicates::storedPoint(v); };
    const paacs::Vec2 A{0, 0}, B{2, 2}, P{2, 0}, Q{0, 2};
    const AP sMid{A, B, 1, 1};   // exact midpoint (1,1) on A-B
    auto fan = [&](const AP& s) {
        return std::vector<std::array<AP, 3>>{{sp(A), sp(P), s},
                                              {sp(P), sp(B), s},
                                              {sp(B), sp(Q), s},
                                              {sp(Q), sp(A), s}};
    };
    // Ok
    CHECK(bnd::steinerSixGate(A, B, P, Q, sMid, fan(sMid), 3, false, true,
                              false) == SteinerGate::Ok);
    // gate 1 nondegeneracy: a collinear M triangle
    {
        auto m = fan(sMid);
        m[0] = {sp(A), sp(P), sp({4, 0})};   // A,P,(4,0) collinear
        CHECK(bnd::steinerSixGate(A, B, P, Q, sMid, m, 3, false, true, false) ==
              SteinerGate::Degenerate);
    }
    // gate 2 orientation: one triangle wound the other way
    {
        auto m = fan(sMid);
        m[1] = {sp(B), sp(P), sMid};   // reversed winding of {P,B,S}
        CHECK(bnd::steinerSixGate(A, B, P, Q, sMid, m, 3, false, true, false) ==
              SteinerGate::Orientation);
    }
    // gate 3 interiority: S collinear but BEYOND B (not strictly between)
    {
        const AP sOut = sp({5, 5});   // on line y=x, past B
        CHECK(bnd::steinerSixGate(A, B, P, Q, sOut, fan(sMid), 3, false, true,
                                  false) == SteinerGate::Interiority);
    }
    // gate 4 boundary simplicity: P,Q straddle A-B, but their segment meets
    // the supporting line outside the finite A-B segment.
    {
        const paacs::Vec2 P2{10, 0}, Q2{10, 20};
        CHECK(bnd::steinerSixGate(A, B, P2, Q2, sMid, fan(sMid), 3, false, true,
                                  false) == SteinerGate::BoundaryNonSimple);
    }
    // gate 5 overlap/T-junction: both boundary vertices lie on the same side
    // of A-B. Check this prerequisite BEFORE the stronger diagonal-cross test
    // so the gate is independently reachable and its proof is non-vacuous.
    {
        const paacs::Vec2 P2{1, -1}, Q2{2, 0};
        CHECK(bnd::steinerSixGate(A, B, P2, Q2, sMid, fan(sMid), 3, false, true,
                                  false) == SteinerGate::OverlapTJunction);
    }
    // gate 6 chain conformity: chain not Existing,Staged,Existing
    {
        CHECK(bnd::steinerSixGate(A, B, P, Q, sMid, fan(sMid), 2, false, true,
                                  false) == SteinerGate::ChainNonConform);
        CHECK(bnd::steinerSixGate(A, B, P, Q, sMid, fan(sMid), 3, true, true,
                                  false) == SteinerGate::ChainNonConform);
    }
}

void testSteinerDyadicSubdivision() {
    TEST("bnd2.steiner_dyadic_subdivision (canonical reversal-symmetric)");
    using bnd::canonicalDyadicParams;
    using bnd::AffineParam;
    auto toks = [](const std::vector<AffineParam>& v) {
        std::vector<std::pair<std::uint64_t, std::uint32_t>> t;
        for (const auto& p : v) t.push_back({p.num, p.level});
        return t;
    };
    CHECK(canonicalDyadicParams(1, 2, 0).empty());
    CHECK(throwsInvalid(
        [&] { (void)canonicalDyadicParams(1, 1, 0); }));
    // exact expected reduced tokens, ascending by t.
    {
        auto d1 = canonicalDyadicParams(1, 2, 1);
        CHECK((toks(d1) ==
               std::vector<std::pair<std::uint64_t, std::uint32_t>>{{1, 1}}));
    }
    {
        auto d2 = canonicalDyadicParams(1, 2, 2);   // 1/4, 3/4
        CHECK((toks(d2) ==
               std::vector<std::pair<std::uint64_t, std::uint32_t>>{{1, 2},
                                                                    {3, 2}}));
    }
    {
        auto d3 = canonicalDyadicParams(1, 2, 3);   // 1/4, 1/2, 3/4
        CHECK((toks(d3) ==
               std::vector<std::pair<std::uint64_t, std::uint32_t>>{
                   {1, 2}, {1, 1}, {3, 2}}));
    }
    {
        auto d4 = canonicalDyadicParams(1, 2, 4);   // 1/4,3/8,5/8,3/4
        CHECK((toks(d4) ==
               std::vector<std::pair<std::uint64_t, std::uint32_t>>{
                   {1, 2}, {3, 3}, {5, 3}, {3, 2}}));
    }
    // reversal symmetry: reversing endpoints yields the SAME reduced token set.
    for (std::uint32_t D = 1; D <= 6; ++D) {
        auto fwd = toks(canonicalDyadicParams(1, 2, D));
        auto rev = toks(canonicalDyadicParams(2, 1, D));
        // The public contract says both calls return the same canonical,
        // ascending token vector; sorting here alone would hide an
        // orientation-dependent output-order bug.
        CHECK(fwd == rev);
        std::sort(fwd.begin(), fwd.end());
        std::sort(rev.begin(), rev.end());
        CHECK(fwd == rev);
        // determinism: identical repeated call.
        CHECK(toks(canonicalDyadicParams(1, 2, D)) ==
              toks(canonicalDyadicParams(1, 2, D)));
        // all D distinct.
        auto s = toks(canonicalDyadicParams(1, 2, D));
        std::sort(s.begin(), s.end());
        CHECK(std::adjacent_find(s.begin(), s.end()) == s.end());
        CHECK(s.size() == D);
    }
    // exact minimum-feature: in the D=3 set {1/4,1/2,3/4} adjacent gaps are 1/4.
    {
        auto d3 = canonicalDyadicParams(1, 2, 3);
        CHECK(bnd::dyadicSeparated(d3[0], d3[1], 2));    // gap 1/4 >= 2^-2
        CHECK(!bnd::dyadicSeparated(d3[0], d3[1], 1));   // gap 1/4 <  2^-1
        CHECK(bnd::dyadicSeparated(d3[0], d3[2], 1));    // gap 1/2 >= 2^-1
    }
    // The separation helper is public predicate-adjacent infrastructure:
    // reject malformed tokens, an invalid policy level, and parameters from
    // different undirected segments BEFORE performing any shift.
    {
        const auto p = bnd::makeAffine(1, 2, 1, 2);
        const auto q = bnd::makeAffine(1, 2, 3, 2);
        CHECK(throwsInvalid([&] { (void)bnd::dyadicSeparated(p, q, 0); }));
        CHECK(throwsInvalid(
            [&] { (void)bnd::dyadicSeparated(p, q, bnd::kMaxSteinerLevel + 1); }));

        auto malformed = p;
        malformed.level = 0;
        CHECK(throwsInvalid(
            [&] { (void)bnd::dyadicSeparated(malformed, q, 2); }));

        const auto foreignSegment = bnd::makeAffine(1, 3, 3, 2);
        CHECK(throwsInvalid(
            [&] { (void)bnd::dyadicSeparated(p, foreignSegment, 2); }));
    }
}

void testSteinerRouteRecovered() {
    TEST("bnd2.steiner_route_recovered (one interior node; chain a-S-b)");
    auto q = buildQuadCS();
    q.f.ws->rebuildAdjacency();
    const auto rev0 = q.f.ws->topologyRevision();
    MockProvider prov;
    auto res = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st);
    CHECK(res.status == bnd::SteinerStatus::Recovered);
    CHECK(res.recoveredChain.size() == 3);
    CHECK(res.recoveredChain[0] == q.f.n[1]);
    CHECK(res.recoveredChain[2] == q.f.n[3]);
    CHECK(res.createdParams.size() == 1);
    CHECK(res.revision == rev0 + 1);
    q.f.ws->rebuildAdjacency();
    CHECK(q.f.ws->adjacencyReady());
    const topo::NodeHandle S = res.recoveredChain[1];
    // the segment is now present as the chain n1-S-n3, both sub-edges required
    const auto e1 = q.f.ws->edges().find(q.f.ws->nodes(), q.f.n[1], S);
    const auto e2 = q.f.ws->edges().find(q.f.ws->nodes(), S, q.f.n[3]);
    CHECK(e1.has_value() && e2.has_value());
    if (e1) CHECK(q.f.ws->edges().get(*e1).requiredConstraint);
    if (e2) CHECK(q.f.ws->edges().get(*e2).requiredConstraint);
    // the crossed diagonal n0-n2 is gone
    CHECK(!q.f.ws->edges().find(q.f.ws->nodes(), q.f.n[0], q.f.n[2]).has_value());
    CHECK(q.f.ws->nodes().activeBySequence().size() == 5);
    CHECK(q.f.ws->cells().activeBySequence().size() == 4);
}

void testSteinerRouteReversedEquivalent() {
    TEST("bnd2.steiner_route_reversed (reversed segment -> identical mesh)");
    auto a = buildQuadCS(false);
    auto b = buildQuadCS(true);
    a.f.ws->rebuildAdjacency();
    b.f.ws->rebuildAdjacency();
    MockProvider prov;
    const auto ra = bnd::recoverSteinerSegment(*a.f.ws, a.cs, a.seg, prov, prov.st);
    const auto rb = bnd::recoverSteinerSegment(*b.f.ws, b.cs, b.seg, prov, prov.st);
    CHECK(ra.status == bnd::SteinerStatus::Recovered);
    CHECK(rb.status == bnd::SteinerStatus::Recovered);
    a.f.ws->rebuildAdjacency();
    b.f.ws->rebuildAdjacency();
    CHECK(a.f.ws->semanticTopologyHash() == b.f.ws->semanticTopologyHash());
}

void testSteinerRouteAlreadyPresent() {
    TEST("bnd2.steiner_route_already_present (required matching edge exists)");
    auto q = buildQuadCS();
    q.f.ws->rebuildAdjacency();
    MockProvider prov;
    bnd::ConstraintSegment diag = q.seg;
    diag.a = 100; diag.b = 102;   // n0-n2 IS present
    q.cs.segments = {diag};
    const auto eh =
        q.f.ws->edges().find(q.f.ws->nodes(), q.f.n[0], q.f.n[2]);
    CHECK(eh.has_value());
    if (eh) {
        q.f.ws->testEdgesMutable().setRequiredConstraint(*eh, true);
        q.f.ws->testEdgesMutable().setClassification(*eh, diag.cls);
    }
    const auto d0 = paacs_tests::observableDigest(*q.f.ws);
    auto res = bnd::recoverSteinerSegment(*q.f.ws, q.cs, diag, prov, prov.st);
    CHECK(res.status == bnd::SteinerStatus::AlreadyPresent);
    CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);   // untouched
}

void testSteinerRouteUnmarkedExisting() {
    TEST("bnd2.steiner_route_unmarked_existing (present is not AlreadyPresent)");
    auto q = buildQuadCS();
    q.f.ws->rebuildAdjacency();
    MockProvider prov;
    bnd::ConstraintSegment diag = q.seg;
    diag.a = 100;
    diag.b = 102;   // direct edge exists, but is not a required constraint
    q.cs.segments = {diag};
    const auto eh =
        q.f.ws->edges().find(q.f.ws->nodes(), q.f.n[0], q.f.n[2]);
    CHECK(eh.has_value());
    if (eh) CHECK(!q.f.ws->edges().get(*eh).requiredConstraint);

    const auto res =
        bnd::recoverSteinerSegment(*q.f.ws, q.cs, diag, prov, prov.st);
    CHECK(res.status == bnd::SteinerStatus::Recovered);
    const auto after =
        q.f.ws->edges().find(q.f.ws->nodes(), q.f.n[0], q.f.n[2]);
    CHECK(after.has_value());
    if (after) {
        const auto d = q.f.ws->edges().get(*after);
        CHECK(d.requiredConstraint);
        CHECK(d.cls.dim == diag.cls.dim && d.cls.tag == diag.cls.tag);
    }
}

void testSteinerRouteConstraintPreflight() {
    TEST("bnd2.steiner_route_constraint_preflight (stamp/params/classification)");
    MockProvider prov;
    auto runRejected = [&](auto corrupt, bnd::SteinerReject why) {
        auto q = buildQuadCS();
        q.f.ws->rebuildAdjacency();
        q.cs.stamp = bnd::OwnershipStamp::of(*q.f.ws);
        corrupt(q);
        q.cs.segments = {q.seg};
        const auto d0 = paacs_tests::observableDigest(*q.f.ws);
        const auto res =
            bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st);
        CHECK(res.status == bnd::SteinerStatus::Rejected);
        CHECK(res.reject == why);
        CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);
    };

    runRejected(
        [](QuadCS& q) { ++q.cs.stamp.revision; },
        bnd::SteinerReject::StaleCandidateStamp);
    runRejected(
        [](QuadCS& q) { q.cs.params.erase(q.f.n[0]); },
        bnd::SteinerReject::PslgRejected);
    runRejected(
        [](QuadCS& q) { q.cs.params[q.f.n[0]].chart = 8; },
        bnd::SteinerReject::ChartMismatch);
    runRejected(
        [](QuadCS& q) { q.seg.cls.tag = -1; },
        bnd::SteinerReject::InvalidClassification);
}

void testSteinerRouteProviderStampTaxonomy() {
    TEST("bnd2.steiner_route_provider_stamp_taxonomy (stamp throws are structured)");
    auto run = [&](MockProvider::Mode mode) {
        auto q = buildQuadCS();
        q.f.ws->rebuildAdjacency();
        q.cs.stamp = bnd::OwnershipStamp::of(*q.f.ws);
        MockProvider prov;
        prov.mode = mode;
        const auto d0 = paacs_tests::observableDigest(*q.f.ws);
        bnd::SteinerResult res;
        bool escaped = false;
        try {
            res = bnd::recoverSteinerSegment(
                *q.f.ws, q.cs, q.seg, prov, prov.st);
        } catch (...) {
            escaped = true;
        }
        CHECK(!escaped);
        CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);
        return res;
    };

    {
        const auto r = run(MockProvider::StampBadAlloc);
        CHECK(r.status == bnd::SteinerStatus::Failed);
        CHECK(r.fail == bnd::SteinerFail::CandidateAllocFailed);
    }
    {
        const auto r = run(MockProvider::StampOther);
        CHECK(r.status == bnd::SteinerStatus::Failed);
        CHECK(r.fail == bnd::SteinerFail::GeometryProviderError);
    }
}

void testSteinerRouteMetadataPublication() {
    TEST("bnd2.steiner_route_metadata_publication (constraint state is atomic)");
    auto q = buildQuadCS();
    q.f.ws->rebuildAdjacency();
    q.cs.stamp = bnd::OwnershipStamp::of(*q.f.ws);
    MockProvider prov;
    const auto res =
        bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st);
    CHECK(res.status == bnd::SteinerStatus::Recovered);
    CHECK(res.recoveredChain.size() == 3);
    if (res.recoveredChain.size() == 3) {
        const auto* pp = q.cs.paramOf(res.recoveredChain[1]);
        CHECK(pp != nullptr);
        if (pp) CHECK(pp->chart == q.seg.chart);
        const auto* ap = q.cs.authParamOf(res.recoveredChain[1]);
        CHECK(ap != nullptr);
        if (ap) {
            CHECK(ap->affine);
            CHECK(ap->lo == 101 && ap->hi == 103);
            CHECK(ap->num == 1 && ap->level == 1);
            CHECK(ap->chart == q.seg.chart);
        }
    }
    CHECK(q.cs.stamp.matches(*q.f.ws));
}

void testSteinerRouteIdempotentChain() {
    TEST("bnd2.steiner_route_idempotent_chain (recovered chain is AlreadyPresent)");
    auto q = buildQuadCS();
    q.f.ws->rebuildAdjacency();
    q.cs.stamp = bnd::OwnershipStamp::of(*q.f.ws);
    MockProvider prov;
    const auto first =
        bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st);
    CHECK(first.status == bnd::SteinerStatus::Recovered);
    CHECK(first.recoveredChain.size() == 3);
    const auto d1 = paacs_tests::observableDigest(*q.f.ws);
    const auto rev1 = q.f.ws->topologyRevision();

    const auto second =
        bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st);
    CHECK(second.status == bnd::SteinerStatus::AlreadyPresent);
    CHECK(second.recoveredChain == first.recoveredChain);
    CHECK(q.f.ws->topologyRevision() == rev1);
    CHECK(paacs_tests::observableDigest(*q.f.ws) == d1);
}

void testSteinerRouteFullSessionPreflight() {
    TEST("bnd2.steiner_route_full_session_preflight (crossing PSLG rejects)");
    auto q = buildQuadCS();
    q.f.ws->rebuildAdjacency();
    q.cs.stamp = bnd::OwnershipStamp::of(*q.f.ws);
    // The prescribed existing diagonal properly crosses q.seg. A1 accepts
    // both individually; the established whole-set PSLG preflight rejects
    // their crossing before any segment is attempted.
    q.cs.segments.push_back(
        bnd::ConstraintSegment{501, 100, 102, 7, {1, 10}});
    const auto d0 = paacs_tests::observableDigest(*q.f.ws);
    MockProvider prov;
    const auto res =
        bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st);
    CHECK(res.status == bnd::SteinerStatus::Rejected);
    CHECK(res.reject == bnd::SteinerReject::PslgRejected);
    CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);
}

void testSteinerRouteExistingClassificationConflict() {
    TEST("bnd2.steiner_route_existing_class_conflict (not crossed-edge reason)");
    auto q = buildQuadCS();
    q.f.ws->rebuildAdjacency();
    MockProvider prov;
    bnd::ConstraintSegment diag = q.seg;
    diag.a = 100;
    diag.b = 102;
    q.cs.segments = {diag};
    const auto eh =
        q.f.ws->edges().find(q.f.ws->nodes(), q.f.n[0], q.f.n[2]);
    CHECK(eh.has_value());
    if (eh) {
        q.f.ws->testEdgesMutable().setRequiredConstraint(*eh, true);
        q.f.ws->testEdgesMutable().setClassification(*eh, {1, 4});
    }
    const auto d0 = paacs_tests::observableDigest(*q.f.ws);
    const auto res =
        bnd::recoverSteinerSegment(*q.f.ws, q.cs, diag, prov, prov.st);
    CHECK(res.status == bnd::SteinerStatus::Rejected);
    CHECK(res.reject == bnd::SteinerReject::InvalidClassification);
    CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);
}

void testSteinerRouteAuthoritativeParamAudit() {
    TEST("bnd2.steiner_route_auth_param_audit (foreign/stale records reject)");
    auto q = buildQuadCS();
    q.f.ws->rebuildAdjacency();
    q.cs.stamp = bnd::OwnershipStamp::of(*q.f.ws);
    const topo::NodeHandle bogus{
        999, 1, q.f.ws->nodes().poolId()};
    q.cs.authParams.emplace(bogus, bnd::AuthoritativeParamRecord{});
    const auto d0 = paacs_tests::observableDigest(*q.f.ws);
    MockProvider prov;
    const auto res =
        bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st);
    CHECK(res.status == bnd::SteinerStatus::Rejected);
    CHECK(res.reject == bnd::SteinerReject::PslgRejected);
    CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);
}

void testSteinerRouteProviderTaxonomy() {
    TEST("bnd2.steiner_route_provider_taxonomy (stale/nonfinite/badalloc/other/domain)");
    MockProvider prov;
    // (a) stale provider identity/revision -> Rejected(StaleGeometryProvider)
    {
        auto q = buildQuadCS();
        q.f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*q.f.ws);
        bnd::GeometryProviderStamp wrong{42, 2};
        auto res = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, wrong);
        CHECK(res.status == bnd::SteinerStatus::Rejected);
        CHECK(res.reject == bnd::SteinerReject::StaleGeometryProvider);
        CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);
    }
    auto run = [&](MockProvider::Mode m) {
        auto q = buildQuadCS();
        q.f.ws->rebuildAdjacency();
        MockProvider p;
        p.mode = m;
        const auto d0 = paacs_tests::observableDigest(*q.f.ws);
        auto res = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, p, p.st);
        CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);   // nothing committed
        return res;
    };
    {
        auto r = run(MockProvider::Nonfinite);
        CHECK(r.status == bnd::SteinerStatus::Rejected);
        CHECK(r.reject == bnd::SteinerReject::NonfiniteMappedPosition);
    }
    {
        auto r = run(MockProvider::Domain);
        CHECK(r.status == bnd::SteinerStatus::Rejected);
        CHECK(r.reject == bnd::SteinerReject::ProviderDomainRejected);
    }
    {
        auto r = run(MockProvider::BadAlloc);
        CHECK(r.status == bnd::SteinerStatus::Failed);
        CHECK(r.fail == bnd::SteinerFail::CandidateAllocFailed);
    }
    {
        auto r = run(MockProvider::Other);
        CHECK(r.status == bnd::SteinerStatus::Failed);
        CHECK(r.fail == bnd::SteinerFail::GeometryProviderError);
    }
}

void testSteinerRoutePolicyAndOutside() {
    TEST("bnd2.steiner_route_policy_outside (disabled refuse; outside invariance)");
    // Steiner disabled -> Refused
    {
        auto q = buildQuadCS();
        q.f.ws->rebuildAdjacency();
        MockProvider prov;
        bnd::SteinerPolicy pol;
        pol.steinerEnabled = false;
        const auto d0 = paacs_tests::observableDigest(*q.f.ws);
        auto res =
            bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
        CHECK(res.status == bnd::SteinerStatus::Refused);
        CHECK(res.refuse == bnd::SteinerRefuse::SteinerDisabled);
        CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);
    }
    // outside invariance: an extra triangle sharing n1-n2 stays untouched.
    {
        topo::NodePool np;
        std::vector<NodeHandle> n;
        n.push_back(np.add({0, 0, 0}, {2, 1}));   // 0
        n.push_back(np.add({2, 0, 0}, {2, 1}));   // 1
        n.push_back(np.add({2, 2, 0}, {2, 1}));   // 2
        n.push_back(np.add({0, 2, 0}, {2, 1}));   // 3
        n.push_back(np.add({4, 1, 0}, {2, 1}));   // 4 (outside, right of n1-n2)
        topo::CellPool cp;
        std::vector<CellHandle> c;
        c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<NodeHandle>{n[0], n[1], n[2]}, {2, 5}, np));
        c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<NodeHandle>{n[0], n[2], n[3]}, {2, 5}, np));
        c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<NodeHandle>{n[1], n[4], n[2]}, {2, 6}, np));  // outside
        topo::EdgePool ep;
        for (const auto& ch : c) ep.insertCellEdges(np, cp, ch);
        Fix f;
        f.n = n; f.c = c;
        f.ws = std::make_unique<topo::TopologyWorkspace>(
            std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
        f.ws->rebuildAdjacency();
        bnd::ConstraintSet cs;
        for (int i = 0; i < 4; ++i) cs.vertices[100 + std::uint64_t(i)] = f.n[std::size_t(i)];
        const paacs::Vec2 uv[4] = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
        for (int i = 0; i < 4; ++i) cs.params[f.n[std::size_t(i)]] = bnd::ParamPoint{uv[i], 7};
        cs.params[f.n[4]] = bnd::ParamPoint{{4, 1}, 7};
        cs.stamp = bnd::OwnershipStamp::of(*f.ws);
        bnd::ConstraintSegment seg;
        seg.id = 500; seg.a = 101; seg.b = 103; seg.chart = 7; seg.cls = {1, 9};
        cs.segments.push_back(seg);
        // snapshot outside cell c[2] + node n4
        const auto out2 = f.ws->cells().get(f.c[2]);
        const auto n4a = f.ws->nodes().get(f.n[4]);
        const auto e14a = f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[2]);
        const auto e14use = e14a ? f.ws->edges().get(*e14a).useCount : 0u;
        MockProvider prov;
        auto res = bnd::recoverSteinerSegment(*f.ws, cs, seg, prov, prov.st);
        CHECK(res.status == bnd::SteinerStatus::Recovered);
        f.ws->rebuildAdjacency();
        CHECK(f.ws->adjacencyReady());
        // outside triangle unchanged (handle current, same connectivity/cls)
        CHECK(f.ws->cells().isCurrent(f.c[2]));
        if (f.ws->cells().isCurrent(f.c[2])) {
            const auto o = f.ws->cells().get(f.c[2]);
            CHECK(o.cls.dim == out2.cls.dim && o.cls.tag == out2.cls.tag);
            for (int i = 0; i < 3; ++i) CHECK(o.nodes[std::size_t(i)] == out2.nodes[std::size_t(i)]);
        }
        const auto n4b = f.ws->nodes().get(f.n[4]);
        CHECK(n4b.coord.x == n4a.coord.x && n4b.coord.y == n4a.coord.y);
        // shared edge n1-n2 preserved (still used by the outside cell)
        const auto e14b = f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[2]);
        CHECK(e14b.has_value());
        if (e14b) CHECK(f.ws->edges().get(*e14b).useCount == e14use);
    }
}

void testSteinerRouteDeterminism() {
    TEST("bnd2.steiner_route_determinism (identical fixtures -> identical hash)");
    auto a = buildQuadCS();
    auto b = buildQuadCS();
    a.f.ws->rebuildAdjacency();
    b.f.ws->rebuildAdjacency();
    MockProvider prov;
    (void)bnd::recoverSteinerSegment(*a.f.ws, a.cs, a.seg, prov, prov.st);
    (void)bnd::recoverSteinerSegment(*b.f.ws, b.cs, b.seg, prov, prov.st);
    a.f.ws->rebuildAdjacency();
    b.f.ws->rebuildAdjacency();
    CHECK(a.f.ws->semanticTopologyHash() == b.f.ws->semanticTopologyHash());
}

void testAscrRejections() {
    TEST("bnd2.ascr_rejections (staged-aware Phase-1 refusals)");
    // (a) no Steiner nodes
    {
        auto f = buildQuad();
        f.ws->rebuildAdjacency();
        auto tx = f.ws->beginSteinerCavityTransaction();
        auto s = fanSpec(f);
        s.steinerNodes.clear();
        CHECK(throwsInvalid([&] { (void)tx.stage(s); }));
        tx.abandon();
    }
    // (b) duplicate localId
    {
        auto f = buildQuad();
        f.ws->rebuildAdjacency();
        auto tx = f.ws->beginSteinerCavityTransaction();
        auto s = fanSpec(f);
        s.steinerNodes.push_back({0, {0.5, 0.5, 0}, {2, 1}});   // dup id 0
        CHECK(throwsInvalid([&] { (void)tx.stage(s); }));
    }
    // (c) staged node declared but never used
    {
        auto f = buildQuad();
        f.ws->rebuildAdjacency();
        auto tx = f.ws->beginSteinerCavityTransaction();
        auto s = fanSpec(f);
        s.steinerNodes.push_back({1, {0.5, 0.5, 0}, {2, 1}});   // id 1 unused
        CHECK(throwsInvalid([&] { (void)tx.stage(s); }));
    }
    // (d) reference to unknown staged localId
    {
        auto f = buildQuad();
        f.ws->rebuildAdjacency();
        auto tx = f.ws->beginSteinerCavityTransaction();
        auto s = fanSpec(f);
        s.added[0].nodes[2] = SteinerNodeRef::Staged(7);   // no such id
        CHECK(throwsInvalid([&] { (void)tx.stage(s); }));
    }
    // (e) boundary not preserved (open fan)
    {
        auto f = buildQuad();
        f.ws->rebuildAdjacency();
        auto tx = f.ws->beginSteinerCavityTransaction();
        auto s = fanSpec(f);
        s.added.pop_back();
        CHECK(throwsInvalid([&] { (void)tx.stage(s); }));
    }
    // (f) duplicate triangle in M
    {
        auto f = buildQuad();
        f.ws->rebuildAdjacency();
        auto tx = f.ws->beginSteinerCavityTransaction();
        auto s = fanSpec(f);
        s.added.push_back(s.added[0]);
        CHECK(throwsInvalid([&] { (void)tx.stage(s); }));
    }
    // (g) foreign node (existing node not in R)
    {
        auto f = buildQuad();
        NodeHandle foreign;
        {
            topo::TopologyPolicy pol; pol.permitOrphanNodes = true;
            f.ws->rebuildAdjacency();
            auto tx = f.ws->beginTransaction(pol);
            tx.addNode({9, 9, 0}, {2, 1});
            auto res = tx.commit();
            foreign = res.node;
            f.ws->rebuildAdjacency();
        }
        auto s = fanSpec(f);
        s.stamp = topo::WorkspaceStamp::of(*f.ws);   // refresh after mutation
        s.added[0].nodes[1] = SteinerNodeRef::Existing(foreign);
        auto tx = f.ws->beginSteinerCavityTransaction();
        CHECK(throwsInvalid([&] { (void)tx.stage(s); }));
    }
    // (h) requirement outside the cavity closure
    {
        auto f = buildQuad();
        f.ws->rebuildAdjacency();
        auto tx = f.ws->beginSteinerCavityTransaction();
        auto s = fanSpec(f);
        s.requirements.push_back(
            {SteinerNodeRef::Existing(f.n[1]),
             SteinerNodeRef::Existing(f.n[3]), {1, 9}});   // n1-n3 not in closure
        CHECK(throwsInvalid([&] { (void)tx.stage(s); }));
    }
}

// ======================================================================
// Increment 4 - MULTI-Steiner recovery ROUTE proofs.
//
// Fan strip: A,B on the x-axis; m top nodes P_i=(i,1) (i=1..m) and one bottom
// apex Q. The segment A-B crosses every P_i-Q edge, so the transverse walk
// produces a K=m demand recovered with m canonical dyadic Steiner nodes in a
// SINGLE atomic ASCR commit. Every node carries a real ConstraintVertexId, no
// vertex lies interior to the segment (so the accepted PSLG/A2 preflight
// passes), and there are no crossing constraints.
// ----------------------------------------------------------------------
struct FanCS {
    Fix f;
    bnd::ConstraintSet cs;
    bnd::ConstraintSegment seg;
    std::vector<NodeHandle> outside;   // optional far cells' nodes
    CellHandle outsideCell{};
};
// m>=1 top nodes -> K=m. `rev` swaps the segment endpoints. `churn` retires
// that many scratch nodes first (allocator history). `concaveDip` lowers the
// middle top node to make the top cavity boundary reflex (concave multi-cell).
// `withOutside` adds one far, declared, same-winding triangle.
FanCS buildFanCS(int m, bool rev = false, int churn = 0,
                 bool concaveDip = false, bool withOutside = false) {
    FanCS q;
    topo::NodePool np;
    for (int k = 0; k < churn; ++k) { auto h = np.add({-100.0 - k, 0, 0}, {2, 1}); np.remove(h); }
    std::vector<NodeHandle> P;
    NodeHandle A = np.add({0, 0, 0}, {2, 1});
    std::vector<paacs::Vec2> Puv;
    for (int i = 1; i <= m; ++i) {
        double y = 1.0;
        if (concaveDip && m >= 3 && i == (m + 1) / 2) y = 0.3;   // reflex dip
        P.push_back(np.add({double(i), y, 0}, {2, 1}));
        Puv.push_back({double(i), y});
    }
    NodeHandle B = np.add({double(m + 1), 0, 0}, {2, 1});
    NodeHandle Q = np.add({(m + 1) / 2.0, -1, 0}, {2, 1});
    NodeHandle Z1{}, Z2{}, Z3{};
    if (withOutside) {   // same (CW) winding as the fan, far away
        Z1 = np.add({50, 50, 0}, {2, 1});
        Z2 = np.add({50, 51, 0}, {2, 1});
        Z3 = np.add({51, 50, 0}, {2, 1});
    }
    topo::CellPool cp;
    auto tri = [&](NodeHandle a, NodeHandle b, NodeHandle c) {
        return cp.add(topo::ElementType::Tri3,
                      std::vector<NodeHandle>{a, b, c}, {2, 5}, np);
    };
    tri(A, P[0], Q);
    for (int i = 0; i + 1 < m; ++i) tri(P[std::size_t(i)], P[std::size_t(i + 1)], Q);
    tri(P[std::size_t(m - 1)], B, Q);
    if (withOutside) q.outsideCell = tri(Z1, Z2, Z3);
    topo::EdgePool ep;
    for (const auto& ch : cp.activeBySequence()) ep.insertCellEdges(np, cp, ch);
    q.f.n = P; q.f.n.insert(q.f.n.begin(), A); q.f.n.push_back(B); q.f.n.push_back(Q);
    q.f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    q.f.ws->rebuildAdjacency();
    auto V = [&](std::uint64_t id, NodeHandle h, paacs::Vec2 uv) {
        q.cs.vertices[id] = h; q.cs.params[h] = bnd::ParamPoint{uv, 7};
    };
    V(100, A, {0, 0}); V(101, B, {double(m + 1), 0});
    for (int i = 0; i < m; ++i) V(200 + std::uint64_t(i), P[std::size_t(i)], Puv[std::size_t(i)]);
    V(300, Q, {(m + 1) / 2.0, -1});
    if (withOutside) {
        V(400, Z1, {50, 50}); V(401, Z2, {50, 51}); V(402, Z3, {51, 50});
        q.outside = {Z1, Z2, Z3};
    }
    q.seg.id = 500; q.seg.a = rev ? 101 : 100; q.seg.b = rev ? 100 : 101;
    q.seg.chart = 7; q.seg.cls = {1, 9};
    q.cs.segments.push_back(q.seg);
    q.cs.stamp = bnd::OwnershipStamp::of(*q.f.ws);
    return q;
}

// A provider that counts EVERY call (stamp + pointAt), to prove the provider is
// untouched on a policy refusal.
struct TallyProvider : bnd::GeometryProvider {
    bnd::GeometryProviderStamp st{42, 1};
    mutable int stampCalls = 0, pointCalls = 0;
    bnd::GeometryProviderStamp stamp() const override { ++stampCalls; return st; }
    paacs::Vec3 pointAt(bnd::ChartId, const paacs::Vec2& uv) const override {
        ++pointCalls; return {uv.u, uv.v, 0};
    }
};

void testSteinerMultiTwoThree() {
    TEST("bnd2.steiner_multi_two_three (K=2 and K=3 dyadic chains recovered)");
    MockProvider prov;
    bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 52;
    for (int m = 2; m <= 3; ++m) {
        auto q = buildFanCS(m);
        const auto r = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
        CHECK(r.status == bnd::SteinerStatus::Recovered);
        CHECK(r.createdParams.size() == std::size_t(m));
        CHECK(r.recoveredChain.size() == std::size_t(m + 2));   // A + m + B
        // every created node has an affine authoritative param on (100,101).
        for (const auto& pr : r.createdParams) {
            const auto* ap = q.cs.authParamOf(pr.first);
            CHECK(ap != nullptr);
            if (ap) { CHECK(ap->affine); CHECK(ap->lo == 100 && ap->hi == 101); CHECK(ap->chart == 7); }
        }
        // the interior chain nodes are pairwise distinct and strictly ordered.
        CHECK(q.cs.stamp.matches(*q.f.ws));
    }
}

void testSteinerMultiReversed() {
    TEST("bnd2.steiner_multi_reversed (a->b and b->a give identical mesh, K=2,3)");
    MockProvider prov;
    bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 52;
    for (int m = 2; m <= 3; ++m) {
        auto a = buildFanCS(m, false);
        auto b = buildFanCS(m, true);
        const auto ra = bnd::recoverSteinerSegment(*a.f.ws, a.cs, a.seg, prov, prov.st, pol);
        const auto rb = bnd::recoverSteinerSegment(*b.f.ws, b.cs, b.seg, prov, prov.st, pol);
        CHECK(ra.status == bnd::SteinerStatus::Recovered);
        CHECK(rb.status == bnd::SteinerStatus::Recovered);
        a.f.ws->rebuildAdjacency();
        b.f.ws->rebuildAdjacency();
        CHECK(a.f.ws->semanticTopologyHash() == b.f.ws->semanticTopologyHash());
    }
}

void testSteinerMultiDeterminism() {
    TEST("bnd2.steiner_multi_determinism (allocator history -> identical mesh)");
    MockProvider prov;
    bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 52;
    auto a = buildFanCS(3, false, 0);
    auto b = buildFanCS(3, false, 37);   // 37 retired scratch nodes first
    (void)bnd::recoverSteinerSegment(*a.f.ws, a.cs, a.seg, prov, prov.st, pol);
    (void)bnd::recoverSteinerSegment(*b.f.ws, b.cs, b.seg, prov, prov.st, pol);
    a.f.ws->rebuildAdjacency();
    b.f.ws->rebuildAdjacency();
    CHECK(a.f.ws->semanticTopologyHash() == b.f.ws->semanticTopologyHash());
}

void testSteinerMultiBudgetAndFeature() {
    TEST("bnd2.steiner_multi_budget_feature (budget + min-feature refuse, no mutation)");
    MockProvider prov;
    {   // K=3 but budget 2 -> SplitBudgetExceeded, workspace untouched
        auto q = buildFanCS(3);
        const auto d0 = paacs_tests::observableDigest(*q.f.ws);
        bnd::SteinerPolicy pol; pol.splitBudget = 2; pol.minFeatureLevel = 52;
        const auto r = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
        CHECK(r.status == bnd::SteinerStatus::Refused);
        CHECK(r.refuse == bnd::SteinerRefuse::SplitBudgetExceeded);
        CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);
    }
    {   // K=3 with a coarse minimum feature -> MinimumFeatureReached, untouched
        auto q = buildFanCS(3);
        const auto d0 = paacs_tests::observableDigest(*q.f.ws);
        bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 1;
        const auto r = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
        CHECK(r.status == bnd::SteinerStatus::Refused);
        CHECK(r.refuse == bnd::SteinerRefuse::MinimumFeatureReached);
        CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);
    }
}

void testSteinerMultiInvalidPolicyProviderUntouched() {
    TEST("bnd2.steiner_multi_invalid_policy (Rejected; provider PROVABLY untouched)");
    auto q = buildFanCS(2);
    const auto d0 = paacs_tests::observableDigest(*q.f.ws);
    TallyProvider prov;
    bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 0;   // invalid
    const auto r = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
    CHECK(r.status == bnd::SteinerStatus::Rejected);
    CHECK(r.reject == bnd::SteinerReject::InvalidSplitPolicy);
    CHECK(prov.stampCalls == 0);   // no provider identity read
    CHECK(prov.pointCalls == 0);   // no mapping
    CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);
    // level > kMaxSteinerLevel is likewise Rejected(InvalidSplitPolicy).
    auto q2 = buildFanCS(2);
    TallyProvider prov2;
    bnd::SteinerPolicy pol2; pol2.splitBudget = 8;
    pol2.minFeatureLevel = bnd::kMaxSteinerLevel + 1;
    const auto r2 = bnd::recoverSteinerSegment(*q2.f.ws, q2.cs, q2.seg, prov2, prov2.st, pol2);
    CHECK(r2.status == bnd::SteinerStatus::Rejected);
    CHECK(r2.reject == bnd::SteinerReject::InvalidSplitPolicy);
    CHECK(prov2.stampCalls == 0 && prov2.pointCalls == 0);
}

void testSteinerMultiConcaveAndOutside() {
    TEST("bnd2.steiner_multi_concave_outside (reflex multi-cell; outside invariant)");
    MockProvider prov;
    bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 52;
    auto q = buildFanCS(3, false, 0, /*concaveDip=*/true, /*withOutside=*/true);
    const auto before = snapshotOutsideTriangle(*q.f.ws, q.outsideCell);
    const auto r = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
    CHECK(r.status == bnd::SteinerStatus::Recovered);
    CHECK(r.createdParams.size() == 3);
    q.f.ws->rebuildAdjacency();
    CHECK(outsideTriangleMatches(*q.f.ws, before));
}

void testSteinerMultiIdempotent() {
    TEST("bnd2.steiner_multi_idempotent (re-recover a K=3 chain -> AlreadyPresent)");
    MockProvider prov;
    bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 52;
    auto q = buildFanCS(3);
    const auto first = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
    CHECK(first.status == bnd::SteinerStatus::Recovered);
    q.f.ws->rebuildAdjacency();
    const auto d1 = paacs_tests::observableDigest(*q.f.ws);
    const auto rev1 = q.f.ws->topologyRevision();
    const auto second = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
    CHECK(second.status == bnd::SteinerStatus::AlreadyPresent);
    CHECK(second.recoveredChain == first.recoveredChain);
    CHECK(q.f.ws->topologyRevision() == rev1);
    CHECK(paacs_tests::observableDigest(*q.f.ws) == d1);
}

void testSteinerMultiSetPrefix() {
    TEST("bnd2.steiner_multi_set_prefix (multi-Steiner set: prefix stop + completeness)");
    // Two disjoint fan strips (K=2 each) + a present boundary edge segment. The
    // provider throws on the 3rd pointAt so the SECOND strip fails mid-set; the
    // report must be complete (each entity once) with an ascending-id prefix.
    topo::NodePool np;
    auto mkfan = [&](double ox, std::uint64_t base, std::vector<NodeHandle>& out) {
        NodeHandle A = np.add({ox + 0, 0, 0}, {2, 1});
        NodeHandle P1 = np.add({ox + 1, 1, 0}, {2, 1});
        NodeHandle P2 = np.add({ox + 2, 1, 0}, {2, 1});
        NodeHandle B = np.add({ox + 3, 0, 0}, {2, 1});
        NodeHandle Q = np.add({ox + 1.5, -1, 0}, {2, 1});
        out = {A, P1, P2, B, Q};
        (void)base;
    };
    std::vector<NodeHandle> f1, f2;
    mkfan(0, 100, f1);
    mkfan(20, 110, f2);
    topo::CellPool cp;
    auto tri = [&](NodeHandle a, NodeHandle b, NodeHandle c) {
        cp.add(topo::ElementType::Tri3, std::vector<NodeHandle>{a, b, c}, {2, 5}, np);
    };
    auto strip = [&](std::vector<NodeHandle>& f) {
        tri(f[0], f[1], f[4]); tri(f[1], f[2], f[4]); tri(f[2], f[3], f[4]);
    };
    strip(f1); strip(f2);
    topo::EdgePool ep;
    for (const auto& ch : cp.activeBySequence()) ep.insertCellEdges(np, cp, ch);
    Fix f;
    f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    f.ws->rebuildAdjacency();
    bnd::ConstraintSet cs;
    auto V = [&](std::uint64_t id, NodeHandle h, double x, double y) {
        cs.vertices[id] = h; cs.params[h] = bnd::ParamPoint{{x, y}, 7};
    };
    V(100, f1[0], 0, 0); V(101, f1[1], 1, 1); V(102, f1[2], 2, 1); V(103, f1[3], 3, 0); V(104, f1[4], 1.5, -1);
    V(110, f2[0], 20, 0); V(111, f2[1], 21, 1); V(112, f2[2], 22, 1); V(113, f2[3], 23, 0); V(114, f2[4], 21.5, -1);
    cs.segments = {{500, 100, 103, 7, {1, 9}},    // strip 1 diagonal (K=2)
                   {501, 110, 113, 7, {1, 9}},    // strip 2 diagonal (K=2)
                   {502, 100, 101, 7, {1, 9}}};   // present boundary edge
    cs.stamp = bnd::OwnershipStamp::of(*f.ws);
    CountingProvider prov;
    prov.throwOn = 3;   // strip1 uses pointAt 1,2 (ok); strip2's first (3rd) throws
    bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 52;
    const auto rep = bnd::recoverSteinerSegments(*f.ws, cs, prov, prov.st, pol);
    CHECK(rep.processed.size() == 2);
    CHECK(rep.unattempted.size() == 1);
    if (rep.processed.size() == 2) {
        CHECK(rep.processed[0].entity == 500);
        CHECK(rep.processed[0].result.status == bnd::SteinerStatus::Recovered);
        CHECK(rep.processed[1].entity == 501);
        CHECK(rep.processed[1].result.status == bnd::SteinerStatus::Failed);
    }
    if (rep.unattempted.size() == 1) CHECK(rep.unattempted[0] == 502);
    // completeness: each of the 3 entities appears exactly once.
    std::map<bnd::ConstraintEntityId, int> seen;
    for (const auto& e : rep.processed) ++seen[e.entity];
    for (const auto id : rep.unattempted) ++seen[id];
    CHECK(seen.size() == 3);
    for (const auto& [id, cnt] : seen) { (void)id; CHECK(cnt == 1); }
}

// ======================================================================
// Increment 4 - EXISTING-NODE REUSE proofs (blocker 2). A manifold 6-triangle
// fan around an interior node X that lies exactly on the segment A-B. X is a
// current mesh node with a StoredParam (parameter totality) but is NOT a
// declared constraint vertex, so it passes the accepted PSLG/A2 preflight; the
// Steiner strip walk reuses it as a chain waypoint. The two external (non-
// crossed) cells of X are pulled into the cavity so X becomes a strictly
// interior existing node and the whole recovery is ONE atomic ASCR commit.
// ----------------------------------------------------------------------
struct ReuseCS {
    Fix f;
    bnd::ConstraintSet cs;
    bnd::ConstraintSegment seg;
    NodeHandle X;
    std::vector<NodeHandle> outside;
    CellHandle outsideCell{};
};
// asym=false -> X at t=0.5; asym=true -> X at t=0.3 (isolates staged/reused
// min-feature). extClass sets the class of ONE external cell of X (for the
// classification-conflict proof). withOutside adds a far same-winding triangle.
ReuseCS buildReuseCS(bool rev = false, mesh::Classification extClass = {2, 5},
                     bool asym = false, bool withOutside = false,
                     bool coincidentSplit = false,
                     mesh::Classification reusedClass = {1, 9}) {
    ReuseCS q;
    const double xx = coincidentSplit ? 1.0 : (asym ? 1.2 : 2.0);
    const double lx = coincidentSplit ? 0.5 : (asym ? 0.6 : 1.0);
    const double rx = coincidentSplit ? 2.5 : (asym ? 2.6 : 3.0);
    topo::NodePool np;
    NodeHandle A = np.add({0, 0, 0}, {2, 1});
    NodeHandle U1 = np.add({lx, 1, 0}, {2, 1});
    NodeHandle L1 = np.add({lx, -1, 0}, {2, 1});
    NodeHandle X = np.add({xx, 0, 0}, reusedClass);
    NodeHandle U2 = np.add({rx, 1, 0}, {2, 1});
    NodeHandle L2 = np.add({rx, -1, 0}, {2, 1});
    NodeHandle B = np.add({4, 0, 0}, {2, 1});
    NodeHandle Z1{}, Z2{}, Z3{};
    if (withOutside) {
        Z1 = np.add({50, 50, 0}, {2, 1});
        Z2 = np.add({50, 51, 0}, {2, 1});
        Z3 = np.add({51, 50, 0}, {2, 1});
    }
    topo::CellPool cp;
    auto tri = [&](NodeHandle a, NodeHandle b, NodeHandle c, mesh::Classification cl) {
        return cp.add(topo::ElementType::Tri3,
                      std::vector<NodeHandle>{a, b, c}, cl, np);
    };
    tri(A, U1, L1, {2, 5});
    tri(U1, X, L1, {2, 5});
    tri(U1, U2, X, extClass);      // an EXTERNAL cell of X
    tri(X, U2, L2, {2, 5});
    tri(X, L2, L1, {2, 5});
    tri(U2, B, L2, {2, 5});
    if (withOutside)
        q.outsideCell = tri(Z1, Z2, Z3, {2, 5});
    topo::EdgePool ep;
    for (const auto& ch : cp.activeBySequence()) ep.insertCellEdges(np, cp, ch);
    q.X = X;
    q.f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    q.f.ws->rebuildAdjacency();
    auto V = [&](std::uint64_t id, NodeHandle h, double x, double y) {
        q.cs.vertices[id] = h; q.cs.params[h] = bnd::ParamPoint{{x, y}, 7};
    };
    V(100, A, 0, 0); V(101, B, 4, 0);
    V(200, U1, lx, 1); V(201, L1, lx, -1); V(202, U2, rx, 1); V(203, L2, rx, -1);
    q.cs.params[X] = bnd::ParamPoint{{xx, 0}, 7};   // reused: param only
    if (withOutside) {
        V(400, Z1, 50, 50); V(401, Z2, 50, 51); V(402, Z3, 51, 50);
        q.outside = {Z1, Z2, Z3};
    }
    q.seg.id = 500; q.seg.a = rev ? 101 : 100; q.seg.b = rev ? 100 : 101;
    q.seg.chart = 7; q.seg.cls = {1, 9};
    q.cs.segments.push_back(q.seg);
    q.cs.stamp = bnd::OwnershipStamp::of(*q.f.ws);
    return q;
}

void testSteinerReuse() {
    TEST("bnd2.steiner_reuse (interior existing node reused as chain waypoint)");
    MockProvider prov;
    bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 52;
    auto q = buildReuseCS();
    const auto r = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
    CHECK(r.status == bnd::SteinerStatus::Recovered);
    CHECK(r.createdParams.size() == 2);        // 2 staged, plus the reused X
    CHECK(r.recoveredChain.size() == 5);        // A, s1, X, s2, B
    bool xin = false;
    for (const auto& h : r.recoveredChain) if (h == q.X) xin = true;
    CHECK(xin);                                 // the reused node is in the chain
    // reversed segment -> byte-identical committed mesh.
    auto b = buildReuseCS(true);
    const auto rb = bnd::recoverSteinerSegment(*b.f.ws, b.cs, b.seg, prov, prov.st, pol);
    CHECK(rb.status == bnd::SteinerStatus::Recovered);
    q.f.ws->rebuildAdjacency();
    b.f.ws->rebuildAdjacency();
    CHECK(q.f.ws->semanticTopologyHash() == b.f.ws->semanticTopologyHash());
    // idempotent re-recovery.
    const auto d1 = paacs_tests::observableDigest(*q.f.ws);
    const auto second = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
    CHECK(second.status == bnd::SteinerStatus::AlreadyPresent);
    CHECK(paacs_tests::observableDigest(*q.f.ws) == d1);
}

void testSteinerReuseMinFeature() {
    TEST("bnd2.steiner_reuse_min_feature (exact staged-to-existing separation)");
    MockProvider prov;
    // asymmetric: X at t=0.3; a dyadic staged node lands at t=0.25 (0.05 from X)
    // while staged/staged (0.5) and staged/endpoint (0.25) are larger. At
    // minFeatureLevel 4 (2^-4 = 0.0625) ONLY the staged/reused pair is too close.
    {
        auto q = buildReuseCS(false, {2, 5}, /*asym=*/true);
        const auto d0 = paacs_tests::observableDigest(*q.f.ws);
        bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 4;
        const auto r = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
        CHECK(r.status == bnd::SteinerStatus::Refused);
        CHECK(r.refuse == bnd::SteinerRefuse::MinimumFeatureReached);
        CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);   // untouched
    }
    {   // a finer minimum feature (level 5) admits the same configuration.
        auto q = buildReuseCS(false, {2, 5}, /*asym=*/true);
        bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 5;
        const auto r = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
        CHECK(r.status == bnd::SteinerStatus::Recovered);
    }
}

void testSteinerReuseClassConflict() {
    TEST("bnd2.steiner_reuse_class_conflict (node and cavity classification)");
    MockProvider prov;
    // one external cell of X carries a different owning-face class; pulling it
    // into the cavity to make X interior must be an InvalidClassification reject.
    auto q = buildReuseCS(false, {2, 7}, /*asym=*/true);
    const auto d0 = paacs_tests::observableDigest(*q.f.ws);
    bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 52;
    const auto r = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
    CHECK(r.status == bnd::SteinerStatus::Rejected);
    CHECK(r.reject == bnd::SteinerReject::InvalidClassification);
    CHECK(paacs_tests::observableDigest(*q.f.ws) == d0);

    // The reused node itself becomes an interior node of the required
    // dim-1 chain. A face-classified node may not be silently retained as
    // though it carried the segment classification.
    auto q2 = buildReuseCS(false, {2, 5}, false, false, false, {2, 1});
    const auto d1 = paacs_tests::observableDigest(*q2.f.ws);
    const auto r2 =
        bnd::recoverSteinerSegment(*q2.f.ws, q2.cs, q2.seg, prov, prov.st, pol);
    CHECK(r2.status == bnd::SteinerStatus::Rejected);
    CHECK(r2.reject == bnd::SteinerReject::InvalidClassification);
    CHECK(paacs_tests::observableDigest(*q2.f.ws) == d1);

    // An arbitrary dimension-0 classification is not enough either. The
    // BND-2A design has no vertex-to-curve ownership relation with which to
    // prove that this geometric vertex belongs to the prescribed curve. An
    // interior reused node therefore must already carry the exact dim-1
    // segment classification; accepting any dim-0 tag would silently absorb
    // an unrelated GVertex into the middle of the constraint.
    auto q3 = buildReuseCS(false, {2, 5}, false, false, false, {0, 77});
    const auto d2 = paacs_tests::observableDigest(*q3.f.ws);
    const auto r3 =
        bnd::recoverSteinerSegment(*q3.f.ws, q3.cs, q3.seg, prov, prov.st, pol);
    CHECK(r3.status == bnd::SteinerStatus::Rejected);
    CHECK(r3.reject == bnd::SteinerReject::InvalidClassification);
    CHECK(paacs_tests::observableDigest(*q3.f.ws) == d2);
}

void testSteinerReuseCoincidentSplit() {
    TEST("bnd2.steiner_reuse_coincident_split (reuse replaces dyadic creation)");
    MockProvider prov;
    bnd::SteinerPolicy pol;
    pol.splitBudget = 1;       // two crossing events, but one is satisfied by X
    pol.minFeatureLevel = 52;
    // K=2 proposes {1/4,3/4}; X is the current Stored node at exactly 1/4.
    // Reuse-before-create must remove that staged split, count only the
    // remaining staged node against the budget, and recover without a
    // coincident node.
    auto q = buildReuseCS(false, {2, 5}, false, false, true, {1, 9});
    const auto r =
        bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
    CHECK(r.status == bnd::SteinerStatus::Recovered);
    CHECK(r.createdParams.size() == 1);
    CHECK(r.recoveredChain.size() == 4);   // A, X, staged(3/4), B
    int occurrences = 0;
    for (const auto& h : r.recoveredChain)
        if (h == q.X) ++occurrences;
    CHECK(occurrences == 1);
}

void testSteinerReuseOutsideInvariance() {
    TEST("bnd2.steiner_reuse_outside (reuse route leaves outside cells invariant)");
    MockProvider prov;
    bnd::SteinerPolicy pol; pol.splitBudget = 8; pol.minFeatureLevel = 52;
    auto q = buildReuseCS(false, {2, 5}, false, /*withOutside=*/true);
    const auto before = snapshotOutsideTriangle(*q.f.ws, q.outsideCell);
    const auto r = bnd::recoverSteinerSegment(*q.f.ws, q.cs, q.seg, prov, prov.st, pol);
    CHECK(r.status == bnd::SteinerStatus::Recovered);
    q.f.ws->rebuildAdjacency();
    CHECK(outsideTriangleMatches(*q.f.ws, before));
}

void testSteinerGeneralGeometryValidator() {
    TEST("bnd2.steiner_general_validator (six production-path adversaries)");
    using Issue = bnd::GeneralSteinerGeometryIssue;
    using Vertex = bnd::GeneralSteinerGeometryVertex;
    using paacs::predicates::storedPoint;

    auto validFan = [&]() {
        bnd::GeneralSteinerGeometry g;
        g.vertices = {
            Vertex{storedPoint({0, 0}), false},
            Vertex{storedPoint({2, 2}), false},
            Vertex{storedPoint({4, 0}), false},
            Vertex{storedPoint({2, -2}), false},
            Vertex{storedPoint({2, 0}), true},
        };
        g.triangles = {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};
        return g;
    };
    struct Adversary {
        Issue expected;
        bnd::GeneralSteinerGeometry geometry;
    };
    std::vector<Adversary> cases;
    {
        auto g = validFan();
        g.vertices.push_back(Vertex{storedPoint({1.5, 0.5}), false});
        cases.push_back({Issue::ForeignVertexInside, std::move(g)});
    }
    {
        auto g = validFan();
        g.vertices.push_back(Vertex{storedPoint({1, 1}), false});
        cases.push_back({Issue::UnintendedEdgeTJunction, std::move(g)});
    }
    {
        bnd::GeneralSteinerGeometry g;
        // Two consistently clockwise triangles whose vertices are mutually
        // exterior, but whose nonadjacent edges cross properly.
        g.vertices = {
            Vertex{storedPoint({0, 0}), false},
            Vertex{storedPoint({2, 3}), false},
            Vertex{storedPoint({4, 0}), false},
            Vertex{storedPoint({0, 2}), false},
            Vertex{storedPoint({4, 2}), false},
            Vertex{storedPoint({2, -1}), false},
        };
        g.triangles = {{0, 1, 2}, {3, 4, 5}};
        cases.push_back({Issue::NonadjacentEdgeCrossing, std::move(g)});
    }
    {
        auto g = validFan();
        std::swap(g.triangles[1][1], g.triangles[1][2]);
        cases.push_back({Issue::InconsistentOrientation, std::move(g)});
    }
    {
        auto g = validFan();
        g.vertices[4].staged = false;
        g.vertices.push_back(Vertex{storedPoint({10, 10}), true});
        cases.push_back({Issue::StagedPointOffSegment, std::move(g)});
    }
    {
        auto g = validFan();
        g.vertices[4].staged = false;
        g.vertices.push_back(Vertex{storedPoint({5, 0}), true});
        cases.push_back({Issue::StagedPointBeyondSegment, std::move(g)});
    }

    MockProvider prov;
    bnd::SteinerPolicy policy;
    policy.splitBudget = 8;
    policy.minFeatureLevel = 52;
    for (const auto& adversary : cases) {
        CHECK(bnd::validateGeneralSteinerGeometry(
                  {0, 0}, {4, 0}, storedPoint({0, 0}), storedPoint({4, 0}),
                  adversary.geometry) == adversary.expected);

        // Substitute the same plain adversary at the production validator
        // boundary, then verify the real public outcome and tier-1 equality.
        auto q = buildFanCS(3);
        const auto before = paacs_tests::observableDigest(*q.f.ws);
        const std::size_t paramsBefore = q.cs.params.size();
        const std::size_t authBefore = q.cs.authParams.size();
        bnd::GeneralSteinerGeometryTestScope scope(adversary.geometry);
        const auto r = bnd::recoverSteinerSegment(
            *q.f.ws, q.cs, q.seg, prov, prov.st, policy);
        CHECK(r.status == bnd::SteinerStatus::Failed);
        CHECK(r.fail == bnd::SteinerFail::InvalidSteinerGeometry);
        CHECK(paacs_tests::observableDigest(*q.f.ws) == before);
        CHECK(q.cs.params.size() == paramsBefore);
        CHECK(q.cs.authParams.size() == authBefore);
    }
}

void testSteinerStripWalkGuards() {
    TEST("bnd2.steiner_strip_walk_guards (collinear edge + monotone waypoint)");
    {
        auto q = buildReuseCS();
        const NodeHandle A = q.cs.vertices.at(100);
        const NodeHandle B = q.cs.vertices.at(101);
        const NodeHandle U1 = q.cs.vertices.at(200);
        q.cs.params.at(U1).p = {1, 0};  // edge A-U1 overlaps A-B
        const auto idx = bnd::buildTriAdjacency(*q.f.ws);
        CHECK(idx.has_value());
        if (idx) {
            const auto walk =
                bnd::swalk::steinerStripWalk(*q.f.ws, q.cs, *idx, A, B);
            CHECK(!walk.ok);
            CHECK(walk.why == bnd::RefusedReason::CollinearOverlapWithEdge);
        }
    }
    {
        auto q = buildReuseCS();
        const NodeHandle A = q.cs.vertices.at(100);
        const NodeHandle B = q.cs.vertices.at(101);
        // The topology reaches X after the first crossed edge, but its
        // authoritative parameter lies beyond B. The exact per-leg progress
        // guard must refuse it instead of starting a backward/overshooting leg.
        q.cs.params.at(q.X).p = {5, 0};
        const auto idx = bnd::buildTriAdjacency(*q.f.ws);
        CHECK(idx.has_value());
        if (idx) {
            const auto walk =
                bnd::swalk::steinerStripWalk(*q.f.ws, q.cs, *idx, A, B);
            CHECK(!walk.ok);
            CHECK(walk.why == bnd::RefusedReason::TangentialTouch);
        }
    }
}

void testSteinerExactEarClip() {
    TEST("bnd2.steiner_exact_earclip (exact affine ear-clip + axis compare)");
    using paacs::predicates::AffinePoint2;
    using paacs::predicates::storedPoint;
    using S = paacs::predicates::Sign;
    // exact axis comparison of a stored point and an affine midpoint.
    const AffinePoint2 a = storedPoint({0, 0});
    const AffinePoint2 b = storedPoint({4, 0});
    const AffinePoint2 mid{{0, 0}, {4, 0}, 1, 1};   // t=1/2 -> (2,0)
    CHECK(paacs::predicates::compareOnAxisAffine(a, mid, 0) == S::Negative);
    CHECK(paacs::predicates::compareOnAxisAffine(mid, b, 0) == S::Negative);
    CHECK(paacs::predicates::compareOnAxisAffine(mid, storedPoint({2, 0}), 0) == S::Zero);
    // a CONCAVE (reflex) pentagon ear-clips into 3 consistently-oriented
    // triangles - the ear ranking and validity are exact, not rounded.
    std::vector<AffinePoint2> P = {
        storedPoint({0, 0}), storedPoint({4, 0}), storedPoint({4, 4}),
        storedPoint({2, 1}), storedPoint({0, 4})};   // (2,1) is the reflex apex
    auto ok = paacs::bnd::scdetail::earClipOriented(P, {0, 1, 2, 3, 4}, S::Positive);
    CHECK(ok.has_value());
    if (ok) {
        CHECK(ok->size() == 3);
        for (const auto& t : *ok)
            CHECK(paacs::predicates::orient2dAffine(P[t[0]], P[t[1]], P[t[2]]) == S::Positive);
    }
    // a fully collinear ring is degenerate -> nullopt.
    std::vector<AffinePoint2> R = {storedPoint({0, 0}), storedPoint({1, 0}),
                                   storedPoint({2, 0})};
    CHECK(!paacs::bnd::scdetail::earClipOriented(R, {0, 1, 2}, S::Positive).has_value());
    // exact validator primitives: proper crossing and on-closed-segment, with a
    // STAGED (affine) point exactly on the segment (the case rounding breaks).
    const AffinePoint2 p00 = storedPoint({0, 0}), p40 = storedPoint({4, 0});
    const AffinePoint2 up = storedPoint({2, 2}), dn = storedPoint({2, -2});
    CHECK(paacs::bnd::scdetail::properCross(p00, p40, up, dn));      // cross
    CHECK(!paacs::bnd::scdetail::properCross(p00, p40, up, storedPoint({2, 1})));
    const AffinePoint2 onSeg{{0, 0}, {4, 0}, 1, 2};                  // t=1/4 -> (1,0)
    CHECK(paacs::bnd::scdetail::onClosedSeg(p00, p40, onSeg));       // exactly on
    CHECK(!paacs::bnd::scdetail::onClosedSeg(p00, p40, storedPoint({1, 1})));
    CHECK(!paacs::bnd::scdetail::onClosedSeg(p00, p40, storedPoint({5, 0})));  // beyond
}

}  // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* bnd2Tests(std::size_t& count) {
    static const Entry k[] = {
        {"bnd2.cavaudit_disk_boundary", testCavAuditDiskAndBoundary},
        {"bnd2.cavaudit_rejections", testCavAuditRejections},
        {"bnd2.ascr_commit", testAscrCommit},
        {"bnd2.ascr_multi_steiner", testAscrMultiSteiner},
        {"bnd2.ascr_multi_steiner_order", testAscrOrderIndependence},
        {"bnd2.ascr_mapping_mixed_slots", testAscrMappingAndMixedSlots},
        {"bnd2.ascr_retired_slot_history", testAscrRetiredSlotHistory},
        {"bnd2.ascr_prepared_move", testAscrPreparedMove},
        {"bnd2.ascr_seq_boundaries", testAscrSeqBoundaries},
        {"bnd2.ascr_stale_candidate_stamp", testAscrStampStale},
        {"bnd2.ascr_foreign_stamp", testAscrStampForeign},
        {"bnd2.ascr_stale_and_foreign_handles", testAscrStaleForeignHandles},
        {"bnd2.ascr_outside_snapshot", testAscrOutsideSnapshot},
        {"bnd2.ascr_determinism", testAscrDeterminism},
        {"bnd2.ascr_abandon_invariance", testAscrAbandonInvariance},
        {"bnd2.ascr_validator_reject", testAscrValidatorReject},
        {"bnd2.steiner_symbolic_validation", testSteinerSymbolicValidation},
        {"bnd2.steiner_lowering_guards", testSteinerLoweringUnknownId},
        {"bnd2.steiner_symbolic_positive", testSteinerSymbolicPositive},
        {"bnd2.steiner_lowering_positive", testSteinerLoweringPositive},
        {"bnd2.steiner_gates", testSteinerGates},
        {"bnd2.steiner_dyadic_subdivision", testSteinerDyadicSubdivision},
        {"bnd2.steiner_set_level", testSteinerSetLevel},
        {"bnd2.steiner_set_whole_preflight", testSteinerSetWholePreflight},
        {"bnd2.steiner_route_recovered", testSteinerRouteRecovered},
        {"bnd2.steiner_route_reversed", testSteinerRouteReversedEquivalent},
        {"bnd2.steiner_route_already_present", testSteinerRouteAlreadyPresent},
        {"bnd2.steiner_route_unmarked_existing", testSteinerRouteUnmarkedExisting},
        {"bnd2.steiner_route_constraint_preflight", testSteinerRouteConstraintPreflight},
        {"bnd2.steiner_route_provider_stamp_taxonomy", testSteinerRouteProviderStampTaxonomy},
        {"bnd2.steiner_route_metadata_publication", testSteinerRouteMetadataPublication},
        {"bnd2.steiner_route_idempotent_chain", testSteinerRouteIdempotentChain},
        {"bnd2.steiner_route_full_session_preflight", testSteinerRouteFullSessionPreflight},
        {"bnd2.steiner_route_existing_class_conflict", testSteinerRouteExistingClassificationConflict},
        {"bnd2.steiner_route_auth_param_audit", testSteinerRouteAuthoritativeParamAudit},
        {"bnd2.steiner_route_provider_taxonomy", testSteinerRouteProviderTaxonomy},
        {"bnd2.steiner_route_policy_outside", testSteinerRoutePolicyAndOutside},
        {"bnd2.steiner_route_determinism", testSteinerRouteDeterminism},
        {"bnd2.steiner_multi_two_three", testSteinerMultiTwoThree},
        {"bnd2.steiner_multi_reversed", testSteinerMultiReversed},
        {"bnd2.steiner_multi_determinism", testSteinerMultiDeterminism},
        {"bnd2.steiner_multi_budget_feature", testSteinerMultiBudgetAndFeature},
        {"bnd2.steiner_multi_invalid_policy", testSteinerMultiInvalidPolicyProviderUntouched},
        {"bnd2.steiner_multi_concave_outside", testSteinerMultiConcaveAndOutside},
        {"bnd2.steiner_multi_idempotent", testSteinerMultiIdempotent},
        {"bnd2.steiner_multi_set_prefix", testSteinerMultiSetPrefix},
        {"bnd2.steiner_reuse", testSteinerReuse},
        {"bnd2.steiner_reuse_min_feature", testSteinerReuseMinFeature},
        {"bnd2.steiner_reuse_class_conflict", testSteinerReuseClassConflict},
        {"bnd2.steiner_reuse_coincident_split", testSteinerReuseCoincidentSplit},
        {"bnd2.steiner_reuse_outside", testSteinerReuseOutsideInvariance},
        {"bnd2.steiner_general_validator", testSteinerGeneralGeometryValidator},
        {"bnd2.steiner_strip_walk_guards", testSteinerStripWalkGuards},
        {"bnd2.steiner_exact_earclip", testSteinerExactEarClip},
        {"bnd2.ascr_rejections", testAscrRejections},
    };
    count = sizeof(k) / sizeof(k[0]);
    return k;
}
}  // namespace paacs_tests
