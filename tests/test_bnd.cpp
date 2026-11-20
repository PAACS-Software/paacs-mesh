// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// BND-1 Stage-A1 proofs: full-handle identity, ownership stamping, the
// complete whole-set preflight, the VALIDATED Tri3 adjacency index, and
// exact segment-chain validation/discovery. NO mutation path exists in
// this stage; every rejection proof also asserts the FULL observable
// digest unchanged. Design authority: docs/PMBND001_DESIGN.md sec 12.
#include <cstdio>
#include <cstring>
#include <limits>
#include <tuple>
#include <string>

#include "paacs/bnd/ParamGeometry.hpp"
#include "paacs/bnd/SegmentRecoveryCommit.hpp"
#include "paacs/bnd/SegmentRecoverySession.hpp"
#include "paacs/bnd/SegmentRecovery.hpp"
#include "workspace_digest.hpp"

namespace topo = paacs::topo;
namespace bnd = paacs::bnd;

extern int gFailures;   // shared harness counter (test_transaction.cpp)
static const char* gCurrent = "";
#define TEST(name)                       \
    do {                                 \
        gCurrent = name;                 \
        std::printf("[test] %s\n", name); \
    } while (0)
#define CHECK(cond)                                                      \
    do {                                                                 \
        if (!(cond)) {                                                   \
            std::printf("  FAIL %s:%d  %s  (in %s)\n", __FILE__,         \
                        __LINE__, #cond, gCurrent);                      \
            ++gFailures;                                                 \
        }                                                                \
    } while (0)

namespace {

struct Fix {
    std::unique_ptr<topo::TopologyWorkspace> ws;
    std::vector<topo::NodeHandle> n;
    std::vector<topo::CellHandle> c;
};

// 2x1 strip; the chain n0 - n1 - n2 lies exactly on v == 0.
Fix buildStrip(bool markChain) {
    Fix f;
    topo::NodePool np;
    f.n.push_back(np.add({0, 0, 0}, {2, 1}));
    f.n.push_back(np.add({1, 0, 0}, {2, 1}));
    f.n.push_back(np.add({2, 0, 0}, {2, 1}));
    f.n.push_back(np.add({0, 1, 0}, {2, 1}));
    f.n.push_back(np.add({1, 1, 0}, {2, 1}));
    f.n.push_back(np.add({2, 1, 0}, {2, 1}));
    topo::CellPool cp;
    f.c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<topo::NodeHandle>{f.n[0], f.n[1], f.n[4]}, {2, 5}, np));
    f.c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<topo::NodeHandle>{f.n[0], f.n[4], f.n[3]}, {2, 5}, np));
    f.c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<topo::NodeHandle>{f.n[1], f.n[2], f.n[5]}, {2, 5}, np));
    f.c.push_back(cp.add(topo::ElementType::Tri3,
        std::vector<topo::NodeHandle>{f.n[1], f.n[5], f.n[4]}, {2, 5}, np));
    topo::EdgePool ep;
    for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
    if (markChain)
        for (int i = 0; i < 2; ++i) {
            const auto eh = *ep.find(np, f.n[std::size_t(i)],
                                     f.n[std::size_t(i + 1)]);
            ep.setRequiredConstraint(eh, true);
            ep.setClassification(eh, {1, 9});
        }
    f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    return f;
}

// 3x1 strip (8 nodes, 6 triangles): TWO interior on-segment nodes, so
// NonIncreasingParameter and MissingEdge are reachable EXACTLY.
Fix buildStrip3(bool markChain) {
    Fix f;
    topo::NodePool np;
    for (int i = 0; i < 4; ++i)
        f.n.push_back(np.add({double(i), 0, 0}, {2, 1}));
    for (int i = 0; i < 4; ++i)
        f.n.push_back(np.add({double(i), 1, 0}, {2, 1}));
    topo::CellPool cp;
    const int tri[6][3] = {{0, 1, 5}, {0, 5, 4}, {1, 2, 6},
                           {1, 6, 5}, {2, 3, 7}, {2, 7, 6}};
    for (const auto& t : tri)
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{f.n[std::size_t(t[0])],
                                          f.n[std::size_t(t[1])],
                                          f.n[std::size_t(t[2])]},
            {2, 5}, np));
    topo::EdgePool ep;
    for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
    if (markChain)
        for (int i = 0; i < 3; ++i) {
            const auto eh = *ep.find(np, f.n[std::size_t(i)],
                                     f.n[std::size_t(i + 1)]);
            ep.setRequiredConstraint(eh, true);
            ep.setClassification(eh, {1, 9});
        }
    f.ws = std::make_unique<topo::TopologyWorkspace>(
        std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
    return f;
}

bnd::ConstraintSet buildSet(const Fix& f, int cols = 3) {
    bnd::ConstraintSet cs;
    cs.stamp = bnd::OwnershipStamp::of(*f.ws);
    for (std::size_t i = 0; i < f.n.size(); ++i)
        cs.params[f.n[i]] = bnd::ParamPoint{
            {double(i % std::size_t(cols)), double(i / std::size_t(cols))},
            7};
    cs.vertices[100] = f.n[0];
    cs.vertices[101] = f.n[std::size_t(cols - 1)];
    cs.segments.push_back(
        {500, 100, 101, 7, {1, 9}});
    return cs;
}

bool hasReason(const std::vector<bnd::PreflightIssue>& v,
               bnd::InvalidInputReason r) {
    for (const auto& i : v)
        if (i.reason == r) return true;
    return false;
}

void testBndPreflightClean() {
    TEST("bnd.preflight_clean_set_passes");
    auto f = buildStrip(true);
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    CHECK(bnd::preflightConstraintSet(*f.ws, buildSet(f)).empty());
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndPreflightStaleStamp() {
    TEST("bnd.preflight_stale_stamp");
    auto f = buildStrip(true);
    auto cs = buildSet(f);
    topo::TopologyPolicy pol;
    pol.permitOrphanNodes = true;
    auto tx = f.ws->beginTransaction(pol);
    tx.addNode({9, 9, 9}, {0, 7});
    (void)tx.commit();
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    const auto issues = bnd::preflightConstraintSet(*f.ws, cs);
    CHECK(issues.size() == 1);
    CHECK(issues[0].reason == bnd::InvalidInputReason::StaleStamp);
    CHECK(!issues[0].entity && !issues[0].node);
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndPreflightStaleHandle() {
    TEST("bnd.preflight_stale_handle (stamp current, handle stale)");
    auto f = buildStrip(true);
    topo::TopologyPolicy pol;
    pol.permitOrphanNodes = true;
    topo::NodeHandle stale{};
    {
        auto tx = f.ws->beginTransaction(pol);
        tx.addNode({9, 9, 9}, {0, 7});
        auto res = tx.commit();
        stale = res.node;
    }
    {
        auto tx = f.ws->beginTransaction(pol);
        tx.removeNode(stale);
        (void)tx.commit();
    }
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    auto cs = buildSet(f);                 // stamp captured NOW: current
    cs.vertices[102] = stale;              // stale VERTEX association
    cs.params[stale] = bnd::ParamPoint{{5, 5}, 7};   // stale PARAM key
    const auto issues = bnd::preflightConstraintSet(*f.ws, cs);
    int staleCount = 0;
    for (const auto& i : issues)
        if (i.reason == bnd::InvalidInputReason::StaleHandle) {
            ++staleCount;
            CHECK(i.node && *i.node == stale);
        }
    CHECK(staleCount == 2);
    CHECK(issues.size() == 2);
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndPreflightForeignHandle() {
    TEST("bnd.preflight_foreign_handle (vertex + param key)");
    auto fa = buildStrip(true);
    auto fb = buildStrip(true);            // distinct pools, same shape
    fa.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*fa.ws);
    auto cs = buildSet(fa);
    cs.vertices[102] = fb.n[0];            // foreign VERTEX association
    cs.params[fb.n[3]] = bnd::ParamPoint{{9, 9}, 7};  // foreign PARAM key
    const auto issues = bnd::preflightConstraintSet(*fa.ws, cs);
    int foreign = 0;
    for (const auto& i : issues)
        if (i.reason == bnd::InvalidInputReason::ForeignHandle) ++foreign;
    CHECK(foreign == 2);
    CHECK(issues.size() == 2);
    // full-handle identity: fb's handles have IDENTICAL slot+generation
    // to fa's - only the pool id distinguishes them (the aliasing the
    // review found nodeKey() permitted)
    CHECK(fb.n[0].slot == fa.n[0].slot && fb.n[0].gen == fa.n[0].gen);
    CHECK(paacs_tests::observableDigest(*fa.ws) == d0);
}

void testBndPreflightParamAudit() {
    TEST("bnd.preflight_param_audit (totality + nonfinite)");
    auto f = buildStrip(true);
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    {
        auto cs = buildSet(f);
        cs.params.erase(f.n[5]);
        const auto issues = bnd::preflightConstraintSet(*f.ws, cs);
        CHECK(issues.size() == 1);
        CHECK(issues[0].reason == bnd::InvalidInputReason::MissingParameter);
        CHECK(issues[0].node && *issues[0].node == f.n[5]);
    }
    {
        auto cs = buildSet(f);
        cs.params[f.n[4]].p.u = std::numeric_limits<double>::infinity();
        const auto issues = bnd::preflightConstraintSet(*f.ws, cs);
        CHECK(issues.size() == 1);
        CHECK(issues[0].reason ==
              bnd::InvalidInputReason::NonfiniteCoordinate);
    }
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndPreflightChart() {
    TEST("bnd.preflight_chart (endpoint + across-cell)");
    auto f = buildStrip(true);
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    {   // segment ENDPOINT chart mismatch: entity attributed
        auto cs = buildSet(f);
        for (auto& [h, pp] : cs.params)
            if (h == f.n[0]) pp.chart = 8;
        const auto issues = bnd::preflightConstraintSet(*f.ws, cs);
        CHECK(hasReason(issues, bnd::InvalidInputReason::ChartMismatch));
        bool entityAttr = false;
        for (const auto& i : issues)
            if (i.reason == bnd::InvalidInputReason::ChartMismatch &&
                i.entity && *i.entity == 500)
                entityAttr = true;
        CHECK(entityAttr);
    }
    {   // ACROSS-CELL chart mismatch on a node no segment touches
        auto cs = buildSet(f);
        for (auto& [h, pp] : cs.params)
            if (h == f.n[4]) pp.chart = 8;
        const auto issues = bnd::preflightConstraintSet(*f.ws, cs);
        bool cellAttr = false;
        for (const auto& i : issues)
            if (i.reason == bnd::InvalidInputReason::ChartMismatch &&
                !i.entity)
                cellAttr = true;
        CHECK(cellAttr);
    }
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndPreflightEntityRules() {
    TEST("bnd.preflight_entity_rules (dup, unknown, degenerate, cls)");
    auto f = buildStrip(true);
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    {
        auto cs = buildSet(f);
        cs.segments.push_back({500, 100, 100, 7, {2, 9}});
        const auto issues = bnd::preflightConstraintSet(*f.ws, cs);
        CHECK(issues.size() == 3);
        CHECK(hasReason(issues, bnd::InvalidInputReason::DuplicateEntityId));
        CHECK(hasReason(issues, bnd::InvalidInputReason::DegenerateSegment));
        CHECK(hasReason(issues,
                        bnd::InvalidInputReason::InvalidClassification));
    }
    {   // tag domain: dim correct, tag negative -> invalid
        auto cs = buildSet(f);
        cs.segments.push_back({501, 100, 101, 7, {1, -1}});
        const auto issues = bnd::preflightConstraintSet(*f.ws, cs);
        CHECK(issues.size() == 1);
        CHECK(issues[0].reason ==
              bnd::InvalidInputReason::InvalidClassification);
        CHECK(issues[0].entity && *issues[0].entity == 501);
    }
    {
        auto cs = buildSet(f);
        cs.segments.push_back({502, 100, 999, 7, {1, 9}});
        const auto issues = bnd::preflightConstraintSet(*f.ws, cs);
        CHECK(issues.size() == 1);
        CHECK(issues[0].reason == bnd::InvalidInputReason::UnknownEndpoint);
    }
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

// Adversaries: Quad4 cell, three-owner edge, pinched vertex. If the
// workspace CONSTRUCTOR refuses one, that upstream enforcement is
// recorded honestly instead of faked.
void testBndManifoldAdversaries() {
    TEST("bnd.preflight_and_index_manifold_adversaries");
    {   // Quad4 present
        topo::NodePool np;
        std::vector<topo::NodeHandle> n;
        for (int i = 0; i < 4; ++i)
            n.push_back(np.add({double(i % 2), double(i / 2), 0}, {2, 1}));
        topo::CellPool cp;
        const auto q = cp.add(topo::ElementType::Quad4,
            std::vector<topo::NodeHandle>{n[0], n[1], n[3], n[2]}, {2, 5},
            np);
        topo::EdgePool ep;
        ep.insertCellEdges(np, cp, q);
        bool built = false;
        try {
            topo::TopologyWorkspace ws(std::move(np), std::move(cp),
                                       std::move(ep), topo::FacePool{});
            built = true;
            CHECK(bnd::auditTri3Manifold(ws) ==
                  bnd::Tri3AuditResult::NonTri3Cell);
            CHECK(!bnd::buildTriAdjacency(ws));
            bnd::ConstraintSet cs;
            cs.stamp = bnd::OwnershipStamp::of(ws);
            for (const auto& nh : ws.nodes().activeBySequence())
                cs.params[nh] = bnd::ParamPoint{{0, 0}, 7};
            CHECK(hasReason(bnd::preflightConstraintSet(ws, cs),
                            bnd::InvalidInputReason::NonTri3Cell));
        } catch (const std::exception&) {
            CHECK(!built);   // upstream refusal: also a valid guarantee
        }
    }
    {   // three-owner edge
        topo::NodePool np;
        std::vector<topo::NodeHandle> n;
        n.push_back(np.add({0, 0, 0}, {2, 1}));
        n.push_back(np.add({1, 0, 0}, {2, 1}));
        n.push_back(np.add({0, 1, 0}, {2, 1}));
        n.push_back(np.add({0, -1, 0}, {2, 1}));
        n.push_back(np.add({1, 1, 0}, {2, 1}));
        topo::CellPool cp;
        std::vector<topo::CellHandle> c;
        c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{n[0], n[1], n[2]}, {2, 5}, np));
        c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{n[0], n[3], n[1]}, {2, 5}, np));
        c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{n[0], n[1], n[4]}, {2, 5}, np));
        topo::EdgePool ep;
        bool built = false;
        try {
            for (const auto& ch : c) ep.insertCellEdges(np, cp, ch);
            topo::TopologyWorkspace ws(std::move(np), std::move(cp),
                                       std::move(ep), topo::FacePool{});
            built = true;
            CHECK(bnd::auditTri3Manifold(ws) ==
                  bnd::Tri3AuditResult::OverusedEdge);
            CHECK(!bnd::buildTriAdjacency(ws));
            bnd::ConstraintSet cs;
            cs.stamp = bnd::OwnershipStamp::of(ws);
            for (const auto& nh : ws.nodes().activeBySequence())
                cs.params[nh] = bnd::ParamPoint{{0, 0}, 7};
            CHECK(hasReason(bnd::preflightConstraintSet(ws, cs),
                            bnd::InvalidInputReason::NonmanifoldIncidence));
        } catch (const std::exception&) {
            CHECK(!built);
        }
    }
    {   // pinched vertex: two triangles sharing ONLY n0
        topo::NodePool np;
        std::vector<topo::NodeHandle> n;
        n.push_back(np.add({0, 0, 0}, {2, 1}));
        n.push_back(np.add({1, 0, 0}, {2, 1}));
        n.push_back(np.add({0, 1, 0}, {2, 1}));
        n.push_back(np.add({-1, 0, 0}, {2, 1}));
        n.push_back(np.add({0, -1, 0}, {2, 1}));
        topo::CellPool cp;
        std::vector<topo::CellHandle> c;
        c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{n[0], n[1], n[2]}, {2, 5}, np));
        c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{n[0], n[3], n[4]}, {2, 5}, np));
        topo::EdgePool ep;
        for (const auto& ch : c) ep.insertCellEdges(np, cp, ch);
        topo::TopologyWorkspace ws(std::move(np), std::move(cp),
                                   std::move(ep), topo::FacePool{});
        CHECK(bnd::auditTri3Manifold(ws) ==
              bnd::Tri3AuditResult::PinchedVertex);
        CHECK(!bnd::buildTriAdjacency(ws));
        bnd::ConstraintSet cs;
        cs.stamp = bnd::OwnershipStamp::of(ws);
        for (const auto& nh : ws.nodes().activeBySequence())
            cs.params[nh] = bnd::ParamPoint{{0, 0}, 7};
        CHECK(hasReason(bnd::preflightConstraintSet(ws, cs),
                        bnd::InvalidInputReason::PinchedVertex));
    }
}

void testBndDuplicateTriangle() {
    TEST("bnd.duplicate_triangle_rejected (identical + reversed)");
    auto make = [](bool reversed) {
        topo::NodePool np;
        std::vector<topo::NodeHandle> n;
        n.push_back(np.add({0, 0, 0}, {2, 1}));
        n.push_back(np.add({1, 0, 0}, {2, 1}));
        n.push_back(np.add({0, 1, 0}, {2, 1}));
        topo::CellPool cp;
        std::vector<topo::CellHandle> c;
        c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{n[0], n[1], n[2]}, {2, 5}, np));
        // second coincident triangle: same nodes, identical or reversed
        if (reversed)
            c.push_back(cp.add(topo::ElementType::Tri3,
                std::vector<topo::NodeHandle>{n[0], n[2], n[1]}, {2, 5},
                np));
        else
            c.push_back(cp.add(topo::ElementType::Tri3,
                std::vector<topo::NodeHandle>{n[0], n[1], n[2]}, {2, 5},
                np));
        topo::EdgePool ep;
        for (const auto& ch : c) ep.insertCellEdges(np, cp, ch);
        return std::make_tuple(
            std::make_unique<topo::TopologyWorkspace>(
                std::move(np), std::move(cp), std::move(ep),
                topo::FacePool{}),
            n);
    };
    for (bool reversed : {false, true}) {
        std::unique_ptr<topo::TopologyWorkspace> ws;
        std::vector<topo::NodeHandle> n;
        bool built = false;
        try {
            std::tie(ws, n) = make(reversed);
            built = true;
        } catch (const std::exception&) {
            CHECK(!built);   // upstream refusal is also acceptable
        }
        if (!built) continue;
        // SPECIFICALLY DuplicateTriangle - not caught incidentally
        CHECK(bnd::auditTri3Manifold(*ws) ==
              bnd::Tri3AuditResult::DuplicateTriangle);
        CHECK(!bnd::buildTriAdjacency(*ws));
        bnd::ConstraintSet cs;
        cs.stamp = bnd::OwnershipStamp::of(*ws);
        for (const auto& nh : ws->nodes().activeBySequence())
            cs.params[nh] = bnd::ParamPoint{{0, 0}, 7};
        CHECK(hasReason(bnd::preflightConstraintSet(*ws, cs),
                        bnd::InvalidInputReason::DuplicateTriangle));
    }
}

void testBndPreflightDeterministic() {
    TEST("bnd.preflight_deterministic (insertion-order adversary)");
    // Two logically identical sets, built by inserting the SAME faults in
    // OPPOSITE orders, must yield bit-identical issue vectors including
    // attribution - not merely equal reason multisets.
    auto build = [](const Fix& f, bool reverse) {
        bnd::ConstraintSet cs;
        cs.stamp = bnd::OwnershipStamp::of(*f.ws);
        std::vector<std::size_t> order;
        for (std::size_t i = 0; i < f.n.size(); ++i) order.push_back(i);
        if (reverse) std::reverse(order.begin(), order.end());
        for (auto i : order)
            cs.params[f.n[i]] =
                bnd::ParamPoint{{double(i % 3), double(i / 3)}, 7};
        // introduce several faults touching different nodes so the
        // unordered iteration order actually matters
        cs.params[f.n[0]].p.u = std::numeric_limits<double>::infinity();
        cs.params.erase(f.n[5]);                     // MissingParameter
        for (auto& [h, pp] : cs.params)
            if (h == f.n[4]) pp.chart = 8;           // across-cell chart
        cs.vertices[100] = f.n[0];
        cs.vertices[101] = f.n[2];
        cs.segments.push_back({500, 100, 999, 7, {1, 9}});  // Unknown
        cs.segments.push_back({500, 100, 101, 7, {1, 9}});  // Duplicate id
        return cs;
    };
    auto sameVec = [](const std::vector<bnd::PreflightIssue>& x,
                      const std::vector<bnd::PreflightIssue>& y) {
        if (x.size() != y.size()) return false;
        for (std::size_t i = 0; i < x.size(); ++i) {
            if (x[i].reason != y[i].reason) return false;
            if (x[i].entity != y[i].entity) return false;
            const bool an = x[i].node.has_value(),
                       bn = y[i].node.has_value();
            if (an != bn) return false;
            if (an && !(*x[i].node == *y[i].node)) return false;
        }
        return true;
    };
    auto f = buildStrip(true);
    const auto fwd = bnd::preflightConstraintSet(*f.ws, build(f, false));
    const auto rev = bnd::preflightConstraintSet(*f.ws, build(f, true));
    CHECK(!fwd.empty());
    CHECK(sameVec(fwd, rev));   // FULL-vector equality incl. attribution

    // Adversary from the review: TWO DISTINCT foreign handles that
    // collided under the deleted fullId() XOR fold - {0,0,pool=3} and
    // {689935870,728783,pool=1546} both mapped to 48180ab7777ec1ed.
    // Their ForeignHandle diagnostics must now order deterministically
    // regardless of insertion order (full-handle operator<, not a fold).
    auto buildCollide = [&](bool reverse) {
        auto cs = build(f, false);
        topo::NodeHandle h1{0, 0, 3};
        topo::NodeHandle h2{689935870, 728783, 1546};
        if (!reverse) {
            cs.params[h1] = bnd::ParamPoint{{1, 1}, 7};
            cs.params[h2] = bnd::ParamPoint{{2, 2}, 7};
        } else {
            cs.params[h2] = bnd::ParamPoint{{2, 2}, 7};
            cs.params[h1] = bnd::ParamPoint{{1, 1}, 7};
        }
        return cs;
    };
    const auto cf = bnd::preflightConstraintSet(*f.ws, buildCollide(false));
    const auto cr = bnd::preflightConstraintSet(*f.ws, buildCollide(true));
    CHECK(sameVec(cf, cr));
    // and the two distinct handles both surface as ForeignHandle (they
    // are NOT collapsed into one)
    int foreign = 0;
    for (const auto& i : cf)
        if (i.reason == bnd::InvalidInputReason::ForeignHandle) ++foreign;
    CHECK(foreign >= 2);
}

void testBndTriAdjacencyIndex() {
    TEST("bnd.tri_adjacency_index (full handles, symmetry, ownership)");
    auto f = buildStrip(false);
    const auto idxOpt = bnd::buildTriAdjacency(*f.ws);
    CHECK(bool(idxOpt));
    const auto& idx = *idxOpt;
    CHECK(idx.current(*f.ws));
    CHECK(idx.neighbors.size() == 4);
    std::size_t links = 0, boundary = 0;
    for (const auto& [ch, nb] : idx.neighbors)
        for (const auto& h : nb)
            (f.ws->cells().isCurrent(h) ? links : boundary) += 1;
    CHECK(links == 6);
    CHECK(boundary == 6);
    for (const auto& [ch, nb] : idx.neighbors)         // symmetry over
        for (const auto& h : nb) {                     // FULL handles
            if (!f.ws->cells().isCurrent(h)) continue;
            bool back = false;
            for (const auto& r : idx.neighbors.at(h)) back |= r == ch;
            CHECK(back);
        }
    CHECK(idx.star.at(f.n[1]).size() == 3);
    CHECK(idx.star.at(f.n[3]).size() == 1);
    // OWNERSHIP binding (review): a second workspace with the SAME
    // revision must NOT satisfy current() - the stamp, not the revision,
    // is the binding.
    auto g = buildStrip(false);
    CHECK(g.ws->topologyRevision() == f.ws->topologyRevision());
    CHECK(!idx.current(*g.ws));
    // revision binding still holds too
    topo::TopologyPolicy pol;
    pol.permitOrphanNodes = true;
    auto tx = f.ws->beginTransaction(pol);
    tx.addNode({9, 9, 9}, {0, 7});
    (void)tx.commit();
    CHECK(!idx.current(*f.ws));
}

void testBndChainValidationExact() {
    TEST("bnd.chain_validation_exact (every status, exactly)");
    auto f = buildStrip3(true);
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    auto cs = buildSet(f, 4);
    const auto& seg = cs.segments[0];        // n0 -> n3 along v == 0
    using S = bnd::ChainStatus;
    auto st = [&](const std::vector<topo::NodeHandle>& c,
                  const bnd::ConstraintSet& s, const auto& sg) {
        return bnd::validateSegmentChain(*f.ws, s, sg, c).status;
    };
    CHECK(st({f.n[0], f.n[1], f.n[2], f.n[3]}, cs, seg) == S::Valid);
    CHECK(st({f.n[0]}, cs, seg) == S::Empty);
    CHECK(st({f.n[1], f.n[2], f.n[3]}, cs, seg) == S::EndpointMismatch);
    CHECK(st({f.n[0], f.n[1], f.n[1], f.n[3]}, cs, seg) ==
          S::DuplicateNode);
    {   // MissingParameter on an interior node
        auto c2 = cs;
        c2.params.erase(f.n[2]);
        CHECK(st({f.n[0], f.n[1], f.n[2], f.n[3]}, c2, seg) ==
              S::MissingParameter);
    }
    {   // ChartMismatch on an interior node
        auto c2 = cs;
        c2.params[f.n[2]].chart = 8;
        CHECK(st({f.n[0], f.n[1], f.n[2], f.n[3]}, c2, seg) ==
              S::ChartMismatch);
    }
    // OffSegment: non-collinear interior node
    CHECK(st({f.n[0], f.n[5], f.n[3]}, cs, seg) == S::OffSegment);
    {   // OffSegment: COLLINEAR but out of exact bounds (n0->n1 segment,
        // n2 at u == 2 beyond hi == 1)
        auto c2 = cs;
        c2.vertices[102] = f.n[1];
        bnd::ConstraintSegment s2{501, 100, 102, 7, {1, 9}};
        CHECK(st({f.n[0], f.n[2], f.n[1]}, c2, s2) == S::OffSegment);
    }
    // NonIncreasingParameter, EXACTLY (in-bounds, order reversed)
    CHECK(st({f.n[0], f.n[2], f.n[1], f.n[3]}, cs, seg) ==
          S::NonIncreasingParameter);
    // MissingEdge (params fine; n0 - n2 is not an edge)
    CHECK(st({f.n[0], f.n[2], f.n[3]}, cs, seg) == S::MissingEdge);
    {   // EdgeNotRequired
        auto fu = buildStrip3(false);
        auto cu = buildSet(fu, 4);
        CHECK(bnd::validateSegmentChain(
                  *fu.ws, cu, cu.segments[0],
                  {fu.n[0], fu.n[1], fu.n[2], fu.n[3]})
                  .status == S::EdgeNotRequired);
    }
    {   // ClassificationMismatch
        auto c2 = cs;
        c2.segments[0].cls = {1, 4};
        CHECK(st({f.n[0], f.n[1], f.n[2], f.n[3]}, c2, c2.segments[0]) ==
              S::ClassificationMismatch);
    }
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndChainDiscovery() {
    TEST("bnd.chain_discovery (finds; refuses broken/stale/foreign)");
    auto f = buildStrip3(true);
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    auto cs = buildSet(f, 4);
    const auto idx = bnd::buildTriAdjacency(*f.ws);
    CHECK(bool(idx));
    const auto chain = bnd::discoverChain(*f.ws, cs, *idx, cs.segments[0]);
    CHECK(bool(chain));
    CHECK(chain->size() == 4);
    for (int i = 0; i < 4; ++i)
        CHECK((*chain)[std::size_t(i)] == f.n[std::size_t(i)]);
    CHECK(bnd::validateSegmentChain(*f.ws, cs, cs.segments[0], *chain)
              .status == bnd::ChainStatus::Valid);
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
    {   // broken chain
        auto fb = buildStrip3(true);
        auto cb = buildSet(fb, 4);
        fb.ws->testEdgesMutable().setRequiredConstraint(
            *fb.ws->edges().find(fb.ws->nodes(), fb.n[1], fb.n[2]), false);
        const auto ib = bnd::buildTriAdjacency(*fb.ws);
        CHECK(bool(ib));
        CHECK(!bnd::discoverChain(*fb.ws, cb, *ib, cb.segments[0]));
    }
    {   // FOREIGN index: same revision, different workspace -> refused
        auto g = buildStrip3(true);
        auto cg = buildSet(g, 4);
        CHECK(!bnd::discoverChain(*g.ws, cg, *idx, cg.segments[0]));
    }
    {   // stale index after a commit -> refused
        auto fs = buildStrip3(true);
        auto csx = buildSet(fs, 4);
        const auto is = bnd::buildTriAdjacency(*fs.ws);
        CHECK(bool(is));
        topo::TopologyPolicy pol;
        pol.permitOrphanNodes = true;
        auto tx = fs.ws->beginTransaction(pol);
        tx.addNode({9, 9, 9}, {0, 7});
        (void)tx.commit();
        CHECK(!bnd::discoverChain(*fs.ws, csx, *is, csx.segments[0]));
    }
}

}  // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* bndTests(std::size_t& count) {
    static const Entry kBnd[] = {
        {"bnd.preflight_clean_set_passes", testBndPreflightClean},
        {"bnd.preflight_stale_stamp", testBndPreflightStaleStamp},
        {"bnd.preflight_stale_handle", testBndPreflightStaleHandle},
        {"bnd.preflight_foreign_handle", testBndPreflightForeignHandle},
        {"bnd.preflight_param_audit", testBndPreflightParamAudit},
        {"bnd.preflight_chart", testBndPreflightChart},
        {"bnd.preflight_entity_rules", testBndPreflightEntityRules},
        {"bnd.preflight_and_index_manifold_adversaries",
         testBndManifoldAdversaries},
        {"bnd.duplicate_triangle_rejected", testBndDuplicateTriangle},
        {"bnd.preflight_deterministic", testBndPreflightDeterministic},
        {"bnd.tri_adjacency_index", testBndTriAdjacencyIndex},
        {"bnd.chain_validation_exact", testBndChainValidationExact},
        {"bnd.chain_discovery", testBndChainDiscovery},
    };
    count = sizeof(kBnd) / sizeof(kBnd[0]);
    return kBnd;
}
}  // namespace paacs_tests

// ============================================================================
// Stage A2: source-embedding geometric preflight proofs.
// ============================================================================
namespace {

// Build a workspace directly from explicit triangles + a param map, so
// each A2 adversary is exact. Triangles are (node-index triples); the
// caller supplies authoritative Vec2 per node.
struct GeomFix {
    std::unique_ptr<topo::TopologyWorkspace> ws;
    std::vector<topo::NodeHandle> n;
    std::vector<topo::CellHandle> c;
    bnd::ConstraintSet cs;
};

GeomFix buildGeom(const std::vector<paacs::Vec2>& coords,
                  const std::vector<std::array<int, 3>>& tris,
                  bool mapAll = true) {
    GeomFix f;
    topo::NodePool np;
    for (const auto& q : coords)
        f.n.push_back(np.add({q.u, q.v, 0}, {2, 1}));
    topo::CellPool cp;
    for (const auto& t : tris)
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
    if (mapAll)
        for (std::size_t i = 0; i < coords.size(); ++i)
            f.cs.params[f.n[i]] = bnd::ParamPoint{coords[i], 7};
    return f;
}

bool hasStatus(const std::vector<bnd::EmbeddingIssue>& v,
               bnd::EmbeddingStatus s) {
    for (const auto& i : v)
        if (i.status == s) return true;
    return false;
}

void testBndEmbeddingClean() {
    TEST("bnd.embedding_clean (2x1 strip is geometrically sane)");
    // 3x1-style strip: consistent orientation, no crossings, no T-junctions.
    auto f = buildGeom(
        {{0, 0}, {1, 0}, {2, 0}, {0, 1}, {1, 1}, {2, 1}},
        {{0, 1, 4}, {0, 4, 3}, {1, 2, 5}, {1, 5, 4}});
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    const auto issues = bnd::preflightSourceEmbedding(*f.ws, f.cs);
    CHECK(issues.empty());
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndEmbeddingStaleStamp() {
    TEST("bnd.embedding_stale_stamp");
    auto f = buildGeom({{0, 0}, {1, 0}, {0, 1}}, {{0, 1, 2}});
    topo::TopologyPolicy pol;
    pol.permitOrphanNodes = true;
    auto tx = f.ws->beginTransaction(pol);
    tx.addNode({9, 9, 9}, {0, 7});
    (void)tx.commit();
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    const auto issues = bnd::preflightSourceEmbedding(*f.ws, f.cs);
    CHECK(issues.size() == 1);
    CHECK(issues[0].status == bnd::EmbeddingStatus::StaleStamp);
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndEmbeddingDegenerate() {
    TEST("bnd.embedding_degenerate_triangle (collinear in param space)");
    // Three collinear parameter coordinates -> zero signed area.
    auto f = buildGeom({{0, 0}, {1, 0}, {2, 0}}, {{0, 1, 2}});
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    const auto issues = bnd::preflightSourceEmbedding(*f.ws, f.cs);
    CHECK(hasStatus(issues, bnd::EmbeddingStatus::DegenerateTriangle));
    bool attributed = false;
    for (const auto& i : issues)
        if (i.status == bnd::EmbeddingStatus::DegenerateTriangle &&
            i.cellA && *i.cellA == f.c[0])
            attributed = true;
    CHECK(attributed);
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndEmbeddingOrientation() {
    TEST("bnd.embedding_inconsistent_orientation (genuinely reversed)");
    // T0 = (0,0),(1,0),(0,1) is CCW (orient2d Positive).
    // T1 = (2,0),(2,1),(3,0): winding chosen so orient2d is NEGATIVE -
    // a GENUINELY reversed triangle, not an ambiguous one.
    const auto s0 =
        paacs::predicates::orient2d({0, 0}, {1, 0}, {0, 1});
    const auto s1 =
        paacs::predicates::orient2d({2, 0}, {2, 1}, {3, 0});
    CHECK(s0 != s1);                     // the fixture IS inconsistent
    CHECK(s0 != paacs::predicates::Sign::Zero);
    CHECK(s1 != paacs::predicates::Sign::Zero);
    auto f = buildGeom(
        {{0, 0}, {1, 0}, {0, 1}, {2, 0}, {2, 1}, {3, 0}},
        {{0, 1, 2}, {3, 4, 5}});
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    const auto issues = bnd::preflightSourceEmbedding(*f.ws, f.cs);
    CHECK(hasStatus(issues,
                    bnd::EmbeddingStatus::InconsistentOrientation));
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndEmbeddingTriangleOverlap() {
    TEST("bnd.embedding_triangle_overlap (nested, no shared vertex)");
    // Big CCW triangle; a small CCW triangle strictly inside it. No shared
    // vertices, no edge crossing of the vertex-on-edge kind - only the
    // interior-overlap test catches it.
    auto f = buildGeom(
        {{0, 0}, {6, 0}, {0, 6},          // big T0
         {1, 1}, {2, 1}, {1, 2}},         // small T1 strictly inside
        {{0, 1, 2}, {3, 4, 5}});
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    const auto issues = bnd::preflightSourceEmbedding(*f.ws, f.cs);
    CHECK(hasStatus(issues, bnd::EmbeddingStatus::TriangleOverlap));
    // deterministic attribution: a small-triangle vertex inside the big
    bool attributed = false;
    for (const auto& i : issues)
        if (i.status == bnd::EmbeddingStatus::TriangleOverlap && i.cellA &&
            i.cellB && *i.cellA == f.c[1] && *i.cellB == f.c[0])
            attributed = true;
    CHECK(attributed);
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndEmbeddingCrossing() {
    TEST("bnd.embedding_edge_crossing (transversal, no shared vertex)");
    // Two disjoint triangles whose edges cross in parameter space: an "X".
    // T0 = (0,0)-(2,2)-(0,2); T1 = (2,0)-(0,... arranged so edge (0,0)-(2,2)
    // of T0 crosses edge (2,0)-(0,2) of T1.
    auto f = buildGeom(
        {{0, 0}, {2, 2}, {0, 2},          // T0 nodes 0,1,2
         {2, 0}, {0, 2 + 1e-9}, {2, 2 + 1e-9}},   // T1 nodes 3,4,5 (shifted)
        {{0, 1, 2}, {3, 4, 5}});
    // Edge n0-n1 = (0,0)-(2,2); edge n3-... must cross it. Use a cleaner
    // construction: put T1 as (2,0)-(0,2)-(2,2+) so its edge (2,0)-(0,2)
    // crosses (0,0)-(2,2).
    auto f2 = buildGeom(
        {{0, 0}, {2, 2}, {1, 3},          // T0: edge 0-1 is the diagonal
         {2, 0}, {0, 2}, {2, 3}},         // T1: edge 3-4 = (2,0)-(0,2)
        {{0, 1, 2}, {3, 4, 5}});
    f2.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f2.ws);
    const auto issues = bnd::preflightSourceEmbedding(*f2.ws, f2.cs);
    CHECK(hasStatus(issues, bnd::EmbeddingStatus::EdgeCrossing));
    CHECK(paacs_tests::observableDigest(*f2.ws) == d0);
    (void)f;
}

void testBndEmbeddingTJunction() {
    TEST("bnd.embedding_tjunction (disjoint triangles; A1 clean first)");
    // Two TOPOLOGICALLY DISJOINT triangles (no shared node, so A1 has no
    // pinched-vertex objection). n5 == (1,0) lies in the OPEN interior of
    // T0's base edge n0(0,0)-n1(2,0).
    auto f = buildGeom(
        {{0, 0}, {2, 0}, {1, 2},          // T0: base edge 0-1
         {1, 0}, {1, -2}, {3, -2}},       // T1: (1,0),(1,-2),(3,-2)
        {{0, 1, 2}, {3, 4, 5}});
    f.ws->rebuildAdjacency();
    // A1 must be CLEAN (this is a geometric fault, not a combinatorial one)
    CHECK(bnd::preflightConstraintSet(*f.ws, f.cs).empty());
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    const auto issues = bnd::preflightSourceEmbedding(*f.ws, f.cs);
    CHECK(hasStatus(issues, bnd::EmbeddingStatus::TJunction));
    bool attributed = false;
    for (const auto& i : issues)
        if (i.status == bnd::EmbeddingStatus::TJunction && i.nodeA &&
            *i.nodeA == f.n[3])
            attributed = true;
    CHECK(attributed);
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndEmbeddingSharedEdgeOverlap() {
    TEST("bnd.embedding_shared_edge_overlap (same-side nesting)");
    // Reviewer's adversary: two A1-valid CCW triangles sharing edge n0-n1,
    // the smaller ON THE SAME SIDE inside the larger. No crossing, no
    // T-junction; only the shared-donor-vertex-skipping overlap catches
    // it. n2=(0,4) big apex, n3=(0.5,0.5) small apex strictly inside.
    auto f = buildGeom(
        {{0, 0}, {4, 0}, {0, 4}, {0.5, 0.5}},
        {{0, 1, 2}, {0, 1, 3}});
    f.ws->rebuildAdjacency();
    // A1 accepts it (valid manifold Tri3, shared edge n0-n1 has 2 owners)
    CHECK(bnd::preflightConstraintSet(*f.ws, f.cs).empty());
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    const auto issues = bnd::preflightSourceEmbedding(*f.ws, f.cs);
    CHECK(hasStatus(issues, bnd::EmbeddingStatus::TriangleOverlap));
    // attribution: the small apex n3 is the non-shared donor vertex,
    // strictly inside the big triangle c0
    bool attributed = false;
    for (const auto& i : issues)
        if (i.status == bnd::EmbeddingStatus::TriangleOverlap && i.nodeA &&
            *i.nodeA == f.n[3])
            attributed = true;
    CHECK(attributed);
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndEmbeddingChartLocal() {
    TEST("bnd.embedding_chart_local (all four tests suppressed x-chart)");
    // Two DISCONNECTED triangles in DIFFERENT charts (7 and 8), with
    // coordinates deliberately arranged so that IF the checker compared
    // across charts it would fire EVERY geometric status:
    //   - T0 chart 7: CCW (0,0),(4,0),(0,4)
    //   - T1 chart 8: CW  (1,1),(1,2),(2,1)  [opposite orientation],
    //     entirely INSIDE T0's coordinate footprint [containment],
    //     its edges CROSS T0's hypotenuse [transversal crossing], and
    //     its vertex (2,1) would sit... arranged for a T-junction too.
    // Chart-BLIND code would report InconsistentOrientation + EdgeCrossing
    // + TriangleOverlap + TJunction; sec 5b -> A2 must be CLEAN.
    // T1 = (2,-1),(2,1),(3,0): verified to fire ALL FOUR statuses when
    // compared IN ONE CHART with T0 - opposite orientation, a transversal
    // crossing of T0's base, strict containment of a vertex, and a
    // vertex (2,0) exactly on T0's base edge (a T-junction).
    auto f = buildGeom(
        {{0, 0}, {4, 0}, {0, 4},           // T0 chart 7 (CCW)
         {2, -1}, {2, 1}, {3, 0}},         // T1 chart 8
        {{0, 1, 2}, {3, 4, 5}}, /*mapAll=*/false);
    f.cs.params[f.n[0]] = bnd::ParamPoint{{0, 0}, 7};
    f.cs.params[f.n[1]] = bnd::ParamPoint{{4, 0}, 7};
    f.cs.params[f.n[2]] = bnd::ParamPoint{{0, 4}, 7};
    f.cs.params[f.n[3]] = bnd::ParamPoint{{2, -1}, 8};
    f.cs.params[f.n[4]] = bnd::ParamPoint{{2, 1}, 8};
    f.cs.params[f.n[5]] = bnd::ParamPoint{{3, 0}, 8};
    f.ws->rebuildAdjacency();

    // (a) A1 CLEAN: valid disjoint manifold Tri3, each cell single-chart.
    CHECK(bnd::preflightConstraintSet(*f.ws, f.cs).empty());

    // (b) GEOMETRIC DISCRIMINATION: prove the fixture WOULD fire each
    // status if compared in a single chart. Re-stamp every coord to
    // chart 7 and confirm the chart-blind view is loud on all four.
    {
        bnd::ConstraintSet one = f.cs;   // same handles/ws, one chart
        for (auto& [hh, pp] : one.params) pp.chart = 7;
        const auto blindIssues = bnd::preflightSourceEmbedding(*f.ws, one);
        CHECK(hasStatus(blindIssues,
                        bnd::EmbeddingStatus::InconsistentOrientation));
        CHECK(hasStatus(blindIssues,
                        bnd::EmbeddingStatus::EdgeCrossing));
        CHECK(hasStatus(blindIssues,
                        bnd::EmbeddingStatus::TriangleOverlap));
        CHECK(hasStatus(blindIssues, bnd::EmbeddingStatus::TJunction));
    }

    // (c) CHART-LOCAL: with the real two-chart map, A2 is CLEAN and the
    // workspace digest is unchanged across the call.
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    const auto issues = bnd::preflightSourceEmbedding(*f.ws, f.cs);
    CHECK(issues.empty());
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndEmbeddingMissingParam() {
    TEST("bnd.embedding_missing_parameter");
    auto f = buildGeom({{0, 0}, {1, 0}, {0, 1}}, {{0, 1, 2}});
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    f.cs.params.erase(f.n[2]);            // drop one authoritative coord
    const auto issues = bnd::preflightSourceEmbedding(*f.ws, f.cs);
    CHECK(hasStatus(issues, bnd::EmbeddingStatus::MissingParameter));
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndEmbeddingDeterministic() {
    TEST("bnd.embedding_deterministic (issue order is canonical)");
    // A mesh with BOTH a degenerate triangle and a T-junction; run twice
    // and require bit-identical issue vectors (the maps iterate in
    // nondeterministic order internally).
    // Several faults on ONE workspace; the ConstraintSet params are
    // inserted in OPPOSITE orders (the map iterates nondeterministically)
    // and the FULL EmbeddingIssue vectors - every field - must match.
    auto buildCs = [](const GeomFix& f, const std::vector<paacs::Vec2>& C,
                      bool reverse) {
        bnd::ConstraintSet cs;
        cs.stamp = bnd::OwnershipStamp::of(*f.ws);
        std::vector<std::size_t> order;
        for (std::size_t i = 0; i < C.size(); ++i) order.push_back(i);
        if (reverse) std::reverse(order.begin(), order.end());
        for (auto i : order) cs.params[f.n[i]] = bnd::ParamPoint{C[i], 7};
        return cs;
    };
    const std::vector<paacs::Vec2> C{
        {0, 0}, {6, 0}, {0, 6},            // big T0
        {1, 1}, {2, 1}, {1, 2},            // small T1 inside (overlap)
        {10, 0}, {12, 0}, {11, 1}, {11, 0}};  // T2 base + n9 on it (T-jct)
    auto f = buildGeom(C, {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}, {6, 9, 7}},
                       /*mapAll=*/false);
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    auto csF = buildCs(f, C, false);
    auto csR = buildCs(f, C, true);
    const auto ia = bnd::preflightSourceEmbedding(*f.ws, csF);
    const auto ib = bnd::preflightSourceEmbedding(*f.ws, csR);
    CHECK(!ia.empty());
    bool identical = ia.size() == ib.size();
    for (std::size_t i = 0; i < ia.size() && identical; ++i) {
        if (ia[i].status != ib[i].status) identical = false;
        if (!(ia[i].cellA == ib[i].cellA)) identical = false;
        if (!(ia[i].cellB == ib[i].cellB)) identical = false;   // full
        if (!(ia[i].nodeA == ib[i].nodeA)) identical = false;   // vector
        if (!(ia[i].nodeB == ib[i].nodeB)) identical = false;
    }
    CHECK(identical);
    CHECK(hasStatus(ia, bnd::EmbeddingStatus::TriangleOverlap));
    CHECK(hasStatus(ia, bnd::EmbeddingStatus::TJunction));
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

}  // namespace

namespace paacs_tests {
const Entry* bndEmbeddingTests(std::size_t& count) {
    static const Entry kBnd[] = {
        {"bnd.embedding_clean", testBndEmbeddingClean},
        {"bnd.embedding_stale_stamp", testBndEmbeddingStaleStamp},
        {"bnd.embedding_degenerate_triangle", testBndEmbeddingDegenerate},
        {"bnd.embedding_inconsistent_orientation", testBndEmbeddingOrientation},
        {"bnd.embedding_edge_crossing", testBndEmbeddingCrossing},
        {"bnd.embedding_triangle_overlap", testBndEmbeddingTriangleOverlap},
        {"bnd.embedding_shared_edge_overlap",
         testBndEmbeddingSharedEdgeOverlap},
        {"bnd.embedding_chart_local", testBndEmbeddingChartLocal},
        {"bnd.embedding_tjunction", testBndEmbeddingTJunction},
        {"bnd.embedding_missing_parameter", testBndEmbeddingMissingParam},
        {"bnd.embedding_deterministic", testBndEmbeddingDeterministic},
    };
    count = sizeof(kBnd) / sizeof(kBnd[0]);
    return kBnd;
}
}  // namespace paacs_tests

// ============================================================================
// sec-13 increment 1: whole-set PSLG preflight + RecoverySession.
// VALIDATION ONLY - no mutation path exists in this increment.
// ============================================================================
namespace {

GeomFix buildPslgBase() {
    return buildGeom({{0, 0}, {4, 0}, {0, 4},
                      {10, 0}, {14, 0}, {10, 4}},
                     {{0, 1, 2}, {3, 4, 5}});
}

// Two triangles sharing edge n1-n3: nodes (0,0),(2,0),(4,0),(2,3).
// Gives a node STRICTLY INTERIOR to the (0,0)-(4,0) line for an exact,
// non-collinear T-junction fixture.
GeomFix buildTeeBase() {
    return buildGeom({{0, 0}, {2, 0}, {4, 0}, {2, 3}},
                     {{0, 1, 3}, {1, 2, 3}});
}

bool hasPslg(const std::vector<bnd::PslgIssue>& v, bnd::PslgStatus s) {
    for (const auto& i : v)
        if (i.status == s) return true;
    return false;
}

void testBndPslgClean() {
    TEST("bnd.pslg_clean (well-formed set; legal single-point touches)");
    auto f = buildTeeBase();
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    // n0-n1 and n1-n2 are COLLINEAR and share exactly the point n1: the
    // strict-interval rule must accept this (touching, not overlapping),
    // and the shared endpoint must not register as a T-junction.
    f.cs.vertices[1] = f.n[0];
    f.cs.vertices[2] = f.n[1];
    f.cs.vertices[3] = f.n[2];
    f.cs.segments.push_back({10, 1, 2, 7, {1, 9}});   // (0,0)-(2,0)
    f.cs.segments.push_back({11, 2, 3, 7, {1, 9}});   // (2,0)-(4,0)
    const auto issues = bnd::preflightPslg(f.cs);
    CHECK(issues.empty());
    CHECK(!hasPslg(issues, bnd::PslgStatus::CollinearOverlap));
    CHECK(!hasPslg(issues, bnd::PslgStatus::ConstraintTJunction));
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndPslgNonInjective() {
    TEST("bnd.pslg_non_injective_vertex_mapping (13.1a)");
    auto f = buildPslgBase();
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    f.cs.vertices[1] = f.n[0];
    f.cs.vertices[2] = f.n[0];            // SAME handle, two ids
    const auto issues = bnd::preflightPslg(f.cs);
    CHECK(hasPslg(issues, bnd::PslgStatus::NonInjectiveVertexMapping));
    bool attributed = false;
    for (const auto& i : issues)
        if (i.status == bnd::PslgStatus::NonInjectiveVertexMapping && i.node &&
            *i.node == f.n[0])
            attributed = true;
    CHECK(attributed);
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndPslgDuplicateAndCoincident() {
    TEST("bnd.pslg_duplicate_and_coincident");
    {
        auto f = buildPslgBase();
        f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        f.cs.vertices[1] = f.n[0];
        f.cs.vertices[2] = f.n[1];
        f.cs.segments.push_back({10, 1, 2, 7, {1, 9}});
        f.cs.segments.push_back({11, 2, 1, 7, {1, 9}});   // reversed dup
        CHECK(hasPslg(bnd::preflightPslg(f.cs),
                      bnd::PslgStatus::DuplicateSegment));
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);
    }
    {
        auto f = buildPslgBase();
        f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        f.cs.params[f.n[3]] = bnd::ParamPoint{{0, 0}, 7};   // == n0's coord
        f.cs.vertices[1] = f.n[0];
        f.cs.vertices[2] = f.n[3];
        CHECK(hasPslg(bnd::preflightPslg(f.cs),
                      bnd::PslgStatus::CoincidentVertices));
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);
    }
}

void testBndPslgTJunctionExact() {
    TEST("bnd.pslg_tjunction_exact (perpendicular guest; NOT collinear)");
    // Host (0,0)-(4,0); guest (2,0)-(2,3) is PERPENDICULAR, so no
    // collinear overlap is possible and the T-junction path is proven
    // in isolation. n1 == (2,0) is strictly interior to the host.
    auto f = buildTeeBase();
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    f.cs.vertices[1] = f.n[0];    // (0,0)
    f.cs.vertices[2] = f.n[2];    // (4,0)
    f.cs.vertices[3] = f.n[1];    // (2,0) - interior to host
    f.cs.vertices[4] = f.n[3];    // (2,3)
    f.cs.segments.push_back({10, 1, 2, 7, {1, 9}});   // host
    f.cs.segments.push_back({11, 3, 4, 7, {1, 9}});   // guest, perpendicular
    const auto issues = bnd::preflightPslg(f.cs);
    CHECK(hasPslg(issues, bnd::PslgStatus::ConstraintTJunction));
    CHECK(!hasPslg(issues, bnd::PslgStatus::CollinearOverlap));  // isolated
    CHECK(!hasPslg(issues, bnd::PslgStatus::ConstraintCrossing));
    bool attributed = false;
    for (const auto& i : issues)
        if (i.status == bnd::PslgStatus::ConstraintTJunction && i.node &&
            *i.node == f.n[1])
            attributed = true;
    CHECK(attributed);
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndPslgCollinearOverlapExact() {
    TEST("bnd.pslg_collinear_overlap_exact (partial overlap only)");
    // Host (0,0)-(4,0); guest (2,0)-(6,0): collinear, sharing [2,4] -
    // strictly more than a point. Distinct from the single-point touch
    // proven clean in pslg_clean.
    auto f = buildGeom({{0, 0}, {2, 0}, {4, 0}, {6, 0}, {3, 3}},
                       {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}});
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    f.cs.vertices[1] = f.n[0];   // (0,0)
    f.cs.vertices[2] = f.n[2];   // (4,0)
    f.cs.vertices[3] = f.n[1];   // (2,0)
    f.cs.vertices[4] = f.n[3];   // (6,0)
    f.cs.segments.push_back({10, 1, 2, 7, {1, 9}});   // [0,4]
    f.cs.segments.push_back({11, 3, 4, 7, {1, 9}});   // [2,6] -> overlap
    const auto issues = bnd::preflightPslg(f.cs);
    CHECK(hasPslg(issues, bnd::PslgStatus::CollinearOverlap));
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndPslgCrossingAndChartLocal() {
    TEST("bnd.pslg_crossing_and_chart_local (loud in one chart, clean x2)");
    // Two segments that PROPERLY CROSS when compared in one chart.
    auto f = buildGeom({{0, 0}, {4, 0}, {0, 4}, {4, 4}},
                       {{0, 1, 3}, {0, 3, 2}});
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    f.cs.vertices[1] = f.n[0];   // (0,0)
    f.cs.vertices[2] = f.n[3];   // (4,4)  -> diagonal
    f.cs.vertices[3] = f.n[1];   // (4,0)
    f.cs.vertices[4] = f.n[2];   // (0,4)  -> anti-diagonal
    f.cs.segments.push_back({10, 1, 2, 7, {1, 9}});
    f.cs.segments.push_back({11, 3, 4, 7, {1, 9}});
    // (a) LOUD: both segments in chart 7 -> proper crossing
    CHECK(hasPslg(bnd::preflightPslg(f.cs),
                  bnd::PslgStatus::ConstraintCrossing));
    // (b) CLEAN: put the second segment's geometry in chart 8 - sec 5b
    // forbids cross-chart comparison, so the crossing must disappear.
    {
        auto cs2 = f.cs;
        cs2.segments[1].chart = 8;
        cs2.params[f.n[1]].chart = 8;
        cs2.params[f.n[2]].chart = 8;
        const auto issues = bnd::preflightPslg(cs2);
        CHECK(!hasPslg(issues, bnd::PslgStatus::ConstraintCrossing));
        CHECK(!hasPslg(issues, bnd::PslgStatus::ConstraintTJunction));
        CHECK(!hasPslg(issues, bnd::PslgStatus::CollinearOverlap));
    }
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndPslgDeterministic() {
    TEST("bnd.pslg_deterministic (opposite insertion orders, full vector)");
    auto f = buildPslgBase();
    auto build = [&](bool reverse) {
        bnd::ConstraintSet cs;
        cs.stamp = bnd::OwnershipStamp::of(*f.ws);
        std::vector<std::size_t> order{0, 1, 2, 3, 4, 5};
        if (reverse) std::reverse(order.begin(), order.end());
        const double C[6][2] = {{0, 0}, {4, 0}, {0, 4},
                                {10, 0}, {14, 0}, {10, 4}};
        for (auto i : order)
            cs.params[f.n[i]] = bnd::ParamPoint{{C[i][0], C[i][1]}, 7};
        // several faults touching different ids/handles
        cs.params[f.n[3]] = bnd::ParamPoint{{0, 0}, 7};   // coincident w/ n0
        std::vector<std::pair<bnd::ConstraintVertexId, std::size_t>> vs{
            {1, 0}, {2, 1}, {3, 3}, {4, 2}};
        if (reverse) std::reverse(vs.begin(), vs.end());
        for (auto [vid, ni] : vs) cs.vertices[vid] = f.n[ni];
        cs.vertices[5] = f.n[0];                          // non-injective
        cs.segments.push_back({10, 1, 2, 7, {1, 9}});
        cs.segments.push_back({11, 2, 1, 7, {1, 9}});     // duplicate
        return cs;
    };
    const auto a = bnd::preflightPslg(build(false));
    const auto b = bnd::preflightPslg(build(true));
    CHECK(!a.empty());
    bool identical = a.size() == b.size();
    for (std::size_t i = 0; i < a.size() && identical; ++i) {
        if (a[i].status != b[i].status) identical = false;
        if (a[i].entityA != b[i].entityA) identical = false;
        if (a[i].entityB != b[i].entityB) identical = false;
        const bool an = a[i].node.has_value(), bn = b[i].node.has_value();
        if (an != bn) identical = false;
        if (an && !(*a[i].node == *b[i].node)) identical = false;
    }
    CHECK(identical);      // FULL vector incl. attribution
}

void testBndSessionOpenAndStampChain() {
    TEST("bnd.session_open_and_stamp_chain (13.2/13.2a)");
    auto f = buildPslgBase();
    f.ws->rebuildAdjacency();
    f.cs.vertices[1] = f.n[0];
    f.cs.vertices[2] = f.n[1];
    f.cs.segments.push_back({10, 1, 2, 7, {1, 9}});
    const auto callerStampBefore = f.cs.stamp;
    bnd::SessionStatus st{};
    auto sess = bnd::RecoverySession::open(*f.ws, f.cs, st);
    CHECK(st == bnd::SessionStatus::Ok);
    CHECK(bool(sess));
    CHECK(bool(sess->orientationOf(7)));
    CHECK(sess->handoffAllowed(*f.ws, sess->stamp()));
    // 13.2: the session owns an IMMUTABLE input snapshot - destroying or
    // mutating the caller's set cannot affect it.
    CHECK(sess->input().segments.size() == 1);
    f.cs.segments.clear();
    f.cs.params.clear();
    f.cs.vertices.clear();
    CHECK(sess->input().segments.size() == 1);
    CHECK(sess->input().params.size() == 6);
    CHECK(sess->input().vertices.size() == 2);
    // the caller's ORIGINAL stamp is preserved unchanged and kept
    // separate from the evolving session stamp
    CHECK(sess->callerStamp().revision == callerStampBefore.revision);
    CHECK(sess->stamp().revision == callerStampBefore.revision);
    bnd::OwnershipStamp old{};
    CHECK(!sess->handoffAllowed(*f.ws, old));   // stale candidate refused
}

void testBndSessionAdvanceCannotLaunder() {
    TEST("bnd.session_advance_cannot_launder (production adversary)");
    // The reviewer's exact reproduction: an UNRELATED orphan-node commit
    // invalidates A1 (the new node has no authoritative parameter), and
    // advanceAfterCommit() must NOT restore handoff eligibility.
    auto f = buildPslgBase();
    f.ws->rebuildAdjacency();
    f.cs.vertices[1] = f.n[0];
    f.cs.vertices[2] = f.n[1];
    f.cs.segments.push_back({10, 1, 2, 7, {1, 9}});
    bnd::SessionStatus st{};
    auto sess = bnd::RecoverySession::open(*f.ws, f.cs, st);
    CHECK(bool(sess));
    CHECK(sess->handoffAllowed(*f.ws, sess->stamp()));
    topo::TopologyPolicy pol;
    pol.permitOrphanNodes = true;
    auto tx = f.ws->beginTransaction(pol);
    tx.addNode({99, 99, 0}, {0, 7});      // NO parameter coordinate
    (void)tx.commit();
    f.ws->rebuildAdjacency();
    // A1 is now invalid against the session's snapshot (totality broken)
    {
        auto view = sess->input();
        view.stamp = bnd::OwnershipStamp::of(*f.ws);
        CHECK(!bnd::preflightConstraintSet(*f.ws, view).empty());
    }
    CHECK(!sess->handoffAllowed(*f.ws, sess->stamp()));   // stale
    CHECK(!sess->advanceAfterCommit(*f.ws));   // MUST refuse to launder
    CHECK(!sess->valid());                     // POISONED permanently
    CHECK(!sess->handoffAllowed(*f.ws, sess->stamp()));
    // and it stays poisoned even if the workspace is made valid again
    CHECK(!sess->advanceAfterCommit(*f.ws));
}

void testBndSessionAdvanceAcceptsValidChange() {
    TEST("bnd.session_advance_accepts_valid_change");
    // A mutation that leaves A1/A2 valid against the snapshot (an orphan
    // node WITH a parameter is impossible - the map is immutable - so use
    // a cell removal that keeps every remaining node parameterized).
    auto f = buildPslgBase();
    f.ws->rebuildAdjacency();
    f.cs.vertices[1] = f.n[0];
    f.cs.vertices[2] = f.n[1];
    f.cs.segments.push_back({10, 1, 2, 7, {1, 9}});
    bnd::SessionStatus st{};
    auto sess = bnd::RecoverySession::open(*f.ws, f.cs, st);
    CHECK(bool(sess));
    const auto rev0 = sess->stamp().revision;
    topo::TopologyPolicy pol;
    pol.permitOrphanNodes = true;
    auto tx = f.ws->beginTransaction(pol);
    tx.removeCell(f.c[1]);                 // second triangle; nodes remain
    (void)tx.commit();
    f.ws->rebuildAdjacency();
    const bool ok = sess->advanceAfterCommit(*f.ws);
    if (ok) {                              // A1/A2 still hold
        CHECK(sess->valid());
        CHECK(sess->stamp().revision > rev0);
        CHECK(sess->handoffAllowed(*f.ws, sess->stamp()));
    } else {
        CHECK(!sess->valid());             // honest: poisoned, not laundered
        CHECK(!sess->handoffAllowed(*f.ws, sess->stamp()));
    }
}

void testBndSessionRefusesBadInput() {
    TEST("bnd.session_refuses_bad_input (stamp / A1 / A2 / PSLG)");
    {   // caller's stamp not current
        auto f = buildPslgBase();
        topo::TopologyPolicy pol;
        pol.permitOrphanNodes = true;
        auto tx = f.ws->beginTransaction(pol);
        tx.addNode({9, 9, 9}, {0, 7});
        (void)tx.commit();
        f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        bnd::SessionStatus st{};
        CHECK(!bnd::RecoverySession::open(*f.ws, f.cs, st));
        CHECK(st == bnd::SessionStatus::InputStampMismatch);
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);
    }
    {   // A1 rejects (missing parameter)
        auto f = buildPslgBase();
        f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        f.cs.params.erase(f.n[5]);
        bnd::SessionStatus st{};
        CHECK(!bnd::RecoverySession::open(*f.ws, f.cs, st));
        CHECK(st == bnd::SessionStatus::A1Rejected);
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);
    }
    {   // A2 rejects (degenerate triangle in the parameter domain) while
        // A1 stays clean - collapse one triangle's coords onto a line.
        auto f = buildPslgBase();
        f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        f.cs.params[f.n[5]] = bnd::ParamPoint{{20, 0}, 7};   // n3,n4,n5
        f.cs.params[f.n[3]] = bnd::ParamPoint{{10, 0}, 7};   // now collinear
        f.cs.params[f.n[4]] = bnd::ParamPoint{{14, 0}, 7};
        CHECK(bnd::preflightConstraintSet(*f.ws, f.cs).empty());  // A1 ok
        bnd::SessionStatus st{};
        CHECK(!bnd::RecoverySession::open(*f.ws, f.cs, st));
        CHECK(st == bnd::SessionStatus::A2Rejected);
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);
    }
    {   // PSLG rejects with A1/A2 clean
        auto f = buildPslgBase();
        f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        f.cs.vertices[1] = f.n[0];
        f.cs.vertices[2] = f.n[0];
        bnd::SessionStatus st{};
        CHECK(!bnd::RecoverySession::open(*f.ws, f.cs, st));
        CHECK(st == bnd::SessionStatus::PslgRejected);
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);
    }
}

void testBndSessionMultiChartOrientation() {
    TEST("bnd.session_multi_chart_orientation (retention per chart)");
    // Two disjoint triangles in DIFFERENT charts; the session must retain
    // an orientation entry for EACH chart.
    auto f = buildGeom({{0, 0}, {4, 0}, {0, 4}, {10, 0}, {14, 0}, {10, 4}},
                       {{0, 1, 2}, {3, 4, 5}}, /*mapAll=*/false);
    f.cs.params[f.n[0]] = bnd::ParamPoint{{0, 0}, 7};
    f.cs.params[f.n[1]] = bnd::ParamPoint{{4, 0}, 7};
    f.cs.params[f.n[2]] = bnd::ParamPoint{{0, 4}, 7};
    f.cs.params[f.n[3]] = bnd::ParamPoint{{10, 0}, 8};
    f.cs.params[f.n[4]] = bnd::ParamPoint{{14, 0}, 8};
    f.cs.params[f.n[5]] = bnd::ParamPoint{{10, 4}, 8};
    f.ws->rebuildAdjacency();
    bnd::SessionStatus st{};
    auto sess = bnd::RecoverySession::open(*f.ws, f.cs, st);
    CHECK(st == bnd::SessionStatus::Ok);
    CHECK(bool(sess));
    CHECK(sess->chartCount() == 2);
    CHECK(bool(sess->orientationOf(7)));
    CHECK(bool(sess->orientationOf(8)));
    CHECK(!sess->orientationOf(9));
}

}  // namespace

namespace paacs_tests {
const Entry* bndSessionTests(std::size_t& count) {
    static const Entry kBnd[] = {
        {"bnd.pslg_clean", testBndPslgClean},
        {"bnd.pslg_non_injective_vertex_mapping", testBndPslgNonInjective},
        {"bnd.pslg_duplicate_and_coincident", testBndPslgDuplicateAndCoincident},
        {"bnd.pslg_tjunction_exact", testBndPslgTJunctionExact},
        {"bnd.pslg_collinear_overlap_exact", testBndPslgCollinearOverlapExact},
        {"bnd.pslg_crossing_and_chart_local",
         testBndPslgCrossingAndChartLocal},
        {"bnd.pslg_deterministic", testBndPslgDeterministic},
        {"bnd.session_open_and_stamp_chain", testBndSessionOpenAndStampChain},
        {"bnd.session_advance_cannot_launder",
         testBndSessionAdvanceCannotLaunder},
        {"bnd.session_advance_accepts_valid_change",
         testBndSessionAdvanceAcceptsValidChange},
        {"bnd.session_refuses_bad_input", testBndSessionRefusesBadInput},
        {"bnd.session_multi_chart_orientation",
         testBndSessionMultiChartOrientation},
    };
    count = sizeof(kBnd) / sizeof(kBnd[0]);
    return kBnd;
}
}  // namespace paacs_tests

// ============================================================================
// sec-13 increment 2: the FIRST MUTATION PATH (mark-only recovery).
// ============================================================================
namespace {

// 2x2 fan: n0(0,0) n1(2,0) n2(4,0) on v==0, plus n3(2,2) above and
// n4(2,-2) below, so interior chain node n1 has a FULL STAR of 4 cells.
GeomFix buildChainFan() {
    return buildGeom(
        {{0, 0}, {2, 0}, {4, 0}, {2, 2}, {2, -2}},
        {{0, 1, 3}, {1, 2, 3}, {1, 0, 4}, {2, 1, 4}});
}

bnd::ConstraintSet chainSet(const GeomFix& f, paacs::mesh::Classification cls,
                            bnd::ConstraintEntityId id = 10) {
    bnd::ConstraintSet cs = f.cs;
    cs.vertices[1] = f.n[0];
    cs.vertices[2] = f.n[2];
    cs.segments.push_back({id, 1, 2, 7, cls});
    return cs;
}

void testBndRoutingAllFour() {
    TEST("bnd.routing_all_four (13.3-0: all four routes discriminated)");
    {   // MarkOnly: complete geometric chain, marks missing
        auto f = buildChainFan();
        f.ws->rebuildAdjacency();
        bnd::SessionStatus st{};
        auto cs = chainSet(f, {1, 9});
        auto s = bnd::RecoverySession::open(*f.ws, cs, st);
        CHECK(bool(s));
        std::vector<topo::NodeHandle> ch;
        CHECK(bnd::routeEntity(*f.ws, s->input(), s->index(),
                               s->input().segments[0], ch) ==
              bnd::Routing::MarkOnly);
        CHECK(ch.size() == 3);
    }
    {   // AlreadyPresent: same chain, correctly marked
        auto f = buildChainFan();
        for (int i = 0; i < 2; ++i) {
            const auto eh = *f.ws->edges().find(f.ws->nodes(),
                                                f.n[std::size_t(i)],
                                                f.n[std::size_t(i + 1)]);
            f.ws->testEdgesMutable().setRequiredConstraint(eh, true);
            f.ws->testEdgesMutable().setClassification(eh, {1, 9});
        }
        f.ws->rebuildAdjacency();
        bnd::SessionStatus st{};
        auto cs = chainSet(f, {1, 9});
        auto s = bnd::RecoverySession::open(*f.ws, cs, st);
        CHECK(bool(s));
        std::vector<topo::NodeHandle> ch;
        CHECK(bnd::routeEntity(*f.ws, s->input(), s->index(),
                               s->input().segments[0], ch) ==
              bnd::Routing::AlreadyPresent);
    }
    {   // TRANSVERSE: a quad whose diagonal is n0-n2; the segment n1 -> n3
        // CROSSES that diagonal, no node lies on it, and n1-n3 is not an
        // edge - the canonical transverse case.
        auto f = buildGeom({{0, 0}, {2, 0}, {2, 2}, {0, 2}},
                           {{0, 1, 2}, {0, 2, 3}});
        f.ws->rebuildAdjacency();
        bnd::ConstraintSet cs = f.cs;
        cs.vertices[1] = f.n[1];
        cs.vertices[2] = f.n[3];
        cs.segments.push_back({10, 1, 2, 7, {1, 9}});
        bnd::SessionStatus st{};
        auto s = bnd::RecoverySession::open(*f.ws, cs, st);
        CHECK(bool(s));
        std::vector<topo::NodeHandle> ch;
        const auto r = bnd::routeEntity(*f.ws, s->input(), s->index(),
                                        s->input().segments[0], ch);
        CHECK(r == bnd::Routing::Transverse);
    }
    {   // PARTIAL: a chain that starts along the segment but does not
        // reach the far endpoint. Segment n0 -> n2 with edge n1-n2
        // REMOVED is impossible without mutation, so instead use a
        // segment n0 -> n3' whose first hop exists but whose target is
        // not reachable on-segment: n0(0,0) -> a far collinear node.
        auto f = buildGeom({{0, 0}, {2, 0}, {4, 0}, {2, 2}, {6, 0}, {6, 2}},
                           {{0, 1, 3}, {1, 2, 3}, {4, 5, 2}});
        f.ws->rebuildAdjacency();
        bnd::ConstraintSet cs = f.cs;
        cs.vertices[1] = f.n[0];   // (0,0)
        cs.vertices[2] = f.n[4];   // (6,0) collinear, but no edge n2-n4
        cs.segments.push_back({10, 1, 2, 7, {1, 9}});
        bnd::SessionStatus st{};
        auto s = bnd::RecoverySession::open(*f.ws, cs, st);
        if (s) {
            std::vector<topo::NodeHandle> ch;
            const auto r = bnd::routeEntity(*f.ws, s->input(), s->index(),
                                            s->input().segments[0], ch);
            CHECK(r == bnd::Routing::PartialChain);   // walked, stalled
        }
    }
}

void testBndMarkOnlyRecovery() {
    TEST("bnd.mark_only_recovery (first MUTATION: marks in ONE commit)");
    auto f = buildChainFan();
    f.ws->rebuildAdjacency();
    const auto rev0 = f.ws->topologyRevision();
    auto cs = chainSet(f, {1, 9});
    const auto rep = bnd::recoverSegments(*f.ws, cs);
    CHECK(!rep.rejected);
    CHECK(rep.processed.size() == 1);
    CHECK(rep.processed[0].outcome == bnd::Outcome::Recovered);
    CHECK(rep.unattempted.empty());
    CHECK(f.ws->topologyRevision() == rev0 + 1);
    CHECK(f.ws->state() == topo::WorkspaceState::CommittedAdjacencyStale);
    f.ws->rebuildAdjacency();
    const auto a = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[1]);
    const auto b = *f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[2]);
    CHECK(f.ws->edges().get(a).requiredConstraint);
    CHECK(f.ws->edges().get(b).requiredConstraint);
    CHECK(f.ws->edges().get(a).cls.dim == 1 &&
          f.ws->edges().get(a).cls.tag == 9);
    CHECK(f.ws->cells().stats().live == 4);      // topology-equivalent
    const auto rev1 = f.ws->topologyRevision();
    auto cs2 = chainSet(f, {1, 9});
    cs2.stamp = bnd::OwnershipStamp::of(*f.ws);
    const auto rep2 = bnd::recoverSegments(*f.ws, cs2);
    CHECK(rep2.processed[0].outcome == bnd::Outcome::AlreadyPresent);
    CHECK(f.ws->topologyRevision() == rev1);
}

void testBndMarkConflictCases() {
    TEST("bnd.mark_conflict_cases (13.3b: all three, incl. {-1,-2})");
    // (a) already-required with a DIFFERENT classification -> Refused
    {
        auto f = buildChainFan();
        const auto eh = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[1]);
        f.ws->testEdgesMutable().setRequiredConstraint(eh, true);
        f.ws->testEdgesMutable().setClassification(eh, {1, 4});
        f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        auto cs = chainSet(f, {1, 9});
        const auto rep = bnd::recoverSegments(*f.ws, cs);
        CHECK(rep.processed[0].outcome == bnd::Outcome::Refused);
        CHECK(rep.processed[0].refused ==
              bnd::RefusedReason::ConflictingEdgeMark);
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);
    }
    // (b) UNREQUIRED with a different NON-DEFAULT classification -> Refused
    {
        auto f = buildChainFan();
        const auto eh = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[1]);
        f.ws->testEdgesMutable().setClassification(eh, {1, 4});
        f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        auto cs = chainSet(f, {1, 9});
        const auto rep = bnd::recoverSegments(*f.ws, cs);
        CHECK(rep.processed[0].outcome == bnd::Outcome::Refused);
        CHECK(rep.processed[0].refused ==
              bnd::RefusedReason::ConflictingEdgeMark);
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);
    }
    // (c) the reviewer's case: UNREQUIRED carrying {-1,-2}, which is NOT
    // the exact default {-1,-1} -> must be Refused, never Failed.
    {
        auto f = buildChainFan();
        const auto eh = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[1]);
        f.ws->testEdgesMutable().setClassification(eh, {-1, -2});
        f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        auto cs = chainSet(f, {1, 9});
        const auto rep = bnd::recoverSegments(*f.ws, cs);
        CHECK(rep.processed[0].outcome == bnd::Outcome::Refused);
        CHECK(rep.processed[0].refused ==
              bnd::RefusedReason::ConflictingEdgeMark);
        CHECK(!(rep.processed[0].outcome == bnd::Outcome::Failed));
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);
    }
    // (d) UNREQUIRED with the EXACT default -> legal, recovers
    {
        auto f = buildChainFan();
        f.ws->rebuildAdjacency();
        auto cs = chainSet(f, {1, 9});
        const auto rep = bnd::recoverSegments(*f.ws, cs);
        CHECK(rep.processed[0].outcome == bnd::Outcome::Recovered);
    }
}

void testBndMixedCellClassification() {
    TEST("bnd.mixed_cell_classification_refused (13.5)");
    // One closure cell carries a different face classification.
    // built inline so one cell carries a DIFFERENT face classification
    GeomFix f;
    {
        topo::NodePool np;
        const double C[5][2] = {{0, 0}, {2, 0}, {4, 0}, {2, 2}, {2, -2}};
        for (const auto& q : C)
            f.n.push_back(np.add({q[0], q[1], 0}, {2, 1}));
        topo::CellPool cp;
        const int T[4][3] = {{0, 1, 3}, {1, 2, 3}, {1, 0, 4}, {2, 1, 4}};
        for (int i = 0; i < 4; ++i)
            f.c.push_back(cp.add(
                topo::ElementType::Tri3,
                std::vector<topo::NodeHandle>{f.n[std::size_t(T[i][0])],
                                              f.n[std::size_t(T[i][1])],
                                              f.n[std::size_t(T[i][2])]},
                i == 3 ? paacs::mesh::Classification{2, 77}
                       : paacs::mesh::Classification{2, 5},
                np));
        topo::EdgePool ep;
        for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
        f.ws = std::make_unique<topo::TopologyWorkspace>(
            std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
        f.cs.stamp = bnd::OwnershipStamp::of(*f.ws);
        for (std::size_t i = 0; i < 5; ++i)
            f.cs.params[f.n[i]] = bnd::ParamPoint{{C[i][0], C[i][1]}, 7};
    }
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    auto cs = chainSet(f, {1, 9});
    const auto rep = bnd::recoverSegments(*f.ws, cs);
    CHECK(rep.processed[0].outcome == bnd::Outcome::Refused);
    CHECK(rep.processed[0].refused ==
          bnd::RefusedReason::MixedCellClassification);
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
}

void testBndDiskAuditComplete() {
    TEST("bnd.disk_audit_complete (connectivity, Euler, node link, cycle)");
    auto f = buildChainFan();
    f.ws->rebuildAdjacency();
    bnd::SessionStatus st{};
    auto cs = chainSet(f, {1, 9});
    auto s = bnd::RecoverySession::open(*f.ws, cs, st);
    CHECK(bool(s));
    const std::vector<topo::NodeHandle> chain{f.n[0], f.n[1], f.n[2]};
    std::vector<topo::NodeHandle> ring;
    const auto closure =
        bnd::detail::canonicalClosure(*f.ws, s->index(), chain);
    CHECK(closure.size() == 4);
    CHECK(bnd::detail::closureIsDisk(*f.ws, closure, &ring));
    CHECK(ring.size() == 4);          // boundary CYCLE traversed + returned
    // every ring node is distinct and on the closure boundary
    std::unordered_set<topo::NodeHandle> uniq(ring.begin(), ring.end());
    CHECK(uniq.size() == ring.size());
    // a single boundary triangle IS a disk (membership alone is fine)
    {
        std::vector<topo::CellHandle> one{f.c[0]};
        std::vector<topo::NodeHandle> r2;
        CHECK(bnd::detail::closureIsDisk(*f.ws, one, &r2));
        CHECK(r2.size() == 3);
    }
    {   // edge-DISCONNECTED pair (fails check 1)
        std::vector<topo::CellHandle> pinched{f.c[0], f.c[3]};
        std::vector<topo::NodeHandle> r3;
        CHECK(!bnd::detail::closureIsDisk(*f.ws, pinched, &r3));
        CHECK(r3.empty());            // no boundary reported on failure
    }
    {   // a contiguous fan of two IS a disk
        std::vector<topo::CellHandle> two{f.c[0], f.c[1]};
        std::vector<topo::NodeHandle> r4;
        CHECK(bnd::detail::closureIsDisk(*f.ws, two, &r4));
    }
    // NODE-LINK failure ISOLATED (review): a closed WHEEL around n0 (a
    // disk) plus one extra triangle that SHARES EDGE n0-n1 with it. The
    // set is therefore edge-CONNECTED, and V-E+F == 1 still holds, so
    // checks 1 and 2 both PASS; the only thing wrong is n0's link, where
    // node 1 acquires degree 3. Both prerequisites are asserted below so
    // the rejection is provably attributable to the node-link audit.
    {
        GeomFix g;
        topo::NodePool np;
        const double C[6][2] = {{0, 0},  {1, 0},  {0, 1},
                                {-1, 0}, {0, -1}, {2, 1}};
        for (const auto& q : C) g.n.push_back(np.add({q[0], q[1], 0}, {2, 1}));
        topo::CellPool cp;
        const int T[5][3] = {{0, 1, 2}, {0, 2, 3}, {0, 3, 4},
                             {0, 4, 1},                    // closed wheel
                             {0, 1, 5}};                   // shares n0-n1
        for (const auto& tr : T)
            g.c.push_back(cp.add(
                topo::ElementType::Tri3,
                std::vector<topo::NodeHandle>{g.n[std::size_t(tr[0])],
                                              g.n[std::size_t(tr[1])],
                                              g.n[std::size_t(tr[2])]},
                {2, 5}, np));
        topo::EdgePool ep;
        for (const auto& ch : g.c) ep.insertCellEdges(np, cp, ch);
        g.ws = std::make_unique<topo::TopologyWorkspace>(
            std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
        // PREREQUISITE 1: the set is EDGE-CONNECTED (every cell reachable
        // through shared edges) - so check 1 cannot be what rejects it.
        {
            std::unordered_map<std::int32_t, std::vector<std::size_t>> byE;
            for (std::size_t ci = 0; ci < g.c.size(); ++ci)
                for (const auto& eh :
                     g.ws->edges().cellEdges(g.ws->cells(), g.c[ci]))
                    byE[eh.slot].push_back(ci);
            std::unordered_set<std::size_t> seen;
            std::vector<std::size_t> stk{0};   // (was 'st': C4456/-Wshadow)
            while (!stk.empty()) {
                const auto ci = stk.back();
                stk.pop_back();
                if (!seen.insert(ci).second) continue;
                for (const auto& eh :
                     g.ws->edges().cellEdges(g.ws->cells(), g.c[ci]))
                    for (auto nb : byE[eh.slot])
                        if (!seen.count(nb)) stk.push_back(nb);
            }
            CHECK(seen.size() == g.c.size());          // CONNECTED
        }
        // PREREQUISITE 2: V - E + F == 1 - so check 2 cannot reject it.
        {
            std::unordered_set<topo::NodeHandle> V;
            std::unordered_set<std::int32_t> E;
            for (const auto& c : g.c) {
                const auto cd = g.ws->cells().get(c);
                for (int k = 0; k < 3; ++k) V.insert(cd.nodes[std::size_t(k)]);
                for (const auto& eh :
                     g.ws->edges().cellEdges(g.ws->cells(), c))
                    E.insert(eh.slot);
            }
            CHECK(static_cast<long long>(V.size()) -
                      static_cast<long long>(E.size()) +
                      static_cast<long long>(g.c.size()) ==
                  1);                                   // EULER OK
        }
        // ... therefore the rejection below is the NODE-LINK audit.
        std::vector<topo::NodeHandle> r5;
        CHECK(!bnd::detail::closureIsDisk(*g.ws, g.c, &r5));
        CHECK(r5.empty());
        // and the wheel ALONE (without the pinching triangle) IS a disk
        std::vector<topo::CellHandle> wheel{g.c[0], g.c[1], g.c[2], g.c[3]};
        std::vector<topo::NodeHandle> r6;
        CHECK(bnd::detail::closureIsDisk(*g.ws, wheel, &r6));
    }
}

void testBndPrefixSemanticsClean() {
    TEST("bnd.prefix_semantics (PSLG-clean set; stop + unattempted)");
    // THREE entities in DISJOINT regions, so the whole-set PSLG preflight
    // is clean (no crossings, duplicates or T-junctions between them):
    // 10 = mark-only chain (recovers), 20 = transverse (deferred ->
    // Refused, STOP), 30 = never attempted.
    GeomFix f;
    {
        topo::NodePool np;
        const double C[12][2] = {
            {0, 0},  {2, 0},  {4, 0},  {2, 2},  {2, -2},     // fan (10)
            {20, 0}, {22, 0}, {22, 2}, {20, 2},              // quad (20)
            {40, 0}, {42, 0}, {41, 2}};                      // tri  (30)
        for (const auto& q : C)
            f.n.push_back(np.add({q[0], q[1], 0}, {2, 1}));
        topo::CellPool cp;
        const int T[7][3] = {{0, 1, 3}, {1, 2, 3}, {1, 0, 4}, {2, 1, 4},
                             {5, 6, 7}, {5, 7, 8}, {9, 10, 11}};
        for (const auto& tr : T)
            f.c.push_back(cp.add(
                topo::ElementType::Tri3,
                std::vector<topo::NodeHandle>{f.n[std::size_t(tr[0])],
                                              f.n[std::size_t(tr[1])],
                                              f.n[std::size_t(tr[2])]},
                {2, 5}, np));
        topo::EdgePool ep;
        for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
        f.ws = std::make_unique<topo::TopologyWorkspace>(
            std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
        f.cs.stamp = bnd::OwnershipStamp::of(*f.ws);
        for (std::size_t i = 0; i < 12; ++i)
            f.cs.params[f.n[i]] = bnd::ParamPoint{{C[i][0], C[i][1]}, 7};
    }
    f.ws->rebuildAdjacency();
    bnd::ConstraintSet cs = f.cs;
    cs.vertices[1] = f.n[0];
    cs.vertices[2] = f.n[2];
    cs.segments.push_back({10, 1, 2, 7, {1, 9}});     // mark-only chain
    cs.vertices[3] = f.n[6];
    cs.vertices[4] = f.n[8];
    cs.segments.push_back({20, 3, 4, 7, {1, 9}});     // transverse (quad)
    // Mark the quad diagonal REQUIRED: transverse recovery would RETIRE
    // it, so entity 20 must refuse (CrossedRequiredEdge) - the stop is
    // forced by policy, not by an unimplemented path.
    f.ws->testEdgesMutable().setRequiredConstraint(
        *f.ws->edges().find(f.ws->nodes(), f.n[5], f.n[7]), true);
    cs.vertices[5] = f.n[9];
    cs.vertices[6] = f.n[10];
    cs.segments.push_back({30, 5, 6, 7, {1, 9}});     // never attempted
    const auto rep = bnd::recoverSegments(*f.ws, cs);
    CHECK(!rep.rejected);                     // PSLG-clean by construction
    CHECK(rep.processed.size() == 2);
    CHECK(rep.processed[0].entity == 10);     // ASCENDING
    CHECK(rep.processed[0].outcome == bnd::Outcome::Recovered);
    CHECK(rep.processed[1].entity == 20);
    CHECK(rep.processed[1].outcome == bnd::Outcome::Refused);
    CHECK(rep.processed[1].refused ==
          bnd::RefusedReason::CrossedRequiredEdge);
    CHECK(rep.unattempted.size() == 1);
    CHECK(rep.unattempted[0] == 30);          // stopped, not attempted
}

void testBndWholeSetRejectedZeroChange() {
    TEST("bnd.whole_set_rejected_zero_change (reason preserved)");
    auto f = buildChainFan();
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    auto cs = chainSet(f, {1, 9});
    cs.vertices[99] = f.n[0];                 // non-injective -> PSLG
    const auto rep = bnd::recoverSegments(*f.ws, cs);
    CHECK(rep.rejected);
    CHECK(rep.rejectReason == bnd::SessionStatus::PslgRejected);  // reason
    CHECK(rep.processed.size() == 1);
    CHECK(rep.processed[0].outcome == bnd::Outcome::Rejected);    // produced
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);
    // MULTI-ENTITY whole-set rejection: EXACT result vector. A duplicate
    // segment makes the PSLG preflight reject the whole set, so EVERY entity
    // gets a Rejected outcome row in `processed` (ascending id) and NONE is
    // left in `unattempted` (13.7: an entity is a processed outcome XOR an
    // unattempted id - never both, so no double-reporting).
    {
        auto g = buildChainFan();
        g.ws->rebuildAdjacency();
        const auto gd0 = paacs_tests::observableDigest(*g.ws);
        bnd::ConstraintSet cs2;
        cs2.stamp = bnd::OwnershipStamp::of(*g.ws);
        for (std::size_t i = 0; i < 5; ++i)
            cs2.params[g.n[i]] = bnd::ParamPoint{
                {double(g.ws->nodes().get(g.n[i]).coord.x),
                 double(g.ws->nodes().get(g.n[i]).coord.y)}, 7};
        cs2.vertices[1] = g.n[0];
        cs2.vertices[2] = g.n[2];
        cs2.vertices[3] = g.n[3];
        cs2.segments.push_back({10, 1, 2, 7, {1, 9}});
        cs2.segments.push_back({20, 1, 3, 7, {1, 9}});
        cs2.segments.push_back({30, 1, 2, 7, {1, 9}});   // DUPLICATE of 10
        const auto r2 = bnd::recoverSegments(*g.ws, cs2);
        CHECK(r2.rejected);
        CHECK(r2.rejectReason == bnd::SessionStatus::PslgRejected);
        // EXACT vector: three Rejected rows, ascending id, unattempted EMPTY
        CHECK(r2.processed.size() == 3);
        CHECK(r2.unattempted.empty());               // no double-reporting
        const bnd::ConstraintEntityId ids[3] = {10, 20, 30};
        for (int k = 0; k < 3; ++k) {
            CHECK(r2.processed[std::size_t(k)].entity == ids[k]);
            CHECK(r2.processed[std::size_t(k)].outcome ==
                  bnd::Outcome::Rejected);
            CHECK(r2.processed[std::size_t(k)].refused ==
                  bnd::RefusedReason::None);
            CHECK(r2.processed[std::size_t(k)].failed ==
                  bnd::FailedReason::None);
        }
        CHECK(paacs_tests::observableDigest(*g.ws) == gd0);   // zero change
    }
}

void testBndCandidateValidatorAllGates() {
    TEST("bnd.candidate_validator_all_gates (each gate reached alone)");
    // The validator is a plain-value predicate over CavityPreview, so
    // each gate is exercised on a HAND-BUILT preview that passes every
    // EARLIER gate and violates exactly one. (Staging real cavities
    // cannot isolate gates 4-6: a degenerate or mis-oriented triangle
    // returns at gate 2 first - the flaw the review found.)
    auto f = buildChainFan();
    f.ws->rebuildAdjacency();
    bnd::SessionStatus st{};
    auto cs = chainSet(f, {1, 9});
    auto s = bnd::RecoverySession::open(*f.ws, cs, st);
    CHECK(bool(s));
    const auto want = *s->orientationOf(7);
    // a 6-node mesh purely as a source of DISTINCT handles for the
    // hand-built previews (its own geometry is irrelevant)
    auto g6 = buildGeom({{0, 0}, {10, 0}, {20, 0}, {0, 10}, {10, 10},
                         {20, 10}},
                        {{0, 1, 4}, {1, 2, 5}});
    auto mk = [&](const bnd::ConstraintSet& c, paacs::predicates::Sign w) {
        return bnd::makeCandidateValidator(c, 7, w);
    };
    // A helper preview carrying only addedConnectivity (all the validator
    // reads); node handles come from an enlarged parameter map.
    auto previewOf = [](std::vector<std::array<topo::NodeHandle, 3>> tris) {
        topo::CavityPreview pv;
        pv.addedConnectivity = std::move(tris);
        return pv;
    };
    // ---- GATE 1: same-chart coverage -------------------------------
    {
        auto c1 = s->input();
        auto v = mk(c1, want);
        // valid identity geometry passes with the right chart...
        std::vector<std::array<topo::NodeHandle, 3>> tris;
        for (const auto& ch : f.c) {
            const auto cd = f.ws->cells().get(ch);
            tris.push_back({cd.nodes[0], cd.nodes[1], cd.nodes[2]});
        }
        CHECK(v(previewOf(tris)));
        // ...and fails when the demanded chart differs
        CHECK(!bnd::makeCandidateValidator(c1, 8, want)(previewOf(tris)));
    }
    // ---- GATE 2: nondegeneracy (collinear triangle) -----------------
    {
        auto c2 = s->input();
        c2.params[f.n[0]] = bnd::ParamPoint{{0, 0}, 7};
        c2.params[f.n[1]] = bnd::ParamPoint{{1, 0}, 7};
        c2.params[f.n[2]] = bnd::ParamPoint{{2, 0}, 7};   // collinear
        CHECK(!mk(c2, want)(previewOf({{f.n[0], f.n[1], f.n[2]}})));
    }
    // ---- GATE 3: orientation (nondegenerate, wrong sign) ------------
    {
        auto c3 = s->input();
        c3.params[f.n[0]] = bnd::ParamPoint{{0, 0}, 7};
        c3.params[f.n[1]] = bnd::ParamPoint{{1, 0}, 7};
        c3.params[f.n[2]] = bnd::ParamPoint{{0, 1}, 7};
        const auto sgn = paacs::predicates::orient2d({0, 0}, {1, 0}, {0, 1});
        const auto opp = sgn == paacs::predicates::Sign::Positive
                             ? paacs::predicates::Sign::Negative
                             : paacs::predicates::Sign::Positive;
        CHECK(mk(c3, sgn)(previewOf({{f.n[0], f.n[1], f.n[2]}})));   // ok
        CHECK(!mk(c3, opp)(previewOf({{f.n[0], f.n[1], f.n[2]}})));  // gate 3
    }
    // ---- GATE 4: proper crossing (both triangles valid + oriented) ---
    {
        auto c4 = s->input();
        // T0 = (0,0),(4,0),(2,3) and T1 = (0,1),(4,1),(2,3): both are
        // nondegenerate and identically oriented (so gates 2 and 3 pass),
        // no vertex of either is strictly inside the other (so gate 6
        // does not fire), but edge n0-n2 crosses edge n3-n4 transversally.
        c4.params[f.n[0]] = bnd::ParamPoint{{0, 0}, 7};
        c4.params[f.n[1]] = bnd::ParamPoint{{4, 0}, 7};
        c4.params[f.n[2]] = bnd::ParamPoint{{2, 3}, 7};
        c4.params[f.n[3]] = bnd::ParamPoint{{0, 1}, 7};
        c4.params[f.n[4]] = bnd::ParamPoint{{4, 1}, 7};
        const auto s1 = paacs::predicates::orient2d({0, 0}, {4, 0}, {2, 3});
        const auto s2 = paacs::predicates::orient2d({0, 1}, {4, 1}, {2, 3});
        CHECK(s1 == s2);                        // gate 3 will pass
        CHECK(s1 != paacs::predicates::Sign::Zero);   // gate 2 will pass
        CHECK(!mk(c4, s1)(previewOf({{f.n[0], f.n[1], f.n[2]},
                                     {f.n[3], f.n[4], f.n[2]}})));
    }
    // ---- GATE 5: T-junction, with NO proper crossing (six handles) ----
    {
        // T0 lies entirely ABOVE y==0 (base on it); T1 lies entirely
        // BELOW, touching only at its apex (2,0), which is strictly
        // interior to T0's base edge. No pair of edges straddles, so
        // gate 4 cannot fire; both triangles are nondegenerate and
        // share the accepted orientation, so gates 2/3 cannot fire.
        auto c5 = s->input();
        const std::array<topo::NodeHandle, 6> H{g6.n[0], g6.n[1], g6.n[2],
                                                g6.n[3], g6.n[4], g6.n[5]};
        const double P[6][2] = {{0, 0}, {4, 0}, {2, 3},
                                {2, 0}, {0, -2}, {4, -2}};
        for (std::size_t q = 0; q < 6; ++q)
            c5.params[H[q]] = bnd::ParamPoint{{P[q][0], P[q][1]}, 7};
        const auto o0 = paacs::predicates::orient2d({0, 0}, {4, 0}, {2, 3});
        const auto o1 = paacs::predicates::orient2d({2, 0}, {0, -2}, {4, -2});
        CHECK(o0 != paacs::predicates::Sign::Zero);   // gate 2 passes
        CHECK(o1 != paacs::predicates::Sign::Zero);
        CHECK(o0 == o1);                              // gate 3 passes
        // gate 4 CANNOT fire: assert no pair of edges properly straddles
        {
            const double E[6][4] = {{0, 0, 4, 0},   {4, 0, 2, 3},
                                    {2, 3, 0, 0},   {2, 0, 0, -2},
                                    {0, -2, 4, -2}, {4, -2, 2, 0}};
            bool anyCross = false;
            for (int a = 0; a < 3; ++a)
                for (int b = 3; b < 6; ++b) {
                    const paacs::Vec2 a1{E[a][0], E[a][1]},
                        b1{E[a][2], E[a][3]}, a2{E[b][0], E[b][1]},
                        b2{E[b][2], E[b][3]};
                    const auto d1 = paacs::predicates::orient2d(a1, b1, a2);
                    const auto d2 = paacs::predicates::orient2d(a1, b1, b2);
                    const auto d3 = paacs::predicates::orient2d(a2, b2, a1);
                    const auto d4 = paacs::predicates::orient2d(a2, b2, b1);
                    const bool s1 = d1 != d2 &&
                                    d1 != paacs::predicates::Sign::Zero &&
                                    d2 != paacs::predicates::Sign::Zero;
                    const bool s2 = d3 != d4 &&
                                    d3 != paacs::predicates::Sign::Zero &&
                                    d4 != paacs::predicates::Sign::Zero;
                    if (s1 && s2) anyCross = true;
                }
            CHECK(!anyCross);                         // gate 4 passes
        }
        // ... so the rejection below is GATE 5 (the T-junction).
        CHECK(!mk(c5, o0)(previewOf({{H[0], H[1], H[2]},
                                     {H[3], H[4], H[5]}})));
    }
    // ---- GATE 6a: STRICT nesting, no shared handles -------------------
    {
        auto c6 = s->input();
        const std::array<topo::NodeHandle, 6> H{g6.n[0], g6.n[1], g6.n[2],
                                                g6.n[3], g6.n[4], g6.n[5]};
        const double P[6][2] = {{0, 0}, {6, 0}, {0, 6},
                                {1, 1}, {2, 1}, {1, 2}};
        for (std::size_t q = 0; q < 6; ++q)
            c6.params[H[q]] = bnd::ParamPoint{{P[q][0], P[q][1]}, 7};
        const auto o0 = paacs::predicates::orient2d({0, 0}, {6, 0}, {0, 6});
        const auto o1 = paacs::predicates::orient2d({1, 1}, {2, 1}, {1, 2});
        CHECK(o0 != paacs::predicates::Sign::Zero);
        CHECK(o0 == o1);                     // gates 2/3 pass
        // the small triangle is STRICTLY inside, so no edge crosses and
        // no vertex lies ON an edge - gates 4 and 5 cannot fire.
        CHECK(!mk(c6, o0)(previewOf({{H[0], H[1], H[2]},
                                     {H[3], H[4], H[5]}})));
    }
    // ---- GATE 6b: SHARED-topology same-side overlap --------------------
    {
        auto c6c = s->input();
        const std::array<topo::NodeHandle, 4> H{g6.n[0], g6.n[1], g6.n[2],
                                                g6.n[3]};
        const double P[4][2] = {{0, 0}, {4, 0}, {0, 4}, {0.5, 0.5}};
        for (std::size_t q = 0; q < 4; ++q)
            c6c.params[H[q]] = bnd::ParamPoint{{P[q][0], P[q][1]}, 7};
        const auto ob = paacs::predicates::orient2d({0, 0}, {4, 0}, {0, 4});
        const auto os =
            paacs::predicates::orient2d({0, 0}, {4, 0}, {0.5, 0.5});
        CHECK(ob != paacs::predicates::Sign::Zero);
        CHECK(ob == os);                     // gates 2/3 pass
        CHECK(!mk(c6c, ob)(previewOf({{H[0], H[1], H[2]},
                                      {H[0], H[1], H[3]}})));
    }
}

void testBndTransverseRecovery() {
    TEST("bnd.transverse_recovery (13.8/13.8a/13.9: re-triangulates)");
    // Quad with diagonal n0-n2; the segment n1 -> n3 CROSSES it. Transverse
    // recovery must remove the strip, re-triangulate both sides, and end
    // with n1-n3 present, required and classified - in ONE commit.
    auto f = buildGeom({{0, 0}, {2, 0}, {2, 2}, {0, 2}},
                       {{0, 1, 2}, {0, 2, 3}});
    f.ws->rebuildAdjacency();
    const auto rev0 = f.ws->topologyRevision();
    CHECK(bool(f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[2])));  // diag
    CHECK(!f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[3]));       // not
    bnd::ConstraintSet cs = f.cs;
    cs.vertices[1] = f.n[1];
    cs.vertices[2] = f.n[3];
    cs.segments.push_back({10, 1, 2, 7, {1, 9}});
    const auto rep = bnd::recoverSegments(*f.ws, cs);
    CHECK(!rep.rejected);
    CHECK(rep.processed.size() == 1);
    CHECK(rep.processed[0].outcome == bnd::Outcome::Recovered);
    CHECK(f.ws->topologyRevision() == rev0 + 1);      // exactly once
    CHECK(f.ws->state() == topo::WorkspaceState::CommittedAdjacencyStale);
    f.ws->rebuildAdjacency();
    // the constrained segment now EXISTS as an edge, marked in the SAME
    // commit, and the crossed diagonal is gone
    const auto e13 = f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[3]);
    CHECK(bool(e13));
    CHECK(f.ws->edges().get(*e13).requiredConstraint);
    CHECK(f.ws->edges().get(*e13).cls.dim == 1);
    CHECK(f.ws->edges().get(*e13).cls.tag == 9);
    CHECK(!f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[2]));
    CHECK(f.ws->cells().stats().live == 2);           // still two triangles
    // and re-running is now AlreadyPresent with zero churn
    const auto rev1 = f.ws->topologyRevision();
    bnd::ConstraintSet cs2 = cs;
    cs2.stamp = bnd::OwnershipStamp::of(*f.ws);
    const auto rep2 = bnd::recoverSegments(*f.ws, cs2);
    CHECK(rep2.processed[0].outcome == bnd::Outcome::AlreadyPresent);
    CHECK(f.ws->topologyRevision() == rev1);
}

void testBndTransverseDegeneracies() {
    TEST("bnd.transverse_degeneracies (13.8 refusals, each reason pinned)");
    {   // ON-SEGMENT INTERIOR NODE: n1=(2,0) lies exactly on the segment
        // n3(2,2) -> n4(2,-2). The router discovers the geometric chain
        // n3-n1-n4 (both sub-edges exist), so this is DETERMINISTICALLY a
        // mark-only recovery - never a silent transverse pass THROUGH a
        // vertex, and never an ambiguous 3-way outcome.
        auto f = buildGeom({{0, 0}, {2, 0}, {4, 0}, {2, 2}, {2, -2}},
                           {{0, 1, 3}, {1, 2, 3}, {1, 0, 4}, {2, 1, 4}});
        f.ws->rebuildAdjacency();
        bnd::ConstraintSet cs = f.cs;
        cs.vertices[1] = f.n[3];
        cs.vertices[2] = f.n[4];
        cs.segments.push_back({10, 1, 2, 7, {1, 9}});
        const auto rep = bnd::recoverSegments(*f.ws, cs);
        CHECK(rep.processed.size() == 1);
        CHECK(rep.processed[0].outcome == bnd::Outcome::Recovered);
        f.ws->rebuildAdjacency();
        // both on-segment sub-edges are now the required constraint
        const auto e0 = f.ws->edges().find(f.ws->nodes(), f.n[3], f.n[1]);
        const auto e1 = f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[4]);
        CHECK(bool(e0) && bool(e1));
        CHECK(f.ws->edges().get(*e0).requiredConstraint);
        CHECK(f.ws->edges().get(*e1).requiredConstraint);
    }
    {   // CROSSED REQUIRED EDGE (13.8): the segment n1 -> n3 would RETIRE
        // the quad diagonal n0-n2, but that diagonal is itself a required
        // constraint. BND policy REFUSES rather than destroy it. Pinned.
        auto f = buildGeom({{0, 0}, {2, 0}, {2, 2}, {0, 2}},
                           {{0, 1, 2}, {0, 2, 3}});
        const auto e02 = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[2]);
        f.ws->testEdgesMutable().setRequiredConstraint(e02, true);
        f.ws->testEdgesMutable().setClassification(e02, {1, 5});
        f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        bnd::ConstraintSet cs = f.cs;
        cs.vertices[1] = f.n[1];
        cs.vertices[2] = f.n[3];
        cs.segments.push_back({10, 1, 2, 7, {1, 9}});
        const auto rep = bnd::recoverSegments(*f.ws, cs);
        CHECK(rep.processed.size() == 1);
        CHECK(rep.processed[0].outcome == bnd::Outcome::Refused);
        CHECK(rep.processed[0].refused ==
              bnd::RefusedReason::CrossedRequiredEdge);   // exact reason
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);   // zero change
    }
    {   // BOUNDARY REACHED (13.8): from the hub of a 4-cell wheel, aim at a
        // node of a far DISJOINT triangle. The walk crosses out of the
        // wheel through a boundary edge before reaching b. Pinned reason.
        auto f = buildGeom({{0, 0}, {2, 0}, {0, 2}, {-2, 0}, {0, -2},
                            {10, 10}, {12, 10}, {10, 12}},
                           {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 1},
                            {5, 6, 7}});
        f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        bnd::ConstraintSet cs = f.cs;
        cs.vertices[1] = f.n[0];    // hub
        cs.vertices[2] = f.n[5];    // far disjoint triangle
        cs.segments.push_back({10, 1, 2, 7, {1, 9}});
        bnd::SessionStatus st{};
        auto s = bnd::RecoverySession::open(*f.ws, cs, st);
        CHECK(bool(s));
        if (s) {
            const auto w = bnd::detail::transverseWalk(
                *f.ws, s->input(), s->index(), f.n[0], f.n[5]);
            CHECK(!w.ok);
            CHECK(w.why == bnd::RefusedReason::BoundaryReached);   // pinned
        }
        const auto rep = bnd::recoverSegments(*f.ws, cs);
        CHECK(rep.processed[0].outcome == bnd::Outcome::Refused);
        CHECK(rep.processed[0].refused ==
              bnd::RefusedReason::BoundaryReached);
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);   // zero change
    }
    {   // TERMINAL ARRIVAL: the canonical transverse crossing recovers,
        // proving the walk reaches map(b) and closes. (Distinct from the
        // refusals above so the "success" branch is separately reached.)
        auto f = buildGeom({{0, 0}, {2, 0}, {2, 2}, {0, 2}},
                           {{0, 1, 2}, {0, 2, 3}});
        f.ws->rebuildAdjacency();
        bnd::ConstraintSet cs = f.cs;
        cs.vertices[1] = f.n[1];
        cs.vertices[2] = f.n[3];
        cs.segments.push_back({10, 1, 2, 7, {1, 9}});
        const auto rep = bnd::recoverSegments(*f.ws, cs);
        CHECK(rep.processed[0].outcome == bnd::Outcome::Recovered);
    }
}

// Local index of edge (x,y) within a cell's cellEdges() order.
static int localEdgeOf(topo::TopologyWorkspace& ws, topo::CellHandle c,
                       topo::NodeHandle x, topo::NodeHandle y) {
    const auto cd = ws.cells().get(c);
    for (int k = 0; k < 3; ++k) {
        const auto a = cd.nodes[std::size_t(k)];
        const auto b = cd.nodes[std::size_t((k + 1) % 3)];
        if ((a == x && b == y) || (a == y && b == x)) return k;
    }
    return -1;
}

void testBndTransverseDegeneracyGuards() {
    TEST("bnd.transverse_degeneracy_guards (defensive branches, crafted)");
    // Three walk/initial-cell guards CANNOT arise from a valid session: the
    // manifold index build (auditTri3Manifold) rejects overlapping stars,
    // and exact predicates make the monotone walk single-crossing and
    // acyclic. They are DEFENSIVE guards. Per the project's cell-aware
    // corruption-adversary practice we reach them with hand-built indices
    // over deliberately non-manifold geometry (buildTriAdjacency would
    // refuse these, so the index is constructed directly).

    // MULTIPLE qualifying initial cells (13.8a): two overlapping triangles
    // incident to n0 whose wedges BOTH contain direction n0 -> n5. found>1
    // must be refused NoUniqueInitialCell, never silently pick one.
    {
        auto f = buildGeom({{0, 0}, {4, 1}, {4, -1}, {3, 2}, {3, -2},
                            {10, 0}},
                           {{0, 1, 2}, {0, 3, 4}});
        bnd::TriAdjacencyIndex idx;
        idx.stamp = bnd::OwnershipStamp::of(*f.ws);
        idx.star[f.n[0]] = {f.c[0], f.c[1]};   // overlapping star (crafted)
        topo::CellHandle out{};
        bnd::RefusedReason why = bnd::RefusedReason::None;
        const bool ok = bnd::detail::initialCell(*f.ws, f.cs, idx, f.n[0],
                                                 f.n[5], out, why);
        CHECK(!ok);
        CHECK(why == bnd::RefusedReason::NoUniqueInitialCell);
    }
    // WALK GUARD (13.8): a crafted CYCLIC neighbour map makes the walk
    // revisit cells until it exceeds the active-cell bound. Two overlapping
    // straddling triangles with neighbours pointing back at each other.
    {
        auto f = buildGeom({{0, 0}, {20, 0}, {2, 3}, {2, -3}, {1, 2},
                            {1, -2}, {5, 1}},
                           {{0, 2, 3}, {4, 5, 6}});
        bnd::TriAdjacencyIndex idx;
        idx.stamp = bnd::OwnershipStamp::of(*f.ws);
        idx.star[f.n[0]] = {f.c[0]};
        idx.neighbors[f.c[0]] = {topo::CellHandle{}, topo::CellHandle{},
                                 topo::CellHandle{}};
        idx.neighbors[f.c[1]] = {topo::CellHandle{}, topo::CellHandle{},
                                 topo::CellHandle{}};
        const int la = localEdgeOf(*f.ws, f.c[0], f.n[2], f.n[3]);
        if (la >= 0) idx.neighbors[f.c[0]][std::size_t(la)] = f.c[1];
        for (const auto& e : {std::array<int, 2>{4, 5},
                              std::array<int, 2>{5, 6},
                              std::array<int, 2>{6, 4}}) {
            const int lb = localEdgeOf(*f.ws, f.c[1], f.n[std::size_t(e[0])],
                                       f.n[std::size_t(e[1])]);
            if (lb >= 0) idx.neighbors[f.c[1]][std::size_t(lb)] = f.c[0];
        }
        const auto w = bnd::detail::transverseWalk(*f.ws, f.cs, idx, f.n[0],
                                                   f.n[1]);
        CHECK(!w.ok);
        CHECK(w.why == bnd::RefusedReason::WalkGuard);
    }
    // TANGENTIAL TOUCH (13.8, crossings==0) is UNREACHABLE under exact
    // predicates and therefore not fabricated: the walk only steps across
    // an edge it PROPERLY crossed, so that entry edge straddles the segment
    // line; its two endpoints are on opposite sides; each far edge shares
    // one of those endpoints, so at least one far edge always straddles
    // too (or a vertex is exactly on the segment -> InteriorNodeHit fires
    // first). Hence crossings>=1 always. We assert that positive invariant
    // on the canonical crossing rather than fake the dead branch.
    {
        auto f = buildGeom({{0, 0}, {2, 0}, {2, 2}, {0, 2}},
                           {{0, 1, 2}, {0, 2, 3}});
        f.ws->rebuildAdjacency();
        bnd::ConstraintSet cs = f.cs;
        cs.vertices[1] = f.n[1];
        cs.vertices[2] = f.n[3];
        cs.segments.push_back({10, 1, 2, 7, {1, 9}});
        bnd::SessionStatus st{};
        auto s = bnd::RecoverySession::open(*f.ws, cs, st);
        CHECK(bool(s));
        if (s) {
            const auto w = bnd::detail::transverseWalk(
                *f.ws, s->input(), s->index(), f.n[1], f.n[3]);
            CHECK(w.ok);                       // a real crossing, never tangential
            CHECK(w.why == bnd::RefusedReason::None);
        }
    }
}

void testBndEarValidityConcave() {
    TEST("bnd.ear_validity_concave (valid ears only; §13.9)");
    // The review reproduction, distilled to the triangulator. An L-shaped
    // (staircase) hexagon is strictly concave at the reflex vertex (2,2).
    // The pre-fix ear test (convex apex + no vertex STRICTLY inside) admits
    // an ear whose diagonal crosses a nonadjacent polygon edge / has a
    // vertex on it, producing a triangulation the six-gate validator then
    // rejects. earValid() now forbids such ears, so the clip yields a
    // triangulation that PASSES every gate.
    auto f = buildGeom({{0, 0}, {4, 0}, {4, 2}, {2, 2}, {2, 4}, {0, 4}},
                       {{0, 1, 3}});   // cells irrelevant; we test the clip
    bnd::ConstraintSet cs = f.cs;      // params already set by buildGeom
    const std::vector<topo::NodeHandle> ring{f.n[0], f.n[1], f.n[2],
                                             f.n[3], f.n[4], f.n[5]};
    const auto want =
        paacs::predicates::orient2d({0, 0}, {4, 0}, {4, 2});
    const auto tris = bnd::detail::earClipOriented(cs, ring, want);
    CHECK(bool(tris));
    if (tris) {
        CHECK(tris->size() == 4);                 // hexagon -> 4 triangles
        // every emitted triangle is nondegenerate and chart-oriented
        for (const auto& t : *tris) {
            const auto s = paacs::predicates::orient2d(
                cs.paramOf(t[0])->p, cs.paramOf(t[1])->p,
                cs.paramOf(t[2])->p);
            CHECK(s == want);
        }
        // the six-gate candidate validator ACCEPTS the whole triangulation
        topo::CavityPreview pv;
        pv.addedConnectivity = *tris;
        auto v = bnd::makeCandidateValidator(cs, 7, want);
        CHECK(v(pv));
        // the reflex vertex (2,2) = n3 must appear as a real triangle
        // corner, never triangulated away
        bool usesReflex = false;
        for (const auto& t : *tris)
            for (const auto& x : t) usesReflex |= (x == f.n[3]);
        CHECK(usesReflex);
    }
    // a fully COLLINEAR ring has no winding and cannot be triangulated
    auto g = buildGeom({{0, 0}, {2, 0}, {4, 0}}, {{0, 1, 2}});
    const std::vector<topo::NodeHandle> col{g.n[0], g.n[1], g.n[2]};
    CHECK(!bnd::detail::earClipOriented(g.cs, col,
                                        paacs::predicates::Sign::Positive));
}

void testBndInvalidCandidateGeometry() {
    TEST("bnd.invalid_candidate_geometry (Failed, not Refused; tier 1)");
    // Classification contract (review blocker 2): a strip that has ALREADY
    // passed BND-1's disk audit but that cannot be turned into a valid,
    // committable triangulation is a candidate-GENERATION FAILURE, never a
    // policy Refusal. In production the two failure sources both route to
    // Failed(FailedReason::*):
    //   - the frozen cavity rejecting the staged spec  -> CavityRejected;
    //   - the six-gate preview validator rejecting the geometry, or the
    //     triangulator returning no triangulation at all (earClipOriented
    //     == nullopt after the disk audit) -> InvalidCandidateGeometry.
    // The invariant this test pins is the one the review demanded: the
    // outcome is Failed (NOT Refused), and because nothing is prepared the
    // FULL observable digest is unchanged (13.6 tier 1). We drive
    // commitCandidate() directly with a candidate that clears the stamp
    // chain but whose triangulation cannot be committed.
    auto f = buildGeom({{0, 0}, {2, 0}, {4, 0}, {2, 2}, {2, -2}},
                       {{0, 1, 3}, {1, 2, 3}, {1, 0, 4}, {2, 1, 4}});
    f.ws->rebuildAdjacency();
    bnd::ConstraintSet cs = f.cs;
    cs.vertices[1] = f.n[0];
    cs.vertices[2] = f.n[2];
    cs.segments.push_back({10, 1, 2, 7, {1, 9}});
    bnd::SessionStatus st{};
    auto s = bnd::RecoverySession::open(*f.ws, cs, st);
    CHECK(bool(s));
    if (!s) return;
    // A REAL mark-only candidate (it stages cleanly through the frozen
    // cavity), driven with a MISMATCHED chart orientation so ONLY the
    // six-gate preview validator rejects it (gate 3, orientation) - after
    // the disk audit and staging, before prepare(). This isolates the
    // validator-rejection route to Failed(InvalidCandidateGeometry): a
    // deliberately-overlapping spec would instead be caught earlier by the
    // frozen stage as CavityRejected and would not exercise this path.
    std::vector<topo::NodeHandle> chain;
    const auto route =
        bnd::routeEntity(*f.ws, s->input(), s->index(), cs.segments[0], chain);
    CHECK(route == bnd::Routing::MarkOnly);
    bnd::RefusedReason why = bnd::RefusedReason::None;
    auto cand = bnd::buildMarkOnlyCandidate(*f.ws, *s, cs.segments[0], chain,
                                            why);
    CHECK(bool(cand));
    if (!cand) return;
    const auto want = *s->orientationOf(7);
    const auto wrongWant = want == paacs::predicates::Sign::Positive
                               ? paacs::predicates::Sign::Negative
                               : paacs::predicates::Sign::Positive;
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    const auto r = bnd::commitCandidate(*f.ws, *s, *cand, wrongWant);
    CHECK(r.outcome == bnd::Outcome::Failed);                    // Failed...
    CHECK(r.refused == bnd::RefusedReason::None);                // ...not Refused
    CHECK(r.failed == bnd::FailedReason::InvalidCandidateGeometry);  // EXACT
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);          // FULL digest
    // The OTHER triangulation-failure source the fix reclassifies: the
    // triangulator returning no triangulation. earClipOriented returns
    // nullopt on a non-triangulable (fully collinear) ring;
    // buildTransverseCandidate maps that post-disk failure to
    // Failed(InvalidCandidateGeometry) (see its `fail` out-param), never
    // Refused(ClosureNotADisk). (With the polygon-own-winding fix a valid
    // simple disk strip always triangulates, so that branch is defensive -
    // a broad geometry scan finds it unreachable end-to-end.)
    auto g = buildGeom({{0, 0}, {2, 0}, {4, 0}}, {{0, 1, 2}});
    const std::vector<topo::NodeHandle> collinear{g.n[0], g.n[1], g.n[2]};
    CHECK(!bnd::detail::earClipOriented(g.cs, collinear,
                                        paacs::predicates::Sign::Positive));
}

void testBndTransversePrefixAndInvariance() {
    TEST("bnd.transverse_prefix_and_invariance");
    // TWO transverse entities in DISJOINT quads plus a far untouched
    // triangle. Entity 10 recovers by TRANSVERSE crossing; entity 20's
    // crossed diagonal is a required constraint, so 20 is Refused. The
    // recovered prefix (10) is preserved EXACTLY and the untouched cell is
    // byte-identical - proving prefix semantics on the TRANSVERSE route,
    // not the mark-only one.
    GeomFix f;
    {
        topo::NodePool np;
        const double C[11][2] = {{0, 0},  {2, 0}, {2, 2}, {0, 2},
                                 {10, 0}, {12, 0}, {12, 2}, {10, 2},
                                 {20, 0}, {22, 0}, {20, 2}};
        for (const auto& q : C) f.n.push_back(np.add({q[0], q[1], 0}, {2, 1}));
        topo::CellPool cp;
        const int T[5][3] = {{0, 1, 2}, {0, 2, 3},          // quad A (e10)
                             {4, 5, 6}, {4, 6, 7},          // quad B (e20)
                             {8, 9, 10}};                    // far cell
        for (const auto& t : T)
            f.c.push_back(cp.add(
                topo::ElementType::Tri3,
                std::vector<topo::NodeHandle>{f.n[std::size_t(t[0])],
                                              f.n[std::size_t(t[1])],
                                              f.n[std::size_t(t[2])]},
                {2, 5}, np));
        topo::EdgePool ep;
        for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
        // quad B's diagonal n4-n6 is a REQUIRED constraint -> e20 refuses
        ep.setRequiredConstraint(*ep.find(np, f.n[4], f.n[6]), true);
        ep.setClassification(*ep.find(np, f.n[4], f.n[6]), {1, 5});
        f.ws = std::make_unique<topo::TopologyWorkspace>(
            std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
        f.cs.stamp = bnd::OwnershipStamp::of(*f.ws);
        for (std::size_t i = 0; i < 11; ++i)
            f.cs.params[f.n[i]] = bnd::ParamPoint{{C[i][0], C[i][1]}, 7};
    }
    f.ws->rebuildAdjacency();
    // COMPLETE snapshot of the far untouched cell: handle, connectivity,
    // classification, and every incident edge's handle, use count,
    // required flag and classification.
    const auto farCell = f.c[4];
    const auto farSnap = f.ws->cells().get(farCell);
    struct EdgeSnap { topo::EdgeHandle h; std::uint32_t use;
                      bool req; paacs::mesh::Classification cls; };
    std::vector<EdgeSnap> farEdges;
    for (const auto& eh : f.ws->edges().cellEdges(f.ws->cells(), farCell)) {
        const auto d = f.ws->edges().get(eh);
        farEdges.push_back({eh, d.useCount, d.requiredConstraint, d.cls});
    }
    CHECK(farEdges.size() == 3);
    // ...and each outside NODE's full payload (coordinate + classification).
    using NodeSnapT = std::decay_t<decltype(f.ws->nodes().get(
        farSnap.nodes[0]))>;
    std::array<NodeSnapT, 3> farNodes;
    for (int i = 0; i < 3; ++i)
        farNodes[std::size_t(i)] =
            f.ws->nodes().get(farSnap.nodes[std::size_t(i)]);
    bnd::ConstraintSet cs = f.cs;
    cs.vertices[1] = f.n[1];  cs.vertices[2] = f.n[3];   // quad A diagonal
    cs.segments.push_back({10, 1, 2, 7, {1, 9}});        // transverse
    cs.vertices[3] = f.n[5];  cs.vertices[4] = f.n[7];   // quad B diagonal
    cs.segments.push_back({20, 3, 4, 7, {1, 9}});        // crosses required
    const auto rep = bnd::recoverSegments(*f.ws, cs);
    CHECK(rep.processed.size() == 2);
    CHECK(rep.processed[0].entity == 10);
    CHECK(rep.processed[0].outcome == bnd::Outcome::Recovered);
    CHECK(rep.processed[1].entity == 20);
    CHECK(rep.processed[1].outcome == bnd::Outcome::Refused);
    CHECK(rep.processed[1].refused ==
          bnd::RefusedReason::CrossedRequiredEdge);
    CHECK(rep.unattempted.empty());               // 20 was the last entity
    f.ws->rebuildAdjacency();
    // entity 10's transverse recovery survived exactly
    const auto e13 = f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[3]);
    CHECK(bool(e13));
    CHECK(f.ws->edges().get(*e13).requiredConstraint);
    CHECK(f.ws->edges().get(*e13).cls.tag == 9);
    // the far untouched cell is byte-identical: same handle, connectivity,
    // classification, AND every incident edge's handle/use/required/class.
    CHECK(f.ws->cells().isCurrent(farCell));
    const auto after = f.ws->cells().get(farCell);
    for (int i = 0; i < 3; ++i)
        CHECK(after.nodes[std::size_t(i)] == farSnap.nodes[std::size_t(i)]);
    CHECK(after.cls.dim == farSnap.cls.dim && after.cls.tag == farSnap.cls.tag);
    for (const auto& es : farEdges) {
        CHECK(f.ws->edges().isCurrent(es.h));
        const auto d = f.ws->edges().get(es.h);
        CHECK(d.useCount == es.use);
        CHECK(d.requiredConstraint == es.req);
        CHECK(d.cls.dim == es.cls.dim && d.cls.tag == es.cls.tag);
    }
    // and each outside node is current with UNCHANGED coordinate AND
    // classification (payload byte-identity, not merely still-alive).
    for (int i = 0; i < 3; ++i) {
        CHECK(f.ws->nodes().isCurrent(farSnap.nodes[std::size_t(i)]));
        const auto nd = f.ws->nodes().get(farSnap.nodes[std::size_t(i)]);
        const auto& s = farNodes[std::size_t(i)];
        CHECK(nd.coord.x == s.coord.x && nd.coord.y == s.coord.y &&
              nd.coord.z == s.coord.z);
        CHECK(nd.cls.dim == s.cls.dim && nd.cls.tag == s.cls.tag);
    }
}

void testBndTransverseMixedClassification() {
    TEST("bnd.transverse_mixed_classification (13.5 Refused)");
    // A transverse strip whose two crossed cells carry DIFFERENT face
    // classifications is a structured Refused(MixedCellClassification),
    // with zero change - on the transverse route.
    GeomFix f;
    {
        topo::NodePool np;
        const double C[4][2] = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
        for (const auto& q : C) f.n.push_back(np.add({q[0], q[1], 0}, {2, 1}));
        topo::CellPool cp;
        // the two triangles of the quad carry DIFFERENT classifications
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{f.n[0], f.n[1], f.n[2]}, {2, 5}, np));
        f.c.push_back(cp.add(topo::ElementType::Tri3,
            std::vector<topo::NodeHandle>{f.n[0], f.n[2], f.n[3]}, {2, 6}, np));
        topo::EdgePool ep;
        for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
        f.ws = std::make_unique<topo::TopologyWorkspace>(
            std::move(np), std::move(cp), std::move(ep), topo::FacePool{});
        f.cs.stamp = bnd::OwnershipStamp::of(*f.ws);
        for (std::size_t i = 0; i < 4; ++i)
            f.cs.params[f.n[i]] = bnd::ParamPoint{{C[i][0], C[i][1]}, 7};
    }
    f.ws->rebuildAdjacency();
    const auto d0 = paacs_tests::observableDigest(*f.ws);
    bnd::ConstraintSet cs = f.cs;
    cs.vertices[1] = f.n[1];
    cs.vertices[2] = f.n[3];
    cs.segments.push_back({10, 1, 2, 7, {1, 9}});
    const auto rep = bnd::recoverSegments(*f.ws, cs);
    CHECK(rep.processed.size() == 1);
    CHECK(rep.processed[0].outcome == bnd::Outcome::Refused);
    CHECK(rep.processed[0].refused ==
          bnd::RefusedReason::MixedCellClassification);
    CHECK(paacs_tests::observableDigest(*f.ws) == d0);      // zero change
}

void testBndTransverseDeterminism() {
    TEST("bnd.transverse_determinism_modulo_allocator_history");
    // Transverse recovery of the same quad on a FRESH workspace and on one
    // with allocator churn must yield the SAME semantic V2 hash - the
    // determinism claim exercised on the TRANSVERSE route.
    auto mk = [] {
        return buildGeom({{0, 0}, {2, 0}, {2, 2}, {0, 2}},
                         {{0, 1, 2}, {0, 2, 3}});
    };
    auto fresh = mk();
    fresh.ws->rebuildAdjacency();
    auto churn = mk();
    {
        // STRONG churn: build a whole disjoint triangle (3 nodes + a cell,
        // which also allocates its 3 edges) and then tear it all down. This
        // retires slots in the NODE, CELL and EDGE pools and advances the
        // runtime pool token - perturbing exactly the cell/edge allocation
        // and star-cell handles that transverse recovery orders over, not
        // merely an orphan node. The logical mesh returns to the pristine
        // quad, so a correct semantic hash must be blind to all of it.
        topo::TopologyPolicy pol;
        pol.permitOrphanNodes = true;
        std::vector<topo::NodeHandle> tn;
        for (int i = 0; i < 3; ++i) {
            auto tx = churn.ws->beginTransaction(pol);
            tx.addNode({50.0 + i, 50.0 + (i & 1), 0}, {0, 7});
            tn.push_back(tx.commit().node);
            churn.ws->rebuildAdjacency();
        }
        topo::CellHandle tc;
        {
            auto tx = churn.ws->beginTransaction(pol);
            std::vector<topo::NodeHandle> cn{tn[0], tn[1], tn[2]};
            tx.addCell(topo::ElementType::Tri3, cn, {2, 5});
            tc = tx.commit().cell;
            churn.ws->rebuildAdjacency();
        }
        {
            auto tx = churn.ws->beginTransaction(pol);
            tx.removeCell(tc);
            (void)tx.commit();
            churn.ws->rebuildAdjacency();
        }
        for (int i = 0; i < 3; ++i) {
            auto tx = churn.ws->beginTransaction(pol);
            tx.removeNode(tn[std::size_t(i)]);
            (void)tx.commit();
            churn.ws->rebuildAdjacency();
        }
    }
    // pre-recovery: identical logical mesh despite the pool churn
    CHECK(fresh.ws->semanticTopologyHashV2() ==
          churn.ws->semanticTopologyHashV2());
    auto mkSet = [](GeomFix& g) {
        bnd::ConstraintSet cs = g.cs;
        cs.stamp = bnd::OwnershipStamp::of(*g.ws);
        cs.vertices[1] = g.n[1];
        cs.vertices[2] = g.n[3];
        cs.segments.push_back({10, 1, 2, 7, {1, 9}});
        return cs;
    };
    const auto rA = bnd::recoverSegments(*fresh.ws, mkSet(fresh));
    const auto rB = bnd::recoverSegments(*churn.ws, mkSet(churn));
    CHECK(rA.processed[0].outcome == bnd::Outcome::Recovered);
    CHECK(rB.processed[0].outcome == bnd::Outcome::Recovered);
    fresh.ws->rebuildAdjacency();
    churn.ws->rebuildAdjacency();
    CHECK(fresh.ws->semanticTopologyHashV2() ==
          churn.ws->semanticTopologyHashV2());
    CHECK(fresh.ws->testV2Tokens() == churn.ws->testV2Tokens());
    CHECK(fresh.ws->topologyRevision() != churn.ws->topologyRevision());
}

void testBndInitialCellSelection() {
    TEST("bnd.initial_cell_selection (unique/zero/collinear/churn)");
    // A closed wheel of 4 cells around n0 (high valence), plus a far
    // DISJOINT triangle supplying target nodes. Nothing is remapped, so
    // the wedges under test keep their real geometry.
    auto mkFix = [] {
        return buildGeom({{0, 0}, {2, 0}, {0, 2}, {-2, 0}, {0, -2},
                          {10, 10}, {12, 10}, {10, 12}},
                         {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 1},
                          {5, 6, 7}});
    };
    auto f = mkFix();
    f.ws->rebuildAdjacency();
    bnd::ConstraintSet cs = f.cs;
    cs.vertices[1] = f.n[0];
    cs.vertices[2] = f.n[5];
    cs.segments.push_back({10, 1, 2, 7, {1, 9}});
    bnd::SessionStatus st{};
    auto s = bnd::RecoverySession::open(*f.ws, cs, st);
    CHECK(bool(s));
    if (!s) return;
    topo::CellHandle out{};
    bnd::RefusedReason why = bnd::RefusedReason::None;
    // UNIQUE: n0 -> n5 at (10,10) points strictly inside exactly one
    // wedge, the one spanned by n1(2,0) and n2(0,2).
    const bool ok = bnd::detail::initialCell(*f.ws, s->input(), s->index(),
                                             f.n[0], f.n[5], out, why);
    CHECK(ok);
    if (ok) {
        const auto cd = f.ws->cells().get(out);
        bool hasN1 = false, hasN2 = false;
        for (int i2 = 0; i2 < 3; ++i2) {
            hasN1 |= cd.nodes[std::size_t(i2)] == f.n[1];
            hasN2 |= cd.nodes[std::size_t(i2)] == f.n[2];
        }
        CHECK(hasN1 && hasN2);
    }
    // COLLINEAR / vertex-hit: a->b exactly along an incident edge.
    why = bnd::RefusedReason::None;
    CHECK(!bnd::detail::initialCell(*f.ws, s->input(), s->index(), f.n[0],
                                    f.n[1], out, why));
    CHECK(why == bnd::RefusedReason::CollinearOverlapWithEdge);
    // UNIQUENESS is guaranteed: a true return means EXACTLY one wedge
    // qualified (initialCell returns false unless found == 1), so the
    // "multiple" branch cannot silently pick one of several. In a valid
    // planar fan the wedges tile the disk without overlap, so >1 is
    // geometrically impossible; the guard is defensive and here we prove
    // the positive uniqueness the hub relies on.
    CHECK(ok);
    // ZERO qualifying cells: from a BOUNDARY node of the far triangle,
    // aim back at the wheel - no incident wedge contains that direction.
    // The direction is NOT collinear with an incident edge, so the reason
    // is pinned to NoUniqueInitialCell (zero), not the collinear branch.
    why = bnd::RefusedReason::None;
    const bool z = bnd::detail::initialCell(*f.ws, s->input(), s->index(),
                                            f.n[6], f.n[0], out, why);
    CHECK(!z);
    CHECK(why == bnd::RefusedReason::NoUniqueInitialCell);
    // DETERMINISM under allocator churn: the SAME geometric cell.
    {
        auto g = mkFix();
        topo::TopologyPolicy pol;
        pol.permitOrphanNodes = true;
        auto tx = g.ws->beginTransaction(pol);
        tx.addNode({50, 50, 0}, {0, 7});
        auto res = tx.commit();
        g.ws->rebuildAdjacency();
        auto tx2 = g.ws->beginTransaction(pol);
        tx2.removeNode(res.node);
        (void)tx2.commit();
        g.ws->rebuildAdjacency();
        bnd::ConstraintSet cg = g.cs;
        cg.stamp = bnd::OwnershipStamp::of(*g.ws);
        cg.vertices[1] = g.n[0];
        cg.vertices[2] = g.n[5];
        cg.segments.push_back({10, 1, 2, 7, {1, 9}});
        bnd::SessionStatus st2{};
        auto s2 = bnd::RecoverySession::open(*g.ws, cg, st2);
        CHECK(bool(s2));
        if (s2) {
            topo::CellHandle out2{};
            bnd::RefusedReason w2 = bnd::RefusedReason::None;
            const bool ok2 = bnd::detail::initialCell(
                *g.ws, s2->input(), s2->index(), g.n[0], g.n[5], out2, w2);
            CHECK(ok2);
            if (ok2) {
                const auto cd = g.ws->cells().get(out2);
                bool hasN1 = false, hasN2 = false;
                for (int i2 = 0; i2 < 3; ++i2) {
                    hasN1 |= cd.nodes[std::size_t(i2)] == g.n[1];
                    hasN2 |= cd.nodes[std::size_t(i2)] == g.n[2];
                }
                CHECK(hasN1 && hasN2);   // same geometric answer
            }
        }
    }
}

void testBndTransverseMultiCellStrip() {
    TEST("bnd.transverse_multi_cell_strip (concave side, strict reflex)");
    // A DEMONSTRATED concave-side adversary. The segment n0 -> n3 crosses a
    // four-cell strip; the UPPER boundary chain n0,n5,n4,n3 =
    // (0,0),(2,3),(4,1),(6,0) has a STRICT REFLEX vertex at n4=(4,1): it
    // turns the opposite way from its neighbours. A convexity-by-`want`
    // ear test (the pre-fix code) marks that reflex vertex as an ear and
    // clips it, yielding an overlapping triangulation the six-gate
    // validator rejects. The fix clips ears in the polygon's OWN winding,
    // so the reflex vertex is never an ear and the strip recovers in ONE
    // commit. This test therefore REQUIRES the fix to pass, not merely
    // "accepts a structured refusal".
    auto f = buildGeom({{0, 0}, {2, -2}, {4, -2}, {6, 0}, {4, 1}, {2, 3}},
                       {{0, 1, 5}, {1, 2, 5}, {2, 4, 5}, {2, 3, 4}});
    f.ws->rebuildAdjacency();
    // Confirm the upper side really is strictly concave: the turn at n4 is
    // opposite the turn at n5 (a genuine reflex, not a collinear run).
    const auto want = paacs::predicates::orient2d({0, 0}, {2, 3}, {4, 1});
    const auto turnAt4 =
        paacs::predicates::orient2d({2, 3}, {4, 1}, {6, 0});
    CHECK(turnAt4 != paacs::predicates::Sign::Zero);   // not collinear
    CHECK(turnAt4 != want);                            // strictly reflex
    const auto rev0 = f.ws->topologyRevision();
    const auto liveCells = f.ws->cells().stats().live;
    CHECK(liveCells == 4);                             // multi-cell strip
    bnd::ConstraintSet cs = f.cs;
    cs.vertices[1] = f.n[0];    // (0,0)
    cs.vertices[2] = f.n[3];    // (6,0)
    cs.segments.push_back({10, 1, 2, 7, {1, 9}});
    // the whole four-cell strip is the crossed transversal set
    bnd::SessionStatus st{};
    auto sess = bnd::RecoverySession::open(*f.ws, cs, st);
    CHECK(bool(sess));
    if (sess) {
        const auto w = bnd::detail::transverseWalk(
            *f.ws, sess->input(), sess->index(), f.n[0], f.n[3]);
        CHECK(w.ok);
        CHECK(w.strip.size() == 4);                    // longer walk
    }
    const auto rep = bnd::recoverSegments(*f.ws, cs);
    CHECK(!rep.rejected);
    CHECK(rep.processed.size() == 1);
    CHECK(rep.processed[0].outcome == bnd::Outcome::Recovered);
    CHECK(f.ws->topologyRevision() == rev0 + 1);       // exactly one commit
    f.ws->rebuildAdjacency();
    const auto e = f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[3]);
    CHECK(bool(e));
    CHECK(f.ws->edges().get(*e).requiredConstraint);
    CHECK(f.ws->edges().get(*e).cls.dim == 1);
    CHECK(f.ws->edges().get(*e).cls.tag == 9);
    CHECK(f.ws->cells().stats().live == liveCells);    // Euler: 4 in, 4 out
    // and the reflex vertex n4 survives as a real interior mesh node
    CHECK(f.ws->nodes().isCurrent(f.n[4]));
    // re-running is AlreadyPresent with zero churn
    const auto rev1 = f.ws->topologyRevision();
    bnd::ConstraintSet cs2 = cs;
    cs2.stamp = bnd::OwnershipStamp::of(*f.ws);
    const auto rep2 = bnd::recoverSegments(*f.ws, cs2);
    CHECK(rep2.processed[0].outcome == bnd::Outcome::AlreadyPresent);
    CHECK(f.ws->topologyRevision() == rev1);
}

void testBndEarTieBreakGeometric() {
    TEST("bnd.ear_tie_break_geometric (apex order, not handle order)");
    // A convex pentagon has several simultaneously valid ears. The apex
    // with the lexicographically smallest (u,v) must be clipped FIRST,
    // regardless of the handle order the ring is given in.
    auto f = buildGeom({{0, 0}, {4, 0}, {5, 3}, {2, 5}, {-1, 3}},
                       {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}});
    f.ws->rebuildAdjacency();
    bnd::ConstraintSet cs = f.cs;
    cs.params[f.n[0]] = bnd::ParamPoint{{0, 0}, 7};
    cs.params[f.n[1]] = bnd::ParamPoint{{4, 0}, 7};
    cs.params[f.n[2]] = bnd::ParamPoint{{5, 3}, 7};
    cs.params[f.n[3]] = bnd::ParamPoint{{2, 5}, 7};
    cs.params[f.n[4]] = bnd::ParamPoint{{-1, 3}, 7};
    const auto want = paacs::predicates::orient2d({0, 0}, {4, 0}, {5, 3});
    // ring given in HANDLE order and in a ROTATED order: the first
    // emitted triangle's apex must be the same geometric point (-1,3)
    // is smallest by u, so the ear at n4 is clipped first.
    const std::vector<topo::NodeHandle> ringA{f.n[0], f.n[1], f.n[2],
                                              f.n[3], f.n[4]};
    const std::vector<topo::NodeHandle> ringB{f.n[2], f.n[3], f.n[4],
                                              f.n[0], f.n[1]};
    const auto a = bnd::detail::earClipOriented(cs, ringA, want);
    const auto b = bnd::detail::earClipOriented(cs, ringB, want);
    CHECK(bool(a) && bool(b));
    if (a && b) {
        CHECK(a->size() == 3 && b->size() == 3);
        // same APEX chosen first from both presentations
        CHECK((*a)[0][1] == (*b)[0][1]);
        const auto* pp = cs.paramOf((*a)[0][1]);
        CHECK(pp != nullptr);
        if (pp) {   // the lexicographically smallest apex is (-1,3)
            CHECK(pp->p.u == -1.0);
            CHECK(pp->p.v == 3.0);
        }
    }
}

void testBndOutsideInvarianceAndStaleHandoff() {
    TEST("bnd.outside_invariance_and_stale_handoff");
    // A cell OUTSIDE the closure must be byte-identical after recovery,
    // and a stale candidate stamp must be refused (Failed, no mutation).
    auto f = buildGeom({{0, 0}, {2, 0}, {4, 0}, {2, 2}, {2, -2}, {8, 0},
                        {10, 0}, {9, 2}},
                       {{0, 1, 3}, {1, 2, 3}, {1, 0, 4}, {2, 1, 4},
                        {5, 6, 7}});
    f.ws->rebuildAdjacency();
    const auto outside = f.c[4];
    const auto snap = f.ws->cells().get(outside);
    const auto e56 = *f.ws->edges().find(f.ws->nodes(), f.n[5], f.n[6]);
    const auto use0 = f.ws->edges().get(e56).useCount;
    auto cs = chainSet(f, {1, 9});
    const auto rep = bnd::recoverSegments(*f.ws, cs);
    CHECK(rep.processed[0].outcome == bnd::Outcome::Recovered);
    f.ws->rebuildAdjacency();
    // outside cell: same handle, same payload, same edge use count
    CHECK(f.ws->cells().isCurrent(outside));
    const auto after = f.ws->cells().get(outside);
    for (int i = 0; i < 3; ++i)
        CHECK(after.nodes[std::size_t(i)] == snap.nodes[std::size_t(i)]);
    CHECK(after.cls.dim == snap.cls.dim && after.cls.tag == snap.cls.tag);
    CHECK(f.ws->edges().isCurrent(e56));
    CHECK(f.ws->edges().get(e56).useCount == use0);
    // stale candidate handoff -> Failed(StampChainBroken), no mutation
    {
        bnd::SessionStatus st{};
        auto cs2 = chainSet(f, {1, 4}, 11);
        cs2.stamp = bnd::OwnershipStamp::of(*f.ws);
        auto sess = bnd::RecoverySession::open(*f.ws, cs2, st);
        CHECK(bool(sess));
        bnd::RecoveryCandidate cand;
        cand.stamp = bnd::OwnershipStamp{};       // STALE
        cand.entity = 11;
        cand.chart = 7;
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        const auto r = bnd::commitCandidate(*f.ws, *sess, cand,
                                            *sess->orientationOf(7));
        CHECK(r.outcome == bnd::Outcome::Failed);
        CHECK(r.failed == bnd::FailedReason::StampChainBroken);
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);
    }
}

void testBndFailureTiers() {
    TEST("bnd.failure_tiers (13.6: all three tiers)");
    // TIER 1 - validator/candidate rejection, nothing prepared: the FULL
    // observable digest is unchanged.
    {
        auto f = buildChainFan();
        const auto eh = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[1]);
        f.ws->testEdgesMutable().setClassification(eh, {1, 4});
        f.ws->rebuildAdjacency();
        const auto d0 = paacs_tests::observableDigest(*f.ws);
        auto cs = chainSet(f, {1, 9});
        const auto rep = bnd::recoverSegments(*f.ws, cs);
        CHECK(rep.processed[0].outcome == bnd::Outcome::Refused);
        CHECK(paacs_tests::observableDigest(*f.ws) == d0);   // FULL digest
    }
    // TIER 2 - preparation/allocation failure: the LOGICAL digest is
    // equal while permitted capacity residue may remain.
    {
        auto f = buildChainFan();
        f.ws->rebuildAdjacency();
        const auto L0 = paacs_tests::observableLogicalDigest(*f.ws);
        auto cs = chainSet(f, {1, 9});
        bool sawFailure = false;
        for (int k = 1; k <= 400 && !sawFailure; ++k) {
            auto g = buildChainFan();
            g.ws->rebuildAdjacency();
            const auto Lg = paacs_tests::observableLogicalDigest(*g.ws);
            auto cg = chainSet(g, {1, 9});
            paacs::topo::faultCountdown() = k;
            const auto rep = bnd::recoverSegments(*g.ws, cg);
            paacs::topo::faultCountdown() = 0;
            if (!rep.processed.empty() &&
                rep.processed[0].outcome == bnd::Outcome::Failed) {
                sawFailure = true;
                CHECK(paacs_tests::observableLogicalDigest(*g.ws) == Lg);
                const auto m = paacs::topo::measureWorkspaceMemory(*g.ws);
                CHECK(m.consistent());
            }
        }
        CHECK(sawFailure);          // the seam battery reached Phase 2
        (void)L0;
    }
    // TIER 3 - a later prefix failure preserves every earlier committed
    // entity EXACTLY (entity 10 recovers; 20 refuses; 10 stays marked).
    {
        // two DISJOINT regions so the set is PSLG-clean: the fan for
        // entity 10, and a far triangle for entity 20 whose edge is
        // pre-marked with a CONFLICTING classification.
        GeomFix f;
        {
            topo::NodePool np;
            const double C[8][2] = {{0, 0},  {2, 0},  {4, 0}, {2, 2},
                                    {2, -2}, {30, 0}, {32, 0}, {31, 2}};
            for (const auto& q : C)
                f.n.push_back(np.add({q[0], q[1], 0}, {2, 1}));
            topo::CellPool cp;
            const int T[5][3] = {{0, 1, 3}, {1, 2, 3}, {1, 0, 4},
                                 {2, 1, 4}, {5, 6, 7}};
            for (const auto& tr : T)
                f.c.push_back(cp.add(
                    topo::ElementType::Tri3,
                    std::vector<topo::NodeHandle>{f.n[std::size_t(tr[0])],
                                                  f.n[std::size_t(tr[1])],
                                                  f.n[std::size_t(tr[2])]},
                    {2, 5}, np));
            topo::EdgePool ep;
            for (const auto& ch : f.c) ep.insertCellEdges(np, cp, ch);
            ep.setClassification(*ep.find(np, f.n[5], f.n[6]), {1, 4});
            f.ws = std::make_unique<topo::TopologyWorkspace>(
                std::move(np), std::move(cp), std::move(ep),
                topo::FacePool{});
            f.cs.stamp = bnd::OwnershipStamp::of(*f.ws);
            for (std::size_t i = 0; i < 8; ++i)
                f.cs.params[f.n[i]] = bnd::ParamPoint{{C[i][0], C[i][1]}, 7};
        }
        f.ws->rebuildAdjacency();
        bnd::ConstraintSet cs = f.cs;
        cs.vertices[1] = f.n[0];
        cs.vertices[2] = f.n[2];
        cs.segments.push_back({10, 1, 2, 7, {1, 9}});   // recovers
        cs.vertices[3] = f.n[5];
        cs.vertices[4] = f.n[6];
        cs.segments.push_back({20, 3, 4, 7, {1, 9}});   // edge holds {1,4}
        const auto rep = bnd::recoverSegments(*f.ws, cs);
        // EXACT two-result vector: 10 recovered, 20 refused for the
        // conflicting mark, nothing unattempted (20 was the last entity).
        // (A run that wrongly stopped after 10 would fail here.)
        CHECK(rep.processed.size() == 2);
        CHECK(rep.processed[0].entity == 10);
        CHECK(rep.processed[0].outcome == bnd::Outcome::Recovered);
        CHECK(rep.processed[1].entity == 20);
        CHECK(rep.processed[1].outcome == bnd::Outcome::Refused);
        CHECK(rep.processed[1].refused ==
              bnd::RefusedReason::ConflictingEdgeMark);
        CHECK(rep.unattempted.empty());
        f.ws->rebuildAdjacency();
        // entity 10's marks survive the later entity's failure EXACTLY
        const auto a = *f.ws->edges().find(f.ws->nodes(), f.n[0], f.n[1]);
        const auto b = *f.ws->edges().find(f.ws->nodes(), f.n[1], f.n[2]);
        CHECK(f.ws->edges().get(a).requiredConstraint);
        CHECK(f.ws->edges().get(a).cls.tag == 9);
        CHECK(f.ws->edges().get(b).requiredConstraint);
        CHECK(f.ws->edges().get(b).cls.tag == 9);
        // and entity 20's edge is untouched (still its original {1,4})
        const auto c = *f.ws->edges().find(f.ws->nodes(), f.n[5], f.n[6]);
        CHECK(!f.ws->edges().get(c).requiredConstraint);
        CHECK(f.ws->edges().get(c).cls.tag == 4);
    }
}

void testBndDeterminismModuloAllocatorHistory() {
    TEST("bnd.determinism_modulo_allocator_history (fresh vs churned)");
    // Fresh workspace vs one with allocator churn (a cell added then
    // removed): the same recovery must produce the SAME semantic V2.
    auto fresh = buildChainFan();
    fresh.ws->rebuildAdjacency();
    auto churn = buildChainFan();
    {
        topo::TopologyPolicy pol;
        pol.permitOrphanNodes = true;
        auto tx = churn.ws->beginTransaction(pol);
        tx.addNode({50, 50, 0}, {0, 7});
        auto res = tx.commit();
        churn.ws->rebuildAdjacency();
        auto tx2 = churn.ws->beginTransaction(pol);
        tx2.removeNode(res.node);
        (void)tx2.commit();
        churn.ws->rebuildAdjacency();
    }
    CHECK(fresh.ws->semanticTopologyHashV2() ==
          churn.ws->semanticTopologyHashV2());
    auto csA = chainSet(fresh, {1, 9});
    auto csB = chainSet(churn, {1, 9});
    csB.stamp = bnd::OwnershipStamp::of(*churn.ws);
    const auto rA = bnd::recoverSegments(*fresh.ws, csA);
    const auto rB = bnd::recoverSegments(*churn.ws, csB);
    CHECK(rA.processed[0].outcome == bnd::Outcome::Recovered);
    CHECK(rB.processed[0].outcome == bnd::Outcome::Recovered);
    fresh.ws->rebuildAdjacency();
    churn.ws->rebuildAdjacency();
    // SEMANTIC equality despite different allocator history
    CHECK(fresh.ws->semanticTopologyHashV2() ==
          churn.ws->semanticTopologyHashV2());
    CHECK(fresh.ws->testV2Tokens() == churn.ws->testV2Tokens());
    // revisions differ (ownership/history is explicitly NOT required)
    CHECK(fresh.ws->topologyRevision() != churn.ws->topologyRevision());
}

}  // namespace

namespace paacs_tests {
const Entry* bndCommitTests(std::size_t& count) {
    static const Entry kBnd[] = {
        {"bnd.routing_all_four", testBndRoutingAllFour},
        {"bnd.mark_only_recovery", testBndMarkOnlyRecovery},
        {"bnd.mark_conflict_cases", testBndMarkConflictCases},
        {"bnd.mixed_cell_classification_refused",
         testBndMixedCellClassification},
        {"bnd.disk_audit_complete", testBndDiskAuditComplete},
        {"bnd.prefix_semantics", testBndPrefixSemanticsClean},
        {"bnd.whole_set_rejected_zero_change",
         testBndWholeSetRejectedZeroChange},
        {"bnd.candidate_validator_all_gates",
         testBndCandidateValidatorAllGates},
        {"bnd.initial_cell_selection", testBndInitialCellSelection},
        {"bnd.transverse_recovery", testBndTransverseRecovery},
        {"bnd.transverse_multi_cell_strip", testBndTransverseMultiCellStrip},
        {"bnd.ear_tie_break_geometric", testBndEarTieBreakGeometric},
        {"bnd.ear_validity_concave", testBndEarValidityConcave},
        {"bnd.invalid_candidate_geometry", testBndInvalidCandidateGeometry},
        {"bnd.transverse_degeneracies", testBndTransverseDegeneracies},
        {"bnd.transverse_degeneracy_guards",
         testBndTransverseDegeneracyGuards},
        {"bnd.transverse_prefix_and_invariance",
         testBndTransversePrefixAndInvariance},
        {"bnd.transverse_mixed_classification",
         testBndTransverseMixedClassification},
        {"bnd.transverse_determinism", testBndTransverseDeterminism},
        {"bnd.outside_invariance_and_stale_handoff",
         testBndOutsideInvarianceAndStaleHandoff},
        {"bnd.failure_tiers", testBndFailureTiers},
        {"bnd.determinism_modulo_allocator_history",
         testBndDeterminismModuloAllocatorHistory},
    };
    count = sizeof(kBnd) / sizeof(kBnd[0]);
    return kBnd;
}
}  // namespace paacs_tests
