// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// PAACS Mesh unit tests - topology layer (M4-B.2 storage, M4-B.3 edges,
// M4-B.4 faces). Split from test_main.cpp so the monolithic TU stays
// compilable in constrained review environments (reviewer request).
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <random>
#include <span>
#include <type_traits>
#include <vector>

#include "paacs/topo/ElementTraits.hpp"
#include "paacs/topo/Handles.hpp"
#include "paacs/topo/Storage.hpp"
#include "paacs/topo/EdgeTopology.hpp"
#include "paacs/topo/FaceTopology.hpp"
#include "paacs/topo/Workspace.hpp"

#include "test_harness.hpp"

using namespace paacs;

// ---------------------------------------------------------------------------
static void testStorageNodeCell() {
    TEST("topo: M4-B.2 generic node/cell storage");
    using topo::NodePool;
    using topo::CellPool;
    using topo::ElementType;
    using topo::NodeHandle;
    using topo::CellHandle;

    // (1) Node allocation, access, free, reuse.
    {
        NodePool np;
        const auto a = np.add({1, 2, 3}, {0, 5});
        CHECK(np.isCurrent(a));
        CHECK(np.get(a).coord.x == 1 && np.get(a).coord.z == 3);
        CHECK(np.get(a).cls.dim == 0 && np.get(a).cls.tag == 5);
        np.setCoord(a, {4, 5, 6});
        CHECK(np.get(a).coord.y == 5);
        np.remove(a);
        CHECK(!np.isCurrent(a));
        const auto b = np.add({7, 8, 9}, {1, 2});
        CHECK(np.isCurrent(b) && !np.isCurrent(a));
        CHECK(np.stats().consistent());
    }

    // (2) Mixed-cell allocation for all ten types (handle-owned).
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 12; ++i)
            n.push_back(np.add({double(i), double(i * 2), double(i * 3)},
                               {3, 1}));
        CellPool cp;
        const std::pair<ElementType, int> types[] = {
            {ElementType::Segment2, 2}, {ElementType::Segment3, 3},
            {ElementType::Tri3, 3},     {ElementType::Tri6, 6},
            {ElementType::Quad4, 4},    {ElementType::Tet4, 4},
            {ElementType::Tet10, 10},   {ElementType::Hex8, 8},
            {ElementType::Wedge6, 6},   {ElementType::Pyramid5, 5},
        };
        for (const auto& [ty, k] : types) {
            std::vector<NodeHandle> conn(n.begin(), n.begin() + k);
            const auto ch = cp.add(ty, conn, {3, 1}, np);
            CHECK(cp.isCurrent(ch));
            CHECK(cp.get(ch).type == ty);
            CHECK(int(cp.get(ch).count) == k);
        }
        CHECK(cp.stats().live == 10 && cp.stats().consistent());
        CHECK(cp.boundNodePool() == np.poolId());
    }

    // (3) Wrong arity rejected BEFORE mutation.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 3; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        const auto before = cp.stats();
        bool threw = false;
        try { cp.add(ElementType::Tet4, n, {3,1}, np); }   // 3 handles, tet=4
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
        const auto after = cp.stats();
        CHECK(after.materializedSlots == before.materializedSlots &&
              after.totalCreated == before.totalCreated);
    }

    // (4) Foreign / stale node handles rejected.
    {
        NodePool owner, other;
        std::vector<NodeHandle> good;
        for (int i = 0; i < 4; ++i) good.push_back(owner.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto foreign = other.add({9,9,9}, {3,1});
        auto bad = good; bad[3] = foreign;
        bool f1 = false;
        try { cp.add(ElementType::Tet4, bad, {3,1}, owner); }
        catch (const std::invalid_argument&) { f1 = true; }
        CHECK(f1);
        auto stale = good[0];
        owner.remove(good[0]);
        auto bad2 = good; bad2[0] = stale;
        bool f2 = false;
        try { cp.add(ElementType::Tet4, bad2, {3,1}, owner); }
        catch (const std::invalid_argument&) { f2 = true; }
        CHECK(f2);
        CHECK(cp.stats().live == 0);
    }

    // (5) Double free deterministic.
    {
        NodePool np;
        auto a = np.add({0,0,0}, {0,0});
        np.remove(a);
        bool threw = false;
        try { np.remove(a); } catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
    }

    // (6) Retirement THROUGH the pool: retired slot unreadable + unreused.
    {
        NodePool np;
        auto a = np.add({1,1,1}, {0,7});
        const auto slot = a.slot;
        np.testRetire(a);                       // generation-exhaust + free
        CHECK(np.stats().retired == 1 && np.stats().consistent());
        bool badRead = false;
        try { (void)np.get(a); } catch (const std::logic_error&) { badRead = true; }
        CHECK(badRead);
        auto b = np.add({2,2,2}, {0,8});         // must not reuse retired slot
        CHECK(b.slot != slot && np.isCurrent(b));
    }

    // (7) Sequence-ordered iteration survives heavy churn.
    {
        NodePool np;
        std::vector<NodeHandle> live;
        std::mt19937 rng(123);
        for (int i = 0; i < 3000; ++i) {
            if (!live.empty() && (rng() & 3) == 0) {
                const std::size_t k = rng() % live.size();
                np.remove(live[k]);
                live.erase(live.begin() + long(k));
            } else {
                live.push_back(np.add({double(i),0,0}, {0,0}));
            }
        }
        const auto order = np.activeBySequence();
        CHECK(order.size() == live.size());
        for (std::size_t i = 1; i < order.size(); ++i)
            CHECK(np.seqOf(order[i]) > np.seqOf(order[i-1]));
        CHECK(np.stats().consistent());
    }

    // (8) Whole-pool move construction preserves destination handles.
    {
        NodePool src;
        auto a = src.add({1,1,1}, {0,1});
        auto b = src.add({2,2,2}, {0,2});
        NodePool dst(std::move(src));
        CHECK(dst.isCurrent(a) && dst.isCurrent(b));
        CHECK(dst.get(a).coord.x == 1 && dst.get(b).coord.y == 2);
        CHECK(dst.activeBySequence().size() == 2);
    }

    // (9) Moved-from pool explicitly empty and unusable.
    {
        NodePool src;
        (void)src.add({0,0,0}, {0,0});
        NodePool dst(std::move(src));
        CHECK(src.activeBySequence().empty());
        CHECK(src.stats().materializedSlots == 0 &&
              src.stats().live == 0 && src.stats().consistent());
        bool threw = false;
        try { (void)src.add({1,1,1}, {0,0}); }
        catch (const std::logic_error&) { threw = true; }
        CHECK(threw);
    }

    // (10) Counters exact after every op incl. validation exception.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 3; ++i) n.push_back(np.add({double(i),0,0}, {2,1}));
        CellPool cp;
        CHECK(cp.stats().consistent() && cp.stats().live == 0);
        cp.add(ElementType::Tri3, n, {2,1}, np);
        CHECK(cp.stats().live == 1 && cp.stats().peakLive == 1 &&
              cp.stats().consistent());
        bool threw = false;
        try {
            std::vector<NodeHandle> two(n.begin(), n.begin() + 2);
            cp.add(ElementType::Tri3, two, {2,1}, np);   // bad arity
        } catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
        CHECK(cp.stats().live == 1 && cp.stats().consistent());
    }

    // (11) ALLOCATION-FAULT injection: a staging throw before alloc leaves
    // the pool exactly consistent - no unreachable live entity (review 1).
    {
        for (int stage = 1; stage <= 2; ++stage) {
            NodePool np;
            (void)np.add({0,0,0}, {0,0});        // one live to start
            const auto before = np.stats();
            topo::faultCountdown() = stage;
            bool threw = false;
            try { (void)np.add({1,1,1}, {0,0}); }
            catch (const std::bad_alloc&) { threw = true; }
            topo::faultCountdown() = 0;
            CHECK(threw);
            const auto after = np.stats();
            CHECK(after.live == before.live);            // no leaked entity
            CHECK(after.live == np.activeBySequence().size());
            CHECK(after.materializedSlots == before.materializedSlots);
            CHECK(after.consistent());
        }
        for (int stage = 1; stage <= 2; ++stage) {
            NodePool np;
            std::vector<NodeHandle> n;
            for (int i = 0; i < 3; ++i) n.push_back(np.add({double(i),0,0}, {2,1}));
            CellPool cp;
            cp.add(ElementType::Tri3, n, {2,1}, np);
            const auto before = cp.stats();
            topo::faultCountdown() = stage;
            bool threw = false;
            try { cp.add(ElementType::Tri3, n, {2,1}, np); }
            catch (const std::bad_alloc&) { threw = true; }
            topo::faultCountdown() = 0;
            CHECK(threw);
            const auto after = cp.stats();
            CHECK(after.live == before.live);
            CHECK(after.live == cp.activeBySequence().size());
            CHECK(after.consistent());
        }
    }

    // (12) Dangling connectivity DETECTABLE + hash FAILS (review 2, 3).
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto cell = cp.add(ElementType::Tet4, n, {3,1}, np);
        CHECK(cp.connectivityIntact(np));
        np.remove(n[0]);                         // delete a referenced node
        CHECK(!cp.connectivityIntact(np));       // detectable now
        CHECK(cp.isCurrent(cell));               // cell still exists...
        bool hashThrew = false;
        try { (void)topo::contentHash(np, cp); } // ...but hashing rejects it
        catch (const std::invalid_argument&) { hashThrew = true; }
        CHECK(hashThrew);
    }

    // (13) contentHash: interleaving-invariant, but ORDER-sensitive (the
    // documented ordered-storage contract, review 4).
    {
        auto build = [](bool churn) {
            auto np = std::make_unique<NodePool>();
            auto cp = std::make_unique<CellPool>();
            if (churn) {
                std::vector<NodeHandle> junk;
                for (int i = 0; i < 10; ++i) junk.push_back(np->add({99,99,99},{9,9}));
                for (auto& h : junk) np->remove(h);
            }
            std::vector<NodeHandle> n;
            for (int i = 0; i < 4; ++i) n.push_back(np->add({double(i),0,0},{3,1}));
            cp->add(ElementType::Tet4, n, {3,1}, *np);
            return std::pair{std::move(np), std::move(cp)};
        };
        auto [na, ca] = build(false);
        auto [nb, cb] = build(true);
        CHECK(topo::contentHash(*na,*ca).value ==
              topo::contentHash(*nb,*cb).value);   // interleaving-invariant
        // Different creation order -> different ordered-storage hash.
        auto nc = std::make_unique<NodePool>();
        auto cc = std::make_unique<CellPool>();
        std::vector<NodeHandle> rev;
        for (int i = 3; i >= 0; --i) rev.insert(rev.begin(),
                                                nc->add({double(i),0,0},{3,1}));
        // rev now holds handles for coords 0..3 but created in reverse.
        cc->add(ElementType::Tet4, rev, {3,1}, *nc);
        CHECK(topo::contentHash(*na,*ca).value !=
              topo::contentHash(*nc,*cc).value);    // order is semantic
    }

    // (14) Pool tokens absent from content: empty pools with distinct pool
    // ids hash identically.
    {
        NodePool n1, n2;
        CellPool c1, c2;
        CHECK(n1.poolId() != n2.poolId());
        CHECK(topo::contentHash(n1,c1).value == topo::contentHash(n2,c2).value);
    }

    // (15) Memory counters: peaks + allocator memory included (review 6).
    {
        NodePool np;
        std::vector<NodeHandle> hs;
        for (int i = 0; i < 100; ++i) hs.push_back(np.add({double(i),0,0},{0,0}));
        const auto peak = np.stats();
        CHECK(peak.peakLive == 100 && peak.peakMaterialized >= 100);
        for (int i = 0; i < 50; ++i) np.remove(hs[size_t(i)]);
        const auto after = np.stats();
        CHECK(after.live == 50);
        CHECK(after.peakLive == 100);            // peak retained
        CHECK(after.allocatorBytes > 0);
        CHECK(after.totalBytes >= after.payloadBytes + after.allocatorBytes -
              1);                                // total includes allocator
        CHECK(after.consistent());
    }
}


// ---------------------------------------------------------------------------
static void testStorageAuditAndSafety() {
#define CHECK_NOTHROW_HASH(NP, CP) do { \
    bool _ok = true; try { (void)topo::contentHash((NP),(CP)); } \
    catch (...) { _ok = false; } CHECK(_ok); } while (0)
    TEST("topo: M4-B.2 storage audit, invariants, exception contract");
    using topo::NodePool;
    using topo::CellPool;
    using topo::ElementType;
    using topo::NodeHandle;
    using topo::CellHandle;

    // (A) Invalid ElementType is rejected structurally, not via overflow.
    {
        CHECK(!topo::isValidElementType(static_cast<ElementType>(255)));
        CHECK(topo::isValidElementType(ElementType::Tet4));
        bool threw = false;
        try { (void)topo::checkedTraits(static_cast<ElementType>(255)); }
        catch (const std::out_of_range&) { threw = true; }
        CHECK(threw);
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        bool addThrew = false;
        try { cp.add(static_cast<ElementType>(200), n, {3,1}, np); }
        catch (const std::out_of_range&) { addThrew = true; }
        CHECK(addThrew);
        CHECK(cp.stats().live == 0 && cp.checkInvariants());
    }

    // (B) Pool-level invariants: retirement THROUGH the pool keeps the pool
    // internally consistent (review item 2 - the seam previously left the
    // live-list stale while allocator said live=0).
    {
        NodePool np;
        auto a = np.add({1,1,1}, {0,7});
        auto b = np.add({2,2,2}, {0,8});
        CHECK(np.checkInvariants());
        np.testRetire(a);
        CHECK(np.stats().retired == 1);
        CHECK(np.checkInvariants());            // <-- the real fix
        CHECK(np.activeBySequence().size() == np.stats().live);
        CHECK(np.stats().live == 1);            // only b remains
        CHECK(np.isCurrent(b) && !np.isCurrent(a));
    }

    // (C) checkInvariants is stronger than allocator consistent(): a
    // healthy pool passes both; churn keeps both true.
    {
        NodePool np;
        std::vector<NodeHandle> live;
        std::mt19937 rng(9);
        for (int i = 0; i < 1500; ++i) {
            if (!live.empty() && (rng() & 3) == 0) {
                const std::size_t k = rng() % live.size();
                np.remove(live[k]);
                live.erase(live.begin() + long(k));
            } else live.push_back(np.add({double(i),0,0}, {0,0}));
            if ((i & 255) == 0) CHECK(np.checkInvariants());
        }
        CHECK(np.checkInvariants());
    }

    // (D) Honest failure contract: a failed insert preserves LOGICAL state
    // (handles, order, hash, accounting) though capacity may grow (item 3).
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        cp.add(ElementType::Tet4, n, {3,1}, np);
        const auto hashBefore = topo::contentHash(np, cp).value;
        const auto liveBefore = cp.stats().live;
        const auto orderBefore = cp.activeBySequence();
        topo::faultCountdown() = 1;              // fail staging stage 1
        bool threw = false;
        try { cp.add(ElementType::Tet4, n, {3,1}, np); }
        catch (const std::bad_alloc&) { threw = true; }
        topo::faultCountdown() = 0;
        CHECK(threw);
        CHECK(cp.stats().live == liveBefore);        // logical state intact
        CHECK(cp.activeBySequence() == orderBefore);
        CHECK(topo::contentHash(np, cp).value == hashBefore);
        CHECK(cp.checkInvariants());
    }

    // (E) Allocator-originated insertion failure (review item 3/5): drive
    // the allocator's creation sequence to max so alloc_.alloc() itself
    // throws; the pool must remain consistent with no leaked entity.
    {
        NodePool np;
        auto keep = np.add({0,0,0}, {0,0});
        (void)keep;
        const auto before = np.stats();
        np.testForceSeqMax();                    // next alloc() throws
        bool threw = false;
        try { (void)np.add({1,1,1}, {0,0}); }
        catch (const std::length_error&) { threw = true; }
        CHECK(threw);
        const auto after = np.stats();
        CHECK(after.live == before.live);
        CHECK(after.live == np.activeBySequence().size());
        CHECK(np.checkInvariants());
    }

    // (F) Storage audit gate: clean pair is edge-ready; a dangling
    // reference makes it NOT edge-ready with a structured diagnostic
    // (review items 4 & 8).
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        cp.add(ElementType::Tet4, n, {3,1}, np);
        CHECK(topo::storageEdgeReady(np, cp));
        CHECK(topo::auditStorage(np, cp).items.empty());
        np.remove(n[0]);                         // strand the cell
        CHECK(!topo::storageEdgeReady(np, cp));
        const auto rep = topo::auditStorage(np, cp);
        bool sawDangling = false;
        for (const auto& d : rep.items)
            sawDangling |= d.code == diag::Code::DanglingCellConnectivity;
        CHECK(sawDangling);
    }

    // (G) Cross-pool binding mismatch is caught by the audit.
    {
        NodePool a, b;
        std::vector<NodeHandle> na;
        for (int i = 0; i < 4; ++i) na.push_back(a.add({double(i),0,0}, {3,1}));
        CellPool cp;
        cp.add(ElementType::Tet4, na, {3,1}, a);   // bound to a
        CHECK(topo::storageEdgeReady(a, cp));
        CHECK(!topo::storageEdgeReady(b, cp));     // wrong node pool
        bool mism = false;
        for (const auto& d : topo::auditStorage(b, cp).items)
            mism |= d.code == diag::Code::CellPoolNodePoolMismatch;
        CHECK(mism);
    }

    // (I) CellPool retirement through the pool keeps it consistent
    // (review item 7 - retirement was only tested on NodePool before).
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto a = cp.add(ElementType::Tet4, n, {3,1}, np);
        auto b = cp.add(ElementType::Tri3,
                        {n[0], n[1], n[2]}, {2,1}, np);
        CHECK(cp.checkInvariants());
        cp.testRetire(a);
        CHECK(cp.stats().retired == 1);
        CHECK(cp.checkInvariants());
        CHECK(cp.isCurrent(b) && !cp.isCurrent(a));
        CHECK(cp.activeBySequence().size() == cp.stats().live);
    }

    // (J) contentHash REQUIRES edge-ready storage (review item): a dangling
    // reference makes it throw up-front, not discover corruption mid-hash.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        cp.add(ElementType::Tet4, n, {3,1}, np);
        CHECK_NOTHROW_HASH(np, cp);              // ok while edge-ready
        np.remove(n[0]);
        bool threw = false;
        try { (void)topo::contentHash(np, cp); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }

    // (H) CellPool move + combined move preserve the runtime relationship
    // and hashing (review item 5).
    {
        auto np = std::make_unique<NodePool>();
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np->add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto cell = cp.add(ElementType::Tet4, n, {3,1}, *np);
        const auto hashBefore = topo::contentHash(*np, cp).value;
        // Move the CellPool: the live CellHandle stays current.
        CellPool cp2(std::move(cp));
        CHECK(cp2.isCurrent(cell));
        CHECK(cp2.get(cell).type == ElementType::Tet4);
        CHECK(cp.activeBySequence().empty());      // moved-from empty
        CHECK(cp.boundNodePool() == 0);            // and unbound
        bool movedAddThrew = false;
        try { cp.add(ElementType::Tet4, n, {3,1}, *np); }
        catch (const std::logic_error&) { movedAddThrew = true; }
        CHECK(movedAddThrew);                      // moved-from unusable
        // Move the NodePool: connectivity remains intact against the moved
        // destination, and the combined hash is unchanged.
        auto np2 = std::make_unique<NodePool>(std::move(*np));
        CHECK(cp2.connectivityIntact(*np2));
        CHECK(topo::contentHash(*np2, cp2).value == hashBefore);
        CHECK(topo::storageEdgeReady(*np2, cp2));
    }
}


// ---------------------------------------------------------------------------
static void testStorageScale() {
    TEST("topo: M4-B.2 storage scale (linear construction, O(n) audit, O(1) remove)");
    using topo::NodePool;
    using topo::CellPool;
    using topo::ElementType;
    using topo::NodeHandle;
    using topo::CellHandle;
    using Clock = std::chrono::steady_clock;
    auto secs = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };

    // Two sizes to check scaling is linear, not quadratic: if T(2n)/T(n)
    // were ~4 we'd have O(n^2). We allow generous slack for a shared CI
    // core but reject quadratic blow-up.
    auto buildNodes = [&](int n) {
        auto t0 = Clock::now();
        NodePool np;
        for (int i = 0; i < n; ++i) np.add({double(i), double(i % 5), 0}, {0, 0});
        return std::pair{secs(t0, Clock::now()), std::move(np)};
    };
    auto [t50k, np50] = buildNodes(50000);
    auto [t100k, np100] = buildNodes(100000);
    (void)np50;
    // 100k must finish quickly (the old code did NOT finish 100k in 120s).
    CHECK(t100k < 5.0);
    // Sub-quadratic: doubling work must not quadruple time (allow 3x slack).
    if (t50k > 1e-4) CHECK(t100k < t50k * 3.0);

    // 200k mixed cells; audit; O(1) removal.
    const int N = 200000;
    NodePool np;
    std::vector<NodeHandle> nh; nh.reserve(size_t(N));
    for (int i = 0; i < N; ++i)
        nh.push_back(np.add({double(i), double(i % 7), 0}, {3, 1}));
    CellPool cp;
    std::vector<CellHandle> ch; ch.reserve(size_t(N));
    const ElementType cyc[3] = {ElementType::Tet4, ElementType::Tri3,
                                ElementType::Segment2};
    const int arity[3] = {4, 3, 2};
    auto tc0 = Clock::now();
    for (int i = 0; i < N; ++i) {
        const int k = i % 3;
        std::vector<NodeHandle> conn;
        for (int j = 0; j < arity[k]; ++j) conn.push_back(nh[size_t((i + j) % N)]);
        ch.push_back(cp.add(cyc[k], conn, {3, 1}, np));
    }
    const double cellTime = secs(tc0, Clock::now());
    CHECK(cp.stats().live == std::size_t(N));
    CHECK(cellTime < 5.0);

    auto ta0 = Clock::now();
    const bool ready = topo::storageEdgeReady(np, cp);
    const double auditTime = secs(ta0, Clock::now());
    CHECK(ready);
    CHECK(auditTime < 2.0);                 // O(n), fraction of a second

    auto tr0 = Clock::now();
    for (int i = 0; i < 50000; ++i) cp.remove(ch[size_t(i)]);
    const double removeTime = secs(tr0, Clock::now());
    CHECK(removeTime < 1.0);                // O(1) each
    CHECK(cp.checkInvariants() && np.checkInvariants());

    // Documented per-entity memory budget (payload + allocator). Inline
    // connectivity keeps cells allocation-free. NOTE on ceilings: the
    // ENGINEERING target is ~139 B/cell at steady state; at a 200k
    // high-water point geometric vector slack (1.5x growth) pushes the
    // /materialized figure to ~198 B/cell on Linux/libstdc++. The test
    // ceiling below is deliberately looser than that measurement so it
    // stays portable across allocators and slack policies - it guards
    // against an order-of-magnitude regression (e.g. arena retention),
    // not against normal growth slack. Actual measured values are recorded
    // in docs/BASELINES and the round-7 note.
    const double nodeBytes = double(np.stats().totalBytes) / double(N);
    const double cellBytes = double(cp.stats().totalBytes) /
                             double(std::max<std::size_t>(cp.stats().materializedSlots, 1));
    CHECK(nodeBytes < 256.0);               // ~115 observed
    CHECK(cellBytes < 320.0);               // ~198 at 200k high-water

    // Million-operation churn retention: with steady-state ~N live nodes,
    // materialized slots must stay bounded near the live set (free-list
    // reuse), not grow with total creation history.
    {
        NodePool cn;
        std::vector<NodeHandle> livev;
        std::mt19937 rng(7);
        const std::size_t target = 20000;
        for (int i = 0; i < 1000000; ++i) {
            if (livev.size() >= target ||
                (!livev.empty() && (rng() & 1))) {
                const std::size_t k = rng() % livev.size();
                cn.remove(livev[k]);
                livev[k] = livev.back();
                livev.pop_back();
            } else {
                livev.push_back(cn.add({double(i),0,0}, {0,0}));
            }
        }
        // Materialized slots bounded by the high-water live set, NOT by the
        // 1,000,000 total operations.
        CHECK(cn.stats().materializedSlots <= target + 1);
        CHECK(cn.stats().totalCreated > 400000);   // lots of churn happened
        CHECK(cn.checkInvariants());
    }

    // CELL-connectivity churn retention (review defect 1): mixed-arity
    // create/remove must keep memory bounded by PEAK live cells, NOT by
    // totalCreated. The old append-only arena retained ~40 MB after 1M
    // create/remove of a single cell; inline connectivity retains nothing.
    {
        NodePool npB;   // (were np/nh/cp: C4456/-Wshadow of an earlier block)
        std::vector<NodeHandle> nhB;
        for (int i = 0; i < 32; ++i) nhB.push_back(npB.add({double(i),0,0}, {3,1}));
        CellPool cpB;
        const topo::ElementType tys[3] = {topo::ElementType::Tet4,
                                          topo::ElementType::Tri3,
                                          topo::ElementType::Segment2};
        const int ar[3] = {4, 3, 2};
        std::vector<CellHandle> livec;
        std::mt19937 rng(11);
        const std::size_t target = 4000;
        std::size_t peakMaterialized = 0;
        for (int i = 0; i < 1000000; ++i) {
            if (livec.size() >= target ||
                (!livec.empty() && (rng() & 1))) {
                const std::size_t k = rng() % livec.size();
                cpB.remove(livec[k]);
                livec[k] = livec.back();
                livec.pop_back();
            } else {
                const int ti = i % 3;
                std::vector<NodeHandle> conn;
                for (int j = 0; j < ar[ti]; ++j)
                    conn.push_back(nhB[size_t((i + j) % 32)]);
                livec.push_back(cpB.add(tys[ti], conn, {3,1}, npB));
            }
            peakMaterialized = std::max(peakMaterialized,
                                        cpB.stats().materializedSlots);
        }
        const auto st = cpB.stats();
        // Bounded by peak live set, NOT by 1,000,000 operations.
        CHECK(st.materializedSlots <= target + 1);
        CHECK(peakMaterialized <= target + 1);
        CHECK(st.totalCreated > 400000);
        // payloadBytes proportional to materialized slots, not history.
        CHECK(st.payloadBytes <=
              (target + 1) * sizeof(std::array<topo::CellPool::LocalNodeRef,
                                               10>) + 4096);
        CHECK(cpB.checkInvariants());
    }
}


// ---------------------------------------------------------------------------
// Detection traits for the ref-qualified connectivity() (rvalue deleted).
// A deleted overload is still selected by overload resolution, so we probe
// with a template that SFINAE-checks invocability WITHOUT naming the
// deleted function in a hard-error context.
template <class S, class = void>
struct ConnLvalueOk : std::false_type {};
template <class S>
struct ConnLvalueOk<S, std::void_t<decltype(
    std::declval<S&>().connectivity())>> : std::true_type {};
template <class S, class = void>
struct ConnRvalueOk : std::false_type {};
template <class S>
struct ConnRvalueOk<S, std::void_t<decltype(
    std::declval<S&&>().connectivity())>> : std::true_type {};
static constexpr bool kConnectivityLvalueOk =
    ConnLvalueOk<paacs::topo::CellPool::CellSnapshot>::value;
static constexpr bool kConnectivityRvalueOk =
    ConnRvalueOk<paacs::topo::CellPool::CellSnapshot>::value;

static void testStorageApiAndListSafety() {
    TEST("topo: M4-B.2 value-API safety + complete list-invariant checks");
    using topo::NodePool;
    using topo::CellPool;
    using topo::ElementType;
    using topo::NodeHandle;

    // (A) By-value snapshots survive reallocating insertion (review item 2 -
    // was heap-use-after-free with reference/span returns).
    {
        NodePool np;
        std::vector<NodeHandle> nh;
        for (int i = 0; i < 4; ++i) nh.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto cell = cp.add(ElementType::Tet4, nh, {3,1}, np);
        const auto snap = cp.get(cell);            // by value
        const auto nd = np.get(nh[0]);             // by value
        for (int i = 0; i < 5000; ++i) np.add({double(i),9,9}, {0,0});
        for (int i = 0; i < 5000; ++i) {
            std::vector<NodeHandle> c = {nh[0],nh[1],nh[2],nh[3]};
            cp.add(ElementType::Tet4, c, {3,1}, np);
        }
        CHECK(snap.count == 4);                    // still valid
        CHECK(nd.coord.x == 0.0);
        CHECK(snap.nodes[0] == nh[0]);
    }

    // (A2) rvalue-snapshot dangling span is now a COMPILE error (review
    // defect 2). The unsafe pattern `cp.get(h).connectivity()` no longer
    // compiles (connectivity() const&& = deleted). The safe pattern below
    // keeps the snapshot as a named lvalue; the by-value nodes/count path
    // is always safe on a temporary.
    {
        NodePool np;
        std::vector<NodeHandle> nh;
        for (int i = 0; i < 4; ++i) nh.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto h = cp.add(ElementType::Tet4, nh, {3,1}, np);
        // Safe: named snapshot outlives the span.
        auto snap = cp.get(h);
        auto sp = snap.connectivity();
        CHECK(sp.size() == 4 && sp[0] == nh[0]);
        // Safe on a temporary: by-value array + count (no borrow).
        CHECK(cp.get(h).count == 4);
        CHECK(cp.get(h).nodes[3] == nh[3]);
        // Compile-time proof the rvalue overload is gone. A deleted
        // function still participates in overload resolution, so a plain
        // requires-expression would name it and hard-error; instead we
        // check that connectivity() is invocable on an lvalue but not an
        // rvalue via a detached trait.
        static_assert(kConnectivityLvalueOk,
                      "connectivity() must work on an lvalue snapshot");
        static_assert(!kConnectivityRvalueOk,
                      "connectivity() must be deleted on rvalue snapshots");
    }

    // (B) span overload avoids per-cell heap vector (review carry-item).
    {
        NodePool np;
        std::array<NodeHandle, 4> a{np.add({0,0,0},{3,1}), np.add({1,0,0},{3,1}),
                                    np.add({0,1,0},{3,1}), np.add({0,0,1},{3,1})};
        CellPool cp;
        auto h = cp.add(ElementType::Tet4,
                        std::span<const NodeHandle>(a.data(), a.size()),
                        {3,1}, np);
        CHECK(cp.isCurrent(h) && cp.get(h).count == 4);
    }

    // (C) forEachActive borrow iteration visits all live cells in order.
    {
        NodePool np;
        std::vector<NodeHandle> nh;
        for (int i = 0; i < 4; ++i) nh.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        for (int i = 0; i < 5; ++i)
            cp.add(ElementType::Tet4, nh, {3,1}, np);
        int n = 0;
        std::uint64_t prev = 0;
        bool ok = true;
        cp.forEachActive([&](topo::CellHandle h, ElementType t,
                             std::span<const NodeHandle> conn,
                             mesh::Classification) {
            ++n;
            if (t != ElementType::Tet4 || conn.size() != 4) ok = false;
            const auto s = cp.seqOf(h);
            if (n > 1 && !(s > prev)) ok = false;
            prev = s;
        });
        CHECK(n == 5 && ok);
    }

    // (D) COMPLETE list-invariant checks catch corruptions the old checker
    // accepted (review item 1): bad tail, bad middle prev, bad stored seq.
    {
        NodePool np;
        std::vector<NodeHandle> nh;
        for (int i = 0; i < 5; ++i) nh.push_back(np.add({double(i),0,0}, {0,0}));
        CHECK(np.checkInvariants());
        // Corrupt tail pointer.
        {
            NodePool p2;
            std::vector<NodeHandle> h2;
            for (int i = 0; i < 5; ++i) h2.push_back(p2.add({double(i),0,0}, {0,0}));
            p2.testSetTail(p2.testHead());         // tail wrong
            CHECK(!p2.checkInvariants());
        }
        // Corrupt a middle prev link.
        {
            NodePool p3;
            std::vector<NodeHandle> h3;
            for (int i = 0; i < 5; ++i) h3.push_back(p3.add({double(i),0,0}, {0,0}));
            p3.testSetPrev(h3[2].slot, 99);        // wrong back-link
            CHECK(!p3.checkInvariants());
        }
        // Corrupt a stored sequence number.
        {
            NodePool p4;
            std::vector<NodeHandle> h4;
            for (int i = 0; i < 5; ++i) h4.push_back(p4.add({double(i),0,0}, {0,0}));
            p4.testSetSeq(h4[2], 999999);          // seq != allocator seq
            CHECK(!p4.checkInvariants());
        }

        // Bounds-safety (review blocker): corrupt ENDPOINTS to
        // out-of-range values. The auditor must return false and
        // auditStorage must report StorageInvariantViolation, never
        // dereference out of bounds (ASan-verified in topo.storage group).
        auto nodeAuditFails = [](NodePool&& p) {
            const auto rep = [&]{
                // build a matching empty cellpool for the gate
                return true;
            }();
            (void)rep;
            return !p.checkInvariants();
        };
        auto mkNodes = [](int n) {
            NodePool p;
            for (int i = 0; i < n; ++i) p.add({double(i),0,0}, {0,0});
            return p;
        };
        { auto p = mkNodes(5); p.testSetTail(2147483647);   // INT_MAX
          CHECK(!p.checkInvariants()); }
        { auto p = mkNodes(5); p.testSetHead(p.testMaterialized());
          CHECK(!p.checkInvariants()); }                    // head == materialized
        { auto p = mkNodes(5); p.testSetTail(-7);           // negative != -1
          CHECK(!p.checkInvariants()); }
        { auto p = mkNodes(5); p.testSetNext(p.testHead(), 99999);
          CHECK(!p.checkInvariants()); }                    // OOB middle next
        { auto p = mkNodes(5); p.testSetPrev(p.testHead(), 88888);
          CHECK(!p.checkInvariants()); }                    // OOB middle prev
        // And auditStorage reports the violation without crashing.
        {
            auto p = mkNodes(5); p.testSetTail(2147483647);
            CellPool empty;
            const auto rep = topo::auditStorage(p, empty);
            bool sawViolation = false;
            for (const auto& d : rep.items)
                sawViolation |= d.code ==
                    diag::Code::StorageInvariantViolation;
            CHECK(sawViolation);
        }
        (void)nodeAuditFails;
    }

    // (D2) CellPool malformed endpoints (review blocker: both pools).
    {
        NodePool np;
        std::vector<NodeHandle> nh;
        for (int i = 0; i < 4; ++i) nh.push_back(np.add({double(i),0,0}, {3,1}));
        auto mkCells = [&](int n) {
            CellPool p;
            for (int i = 0; i < n; ++i) p.add(ElementType::Tet4, nh, {3,1}, np);
            return p;
        };
        { auto p = mkCells(5); p.testSetTail(2147483647);
          CHECK(!p.checkInvariants()); }
        { auto p = mkCells(5); p.testSetHead(p.testMaterialized());
          CHECK(!p.checkInvariants()); }
        { auto p = mkCells(5); p.testSetTail(-7);
          CHECK(!p.checkInvariants()); }
        { auto p = mkCells(5); p.testSetNext(p.testHead(), 99999);
          CHECK(!p.checkInvariants()); }
        { auto p = mkCells(5); p.testSetPrev(p.testHead(), 88888);
          CHECK(!p.checkInvariants()); }
        {
            auto p = mkCells(5); p.testSetTail(2147483647);
            const auto rep = topo::auditStorage(np, p);
            bool sawViolation = false;
            for (const auto& d : rep.items)
                sawViolation |= d.code ==
                    diag::Code::StorageInvariantViolation;
            CHECK(sawViolation);
        }
    }

    // (E) Memory budget: guards against order-of-magnitude regression
    // (e.g. arena retention), NOT against normal geometric vector slack.
    // Steady-state ~139 B/cell; ceiling deliberately portable-loose
    // (review carry-item: the old <200 assertion had almost no margin at a
    // high-water point). Actual Linux measurement recorded in the round-7
    // note.
    {
        NodePool np;
        std::vector<NodeHandle> nh;
        for (int i = 0; i < 5000; ++i) nh.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        for (int i = 0; i < 5000; ++i) {
            std::vector<NodeHandle> c = {nh[i], nh[(i+1)%5000],
                                        nh[(i+2)%5000], nh[(i+3)%5000]};
            cp.add(ElementType::Tet4, c, {3,1}, np);
        }
        const double perCell = double(cp.stats().totalBytes) /
                               double(cp.stats().materializedSlots);
        CHECK(perCell < 320.0);                    // ~140 observed steady-state
    }
}


// ---------------------------------------------------------------------------
static void testEdgeTopology() {
    TEST("topo: M4-B.3 explicit reusable edge topology");
    using topo::NodePool;
    using topo::CellPool;
    using topo::EdgePool;
    using topo::ElementType;
    using topo::NodeHandle;
    using topo::CellHandle;

    // (A) Shared-edge reuse across cells: two tets sharing a face share 3
    // edges, so total unique = 6 + 6 - 3 = 9.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 5; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c1 = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        auto c2 = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[4]}, {3,1}, np);
        EdgePool ep;
        auto e1 = ep.insertCellEdges(np, cp, c1);
        auto e2 = ep.insertCellEdges(np, cp, c2);
        CHECK(e1.size() == 6 && e2.size() == 6);
        CHECK(ep.liveCount() == 9);
        CHECK(ep.checkInvariants(np));
        // Shared edge (0,1) has use count 2.
        auto sh = ep.find(np, n[0], n[1]);
        CHECK(sh.has_value());
        CHECK(ep.testUseCount(*sh) == 2);
        // A non-shared edge (0,3) has use count 1.
        auto ns = ep.find(np, n[0], n[3]);
        CHECK(ns.has_value() && ep.testUseCount(*ns) == 1);
    }

    // (B) Use-count maintenance on removal: removing one cell decrements
    // shared edges and retires cell-exclusive ones (edge exists iff used).
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 5; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c1 = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        auto c2 = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[4]}, {3,1}, np);
        EdgePool ep;
        ep.insertCellEdges(np, cp, c1);
        ep.insertCellEdges(np, cp, c2);
        CHECK(ep.liveCount() == 9);
        ep.removeCellEdges(cp, c2);
        CHECK(ep.liveCount() == 6);            // c2-exclusive edges retired
        CHECK(ep.checkInvariants(np));
        auto sh = ep.find(np, n[0], n[1]);
        CHECK(sh.has_value() && ep.testUseCount(*sh) == 1);
        // A c2-exclusive edge (0,4) is gone.
        CHECK(!ep.find(np, n[0], n[4]).has_value());
        ep.removeCellEdges(cp, c1);
        CHECK(ep.liveCount() == 0);            // all edges retired
        CHECK(ep.checkInvariants(np));
    }

    // (C) Reuse of retired slots: re-inserting after full removal reuses
    // slots (three-state accounting), not unbounded growth.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        EdgePool ep;
        for (int rep = 0; rep < 1000; ++rep) {
            ep.insertCellEdges(np, cp, c);
            CHECK(ep.liveCount() == 6);
            ep.removeCellEdges(cp, c);
            CHECK(ep.liveCount() == 0);
        }
        // Materialized bounded by peak (6), not by 1000*6 insertions.
        CHECK(ep.stats().materializedSlots <= 6);
        CHECK(ep.checkInvariants(np));
        const auto statsBefore = ep.stats();
        EdgePool moved(std::move(ep));
        CHECK(moved.stats().materializedSlots == statsBefore.materializedSlots);
        CHECK(moved.stats().reusableFree == statsBefore.reusableFree);
        CHECK(moved.stats().totalCreated == statsBefore.totalCreated);
        CHECK(ep.poolId() == 0 && ep.stats().materializedSlots == 0);
        const auto reused = moved.insertCellEdges(np, cp, c);
        CHECK(reused.size() == 6 && moved.stats().materializedSlots <= 6);
        CHECK(moved.checkInvariants(np, cp));
    }

    // (D) Self-edge rejected: find with identical endpoints is nullopt; a
    // degenerate pair never creates a record.
    {
        NodePool np;
        auto a = np.add({0,0,0}, {3,1});
        EdgePool ep;
        CHECK(!ep.find(np, a, a).has_value());
    }

    // (E) Duplicate insertion of the same cell's edges is rejected (keeps
    // use counts honest).
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        EdgePool ep;
        ep.insertCellEdges(np, cp, c);
        bool threw = false;
        try { ep.insertCellEdges(np, cp, c); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
        CHECK(ep.liveCount() == 6);            // unchanged
    }

    // (F) High-valence degeneracy: many cells sharing one edge. Build a fan
    // of tets all sharing edge (0,1); that edge's use count == fan size.
    {
        NodePool np;
        auto p0 = np.add({0,0,0}, {3,1});
        auto p1 = np.add({1,0,0}, {3,1});
        std::vector<NodeHandle> ring;
        const int F = 40;
        for (int i = 0; i < F + 1; ++i)
            ring.push_back(np.add({0.5, std::cos(0.1*i), std::sin(0.1*i)}, {3,1}));
        CellPool cp;
        EdgePool ep;
        for (int i = 0; i < F; ++i) {
            auto c = cp.add(ElementType::Tet4,
                            {p0, p1, ring[i], ring[i+1]}, {3,1}, np);
            ep.insertCellEdges(np, cp, c);
        }
        auto shared = ep.find(np, p0, p1);
        CHECK(shared.has_value());
        CHECK(ep.testUseCount(*shared) == std::uint32_t(F));  // high valence
        CHECK(ep.checkInvariants(np));
    }

    // (G) Constraint flags and classification survive edits to other edges.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 5; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c1 = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        EdgePool ep;
        ep.insertCellEdges(np, cp, c1);
        auto e01 = ep.find(np, n[0], n[1]);
        CHECK(e01.has_value());
        ep.setRequiredConstraint(*e01, true);
        ep.setClassification(*e01, {1, 42});
        // Add another cell (does not touch edge 0-1's flags).
        auto c2 = cp.add(ElementType::Tet4, {n[0],n[2],n[3],n[4]}, {3,1}, np);
        ep.insertCellEdges(np, cp, c2);
        auto again = ep.find(np, n[0], n[1]);
        CHECK(again.has_value());
        CHECK(ep.get(*again).requiredConstraint == true);
        CHECK(ep.get(*again).cls.dim == 1 && ep.get(*again).cls.tag == 42);
        CHECK(ep.checkInvariants(np));
    }

    // (H) Stale cell handle rejected.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        EdgePool ep;
        ep.insertCellEdges(np, cp, c);
        ep.removeCellEdges(cp, c);
        cp.remove(c);
        bool threw = false;
        try { ep.insertCellEdges(np, cp, c); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }

    // (I) Mixed element types extract the right edge counts (traits-driven).
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 8; ++i) n.push_back(np.add({double(i),double(i%3),0}, {3,1}));
        CellPool cp;
        EdgePool ep;
        // Hex8 has 12 edges; a single hex => 12 unique edges.
        auto hx = cp.add(ElementType::Hex8,
                         {n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]}, {3,1}, np);
        ep.insertCellEdges(np, cp, hx);
        CHECK(ep.liveCount() == 12);
        CHECK(ep.checkInvariants(np));
    }
}


// ---------------------------------------------------------------------------
static void testEdgeLifecycle() {
    TEST("topo: M4-B.3 edge lifecycle integrity (ghost detection, rollback)");
    using topo::NodePool;
    using topo::CellPool;
    using topo::EdgePool;
    using topo::ElementType;
    using topo::NodeHandle;
    using topo::CellHandle;

    // (A) Ghost edge: cell removed WITHOUT removeCellEdges. The cell-aware
    // audit must catch it even though the nodes-only record check passes.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        EdgePool ep;
        ep.insertCellEdges(np, cp, c);
        CHECK(ep.checkInvariants(np));              // record-level: ok
        CHECK(ep.checkInvariants(np, cp));          // cell-aware: ok
        cp.remove(c);                               // uncoordinated removal
        CHECK(ep.checkInvariants(np));              // record-level STILL ok
        CHECK(!ep.checkInvariants(np, cp));         // cell-aware CATCHES ghost
        const auto rep = topo::auditEdgeTopology(np, cp, ep);
        bool ghost = false;
        for (const auto& d : rep.items)
            ghost |= d.code == diag::Code::GhostEdge;
        CHECK(ghost);
        CHECK(!topo::edgeTopologyReady(np, cp, ep));
    }

    // (B) Reused cell slot with stale incidence: after the ghost above,
    // reusing the slot for a new cell and inserting its edges must be
    // REJECTED (stale incidence), not silently overwrite.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 5; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c1 = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        EdgePool ep;
        ep.insertCellEdges(np, cp, c1);
        cp.remove(c1);                              // ghost incidence at slot
        auto c2 = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[4]}, {3,1}, np);
        // c2 likely reuses c1's slot; inserting must reject stale incidence.
        if (c2.slot == c1.slot) {
            bool threw = false;
            try { ep.insertCellEdges(np, cp, c2); }
            catch (const std::invalid_argument&) { threw = true; }
            CHECK(threw);
        }
    }

    // (C) Node removed before edge removal: cell-aware audit catches the
    // dangling endpoint via storageEdgeReady short-circuit.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        EdgePool ep;
        ep.insertCellEdges(np, cp, c);
        np.remove(n[0]);                            // strands cell + edges
        CHECK(!topo::edgeTopologyReady(np, cp, ep));
        const auto rep = topo::auditEdgeTopology(np, cp, ep);
        CHECK(!rep.items.empty());
    }

    // (D) Allocation failure during insertion => strong exception safety:
    // pool unchanged, no partial edges, audit still clean.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        EdgePool ep;
        // Inject a fault at successive commit steps and verify rollback.
        for (int stage = 1; stage <= 6; ++stage) {
            topo::faultCountdown() = stage;
            bool threw = false;
            try { ep.insertCellEdges(np, cp, c); }
            catch (const std::bad_alloc&) { threw = true; }
            topo::faultCountdown() = 0;
            if (!threw) {
                // Insertion succeeded (fault landed past the mutation);
                // undo so the next stage starts clean.
                ep.removeCellEdges(cp, c);
                continue;
            }
            // Rolled back: no live edges, audit clean.
            CHECK(ep.liveCount() == 0);
            CHECK(ep.checkInvariants(np, cp));
        }
    }

    // (E) Partial-insertion rollback leaves NO orphan records and the edge
    // table empty.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        EdgePool ep;
        topo::faultCountdown() = 3;                 // fail mid-commit
        bool threw = false;
        try { ep.insertCellEdges(np, cp, c); }
        catch (const std::bad_alloc&) { threw = true; }
        topo::faultCountdown() = 0;
        CHECK(threw);
        CHECK(ep.liveCount() == 0);
        CHECK(ep.stats().materializedSlots == ep.stats().reusableFree +
              ep.stats().retired);                  // no live orphan
        CHECK(ep.checkInvariants(np, cp));
        // A subsequent clean insert works normally.
        auto e = ep.insertCellEdges(np, cp, c);
        CHECK(e.size() == 6 && ep.liveCount() == 6);
        CHECK(ep.checkInvariants(np, cp));
    }

    // (F) removeCellEdges validates against the CellPool: foreign pool and
    // never-inserted cells are rejected.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp, other;
        auto c = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        EdgePool ep;
        ep.insertCellEdges(np, cp, c);
        bool foreignThrew = false;
        try { ep.removeCellEdges(other, c); }
        catch (const std::invalid_argument&) { foreignThrew = true; }
        CHECK(foreignThrew);
        ep.removeCellEdges(cp, c);                  // correct pool: ok
        CHECK(ep.liveCount() == 0);
    }

    // (G) Coordinated mixed-cell churn keeps the cell-aware audit clean
    // throughout (edges added and removed with cells, in order).
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 16; ++i) n.push_back(np.add({double(i),double(i%3),0}, {3,1}));
        CellPool cp;
        EdgePool ep;
        std::mt19937 rng(5);
        std::vector<CellHandle> live;
        const ElementType tys[2] = {ElementType::Tet4, ElementType::Tri3};
        const int ar[2] = {4, 3};
        for (int i = 0; i < 4000; ++i) {
            if (!live.empty() && (rng() & 1)) {
                const std::size_t k = rng() % live.size();
                ep.removeCellEdges(cp, live[k]);   // edges first
                cp.remove(live[k]);                 // then cell
                live[k] = live.back();
                live.pop_back();
            } else {
                const int t = i % 2;
                std::vector<NodeHandle> conn;
                for (int j = 0; j < ar[t]; ++j)
                    conn.push_back(n[std::size_t((i + j) % 16)]);
                // skip degenerate (repeated node) combos
                std::sort(conn.begin(), conn.end());
                if (std::adjacent_find(conn.begin(), conn.end()) != conn.end())
                    continue;
                auto c = cp.add(tys[t], conn, {3,1}, np);
                ep.insertCellEdges(np, cp, c);
                live.push_back(c);
            }
            if ((i & 511) == 0) CHECK(ep.checkInvariants(np, cp));
        }
        CHECK(ep.checkInvariants(np, cp));
        CHECK(topo::edgeTopologyReady(np, cp, ep));
    }
}


// ---------------------------------------------------------------------------
static void testEdgeHardening() {
    TEST("topo: M4-B.3 edge hardening (stale-cell, foreign-pool, corruption)");
    using topo::NodePool;
    using topo::CellPool;
    using topo::EdgePool;
    using topo::ElementType;
    using topo::NodeHandle;
    using topo::CellHandle;
    using topo::EdgeHandle;

    // (A) removeCellEdges rejects a stale (already-removed) cell.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        EdgePool ep;
        ep.insertCellEdges(np, cp, c);
        cp.remove(c);
        bool threw = false;
        try { ep.removeCellEdges(cp, c); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
        CHECK(ep.liveCount() == 6);             // edges NOT silently dropped
    }

    // (B) find() rejects a foreign NodePool with matching sequences.
    {
        NodePool A;
        std::vector<NodeHandle> na;
        for (int i = 0; i < 4; ++i) na.push_back(A.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4, {na[0],na[1],na[2],na[3]}, {3,1}, A);
        EdgePool ep;
        ep.insertCellEdges(A, cp, c);
        NodePool B;                              // same sequences, diff pool
        std::vector<NodeHandle> nb;
        for (int i = 0; i < 4; ++i) nb.push_back(B.add({double(i)*7,0,0}, {3,1}));
        CHECK(!ep.find(B, nb[0], nb[1]).has_value());   // no foreign leak
        CHECK(ep.find(A, na[0], na[1]).has_value());    // own pool still ok
    }

    // (C) Real allocation-failure rollback via the fault seam at every
    // commit step: pool fully consistent, no partial edge, on each throw.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 8; ++i) n.push_back(np.add({double(i),double(i%3),0}, {3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Hex8,
                        {n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]}, {3,1}, np);
        for (int stage = 1; stage <= 30; ++stage) {
            EdgePool ep;
            topo::faultCountdown() = stage;
            bool threw = false;
            try { ep.insertCellEdges(np, cp, c); }
            catch (const std::bad_alloc&) { threw = true; }
            topo::faultCountdown() = 0;
            if (threw) {
                CHECK(ep.liveCount() == 0);
                CHECK(ep.checkInvariants(np, cp));
                CHECK(ep.stats().materializedSlots ==
                      ep.stats().reusableFree + ep.stats().retired);
            } else {
                CHECK(ep.liveCount() == 12);   // Hex8 => 12 edges committed
                CHECK(ep.checkInvariants(np, cp));
            }
        }
    }

    // (D) Edge-level corruption probes: malformed endpoints/links and table
    // corruption are caught by checkInvariants without OOB access.
    {
        auto mk = [&](NodePool& np, CellPool& cp, EdgePool& ep) {
            std::vector<NodeHandle> n;
            for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
            auto c = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
            ep.insertCellEdges(np, cp, c);
        };
        { NodePool np; CellPool cp; EdgePool ep; mk(np,cp,ep);
          ep.testSetTail(2147483647); CHECK(!ep.checkInvariants(np)); }
        { NodePool np; CellPool cp; EdgePool ep; mk(np,cp,ep);
          ep.testSetHead(ep.testMaterialized()); CHECK(!ep.checkInvariants(np)); }
        { NodePool np; CellPool cp; EdgePool ep; mk(np,cp,ep);
          ep.testSetTail(-7); CHECK(!ep.checkInvariants(np)); }
        { NodePool np; CellPool cp; EdgePool ep; mk(np,cp,ep);
          ep.testSetNext(ep.testHead(), 99999); CHECK(!ep.checkInvariants(np)); }
        { NodePool np; CellPool cp; EdgePool ep; mk(np,cp,ep);
          ep.testSetPrev(ep.testHead(), 88888); CHECK(!ep.checkInvariants(np)); }
        // Table entry pointing to the wrong handle.
        { NodePool np; CellPool cp; EdgePool ep; mk(np,cp,ep);
          auto live = ep.activeBySequence();
          ep.testCorruptTableEntry(live[0], live[1]);
          CHECK(!ep.checkInvariants(np)); }
        // Missing table entry for a live edge.
        { NodePool np; CellPool cp; EdgePool ep; mk(np,cp,ep);
          auto live = ep.activeBySequence();
          ep.testDropTableEntry(live[0]);
          CHECK(!ep.checkInvariants(np)); }
        // Use count disagreeing with incidence (needs the cells form).
        { NodePool np; CellPool cp; EdgePool ep;
          std::vector<NodeHandle> n;
          for (int i = 0; i < 4; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
          auto c = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
          ep.insertCellEdges(np, cp, c);
          auto live = ep.activeBySequence();
          ep.testSetUseCount(live[0], 7);
          CHECK(!ep.checkInvariants(np, cp)); }
    }

    // (E) Ordinary REMOVAL under allocation pressure never terminates and
    // never leaves inconsistent state: removeCellEdges pre-reserves
    // free-list capacity before mutation (review blocker). Uses the
    // maybeFault seam to fail the reservation step; on throw, logical state
    // is preserved; on success, removal completes with no allocation in the
    // retirement loop.
    {
        for (int stage = 1; stage <= 8; ++stage) {
            NodePool np;
            std::vector<NodeHandle> n;
            for (int i = 0; i < 8; ++i)
                n.push_back(np.add({double(i),double(i%3),0}, {3,1}));
            CellPool cp;
            auto c1 = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
            auto c2 = cp.add(ElementType::Tet4, {n[4],n[5],n[6],n[7]}, {3,1}, np);
            EdgePool ep;
            ep.insertCellEdges(np, cp, c1);
            ep.insertCellEdges(np, cp, c2);
            const auto liveBefore = ep.liveCount();
            topo::faultCountdown() = stage;
            bool threw = false;
            try { ep.removeCellEdges(cp, c2); }
            catch (const std::bad_alloc&) { threw = true; }
            topo::faultCountdown() = 0;
            // Consistent either way.
            CHECK(ep.checkInvariants(np, cp));
            if (threw) {
                CHECK(ep.liveCount() == liveBefore);   // logical state kept
                // A clean retry succeeds.
                ep.removeCellEdges(cp, c2);
            }
            CHECK(ep.checkInvariants(np, cp));
        }
    }

    static_assert(!std::is_nothrow_move_constructible_v<EdgePool>);
    static_assert(std::is_nothrow_move_constructible_v<
                  topo::HandleAllocator<topo::EdgeTag>>);

    // (F0) Empty public move is tokenless: it preserves the source identity
    // in the target and does not consume the next EdgeTag pool token.
    {
        EdgePool src;
        const auto sourcePool = src.poolId();
        EdgePool dst(std::move(src));
        CHECK(dst.poolId() == sourcePool && dst.liveCount() == 0);
        CHECK(src.poolId() == 0 && src.liveCount() == 0);
        CHECK(src.stats().materializedSlots == 0 && src.stats().consistent());
        EdgePool nextOrdinary;
        CHECK(nextOrdinary.poolId() == sourcePool + 1);
    }

    // (F) EdgePool public move construction qualification.
    {
        NodePool np;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 5; ++i) n.push_back(np.add({double(i),0,0}, {3,1}));
        CellPool cp;
        auto c1 = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[3]}, {3,1}, np);
        auto c2 = cp.add(ElementType::Tet4, {n[0],n[1],n[2],n[4]}, {3,1}, np);
        EdgePool src;
        src.insertCellEdges(np, cp, c1);
        src.insertCellEdges(np, cp, c2);
        auto shared = src.find(np, n[0], n[1]);
        CHECK(shared.has_value());
        const auto hBefore = *shared;
        src.setRequiredConstraint(hBefore, true);
        src.setClassification(hBefore, {1, 91});
        const auto poolBefore = src.poolId();
        const auto orderBefore = src.activeBySequence();
        const std::vector<topo::EdgeHandle> c1Before(
            src.cellEdges(cp, c1).begin(), src.cellEdges(cp, c1).end());
        const std::vector<topo::EdgeHandle> c2Before(
            src.cellEdges(cp, c2).begin(), src.cellEdges(cp, c2).end());
        std::vector<std::uint64_t> seqBefore;
        for (const auto& h : orderBefore) seqBefore.push_back(src.seqOf(h));
        const auto statsBefore = src.stats();
        EdgePool dst(std::move(src));
        // Existing handles current in the moved-to pool.
        CHECK(dst.poolId() == poolBefore);
        CHECK(dst.isCurrent(hBefore));
        CHECK(dst.find(np, n[0], n[1]).has_value());
        CHECK(*dst.find(np, n[0], n[1]) == hBefore);
        // Table + incidence agree; audit passes.
        CHECK(dst.checkInvariants(np, cp));
        CHECK(topo::edgeTopologyReady(np, cp, dst));
        CHECK(dst.liveCount() == 9);
        CHECK(dst.activeBySequence() == orderBefore);
        CHECK(std::vector<topo::EdgeHandle>(dst.cellEdges(cp, c1).begin(),
                  dst.cellEdges(cp, c1).end()) == c1Before);
        CHECK(std::vector<topo::EdgeHandle>(dst.cellEdges(cp, c2).begin(),
                  dst.cellEdges(cp, c2).end()) == c2Before);
        for (std::size_t i = 0; i < orderBefore.size(); ++i)
            CHECK(dst.seqOf(orderBefore[i]) == seqBefore[i]);
        CHECK(dst.stats().materializedSlots == statsBefore.materializedSlots);
        CHECK(dst.stats().live == statsBefore.live);
        CHECK(dst.stats().reusableFree == statsBefore.reusableFree);
        CHECK(dst.stats().retired == statsBefore.retired);
        CHECK(dst.stats().totalCreated == statsBefore.totalCreated);
        CHECK(dst.get(hBefore).requiredConstraint);
        CHECK(dst.get(hBefore).cls.dim == 1 && dst.get(hBefore).cls.tag == 91);
        // Moved-from is empty/unbound/unusable.
        CHECK(src.poolId() == 0 && src.liveCount() == 0);
        CHECK(src.boundNodePool() == 0 && src.boundCellPool() == 0);
        CHECK(src.activeBySequence().empty());
        CHECK(src.stats().materializedSlots == 0 && src.stats().live == 0 &&
              src.stats().reusableFree == 0 && src.stats().retired == 0 &&
              src.stats().totalCreated == 0 && src.stats().consistent());
        bool unusable = false;
        try { (void)src.insertCellEdges(np, cp, c1); }
        catch (const std::logic_error&) { unusable = true; }
        CHECK(unusable);
        // Removal still works after move.
        dst.removeCellEdges(cp, c2);
        CHECK(dst.liveCount() == 6);
        CHECK(dst.checkInvariants(np, cp));
    }
}


// ---------------------------------------------------------------------------
static void testFaceTopology() {
    TEST("topo: M4-B.4 explicit face topology");
    using topo::NodePool; using topo::CellPool; using topo::EdgePool;
    using topo::FacePool; using topo::ElementType; using topo::NodeHandle;
    using topo::CellHandle; using topo::FaceHandle; using topo::FaceClass;

    auto edgesReady = [](NodePool& np, CellPool& cp, EdgePool& ep,
                         const std::vector<CellHandle>& cells) {
        for (auto c : cells) ep.insertCellEdges(np, cp, c);
    };

    // (A) Face sharing + boundary/interior: two tets sharing a face =>
    // 7 unique faces, 1 interior, 6 boundary.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<5;++i) n.push_back(np.add({double(i),double(i%2),0},{3,1}));
        CellPool cp;
        auto c1 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        auto c2 = cp.add(ElementType::Tet4,{n[0],n[2],n[1],n[4]},{3,1},np);
        EdgePool ep; edgesReady(np,cp,ep,{c1,c2});
        FacePool fp;
        auto f1 = fp.insertCellFaces(np,cp,ep,c1);
        auto f2 = fp.insertCellFaces(np,cp,ep,c2);
        CHECK(f1.size()==4 && f2.size()==4);
        CHECK(fp.liveCount()==7);
        CHECK(fp.boundaryFaces().size()==6);
        CHECK(fp.checkInvariants(np,cp));
        CHECK(topo::faceTopologyReady(np,cp,ep,fp));
        std::array<NodeHandle,3> sh{n[0],n[1],n[2]};
        auto shared = fp.find(np, std::span<const NodeHandle>(sh.data(),3));
        CHECK(shared.has_value());
        CHECK(fp.classify(*shared)==FaceClass::Interior);
    }

    // (B) Orientation sense: a valid interior face has opposite senses.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<5;++i) n.push_back(np.add({double(i),0,0},{3,1}));
        CellPool cp;
        auto cA = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        auto cB = cp.add(ElementType::Tet4,{n[0],n[2],n[1],n[4]},{3,1},np);
        EdgePool ep; edgesReady(np,cp,ep,{cA,cB});
        FacePool fp;
        fp.insertCellFaces(np,cp,ep,cA);
        fp.insertCellFaces(np,cp,ep,cB);
        std::array<NodeHandle,3> sh{n[0],n[1],n[2]};
        auto shared = fp.find(np, std::span<const NodeHandle>(sh.data(),3));
        CHECK(shared.has_value());
        auto fd = fp.get(*shared);
        CHECK(fd.incident.size()==2);
        CHECK(fd.incident[0].sense + fd.incident[1].sense == 0);  // opposite
    }

    // (C) Nonmanifold detection: three tets sharing one face.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<6;++i) n.push_back(np.add({double(i),double(i%3),0},{3,1}));
        CellPool cp;
        auto c1 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        auto c2 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[4]},{3,1},np);
        auto c3 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[5]},{3,1},np);
        EdgePool ep; edgesReady(np,cp,ep,{c1,c2,c3});
        FacePool fp;
        fp.insertCellFaces(np,cp,ep,c1);
        fp.insertCellFaces(np,cp,ep,c2);
        fp.insertCellFaces(np,cp,ep,c3);
        std::array<NodeHandle,3> sh{n[0],n[1],n[2]};
        auto shared = fp.find(np, std::span<const NodeHandle>(sh.data(),3));
        CHECK(shared.has_value());
        CHECK(fp.classify(*shared)==FaceClass::NonManifold);
        CHECK(fp.testIncidenceCount(*shared)==3);
        const auto sharedBefore = fp.get(*shared);
        const auto orderBefore = fp.activeBySequence();
        FacePool moved(std::move(fp));
        CHECK(moved.activeBySequence() == orderBefore);
        CHECK(moved.classify(*shared)==FaceClass::NonManifold);
        const auto sharedAfter = moved.get(*shared);
        CHECK(sharedAfter.incident.size() == sharedBefore.incident.size());
        for (std::size_t i = 0; i < sharedBefore.incident.size(); ++i) {
            CHECK(sharedAfter.incident[i].cellSlot ==
                  sharedBefore.incident[i].cellSlot);
            CHECK(sharedAfter.incident[i].cellGen ==
                  sharedBefore.incident[i].cellGen);
            CHECK(sharedAfter.incident[i].sense ==
                  sharedBefore.incident[i].sense);
        }
        CHECK(moved.checkInvariants(np,cp));
        CHECK(fp.poolId() == 0 && fp.liveCount() == 0);
    }

    // (D) Removal: removing a cell drops its incidence; shared face reverts
    // interior->boundary; cell-exclusive faces retire.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<5;++i) n.push_back(np.add({double(i),double(i%2),0},{3,1}));
        CellPool cp;
        auto c1 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        auto c2 = cp.add(ElementType::Tet4,{n[0],n[2],n[1],n[4]},{3,1},np);
        EdgePool ep; edgesReady(np,cp,ep,{c1,c2});
        FacePool fp;
        fp.insertCellFaces(np,cp,ep,c1);
        fp.insertCellFaces(np,cp,ep,c2);
        CHECK(fp.liveCount()==7);
        fp.removeCellFaces(cp,c2);
        CHECK(fp.liveCount()==4);            // c2-exclusive 3 faces retired
        std::array<NodeHandle,3> sh{n[0],n[1],n[2]};
        auto shared = fp.find(np, std::span<const NodeHandle>(sh.data(),3));
        CHECK(shared.has_value());
        CHECK(fp.classify(*shared)==FaceClass::Boundary);  // reverted
        CHECK(fp.checkInvariants(np,cp));
        fp.removeCellFaces(cp,c1);
        CHECK(fp.liveCount()==0);
        CHECK(fp.checkInvariants(np,cp));
        const auto statsBefore = fp.stats();
        FacePool moved(std::move(fp));
        CHECK(moved.stats().materializedSlots == statsBefore.materializedSlots);
        CHECK(moved.stats().reusableFree == statsBefore.reusableFree);
        CHECK(moved.stats().totalCreated == statsBefore.totalCreated);
        CHECK(fp.poolId() == 0 && fp.stats().materializedSlots == 0);
        const auto reused = moved.insertCellFaces(np,cp,ep,c1);
        CHECK(reused.size()==4);
        CHECK(moved.stats().materializedSlots == statsBefore.materializedSlots);
        CHECK(moved.checkInvariants(np,cp));
    }

    // (E) Mixed element: Hex8 has 6 quad faces.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<8;++i) n.push_back(np.add({double(i),double(i%4),double(i%2)},{3,1}));
        CellPool cp;
        auto h = cp.add(ElementType::Hex8,
                        {n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]},{3,1},np);
        EdgePool ep; edgesReady(np,cp,ep,{h});
        FacePool fp;
        auto ff = fp.insertCellFaces(np,cp,ep,h);
        CHECK(ff.size()==6);
        CHECK(fp.liveCount()==6);
        CHECK(fp.boundaryFaces().size()==6);
        CHECK(fp.checkInvariants(np,cp));
    }

    // (F) Rejections: no-face element (Tri3), stale cell, duplicate, foreign.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<4;++i) n.push_back(np.add({double(i),0,0},{3,1}));
        CellPool cp;
        auto tri = cp.add(ElementType::Tri3,{n[0],n[1],n[2]},{2,1},np);
        EdgePool ep; ep.insertCellEdges(np,cp,tri);
        FacePool fp;
        bool noFace=false;
        try { fp.insertCellFaces(np,cp,ep,tri); }
        catch(const std::invalid_argument&){ noFace=true; }
        CHECK(noFace);                        // 2D element has no faces
        auto tet = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        ep.insertCellEdges(np,cp,tet);
        fp.insertCellFaces(np,cp,ep,tet);
        bool dup=false;
        try { fp.insertCellFaces(np,cp,ep,tet); }
        catch(const std::invalid_argument&){ dup=true; }
        CHECK(dup);
        fp.removeCellFaces(cp,tet);
        cp.remove(tet);
        bool stale=false;
        try { fp.insertCellFaces(np,cp,ep,tet); }
        catch(const std::invalid_argument&){ stale=true; }
        CHECK(stale);
    }

    // (G) Transactional insert: fault injection => full rollback, clean audit.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<8;++i) n.push_back(np.add({double(i),double(i%4),double(i%2)},{3,1}));
        CellPool cp;
        auto h = cp.add(ElementType::Hex8,
                        {n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]},{3,1},np);
        EdgePool ep; ep.insertCellEdges(np,cp,h);
        for (int stage=1; stage<=14; ++stage) {
            FacePool fp;
            topo::faultCountdown()=stage;
            bool threw=false;
            try { fp.insertCellFaces(np,cp,ep,h); }
            catch(const std::bad_alloc&){ threw=true; }
            topo::faultCountdown()=0;
            if (threw) {
                CHECK(fp.liveCount()==0);
                CHECK(fp.checkInvariants(np,cp));
            } else {
                CHECK(fp.liveCount()==6);
                CHECK(fp.checkInvariants(np,cp));
            }
        }
    }

    // (H) Ghost face: cell removed without removeCellFaces => cell-aware
    // audit catches it.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<4;++i) n.push_back(np.add({double(i),0,0},{3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        EdgePool ep; ep.insertCellEdges(np,cp,c);
        FacePool fp; fp.insertCellFaces(np,cp,ep,c);
        CHECK(fp.checkInvariants(np,cp));
        cp.remove(c);
        CHECK(fp.checkInvariants(np));        // record-only: still ok
        CHECK(!fp.checkInvariants(np,cp));    // cell-aware: catches ghost
        // The top-level audit short-circuits on the (now also ghosted) edge
        // layer first; the face-level cell-aware check is what proves the
        // ghost face directly. Confirm the face pool itself is not ready.
        CHECK(!topo::faceTopologyReady(np,cp,ep,fp));
    }

    // (I) Isolated ghost face reaching the GhostFace diagnostic: keep edges
    // valid (remove the cell's edges too, in order) but DON'T remove faces,
    // then delete the cell -> only the face layer is ghosted.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<4;++i) n.push_back(np.add({double(i),0,0},{3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        EdgePool ep; ep.insertCellEdges(np,cp,c);
        FacePool fp; fp.insertCellFaces(np,cp,ep,c);
        ep.removeCellEdges(cp,c);             // edges cleaned up properly
        cp.remove(c);                         // now only faces are ghosted
        CHECK(topo::edgeTopologyReady(np,cp,ep));   // edges fine
        const auto rep = topo::auditFaceTopology(np,cp,ep,fp);
        bool ghost=false;
        for (const auto& d : rep.items) ghost |= d.code==diag::Code::GhostFace;
        CHECK(ghost);
    }
}


// ---------------------------------------------------------------------------
static void testFaceHardening() {
    TEST("topo: M4-B.4 face hardening (gate, orientation, corruption, move)");
    using topo::NodePool; using topo::CellPool; using topo::EdgePool;
    using topo::FacePool; using topo::ElementType; using topo::NodeHandle;
    using topo::CellHandle; using topo::FaceHandle; using topo::FaceClass;

    // (A) Gate: building faces without the cell's edges present is rejected.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<4;++i) n.push_back(np.add({double(i),0,0},{3,1}));
        CellPool cp;
        auto c = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        EdgePool ep;                          // no edges inserted
        FacePool fp;
        bool threw=false;
        try { fp.insertCellFaces(np,cp,ep,c); }
        catch(const std::invalid_argument&){ threw=true; }
        CHECK(threw);
        CHECK(fp.liveCount()==0);
        // With edges present it works.
        ep.insertCellEdges(np,cp,c);
        auto ff = fp.insertCellFaces(np,cp,ep,c);
        CHECK(ff.size()==4 && fp.checkInvariants(np,cp));
    }

    // (B) Same-sense interior face FAILS the cell-aware audit.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<5;++i) n.push_back(np.add({double(i),double(i%2),0},{3,1}));
        CellPool cp;
        // Both list shared face {0,1,2} in the SAME winding.
        auto c1 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        auto c2 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[4]},{3,1},np);
        EdgePool ep; ep.insertCellEdges(np,cp,c1); ep.insertCellEdges(np,cp,c2);
        FacePool fp;
        fp.insertCellFaces(np,cp,ep,c1);
        fp.insertCellFaces(np,cp,ep,c2);
        std::array<NodeHandle,3> sh{n[0],n[1],n[2]};
        auto shared = fp.find(np, std::span<const NodeHandle>(sh.data(),3));
        CHECK(shared.has_value());
        auto fd = fp.get(*shared);
        CHECK(fd.incident[0].sense + fd.incident[1].sense != 0);  // same sense
        CHECK(!fp.checkInvariants(np,cp));    // audit rejects it
    }

    // (C) Opposite-sense interior face PASSES.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<5;++i) n.push_back(np.add({double(i),0,0},{3,1}));
        CellPool cp;
        auto cA = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        auto cB = cp.add(ElementType::Tet4,{n[0],n[2],n[1],n[4]},{3,1},np);
        EdgePool ep; ep.insertCellEdges(np,cp,cA); ep.insertCellEdges(np,cp,cB);
        FacePool fp;
        fp.insertCellFaces(np,cp,ep,cA);
        fp.insertCellFaces(np,cp,ep,cB);
        CHECK(fp.checkInvariants(np,cp));
    }

    // (D) Corruption probes: bad endpoints/links, table entry corruption,
    // dropped table entry, flipped sense.
    {
        auto mk = [&](NodePool& np, CellPool& cp, EdgePool& ep, FacePool& fp,
                      std::vector<NodeHandle>& n) {
            for (int i=0;i<4;++i) n.push_back(np.add({double(i),0,0},{3,1}));
            auto c = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
            ep.insertCellEdges(np,cp,c);
            fp.insertCellFaces(np,cp,ep,c);
        };
        { NodePool np; CellPool cp; EdgePool ep; FacePool fp; std::vector<NodeHandle> n; mk(np,cp,ep,fp,n);
          fp.testSetTail(2147483647); CHECK(!fp.checkInvariants(np)); }
        { NodePool np; CellPool cp; EdgePool ep; FacePool fp; std::vector<NodeHandle> n; mk(np,cp,ep,fp,n);
          fp.testSetHead(fp.testMaterialized()); CHECK(!fp.checkInvariants(np)); }
        { NodePool np; CellPool cp; EdgePool ep; FacePool fp; std::vector<NodeHandle> n; mk(np,cp,ep,fp,n);
          fp.testSetNext(fp.testHead(), 99999); CHECK(!fp.checkInvariants(np)); }
        { NodePool np; CellPool cp; EdgePool ep; FacePool fp; std::vector<NodeHandle> n; mk(np,cp,ep,fp,n);
          fp.testSetPrev(fp.testHead(), 88888); CHECK(!fp.checkInvariants(np)); }
        { NodePool np; CellPool cp; EdgePool ep; FacePool fp; std::vector<NodeHandle> n; mk(np,cp,ep,fp,n);
          auto live = fp.activeBySequence();
          fp.testCorruptTableEntry(live[0], live[1]);
          CHECK(!fp.checkInvariants(np)); }
        { NodePool np; CellPool cp; EdgePool ep; FacePool fp; std::vector<NodeHandle> n; mk(np,cp,ep,fp,n);
          auto live = fp.activeBySequence();
          fp.testDropTableEntry(live[0]);
          CHECK(!fp.checkInvariants(np)); }
        { NodePool np; CellPool cp; EdgePool ep; FacePool fp; std::vector<NodeHandle> n; mk(np,cp,ep,fp,n);
          auto live = fp.activeBySequence();
          fp.testFlipSense(live[0], 0);       // boundary face sense flip
          CHECK(!fp.checkInvariants(np,cp)); } // recomputed sense disagrees
    }

    // (E) Real allocation-failure rollback across commit steps.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<8;++i) n.push_back(np.add({double(i),double(i%4),double(i%2)},{3,1}));
        CellPool cp;
        auto h = cp.add(ElementType::Hex8,
                        {n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]},{3,1},np);
        EdgePool ep; ep.insertCellEdges(np,cp,h);
        for (int stage=1; stage<=16; ++stage) {
            FacePool fp;
            topo::faultCountdown()=stage;
            bool threw=false;
            try { fp.insertCellFaces(np,cp,ep,h); }
            catch(const std::bad_alloc&){ threw=true; }
            topo::faultCountdown()=0;
            CHECK(fp.checkInvariants(np,cp));
            if (threw) CHECK(fp.liveCount()==0);
            else CHECK(fp.liveCount()==6);
        }
    }

    // (F) Removal under allocation pressure: pre-reserves before mutation.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<8;++i) n.push_back(np.add({double(i),double(i%4),double(i%2)},{3,1}));
        CellPool cp;
        auto h1 = cp.add(ElementType::Hex8,{n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]},{3,1},np);
        EdgePool ep; ep.insertCellEdges(np,cp,h1);
        FacePool fp; fp.insertCellFaces(np,cp,ep,h1);
        for (int stage=1; stage<=6; ++stage) {
            // rebuild fresh each stage
            NodePool np2; std::vector<NodeHandle> m;
            for (int i=0;i<8;++i) m.push_back(np2.add({double(i),double(i%4),double(i%2)},{3,1}));
            CellPool cp2; auto h2 = cp2.add(ElementType::Hex8,{m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7]},{3,1},np2);
            EdgePool ep2; ep2.insertCellEdges(np2,cp2,h2);
            FacePool fp2; fp2.insertCellFaces(np2,cp2,ep2,h2);
            const auto before = fp2.liveCount();
            topo::faultCountdown()=stage;
            bool threw=false;
            try { fp2.removeCellFaces(cp2,h2); }
            catch(const std::bad_alloc&){ threw=true; }
            topo::faultCountdown()=0;
            CHECK(fp2.checkInvariants(np2,cp2));
            if (threw) { CHECK(fp2.liveCount()==before); fp2.removeCellFaces(cp2,h2); }
            CHECK(fp2.checkInvariants(np2,cp2));
        }
    }

    static_assert(!std::is_nothrow_move_constructible_v<FacePool>);
    static_assert(std::is_nothrow_move_constructible_v<
                  topo::HandleAllocator<topo::FaceTag>>);

    // (G0) Empty public move uses a tokenless destination allocator.
    {
        FacePool src;
        const auto sourcePool = src.poolId();
        FacePool dst(std::move(src));
        CHECK(dst.poolId() == sourcePool && dst.liveCount() == 0);
        CHECK(src.poolId() == 0 && src.liveCount() == 0);
        CHECK(src.stats().materializedSlots == 0 && src.stats().consistent());
        FacePool nextOrdinary;
        CHECK(nextOrdinary.poolId() == sourcePool + 1);
    }

    // (G) FacePool public move construction qualification.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<5;++i) n.push_back(np.add({double(i),0,0},{3,1}));
        CellPool cp;
        auto cA = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        auto cB = cp.add(ElementType::Tet4,{n[0],n[2],n[1],n[4]},{3,1},np);
        EdgePool ep; ep.insertCellEdges(np,cp,cA); ep.insertCellEdges(np,cp,cB);
        FacePool src;
        src.insertCellFaces(np,cp,ep,cA);
        src.insertCellFaces(np,cp,ep,cB);
        std::array<NodeHandle,3> sh{n[0],n[1],n[2]};
        auto shared = src.find(np, std::span<const NodeHandle>(sh.data(),3));
        CHECK(shared.has_value());
        const auto hBefore = *shared;
        const auto poolBefore = src.poolId();
        const auto orderBefore = src.activeBySequence();
        const std::vector<FaceHandle> cABefore(
            src.cellFaces(cp, cA).begin(), src.cellFaces(cp, cA).end());
        const std::vector<FaceHandle> cBBefore(
            src.cellFaces(cp, cB).begin(), src.cellFaces(cp, cB).end());
        const auto sharedBefore = src.get(hBefore);
        std::vector<std::uint64_t> seqBefore;
        for (const auto& h : orderBefore) seqBefore.push_back(src.seqOf(h));
        const auto statsBefore = src.stats();
        FacePool dst(std::move(src));
        CHECK(dst.poolId() == poolBefore);
        CHECK(dst.isCurrent(hBefore));
        CHECK(dst.checkInvariants(np,cp));
        CHECK(dst.liveCount()==7);
        CHECK(dst.activeBySequence() == orderBefore);
        CHECK(std::vector<FaceHandle>(dst.cellFaces(cp, cA).begin(),
                  dst.cellFaces(cp, cA).end()) == cABefore);
        CHECK(std::vector<FaceHandle>(dst.cellFaces(cp, cB).begin(),
                  dst.cellFaces(cp, cB).end()) == cBBefore);
        for (std::size_t i = 0; i < orderBefore.size(); ++i)
            CHECK(dst.seqOf(orderBefore[i]) == seqBefore[i]);
        const auto sharedAfter = dst.get(hBefore);
        CHECK(sharedAfter.incident.size() == sharedBefore.incident.size());
        for (std::size_t i = 0; i < sharedBefore.incident.size(); ++i) {
            CHECK(sharedAfter.incident[i].cellSlot ==
                  sharedBefore.incident[i].cellSlot);
            CHECK(sharedAfter.incident[i].cellGen ==
                  sharedBefore.incident[i].cellGen);
            CHECK(sharedAfter.incident[i].sense ==
                  sharedBefore.incident[i].sense);
        }
        CHECK(dst.stats().materializedSlots == statsBefore.materializedSlots);
        CHECK(dst.stats().live == statsBefore.live);
        CHECK(dst.stats().reusableFree == statsBefore.reusableFree);
        CHECK(dst.stats().retired == statsBefore.retired);
        CHECK(dst.stats().totalCreated == statsBefore.totalCreated);
        CHECK(src.poolId()==0 && src.liveCount()==0);
        CHECK(src.boundNodePool()==0 && src.boundCellPool()==0);
        CHECK(src.activeBySequence().empty());
        CHECK(src.stats().materializedSlots == 0 && src.stats().live == 0 &&
              src.stats().reusableFree == 0 && src.stats().retired == 0 &&
              src.stats().totalCreated == 0 && src.stats().consistent());
        bool unusable = false;
        try { (void)src.insertCellFaces(np,cp,ep,cA); }
        catch (const std::logic_error&) { unusable = true; }
        CHECK(unusable);
        dst.removeCellFaces(cp,cB);           // removal works after move
        CHECK(dst.liveCount()==4);
        CHECK(dst.checkInvariants(np,cp));
    }

    // (H) Cell-EDGE incidence gate (review blocker): a duplicate cell whose
    // OWN edges were never inserted is rejected, even though an identical
    // cell already created the edge records.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<4;++i) n.push_back(np.add({double(i),0,0},{3,1}));
        CellPool cp;
        auto c1 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        auto c2 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        EdgePool ep; ep.insertCellEdges(np,cp,c1);   // only c1's edges
        CHECK(ep.hasCellEdges(cp,c1));
        CHECK(!ep.hasCellEdges(cp,c2));              // c2 has no edge incidence
        FacePool fp;
        bool threw=false;
        try { fp.insertCellFaces(np,cp,ep,c2); }
        catch(const std::invalid_argument&){ threw=true; }
        CHECK(threw);
        CHECK(fp.liveCount()==0);
        // c1 works (its edges are present).
        auto ff = fp.insertCellFaces(np,cp,ep,c1);
        CHECK(ff.size()==4 && fp.checkInvariants(np,cp));
        CHECK(topo::faceTopologyReady(np,cp,ep,fp));
    }

    // (I) Specific diagnostics: orientation vs incidence codes are distinct.
    {
        NodePool np; std::vector<NodeHandle> n;
        for (int i=0;i<5;++i) n.push_back(np.add({double(i),double(i%2),0},{3,1}));
        CellPool cp;
        // Same-winding shared face => orientation mismatch.
        auto c1 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
        auto c2 = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[4]},{3,1},np);
        EdgePool ep; ep.insertCellEdges(np,cp,c1); ep.insertCellEdges(np,cp,c2);
        FacePool fp;
        fp.insertCellFaces(np,cp,ep,c1);
        fp.insertCellFaces(np,cp,ep,c2);
        const auto rep = topo::auditFaceTopology(np,cp,ep,fp);
        bool orient=false;
        for (const auto& d : rep.items)
            orient |= d.code==diag::Code::FaceOrientationMismatch;
        CHECK(orient);                         // specific code, not GhostFace
    }

    // (J) CELL-AWARE corruption probes (review blocker: the cell-aware
    // invariant form previously segfaulted under ASan on a corrupt head_
    // because incidence traversal preceded structural validation). Every
    // corruption must be REPORTED - by checkInvariants(nodes),
    // checkInvariants(nodes,cells) AND auditFaceTopology - never crash.
    {
        auto scenario = [&](auto&& corrupt, bool senseOnly) {
            NodePool np; std::vector<NodeHandle> n;
            for (int i=0;i<4;++i) n.push_back(np.add({double(i),0,0},{3,1}));
            CellPool cp;
            auto c = cp.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},np);
            EdgePool ep; ep.insertCellEdges(np,cp,c);
            FacePool fp; fp.insertCellFaces(np,cp,ep,c);
            corrupt(fp);
            if (senseOnly) CHECK(fp.checkInvariants(np));   // not structural
            else           CHECK(!fp.checkInvariants(np));
            CHECK(!fp.checkInvariants(np,cp));              // no crash, false
            const auto rep = topo::auditFaceTopology(np,cp,ep,fp);
            CHECK(!rep.acceptable(diag::Severity::Warning));
        };
        scenario([](FacePool& f){ f.testSetHead(f.testMaterialized()); }, false);
        scenario([](FacePool& f){ f.testSetHead(2147483647); },          false);
        scenario([](FacePool& f){ f.testSetTail(-7); },                  false);
        scenario([](FacePool& f){ f.testSetNext(f.testHead(), 99999); }, false);
        scenario([](FacePool& f){ f.testSetPrev(f.testHead(), 88888); }, false);
        scenario([](FacePool& f){
            auto v=f.activeBySequence(); f.testDropTableEntry(v[0]); },  false);
        scenario([](FacePool& f){
            auto v=f.activeBySequence(); f.testCorruptTableEntry(v[0], v[1]); }, false);
        scenario([](FacePool& f){
            auto v=f.activeBySequence(); f.testFlipSense(v[0], 0); },    true);
    }
}


// ---------------------------------------------------------------------------
static void testCompleteTopologyFactories() {
    TEST("topo: complete edge and face factories preserve incremental identity");
    using topo::CellPool;
    using topo::EdgePool;
    using topo::ElementType;
    using topo::FacePool;
    using topo::NodeHandle;
    using topo::NodePool;

    auto requireSameHandle = [](const auto& a,const auto& b) {
        CHECK(a.slot == b.slot);
        CHECK(a.gen == b.gen);
    };
    auto requireEquivalent = [&](NodePool& nodes,CellPool& cells,
                                 bool markProtected) {
        EdgePool incrementalEdges;
        for (const auto& cell : cells.activeBySequence())
            incrementalEdges.insertCellEdges(nodes,cells,cell);
        topo::EdgePoolCompleteBuildMetrics edgeMetrics;
        auto completeEdges = EdgePool::buildComplete(
            nodes,cells,&edgeMetrics);

        FacePool incrementalFaces;
        for (const auto& cell : cells.activeBySequence())
            incrementalFaces.insertCellFaces(
                nodes,cells,incrementalEdges,cell);
        topo::FacePoolCompleteBuildMetrics faceMetrics;
        auto completeFaces = FacePool::buildComplete(
            nodes,cells,completeEdges,&faceMetrics);

        const auto cellsInOrder = cells.activeBySequence();
        if (markProtected && incrementalEdges.liveCount() != 0) {
            const auto ie = incrementalEdges.activeBySequence().front();
            const auto be = completeEdges.activeBySequence().front();
            incrementalEdges.setRequiredConstraint(ie,true);
            completeEdges.setRequiredConstraint(be,true);
            incrementalEdges.setClassification(ie,{1,77});
            completeEdges.setClassification(be,{1,77});
            const auto iff = incrementalFaces.activeBySequence().front();
            const auto bff = completeFaces.activeBySequence().front();
            incrementalFaces.setClassification(iff,{2,88});
            completeFaces.setClassification(bff,{2,88});
        }

        const auto oldEdges = incrementalEdges.activeBySequence();
        const auto bulkEdges = completeEdges.activeBySequence();
        CHECK(oldEdges.size() == bulkEdges.size());
        for (std::size_t i = 0; i < oldEdges.size(); ++i) {
            requireSameHandle(oldEdges[i],bulkEdges[i]);
            CHECK(incrementalEdges.seqOf(oldEdges[i]) ==
                  completeEdges.seqOf(bulkEdges[i]));
            const auto a = incrementalEdges.get(oldEdges[i]);
            const auto b = completeEdges.get(bulkEdges[i]);
            CHECK(a.a == b.a && a.b == b.b);
            CHECK(a.cls.dim == b.cls.dim && a.cls.tag == b.cls.tag);
            CHECK(a.requiredConstraint == b.requiredConstraint);
            CHECK(a.useCount == b.useCount);
        }
        for (const auto& cell : cellsInOrder) {
            const auto a = incrementalEdges.cellEdges(cells,cell);
            const auto b = completeEdges.cellEdges(cells,cell);
            CHECK(a.size() == b.size());
            for (std::size_t i = 0; i < a.size(); ++i)
                requireSameHandle(a[i],b[i]);
        }

        const auto oldFaces = incrementalFaces.activeBySequence();
        const auto bulkFaces = completeFaces.activeBySequence();
        CHECK(oldFaces.size() == bulkFaces.size());
        for (std::size_t i = 0; i < oldFaces.size(); ++i) {
            requireSameHandle(oldFaces[i],bulkFaces[i]);
            CHECK(incrementalFaces.seqOf(oldFaces[i]) ==
                  completeFaces.seqOf(bulkFaces[i]));
            const auto a = incrementalFaces.get(oldFaces[i]);
            const auto b = completeFaces.get(bulkFaces[i]);
            CHECK(a.shape == b.shape && a.cornerCount == b.cornerCount);
            CHECK(a.corners == b.corners);
            CHECK(a.cls.dim == b.cls.dim && a.cls.tag == b.cls.tag);
            CHECK(a.incident.size() == b.incident.size());
            for (std::size_t j = 0; j < a.incident.size(); ++j) {
                CHECK(a.incident[j].cellSlot == b.incident[j].cellSlot);
                CHECK(a.incident[j].cellGen == b.incident[j].cellGen);
                CHECK(a.incident[j].sense == b.incident[j].sense);
            }
        }
        for (const auto& cell : cellsInOrder) {
            const auto a = incrementalFaces.cellFaces(cells,cell);
            const auto b = completeFaces.cellFaces(cells,cell);
            CHECK(a.size() == b.size());
            for (std::size_t i = 0; i < a.size(); ++i)
                requireSameHandle(a[i],b[i]);
        }

        const auto oldEdgeStats = incrementalEdges.stats();
        const auto bulkEdgeStats = completeEdges.stats();
        CHECK(oldEdgeStats.live == bulkEdgeStats.live);
        CHECK(oldEdgeStats.materializedSlots == bulkEdgeStats.materializedSlots);
        CHECK(oldEdgeStats.totalCreated == bulkEdgeStats.totalCreated);
        const auto oldFaceStats = incrementalFaces.stats();
        const auto bulkFaceStats = completeFaces.stats();
        CHECK(oldFaceStats.live == bulkFaceStats.live);
        CHECK(oldFaceStats.materializedSlots == bulkFaceStats.materializedSlots);
        CHECK(oldFaceStats.totalCreated == bulkFaceStats.totalCreated);
        // Capacity is intentionally not compared: the direct factories make
        // one exact reserve while the legacy path grows geometrically.
        CHECK(topo::edgeTopologyReady(nodes,cells,completeEdges));
        CHECK(topo::faceTopologyReady(
            nodes,cells,completeEdges,completeFaces));
        CHECK(edgeMetrics.validationCellVisits == cellsInOrder.size());
        CHECK(edgeMetrics.buildCellVisits == cellsInOrder.size());
        CHECK(edgeMetrics.edgeUses == cellsInOrder.size() * 6U);
        CHECK(edgeMetrics.uniqueEdges == completeEdges.liveCount());
        CHECK(edgeMetrics.finalEdgeAudits == 1);
        CHECK(edgeMetrics.perCellGlobalAudits == 0);
        CHECK(edgeMetrics.runComplete);
        CHECK(faceMetrics.validationCellVisits == cellsInOrder.size());
        CHECK(faceMetrics.buildCellVisits == cellsInOrder.size());
        CHECK(faceMetrics.faceUses == cellsInOrder.size() * 4U);
        CHECK(faceMetrics.uniqueFaces == completeFaces.liveCount());
        CHECK(faceMetrics.suppliedEdgeAudits == 1);
        CHECK(faceMetrics.finalFaceAudits == 1);
        CHECK(faceMetrics.perCellGlobalAudits == 0);
        CHECK(faceMetrics.runComplete);
    };

    // Empty pools retain the established tokenful/unbound identity.
    {
        NodePool nodes;
        CellPool cells;
        requireEquivalent(nodes,cells,false);
    }
    // Boundary-only, two-sided interior, nonmanifold, and protected rows.
    for (int owners = 1; owners <= 3; ++owners) {
        NodePool nodes;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 6; ++i)
            n.push_back(nodes.add({double(i),double(i % 2),double(i % 3)},
                                  {3,1}));
        CellPool cells;
        cells.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},
                  {3,1},nodes);
        if (owners >= 2)
            cells.add(ElementType::Tet4,{n[0],n[2],n[1],n[4]},
                      {3,1},nodes);
        if (owners >= 3)
            cells.add(ElementType::Tet4,{n[0],n[1],n[2],n[5]},
                      {3,1},nodes);
        requireEquivalent(nodes,cells,owners == 2);
    }
    // Source allocator churn: the current cell reuses a retired slot with a
    // new generation, and both construction routes must retain that identity.
    {
        NodePool nodes;
        std::vector<NodeHandle> n;
        for (int i = 0; i < 6; ++i)
            n.push_back(nodes.add({double(i),double(i % 2),double(i % 3)},
                                  {3,1}));
        CellPool cells;
        const auto retired = cells.add(
            ElementType::Tet4,{n[0],n[1],n[2],n[3]},{3,1},nodes);
        cells.remove(retired);
        const auto reused = cells.add(
            ElementType::Tet4,{n[0],n[2],n[1],n[4]},{3,1},nodes);
        CHECK(reused.slot == retired.slot && reused.gen != retired.gen);
        cells.add(ElementType::Tet4,{n[0],n[1],n[2],n[5]},
                  {3,1},nodes);
        requireEquivalent(nodes,cells,false);
    }

    // The direct factories must also preserve the canonical semantic stream,
    // not merely agree row-by-row with the incremental construction path.
    // Separate source pools deliberately give the workspaces different pool
    // identities; V2 normalizes those identities to persistent ranks.
    auto makeWorkspace = [](bool complete) {
        NodePool nodes;
        std::vector<NodeHandle> n;
        n.push_back(nodes.add({0.0,0.0,0.0},{3,1}));
        n.push_back(nodes.add({1.0,0.0,0.0},{3,1}));
        n.push_back(nodes.add({0.0,1.0,0.0},{3,1}));
        n.push_back(nodes.add({0.0,0.0,1.0},{3,1}));
        n.push_back(nodes.add({0.0,0.0,-1.0},{3,1}));
        CellPool cells;
        cells.add(ElementType::Tet4,{n[0],n[1],n[2],n[3]},
                  {3,1},nodes);
        cells.add(ElementType::Tet4,{n[0],n[2],n[1],n[4]},
                  {3,1},nodes);

        auto edges = [&]() -> EdgePool {
            if (complete)
                return EdgePool::buildComplete(nodes,cells);
            EdgePool result;
            for (const auto& cell : cells.activeBySequence())
                result.insertCellEdges(nodes,cells,cell);
            return result;
        }();
        auto faces = [&]() -> FacePool {
            if (complete)
                return FacePool::buildComplete(nodes,cells,edges);
            FacePool result;
            for (const auto& cell : cells.activeBySequence())
                result.insertCellFaces(nodes,cells,edges,cell);
            return result;
        }();
        const auto edge = edges.activeBySequence().front();
        edges.setRequiredConstraint(edge,true);
        edges.setClassification(edge,{1,77});
        const auto face = faces.activeBySequence().front();
        faces.setClassification(face,{2,88});
        return std::make_unique<topo::TopologyWorkspace>(
            std::move(nodes),std::move(cells),
            std::move(edges),std::move(faces));
    };
    const auto incrementalWorkspace = makeWorkspace(false);
    const auto completeWorkspace = makeWorkspace(true);
    CHECK(incrementalWorkspace->semanticTopologyHashV2() ==
          completeWorkspace->semanticTopologyHashV2());
    CHECK(incrementalWorkspace->semanticTopologyHash() ==
          completeWorkspace->semanticTopologyHash());
#ifdef PAACS_TOPO_TESTING
    CHECK(incrementalWorkspace->testV2Tokens() ==
          completeWorkspace->testV2Tokens());
#endif
}



// Registration hook consumed by test_main.cpp.
namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* topologyTests(std::size_t& count) {
    static const Entry kTopo[] = {
        {"topo.storage",              testStorageNodeCell},
        {"topo.storage.audit",        testStorageAuditAndSafety},
        {"topo.storage.scale",        testStorageScale},
        {"topo.storage.api",          testStorageApiAndListSafety},
        {"topo.edges",                testEdgeTopology},
        {"topo.edges.lifecycle",      testEdgeLifecycle},
        {"topo.edges.hardening",      testEdgeHardening},
        {"topo.faces",                testFaceTopology},
        {"topo.faces.hardening",      testFaceHardening},
        {"topo.complete_factories",   testCompleteTopologyFactories},
    };
    count = sizeof(kTopo) / sizeof(kTopo[0]);
    return kTopo;
}
} // namespace paacs_tests
