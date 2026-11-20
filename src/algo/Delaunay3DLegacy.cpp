// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Delaunay3D.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/topo/TetArena.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace paacs::algo::detail {
namespace {

using predicates::Sign;

// Heuristic-only double orientation for the walk's steepest-descent face
// choice; every DECISION uses predicates::orient3d.
inline double orientHeur(const Vec3& a, const Vec3& b, const Vec3& c,
                         const Vec3& d) {
    return (b - a).cross(c - a).dot(d - a);
}

// Face k is opposite local vertex k; for a positively oriented FINITE tet
// the kFace winding points outward: orient3d(f0,f1,f2, v[k]) < 0.
constexpr int kFace[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};

using ACell = topo::TetArena::Cell;

std::uint64_t morton3(std::uint32_t x, std::uint32_t y, std::uint32_t z) {
    auto part = [](std::uint64_t v) {
        v &= 0x1FFFFF;
        v = (v | v << 32) & 0x1F00000000FFFFULL;
        v = (v | v << 16) & 0x1F0000FF0000FFULL;
        v = (v | v << 8)  & 0x100F00F00F00F00FULL;
        v = (v | v << 4)  & 0x10C30C30C30C30C3ULL;
        v = (v | v << 2)  & 0x1249249249249249ULL;
        return v;
    };
    return part(x) | (part(y) << 1) | (part(z) << 2);
}

} // namespace

// Bowyer-Watson with ONE symbolic vertex at infinity (index n). The old
// big-enclosing-simplex approach fails on flat point clouds: near-planar
// real tets have circumspheres so large they contain the finite-coordinate
// super vertices, so those tets can never survive - measured as zero output
// tets on clouds flattened to 1e-9. Infinite cells are handled by witness-
// based symbolic predicates instead of geometry:
//   an infinite cell (finite facet F + INF) conflicts with p  iff  p lies
//   strictly on the opposite side of F from the adjacent finite cell's
//   witness vertex, or exactly on F's plane (coplanar hull extension).
std::vector<Tet> delaunay3DLegacy64(const std::vector<Vec3>& pts,
                                    Delaunay3DStats* stats) {
    Delaunay3DStats localStats;
    Delaunay3DStats& st = stats ? *stats : localStats;
    st = Delaunay3DStats{};
    std::vector<Tet> out;
    const Index n = static_cast<Index>(pts.size());
    st.inputPoints = static_cast<std::uint64_t>(pts.size());
    // Structured rejection of non-finite input (review finding: a NaN point
    // previously crashed inside Morton quantization). Meshers guarantee
    // finite nodes, so this is a caller-contract violation, reported loudly
    // with a count rather than undefined behavior.
    {
        std::uint64_t nonFinite = 0;
        for (const auto& p : pts)
            if (!std::isfinite(p.x) || !std::isfinite(p.y) ||
                !std::isfinite(p.z))
                ++nonFinite;
        if (nonFinite > 0)
            throw std::invalid_argument(
                "delaunay3D: " + std::to_string(nonFinite) +
                " input point(s) have non-finite coordinates");
    }
    if (n < 4) {
        st.pointsSkipped = static_cast<std::uint64_t>(n);
        st.runComplete = true;
        return out;
    }
    const Index INF = n;

    // Normalized working frame (PM-NUM-001): predicates and heuristics see
    // unit-scale coordinates regardless of model position/size. The frame is
    // a uniform similarity, so Delaunay topology is unchanged.
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const auto& p : pts) {
        lo.x = std::min(lo.x, p.x); hi.x = std::max(hi.x, p.x);
        lo.y = std::min(lo.y, p.y); hi.y = std::max(hi.y, p.y);
        lo.z = std::min(lo.z, p.z); hi.z = std::max(hi.z, p.z);
    }
    const Vec3 mid = (lo + hi) * 0.5;
    const double span =
        std::max({hi.x - lo.x, hi.y - lo.y, hi.z - lo.z, 1e-300});
    const double inv = 1.0 / span;

    std::vector<Vec3> P(static_cast<size_t>(n) + 1);
    for (Index i = 0; i < n; ++i) P[size_t(i)] = (pts[size_t(i)] - mid) * inv;
    const double qn = std::numeric_limits<double>::quiet_NaN();
    P[size_t(INF)] = {qn, qn, qn};   // must never be used geometrically

    // Morton insertion order (spatial coherence => short walks).
    std::vector<Index> order(static_cast<size_t>(n));
    {
        std::vector<std::uint64_t> key(static_cast<size_t>(n));
        for (Index i = 0; i < n; ++i) {
            const Vec3& q = P[size_t(i)];
            auto quant = [](double x) {
                return static_cast<std::uint32_t>(
                    std::clamp((x + 0.5) * 2097151.0, 0.0, 2097151.0));
            };
            key[size_t(i)] = morton3(quant(q.x), quant(q.y), quant(q.z));
            order[size_t(i)] = i;
        }
        std::stable_sort(order.begin(), order.end(),
                         [&](Index a, Index b) {
                             return key[size_t(a)] < key[size_t(b)];
                         });
    }

    // M4: reusable arena. Cavity commits free their dead cells and the
    // free list feeds the replacement allocations, so construction memory
    // tracks the live complex, not construction history (M3-A retained
    // every dead cell: ~145-208 MB at 192k output tets).
    topo::TetArena arena;
    auto isInf = [&](const ACell& t) { return t.v[3] == INF; };

    auto makeTet = [&](Index va, Index vb, Index vc, Index vd) -> int {
        std::array<Index, 4> v{va, vb, vc, vd};
        // Canonicalize the infinite vertex to slot 3; infinite cells carry
        // no geometric orientation.
        for (int k = 0; k < 3; ++k)
            if (v[size_t(k)] == INF) std::swap(v[size_t(k)], v[3]);
        if (v[3] != INF) {
            const Sign o = predicates::orient3d(
                P[size_t(v[0])], P[size_t(v[1])], P[size_t(v[2])],
                P[size_t(v[3])]);
            if (o == Sign::Negative) std::swap(v[2], v[3]);
            if (o == Sign::Zero)
                throw std::logic_error("delaunay3D: degenerate tet requested");
        }
        return arena.alloc(v);
    };

    // ---- Bootstrap: first non-degenerate quadruple in Morton order --------
    // (exact test, so 1e-9-flat clouds seed correctly), then one finite tet
    // plus four infinite cells over its facets.
    {
        Index a = order[0], b = kInvalidIndex, c = kInvalidIndex,
              d = kInvalidIndex;
        size_t bi = 0;
        for (size_t i = 1; i < order.size() && b == kInvalidIndex; ++i) {
            ++st.bootstrapDistinctTests;
            if (!(P[size_t(order[i])].x == P[size_t(a)].x &&
                  P[size_t(order[i])].y == P[size_t(a)].y &&
                  P[size_t(order[i])].z == P[size_t(a)].z)) {
                b = order[i];
                bi = i;
            }
        }
        for (size_t i = bi + 1; i < order.size() && d == kInvalidIndex; ++i) {
            for (size_t j = i + 1; j < order.size() && d == kInvalidIndex;
                 ++j) {
                ++st.bootstrapOrientationTests;
                if (predicates::orient3d(P[size_t(a)], P[size_t(b)],
                                         P[size_t(order[i])],
                                         P[size_t(order[j])]) != Sign::Zero) {
                    c = order[i];
                    d = order[j];
                }
            }
        }
        if (d == kInvalidIndex) {
            // Structured limitation: the entire input is exactly collinear /
            // coplanar - no 3D triangulation exists.
            st.pointsSkipped = static_cast<std::uint64_t>(n);
            st.runComplete = true;
            return out;
        }
        const int T = makeTet(a, b, c, d);
        int inf[4];
        for (int k = 0; k < 4; ++k) {
            inf[k] = makeTet(arena.cell(T).v[size_t(kFace[k][0])],
                             arena.cell(T).v[size_t(kFace[k][1])],
                             arena.cell(T).v[size_t(kFace[k][2])], INF);
            arena.cell(inf[k]).nb[3] = T;
            arena.cell(T).nb[size_t(k)] = inf[k];
        }
        // Infinite-infinite adjacency: cells over facets j,k share the edge
        // opposite {j,k} plus INF.
        auto sideFace = [&](int t, Index e0, Index e1) {
            for (int k = 0; k < 3; ++k) {
                const Index x = arena.cell(t).v[size_t(kFace[k][0])],
                            y = arena.cell(t).v[size_t(kFace[k][1])],
                            z = arena.cell(t).v[size_t(kFace[k][2])];
                const bool has0 = x == e0 || y == e0 || z == e0;
                const bool has1 = x == e1 || y == e1 || z == e1;
                const bool hasI = x == INF || y == INF || z == INF;
                if (has0 && has1 && hasI) return k;
            }
            return -1;
        };
        for (int j = 0; j < 4; ++j)
            for (int k = j + 1; k < 4; ++k) {
                Index shared[2];
                int m = 0;
                for (Index vv : arena.cell(T).v)
                    if (vv != arena.cell(T).v[size_t(j)] &&
                        vv != arena.cell(T).v[size_t(k)])
                        shared[m++] = vv;
                const int fj = sideFace(inf[j], shared[0], shared[1]);
                const int fk = sideFace(inf[k], shared[0], shared[1]);
                arena.cell(inf[j]).nb[size_t(fj)] = inf[k];
                arena.cell(inf[k]).nb[size_t(fk)] = inf[j];
            }
        st.pointsInserted += 4;
    }
    const std::array<Index, 4> seeded{
        arena.cell(0).v[0], arena.cell(0).v[1],
        arena.cell(0).v[2], arena.cell(0).v[3]};

    // Conflict predicate. Finite: exact insphere, include-on-zero (merges
    // cospherical shells - required by structured lattices). Infinite:
    // witness rule described above, include-on-coplanar.
    auto conflicts = [&](const ACell& t, const Vec3& p) {
        ++st.conflictQueries;
        if (!isInf(t)) {
            ++st.conflictInsphereTests;
            return predicates::insphere(P[size_t(t.v[0])], P[size_t(t.v[1])],
                                        P[size_t(t.v[2])], P[size_t(t.v[3])],
                                        p) != Sign::Negative;
        }
        const int fin = t.nb[3];
        const ACell& ft = arena.cell(fin);
        Index w = kInvalidIndex;
        for (Index vv : ft.v)
            if (vv != t.v[0] && vv != t.v[1] && vv != t.v[2]) w = vv;
        st.conflictOrientationTests += 2;
        const Sign ow = predicates::orient3d(
            P[size_t(t.v[0])], P[size_t(t.v[1])], P[size_t(t.v[2])],
            P[size_t(w)]);
        const Sign op = predicates::orient3d(
            P[size_t(t.v[0])], P[size_t(t.v[1])], P[size_t(t.v[2])], p);
        if (op == ow) return false;        // strictly on the hull's inside
        if (op != Sign::Zero) return true; // strictly outside this facet
        // p exactly on the facet plane (ubiquitous on planar boundary
        // patches such as cylinder caps). Blanket inclusion floods the
        // whole coplanar hull region (measured: ~178 resurrections per
        // insertion). The correct rule is the in-plane Delaunay condition:
        // conflict iff p lies inside the facet's circumCIRCLE, evaluated
        // exactly as the circumsphere through the facet and the finite
        // witness w (that sphere meets the facet plane in precisely the
        // circumcircle). Consistency bonus: p inside that sphere means the
        // finite neighbor conflicts too, so the facet becomes an internal
        // cavity face and can never spawn a degenerate coplanar tet.
        if (ow == Sign::Positive)
        {
            ++st.conflictInsphereTests;
            return predicates::insphere(P[size_t(t.v[0])], P[size_t(t.v[1])],
                                        P[size_t(t.v[2])], P[size_t(w)],
                                        p) != Sign::Negative;
        }
        ++st.conflictInsphereTests;
        return predicates::insphere(P[size_t(t.v[0])], P[size_t(t.v[2])],
                                    P[size_t(t.v[1])], P[size_t(w)],
                                    p) != Sign::Negative;
    };

    int lastAlive = 0;
    std::vector<int> bad, stack, work, stillBad;
    struct Pending { int newTet = -1; int side = -1; Index e0 = 0, e1 = 0; };
    struct ShellEdgePairTable {
        struct Slot {
            Index e0 = 0, e1 = 0;
            int newTet = -1, side = -1;
            std::uint32_t epoch = 0;
            std::uint8_t uses = 0;
        };
        enum class Result : std::uint8_t { First, Second, Third };
        std::vector<Slot> slots;
        std::uint32_t epoch = 0;
        std::size_t mask = 0;
        std::size_t unmatched = 0;
        std::uint64_t* probes = nullptr;

        static std::uint64_t mix(std::uint64_t value) noexcept {
            value ^= value >> 30;
            value *= 0xbf58476d1ce4e5b9ULL;
            value ^= value >> 27;
            value *= 0x94d049bb133111ebULL;
            return value ^ (value >> 31);
        }
        void begin(std::size_t expected, std::uint64_t& probeCounter) {
            std::size_t capacity = 8;
            while (capacity < expected * 2U) capacity <<= 1U;
            if (slots.size() < capacity) slots.resize(capacity);
            mask = slots.size() - 1U;
            probes = &probeCounter;
            unmatched = 0;
            if (++epoch == 0) {
                for (auto& slot : slots) slot.epoch = 0;
                epoch = 1;
            }
        }
        Result insert(const Pending& value, Pending& first) noexcept {
            const std::uint64_t a = static_cast<std::uint64_t>(value.e0);
            const std::uint64_t b = static_cast<std::uint64_t>(value.e1);
            std::size_t position = static_cast<std::size_t>(
                mix(a) ^ (mix(b) + 0x9e3779b97f4a7c15ULL)) & mask;
            for (;;) {
                if (probes) ++*probes;
                Slot& slot = slots[position];
                if (slot.epoch != epoch) {
                    slot = Slot{value.e0, value.e1, value.newTet,
                                value.side, epoch, 1};
                    ++unmatched;
                    return Result::First;
                }
                if (slot.e0 == value.e0 && slot.e1 == value.e1) {
                    first = Pending{slot.newTet, slot.side, slot.e0, slot.e1};
                    if (slot.uses != 1) return Result::Third;
                    slot.uses = 2;
                    --unmatched;
                    return Result::Second;
                }
                position = (position + 1U) & mask;
            }
        }
    } pairTable;
    struct ShellFaceT { Index f0, f1, f2; int outer, outerFace; };
    std::vector<ShellFaceT> shell;

    for (Index oi = 0; oi < n; ++oi) {
        const Index pi = order[size_t(oi)];
        if (std::find(seeded.begin(), seeded.end(), pi) != seeded.end())
            continue;  // bootstrap vertices
        const Vec3 p = P[size_t(pi)];

        // ---- 1. Locate ------------------------------------------------------
        int cur = lastAlive;
        if (!arena.cell(cur).alive) {
            cur = -1;
            for (int t = arena.slotCount() - 1; t >= 0; --t)
                if (arena.cell(t).alive) { cur = t; break; }
        }
        int located = -1;
        const int maxSteps = 4 * arena.slotCount() + 64;
        for (int step = 0; cur >= 0 && step < maxSteps; ++step) {
            ++st.walkCellVisits;
            const ACell& t = arena.cell(cur);
            if (isInf(t)) {
                // p conflicts with an infinite cell it walked into (see
                // note below); otherwise continue inward.
                if (conflicts(t, p)) { located = cur; break; }
                cur = t.nb[3];
                continue;
            }
            int bestFace = -1;
            double bestVal = -1e300;
            for (int k = 0; k < 4; ++k) {
                ++st.walkOrientationTests;
                const Vec3 &f0 = P[size_t(t.v[size_t(kFace[k][0])])],
                           &f1 = P[size_t(t.v[size_t(kFace[k][1])])],
                           &f2 = P[size_t(t.v[size_t(kFace[k][2])])];
                if (predicates::orient3d(f0, f1, f2, p) != Sign::Positive)
                    continue;   // not strictly outside across this face
                const double h = orientHeur(f0, f1, f2, p);
                if (h > bestVal) { bestVal = h; bestFace = k; }
            }
            if (bestFace < 0) { located = cur; break; }
            cur = t.nb[size_t(bestFace)];
        }

        // Duplicate input: bitwise-equal to a vertex of the located tet.
        if (located >= 0) {
            bool dup = false;
            for (Index vv : arena.cell(located).v) {
                if (vv == INF) continue;
                dup = dup || (P[size_t(vv)].x == p.x &&
                              P[size_t(vv)].y == p.y &&
                              P[size_t(vv)].z == p.z);
            }
            if (dup) { ++st.pointsDuplicate; continue; }
        }

        // ---- 2. Seed the conflict cavity ------------------------------------
        int seed = -1;
        if (located >= 0 && conflicts(arena.cell(located), p)) {
            seed = located;
        } else if (located >= 0) {
            for (int k = 0; k < 4 && seed < 0; ++k) {
                const int nb = arena.cell(located).nb[size_t(k)];
                if (nb >= 0 && arena.cell(nb).alive &&
                    conflicts(arena.cell(nb), p))
                    seed = nb;
            }
        }
        if (seed < 0) {
            ++st.walkFallbackScans;
            for (int t = 0; t < arena.slotCount(); ++t)
                if (arena.cell(t).alive) {
                    ++st.fallbackCellTests;
                    if (conflicts(arena.cell(t), p)) {
                    seed = t;
                    break;
                    }
                }
        }
        if (seed < 0) { ++st.pointsSkipped; continue; }

        // ---- 3. Connected conflict region (BFS over adjacency) --------------
        bad.clear();
        stack.clear();
        stack.push_back(seed);
        arena.beginEpoch();               // O(1) invalidation of old marks
        arena.markVisited(seed);          // marked == in the candidate cavity
        while (!stack.empty()) {
            const int t = stack.back();
            stack.pop_back();
            bad.push_back(t);
            ++st.cavityCellVisits;
            for (int k = 0; k < 4; ++k) {
                ++st.cavityNeighbourVisits;
                const int nb = arena.cell(t).nb[size_t(k)];
                if (nb < 0 || arena.visited(nb)) continue;
                if (!conflicts(arena.cell(nb), p)) continue;
                arena.markVisited(nb);
                stack.push_back(nb);
            }
        }

        // ---- 3b. Trim to a valid (star-shaped) cavity ------------------------
        // Every FINITE shell face must see p strictly from the correct side,
        // or its owning bad cell is resurrected. Faces containing INF need
        // no geometric check. For a finite shell face owned by an INFINITE
        // bad cell the required side is opposite the surviving finite
        // neighbor's witness vertex.
        //
        // Worklist algorithm: re-test a bad cell only when its shell can
        // have changed (a neighbor was resurrected). The earlier
        // pass-until-stable loop re-tested every face of every bad cell per
        // pass; on planar boundary patches (points exactly coplanar with
        // shell faces are common there) that produced ~178 resurrections
        // and >300 predicate tests PER INSERTION - measured as a 41% exact
        // orient3d rate and a superlinear series. Each cell is resurrected
        // at most once and each pop costs <= 4 predicate tests.
        auto shellOk = [&](int bt) {
            ++st.trimCellVisits;
            const ACell& t = arena.cell(bt);
            for (int k = 0; k < 4; ++k) {
                const int nb = t.nb[size_t(k)];
                if (nb >= 0 && arena.visited(nb)) continue;   // internal
                const Index f0 = t.v[size_t(kFace[k][0])],
                            f1 = t.v[size_t(kFace[k][1])],
                            f2 = t.v[size_t(kFace[k][2])];
                if (f0 == INF || f1 == INF || f2 == INF) continue;
                ++st.trimOrientationTests;
                const Sign op = predicates::orient3d(
                    P[size_t(f0)], P[size_t(f1)], P[size_t(f2)], p);
                if (t.v[3] != INF) {
                    if (op != Sign::Negative) return false;
                } else {
                    const ACell& ft = arena.cell(nb);
                    Index w = kInvalidIndex;
                    for (Index vv : ft.v)
                        if (vv != f0 && vv != f1 && vv != f2) w = vv;
                    const Sign ow = predicates::orient3d(
                        P[size_t(f0)], P[size_t(f1)], P[size_t(f2)],
                        P[size_t(w)]);
                    ++st.trimOrientationTests;
                    if (op == Sign::Zero || op == ow) return false;
                }
            }
            return true;
        };
        {
            work = bad;
            while (!work.empty()) {
                const int bt = work.back();
                work.pop_back();
                if (!arena.visited(bt)) continue;   // already resurrected
                if (shellOk(bt)) continue;
                arena.unmark(bt);                   // resurrect
                ++st.cavityTrimResurrections;
                // Neighbors still in the cavity gained a shell face.
                for (int k = 0; k < 4; ++k) {
                    const int nb = arena.cell(bt).nb[size_t(k)];
                    if (nb >= 0 && arena.visited(nb))
                        work.push_back(nb);
                }
            }
            stillBad.clear();
            stillBad.reserve(bad.size());
            for (const int bt : bad)
                if (arena.visited(bt)) stillBad.push_back(bt);
            bad.swap(stillBad);
        }
        st.peakCavityCells = std::max(
            st.peakCavityCells, static_cast<std::uint64_t>(bad.size()));
        if (bad.empty()) { ++st.pointsSkipped; continue; }

        // ---- 4. Commit: snapshot shell, free cavity, allocate from the
        //         freed slots, wire adjacency --------------------------------
        // The M3-A commit read cavity cells while creating replacements,
        // which is impossible once freed slots are immediately reused. So:
        // (a) snapshot every shell face (vertices + surviving neighbor and
        // its back-face index) in deterministic BFS order, (b) free the
        // cavity - its slots go to the free list, (c) allocate replacements
        // (drawing from exactly those slots: near-zero net growth per
        // insertion), (d) wire. Creation SEQUENCE numbers, not slot
        // numbers, provide deterministic output order, so the emitted
        // triangulation is bit-identical to M3-A (delaunay.golden).
        using ShellFace = ShellFaceT;
        shell.clear();
        for (const int bt : bad) {
            for (int k = 0; k < 4; ++k) {
                const int nb = arena.cell(bt).nb[size_t(k)];
                if (nb >= 0 && arena.visited(nb)) continue;   // internal
                int outerFace = -1;
                if (nb >= 0)
                    for (int kk = 0; kk < 4; ++kk)
                        if (arena.cell(nb).nb[size_t(kk)] == bt) outerFace = kk;
                shell.push_back({arena.cell(bt).v[size_t(kFace[k][0])],
                                 arena.cell(bt).v[size_t(kFace[k][1])],
                                 arena.cell(bt).v[size_t(kFace[k][2])],
                                 nb, outerFace});
                ++st.shellFacesEmitted;
            }
        }
        st.peakShellFaces = std::max(
            st.peakShellFaces, static_cast<std::uint64_t>(shell.size()));
        for (const int bt : bad) arena.free(bt);

        pairTable.begin(shell.size() * 3U, st.shellEdgePairProbes);
        st.peakShellPairTableSlots = std::max(
            st.peakShellPairTableSlots,
            static_cast<std::uint64_t>(pairTable.slots.size()));
        int firstNew = -1;
        for (const ShellFace& f : shell) {
            const int nt = makeTet(f.f0, f.f1, f.f2, pi);
            if (firstNew < 0) firstNew = nt;
            const int pFace = [&] {
                for (int kk = 0; kk < 4; ++kk)
                    if (arena.cell(nt).v[size_t(kk)] == pi) return kk;
                return -1;
            }();
            arena.cell(nt).nb[size_t(pFace)] = f.outer;
            if (f.outer >= 0)
                arena.cell(f.outer).nb[size_t(f.outerFace)] = nt;
            for (int kk = 0; kk < 4; ++kk) {
                if (kk == pFace) continue;
                Index e[2];
                int m = 0;
                for (int j = 0; j < 3; ++j) {
                    const Index vv = arena.cell(nt).v[size_t(kFace[kk][j])];
                    if (vv != pi) e[m++] = vv;
                }
                if (e[0] > e[1]) std::swap(e[0], e[1]);
                const Pending current{nt, kk, e[0], e[1]};
                Pending first{};
                ++st.shellEdgePairInsertions;
                const auto result = pairTable.insert(current, first);
                if (result == ShellEdgePairTable::Result::Third)
                    throw std::logic_error(
                        "delaunay3D: cavity shell edge used more than twice");
                if (result == ShellEdgePairTable::Result::Second) {
                    arena.cell(first.newTet).nb[size_t(first.side)] = nt;
                    arena.cell(nt).nb[size_t(kk)] = first.newTet;
                }
            }
        }
        if (pairTable.unmatched != 0)
            throw std::logic_error("delaunay3D: cavity shell is not closed");
        if (firstNew >= 0)
            lastAlive = firstNew;
        ++st.pointsInserted;
    }

    // ---- Audit: reciprocal adjacency (construction invariant) --------------
    {
        // Review finding: "some back-reference" is weaker than required.
        // Verify the neighbor's back-pointing face carries EXACTLY the
        // shared vertex set of this face.
        bool ok = true;
        for (int t = 0; t < arena.slotCount() && ok; ++t) {
            if (!arena.cell(t).alive) continue;
            for (int k = 0; k < 4 && ok; ++k) {
                ++st.adjacencyAuditFaceVisits;
                const int nb = arena.cell(t).nb[size_t(k)];
                if (nb < 0) continue;
                if (!arena.cell(nb).alive) { ok = false; break; }
                std::array<Index, 3> mine{
                    arena.cell(t).v[size_t(kFace[k][0])],
                    arena.cell(t).v[size_t(kFace[k][1])],
                    arena.cell(t).v[size_t(kFace[k][2])]};
                std::sort(mine.begin(), mine.end());
                bool matched = false;
                for (int kk = 0; kk < 4; ++kk) {
                    if (arena.cell(nb).nb[size_t(kk)] != t) continue;
                    std::array<Index, 3> theirs{
                        arena.cell(nb).v[size_t(kFace[kk][0])],
                        arena.cell(nb).v[size_t(kFace[kk][1])],
                        arena.cell(nb).v[size_t(kFace[kk][2])]};
                    std::sort(theirs.begin(), theirs.end());
                    matched = matched || theirs == mine;
                }
                ok = ok && matched;
            }
        }
        st.auditPassed = ok;
        if (!ok)
            throw std::logic_error(
                "delaunay3D: adjacency audit failed (shared-face identity)");
    }

    // ---- Output: finite live cells in CREATION-SEQUENCE order -------------
    // Slot numbers are reuse-scrambled; the sequence stamp reproduces the
    // M3-A append order exactly (verified by delaunay.golden and the
    // cross-thread mesh hash). Sort-based dedup retained as insurance.
    std::vector<std::pair<std::uint64_t, Tet>> bySeq;
    for (std::int32_t sSlot = 0; sSlot < arena.slotCount(); ++sSlot) {
        if (!arena.alive(sSlot)) continue;
        const ACell& c = arena.cell(sSlot);
        if (c.v[3] == INF) continue;
        bySeq.push_back({c.seq, Tet{c.v}});
    }
    std::sort(bySeq.begin(), bySeq.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<std::pair<std::array<Index, 4>, Index>> keyed;
    for (const auto& [seq, tet] : bySeq) {
        std::array<Index, 4> key = tet.v;
        std::sort(key.begin(), key.end());
        keyed.push_back({key, static_cast<Index>(out.size())});
        out.push_back(tet);
    }
    std::sort(keyed.begin(), keyed.end());
    std::vector<Index> keep;
    keep.reserve(keyed.size());
    for (size_t i = 0; i < keyed.size(); ++i) {
        if (i > 0 && keyed[i].first == keyed[i - 1].first) continue;
        keep.push_back(keyed[i].second);
    }
    std::sort(keep.begin(), keep.end());   // creation order (determinism)
    std::vector<Tet> dedup;
    dedup.reserve(keep.size());
    for (const Index id : keep) dedup.push_back(out[size_t(id)]);
    st.outputCells = static_cast<std::uint64_t>(dedup.size());

    // Output audit: every point accounted as inserted must appear in at
    // least one finite tetrahedron (no silent loss can hide in the
    // super-structure). Duplicates are exempt by definition.
    if (!dedup.empty()) {
        std::vector<char> used(static_cast<size_t>(n), 0);
        for (const auto& t : dedup)
            for (const Index vv : t.v) used[size_t(vv)] = 1;
        std::uint64_t usedCount = 0;
        for (const char u : used) usedCount += static_cast<std::uint64_t>(u);
        if (usedCount < st.pointsInserted) {
            st.auditPassed = false;
            throw std::logic_error(
                "delaunay3D: output audit failed (" +
                std::to_string(st.pointsInserted - usedCount) +
                " inserted point(s) not used by any finite tetrahedron)");
        }
    }

    // M4 instrumentation: expose construction memory.
    {
        const auto mem = arena.memory();
        st.arenaPeakLiveCells = mem.peakLiveCells;
        st.arenaMaterializedSlots = mem.materializedSlots;
        // For the record vs the append-only M3-A arena:
        // mem.totalCellsCreated is what it would have retained.
        st.arenaBytes = mem.bytes;
        st.arenaTotalCellsCreated = mem.totalCellsCreated;
    }
    st.runComplete = true;
    return dedup;
}

} // namespace paacs::algo::detail
