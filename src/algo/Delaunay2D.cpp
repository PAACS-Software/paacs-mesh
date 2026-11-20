// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Delaunay2D.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace paacs::algo {

static double orient2d(const Vec2& a, const Vec2& b, const Vec2& c) {
    return (b.u - a.u) * (c.v - a.v) - (b.v - a.v) * (c.u - a.u);
}

bool inCircumcircle(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& p) {
    // Standard lifted determinant; assumes (a,b,c) CCW.
    const double ax = a.u - p.u, ay = a.v - p.v;
    const double bx = b.u - p.u, by = b.v - p.v;
    const double cx = c.u - p.u, cy = c.v - p.v;
    const double det =
        (ax * ax + ay * ay) * (bx * cy - cx * by) -
        (bx * bx + by * by) * (ax * cy - cx * ay) +
        (cx * cx + cy * cy) * (ax * by - bx * ay);
    return det > 1e-12;
}

#if defined(_MSC_VER)
#define PAACS_CK0_FORCEINLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define PAACS_CK0_FORCEINLINE inline __attribute__((always_inline))
#else
#define PAACS_CK0_FORCEINLINE inline
#endif

template <bool Measure>
PAACS_CK0_FORCEINLINE std::vector<Tri> delaunay2DImpl(
    const std::vector<Vec2>& pts,
    ConstructionMetrics* metrics,
    std::span<InsertionRecord> perInsertion) {
    std::vector<Tri> out;
    const Index n = static_cast<Index>(pts.size());
    if (n < 3) return out;

    if constexpr (Measure)
        metrics->detailComplete = perInsertion.size() >= static_cast<std::size_t>(n);

    // Normalized working frame (PM-NUM-001/002): predicates run in a
    // unit-scale frame so absolute thresholds and the super-triangle are
    // conditioned identically for models of any size or position. The old
    // "span + 1.0" super-triangle was absolute-scale and destroyed tiny
    // models. Connectivity is index-based, so output is frame-independent.
    double lo_u = 1e300, lo_v = 1e300, hi_u = -1e300, hi_v = -1e300;
    for (const auto& p : pts) {
        lo_u = std::min(lo_u, p.u); hi_u = std::max(hi_u, p.u);
        lo_v = std::min(lo_v, p.v); hi_v = std::max(hi_v, p.v);
    }
    const double ext = std::max(hi_u - lo_u, hi_v - lo_v);
    const double scale = ext > 0 ? ext : 1.0;
    const double invS = 1.0 / scale;
    const Vec2 midW{(lo_u + hi_u) / 2, (lo_v + hi_v) / 2};

    std::vector<Vec2> P(pts.size());
    for (size_t i = 0; i < pts.size(); ++i)
        P[i] = {(pts[i].u - midW.u) * invS, (pts[i].v - midW.v) * invS};
    const double span = 1.0;
    const Vec2 mid{0, 0};
    P.push_back({mid.u - 20 * span, mid.v - 10 * span});   // n
    P.push_back({mid.u + 20 * span, mid.v - 10 * span});   // n+1
    P.push_back({mid.u,             mid.v + 20 * span});   // n+2

    std::vector<Tri> tris{{{n, n + 1, n + 2}}};

    if constexpr (Measure) metrics->peakLiveTriangles = 1;

    for (Index ip = 0; ip < n; ++ip) {
        const Vec2& p = P[static_cast<size_t>(ip)];

        InsertionRecord record{};
        if constexpr (Measure) {
            record.inputOrdinal = static_cast<std::uint64_t>(ip);
            record.usedFullScan = true;
        }

        // Cavity: all triangles whose circumcircle contains p.
        std::vector<Tri> bad, keep;
        for (const auto& t : tris) {
            const Vec2 &a = P[size_t(t.v[0])], &b = P[size_t(t.v[1])],
                       &c = P[size_t(t.v[2])];
            if constexpr (Measure) {
                ++metrics->trianglesScanned;
                ++record.trianglesScanned;
                const bool inCavity = inCircumcircle(a, b, c, p);
                if (inCavity) {
                    ++metrics->cavityMembers;
                    ++record.cavityMembers;
                    bad.push_back(t);
                } else keep.push_back(t);
            } else {
                (inCircumcircle(a, b, c, p) ? bad : keep).push_back(t);
            }
        }

        // Cavity boundary = edges appearing exactly once among bad triangles.
        std::map<std::pair<Index, Index>, int> edgeCount;
        auto ekey = [](Index a, Index b) {
            return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
        };
        for (const auto& t : bad)
            for (int e = 0; e < 3; ++e) {
                if constexpr (Measure) {
                    const auto key = ekey(t.v[e], t.v[(e + 1) % 3]);
                    const bool inserted = edgeCount.find(key) == edgeCount.end();
                    edgeCount[key]++;
                    if (inserted) {
                        ++metrics->boundaryMapInsertions;
                        ++record.boundaryMapInsertions;
                    }
                } else edgeCount[ekey(t.v[e], t.v[(e + 1) % 3])]++;
            }

        if constexpr (Measure) {
            const auto mapSize = static_cast<std::uint32_t>(std::min<std::size_t>(
                edgeCount.size(), (std::numeric_limits<std::uint32_t>::max)()));
            metrics->peakBoundaryMapEntries =
                std::max(metrics->peakBoundaryMapEntries, mapSize);
        }

        tris = std::move(keep);
        for (const auto& t : bad)
            for (int e = 0; e < 3; ++e) {
                Index a = t.v[e], b = t.v[(e + 1) % 3];
                if (edgeCount[ekey(a, b)] != 1) continue;
                Tri nt{{a, b, ip}};
                // Enforce CCW.
                if (orient2d(P[size_t(a)], P[size_t(b)], p) < 0)
                    std::swap(nt.v[0], nt.v[1]);
                tris.push_back(nt);
                if constexpr (Measure) {
                    ++metrics->boundaryEdges;
                    ++record.boundaryEdgesEmitted;
                    ++metrics->cellsCommitted;
                    ++record.cellsCommitted;
                }
            }

        if constexpr (Measure) {
            ++metrics->insertions;
            ++metrics->fullScanInsertions;
            const auto live = static_cast<std::uint32_t>(std::min<std::size_t>(
                tris.size(), (std::numeric_limits<std::uint32_t>::max)()));
            record.liveTrianglesAfterCommit = live;
            metrics->peakLiveTriangles = std::max(metrics->peakLiveTriangles, live);
            const auto cavity = static_cast<std::uint32_t>(std::min<std::size_t>(
                bad.size(), (std::numeric_limits<std::uint32_t>::max)()));
            metrics->peakCavitySize = std::max(metrics->peakCavitySize, cavity);
            ++metrics->cavityHistogram[std::min<std::size_t>(bad.size(), 15)];
            ++metrics->walkHistogram[0];
            if (static_cast<std::size_t>(ip) < perInsertion.size()) {
                perInsertion[static_cast<std::size_t>(ip)] = record;
                ++metrics->rowsWritten;
            } else {
                ++metrics->rowsDropped;
                metrics->detailComplete = false;
            }
        }
    }

    // Drop triangles touching the super-triangle.
    for (const auto& t : tris)
        if (t.v[0] < n && t.v[1] < n && t.v[2] < n) out.push_back(t);
    return out;
}

std::vector<Tri> delaunay2D(const std::vector<Vec2>& pts) {
    return delaunay2DImpl<false>(pts, nullptr, {});
}

std::vector<Tri> delaunay2DMeasured(
    const std::vector<Vec2>& pts,
    ConstructionMetrics& metrics,
    std::span<InsertionRecord> perInsertion) {
    metrics.reset();
    std::vector<Tri> out =
        delaunay2DImpl<true>(pts, &metrics, perInsertion);
    metrics.runComplete = true;
    return out;
}

#undef PAACS_CK0_FORCEINLINE

} // namespace paacs::algo
