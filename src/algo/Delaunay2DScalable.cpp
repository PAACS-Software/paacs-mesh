// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Delaunay2DScalable.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "paacs/predicates/Predicates.hpp"

namespace paacs::algo {
namespace {

using predicates::Sign;
constexpr int kNoTriangle = -1;
constexpr std::array<int, 3> kEdgeA{1, 2, 0};
constexpr std::array<int, 3> kEdgeB{2, 0, 1};

struct TriangleSlot {
    std::array<Index, 3> v{kInvalidIndex, kInvalidIndex, kInvalidIndex};
    std::array<int, 3> nb{kNoTriangle, kNoTriangle, kNoTriangle};
    bool alive = false;
};

struct BoundaryEdge {
    Index a = kInvalidIndex;
    Index b = kInvalidIndex;
    int outside = kNoTriangle;
    int outsideEdge = -1;
};

std::uint64_t spreadBits(std::uint32_t value) {
    std::uint64_t x = value & 0x1fffffu;
    x = (x | (x << 16)) & 0x0000ffff0000ffffULL;
    x = (x | (x << 8)) & 0x00ff00ff00ff00ffULL;
    x = (x | (x << 4)) & 0x0f0f0f0f0f0f0f0fULL;
    x = (x | (x << 2)) & 0x3333333333333333ULL;
    x = (x | (x << 1)) & 0x5555555555555555ULL;
    return x;
}

std::uint32_t quantize(double value) {
    const double clamped = std::clamp(value + 0.5, 0.0, 1.0);
    constexpr double scale = static_cast<double>((1u << 21) - 1u);
    return static_cast<std::uint32_t>(clamped * scale);
}

std::uint64_t morton(const Vec2& point) {
    return spreadBits(quantize(point.u)) |
           (spreadBits(quantize(point.v)) << 1);
}

bool samePoint(const Vec2& lhs, const Vec2& rhs) {
    return lhs.u == rhs.u && lhs.v == rhs.v;
}

int edgeIndex(const TriangleSlot& triangle, Index a, Index b) {
    const auto lo = std::min(a, b);
    const auto hi = std::max(a, b);
    for (int edge = 0; edge < 3; ++edge) {
        const Index ea = triangle.v[static_cast<std::size_t>(kEdgeA[edge])];
        const Index eb = triangle.v[static_cast<std::size_t>(kEdgeB[edge])];
        if (std::min(ea, eb) == lo && std::max(ea, eb) == hi) return edge;
    }
    return -1;
}

void canonicalize(Tri& triangle) {
    if (triangle.v[1] < triangle.v[0] && triangle.v[1] < triangle.v[2])
        triangle.v = {triangle.v[1], triangle.v[2], triangle.v[0]};
    else if (triangle.v[2] < triangle.v[0] && triangle.v[2] < triangle.v[1])
        triangle.v = {triangle.v[2], triangle.v[0], triangle.v[1]};
}

} // namespace

std::vector<Tri> delaunay2DScalable(
    const std::vector<Vec2>& points,
    Delaunay2DScalableMetrics* metricsOut,
    Delaunay2DScalableOptions options) {
    Delaunay2DScalableMetrics localMetrics{};
    Delaunay2DScalableMetrics& metrics =
        metricsOut == nullptr ? localMetrics : *metricsOut;
    metrics.reset();
    metrics.inputPoints = points.size();

    for (std::size_t i = 0; i < points.size(); ++i)
        if (!std::isfinite(points[i].u) || !std::isfinite(points[i].v))
            throw std::invalid_argument(
                "delaunay2DScalable: non-finite input at ordinal " +
                std::to_string(i));
    constexpr std::size_t kMaxDensePoints =
        (static_cast<std::size_t>((std::numeric_limits<int>::max)()) - 8) / 2;
    if (points.size() > kMaxDensePoints)
        throw std::length_error("delaunay2DScalable: dense index limit exceeded");

    std::vector<Index> duplicateOrder(points.size());
    for (std::size_t i = 0; i < points.size(); ++i)
        duplicateOrder[i] = static_cast<Index>(i);
    std::sort(duplicateOrder.begin(), duplicateOrder.end(),
              [&](Index lhs, Index rhs) {
                  ++metrics.coordinateSortComparisons;
                  const Vec2& a = points[static_cast<std::size_t>(lhs)];
                  const Vec2& b = points[static_cast<std::size_t>(rhs)];
                  return std::tie(a.u, a.v, lhs) < std::tie(b.u, b.v, rhs);
              });
    for (std::size_t i = 1; i < duplicateOrder.size(); ++i) {
        const Vec2& previous =
            points[static_cast<std::size_t>(duplicateOrder[i - 1])];
        const Vec2& current =
            points[static_cast<std::size_t>(duplicateOrder[i])];
        if (samePoint(previous, current))
            throw std::invalid_argument("delaunay2DScalable: duplicate input");
    }
    if (points.size() < 3) {
        metrics.pointsInserted = points.size();
        metrics.adjacencyAuditPassed = true;
        metrics.runComplete = true;
        return {};
    }

    const Index n = static_cast<Index>(points.size());
    long double loU = points[0].u, hiU = points[0].u;
    long double loV = points[0].v, hiV = points[0].v;
    for (const Vec2& point : points) {
        loU = std::min(loU, static_cast<long double>(point.u));
        hiU = std::max(hiU, static_cast<long double>(point.u));
        loV = std::min(loV, static_cast<long double>(point.v));
        hiV = std::max(hiV, static_cast<long double>(point.v));
    }
    const long double extent = std::max(hiU - loU, hiV - loV);
    const long double scale = extent > 0.0L ? extent : 1.0L;
    const long double midpointU = (loU + hiU) * 0.5L;
    const long double midpointV = (loV + hiV) * 0.5L;

    std::vector<Vec2> spatial;
    spatial.reserve(points.size());
    for (const Vec2& point : points)
        spatial.push_back(
            {static_cast<double>((static_cast<long double>(point.u) -
                                  midpointU) /
                                 scale),
             static_cast<double>((static_cast<long double>(point.v) -
                                  midpointV) /
                                 scale)});

    std::vector<Vec2> working = points;

    std::vector<Index> insertionOrder(points.size());
    for (Index i = 0; i < n; ++i) insertionOrder[static_cast<std::size_t>(i)] = i;
    std::sort(insertionOrder.begin(), insertionOrder.end(),
              [&](Index lhs, Index rhs) {
                  ++metrics.insertionOrderComparisons;
                  const Vec2& a = spatial[static_cast<std::size_t>(lhs)];
                  const Vec2& b = spatial[static_cast<std::size_t>(rhs)];
                  const Vec2& originalA = points[static_cast<std::size_t>(lhs)];
                  const Vec2& originalB = points[static_cast<std::size_t>(rhs)];
                  return std::tuple{morton(a), originalA.u, originalA.v, lhs} <
                         std::tuple{morton(b), originalB.u, originalB.v, rhs};
              });

    // Exact finite-point bootstrap: the strict convex hull supplies the only
    // initial vertices. Collinear hull points are inserted later on the
    // corresponding boundary edge, so every declared identity is retained.
    std::vector<Index> lower;
    std::vector<Index> upper;
    lower.reserve(points.size());
    upper.reserve(points.size());
    auto hullOrientation = [&](Index a, Index b, Index c) {
        ++metrics.hullOrientationTests;
        return predicates::orient2d(
            working[static_cast<std::size_t>(a)],
            working[static_cast<std::size_t>(b)],
            working[static_cast<std::size_t>(c)]);
    };
    for (const Index index : duplicateOrder) {
        while (lower.size() >= 2 &&
               hullOrientation(lower[lower.size() - 2], lower.back(), index) !=
                   Sign::Positive)
            lower.pop_back();
        lower.push_back(index);
    }
    for (auto it = duplicateOrder.rbegin(); it != duplicateOrder.rend(); ++it) {
        const Index index = *it;
        while (upper.size() >= 2 &&
               hullOrientation(upper[upper.size() - 2], upper.back(), index) !=
                   Sign::Positive)
            upper.pop_back();
        upper.push_back(index);
    }
    std::vector<Index> hull;
    hull.reserve(lower.size() + upper.size());
    for (std::size_t i = 0; i + 1 < lower.size(); ++i) hull.push_back(lower[i]);
    for (std::size_t i = 0; i + 1 < upper.size(); ++i) hull.push_back(upper[i]);
    if (hull.size() < 3) {
        metrics.pointsInserted = points.size();
        metrics.adjacencyAuditPassed = true;
        metrics.runComplete = true;
        return {};
    }

    std::vector<bool> bootstrap(points.size(), false);
    for (const Index index : hull)
        bootstrap[static_cast<std::size_t>(index)] = true;

    std::vector<TriangleSlot> triangles;
    triangles.reserve(points.size() * 2 + 8);
    for (std::size_t i = 1; i + 1 < hull.size(); ++i) {
        TriangleSlot triangle{{hull[0], hull[i], hull[i + 1]},
                              {kNoTriangle, kNoTriangle, kNoTriangle}, true};
        if (i + 2 < hull.size()) triangle.nb[1] = static_cast<int>(i);
        if (i > 1) triangle.nb[2] = static_cast<int>(i - 2);
        triangles.push_back(triangle);
    }

    metrics.bootstrapPoints = hull.size();
    auto pointKey = [&](Index index) {
        const Vec2& point = working[static_cast<std::size_t>(index)];
        return std::tuple{point.u, point.v, index};
    };
    auto geometricEdgeLess = [&](Index a, Index b, Index c, Index d) {
        if (pointKey(b) < pointKey(a)) std::swap(a, b);
        if (pointKey(d) < pointKey(c)) std::swap(c, d);
        return std::tuple{pointKey(a), pointKey(b)} <
               std::tuple{pointKey(c), pointKey(d)};
    };
    auto replaceOutsideOwner = [&](int outside, Index a, Index b,
                                   int oldOwner, int newOwner) {
        if (outside < 0) return;
        TriangleSlot& triangle = triangles[static_cast<std::size_t>(outside)];
        const int edge = edgeIndex(triangle, a, b);
        if (edge < 0 || triangle.nb[static_cast<std::size_t>(edge)] != oldOwner)
            throw std::logic_error(
                "delaunay2DScalable: bootstrap reciprocal edge missing");
        triangle.nb[static_cast<std::size_t>(edge)] = newOwner;
    };
    auto neighborAcross = [](const TriangleSlot& triangle, Index a, Index b) {
        const int edge = edgeIndex(triangle, a, b);
        if (edge < 0)
            throw std::logic_error(
                "delaunay2DScalable: bootstrap boundary edge missing");
        return triangle.nb[static_cast<std::size_t>(edge)];
    };
    // Lawson-legalize the finite hull fan through a deterministic local edge
    // worklist. Geometry-keyed map membership both orders work independently
    // of caller indices and updates an already-pending edge instead of
    // duplicating it. Exact-zero diagonals monotonically select the
    // geometry-lexicographically smaller undirected diagonal.
    using BootstrapEdgeKey = std::pair<Index, Index>;
    auto canonicalBootstrapEdge = [&](Index a, Index b) {
        if (pointKey(b) < pointKey(a)) std::swap(a, b);
        return BootstrapEdgeKey{a, b};
    };
    auto bootstrapEdgeLess = [&](const BootstrapEdgeKey& lhs,
                                 const BootstrapEdgeKey& rhs) {
        return std::tuple{pointKey(lhs.first), pointKey(lhs.second)} <
               std::tuple{pointKey(rhs.first), pointKey(rhs.second)};
    };
    std::map<BootstrapEdgeKey, std::array<int, 2>,
             decltype(bootstrapEdgeLess)>
        bootstrapWorklist(bootstrapEdgeLess);
    auto enqueueBootstrapEdge = [&](int owner, Index a, Index b) {
        if (owner < 0) return;
        const TriangleSlot& triangle =
            triangles[static_cast<std::size_t>(owner)];
        const int edge = edgeIndex(triangle, a, b);
        if (edge < 0)
            throw std::logic_error(
                "delaunay2DScalable: bootstrap enqueue edge missing");
        const int neighbor = triangle.nb[static_cast<std::size_t>(edge)];
        if (neighbor < 0) return;
        const BootstrapEdgeKey key = canonicalBootstrapEdge(a, b);
        const std::array<int, 2> owners{std::min(owner, neighbor),
                                        std::max(owner, neighbor)};
        const auto [it, inserted] = bootstrapWorklist.emplace(key, owners);
        if (inserted) {
            ++metrics.bootstrapWorklistInsertions;
        } else {
            it->second = owners;
            ++metrics.bootstrapWorklistUpdates;
        }
        metrics.bootstrapWorklistHighWater = std::max<std::uint64_t>(
            metrics.bootstrapWorklistHighWater, bootstrapWorklist.size());
    };
    for (std::size_t slot = 0; slot < triangles.size(); ++slot)
        for (int edge = 0; edge < 3; ++edge) {
            const int neighbor =
                triangles[slot].nb[static_cast<std::size_t>(edge)];
            if (neighbor > static_cast<int>(slot))
                enqueueBootstrapEdge(
                    static_cast<int>(slot),
                    triangles[slot].v[static_cast<std::size_t>(kEdgeA[edge])],
                    triangles[slot].v[static_cast<std::size_t>(kEdgeB[edge])]);
        }

    while (!bootstrapWorklist.empty()) {
        const auto item = *bootstrapWorklist.begin();
        bootstrapWorklist.erase(bootstrapWorklist.begin());
        ++metrics.bootstrapLegalizationEdgeVisits;
        const Index keyA = item.first.first;
        const Index keyB = item.first.second;
        const int firstSlot = item.second[0];
        const int secondSlot = item.second[1];
        if (firstSlot < 0 || secondSlot < 0 ||
            static_cast<std::size_t>(firstSlot) >= triangles.size() ||
            static_cast<std::size_t>(secondSlot) >= triangles.size()) {
            ++metrics.bootstrapStaleEdgeVisits;
            continue;
        }
        const TriangleSlot oldFirst =
            triangles[static_cast<std::size_t>(firstSlot)];
        const TriangleSlot oldSecond =
            triangles[static_cast<std::size_t>(secondSlot)];
        const int edge = edgeIndex(oldFirst, keyA, keyB);
        const int secondEdge = edgeIndex(oldSecond, keyA, keyB);
        if (!oldFirst.alive || !oldSecond.alive || edge < 0 ||
            secondEdge < 0 ||
            oldFirst.nb[static_cast<std::size_t>(edge)] != secondSlot ||
            oldSecond.nb[static_cast<std::size_t>(secondEdge)] != firstSlot) {
            ++metrics.bootstrapStaleEdgeVisits;
            continue;
        }

        const Index a =
            oldFirst.v[static_cast<std::size_t>(kEdgeA[edge])];
        const Index b =
            oldFirst.v[static_cast<std::size_t>(kEdgeB[edge])];
        const Index c = oldFirst.v[static_cast<std::size_t>(edge)];
        const Index d = oldSecond.v[static_cast<std::size_t>(secondEdge)];
        ++metrics.bootstrapLegalizationOrientationTests;
        if (predicates::orient2d(
                working[static_cast<std::size_t>(c)],
                working[static_cast<std::size_t>(a)],
                working[static_cast<std::size_t>(d)]) != Sign::Positive)
            continue;
        ++metrics.bootstrapLegalizationOrientationTests;
        if (predicates::orient2d(
                working[static_cast<std::size_t>(c)],
                working[static_cast<std::size_t>(d)],
                working[static_cast<std::size_t>(b)]) != Sign::Positive)
            continue;
        ++metrics.bootstrapLegalizationIncircleTests;
        const Sign circle = predicates::incircle(
            working[static_cast<std::size_t>(oldFirst.v[0])],
            working[static_cast<std::size_t>(oldFirst.v[1])],
            working[static_cast<std::size_t>(oldFirst.v[2])],
            working[static_cast<std::size_t>(d)]);
        if (circle == Sign::Zero) ++metrics.bootstrapCocircularTieTests;
        if (circle == Sign::Negative ||
            (circle == Sign::Zero && !geometricEdgeLess(c, d, a, b)))
            continue;

        const int tBC = neighborAcross(oldFirst, b, c);
        const int tCA = neighborAcross(oldFirst, c, a);
        const int uAD = neighborAcross(oldSecond, a, d);
        const int uDB = neighborAcross(oldSecond, d, b);
        TriangleSlot first{{c, a, d},
                           {kNoTriangle, kNoTriangle, kNoTriangle}, true};
        TriangleSlot second{{c, d, b},
                            {kNoTriangle, kNoTriangle, kNoTriangle}, true};
        first.nb[static_cast<std::size_t>(edgeIndex(first, c, a))] = tCA;
        first.nb[static_cast<std::size_t>(edgeIndex(first, a, d))] = uAD;
        first.nb[static_cast<std::size_t>(edgeIndex(first, c, d))] = secondSlot;
        second.nb[static_cast<std::size_t>(edgeIndex(second, c, d))] = firstSlot;
        second.nb[static_cast<std::size_t>(edgeIndex(second, d, b))] = uDB;
        second.nb[static_cast<std::size_t>(edgeIndex(second, b, c))] = tBC;
        triangles[static_cast<std::size_t>(firstSlot)] = first;
        triangles[static_cast<std::size_t>(secondSlot)] = second;
        replaceOutsideOwner(tCA, c, a, firstSlot, firstSlot);
        replaceOutsideOwner(uAD, a, d, secondSlot, firstSlot);
        replaceOutsideOwner(uDB, d, b, secondSlot, secondSlot);
        replaceOutsideOwner(tBC, b, c, firstSlot, secondSlot);
        ++metrics.bootstrapLegalizationFlips;

        for (const int owner : {firstSlot, secondSlot}) {
            const TriangleSlot& changed =
                triangles[static_cast<std::size_t>(owner)];
            for (int changedEdge = 0; changedEdge < 3; ++changedEdge)
                enqueueBootstrapEdge(
                    owner,
                    changed.v[static_cast<std::size_t>(kEdgeA[changedEdge])],
                    changed.v[static_cast<std::size_t>(kEdgeB[changedEdge])]);
        }
    }
    std::vector<int> freeSlots;
    freeSlots.reserve(points.size());
    std::vector<std::uint32_t> visited(triangles.size(), 0);
    std::vector<std::uint32_t> inCavity(triangles.size(), 0);
    std::uint32_t epoch = 0;
    std::uint64_t liveTriangles = triangles.size();
    int lastLocated = 0;
    metrics.pointsInserted = hull.size();
    metrics.trianglesCommitted = triangles.size();

    std::vector<int> stack;
    std::vector<int> cavity;
    std::vector<BoundaryEdge> boundary;
    std::vector<int> boundaryOrder;
    std::vector<int> newSlots;
    std::vector<int> edgeByStart(points.size() + 3, -1);
    std::vector<std::uint32_t> edgeStamp(points.size() + 3, 0);
    stack.reserve(64);
    cavity.reserve(64);
    boundary.reserve(64);
    boundaryOrder.reserve(64);
    newSlots.reserve(64);

    auto updateMemory = [&] {
        const std::uint64_t scratchBytes =
            stack.capacity() * sizeof(int) + cavity.capacity() * sizeof(int) +
            boundary.capacity() * sizeof(BoundaryEdge) +
            boundaryOrder.capacity() * sizeof(int) +
            newSlots.capacity() * sizeof(int) +
            edgeByStart.capacity() * sizeof(int) +
            edgeStamp.capacity() * sizeof(std::uint32_t) +
            visited.capacity() * sizeof(std::uint32_t) +
            inCavity.capacity() * sizeof(std::uint32_t) +
            freeSlots.capacity() * sizeof(int);
        metrics.scratchHighWaterBytes =
            std::max(metrics.scratchHighWaterBytes, scratchBytes);
        const std::uint64_t persistentBytes =
            spatial.capacity() * sizeof(Vec2) +
            working.capacity() * sizeof(Vec2) +
            duplicateOrder.capacity() * sizeof(Index) +
            insertionOrder.capacity() * sizeof(Index) +
            triangles.capacity() * sizeof(TriangleSlot);
        metrics.accountedCapacityHighWaterBytes = std::max(
            metrics.accountedCapacityHighWaterBytes,
            scratchBytes + persistentBytes);
        metrics.peakLiveTriangles =
            std::max(metrics.peakLiveTriangles, liveTriangles);
        metrics.materializedTriangleSlots =
            std::max<std::uint64_t>(metrics.materializedTriangleSlots,
                                    triangles.size());
    };

    auto addTriangle = [&](const TriangleSlot& triangle) {
        int slot = kNoTriangle;
        if (!freeSlots.empty()) {
            slot = freeSlots.back();
            freeSlots.pop_back();
            triangles[static_cast<std::size_t>(slot)] = triangle;
        } else {
            if (triangles.size() >=
                static_cast<std::size_t>((std::numeric_limits<int>::max)()))
                throw std::length_error(
                    "delaunay2DScalable: triangle slot limit exceeded");
            slot = static_cast<int>(triangles.size());
            triangles.push_back(triangle);
            visited.push_back(0);
            inCavity.push_back(0);
        }
        ++liveTriangles;
        return slot;
    };

    auto conflicts = [&](int slot, const Vec2& point) {
        ++metrics.cavityPredicateTests;
        const TriangleSlot& triangle = triangles[static_cast<std::size_t>(slot)];
        return predicates::incircle(
                   working[static_cast<std::size_t>(triangle.v[0])],
                   working[static_cast<std::size_t>(triangle.v[1])],
                   working[static_cast<std::size_t>(triangle.v[2])], point) !=
               Sign::Negative;
    };

    updateMemory();
    for (const Index pointIndex : insertionOrder) {
        if (bootstrap[static_cast<std::size_t>(pointIndex)]) continue;
        const Vec2& point = working[static_cast<std::size_t>(pointIndex)];
        if (++epoch == 0) {
            std::fill(visited.begin(), visited.end(), 0);
            std::fill(inCavity.begin(), inCavity.end(), 0);
            std::fill(edgeStamp.begin(), edgeStamp.end(), 0);
            epoch = 1;
        }

        int current = lastLocated;
        if (current < 0 ||
            !triangles[static_cast<std::size_t>(current)].alive) {
            current = kNoTriangle;
            for (std::size_t i = 0; i < triangles.size(); ++i)
                if (triangles[i].alive) {
                    current = static_cast<int>(i);
                    break;
                }
        }

        int located = kNoTriangle;
        const std::uint64_t topologyWalkGuard = 4 * liveTriangles + 32;
        const std::uint64_t walkGuard =
            std::min(topologyWalkGuard, options.maxWalkSteps);
        for (std::uint64_t step = 0;
             current >= 0 && step < walkGuard; ++step) {
            const TriangleSlot& triangle =
                triangles[static_cast<std::size_t>(current)];
            ++metrics.locationTriangleTests;
            int exitEdge = -1;
            for (int edge = 0; edge < 3; ++edge) {
                const Vec2& a = working[static_cast<std::size_t>(
                    triangle.v[static_cast<std::size_t>(kEdgeA[edge])])];
                const Vec2& b = working[static_cast<std::size_t>(
                    triangle.v[static_cast<std::size_t>(kEdgeB[edge])])];
                if (predicates::orient2d(a, b, point) == Sign::Negative) {
                    exitEdge = edge;
                    break;
                }
            }
            if (exitEdge < 0) {
                located = current;
                break;
            }
            current = triangle.nb[static_cast<std::size_t>(exitEdge)];
            ++metrics.walkSteps;
        }

        int seed = kNoTriangle;
        if (located >= 0 && conflicts(located, point)) seed = located;
        if (seed < 0 && located >= 0) {
            const TriangleSlot& triangle =
                triangles[static_cast<std::size_t>(located)];
            for (const int neighbor : triangle.nb)
                if (neighbor >= 0 &&
                    triangles[static_cast<std::size_t>(neighbor)].alive &&
                    conflicts(neighbor, point)) {
                    seed = neighbor;
                    break;
                }
        }
        if (seed < 0) {
            ++metrics.fullScanFallbacks;
            for (std::size_t slot = 0; slot < triangles.size(); ++slot) {
                if (!triangles[slot].alive) continue;
                ++metrics.fullScanTriangleTests;
                if (conflicts(static_cast<int>(slot), point)) {
                    seed = static_cast<int>(slot);
                    break;
                }
            }
        }
        if (seed < 0)
            throw std::logic_error(
                "delaunay2DScalable: no conflict seed after guarded fallback");

        stack.clear();
        cavity.clear();
        stack.push_back(seed);
        visited[static_cast<std::size_t>(seed)] = epoch;
        inCavity[static_cast<std::size_t>(seed)] = epoch;
        while (!stack.empty()) {
            const int slot = stack.back();
            stack.pop_back();
            cavity.push_back(slot);
            const TriangleSlot& triangle =
                triangles[static_cast<std::size_t>(slot)];
            for (const int neighbor : triangle.nb) {
                if (neighbor < 0 ||
                    visited[static_cast<std::size_t>(neighbor)] == epoch)
                    continue;
                visited[static_cast<std::size_t>(neighbor)] = epoch;
                if (!triangles[static_cast<std::size_t>(neighbor)].alive ||
                    !conflicts(neighbor, point))
                    continue;
                inCavity[static_cast<std::size_t>(neighbor)] = epoch;
                stack.push_back(neighbor);
            }
        }
        metrics.cavityMembers += cavity.size();
        metrics.peakCavitySize =
            std::max<std::uint64_t>(metrics.peakCavitySize, cavity.size());

        boundary.clear();
        for (const int slot : cavity) {
            const TriangleSlot& triangle =
                triangles[static_cast<std::size_t>(slot)];
            for (int edge = 0; edge < 3; ++edge) {
                const int outside = triangle.nb[static_cast<std::size_t>(edge)];
                if (outside >= 0 &&
                    inCavity[static_cast<std::size_t>(outside)] == epoch)
                    continue;
                const Index a = triangle.v[static_cast<std::size_t>(kEdgeA[edge])];
                const Index b = triangle.v[static_cast<std::size_t>(kEdgeB[edge])];
                const int outsideEdge = outside < 0 ? -1 :
                    edgeIndex(triangles[static_cast<std::size_t>(outside)], a, b);
                if (outside >= 0 && outsideEdge < 0)
                    throw std::logic_error(
                        "delaunay2DScalable: reciprocal boundary edge missing");
                boundary.push_back({a, b, outside, outsideEdge});
            }
        }
        if (boundary.empty())
            throw std::logic_error("delaunay2DScalable: empty cavity boundary");
        metrics.boundaryEdges += boundary.size();
        metrics.peakBoundarySize =
            std::max<std::uint64_t>(metrics.peakBoundarySize, boundary.size());

        int start = 0;
        for (std::size_t i = 0; i < boundary.size(); ++i) {
            const std::size_t vertex = static_cast<std::size_t>(boundary[i].a);
            if (edgeStamp[vertex] == epoch)
                throw std::logic_error(
                    "delaunay2DScalable: cavity boundary is not a simple cycle");
            edgeStamp[vertex] = epoch;
            edgeByStart[vertex] = static_cast<int>(i);
            if (std::tie(boundary[i].a, boundary[i].b) <
                std::tie(boundary[static_cast<std::size_t>(start)].a,
                         boundary[static_cast<std::size_t>(start)].b))
                start = static_cast<int>(i);
        }
        boundaryOrder.clear();
        int edge = start;
        for (std::size_t count = 0; count < boundary.size(); ++count) {
            boundaryOrder.push_back(edge);
            const Index nextVertex = boundary[static_cast<std::size_t>(edge)].b;
            if (nextVertex < 0 ||
                edgeStamp[static_cast<std::size_t>(nextVertex)] != epoch)
                throw std::logic_error(
                    "delaunay2DScalable: open cavity boundary");
            edge = edgeByStart[static_cast<std::size_t>(nextVertex)];
        }
        if (edge != start)
            throw std::logic_error(
                "delaunay2DScalable: multiple cavity boundary cycles");

        for (const int boundaryIndex : boundaryOrder) {
            const BoundaryEdge& item =
                boundary[static_cast<std::size_t>(boundaryIndex)];
            const Sign orientation = predicates::orient2d(
                working[static_cast<std::size_t>(item.a)],
                working[static_cast<std::size_t>(item.b)], point);
            if (orientation == Sign::Negative)
                throw std::logic_error(
                    "delaunay2DScalable: non-star-shaped conflict cavity");
            if (orientation == Sign::Zero && item.outside >= 0)
                throw std::logic_error(
                    "delaunay2DScalable: degenerate interior boundary edge");
        }

        for (const int slot : cavity) {
            TriangleSlot& triangle = triangles[static_cast<std::size_t>(slot)];
            triangle.alive = false;
            triangle.nb = {kNoTriangle, kNoTriangle, kNoTriangle};
            freeSlots.push_back(slot);
            --liveTriangles;
        }

        newSlots.assign(boundaryOrder.size(), kNoTriangle);
        for (std::size_t i = 0; i < boundaryOrder.size(); ++i) {
            const BoundaryEdge& item = boundary[static_cast<std::size_t>(
                boundaryOrder[i])];
            if (predicates::orient2d(
                    working[static_cast<std::size_t>(item.a)],
                    working[static_cast<std::size_t>(item.b)], point) ==
                Sign::Zero)
                continue;
            TriangleSlot replacement{{item.a, item.b, pointIndex},
                                     {kNoTriangle, kNoTriangle, item.outside},
                                     true};
            newSlots[i] = addTriangle(replacement);
            ++metrics.trianglesCommitted;
        }

        for (std::size_t i = 0; i < boundaryOrder.size(); ++i) {
            const int slot = newSlots[i];
            if (slot < 0) continue;
            const std::size_t previous =
                (i + boundaryOrder.size() - 1) % boundaryOrder.size();
            const std::size_t next = (i + 1) % boundaryOrder.size();
            if (newSlots[previous] >= 0)
                triangles[static_cast<std::size_t>(slot)].nb[1] =
                    newSlots[previous];
            if (newSlots[next] >= 0)
                triangles[static_cast<std::size_t>(slot)].nb[0] = newSlots[next];

            const BoundaryEdge& item = boundary[static_cast<std::size_t>(
                boundaryOrder[i])];
            if (item.outside >= 0)
                triangles[static_cast<std::size_t>(item.outside)]
                    .nb[static_cast<std::size_t>(item.outsideEdge)] = slot;
        }

        lastLocated = kNoTriangle;
        for (const int slot : newSlots)
            if (slot >= 0) {
                lastLocated = slot;
                break;
            }
        if (lastLocated < 0)
            throw std::logic_error(
                "delaunay2DScalable: insertion created no live triangle");
        ++metrics.pointsInserted;
        updateMemory();
    }

    bool adjacencyValid = true;
    for (std::size_t slot = 0; slot < triangles.size(); ++slot) {
        const TriangleSlot& triangle = triangles[slot];
        if (!triangle.alive) continue;
        for (int edge = 0; edge < 3; ++edge) {
            const int neighbor = triangle.nb[static_cast<std::size_t>(edge)];
            if (neighbor < 0) continue;
            if (static_cast<std::size_t>(neighbor) >= triangles.size() ||
                !triangles[static_cast<std::size_t>(neighbor)].alive) {
                adjacencyValid = false;
                continue;
            }
            const Index a = triangle.v[static_cast<std::size_t>(kEdgeA[edge])];
            const Index b = triangle.v[static_cast<std::size_t>(kEdgeB[edge])];
            const int back = edgeIndex(
                triangles[static_cast<std::size_t>(neighbor)], a, b);
            if (back < 0 || triangles[static_cast<std::size_t>(neighbor)]
                                    .nb[static_cast<std::size_t>(back)] !=
                                static_cast<int>(slot))
                adjacencyValid = false;
        }
    }
    metrics.adjacencyAuditPassed = adjacencyValid;
    if (!adjacencyValid)
        throw std::logic_error("delaunay2DScalable: adjacency audit failed");

    std::vector<Tri> output;
    output.reserve(static_cast<std::size_t>(liveTriangles));
    for (const TriangleSlot& triangle : triangles) {
        if (!triangle.alive || triangle.v[0] >= n || triangle.v[1] >= n ||
            triangle.v[2] >= n)
            continue;
        Tri item{triangle.v};
        canonicalize(item);
        output.push_back(item);
    }
    std::sort(output.begin(), output.end(),
              [](const Tri& lhs, const Tri& rhs) { return lhs.v < rhs.v; });
    metrics.runComplete = true;
    return output;
}

} // namespace paacs::algo
