// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <future>
#include <limits>
#include <stdexcept>
#include <vector>

#include "paacs/algo/Delaunay2DScalable.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "test_harness.hpp"

namespace {
using paacs::Index;
using paacs::Vec2;
using paacs::algo::Delaunay2DScalableMetrics;
using paacs::algo::Delaunay2DScalableOptions;
using paacs::algo::Tri;
using paacs::predicates::Sign;

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

double unit(std::uint64_t value) {
    return static_cast<double>(value >> 11) *
           (1.0 / 9007199254740992.0);
}

std::vector<Vec2> fixture(std::size_t count) {
    std::vector<Vec2> points;
    points.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const auto ordinal = static_cast<std::uint64_t>(i);
        points.push_back({unit(splitmix64(0x12345678ULL + 2 * ordinal)),
                          unit(splitmix64(0x12345678ULL + 2 * ordinal + 1))});
    }
    return points;
}

std::vector<Vec2> nearDegeneratePrefix(std::size_t count) {
    constexpr std::uint64_t seed = 0xC0FFEEULL + 1000003ULL * 4ULL;
    std::vector<Vec2> points;
    points.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const std::uint64_t ordinal = static_cast<std::uint64_t>(i);
        const double u = unit(splitmix64(seed + 2 * ordinal));
        const double b = unit(splitmix64(seed + 2 * ordinal + 1));
        points.push_back({u, 0.3 + 0.4 * u + (b - 0.5) * 1.0e-11});
    }
    return points;
}

bool retainsEveryPoint(std::size_t pointCount,
                       const std::vector<Tri>& triangles) {
    std::vector<bool> used(pointCount, false);
    for (const Tri& triangle : triangles)
        for (const Index index : triangle.v) {
            if (index < 0 || static_cast<std::size_t>(index) >= pointCount)
                return false;
            used[static_cast<std::size_t>(index)] = true;
        }
    return std::all_of(used.begin(), used.end(), [](bool value) { return value; });
}

std::vector<std::array<Index, 3>> unoriented(std::vector<Tri> triangles) {
    std::vector<std::array<Index, 3>> result;
    result.reserve(triangles.size());
    for (auto& triangle : triangles) {
        std::sort(triangle.v.begin(), triangle.v.end());
        result.push_back(triangle.v);
    }
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<std::array<Index, 3>> geometricLabels(
    const std::vector<Vec2>& reference,
    const std::vector<Vec2>& permuted,
    const std::vector<Tri>& triangles) {
    std::vector<Index> label(permuted.size(), -1);
    for (std::size_t i = 0; i < permuted.size(); ++i)
        for (std::size_t j = 0; j < reference.size(); ++j)
            if (permuted[i].u == reference[j].u &&
                permuted[i].v == reference[j].v) {
                label[i] = static_cast<Index>(j);
                break;
            }
    std::vector<std::array<Index, 3>> result;
    result.reserve(triangles.size());
    for (const auto& triangle : triangles) {
        std::array<Index, 3> ids{
            label[static_cast<std::size_t>(triangle.v[0])],
            label[static_cast<std::size_t>(triangle.v[1])],
            label[static_cast<std::size_t>(triangle.v[2])]};
        std::sort(ids.begin(), ids.end());
        result.push_back(ids);
    }
    std::sort(result.begin(), result.end());
    return result;
}

bool sameMetrics(const Delaunay2DScalableMetrics& lhs,
                 const Delaunay2DScalableMetrics& rhs) {
    return lhs.inputPoints == rhs.inputPoints &&
           lhs.coordinateSortComparisons == rhs.coordinateSortComparisons &&
           lhs.insertionOrderComparisons == rhs.insertionOrderComparisons &&
           lhs.bootstrapPoints == rhs.bootstrapPoints &&
           lhs.hullOrientationTests == rhs.hullOrientationTests &&
           lhs.bootstrapLegalizationEdgeVisits ==
               rhs.bootstrapLegalizationEdgeVisits &&
           lhs.bootstrapLegalizationOrientationTests ==
               rhs.bootstrapLegalizationOrientationTests &&
           lhs.bootstrapLegalizationIncircleTests ==
               rhs.bootstrapLegalizationIncircleTests &&
           lhs.bootstrapCocircularTieTests ==
               rhs.bootstrapCocircularTieTests &&
           lhs.bootstrapLegalizationFlips ==
               rhs.bootstrapLegalizationFlips &&
           lhs.bootstrapWorklistInsertions ==
               rhs.bootstrapWorklistInsertions &&
           lhs.bootstrapWorklistUpdates == rhs.bootstrapWorklistUpdates &&
           lhs.bootstrapStaleEdgeVisits == rhs.bootstrapStaleEdgeVisits &&
           lhs.bootstrapWorklistHighWater == rhs.bootstrapWorklistHighWater &&
           lhs.pointsInserted == rhs.pointsInserted &&
           lhs.locationTriangleTests == rhs.locationTriangleTests &&
           lhs.walkSteps == rhs.walkSteps &&
           lhs.cavityPredicateTests == rhs.cavityPredicateTests &&
           lhs.cavityMembers == rhs.cavityMembers &&
           lhs.boundaryEdges == rhs.boundaryEdges &&
           lhs.trianglesCommitted == rhs.trianglesCommitted &&
           lhs.fullScanFallbacks == rhs.fullScanFallbacks &&
           lhs.fullScanTriangleTests == rhs.fullScanTriangleTests &&
           lhs.peakLiveTriangles == rhs.peakLiveTriangles &&
           lhs.materializedTriangleSlots == rhs.materializedTriangleSlots &&
           lhs.peakCavitySize == rhs.peakCavitySize &&
           lhs.peakBoundarySize == rhs.peakBoundarySize &&
           lhs.scratchHighWaterBytes == rhs.scratchHighWaterBytes &&
           lhs.accountedCapacityHighWaterBytes ==
               rhs.accountedCapacityHighWaterBytes &&
           lhs.adjacencyAuditPassed == rhs.adjacencyAuditPassed &&
           lhs.runComplete == rhs.runComplete;
}

std::uint64_t totalConstructionWork(
    const Delaunay2DScalableMetrics& metrics) {
    return metrics.coordinateSortComparisons +
           metrics.insertionOrderComparisons +
           metrics.hullOrientationTests +
           metrics.bootstrapLegalizationEdgeVisits +
           metrics.bootstrapLegalizationOrientationTests +
           metrics.bootstrapLegalizationIncircleTests +
           metrics.bootstrapLegalizationFlips +
           metrics.bootstrapWorklistInsertions +
           metrics.bootstrapWorklistUpdates +
           metrics.locationTriangleTests +
           metrics.cavityPredicateTests + metrics.cavityMembers +
           metrics.boundaryEdges + metrics.trianglesCommitted +
           metrics.fullScanTriangleTests;
}

std::vector<Vec2> allHullRationalCircle(std::size_t count) {
    std::vector<Vec2> points;
    points.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const double t = count < 2 ? 0.0 :
            -4.0 + 8.0 * static_cast<double>(i) /
                       static_cast<double>(count - 1);
        const double denominator = 1.0 + t * t;
        points.push_back({(1.0 - t * t) / denominator,
                          (2.0 * t) / denominator});
    }
    return points;
}

std::vector<Vec2> hullHeavyMixed(std::size_t count) {
    if (count < 4) return allHullRationalCircle(count);
    const std::size_t hullCount = std::max<std::size_t>(3, (3 * count) / 4);
    std::vector<Vec2> points;
    points.reserve(count);
    for (std::size_t i = 0; i < hullCount; ++i) {
        const double x = static_cast<double>(i);
        points.push_back({x, x * x});
    }
    const double spacing = static_cast<double>(hullCount) /
                           static_cast<double>(count - hullCount);
    for (std::size_t i = hullCount; i < count; ++i) {
        const double x =
            (static_cast<double>(i - hullCount) + 0.5) * spacing;
        points.push_back({x, x * x + 0.5});
    }
    return points;
}

void scalableExactLocalConstruction() {
    TEST("ck1.delaunay2d scalable exact local construction");
    const auto points = fixture(256);
    Delaunay2DScalableMetrics firstMetrics{};
    Delaunay2DScalableMetrics secondMetrics{};
    const auto first = paacs::algo::delaunay2DScalable(points, &firstMetrics);
    const auto second = paacs::algo::delaunay2DScalable(points, &secondMetrics);

    CHECK(!first.empty());
    CHECK(unoriented(first) == unoriented(second));
    CHECK(sameMetrics(firstMetrics, secondMetrics));
    CHECK(firstMetrics.runComplete);
    CHECK(firstMetrics.adjacencyAuditPassed);
    CHECK(firstMetrics.pointsInserted == points.size());
    CHECK(firstMetrics.fullScanFallbacks == 0);
    CHECK(firstMetrics.fullScanTriangleTests == 0);
    CHECK(firstMetrics.bootstrapWorklistUpdates > 0);
    CHECK(firstMetrics.bootstrapStaleEdgeVisits == 0);
    CHECK(firstMetrics.locationTriangleTests < 64 * points.size());
    CHECK(firstMetrics.materializedTriangleSlots < 4 * points.size());

    std::vector<bool> used(points.size(), false);
    std::array<Index, 3> previous{-1, -1, -1};
    bool firstRow = true;
    for (const auto& triangle : first) {
        CHECK(triangle.v[0] < triangle.v[1]);
        CHECK(triangle.v[0] < triangle.v[2]);
        CHECK(paacs::predicates::orient2d(
                  points[static_cast<std::size_t>(triangle.v[0])],
                  points[static_cast<std::size_t>(triangle.v[1])],
                  points[static_cast<std::size_t>(triangle.v[2])]) ==
              Sign::Positive);
        if (!firstRow) CHECK(previous < triangle.v);
        previous = triangle.v;
        firstRow = false;
        for (const Index vertex : triangle.v)
            used[static_cast<std::size_t>(vertex)] = true;
        for (std::size_t point = 0; point < points.size(); ++point) {
            if (point == static_cast<std::size_t>(triangle.v[0]) ||
                point == static_cast<std::size_t>(triangle.v[1]) ||
                point == static_cast<std::size_t>(triangle.v[2]))
                continue;
            CHECK(paacs::predicates::incircle(
                      points[static_cast<std::size_t>(triangle.v[0])],
                      points[static_cast<std::size_t>(triangle.v[1])],
                      points[static_cast<std::size_t>(triangle.v[2])],
                      points[point]) != Sign::Positive);
        }
    }
    CHECK(std::all_of(used.begin(), used.end(), [](bool value) { return value; }));

    const auto legacy = paacs::algo::delaunay2D(points);
    CHECK(unoriented(legacy) == unoriented(first));
}

void scalableOrderingAndZeroPolicy() {
    TEST("ck1.delaunay2d spatial order, original identity and zero policy");
    const auto points = fixture(96);
    auto reversed = points;
    std::reverse(reversed.begin(), reversed.end());
    const auto forwardTriangles = paacs::algo::delaunay2DScalable(points);
    const auto reversedTriangles = paacs::algo::delaunay2DScalable(reversed);
    CHECK(geometricLabels(points, points, forwardTriangles) ==
          geometricLabels(points, reversed, reversedTriangles));

    const std::vector<Vec2> cocircular{{0.0, 0.0}, {1.0, 0.0},
                                       {1.0, 1.0}, {0.0, 1.0}};
    const auto first = paacs::algo::delaunay2DScalable(cocircular);
    const auto second = paacs::algo::delaunay2DScalable(cocircular);
    CHECK(unoriented(first) == unoriented(second));
    CHECK(first.size() == 2);

    const std::vector<Vec2> collinear{{0.0, 0.0}, {0.25, 0.0},
                                      {0.5, 0.0}, {0.75, 0.0}, {1.0, 0.0}};
    Delaunay2DScalableMetrics collinearMetrics{};
    CHECK(paacs::algo::delaunay2DScalable(collinear, &collinearMetrics).empty());
    CHECK(collinearMetrics.runComplete);
    CHECK(collinearMetrics.pointsInserted == collinear.size());

    auto rejected = [](const std::vector<Vec2>& input) {
        try {
            (void)paacs::algo::delaunay2DScalable(input);
        } catch (const std::invalid_argument&) {
            return true;
        } catch (...) {
            return false;
        }
        return false;
    };
    CHECK(rejected({{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}}));
    CHECK(rejected({{0.0, 0.0}, {1.0, 0.0},
                    {(std::numeric_limits<double>::quiet_NaN)(), 1.0}}));
}

void scalableOperationCountScaling() {
    TEST("ck1.delaunay2d operation-count and memory scaling");
    std::array<Delaunay2DScalableMetrics, 3> metrics{};
    constexpr std::array<std::size_t, 3> sizes{256, 512, 1024};
    for (std::size_t i = 0; i < sizes.size(); ++i) {
        const auto points = fixture(sizes[i]);
        const auto triangles =
            paacs::algo::delaunay2DScalable(points, &metrics[i]);
        CHECK(!triangles.empty());
        CHECK(metrics[i].runComplete);
        CHECK(metrics[i].adjacencyAuditPassed);
        CHECK(metrics[i].fullScanFallbacks == 0);
        CHECK(metrics[i].materializedTriangleSlots < 4 * sizes[i]);
        CHECK(metrics[i].scratchHighWaterBytes < 256 * sizes[i]);
        CHECK(metrics[i].accountedCapacityHighWaterBytes < 512 * sizes[i]);
    }
    CHECK(metrics[1].locationTriangleTests <
          3 * metrics[0].locationTriangleTests);
    CHECK(metrics[2].locationTriangleTests <
          3 * metrics[1].locationTriangleTests);
    CHECK(metrics[1].cavityPredicateTests <
          3 * metrics[0].cavityPredicateTests);
    CHECK(metrics[2].cavityPredicateTests <
          3 * metrics[1].cavityPredicateTests);
}

void scalableAdversarialExactnessAndConcurrency() {
    TEST("ck1.delaunay2d adversarial exactness and independent concurrency");
    std::vector<Vec2> narrow;
    narrow.reserve(66);
    for (std::size_t i = 0; i < 33; ++i) {
        const double u = static_cast<double>(i) / 32.0;
        narrow.push_back({u, 0.0});
        narrow.push_back({u, 1.0e-12 + (i % 2 == 0 ? 0.0 : 1.0e-15)});
    }
    Delaunay2DScalableMetrics serialMetrics{};
    const auto serial =
        paacs::algo::delaunay2DScalable(narrow, &serialMetrics);
    CHECK(serialMetrics.runComplete);
    CHECK(serialMetrics.adjacencyAuditPassed);
    CHECK(serialMetrics.fullScanFallbacks == 0);
    CHECK(!serial.empty());

    struct ConcurrentResult {
        std::vector<Tri> output;
        Delaunay2DScalableMetrics metrics{};
    };
    std::array<std::future<ConcurrentResult>, 4> futures;
    for (auto& future : futures)
        future = std::async(std::launch::async, [&narrow] {
            ConcurrentResult result{};
            result.output = paacs::algo::delaunay2DScalable(
                narrow, &result.metrics);
            return result;
        });
    for (auto& future : futures) {
        const ConcurrentResult result = future.get();
        CHECK(unoriented(result.output) == unoriented(serial));
        CHECK(sameMetrics(result.metrics, serialMetrics));
    }

    const double base = 1.0e12;
    const std::vector<Vec2> translated{
        {base, base}, {base + 1.0, base},
        {base + 2.0, base + 0.0001220703125},
        {base, base + 2.0}, {base + 1.0, base + 1.0001220703125},
        {base + 2.0, base + 2.0}};
    const auto translatedTriangles =
        paacs::algo::delaunay2DScalable(translated);
    CHECK(!translatedTriangles.empty());
    for (const Tri& triangle : translatedTriangles)
        CHECK(paacs::predicates::orient2d(
                  translated[static_cast<std::size_t>(triangle.v[0])],
                  translated[static_cast<std::size_t>(triangle.v[1])],
                  translated[static_cast<std::size_t>(triangle.v[2])]) ==
              Sign::Positive);
}

void scalableGuardedFallback() {
    TEST("ck1.delaunay2d guarded full-scan fallback");
    const auto points = fixture(96);
    Delaunay2DScalableMetrics localMetrics{};
    Delaunay2DScalableMetrics fallbackMetrics{};
    const auto local = paacs::algo::delaunay2DScalable(points, &localMetrics);
    const auto fallback = paacs::algo::delaunay2DScalable(
        points, &fallbackMetrics, Delaunay2DScalableOptions{0});
    CHECK(unoriented(local) == unoriented(fallback));
    CHECK(localMetrics.fullScanFallbacks == 0);
    CHECK(fallbackMetrics.fullScanFallbacks +
              fallbackMetrics.bootstrapPoints == points.size());
    CHECK(fallbackMetrics.fullScanTriangleTests > points.size());
    CHECK(fallbackMetrics.runComplete);
    CHECK(fallbackMetrics.adjacencyAuditPassed);
}

void scalableFiniteBootstrapRegressions() {
    TEST("ck1.delaunay2d finite bootstrap retains difficult identities");
    for (const std::size_t count : {256u, 1024u, 4096u}) {
        const auto points = nearDegeneratePrefix(count);
        Delaunay2DScalableMetrics metrics{};
        const auto triangles =
            paacs::algo::delaunay2DScalable(points, &metrics);
        CHECK(!triangles.empty());
        CHECK(retainsEveryPoint(points.size(), triangles));
        CHECK(metrics.pointsInserted == points.size());
        CHECK(metrics.fullScanFallbacks == 0);
        CHECK(metrics.adjacencyAuditPassed);
        CHECK(metrics.runComplete);
    }

    const std::vector<Vec2> hullCollinear{
        {0, 0}, {1, 0}, {2, 0}, {3, 0}, {3, 1}, {3, 2},
        {2, 2}, {1, 2}, {0, 2}, {0, 1}, {1, 1}, {2, 1}};
    const auto hullTriangles = paacs::algo::delaunay2DScalable(hullCollinear);
    CHECK(retainsEveryPoint(hullCollinear.size(), hullTriangles));

    std::vector<Vec2> grid;
    for (int v = 0; v < 8; ++v)
        for (int u = 0; u < 8; ++u)
            grid.push_back({static_cast<double>(u), static_cast<double>(v)});
    const auto gridTriangles = paacs::algo::delaunay2DScalable(grid);
    CHECK(retainsEveryPoint(grid.size(), gridTriangles));

    const std::vector<Vec2> circle{{5, 0}, {3, 4}, {0, 5}, {-3, 4},
                                   {-5, 0}, {-3, -4}, {0, -5}, {3, -4}};
    auto permuted = circle;
    std::rotate(permuted.begin(), permuted.begin() + 3, permuted.end());
    std::reverse(permuted.begin(), permuted.end());
    Delaunay2DScalableMetrics circleMetrics{};
    const auto circleTriangles =
        paacs::algo::delaunay2DScalable(circle, &circleMetrics);
    const auto permutedTriangles = paacs::algo::delaunay2DScalable(permuted);
    CHECK(retainsEveryPoint(circle.size(), circleTriangles));
    CHECK(circleMetrics.bootstrapCocircularTieTests > 0);
    CHECK(geometricLabels(circle, circle, circleTriangles) ==
          geometricLabels(circle, permuted, permutedTriangles));
}

void scalableBootstrapWorklistScaling() {
    TEST("ck1.delaunay2d local bootstrap worklist scaling");
    constexpr std::array<std::size_t, 4> sizes{128, 256, 512, 1024};
    for (const bool mixed : {false, true}) {
        std::array<Delaunay2DScalableMetrics, sizes.size()> metrics{};
        for (std::size_t i = 0; i < sizes.size(); ++i) {
            const auto points = mixed ? hullHeavyMixed(sizes[i])
                                      : allHullRationalCircle(sizes[i]);
            const auto first =
                paacs::algo::delaunay2DScalable(points, &metrics[i]);
            Delaunay2DScalableMetrics repeatedMetrics{};
            const auto repeated =
                paacs::algo::delaunay2DScalable(points, &repeatedMetrics);
            CHECK(unoriented(first) == unoriented(repeated));
            CHECK(sameMetrics(metrics[i], repeatedMetrics));
            CHECK(retainsEveryPoint(points.size(), first));
            CHECK(metrics[i].runComplete);
            CHECK(metrics[i].adjacencyAuditPassed);
            CHECK(metrics[i].fullScanFallbacks == 0);
            CHECK(metrics[i].bootstrapLegalizationEdgeVisits > 0);
            CHECK(metrics[i].bootstrapLegalizationIncircleTests > 0);
            CHECK(metrics[i].bootstrapWorklistHighWater > 0);
            if (mixed) {
                CHECK(metrics[i].bootstrapPoints >= sizes[i] / 2);
                CHECK(metrics[i].bootstrapPoints < sizes[i]);
                CHECK(metrics[i].locationTriangleTests > 0);
                CHECK(metrics[i].cavityPredicateTests > 0);
            } else {
                CHECK(metrics[i].bootstrapPoints == sizes[i]);
                CHECK(metrics[i].locationTriangleTests == 0);
                CHECK(metrics[i].cavityPredicateTests == 0);
                CHECK(metrics[i].bootstrapLegalizationFlips > 0);
            }
            if (i != 0)
                CHECK(totalConstructionWork(metrics[i]) <
                      3 * totalConstructionWork(metrics[i - 1]));
        }
    }

    const auto reference = allHullRationalCircle(128);
    auto permuted = reference;
    std::rotate(permuted.begin(), permuted.begin() + 37, permuted.end());
    std::reverse(permuted.begin(), permuted.end());
    const auto referenceTriangles =
        paacs::algo::delaunay2DScalable(reference);
    const auto permutedTriangles =
        paacs::algo::delaunay2DScalable(permuted);
    CHECK(geometricLabels(reference, reference, referenceTriangles) ==
          geometricLabels(reference, permuted, permutedTriangles));
}
} // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };

const Entry* delaunayScalableTests(std::size_t& count) {
    static const Entry tests[] = {
        {"ck1.delaunay2d_local", scalableExactLocalConstruction},
        {"ck1.delaunay2d_zero_order", scalableOrderingAndZeroPolicy},
        {"ck1.delaunay2d_scaling", scalableOperationCountScaling},
        {"ck1.delaunay2d_fallback", scalableGuardedFallback},
        {"ck1.delaunay2d_adversarial", scalableAdversarialExactnessAndConcurrency},
        {"ck1.delaunay2d_bootstrap", scalableFiniteBootstrapRegressions},
        {"ck1.delaunay2d_bootstrap_worklist", scalableBootstrapWorklistScaling},
    };
    count = sizeof(tests) / sizeof(tests[0]);
    return tests;
}
} // namespace paacs_tests
