// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// CK-0 approved measured-entry functional and semantic-equivalence proofs.
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <thread>
#include <vector>

#include "paacs/algo/Delaunay2D.hpp"
#include "test_harness.hpp"

namespace {
using paacs::Vec2;
using paacs::algo::ConstructionMetrics;
using paacs::algo::InsertionRecord;
using paacs::algo::Tri;

bool sameTriangles(const std::vector<Tri>& a, const std::vector<Tri>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i].v != b[i].v) return false;
    return true;
}

bool sameMetrics(const ConstructionMetrics& lhs,
                 const ConstructionMetrics& rhs) {
    if (lhs.insertions != rhs.insertions ||
        lhs.trianglesScanned != rhs.trianglesScanned ||
        lhs.cavityMembers != rhs.cavityMembers ||
        lhs.boundaryEdges != rhs.boundaryEdges ||
        lhs.boundaryMapInsertions != rhs.boundaryMapInsertions ||
        lhs.boundaryMapErasures != rhs.boundaryMapErasures ||
        lhs.cellsCommitted != rhs.cellsCommitted ||
        lhs.fullScanInsertions != rhs.fullScanInsertions ||
        lhs.rowsWritten != rhs.rowsWritten ||
        lhs.rowsDropped != rhs.rowsDropped ||
        lhs.peakLiveTriangles != rhs.peakLiveTriangles ||
        lhs.peakBoundaryMapEntries != rhs.peakBoundaryMapEntries ||
        lhs.peakCavitySize != rhs.peakCavitySize ||
        lhs.runComplete != rhs.runComplete ||
        lhs.detailComplete != rhs.detailComplete)
        return false;
    for (std::size_t i = 0; i < 16; ++i)
        if (lhs.cavityHistogram[i] != rhs.cavityHistogram[i] ||
            lhs.walkHistogram[i] != rhs.walkHistogram[i])
            return false;
    return true;
}

bool sameInsertionRecord(const InsertionRecord& lhs,
                         const InsertionRecord& rhs) {
    return lhs.inputOrdinal == rhs.inputOrdinal &&
           lhs.trianglesScanned == rhs.trianglesScanned &&
           lhs.cavityMembers == rhs.cavityMembers &&
           lhs.boundaryMapInsertions == rhs.boundaryMapInsertions &&
           lhs.boundaryMapErasures == rhs.boundaryMapErasures &&
           lhs.boundaryEdgesEmitted == rhs.boundaryEdgesEmitted &&
           lhs.cellsCommitted == rhs.cellsCommitted &&
           lhs.liveTrianglesAfterCommit == rhs.liveTrianglesAfterCommit &&
           lhs.walkSteps == rhs.walkSteps &&
           lhs.usedFullScan == rhs.usedFullScan;
}

enum class Family : std::uint64_t {
    Uniform,
    Jittered,
    Clustered,
    StronglyNonuniform,
    NearDegenerate
};

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

double u01(std::uint64_t h) {
    return static_cast<double>(h >> 11) * (1.0 / 9007199254740992.0);
}

std::uint32_t bitreverse(std::uint32_t x, int bits) {
    std::uint32_t result = 0;
    for (int i = 0; i < bits; ++i) {
        result = (result << 1) | (x & 1u);
        x >>= 1;
    }
    return result;
}

std::uint32_t deinterleave(std::uint32_t code) {
    std::uint32_t result = 0;
    for (int i = 0; i < 12; ++i)
        result |= ((code >> (2 * i)) & 1u) << i;
    return result;
}

std::vector<Vec2> familyPrefix(Family family, std::size_t count) {
    const std::uint64_t seed =
        0xC0FFEEULL + 1000003ULL * static_cast<std::uint64_t>(family);
    std::vector<Vec2> points;
    points.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const std::uint64_t ordinal = static_cast<std::uint64_t>(i);
        const double a = u01(splitmix64(seed + 2 * ordinal));
        const double b = u01(splitmix64(seed + 2 * ordinal + 1));
        switch (family) {
            case Family::Uniform:
                points.push_back({a, b});
                break;
            case Family::Jittered: {
                const std::uint32_t code =
                    bitreverse(static_cast<std::uint32_t>(i), 24);
                const std::uint32_t x = deinterleave(code);
                const std::uint32_t y = deinterleave(code >> 1);
                constexpr double range = 4096.0;
                points.push_back({(x + 0.5 + 0.5 * (a - 0.5)) / range,
                                  (y + 0.5 + 0.5 * (b - 0.5)) / range});
                break;
            }
            case Family::Clustered: {
                const int cluster =
                    static_cast<int>(splitmix64(seed + ordinal) % 8);
                const double x = 0.125 + 0.75 * (cluster % 4) / 3.0;
                const double y = 0.2 + 0.6 * (cluster / 4);
                points.push_back(
                    {x + 0.02 * (a - 0.5), y + 0.02 * (b - 0.5)});
                break;
            }
            case Family::StronglyNonuniform:
                points.push_back({a * a * a, b * b * b});
                break;
            case Family::NearDegenerate:
                points.push_back({a, 0.3 + 0.4 * a + (b - 0.5) * 1e-11});
                break;
        }
    }
    return points;
}

std::uint64_t semanticHash(std::vector<Tri> triangles) {
    for (auto& triangle : triangles)
        std::sort(triangle.v.begin(), triangle.v.end());
    std::sort(triangles.begin(), triangles.end(),
              [](const Tri& lhs, const Tri& rhs) { return lhs.v < rhs.v; });
    std::uint64_t hash = 0xCBF29CE484222325ULL;
    const auto put = [&hash](std::uint64_t value) {
        for (int byte = 0; byte < 8; ++byte) {
            hash ^= (value >> (8 * byte)) & 0xFFu;
            hash *= 0x100000001B3ULL;
        }
    };
    put(triangles.size());
    for (const auto& triangle : triangles)
        for (const auto vertex : triangle.v)
            put(static_cast<std::uint64_t>(vertex));
    return hash;
}

void checkAllFamilies() {
    constexpr std::array families{
        Family::Uniform, Family::Jittered, Family::Clustered,
        Family::StronglyNonuniform, Family::NearDegenerate};
    constexpr std::array<std::size_t, 2> sizes{16, 64};
    for (const Family family : families) {
        for (const std::size_t size : sizes) {
            const auto points = familyPrefix(family, size);
            const auto ordinary = paacs::algo::delaunay2D(points);
            std::vector<InsertionRecord> rows(size);
            ConstructionMetrics firstMetrics{};
            ConstructionMetrics secondMetrics{};
            const auto first =
                paacs::algo::delaunay2DMeasured(points, firstMetrics, rows);
            const auto second =
                paacs::algo::delaunay2DMeasured(points, secondMetrics);
            CHECK(sameTriangles(ordinary, first));
            CHECK(sameTriangles(first, second));
            CHECK(semanticHash(ordinary) == semanticHash(first));
            CHECK(firstMetrics.runComplete);
            CHECK(firstMetrics.detailComplete);
            CHECK(secondMetrics.runComplete);
            CHECK(!secondMetrics.detailComplete);
            CHECK(secondMetrics.rowsWritten == 0);
            CHECK(secondMetrics.rowsDropped == size);
            CHECK(firstMetrics.trianglesScanned == secondMetrics.trianglesScanned);
            CHECK(firstMetrics.cavityMembers == secondMetrics.cavityMembers);
            CHECK(firstMetrics.cellsCommitted == secondMetrics.cellsCommitted);
        }
    }
}

std::vector<Vec2> fixture() {
    return {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0},
            {0.2, 0.3}, {0.8, 0.25}, {0.65, 0.8}, {0.35, 0.7}};
}

void metricsEquivalentAndAggregated() {
    TEST("delaunay.metrics equivalent output and complete aggregates");
    const auto points = fixture();
    const auto ordinary = paacs::algo::delaunay2D(points);
    std::array<InsertionRecord, 8> detail{};
    ConstructionMetrics metrics{};
    const auto measured = paacs::algo::delaunay2DMeasured(points, metrics, detail);

    CHECK(sameTriangles(ordinary, measured));
    CHECK(metrics.runComplete);
    CHECK(metrics.detailComplete);
    CHECK(metrics.insertions == points.size());
    CHECK(metrics.fullScanInsertions == points.size());
    CHECK(metrics.rowsWritten == points.size());
    CHECK(metrics.rowsDropped == 0);
    CHECK(metrics.trianglesScanned > 0);
    CHECK(metrics.cavityMembers > 0);
    CHECK(metrics.boundaryMapInsertions > 0);
    CHECK(metrics.boundaryMapErasures == 0);
    CHECK(metrics.boundaryEdges > 0);
    CHECK(metrics.cellsCommitted == metrics.boundaryEdges);
    CHECK(metrics.peakLiveTriangles > 0);
    CHECK(metrics.peakBoundaryMapEntries > 0);
    CHECK(metrics.peakCavitySize > 0);

    std::uint64_t scans = 0, cavities = 0, insertions = 0;
    std::uint64_t erasures = 0, edges = 0, committed = 0;
    for (std::size_t i = 0; i < detail.size(); ++i) {
        CHECK(detail[i].inputOrdinal == i);
        CHECK(detail[i].usedFullScan);
        CHECK(detail[i].walkSteps == 0);
        CHECK(detail[i].liveTrianglesAfterCommit > 0);
        scans += detail[i].trianglesScanned;
        cavities += detail[i].cavityMembers;
        insertions += detail[i].boundaryMapInsertions;
        erasures += detail[i].boundaryMapErasures;
        edges += detail[i].boundaryEdgesEmitted;
        committed += detail[i].cellsCommitted;
    }
    CHECK(scans == metrics.trianglesScanned);
    CHECK(cavities == metrics.cavityMembers);
    CHECK(insertions == metrics.boundaryMapInsertions);
    CHECK(erasures == metrics.boundaryMapErasures);
    CHECK(edges == metrics.boundaryEdges);
    CHECK(committed == metrics.cellsCommitted);

    std::uint64_t cavityHistogramTotal = 0, walkHistogramTotal = 0;
    for (std::size_t i = 0; i < 16; ++i) {
        cavityHistogramTotal += metrics.cavityHistogram[i];
        walkHistogramTotal += metrics.walkHistogram[i];
    }
    CHECK(cavityHistogramTotal == points.size());
    CHECK(walkHistogramTotal == points.size());
    CHECK(metrics.walkHistogram[0] == points.size());
    checkAllFamilies();
}

void metricsTruncationAndReset() {
    TEST("delaunay.metrics run/detail completeness and reset semantics");
    const auto points = fixture();
    std::array<InsertionRecord, 3> detail{};
    ConstructionMetrics metrics{};
    metrics.insertions = 999;
    metrics.runComplete = true;
    metrics.detailComplete = false;
    const auto ordinary = paacs::algo::delaunay2D(points);
    const auto measured = paacs::algo::delaunay2DMeasured(points, metrics, detail);
    CHECK(sameTriangles(ordinary, measured));
    CHECK(metrics.runComplete);
    CHECK(!metrics.detailComplete);
    CHECK(metrics.rowsWritten == detail.size());
    CHECK(metrics.rowsDropped == points.size() - detail.size());
    CHECK(metrics.insertions == points.size());
    for (std::size_t i = 0; i < detail.size(); ++i)
        CHECK(detail[i].inputOrdinal == i);

    metrics.insertions = 123;
    metrics.rowsDropped = 456;
    metrics.runComplete = false;
    metrics.detailComplete = false;
    const std::vector<Vec2> tooSmall{{0.0, 0.0}, {1.0, 0.0}};
    const auto empty = paacs::algo::delaunay2DMeasured(tooSmall, metrics);
    CHECK(empty.empty());
    CHECK(metrics.runComplete);
    CHECK(metrics.detailComplete);
    CHECK(metrics.insertions == 0);
    CHECK(metrics.rowsWritten == 0);
    CHECK(metrics.rowsDropped == 0);
}

void metricsConcurrentIndependentCalls() {
    TEST("delaunay.metrics concurrent independent calls match serial references");
    const auto firstPoints = familyPrefix(Family::Uniform, 128);
    const auto secondPoints = familyPrefix(Family::NearDegenerate, 128);

    std::vector<InsertionRecord> firstSerialRows(firstPoints.size());
    std::vector<InsertionRecord> secondSerialRows(secondPoints.size());
    ConstructionMetrics firstSerialMetrics{};
    ConstructionMetrics secondSerialMetrics{};
    const auto firstSerial = paacs::algo::delaunay2DMeasured(
        firstPoints, firstSerialMetrics, firstSerialRows);
    const auto secondSerial = paacs::algo::delaunay2DMeasured(
        secondPoints, secondSerialMetrics, secondSerialRows);

    std::vector<InsertionRecord> firstConcurrentRows(firstPoints.size());
    std::vector<InsertionRecord> secondConcurrentRows(secondPoints.size());
    ConstructionMetrics firstConcurrentMetrics{};
    ConstructionMetrics secondConcurrentMetrics{};
    std::vector<Tri> firstConcurrent;
    std::vector<Tri> secondConcurrent;

    std::thread firstThread([&] {
        firstConcurrent = paacs::algo::delaunay2DMeasured(
            firstPoints, firstConcurrentMetrics, firstConcurrentRows);
    });
    std::thread secondThread([&] {
        secondConcurrent = paacs::algo::delaunay2DMeasured(
            secondPoints, secondConcurrentMetrics, secondConcurrentRows);
    });
    firstThread.join();
    secondThread.join();

    CHECK(sameTriangles(firstSerial, firstConcurrent));
    CHECK(sameTriangles(secondSerial, secondConcurrent));
    CHECK(sameMetrics(firstSerialMetrics, firstConcurrentMetrics));
    CHECK(sameMetrics(secondSerialMetrics, secondConcurrentMetrics));
    CHECK(firstSerialRows.size() == firstConcurrentRows.size());
    CHECK(secondSerialRows.size() == secondConcurrentRows.size());
    for (std::size_t i = 0; i < firstSerialRows.size(); ++i)
        CHECK(sameInsertionRecord(firstSerialRows[i], firstConcurrentRows[i]));
    for (std::size_t i = 0; i < secondSerialRows.size(); ++i)
        CHECK(sameInsertionRecord(secondSerialRows[i],
                                  secondConcurrentRows[i]));
}
} // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };

const Entry* delaunayMetricTests(std::size_t& count) {
    static const Entry tests[] = {
        {"delaunay.metrics_equivalence", metricsEquivalentAndAggregated},
        {"delaunay.metrics_completeness", metricsTruncationAndReset},
        {"delaunay.metrics_concurrency", metricsConcurrentIndependentCalls},
    };
    count = sizeof(tests) / sizeof(tests[0]);
    return tests;
}
} // namespace paacs_tests
