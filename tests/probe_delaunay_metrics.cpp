// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// CK-0 measured/unmeasured allocation-ordinal equivalence probe.
// Run without sanitizers: this executable replaces global operator new.
#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <vector>

#include "paacs/algo/Delaunay2D.hpp"

namespace {
bool gArmed = false;
std::size_t gCalls = 0;
std::size_t gFailAt = 0;
int gBad = 0;
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t n) {
    if (gArmed) {
        ++gCalls;
        if (gFailAt != 0 && gCalls == gFailAt) throw std::bad_alloc();
    }
    if (void* p = std::malloc(n == 0 ? 1 : n)) return p;
    throw std::bad_alloc();
}
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
using paacs::Vec2;
using paacs::algo::ConstructionMetrics;
using paacs::algo::InsertionRecord;
using paacs::algo::Tri;

#define REQUIRE(c)                                                        \
    do {                                                                  \
        if (!(c)) {                                                       \
            std::printf("FAIL %d: %s\n", __LINE__, #c);                  \
            ++gBad;                                                       \
        }                                                                 \
    } while (0)

bool sameTriangles(const std::vector<Tri>& a, const std::vector<Tri>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i].v != b[i].v) return false;
    return true;
}

struct Attempt {
    bool threw = false;
    std::size_t calls = 0;
    std::vector<Tri> output;
    ConstructionMetrics metrics{};
};

Attempt ordinaryAttempt(const std::vector<Vec2>& points, std::size_t failAt) {
    Attempt result{};
    gCalls = 0;
    gFailAt = failAt;
    gArmed = true;
    try {
        result.output = paacs::algo::delaunay2D(points);
    } catch (const std::bad_alloc&) {
        result.threw = true;
    }
    gArmed = false;
    result.calls = gCalls;
    return result;
}

Attempt measuredAttempt(const std::vector<Vec2>& points, std::size_t failAt) {
    Attempt result{};
    std::array<InsertionRecord, 32> rows{};
    gCalls = 0;
    gFailAt = failAt;
    gArmed = true;
    try {
        result.output =
            paacs::algo::delaunay2DMeasured(points, result.metrics, rows);
    } catch (const std::bad_alloc&) {
        result.threw = true;
    }
    gArmed = false;
    result.calls = gCalls;
    return result;
}

std::vector<Vec2> fixture() {
    return {{0.00, 0.00}, {1.00, 0.00}, {1.00, 1.00}, {0.00, 1.00},
            {0.15, 0.20}, {0.35, 0.12}, {0.62, 0.18}, {0.84, 0.33},
            {0.78, 0.67}, {0.61, 0.86}, {0.37, 0.82}, {0.16, 0.61},
            {0.31, 0.42}, {0.52, 0.36}, {0.68, 0.51}, {0.49, 0.64}};
}

std::vector<Vec2> nearDegenerateFixture() {
    std::vector<Vec2> points;
    points.reserve(16);
    for (std::size_t i = 0; i < 16; ++i) {
        const double u = static_cast<double>(i) / 15.0;
        const double noise = (i % 2 == 0 ? 1.0 : -1.0) * 1e-11;
        points.push_back({u, 0.3 + 0.4 * u + noise});
    }
    return points;
}

std::size_t sweepFixture(const std::vector<Vec2>& points) {
    const Attempt ordinaryClean = ordinaryAttempt(points, 0);
    const Attempt measuredClean = measuredAttempt(points, 0);
    REQUIRE(!ordinaryClean.threw);
    REQUIRE(!measuredClean.threw);
    REQUIRE(ordinaryClean.calls == measuredClean.calls);
    REQUIRE(sameTriangles(ordinaryClean.output, measuredClean.output));
    REQUIRE(measuredClean.metrics.runComplete);
    REQUIRE(measuredClean.metrics.detailComplete);

    const std::size_t allocationCount = ordinaryClean.calls;
    std::size_t caught = 0;
    for (std::size_t fault = 1; fault <= allocationCount; ++fault) {
        const Attempt ordinary = ordinaryAttempt(points, fault);
        const Attempt measured = measuredAttempt(points, fault);
        REQUIRE(ordinary.threw);
        REQUIRE(measured.threw);
        REQUIRE(ordinary.calls == fault);
        REQUIRE(measured.calls == fault);
        REQUIRE(!measured.metrics.runComplete);
        ++caught;
    }

    const Attempt ordinaryPastEnd = ordinaryAttempt(points, allocationCount + 1);
    const Attempt measuredPastEnd = measuredAttempt(points, allocationCount + 1);
    REQUIRE(!ordinaryPastEnd.threw);
    REQUIRE(!measuredPastEnd.threw);
    REQUIRE(ordinaryPastEnd.calls == allocationCount);
    REQUIRE(measuredPastEnd.calls == allocationCount);
    REQUIRE(sameTriangles(ordinaryPastEnd.output, measuredPastEnd.output));
    REQUIRE(measuredPastEnd.metrics.runComplete);

    const Attempt ordinaryRetry = ordinaryAttempt(points, 0);
    const Attempt measuredRetry = measuredAttempt(points, 0);
    REQUIRE(!ordinaryRetry.threw);
    REQUIRE(!measuredRetry.threw);
    REQUIRE(sameTriangles(ordinaryClean.output, ordinaryRetry.output));
    REQUIRE(sameTriangles(ordinaryRetry.output, measuredRetry.output));
    REQUIRE(ordinaryRetry.calls == measuredRetry.calls);
    return caught;
}
} // namespace

int main() {
    std::size_t caught = 0;
    caught += sweepFixture({});
    caught += sweepFixture({{0.0, 0.0}, {1.0, 0.0}});
    caught += sweepFixture({{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}});
    caught += sweepFixture(fixture());
    caught += sweepFixture(nearDegenerateFixture());
    REQUIRE(caught > 0);

    if (gBad != 0) {
        std::printf("probe_delaunay_metrics FAILED: %d assertions\n", gBad);
        return 1;
    }
    std::printf("probe_delaunay_metrics PASSED: %zu allocation faults, "
                "zero ordinal differences, clean retry\n", caught);
    return 0;
}
