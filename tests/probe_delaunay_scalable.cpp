// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// CK-1 scalable-constructor allocation-failure containment probe.
// Run without sanitizers: this executable replaces global operator new.
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <vector>

#include "paacs/algo/Delaunay2DScalable.hpp"

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
void* operator new(std::size_t size) {
    if (gArmed) {
        ++gCalls;
        if (gFailAt != 0 && gCalls == gFailAt) throw std::bad_alloc();
    }
    if (void* memory = std::malloc(size == 0 ? 1 : size)) return memory;
    throw std::bad_alloc();
}
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete(void* memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void* memory, std::size_t) noexcept { std::free(memory); }
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
using paacs::Vec2;
using paacs::algo::Delaunay2DScalableMetrics;
using paacs::algo::Tri;

#define REQUIRE(condition)                                                \
    do {                                                                  \
        if (!(condition)) {                                               \
            std::printf("FAIL %d: %s\n", __LINE__, #condition);          \
            ++gBad;                                                       \
        }                                                                 \
    } while (0)

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9E3779B97F4A7C15ULL;
    value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9ULL;
    value = (value ^ (value >> 27)) * 0x94D049BB133111EBULL;
    return value ^ (value >> 31);
}

double unit(std::uint64_t value) {
    return static_cast<double>(value >> 11) *
           (1.0 / 9007199254740992.0);
}

std::vector<Vec2> randomFixture(std::size_t count) {
    std::vector<Vec2> points;
    points.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const std::uint64_t ordinal = static_cast<std::uint64_t>(i);
        points.push_back({unit(splitmix64(0xCAFE1234ULL + 2 * ordinal)),
                          unit(splitmix64(0xCAFE1234ULL + 2 * ordinal + 1))});
    }
    return points;
}

std::vector<Vec2> allHullFixture(std::size_t count) {
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

bool sameTriangles(const std::vector<Tri>& lhs,
                   const std::vector<Tri>& rhs) {
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i)
        if (lhs[i].v != rhs[i].v) return false;
    return true;
}

struct Attempt {
    bool threw = false;
    std::size_t calls = 0;
    std::vector<Tri> destination{{{{71, 72, 73}}}};
    Delaunay2DScalableMetrics metrics{};
};

Attempt run(const std::vector<Vec2>& points, std::size_t failAt) {
    Attempt result{};
    gCalls = 0;
    gFailAt = failAt;
    gArmed = true;
    try {
        result.destination =
            paacs::algo::delaunay2DScalable(points, &result.metrics);
    } catch (const std::bad_alloc&) {
        result.threw = true;
    } catch (...) {
        gArmed = false;
        std::printf("FAIL %d: unexpected exception escaped fault injection\n",
                    __LINE__);
        ++gBad;
        result.threw = true;
    }
    gArmed = false;
    result.calls = gCalls;
    return result;
}

std::size_t sweep(const std::vector<Vec2>& points) {
    // Exact predicates own process-local initialization. Warm it outside the
    // measured constructor interval so every injected ordinal belongs to the
    // constructor itself and is reproducible on the disarmed retry.
    (void)paacs::algo::delaunay2DScalable(points);
    const Attempt clean = run(points, 0);
    REQUIRE(!clean.threw);
    REQUIRE(clean.metrics.runComplete);
    REQUIRE(clean.metrics.adjacencyAuditPassed);
    const std::size_t allocationCount = clean.calls;

    for (std::size_t fault = 1; fault <= allocationCount; ++fault) {
        const Attempt failed = run(points, fault);
        REQUIRE(failed.threw);
        REQUIRE(failed.calls == fault);
        REQUIRE(!failed.metrics.runComplete);
        REQUIRE(failed.destination.size() == 1);
        REQUIRE(failed.destination[0].v ==
                (std::array<paacs::Index, 3>{71, 72, 73}));

        const Attempt retry = run(points, 0);
        REQUIRE(!retry.threw);
        REQUIRE(retry.metrics.runComplete);
        REQUIRE(sameTriangles(clean.destination, retry.destination));
    }

    const Attempt pastEnd = run(points, allocationCount + 1);
    REQUIRE(!pastEnd.threw);
    REQUIRE(pastEnd.calls == allocationCount);
    REQUIRE(sameTriangles(clean.destination, pastEnd.destination));
    return allocationCount;
}
} // namespace

int main() {
    std::size_t caught = 0;
    caught += sweep({});
    caught += sweep({{0.0, 0.0}, {1.0, 0.0}});
    caught += sweep({{0.0, 0.0}, {1.0, 0.0},
                     {1.0, 1.0}, {0.0, 1.0}});
    caught += sweep(randomFixture(32));
    caught += sweep(randomFixture(96));
    const std::size_t worklistFaults = sweep(allHullFixture(128));
    caught += worklistFaults;
    REQUIRE(worklistFaults > 0);
    REQUIRE(caught > 0);

    if (gBad != 0) {
        std::printf("probe_delaunay_scalable FAILED: %d assertions\n", gBad);
        return 1;
    }
    std::printf("probe_delaunay_scalable PASSED: %zu allocation faults "
                "(%zu all-hull worklist), zero escapes, destination "
                "preserved, clean retries\n",
                caught, worklistFaults);
    return 0;
}
