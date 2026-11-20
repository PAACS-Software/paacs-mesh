// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// Incremental Bowyer-Watson Delaunay triangulation in the plane.
// Prototype-grade predicates (double precision with tolerance); production
// will swap in Shewchuk's adaptive exact predicates behind the same API.
//
#include <array>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>
#include "paacs/core/Types.hpp"

namespace paacs::algo {

struct Tri { std::array<Index, 3> v; };

// CK-0 approved fixed-storage instrumentation. These are observations only;
// the production entry point below never receives or reads them.
struct InsertionRecord {
    std::uint64_t inputOrdinal = 0;
    std::uint64_t trianglesScanned = 0;
    std::uint64_t cavityMembers = 0;
    std::uint64_t boundaryMapInsertions = 0;
    std::uint64_t boundaryMapErasures = 0;
    std::uint64_t boundaryEdgesEmitted = 0;
    std::uint64_t cellsCommitted = 0;
    std::uint32_t liveTrianglesAfterCommit = 0;
    std::uint32_t walkSteps = 0;
    bool usedFullScan = false;
};

struct ConstructionMetrics {
    std::uint64_t insertions = 0;
    std::uint64_t trianglesScanned = 0; // conflict/cavity membership scans
    std::uint64_t cavityMembers = 0;
    std::uint64_t boundaryEdges = 0;
    std::uint64_t boundaryMapInsertions = 0;
    std::uint64_t boundaryMapErasures = 0;
    std::uint64_t cellsCommitted = 0;
    std::uint64_t fullScanInsertions = 0;
    std::uint64_t rowsWritten = 0;
    std::uint64_t rowsDropped = 0;
    std::uint32_t peakLiveTriangles = 0;
    std::uint32_t peakBoundaryMapEntries = 0;
    std::uint32_t peakCavitySize = 0;
    std::uint32_t cavityHistogram[16] = {};
    std::uint32_t walkHistogram[16] = {};
    bool runComplete = false;
    bool detailComplete = true;

    void reset() noexcept { *this = ConstructionMetrics{}; }
};

static_assert(std::is_standard_layout_v<InsertionRecord>);
static_assert(std::is_trivially_copyable_v<InsertionRecord>);
static_assert(std::is_nothrow_copy_constructible_v<InsertionRecord>);
static_assert(std::is_nothrow_copy_assignable_v<InsertionRecord>);
static_assert(std::is_nothrow_move_constructible_v<InsertionRecord>);
static_assert(std::is_nothrow_move_assignable_v<InsertionRecord>);
static_assert(std::is_standard_layout_v<ConstructionMetrics>);
static_assert(std::is_trivially_copyable_v<ConstructionMetrics>);
static_assert(std::is_nothrow_copy_constructible_v<ConstructionMetrics>);
static_assert(std::is_nothrow_copy_assignable_v<ConstructionMetrics>);
static_assert(std::is_nothrow_move_constructible_v<ConstructionMetrics>);
static_assert(std::is_nothrow_move_assignable_v<ConstructionMetrics>);

// Triangulates the given points. Returns triangles indexing into `pts`.
std::vector<Tri> delaunay2D(const std::vector<Vec2>& pts);

// Additive measured entry point. Instrumentation introduces no allocation;
// ordinary constructor/result allocations and exception behavior are retained.
std::vector<Tri> delaunay2DMeasured(
    const std::vector<Vec2>& pts,
    ConstructionMetrics& metrics,
    std::span<InsertionRecord> perInsertion = {});

// True if p lies strictly inside the circumcircle of (a,b,c) (CCW).
bool inCircumcircle(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& p);

} // namespace paacs::algo
