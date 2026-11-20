// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// CK-1 bounded scalable 2D candidate constructor.
//
// This is additive: legacy delaunay2D() and its accepted goldens remain the
// comparison oracle. The transient kernel uses dense local triangle slots,
// reciprocal neighbours, deterministic Morton insertion, neighbour walking,
// connected conflict-cavity traversal and reusable generation-marked scratch.
// Every topology-affecting sign routes through paacs::predicates.
//
// Exact zero policy (part of this entry's contract):
//   * incircle Zero is included in the conflict cavity (closed-cavity rule);
//   * during finite-point bootstrap legalization, incircle Zero selects the
//     geometry-lexicographically smaller undirected diagonal;
//   * orient2d Zero is never committed as a triangle;
//   * exact duplicate input coordinates are rejected before construction.
// Spatial normalization affects insertion order only; predicates always see
// the caller's original coordinates.
// This policy is deterministic but does not change the legacy constructor's
// accepted degeneracy behavior or goldens.

#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

#include "paacs/algo/Delaunay2D.hpp"

namespace paacs::algo {

struct Delaunay2DScalableMetrics {
    std::uint64_t inputPoints = 0;
    std::uint64_t coordinateSortComparisons = 0;
    std::uint64_t insertionOrderComparisons = 0;
    std::uint64_t bootstrapPoints = 0;
    std::uint64_t hullOrientationTests = 0;
    std::uint64_t bootstrapLegalizationEdgeVisits = 0;
    std::uint64_t bootstrapLegalizationOrientationTests = 0;
    std::uint64_t bootstrapLegalizationIncircleTests = 0;
    std::uint64_t bootstrapCocircularTieTests = 0;
    std::uint64_t bootstrapLegalizationFlips = 0;
    std::uint64_t bootstrapWorklistInsertions = 0;
    std::uint64_t bootstrapWorklistUpdates = 0;
    std::uint64_t bootstrapStaleEdgeVisits = 0;
    std::uint64_t bootstrapWorklistHighWater = 0;
    std::uint64_t pointsInserted = 0;
    std::uint64_t locationTriangleTests = 0;
    std::uint64_t walkSteps = 0;
    std::uint64_t cavityPredicateTests = 0;
    std::uint64_t cavityMembers = 0;
    std::uint64_t boundaryEdges = 0;
    std::uint64_t trianglesCommitted = 0;
    std::uint64_t fullScanFallbacks = 0;
    std::uint64_t fullScanTriangleTests = 0;
    std::uint64_t peakLiveTriangles = 0;
    std::uint64_t materializedTriangleSlots = 0;
    std::uint64_t peakCavitySize = 0;
    std::uint64_t peakBoundarySize = 0;
    std::uint64_t scratchHighWaterBytes = 0;
    // Sum of retained vector capacities at observation points. This is not a
    // process peak-live measurement and excludes allocator reallocation
    // overlap and the returned output buffer.
    std::uint64_t accountedCapacityHighWaterBytes = 0;
    bool adjacencyAuditPassed = false;
    bool runComplete = false;

    void reset() noexcept { *this = Delaunay2DScalableMetrics{}; }
};

struct Delaunay2DScalableOptions {
    // Upper bound for one neighbour walk. The default retains the internal
    // topology-sized guard. Zero deliberately bypasses walking and exercises
    // the counted full-complex fallback without changing accepted output.
    std::uint64_t maxWalkSteps = (std::numeric_limits<std::uint64_t>::max)();
};

static_assert(std::is_standard_layout_v<Delaunay2DScalableMetrics>);
static_assert(std::is_trivially_copyable_v<Delaunay2DScalableMetrics>);
static_assert(std::is_nothrow_copy_constructible_v<Delaunay2DScalableMetrics>);
static_assert(std::is_nothrow_move_constructible_v<Delaunay2DScalableMetrics>);
static_assert(std::is_standard_layout_v<Delaunay2DScalableOptions>);
static_assert(std::is_trivially_copyable_v<Delaunay2DScalableOptions>);

// Throws std::invalid_argument for non-finite or duplicate input,
// std::length_error when dense 32-bit local indexing cannot represent the
// request, and propagates allocation failure. The caller's point vector is
// never mutated. When metrics is non-null it is reset on entry and remains a
// partial diagnostic snapshot if construction throws; runComplete is true
// only when the returned triangulation and final audit are complete.
std::vector<Tri> delaunay2DScalable(
    const std::vector<Vec2>& points,
    Delaunay2DScalableMetrics* metrics = nullptr,
    Delaunay2DScalableOptions options = {});

} // namespace paacs::algo
