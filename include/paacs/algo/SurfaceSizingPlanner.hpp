// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/mesh/Mesh.hpp"
#include "paacs/sizing/SizeField.hpp"

#include <cstdint>
#include <functional>
#include <limits>
#include <vector>

namespace paacs::algo {

struct SurfaceSizingPlacement {
    int faceTag = -1;
    Index firstNode = kInvalidIndex;
    Index secondNode = kInvalidIndex;
    Vec3 point{};
    double edgeLength = 0.0;
    double targetSize = 0.0;
    double violationRatio = 0.0;
};
using SurfaceProjectionCallback = std::function<bool(
    int faceTag, const Vec3& candidate, Vec3& projected)>;

struct SurfaceSizingPlanRequest {
    const sizing::SizeField* field = nullptr;
    double maximumEdgeToTargetRatio = 1.35;
    std::uint64_t maximumPlacements = 1'000'000;
    int executionThreads = 0;
    SurfaceProjectionCallback project;
};
enum class SurfaceSizingPlanReason : std::uint8_t {
    None, EmptySurface, MissingField, UnsupportedCell,
    InvalidConnectivity, InvalidClassification, InvalidTargetSize,
    ProjectionFailed, BudgetExceeded, AllocationFailed, InternalError
};
struct SurfaceSizingPlanMetrics {
    std::uint64_t sourceEdges = 0;
    std::uint64_t violatingEdges = 0;
    std::uint64_t projectedPlacements = 0;
    double maximumViolationRatio = 0.0;
};
struct SurfaceSizingPlanOutcome {
    bool planned = false;
    SurfaceSizingPlanReason reason = SurfaceSizingPlanReason::InternalError;
    std::vector<SurfaceSizingPlacement> placements;
    SurfaceSizingPlanMetrics metrics{};
};

// One deterministic incremental placement pass. The caller can feed the
// result into the existing persistent constrained-surface queue, update only
// the affected cavities, and invoke the planner again if further refinement is
// required. No whole-surface reconstruction is performed here.
SurfaceSizingPlanOutcome planSurfaceSizingPass(
    const mesh::Mesh& surface,
    const SurfaceSizingPlanRequest& request) noexcept;

}  // namespace paacs::algo
