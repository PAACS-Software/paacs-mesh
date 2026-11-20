// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/PreparedSweepMesher.hpp"

#include <cstdint>
#include <limits>

namespace paacs::algo {

struct CurvedPrismLayerRequest {
    std::uint32_t layers = 5;
    double firstLayerThickness = 1.0e-3;
    double cellToCellGrowth = 1.2;
    double maximumAllowedGrowth = 1.3;
    double maximumTotalThickness = 1.0;
    double minimumTotalThickness = 1.0e-8;
    double featureAngleDegrees = 45.0;
    std::uint32_t normalSmoothingIterations = 4;
    std::uint32_t thicknessSmoothingIterations = 4;
    double maximumAdjacentThicknessRatio = 1.5;
    double collisionSafetyFactor = 0.75;
    std::uint32_t maximumGlobalBackoffSteps = 6;
    double globalBackoffFactor = 0.75;
    const mesh::Mesh* collisionSurface = nullptr;
    mesh::Classification volumeClassification{3,1};
    mesh::Classification targetClassification{2,2};
    mesh::Classification sideClassification{2,3};
    bool emitBoundaryFaces = true;
    int executionThreads = 0;
    std::uint64_t maximumCollisionTriangles = 10'000'000;
};

enum class CurvedPrismReason : std::uint8_t {
    None, EmptySource, UnsupportedSource, InvalidConnectivity,
    InvalidLayerControl, InvalidClassification, InvalidCollisionSurface,
    CollisionBudgetExceeded, LayerCollapsed, SweepFailed,
    AllocationFailed, InternalError
};
struct CurvedPrismMetrics {
    std::uint64_t sourceNodes = 0;
    std::uint64_t sourceTriangles = 0;
    std::uint64_t collisionTriangles = 0;
    std::uint64_t rayQueries = 0;
    std::uint64_t clippedNodes = 0;
    std::uint64_t featureNodes = 0;
    std::uint32_t backoffSteps = 0;
    double requestedThickness = 0.0;
    double minimumThickness = 0.0;
    double maximumThickness = 0.0;
    PreparedSweepMetrics sweep{};
};
struct CurvedPrismOutcome {
    bool published = false;
    CurvedPrismReason reason = CurvedPrismReason::InternalError;
    bool destinationUnchanged = true;
    CurvedPrismMetrics metrics{};
};

// Geometry-kernel-independent collision-aware layer front. It accepts a
// prepared oriented Tri3 patch, computes smoothed feature-aware node normals,
// clips them against a source or supplied obstacle surface, and delegates the
// final direct Wedge6 emission to PreparedMixedSweep.
CurvedPrismOutcome generateCurvedPrismLayers(
    const mesh::Mesh& sourceSurface,
    const CurvedPrismLayerRequest& request,
    mesh::Mesh& destination) noexcept;

}  // namespace paacs::algo
