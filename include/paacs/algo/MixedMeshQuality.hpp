// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/mesh/Mesh.hpp"

#include <cstdint>
#include <limits>

namespace paacs::algo {

enum class MixedQualityMetric : std::uint8_t { SubTetMeanRatio };
struct MixedMeshQualityReport {
    bool complete = true;
    bool valid = false;
    std::uint64_t volumeCells = 0;
    std::uint64_t tet4 = 0;
    std::uint64_t pyramid5 = 0;
    std::uint64_t prism6 = 0;
    std::uint64_t hex8 = 0;
    std::uint64_t unsupportedCells = 0;
    std::uint64_t invalidConnectivity = 0;
    std::uint64_t invertedCells = 0;
    std::uint64_t nonfiniteQuality = 0;
    double minimumQuality = 0.0;
    double percentile01 = 0.0;
    double percentile1 = 0.0;
    double percentile5 = 0.0;
    double medianQuality = 0.0;
    double meanQuality = 0.0;
    double minimumVolume = 0.0;
    double totalVolume = 0.0;
    // Distinguishes an exhausted allocation from other incomplete analyses.
    bool allocationFailed = false;
};
MixedMeshQualityReport analyzeMixedMeshQuality(
    const mesh::Mesh& value,
    MixedQualityMetric metric = MixedQualityMetric::SubTetMeanRatio) noexcept;

struct MixedMeshOptimizationRequest {
    std::uint32_t maximumIterations = 4;
    double relaxation = 0.5;
    double minimumAcceptedImprovement = 1.0e-6;
    std::uint64_t maximumNodeTrials =
        (std::numeric_limits<std::uint64_t>::max)();
};
enum class MixedMeshOptimizationReason : std::uint8_t {
    None, InvalidRequest, InvalidInput, BudgetExceeded,
    AllocationFailed, InternalError
};
struct MixedMeshOptimizationOutcome {
    bool completed = false;
    MixedMeshOptimizationReason reason = MixedMeshOptimizationReason::InternalError;
    std::uint64_t iterations = 0;
    std::uint64_t nodeTrials = 0;
    std::uint64_t acceptedMoves = 0;
    MixedMeshQualityReport before{};
    MixedMeshQualityReport after{};
};

// Deterministic smart-Laplacian repair. Only nodes classified on a volume
// entity (dim=3) are eligible. Every move is accepted only when all incident
// mixed cells remain positive and the incident worst quality improves.
MixedMeshOptimizationOutcome optimizeMixedMeshQuality(
    mesh::Mesh& value,
    const MixedMeshOptimizationRequest& request = {}) noexcept;

}  // namespace paacs::algo
