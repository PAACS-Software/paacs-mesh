// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/sizing/RegionSizeField.hpp"
#include <atomic>
#include <cstdint>

namespace paacs::algo {
struct AdaptiveBoxRequest {
    Vec3 lower{}, upper{1, 1, 1};
    double targetSize = .25;
    const sizing::BoundedSizeField *field = nullptr;
    std::uint32_t maximumDepth = 12;
    std::uint64_t maximumLeaves = 250'000, maximumNodes = 3'000'000, maximumCells = 5'000'000;
    std::uint64_t maximumOutputBytes = 512ULL * 1024 * 1024;
    std::uint64_t maximumWorkingBytes = 512ULL * 1024 * 1024;
    double minimumQuality = 0;
    int threads = 1;
    const std::atomic<bool> *cancellation = nullptr;
    mesh::LinearVolumeAuditLimits auditLimits{};
    // Zero reports only. A positive limit >=1 rejects a mesh whose measured
    // face-neighbor cbrt(volume) ratio exceeds it (up to roundoff).
    double maximumAdjacentSizeRatio = 0;
    // Zero disables this optional query-work gate for legacy callers.
    std::uint64_t maximumSizingWork = 0;
};
enum class AdaptiveBoxReason : std::uint8_t {
    None,
    InvalidRequest,
    DepthLimit,
    LeafBudget,
    OutputBudget,
    WorkingBudget,
    Cancelled,
    NumericallyUnresolved,
    InvalidCell,
    QualityRequirement,
    BoundaryAuditFailed,
    AllocationFailed,
    InternalError,
    GrowthRequirement,
    SizingBudget
};
struct AdaptiveBoxResult {
    bool published = false, destinationUnchanged = true;
    AdaptiveBoxReason reason = AdaptiveBoxReason::InternalError;
    std::uint64_t leaves = 0, treeNodes = 0, nodes = 0, tetrahedra = 0, boundaryTriangles = 0;
    std::uint64_t balanceSplits = 0, sizingQueries = 0, outputBytes = 0, workingBytes = 0;
    std::uint32_t deepestLevel = 0, balancePasses = 0;
    double smallestBlockWidth = 0, largestBlockWidth = 0;
    double minimumQuality = 0, minimumCellVolume = 0, maximumCellVolume = 0;
    double analyticVolume = 0, meshVolume = 0;
    double planningSeconds = 0, generationSeconds = 0, validationSeconds = 0, auditSeconds = 0;
    mesh::LinearVolumeAudit audit{};
    std::uint64_t sizingWork = 0;
};
// Axis-aligned box, isotropic octree refinement, strong 2:1 balance, conforming
// Tet4 cells and classified Tri3 exterior. Region controls are influence zones,
// not material interfaces. Field minima, not sparse samples, choose refinement.
// The target caps the maximum leaf-block width AND all tetrahedron edges
// throughout that leaf (up to coordinate roundoff). Dyadic refinement
// quantizes actual cell sizes.
// Face tags: 1=xmin, 2=xmax, 3=ymin, 4=ymax, 5=zmin, 6=zmax; volume tag 1.
// Every failure leaves destination unchanged. Full boundary audit always runs.
AdaptiveBoxResult meshAdaptiveBox(const AdaptiveBoxRequest &, mesh::Mesh &) noexcept;
const char *adaptiveBoxReasonName(AdaptiveBoxReason) noexcept;
} // namespace paacs::algo
