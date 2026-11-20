// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/algo/PreparedSweepMesher.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/sizing/CappedSizeProfile.hpp"
#include <atomic>

namespace paacs::algo {
struct LayerScheduleRequest {
    std::uint32_t wallLayers = 5;
    double firstLayerThickness = .01;
    double growthRatio = 1.2, maximumAllowedGrowth = 1.3;
    // Zero ends at the last wall layer. Otherwise fill exactly this depth,
    // retaining all specified wall layers, then grow towards the core cap.
    double totalDepth = 1;
    double maximumCoreThickness = .15;
    std::uint32_t maximumLayers = 4096;
    std::uint64_t maximumBytes = 1024 * 1024;
};
enum class LayerScheduleReason : std::uint8_t {
    None,
    InvalidRequest,
    LayerBudget,
    ByteBudget,
    WallLayersExceedDepth,
    TransitionDoesNotFit,
    NumericallyUnresolved,
    AllocationFailed,
    InternalError
};
struct LayerSchedule {
    bool planned = false;
    LayerScheduleReason reason = LayerScheduleReason::InternalError;
    std::uint32_t wallLayers = 0, coreLayers = 0;
    double wallThickness = 0, totalDepth = 0, maximumAdjacentGrowth = 1;
    double minimumThickness = 0, maximumThickness = 0;
    std::vector<double> offsets;
};
LayerSchedule planLayerSchedule(const LayerScheduleRequest &) noexcept;

struct LayeredPatchRequest {
    LayerScheduleRequest schedule{};
    Vec3 direction{0, 0, 1};
    mesh::Classification volumeClassification{3, 1}, targetClassification{2, 2},
        sideClassification{2, 3};
    std::uint64_t maximumSourceCells = 250'000, maximumSourceNodes = 500'000;
    // Includes every overlapping-x sweep candidate, even when its y bounds
    // reject it. Exhaustion refuses; intersection checking is never skipped.
    std::uint64_t maximumSourcePairChecks = 20'000'000;
    std::uint64_t maximumNodes = 5'000'000, maximumCells = 5'000'000;
    std::uint64_t maximumOutputBytes = 512ULL * 1024 * 1024;
    std::uint64_t maximumWorkingBytes = 512ULL * 1024 * 1024;
    mesh::LinearVolumeAuditLimits auditLimits{};
    int threads = 1;
    const std::atomic<bool> *cancellation = nullptr;
    // Optional reusable normal-size profile replaces the wall/core schedule.
    // Does not modify tangential source topology. Caller retains the profile.
    const sizing::CappedSizeProfile *sizeProfile = nullptr;
    sizing::ProfileDistributionLimits profileLimits{};
};
enum class LayeredPatchReason : std::uint8_t {
    None,
    InvalidRequest,
    ScheduleRejected,
    InvalidSource,
    NonplanarSource,
    IntersectingSource,
    SourceBudget,
    OutputBudget,
    WorkingBudget,
    Cancelled,
    SweepRejected,
    BoundaryAuditFailed,
    AllocationFailed,
    InternalError
};
struct LayeredPatchResult {
    bool published = false, destinationUnchanged = true;
    LayeredPatchReason reason = LayeredPatchReason::InternalError;
    LayerSchedule schedule{};
    PreparedSweepOutcome sweep{};
    mesh::LinearVolumeAudit audit{};
    std::uint64_t sourcePairChecks = 0, outputBytes = 0, plannedOutputBytes = 0,
                  plannedWorkingBytes = 0;
    bool usedSizeProfile = false;
    std::uint64_t profileWorkingBytes = 0, plateauCellsRemoved = 0;
};
// Prepared planar Tri3/convex Quad4 patches, structured or unstructured.
// Straight normal extrusion: triangles become Prism6 (Wedge6), quads Hex8.
// Validates projected intersections/T-junctions, source winding via the sweep,
// positive subcells and full volume-boundary incidence before publication.
// This does not join a separate tetrahedral core or handle curved wall fronts.
LayeredPatchResult meshLayeredPatch(const mesh::Mesh &source, const LayeredPatchRequest &,
                                    mesh::Mesh &destination) noexcept;
} // namespace paacs::algo
