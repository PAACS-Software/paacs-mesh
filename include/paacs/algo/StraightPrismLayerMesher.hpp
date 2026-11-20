// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/mesh/Mesh.hpp"
#include "paacs/sizing/MeshControls.hpp"

#include <cstdint>
#include <limits>

namespace paacs::algo {

enum class PrismLayerStatus : std::uint8_t {
    Published,
    Rejected,
    Failed
};

enum class PrismLayerReason : std::uint8_t {
    None,
    EmptySource,
    UnsupportedSourceCell,
    InvalidSourceConnectivity,
    NonmanifoldSource,
    InconsistentSourceOrientation,
    NonfiniteDirection,
    ZeroDirection,
    InvalidLayerCount,
    InvalidFirstLayerThickness,
    InvalidGrowthRate,
    GrowthLimitExceeded,
    ThicknessBudgetExceeded,
    InvalidClassification,
    BudgetExceeded,
    CountOverflow,
    AllocationFailed,
    AuditFailed,
    InternalError
};

enum class PrismLayerPhase : std::uint8_t {
    None,
    Preflight,
    Schedule,
    Construction,
    Audit,
    Publication
};

enum class PrismLayerValidationMode : std::uint8_t {
    CertifiedStraightSweep,
    FullTopologyAudit
};

struct StraightPrismLayerRequest {
    Vec3 direction{0.0, 0.0, 1.0};
    std::uint32_t layers = 1;
    double firstLayerThickness = 1.0;
    // Width(k+1) / width(k).  V6 intentionally supports outward-growing
    // layers only. A future curved-normal engine can add compression safely.
    double cellToCellGrowth = 1.0;
    double maximumAllowedGrowth = 1.3;
    double maximumTotalThickness =
        (std::numeric_limits<double>::max)();
    mesh::Classification volumeClassification{3, 1};
    mesh::Classification targetClassification{2, 2};
    mesh::Classification sideClassification{2, 3};
    bool emitBoundaryFaces = true;
    // Zero uses the current request-independent default. Positive values
    // bound worker participation for this operation only.
    int executionThreads = 0;
    PrismLayerValidationMode validation =
        PrismLayerValidationMode::CertifiedStraightSweep;
    std::uint64_t maxNodes =
        (std::numeric_limits<std::uint64_t>::max)();
    std::uint64_t maxCells =
        (std::numeric_limits<std::uint64_t>::max)();
    std::uint64_t maxBoundaryFaces =
        (std::numeric_limits<std::uint64_t>::max)();
};

struct PrismLayerAuditReport {
    bool complete = true;
    bool valid = false;
    bool allocationFailed = false;
    std::uint64_t nodes = 0;
    std::uint64_t prisms = 0;
    std::uint64_t boundaryTriangles = 0;
    std::uint64_t boundaryQuads = 0;
    std::uint64_t unsupportedCells = 0;
    std::uint64_t invalidConnectivity = 0;
    std::uint64_t duplicateCellNodes = 0;
    std::uint64_t invalidClassifications = 0;
    std::uint64_t nonPositiveVolumes = 0;
    std::uint64_t nonmanifoldFaces = 0;
    std::uint64_t missingBoundaryFaces = 0;
    std::uint64_t duplicateBoundaryFaces = 0;
    std::uint64_t unexpectedBoundaryFaces = 0;
    std::uint64_t interiorFaces = 0;
    double minimumPrismVolume = 0.0;
    double maximumPrismVolume = 0.0;
    double totalVolume = 0.0;
};

struct PrismLayerMetrics {
    std::uint64_t sourceNodes = 0;
    std::uint64_t sourceTriangles = 0;
    std::uint64_t sourceBoundaryEdges = 0;
    std::uint64_t nodes = 0;
    std::uint64_t prisms = 0;
    std::uint64_t boundaryTriangles = 0;
    std::uint64_t boundaryQuads = 0;
    std::uint64_t connectivityEntries = 0;
    std::uint64_t constructionWorkUnits = 0;
    std::uint64_t logicalBytes = 0;
    double totalThickness = 0.0;
    double maximumAdjacentGrowth = 1.0;
    PrismLayerAuditReport audit{};
};

struct PrismLayerOutcome {
    PrismLayerStatus status = PrismLayerStatus::Failed;
    PrismLayerReason reason = PrismLayerReason::InternalError;
    PrismLayerPhase phase = PrismLayerPhase::None;
    bool destinationUnchanged = true;
    PrismLayerMetrics metrics{};

    bool published() const noexcept {
        return status == PrismLayerStatus::Published;
    }
};

// Extrudes an already-valid, consistently oriented Tri3 surface patch along a
// single direction.  It is a deliberately bounded first prism engine: no
// curved normals, collision clipping or CAD-side face decomposition is
// claimed.  Those require a separate layer-front engine.
PrismLayerOutcome generateStraightPrismLayers(
    const mesh::Mesh& sourceSurface,
    const StraightPrismLayerRequest& request,
    mesh::Mesh& destination) noexcept;

PrismLayerAuditReport auditPrismLayerMesh(
    const mesh::Mesh& value,
    bool requireBoundaryFaces = false) noexcept;

// Copies one validated compiled face-layer control into the straight-sweep
// request. Geometry direction and classifications remain explicit geometric
// inputs; all scalar layer policy is resolved before entering hot loops.
bool applyCompiledPrismLayerControl(
    int faceTag,
    const sizing::CompiledMeshControls& controls,
    StraightPrismLayerRequest& request) noexcept;

}  // namespace paacs::algo
