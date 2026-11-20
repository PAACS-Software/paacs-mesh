// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/mesh/Mesh.hpp"

#include <cstdint>
#include <limits>
#include <vector>

namespace paacs::algo {

enum class PreparedSweepStatus : std::uint8_t { Published, Rejected, Failed };
enum class PreparedSweepReason : std::uint8_t {
    None,
    EmptySource,
    UnsupportedSourceCell,
    InvalidSourceConnectivity,
    NonmanifoldSource,
    InconsistentSourceOrientation,
    InvalidTarget,
    InvalidLayerCount,
    InvalidGrowth,
    InvalidThickness,
    InvalidClassification,
    CountOverflow,
    BudgetExceeded,
    AllocationFailed,
    NonPositiveCell,
    AuditFailed,
    InternalError
};
enum class PreparedSweepPhase : std::uint8_t {
    None, Preflight, Schedule, Construction, Audit, Publication
};

struct SweepSideClassification {
    Index first = kInvalidIndex, second = kInvalidIndex;
    mesh::Classification classification{2, 3};
};
struct PreparedMixedSweepRequest {
    std::uint32_t layers = 1;
    // Positive ratio; values below one refine towards the target. Both q
    // and 1/q must be within maximumAllowedGrowth.
    double cellToCellGrowth = 1.0;
    double maximumAllowedGrowth = 1.3;
    Vec3 direction{0.0, 0.0, 1.0};
    double totalThickness = 1.0;
    std::vector<Vec3> targetNodes;
    std::vector<mesh::Classification> targetNodeClassifications;
    // Optional topology-preserving maps for geometry-recognized sweeps.
    // Node arrays must have source.nodeCount() entries; side records must
    // cover each undirected source boundary edge exactly once.
    std::vector<mesh::Classification> intermediateNodeClassifications;
    std::vector<SweepSideClassification> sideEdgeClassifications;
    mesh::Classification volumeClassification{3, 1};
    mesh::Classification targetClassification{2, 2};
    mesh::Classification sideClassification{2, 3};
    bool emitBoundaryFaces = true;
    // FullAudit checks every emitted volume cell.  When false, a straight
    // translation uses an analytic source-cell certificate and verifies only
    // the source rows, so the cost is independent of the number of layers.
    // Explicit source-to-target maps and translations near coordinate-resolution
    // limits retain the full emitted-cell audit.
    bool fullAudit = true;
    int executionThreads = 0;
    std::uint64_t maximumNodes =
        (std::numeric_limits<std::uint64_t>::max)();
    std::uint64_t maximumCells =
        (std::numeric_limits<std::uint64_t>::max)();
    std::uint64_t maximumBoundaryFaces =
        (std::numeric_limits<std::uint64_t>::max)();
    // Optional normalized boundaries: exactly layers+1 entries, 0 ... 1,
    // strictly increasing. Overrides geometric cellToCellGrowth spacing.
    // Actual adjacent width ratios must satisfy maximumAllowedGrowth.
    std::vector<double> layerFractions;
};

struct PreparedSweepAuditReport {
    bool complete = true;
    bool valid = false;
    bool allocationFailed = false;
    std::uint64_t nodes = 0;
    std::uint64_t wedges = 0;
    std::uint64_t hexes = 0;
    std::uint64_t boundaryTriangles = 0;
    std::uint64_t boundaryQuads = 0;
    std::uint64_t unsupportedCells = 0;
    std::uint64_t invalidConnectivity = 0;
    std::uint64_t invalidClassifications = 0;
    std::uint64_t nonPositiveCells = 0;
    double minimumCellVolume = 0.0;
    double maximumCellVolume = 0.0;
    double totalVolume = 0.0;
};
struct PreparedSweepMetrics {
    std::uint64_t sourceNodes = 0;
    std::uint64_t sourceTriangles = 0;
    std::uint64_t sourceQuads = 0;
    std::uint64_t boundaryEdges = 0;
    std::uint64_t nodes = 0;
    std::uint64_t wedges = 0;
    std::uint64_t hexes = 0;
    std::uint64_t boundaryFaces = 0;
    std::uint64_t connectivityEntries = 0;
    double maximumAdjacentGrowth = 1.0;
    double minimumCellVolume = 0.0;
    double totalVolume = 0.0;
    PreparedSweepAuditReport audit{};
    bool usedFullGeometricAudit = false;
    bool usedExplicitLayerFractions = false;
};
struct PreparedSweepOutcome {
    PreparedSweepStatus status = PreparedSweepStatus::Failed;
    PreparedSweepReason reason = PreparedSweepReason::InternalError;
    PreparedSweepPhase phase = PreparedSweepPhase::None;
    bool destinationUnchanged = true;
    PreparedSweepMetrics metrics{};
    bool published() const noexcept { return status == PreparedSweepStatus::Published; }
};

PreparedSweepOutcome generatePreparedMixedSweep(
    const mesh::Mesh& sourceSurface,
    const PreparedMixedSweepRequest& request,
    mesh::Mesh& destination) noexcept;
PreparedSweepAuditReport auditPreparedMixedSweep(
    const mesh::Mesh& value) noexcept;

}  // namespace paacs::algo
