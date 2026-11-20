// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/StructuredBlock.hpp"
#include "paacs/algo/StructuredGrowth.hpp"

#include <cstdint>
#include <memory>

namespace paacs::algo {

enum class MappedHexStatus : std::uint8_t {
    Published,
    Rejected,
    Failed
};

enum class MappedHexReason : std::uint8_t {
    None,
    InvalidCellCounts,
    InvalidClassification,
    NonfiniteFrame,
    NonPositiveFrame,
    BudgetExceeded,
    InvalidGrowth,
    CountOverflow,
    AllocationFailed,
    AuditFailed,
    InternalError
};

enum class MappedHexPhase : std::uint8_t {
    None,
    Preflight,
    Construction,
    Audit,
    Publication
};

enum class MappedHexValidationMode : std::uint8_t {
    // Constant-cost certificate derived from the affine map, exact count
    // formulas and deterministic direct-index construction. Automatically
    // falls back to emitted-cell auditing near coordinate-resolution limits.
    CertifiedAffine,
    // General O(cells log cells) sampled-Jacobian and face-incidence audit.
    FullTopologyAudit
};

struct MappedHexBlockRequest {
    StructuredBlockGeometry geometry{};
    StructuredBlockGrowth growth{};
    // Zero uses the current request-independent default. Positive values
    // bound worker participation for this operation only.
    int executionThreads = 0;
    // When true, one outward-oriented Quad4 is emitted for every exterior
    // Hex8 face. Boundary-condition ownership therefore remains available in
    // the ordinary flat mesh format without assigning false CAD identities to
    // edge/corner nodes.
    bool emitBoundaryQuads = true;
    MappedHexValidationMode validation =
        MappedHexValidationMode::CertifiedAffine;
};

struct Hex8AuditReport {
    bool complete = true;
    bool valid = false;
    bool allocationFailed = false;
    std::uint64_t nodes = 0;
    std::uint64_t hexes = 0;
    std::uint64_t boundaryQuads = 0;
    std::uint64_t unsupportedCells = 0;
    std::uint64_t nonfiniteNodes = 0;
    std::uint64_t invalidConnectivity = 0;
    std::uint64_t duplicateCellNodes = 0;
    std::uint64_t invalidClassifications = 0;
    std::uint64_t nonPositiveJacobians = 0;
    std::uint64_t nonmanifoldFaces = 0;
    std::uint64_t missingBoundaryQuads = 0;
    std::uint64_t duplicateBoundaryQuads = 0;
    std::uint64_t unexpectedBoundaryQuads = 0;
    std::uint64_t boundaryFaces = 0;
    std::uint64_t interiorFaces = 0;
    double minimumJacobian = 0.0;
    double maximumJacobian = 0.0;
    double minimumScaledJacobian = 0.0;
    double totalVolume = 0.0;
};

struct MappedHexMetrics {
    std::uint64_t nodes = 0;
    std::uint64_t hexes = 0;
    std::uint64_t boundaryQuads = 0;
    std::uint64_t connectivityEntries = 0;
    std::uint64_t constructionWorkUnits = 0;
    std::uint64_t logicalBytes = 0;
    std::uint64_t growthPreparationWorkUnits = 0;
    double blockDeterminant = 0.0;
    double minimumCellFraction = 0.0;
    double maximumCellFraction = 0.0;
    double maximumAdjacentGrowth = 1.0;
    Hex8AuditReport audit{};
    // True when the reusable private construction buffers already had the
    // exact requested shape, so allocation/value-initialization was skipped.
    bool reusedConstructionStorage = false;
    // Also true when a coordinate-resolution guard requires emitted-cell checks.
    bool usedFullGeometricAudit = false;
};

struct MappedHexOutcome {
    MappedHexStatus status = MappedHexStatus::Failed;
    MappedHexReason reason = MappedHexReason::InternalError;
    MappedHexPhase phase = MappedHexPhase::None;
    bool destinationUnchanged = true;
    MappedHexMetrics metrics{};

    bool published() const noexcept {
        return status == MappedHexStatus::Published;
    }
};


// Optional request-local workspace for repeated mapped-block generation.
// It preserves the ordinary explicit Mesh result and atomic publication, but
// reuses private buffers across remeshes. The first two same-size calls prime
// the destination/workspace double buffer; subsequent calls perform no large
// construction allocation or value-initialization.
class MappedHexWorkspace {
public:
    MappedHexWorkspace();
    ~MappedHexWorkspace();
    MappedHexWorkspace(MappedHexWorkspace&&) noexcept;
    MappedHexWorkspace& operator=(MappedHexWorkspace&&) noexcept;
    MappedHexWorkspace(const MappedHexWorkspace&) = delete;
    MappedHexWorkspace& operator=(const MappedHexWorkspace&) = delete;

    // Releases the spare output-sized buffer retained for the next remesh.
    // The currently published destination Mesh is unaffected.
    void release() noexcept;
    std::uint64_t retainedLogicalBytes() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    friend MappedHexOutcome generateMappedHexBlock(
        const MappedHexBlockRequest&, MappedHexWorkspace&,
        mesh::Mesh&) noexcept;
};

// Generates a deterministic affine mapped Hex8 block directly into the normal
// PAACS flat mesh representation. Construction is private and publication is
// one noexcept swap, so every rejection/failure leaves destination unchanged.
MappedHexOutcome generateMappedHexBlock(
    const MappedHexBlockRequest& request,
    mesh::Mesh& destination) noexcept;

// Repeated-generation overload using a private double-buffer workspace.
MappedHexOutcome generateMappedHexBlock(
    const MappedHexBlockRequest& request,
    MappedHexWorkspace& workspace,
    mesh::Mesh& destination) noexcept;

// Full sampled-Jacobian and face-incidence audit for Hex8 meshes. It accepts
// optional Quad4 boundary cells and can require exact one-to-one boundary
// coverage. The 2x2x2 Gauss volume integration is exact for trilinear Hex8.
Hex8AuditReport auditHex8Mesh(
    const mesh::Mesh& value,
    bool requireBoundaryQuads = false) noexcept;

}  // namespace paacs::algo
