// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/StructuredBlock.hpp"
#include "paacs/algo/StructuredGrowth.hpp"
#include "paacs/mesh/PolyhedralMesh.hpp"

#include <cstdint>

namespace paacs::algo {

enum class StructuredPolyStatus : std::uint8_t {
    Published,
    Rejected,
    Failed
};

enum class StructuredPolyReason : std::uint8_t {
    None,
    InvalidCellCounts,
    InvalidClassification,
    NonfiniteFrame,
    NonPositiveFrame,
    BudgetExceeded,
    InvalidGrowth,
    CountOverflow,
    UnsupportedCellType,
    NonmanifoldInput,
    AllocationFailed,
    AuditFailed,
    InternalError
};

enum class StructuredPolyPhase : std::uint8_t {
    None,
    Preflight,
    Construction,
    Conversion,
    Audit,
    Publication
};

enum class StructuredPolyValidationMode : std::uint8_t {
    // Constant-cost certificate for the deterministic structured face layout.
    CertifiedStructured,
    // Full owner/neighbour, orientation and positive-volume audit.
    FullAudit
};

struct StructuredPolyMetrics {
    std::uint64_t nodes = 0;
    std::uint64_t faces = 0;
    std::uint64_t cells = 0;
    std::uint64_t faceNodeEntries = 0;
    std::uint64_t cellFaceEntries = 0;
    std::uint64_t constructionWorkUnits = 0;
    std::uint64_t logicalBytes = 0;
    std::uint64_t growthPreparationWorkUnits = 0;
    double blockDeterminant = 0.0;
    double minimumCellFraction = 0.0;
    double maximumCellFraction = 0.0;
    double maximumAdjacentGrowth = 1.0;
    mesh::PolyhedralAuditReport audit{};
};

struct StructuredPolyOutcome {
    StructuredPolyStatus status = StructuredPolyStatus::Failed;
    StructuredPolyReason reason = StructuredPolyReason::InternalError;
    StructuredPolyPhase phase = StructuredPolyPhase::None;
    bool destinationUnchanged = true;
    StructuredPolyMetrics metrics{};

    bool published() const noexcept {
        return status == StructuredPolyStatus::Published;
    }
};

struct StructuredPolyBlockRequest {
    StructuredBlockGeometry geometry{};
    StructuredBlockGrowth growth{};
    StructuredPolyValidationMode validation =
        StructuredPolyValidationMode::CertifiedStructured;
};

// Generates the same affine block as the mapped-Hex engine, but publishes a
// true owner/neighbour face-based mesh. Each initial cell is a six-face convex
// polyhedron; the representation itself supports arbitrary polygonal faces.
StructuredPolyOutcome generateStructuredPolyhedralBlock(
    const StructuredPolyBlockRequest& request,
    mesh::PolyhedralMesh& destination) noexcept;

// V5-compatible uniform-spacing convenience overload.
StructuredPolyOutcome generateStructuredPolyhedralBlock(
    const StructuredBlockGeometry& geometry,
    mesh::PolyhedralMesh& destination,
    StructuredPolyValidationMode validation =
        StructuredPolyValidationMode::CertifiedStructured) noexcept;

// Converts normal explicit Hex8 cells into the face-based polyhedral format.
// Quad4 cells, when present, supply boundary classifications. Other cell types
// are rejected rather than silently discarded.
StructuredPolyOutcome convertHex8ToPolyhedral(
    const mesh::Mesh& source,
    mesh::PolyhedralMesh& destination) noexcept;

// V6 mixed-cell conversion. Hex8 and Prism6 volume cells plus Tri3/Quad4
// boundary cells are accepted; every other type is rejected explicitly.
StructuredPolyOutcome convertSupportedVolumeMeshToPolyhedral(
    const mesh::Mesh& source,
    mesh::PolyhedralMesh& destination) noexcept;

}  // namespace paacs::algo
