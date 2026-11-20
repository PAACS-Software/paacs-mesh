// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/mesh/Mesh.hpp"

#include <cstdint>
#include <limits>
#include <span>

namespace paacs::algo {

enum class TransitionTemplateKind : std::uint8_t {
    QuadToApexPyramid5,
    Pyramid5ToTet4,
    Wedge6ToTet4,
    Hex8ToTet4
};
struct TransitionTemplateInfo {
    TransitionTemplateKind kind = TransitionTemplateKind::QuadToApexPyramid5;
    mesh::CellType input = mesh::CellType::Quad4;
    mesh::CellType output = mesh::CellType::Pyramid5;
    std::uint8_t outputCells = 1;
    bool conformingForSharedCanonicalFaces = true;
};
std::span<const TransitionTemplateInfo> transitionTemplateCatalog() noexcept;

enum class TransitionStatus : std::uint8_t { Published, Rejected, Failed };
enum class TransitionReason : std::uint8_t {
    None, EmptyInput, UnsupportedCell, InvalidConnectivity,
    InvalidClassification, NonfiniteGeometry, DegenerateCell,
    NonPositiveCell, CountOverflow, BudgetExceeded, AllocationFailed,
    AuditFailed, InternalError
};
struct PyramidTransitionRequest {
    Vec3 apex{};
    mesh::Classification volumeClassification{3,1};
    mesh::Classification apexClassification{3,1};
    mesh::Classification sideClassification{2,2};
    bool emitBoundaryFaces = true;
    std::uint64_t maximumCells =
        (std::numeric_limits<std::uint64_t>::max)();
};
struct TransitionMetrics {
    std::uint64_t inputCells = 0;
    std::uint64_t outputNodes = 0;
    std::uint64_t outputTet4 = 0;
    std::uint64_t outputPyramid5 = 0;
    std::uint64_t boundaryTriangles = 0;
    std::uint64_t boundaryQuads = 0;
    double minimumVolume = 0.0;
    double totalVolume = 0.0;
};
struct TransitionOutcome {
    TransitionStatus status = TransitionStatus::Failed;
    TransitionReason reason = TransitionReason::InternalError;
    bool destinationUnchanged = true;
    TransitionMetrics metrics{};
    bool published() const noexcept { return status == TransitionStatus::Published; }
};

TransitionOutcome generatePyramidTransitionFan(
    const mesh::Mesh& quadSurface,
    const PyramidTransitionRequest& request,
    mesh::Mesh& destination) noexcept;

struct MixedToTet4Request {
    bool retainBoundary = true;
    bool triangulateBoundaryQuads = true;
    std::uint64_t maximumTet4 =
        (std::numeric_limits<std::uint64_t>::max)();
};
TransitionOutcome convertSupportedMixedMeshToTet4(
    const mesh::Mesh& source,
    const MixedToTet4Request& request,
    mesh::Mesh& destination) noexcept;

}  // namespace paacs::algo
