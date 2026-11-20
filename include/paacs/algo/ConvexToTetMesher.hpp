// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include <atomic>

namespace paacs::algo {
enum class ConvexToTetReason : std::uint8_t {
    None,
    InvalidRequest,
    InvalidSource,
    NonConvexCell,
    NonplanarFace,
    BoundaryAuditFailed,
    OutputBudget,
    WorkingBudget,
    QualityRequirement,
    Cancelled,
    AllocationFailed,
    InternalError
};
struct ConvexToTetOptions {
    std::uint64_t maximumNodes = 5'000'000, maximumCells = 5'000'000;
    std::uint64_t maximumOutputBytes = 512ULL * 1024 * 1024;
    std::uint64_t maximumWorkingBytes = 512ULL * 1024 * 1024;
    mesh::LinearVolumeAuditLimits auditLimits{};
    double relativePlanarityTolerance = 1e-10;
    double minimumQuality = 0, warningQuality = .1;
    const std::atomic<bool> *cancellation = nullptr;
};
struct ConvexToTetResult {
    bool published = false, destinationUnchanged = true;
    ConvexToTetReason reason = ConvexToTetReason::InternalError;
    std::uint64_t tetrahedra = 0, boundaryTriangles = 0, outputBytes = 0, workingBytes = 0;
    std::uint64_t lowQualityCells = 0;
    double minimumQuality = 1, volume = 0;
    mesh::LinearVolumeAudit audit{};
};
// Conforming pulling triangulation of a single-region Tet4/Prism6/Hex8 mesh
// with a complete oriented Tri3/Quad4 exterior. Requires convex cells and
// planar faces; rejects warped/concave cells. No additional nodes: a prism
// gives 3 tets and a hex 6. Quad diagonals depend on shared global node IDs.
// Preserves node/face/body tags and piecewise planar geometry. Does not detect
// global intersections or project curved CAD. Work/output limits cover array
// payload separately; source and caller's old destination are not included.
// Full source and result topology checks; failure leaves destination intact.
ConvexToTetResult convertConvexCellsToTet4(const mesh::Mesh &, const ConvexToTetOptions &,
                                           mesh::Mesh &destination) noexcept;
const char *convexToTetReasonName(ConvexToTetReason) noexcept;
} // namespace paacs::algo
