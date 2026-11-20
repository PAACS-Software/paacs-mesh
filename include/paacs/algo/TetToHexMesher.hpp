// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/algo/PrimitiveMesher.hpp"
namespace paacs::algo {
struct TetToHexOptions {
    std::uint64_t maximumNodes = 5'000'000, maximumCells = 5'000'000;
    std::uint64_t maximumOutputBytes = 512ULL * 1024 * 1024,
                  maximumWorkingBytes = 512ULL * 1024 * 1024;
    mesh::LinearVolumeAuditLimits auditLimits{};
    const std::atomic<bool> *cancellation = nullptr;
    int threads = 1;
    double minimumQuality = 0, warningQuality = .1;
    PrimitiveQualityPolicy qualityPolicy = PrimitiveQualityPolicy::Warn;
};
struct TetToHexResult {
    bool published = false, destinationUnchanged = true;
    PrimitiveReason reason = PrimitiveReason::InternalError;
    std::uint64_t nodes = 0, hexes = 0, boundaryQuads = 0, outputBytes = 0, workingBytes = 0;
    double volume = 0, minimumCellVolume = 0, minimumQuality = 0;
    PrimitiveQualityReport quality{};
    mesh::LinearVolumeAudit audit{};
};
// Single-region, nonoverlapping positive Tet4 source with complete, oriented
// Tri3 exterior. One shared node per edge/face; four Hex8 per Tet4. Each Tri3
// becomes three Quad4. Preserves the piecewise planar geometry; no CAD projection.
// At differently tagged faces, an edge's endpoints must share a curve/vertex
// classification. Unresolved ownership is refused. Does not test global overlap.
// Full Hex8 Jacobian and boundary checks; strong destination rollback on failure.
TetToHexResult convertTet4ToHex8(const mesh::Mesh &source, const TetToHexOptions &,
                                 mesh::Mesh &destination) noexcept;
} // namespace paacs::algo
