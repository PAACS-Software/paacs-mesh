// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/mesh/Mesh.hpp"
#include <atomic>
#include <cstdint>
#include <span>

namespace paacs::mesh {
struct LinearVolumeAuditLimits {
    std::uint64_t maximumFaceUses = 12000000;
    // Requested vector storage, excluding the caller's mesh and allocator
    // bookkeeping. Checked before allocating either temporary array.
    std::uint64_t maximumTemporaryBytes = 512ULL * 1024 * 1024;
    const std::atomic<bool> *cancellation = nullptr;
    // Optional, caller-owned positive volumes indexed by ALL mesh cell IDs;
    // boundary/line entries are ignored. Piggybacks adjacent-volume growth
    // measurement on the same exact face grouping. Empty preserves the compact
    // legacy audit path. The audit does not certify these supplied volumes.
    std::span<const double> cellVolumes{};
};
struct LinearVolumeAudit {
    bool complete = false;
    bool valid = false;
    bool budgetExceeded = false;
    bool cancelled = false, allocationFailed = false;
    std::uint64_t temporaryBytes = 0, maximumBucketSize = 0;
    std::uint64_t volumeCells = 0, boundaryCells = 0, faceUses = 0;
    std::uint64_t exteriorFaces = 0, interiorFaces = 0;
    std::uint64_t nonmanifoldFaces = 0, inconsistentOrientation = 0;
    std::uint64_t missingBoundary = 0, unexpectedBoundary = 0, duplicateBoundary = 0;
    std::uint64_t invalidConnectivity = 0, nonfiniteNodes = 0, unsupportedCells = 0;
    std::uint64_t invalidVolumeData = 0;
    bool volumeGrowthMeasured = false;
    double maximumAdjacentVolumeRatio = 1, maximumAdjacentSizeRatio = 1;
};
// Exact face-incidence audit for Tet4/Hex8/Prism6/Pyramid5 and their Tri3/Quad4
// exterior faces. Groups by the smallest corner, then sorts full cyclic keys
// within each group; no probabilistic hashes or construction receipts are
// trusted. Uses 32-bit temporary corner IDs when possible (Mesh stays 64-bit).
// Line2 cells are ignored. Jacobians, geometric intersections and CAD
// ownership need separate checks. Never modifies the mesh. Cancellation is
// polled during scans and between groups, not inside a single group's sort.
LinearVolumeAudit auditLinearVolumeBoundary(const Mesh &mesh,
                                            const LinearVolumeAuditLimits &limits = {}) noexcept;
} // namespace paacs::mesh
