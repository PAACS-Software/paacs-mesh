// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/mesh/Mesh.hpp"
#include <atomic>
#include <cstdint>
#include <vector>

namespace paacs::algo {
enum class SurfaceFeatureFlag : std::uint8_t { Sharp = 1, OpenBoundary = 2, TagBoundary = 4 };
struct SurfaceFeatureEdge {
    Index first = kInvalidIndex, second = kInvalidIndex;
    Index firstFace = kInvalidIndex, secondFace = kInvalidIndex;
    double length = 0, angleDegrees = 0;
    std::uint8_t flags = 0;
};
struct SurfaceFeatureOptions {
    // Disabled is a constant-time no-op: no input traversal or allocation.
    bool enabled = false;
    double angleDegrees = 30;
    // Mark sharp only when measured angle > angleDegrees + angularToleranceDegrees.
    double angularToleranceDegrees = 1e-8;
    // Filters the returned edges only; never suppresses topology checks or
    // changes the mesh. This is not a minimum element size or defeaturing.
    double minimumLength = 0;
    bool includeOpenBoundaries = true, includeTagBoundaries = false;
    std::uint64_t maximumEdgeUses = 12'000'000, maximumEdges = 1'000'000;
    std::uint64_t maximumWorkingBytes = 256ULL * 1024 * 1024;
    std::uint64_t maximumOutputBytes = 64ULL * 1024 * 1024;
    const std::atomic<bool> *cancellation = nullptr;
};
enum class SurfaceFeatureReason : std::uint8_t {
    None,
    Disabled,
    InvalidRequest,
    EmptySurface,
    InvalidConnectivity,
    NumericallyUnresolved,
    NonManifoldSurface,
    InconsistentWinding,
    WorkingBudget,
    OutputBudget,
    Cancelled,
    AllocationFailed,
    InternalError
};
struct SurfaceFeatureMetrics {
    std::uint64_t surfaceFaces = 0, highOrderFaces = 0, edgeUses = 0, uniqueEdges = 0;
    // Candidate counts precede the length filter. A sharp edge can also be a tag boundary.
    std::uint64_t sharpEdges = 0, openBoundaryEdges = 0, tagBoundaryEdges = 0;
    std::uint64_t filteredShortEdges = 0, selectedEdges = 0;
    std::uint64_t nonManifoldEdges = 0, inconsistentWindingEdges = 0;
    std::uint64_t workingBytes = 0, outputBytes = 0;
};
struct SurfaceFeatureResult {
    bool completed = false;
    SurfaceFeatureReason reason = SurfaceFeatureReason::InternalError;
    SurfaceFeatureMetrics metrics{};
    // Ascending unordered node pair; face IDs are original mesh cell ordinals.
    // Empty on every failure, even if diagnostic counts are partially available.
    std::vector<SurfaceFeatureEdge> edges;
};
// Opt-in analysis of explicit Tri3/Tri6/Quad4 surface cells. Volume and line
// cells are ignored; volume boundaries are not extracted automatically.
// Tri6 uses the linear corner normal/edge length, not a curved-surface test.
// Quad4 uses its normalized Newell normal. Input winding must be consistent.
// Open surfaces are supported. Nonmanifold or inconsistently wound adjacency
// refuses the result. This is not a self-intersection or CAD validity audit.
// No mesh mutation, node welding, smoothing, refinement, or CAD tolerance change.
// O(F + U log U) time, bounded flat arrays; cancellation polls around sort and
// throughout linear passes (the in-place sort itself is not interruptible).
SurfaceFeatureResult detectSurfaceFeatures(const mesh::Mesh &mesh,
                                           const SurfaceFeatureOptions &options = {}) noexcept;
const char *surfaceFeatureReasonName(SurfaceFeatureReason reason) noexcept;
} // namespace paacs::algo
