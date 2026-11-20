// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/algo/Delaunay2D.hpp"
#include <span>

namespace paacs::algo {
struct PlanarSegment {
    Index first = -1, second = -1;
};
enum class PlanarRecoveryReason : std::uint8_t {
    None,
    InvalidInput,
    InvalidTriangulation,
    CrossingConstraints,
    UnsplitSegment,
    WorkBudget,
    NotRecoverable,
    AllocationFailed,
    InternalError
};
struct PlanarRecoveryOutcome {
    bool completed = false;
    PlanarRecoveryReason reason = PlanarRecoveryReason::InternalError;
    std::uint64_t work = 0, flips = 0, recovered = 0;
};
// Bounded deterministic edge-flip insertion into a valid planar triangulation.
// Points and segment endpoints are preserved; no Steiner points are added.
// Uses exact orientation signs. Collinear vertices must already split a
// requested segment. The result is constrained, but is not re-legalized to
// constrained Delaunay. On failure the triangle vector remains unchanged.
PlanarRecoveryOutcome recoverPlanarConstraints(std::span<const Vec2> points,
                                               std::span<const PlanarSegment> segments,
                                               std::vector<Tri> &triangles,
                                               std::uint64_t maximumWork = 50000000) noexcept;
} // namespace paacs::algo
