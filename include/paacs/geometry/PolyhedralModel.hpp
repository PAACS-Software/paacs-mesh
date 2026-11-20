// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/geometry/Kernel.hpp"
#include <cstdint>
#include <span>

namespace paacs::geo {

struct PolyhedralFace {
    // Zero-based vertex indices; no repeated closing vertex. Outer loop first,
    // then oppositely wound holes. Faces must form one consistently oriented
    // connected closed shell. A globally inward shell is reversed on build.
    std::vector<std::vector<std::uint32_t>> loops;
};
struct PolyhedralInput {
    std::vector<Vec3> vertices;
    std::vector<PolyhedralFace> faces;
    int solidTag = 1;
};
struct PolyhedralLimits {
    std::uint64_t maximumVertices = 100000;
    std::uint64_t maximumFaces = 100000;
    std::uint64_t maximumWork = 10000000;
    double relativeTolerance = 1.0e-10;
};
enum class PolyhedralReason : std::uint8_t {
    None,
    Empty,
    InvalidIndex,
    Nonfinite,
    Degenerate,
    Nonplanar,
    InvalidLoop,
    OpenOrNonmanifold,
    InconsistentOrientation,
    DisconnectedShell,
    BudgetExceeded,
    AllocationFailed,
    InternalError
};
struct PolyhedralBuild {
    std::unique_ptr<Model> model;
    PolyhedralReason reason = PolyhedralReason::InternalError;
    std::uint64_t work = 0;
    double volume = 0.0;
    bool built() const noexcept { return model != nullptr; }
};

// Validates planar loops, holes, edge incidence, orientation, shell connectivity
// and finite positive volume. It is not a general inter-face intersection or
// CAD-healing certificate. Geometry is immutable after construction.
PolyhedralBuild makePolyhedralModel(const PolyhedralInput &input,
                                    const PolyhedralLimits &limits = {}) noexcept;

// Convenience construction for a polygon with optional polygonal holes,
// translated by an arbitrary non-tangent sweep vector. Profiles use x/y in the
// plane origin + u*axisU + v*axisV. This still passes the full builder checks.
PolyhedralInput extrudedPolygon(std::span<const std::vector<Vec2>> loops, Vec3 sweep,
                                Vec3 origin = {}, Vec3 axisU = {1, 0, 0}, Vec3 axisV = {0, 1, 0});

const char *polyhedralReasonName(PolyhedralReason reason) noexcept;
} // namespace paacs::geo
