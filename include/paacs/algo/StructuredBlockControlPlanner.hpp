// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/StructuredBlock.hpp"
#include "paacs/algo/StructuredGrowth.hpp"
#include "paacs/sizing/MeshControls.hpp"

#include <array>
#include <cstdint>

namespace paacs::algo {

// Compiles body and boundary-face sizing into the three integer counts needed
// by a structured volume engine.  It runs once per body; Hex/poly hot loops
// receive only final counts and precomputed axis coordinates.
struct StructuredBlockSizingPolicy {
    double baseSize = 1.0;
    double sizingScale = 1.0;
    std::array<std::uint32_t, 3> minimumCells{1, 1, 1};
    std::array<std::uint32_t, 3> maximumCells{4096, 4096, 4096};
    bool includeBoundaryFaceSizes = true;
    StructuredBlockGrowth growth{};
};

enum class StructuredBlockPlanReason : std::uint8_t {
    None,
    InvalidSolidTag,
    InvalidBaseSize,
    InvalidSizingScale,
    InvalidAxis,
    InvalidCellBounds,
    CountOverflow,
    GrowthRejected,
    UnsupportedEngine,
    InternalError
};

struct StructuredBlockControlPlan {
    bool planned = false;
    StructuredBlockPlanReason reason =
        StructuredBlockPlanReason::InternalError;
    sizing::VolumeEngineKind requestedEngine =
        sizing::VolumeEngineKind::Auto;
    StructuredBlockGeometry geometry{};
    StructuredBlockGrowth growth{};
    std::array<double, 3> effectiveTargetSizes{};
    std::uint64_t controlLookups = 0;
    std::uint64_t preparationWorkUnits = 0;
};

// Face-size precedence is directional.  A W face constrains U/V spacing, a V
// face constrains U/W, and a U face constrains V/W.  The body size constrains
// all three directions.  Auto is allowed; orchestration may then select Hex or
// poly only after proving the corresponding structured descriptor exists.
StructuredBlockControlPlan planStructuredBlockFromControls(
    const StructuredBlockGeometry& geometryTemplate,
    int solidTag,
    const sizing::CompiledMeshControls& controls,
    const StructuredBlockSizingPolicy& policy = {}) noexcept;

}  // namespace paacs::algo
