// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/StructuredBlockControlPlanner.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace paacs::algo {
namespace {

bool finitePositive(double value) noexcept {
    return std::isfinite(value) && value > 0.0;
}

std::size_t faceIndex(StructuredBlockFace face) noexcept {
    return static_cast<std::size_t>(face);
}

bool validStructuredEngine(sizing::VolumeEngineKind engine) noexcept {
    return engine == sizing::VolumeEngineKind::Auto ||
        engine == sizing::VolumeEngineKind::MappedHex8 ||
        engine == sizing::VolumeEngineKind::StructuredPolyhedral;
}

}  // namespace

StructuredBlockControlPlan planStructuredBlockFromControls(
        const StructuredBlockGeometry& geometryTemplate,
        int solidTag,
        const sizing::CompiledMeshControls& controls,
        const StructuredBlockSizingPolicy& policy) noexcept {
    StructuredBlockControlPlan out;
    out.geometry = geometryTemplate;
    out.growth = policy.growth;
    if (solidTag < 0) {
        out.reason = StructuredBlockPlanReason::InvalidSolidTag;
        return out;
    }
    if (!finitePositive(policy.baseSize)) {
        out.reason = StructuredBlockPlanReason::InvalidBaseSize;
        return out;
    }
    if (!finitePositive(policy.sizingScale)) {
        out.reason = StructuredBlockPlanReason::InvalidSizingScale;
        return out;
    }
    const std::array<double, 3> lengths{
        geometryTemplate.axisU.norm(), geometryTemplate.axisV.norm(),
        geometryTemplate.axisW.norm()};
    for (double length : lengths) {
        if (!finitePositive(length)) {
            out.reason = StructuredBlockPlanReason::InvalidAxis;
            return out;
        }
    }
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (policy.minimumCells[axis] == 0 ||
            policy.minimumCells[axis] > policy.maximumCells[axis]) {
            out.reason = StructuredBlockPlanReason::InvalidCellBounds;
            return out;
        }
    }

    out.requestedEngine = controls.engineForBody(solidTag);
    ++out.controlLookups;
    if (!validStructuredEngine(out.requestedEngine)) {
        out.reason = StructuredBlockPlanReason::UnsupportedEngine;
        return out;
    }

    std::array<double, 3> target{
        controls.effectiveSize(3, solidTag, policy.baseSize),
        controls.effectiveSize(3, solidTag, policy.baseSize),
        controls.effectiveSize(3, solidTag, policy.baseSize)};
    out.controlLookups += 3;

    if (policy.includeBoundaryFaceSizes) {
        const auto capFace = [&](StructuredBlockFace face,
                                 std::size_t firstAxis,
                                 std::size_t secondAxis) noexcept {
            const auto classification = geometryTemplate.
                boundaryClassifications[faceIndex(face)];
            if (classification.dim != 2 || classification.tag < 0) return;
            target[firstAxis] = controls.effectiveSize(
                2, classification.tag, target[firstAxis]);
            target[secondAxis] = controls.effectiveSize(
                2, classification.tag, target[secondAxis]);
            out.controlLookups += 2;
        };
        capFace(StructuredBlockFace::WMin, 0, 1);
        capFace(StructuredBlockFace::WMax, 0, 1);
        capFace(StructuredBlockFace::VMin, 0, 2);
        capFace(StructuredBlockFace::VMax, 0, 2);
        capFace(StructuredBlockFace::UMin, 1, 2);
        capFace(StructuredBlockFace::UMax, 1, 2);
    }

    std::array<std::uint32_t, 3> counts{};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        target[axis] *= policy.sizingScale;
        if (!finitePositive(target[axis])) {
            out.reason = StructuredBlockPlanReason::InvalidSizingScale;
            return out;
        }
        const double raw = std::ceil(lengths[axis] / target[axis]);
        if (!std::isfinite(raw) || raw > static_cast<double>(
                (std::numeric_limits<std::uint32_t>::max)())) {
            out.reason = StructuredBlockPlanReason::CountOverflow;
            return out;
        }
        const auto count = static_cast<std::uint32_t>(
            std::max(1.0, raw));
        counts[axis] = std::clamp(
            count, policy.minimumCells[axis], policy.maximumCells[axis]);
        out.effectiveTargetSizes[axis] = target[axis];
    }

    const std::array<StructuredAxisGrowth, 3> growth{
        policy.growth.u, policy.growth.v, policy.growth.w};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const auto distribution = buildStructuredAxisDistribution(
            counts[axis], growth[axis]);
        out.preparationWorkUnits +=
            static_cast<std::uint64_t>(counts[axis]) + 1U;
        if (!distribution.valid) {
            out.reason = StructuredBlockPlanReason::GrowthRejected;
            return out;
        }
    }

    out.geometry.cellsU = counts[0];
    out.geometry.cellsV = counts[1];
    out.geometry.cellsW = counts[2];
    out.preparationWorkUnits += out.controlLookups + 3U;
    out.planned = true;
    out.reason = StructuredBlockPlanReason::None;
    return out;
}

}  // namespace paacs::algo
