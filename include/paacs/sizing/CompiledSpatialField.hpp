// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/core/Types.hpp"
#include "paacs/sizing/SizeField.hpp"

#include <cstdint>
#include <vector>

namespace paacs::sizing {

enum class SpatialControlKind : std::uint8_t {
    Box,
    Ball,
    Cylinder,
    DistanceRamp
};

struct SpatialSizeControl {
    SpatialControlKind kind = SpatialControlKind::Ball;
    Vec3 a{};
    Vec3 b{};
    double radius = 0.0;
    double outerRadius = 0.0;
    double targetSize = 0.0;
    double outerSize = 0.0;
};

struct CompiledSpatialFieldRequest {
    Vec3 lower{};
    Vec3 upper{1.0, 1.0, 1.0};
    std::uint32_t samplesX = 32;
    std::uint32_t samplesY = 32;
    std::uint32_t samplesZ = 32;
    double baseSize = 1.0;
    double minimumSize = 0.0;
    double maximumSize = 1.0e300;
    double maximumAdjacentGrowth = 1.3;
    std::uint32_t maximumGrowthSweeps = 8;
    std::uint64_t maximumSamples = 64ULL * 1024ULL * 1024ULL;
    int executionThreads = 0;
    std::vector<SpatialSizeControl> controls;
};

enum class SpatialFieldCompileReason : std::uint8_t {
    None,
    InvalidBounds,
    InvalidResolution,
    InvalidSize,
    InvalidGrowth,
    InvalidControl,
    BudgetExceeded,
    CountOverflow,
    AllocationFailed,
    NumericFailure,
    InternalError
};

struct SpatialFieldCompileMetrics {
    std::uint64_t samples = 0;
    std::uint64_t controlEvaluations = 0;
    std::uint64_t growthRelaxations = 0;
    std::uint32_t growthSweeps = 0;
    std::uint64_t logicalBytes = 0;
    double minimumCompiledSize = 0.0;
    double maximumCompiledSize = 0.0;
    double measuredMaximumAdjacentGrowth = 1.0;
};

struct CompiledSpatialFieldOutcome;

class CompiledSpatialField final : public SizeField {
public:
    double eval(const Vec3& point) const override;

    const std::vector<float>& values() const noexcept { return values_; }
    Vec3 lower() const noexcept { return lower_; }
    Vec3 upper() const noexcept { return upper_; }
    std::uint32_t samplesX() const noexcept { return nx_; }
    std::uint32_t samplesY() const noexcept { return ny_; }
    std::uint32_t samplesZ() const noexcept { return nz_; }
    double baseSize() const noexcept { return baseSize_; }
    std::uint64_t deterministicHash() const noexcept { return hash_; }

private:
    friend CompiledSpatialFieldOutcome compileSpatialField(
        const CompiledSpatialFieldRequest&) noexcept;

    Vec3 lower_{};
    Vec3 upper_{1.0, 1.0, 1.0};
    std::uint32_t nx_ = 0;
    std::uint32_t ny_ = 0;
    std::uint32_t nz_ = 0;
    double baseSize_ = 1.0;
    std::vector<float> values_;
    std::uint64_t hash_ = 0;
};

struct CompiledSpatialFieldOutcome {
    bool compiled = false;
    SpatialFieldCompileReason reason =
        SpatialFieldCompileReason::InternalError;
    CompiledSpatialField field{};
    SpatialFieldCompileMetrics metrics{};
};

CompiledSpatialFieldOutcome compileSpatialField(
    const CompiledSpatialFieldRequest& request) noexcept;

}  // namespace paacs::sizing
