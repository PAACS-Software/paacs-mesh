// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/sizing/BoundedSizeField.hpp"
#include <cstdint>
#include <span>

namespace paacs::sizing {
enum class SizeRegionKind : std::uint8_t { Box, Ball, Rectangle };
enum class BoxFace : std::uint8_t { XMin, XMax, YMin, YMax, ZMin, ZMax };

struct RegionSizeControl {
    SizeRegionKind kind = SizeRegionKind::Box;
    // Ordered bounds for Box/Rectangle; lower is the center for Ball.
    Vec3 lower{}, upper{1, 1, 1};
    double radius = 1;
    double targetSize = .1;
    // Maps to a physical size slope 2*(g-1)/(g+1). This is a field
    // control, NOT a promise about every neighboring element's edge lengths.
    double growthRatio = 1.3;
};
struct RegionSizeSpecification {
    double baseSize = .25;
    double minimumSize = 0;
    std::vector<RegionSizeControl> controls;
    std::uint64_t maximumControls = 1024;
    std::uint64_t maximumBytes = 1024 * 1024;
};
enum class RegionSizeReason : std::uint8_t {
    None,
    InvalidSize,
    InvalidRegion,
    InvalidGrowth,
    BudgetExceeded,
    AllocationFailed,
    InternalError
};
struct RegionSizeOutcome;
class RegionSizeField final : public BoundedSizeField {
  public:
    double eval(const Vec3 &point) const override;
    // Exact minimum of the analytic field over an axis-aligned box.
    // Unlike point sampling, this cannot miss a small enclosed region.
    double minimumOnBox(const Vec3 &lower, const Vec3 &upper) const noexcept override;
    bool valid() const noexcept override { return valid_; }
    std::uint64_t boxQueryWorkBound() const noexcept override { return 1 + controls_.size(); }
    double baseSize() const noexcept { return baseSize_; }
    double minimumSize() const noexcept { return minimumSize_; }
    double gradientBound() const noexcept { return gradientBound_; }
    std::span<const RegionSizeControl> controls() const noexcept { return controls_; }

  private:
    friend RegionSizeOutcome compileRegionSizeField(const RegionSizeSpecification &) noexcept;
    bool valid_ = false;
    double baseSize_ = .25, minimumSize_ = 0, gradientBound_ = 0;
    std::vector<RegionSizeControl> controls_;
};
struct RegionSizeOutcome {
    bool compiled = false;
    RegionSizeReason reason = RegionSizeReason::InternalError;
    std::uint64_t bytes = 0;
    RegionSizeField field;
};
RegionSizeOutcome compileRegionSizeField(const RegionSizeSpecification &) noexcept;
// Invalid face enum produces an invalid control, rejected at compilation.
RegionSizeControl boxFaceSizeControl(Vec3 lower, Vec3 upper, BoxFace face, double size,
                                     double growth = 1.3) noexcept;
} // namespace paacs::sizing
