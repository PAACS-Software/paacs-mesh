// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include <atomic>
#include <cstdint>
#include <span>
#include <vector>

namespace paacs::sizing {
// Physical coordinates along a finite sweep. begin == end selects a plane.
struct ProfileSizeControl {
    double begin = 0, end = 0, size = .1, growth = 1.3;
};
struct CappedSizeSpecification {
    double length = 1, globalSize = .25, minimumSize = 0;
    std::vector<ProfileSizeControl> controls;
    std::uint64_t maximumControls = 1024, maximumKnots = 16384;
    std::uint64_t maximumWork = 4'000'000, maximumBytes = 4ULL * 1024 * 1024;
    const std::atomic<bool> *cancellation = nullptr;
};
enum class ProfileReason : std::uint8_t {
    None,
    InvalidRequest,
    BudgetExceeded,
    Cancelled,
    NumericallyUnresolved,
    AllocationFailed,
    InternalError
};
struct ProfileKnot {
    double position = 0, size = 0;
};
struct CappedSizeOutcome;
class CappedSizeProfile {
  public:
    bool valid() const noexcept { return valid_ && knots_.size() >= 2; }
    double length() const noexcept { return length_; }
    double globalSize() const noexcept { return global_; }
    double minimumSize() const noexcept { return minimum_; }
    double gradientBound() const noexcept { return gradient_; }
    double maximumGrowth() const noexcept { return growth_; }
    std::span<const ProfileKnot> knots() const noexcept { return knots_; }
    // Invalid profiles, nonfinite arguments and out-of-domain queries return NaN.
    double eval(double position) const noexcept;
    double minimumOnInterval(double begin, double end) const noexcept;

  private:
    bool valid_ = false;
    double length_ = 0, global_ = 0, minimum_ = 0, gradient_ = 0, growth_ = 1;
    std::vector<ProfileKnot> knots_;
    friend CappedSizeOutcome compileCappedSizeProfile(const CappedSizeSpecification &) noexcept;
};
struct CappedSizeOutcome {
    bool compiled = false;
    ProfileReason reason = ProfileReason::InternalError;
    std::uint64_t work = 0, workingBytes = 0;
    CappedSizeProfile profile;
};
// h(z) = min(global, min_i(max(minimumSize,size_i) + s_i*distance(z,[a_i,b_i]))),
// s_i = 2*(g_i-1)/(g_i+1). Deterministic lower envelope, independent of a grid.
CappedSizeOutcome compileCappedSizeProfile(const CappedSizeSpecification &) noexcept;
struct ProfileDistributionLimits {
    std::uint64_t maximumCells = 100'000, maximumBytes = 16ULL * 1024 * 1024;
    // Zero uses the largest control growth. One requests a uniform fine grid.
    double maximumAdjacentGrowth = 0;
    const std::atomic<bool> *cancellation = nullptr;
};
struct ProfileDistribution {
    bool planned = false;
    ProfileReason reason = ProfileReason::InternalError;
    std::vector<double> positions; // Strictly increasing, exactly 0 ... length.
    double minimumWidth = 0, maximumWidth = 0, maximumAdjacentGrowth = 1;
    double maximumWidthToTarget = 0;
    std::uint64_t workingBytes = 0, plateauCellsRemoved = 0;
};
// Every interval width <= the minimum target over that entire interval.
// Adjacent physical widths satisfy the requested ratio (relative slack 1e-10).
// Constant runs are compacted without violating either bound. Budget preflight
// is conservative and precedes compaction; no allocation occurs during it.
ProfileDistribution planProfileDistribution(const CappedSizeProfile &,
                                            const ProfileDistributionLimits & = {}) noexcept;
const char *profileReasonName(ProfileReason) noexcept;
} // namespace paacs::sizing
