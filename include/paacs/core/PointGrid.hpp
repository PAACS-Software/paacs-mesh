// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Uniform spatial hash for finite points. Hash collisions are resolved by
// distance checks; extreme radii use a bounded scan of populated buckets.
#include "paacs/core/Types.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace paacs::core {
class PointGrid {
  public:
    explicit PointGrid(double cellSize) : cell_(cellSize), inv_(1.0 / cellSize) {
        if (!(cellSize > 0.0) || !std::isfinite(cellSize))
            throw std::invalid_argument("PointGrid: cell size must be finite and positive");
    }
    void insert(const Vec3 &p) {
        requireFinite(p);
        map_[key(p)].push_back(p);
        ++size_;
    }
    // Strict radius: zero never matches; negative/non-finite radii are errors.
    bool anyWithin(const Vec3 &p, double r) const {
        requireFinite(p);
        if (r < 0.0 || !std::isfinite(r))
            throw std::invalid_argument("PointGrid: radius must be finite and nonnegative");
        if (r == 0.0 || size_ == 0)
            return false;
        const double r2 = r * r;
        const auto within = [&](const Vec3 &q) {
            const Vec3 d = q - p;
            // Retain the fast squared-distance path at ordinary scales.
            // hypot avoids overflow/underflow for extreme binary64 inputs.
            return std::isnormal(r2) ? d.norm2() < r2 : std::hypot(d.x, d.y, d.z) < r;
        };
        const double spanValue = std::ceil(scaled(r));
        // Never cast an unbounded floating-point value or enumerate billions
        // of empty cells. A populated-bucket scan also handles subnormal cells.
        if (!std::isfinite(spanValue) || spanValue > 32.0) {
            for (const auto &[hash, points] : map_) {
                (void)hash;
                for (const auto &q : points)
                    if (within(q))
                        return true;
            }
            return false;
        }
        // One conservative extra cell covers rounding at grid boundaries.
        const long long span = static_cast<long long>(spanValue) + 1;
        const long long cx = coord(p.x), cy = coord(p.y), cz = coord(p.z);
        for (long long dz = -span; dz <= span; ++dz)
            for (long long dy = -span; dy <= span; ++dy)
                for (long long dx = -span; dx <= span; ++dx) {
                    const auto it = map_.find(pack(cx + dx, cy + dy, cz + dz));
                    if (it != map_.end())
                        for (const auto &q : it->second)
                            if (within(q))
                                return true;
                }
        return false;
    }
    std::size_t size() const noexcept { return size_; }

  private:
    static void requireFinite(const Vec3 &p) {
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
            throw std::invalid_argument("PointGrid: non-finite position");
    }
    double scaled(double x) const noexcept { return std::isfinite(inv_) ? x * inv_ : x / cell_; }
    long long coord(double x) const noexcept {
        constexpr double limit = 4.0e18;
        return static_cast<long long>(std::clamp(std::floor(scaled(x)), -limit, limit));
    }
    static std::uint64_t pack(long long x, long long y, long long z) noexcept {
        std::uint64_t h = static_cast<std::uint64_t>(x) * 0x9E3779B97F4A7C15ULL;
        h ^= static_cast<std::uint64_t>(y) * 0xC2B2AE3D27D4EB4FULL;
        h ^= static_cast<std::uint64_t>(z) * 0x165667B19E3779F9ULL;
        return h ^ (h >> 29);
    }
    std::uint64_t key(const Vec3 &p) const noexcept {
        return pack(coord(p.x), coord(p.y), coord(p.z));
    }
    double cell_, inv_;
    std::size_t size_ = 0;
    std::unordered_map<std::uint64_t, std::vector<Vec3>> map_;
};
} // namespace paacs::core
