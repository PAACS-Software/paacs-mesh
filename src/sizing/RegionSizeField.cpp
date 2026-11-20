// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sizing/RegionSizeField.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <tuple>

namespace paacs::sizing {
namespace {
bool finite(Vec3 p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
bool bounds(Vec3 lo, Vec3 hi) noexcept {
    return finite(lo) && finite(hi) && lo.x <= hi.x && lo.y <= hi.y && lo.z <= hi.z;
}
double slope(double g) noexcept {
    return 2 * (1 - 2 / (g + 1));
}
double separation(double a, double b, double c, double d) noexcept {
    return std::max({0., a - d, c - b});
}
double boxDistance(Vec3 lo, Vec3 hi, Vec3 a, Vec3 b) noexcept {
    return std::hypot(separation(lo.x, hi.x, a.x, b.x), separation(lo.y, hi.y, a.y, b.y),
                      separation(lo.z, hi.z, a.z, b.z));
}
double distance(const RegionSizeControl &c, Vec3 lo, Vec3 hi) noexcept {
    if (c.kind == SizeRegionKind::Ball)
        return std::max(0., boxDistance(lo, hi, c.lower, c.lower) - c.radius);
    return boxDistance(lo, hi, c.lower, c.upper);
}
auto key(const RegionSizeControl &c) noexcept {
    return std::tuple{c.kind,    c.lower.x, c.lower.y, c.lower.z,    c.upper.x,
                      c.upper.y, c.upper.z, c.radius,  c.targetSize, c.growthRatio};
}
} // namespace

double RegionSizeField::minimumOnBox(const Vec3 &lo, const Vec3 &hi) const noexcept {
    if (!valid_ || !bounds(lo, hi))
        return std::numeric_limits<double>::quiet_NaN();
    double value = baseSize_;
    for (const auto &c : controls_) {
        const double h = std::max(minimumSize_, c.targetSize), s = slope(c.growthRatio);
        if (s == 0) {
            value = std::min(value, h);
            continue;
        }
        const double d = distance(c, lo, hi);
        // Saturate before multiplication, including overflowed distances.
        if (d < (value - h) / s)
            value = h + s * d;
    }
    return value;
}
double RegionSizeField::eval(const Vec3 &p) const {
    return minimumOnBox(p, p);
}

RegionSizeOutcome compileRegionSizeField(const RegionSizeSpecification &request) noexcept {
    RegionSizeOutcome out;
    if (!std::isfinite(request.baseSize) || !(request.baseSize > 0) ||
        !std::isfinite(request.minimumSize) || request.minimumSize < 0 ||
        request.minimumSize > request.baseSize) {
        out.reason = RegionSizeReason::InvalidSize;
        return out;
    }
    if (request.controls.size() > request.maximumControls ||
        request.controls.size() > request.maximumBytes / sizeof(RegionSizeControl)) {
        out.reason = RegionSizeReason::BudgetExceeded;
        return out;
    }
    for (const auto &c : request.controls) {
        if (!std::isfinite(c.targetSize) || !(c.targetSize > 0) ||
            c.targetSize > request.baseSize) {
            out.reason = RegionSizeReason::InvalidSize;
            return out;
        }
        if (!std::isfinite(c.growthRatio) || c.growthRatio < 1) {
            out.reason = RegionSizeReason::InvalidGrowth;
            return out;
        }
        bool valid = false;
        switch (c.kind) {
        case SizeRegionKind::Ball:
            valid = finite(c.lower) && std::isfinite(c.radius) && c.radius > 0;
            break;
        case SizeRegionKind::Box:
            valid = bounds(c.lower, c.upper) && c.lower.x < c.upper.x && c.lower.y < c.upper.y &&
                    c.lower.z < c.upper.z;
            break;
        case SizeRegionKind::Rectangle:
            valid = bounds(c.lower, c.upper) &&
                    (int(c.lower.x == c.upper.x) + int(c.lower.y == c.upper.y) +
                         int(c.lower.z == c.upper.z) ==
                     1);
            break;
        }
        if (!valid) {
            out.reason = RegionSizeReason::InvalidRegion;
            return out;
        }
    }
    try {
        auto &f = out.field;
        f.baseSize_ = request.baseSize;
        f.minimumSize_ = request.minimumSize;
        f.controls_ = request.controls;
        for (auto &c : f.controls_) {
            if (c.kind == SizeRegionKind::Ball)
                c.upper = c.lower;
            else
                c.radius = 0;
            f.gradientBound_ = std::max(f.gradientBound_, slope(c.growthRatio));
        }
        std::sort(f.controls_.begin(), f.controls_.end(),
                  [](const auto &a, const auto &b) { return key(a) < key(b); });
        f.valid_ = true;
        out.bytes = f.controls_.capacity() * sizeof(RegionSizeControl);
        out.compiled = true;
        out.reason = RegionSizeReason::None;
    } catch (const std::bad_alloc &) {
        out.reason = RegionSizeReason::AllocationFailed;
    } catch (...) {
        out.reason = RegionSizeReason::InternalError;
    }
    return out;
}
RegionSizeControl boxFaceSizeControl(Vec3 lo, Vec3 hi, BoxFace face, double size,
                                     double growth) noexcept {
    RegionSizeControl c;
    c.kind = SizeRegionKind::Rectangle;
    c.lower = lo;
    c.upper = hi;
    c.targetSize = size;
    c.growthRatio = growth;
    if (!bounds(lo, hi) || lo.x == hi.x || lo.y == hi.y || lo.z == hi.z) {
        c.targetSize = -1;
        return c;
    }
    switch (face) {
    case BoxFace::XMin:
        c.upper.x = lo.x;
        break;
    case BoxFace::XMax:
        c.lower.x = hi.x;
        break;
    case BoxFace::YMin:
        c.upper.y = lo.y;
        break;
    case BoxFace::YMax:
        c.lower.y = hi.y;
        break;
    case BoxFace::ZMin:
        c.upper.z = lo.z;
        break;
    case BoxFace::ZMax:
        c.lower.z = hi.z;
        break;
    default:
        c.targetSize = -1;
        break;
    }
    return c;
}
} // namespace paacs::sizing
