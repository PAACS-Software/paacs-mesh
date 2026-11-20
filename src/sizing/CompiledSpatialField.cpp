// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sizing/CompiledSpatialField.hpp"

#include "paacs/core/Parallel.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>

namespace paacs::sizing {
namespace {

bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
bool finitePositive(double value) noexcept {
    return std::isfinite(value) && value > 0.0;
}
bool checkedMul(std::uint64_t a, std::uint64_t b,
                std::uint64_t& out) noexcept {
    if (a != 0 && b > (std::numeric_limits<std::uint64_t>::max)() / a)
        return false;
    out = a * b;
    return true;
}
std::size_t index(std::uint32_t i, std::uint32_t j, std::uint32_t k,
                  std::uint32_t nx, std::uint32_t ny) noexcept {
    return (static_cast<std::size_t>(k) * ny + j) * nx + i;
}
double segmentDistance(const Vec3& point, const Vec3& a,
                       const Vec3& b) noexcept {
    const Vec3 ab = b - a;
    const double length2 = ab.norm2();
    if (!(length2 > 0.0)) return (point - a).norm();
    const double t = std::clamp((point - a).dot(ab) / length2, 0.0, 1.0);
    return (point - (a + ab * t)).norm();
}

bool validControl(const SpatialSizeControl& c) noexcept {
    if (!finite(c.a) || !finite(c.b) || !finitePositive(c.targetSize))
        return false;
    switch (c.kind) {
        case SpatialControlKind::Box:
            return c.a.x != c.b.x && c.a.y != c.b.y && c.a.z != c.b.z;
        case SpatialControlKind::Ball:
            return finitePositive(c.radius);
        case SpatialControlKind::Cylinder:
            return finitePositive(c.radius) && (c.b - c.a).norm2() > 0.0;
        case SpatialControlKind::DistanceRamp:
            return c.radius >= 0.0 && std::isfinite(c.radius) &&
                   c.outerRadius > c.radius &&
                   std::isfinite(c.outerRadius) &&
                   finitePositive(c.outerSize);
    }
    return false;
}

double controlValue(const SpatialSizeControl& c,
                    const Vec3& point) noexcept {
    switch (c.kind) {
        case SpatialControlKind::Box: {
            const Vec3 lo{std::min(c.a.x, c.b.x),
                          std::min(c.a.y, c.b.y),
                          std::min(c.a.z, c.b.z)};
            const Vec3 hi{std::max(c.a.x, c.b.x),
                          std::max(c.a.y, c.b.y),
                          std::max(c.a.z, c.b.z)};
            return point.x >= lo.x && point.x <= hi.x &&
                   point.y >= lo.y && point.y <= hi.y &&
                   point.z >= lo.z && point.z <= hi.z
                ? c.targetSize
                : (std::numeric_limits<double>::infinity)();
        }
        case SpatialControlKind::Ball:
            return (point - c.a).norm() <= c.radius
                ? c.targetSize
                : (std::numeric_limits<double>::infinity)();
        case SpatialControlKind::Cylinder:
            return segmentDistance(point, c.a, c.b) <= c.radius
                ? c.targetSize
                : (std::numeric_limits<double>::infinity)();
        case SpatialControlKind::DistanceRamp: {
            const double distance = (point - c.a).norm();
            if (distance <= c.radius) return c.targetSize;
            if (distance >= c.outerRadius) return c.outerSize;
            const double t = (distance - c.radius) /
                             (c.outerRadius - c.radius);
            return c.targetSize + t * (c.outerSize - c.targetSize);
        }
    }
    return (std::numeric_limits<double>::infinity)();
}

void hashBytes(std::uint64_t& hash, const void* data,
               std::size_t size) noexcept {
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
        hash ^= bytes[i];
        hash *= 1099511628211ULL;
    }
}
template <class T>
void hashValue(std::uint64_t& hash, const T& value) noexcept {
    hashBytes(hash, &value, sizeof(value));
}
std::uint64_t fieldHash(const CompiledSpatialField& field) noexcept {
    std::uint64_t hash = 1469598103934665603ULL;
    const std::uint32_t version = 1;
    hashValue(hash, version);
    hashValue(hash, field.lower().x); hashValue(hash, field.lower().y);
    hashValue(hash, field.lower().z); hashValue(hash, field.upper().x);
    hashValue(hash, field.upper().y); hashValue(hash, field.upper().z);
    hashValue(hash, field.samplesX()); hashValue(hash, field.samplesY());
    hashValue(hash, field.samplesZ()); hashValue(hash, field.baseSize());
    if (!field.values().empty())
        hashBytes(hash, field.values().data(),
                  field.values().size() * sizeof(float));
    return hash;
}

}  // namespace

double CompiledSpatialField::eval(const Vec3& point) const {
    if (values_.empty() || nx_ < 2 || ny_ < 2 || nz_ < 2)
        return baseSize_;
    const Vec3 span = upper_ - lower_;
    const auto coordinate = [](double value, double lo, double length,
                               std::uint32_t samples) {
        const double normalized = std::clamp((value - lo) / length, 0.0, 1.0);
        const double grid = normalized * static_cast<double>(samples - 1U);
        const auto lowerIndex = static_cast<std::uint32_t>(
            std::min(grid, static_cast<double>(samples - 2U)));
        return std::pair{lowerIndex, grid - lowerIndex};
    };
    const auto [i, u] = coordinate(point.x, lower_.x, span.x, nx_);
    const auto [j, v] = coordinate(point.y, lower_.y, span.y, ny_);
    const auto [k, w] = coordinate(point.z, lower_.z, span.z, nz_);
    const auto sample = [&](std::uint32_t di, std::uint32_t dj,
                            std::uint32_t dk) {
        return static_cast<double>(values_[index(
            i + di, j + dj, k + dk, nx_, ny_)]);
    };
    const double c00 = sample(0,0,0) * (1.0-u) + sample(1,0,0) * u;
    const double c10 = sample(0,1,0) * (1.0-u) + sample(1,1,0) * u;
    const double c01 = sample(0,0,1) * (1.0-u) + sample(1,0,1) * u;
    const double c11 = sample(0,1,1) * (1.0-u) + sample(1,1,1) * u;
    const double c0 = c00 * (1.0-v) + c10 * v;
    const double c1 = c01 * (1.0-v) + c11 * v;
    return c0 * (1.0-w) + c1 * w;
}

CompiledSpatialFieldOutcome compileSpatialField(
        const CompiledSpatialFieldRequest& request) noexcept {
    CompiledSpatialFieldOutcome out;
    if (!finite(request.lower) || !finite(request.upper) ||
        !(request.upper.x > request.lower.x) ||
        !(request.upper.y > request.lower.y) ||
        !(request.upper.z > request.lower.z)) {
        out.reason = SpatialFieldCompileReason::InvalidBounds;
        return out;
    }
    if (request.samplesX < 2 || request.samplesY < 2 ||
        request.samplesZ < 2 || request.maximumGrowthSweeps == 0) {
        out.reason = SpatialFieldCompileReason::InvalidResolution;
        return out;
    }
    if (!finitePositive(request.baseSize) ||
        !std::isfinite(request.minimumSize) || request.minimumSize < 0.0 ||
        !finitePositive(request.maximumSize) ||
        request.minimumSize > request.maximumSize) {
        out.reason = SpatialFieldCompileReason::InvalidSize;
        return out;
    }
    if (!std::isfinite(request.maximumAdjacentGrowth) ||
        request.maximumAdjacentGrowth < 1.0) {
        out.reason = SpatialFieldCompileReason::InvalidGrowth;
        return out;
    }
    for (const auto& control : request.controls)
        if (!validControl(control)) {
            out.reason = SpatialFieldCompileReason::InvalidControl;
            return out;
        }
    std::uint64_t xy = 0, sampleCount = 0;
    if (!checkedMul(request.samplesX, request.samplesY, xy) ||
        !checkedMul(xy, request.samplesZ, sampleCount)) {
        out.reason = SpatialFieldCompileReason::CountOverflow;
        return out;
    }
    if (sampleCount > request.maximumSamples ||
        sampleCount > (std::numeric_limits<std::size_t>::max)()) {
        out.reason = SpatialFieldCompileReason::BudgetExceeded;
        return out;
    }

    try {
        auto& field = out.field;
        field.lower_ = request.lower;
        field.upper_ = request.upper;
        field.nx_ = request.samplesX;
        field.ny_ = request.samplesY;
        field.nz_ = request.samplesZ;
        field.baseSize_ = std::clamp(
            request.baseSize, request.minimumSize, request.maximumSize);
        field.values_.resize(static_cast<std::size_t>(sampleCount));
        const Vec3 span = request.upper - request.lower;
        core::parallelForWithThreads(
            0, field.values_.size(), request.executionThreads,
            [&](std::size_t ordinal) {
                const std::uint64_t row = ordinal / request.samplesX;
                const std::uint32_t i = static_cast<std::uint32_t>(
                    ordinal % request.samplesX);
                const std::uint32_t j = static_cast<std::uint32_t>(
                    row % request.samplesY);
                const std::uint32_t k = static_cast<std::uint32_t>(
                    row / request.samplesY);
                const Vec3 point{
                    request.lower.x + span.x * i / (request.samplesX - 1U),
                    request.lower.y + span.y * j / (request.samplesY - 1U),
                    request.lower.z + span.z * k / (request.samplesZ - 1U)};
                double value = field.baseSize_;
                for (const auto& control : request.controls)
                    value = std::min(value, controlValue(control, point));
                value = std::clamp(
                    value, request.minimumSize, request.maximumSize);
                field.values_[ordinal] = static_cast<float>(value);
            });

        const float ratio = static_cast<float>(
            request.maximumAdjacentGrowth);
        std::uint64_t relaxations = 0;
        std::uint32_t sweeps = 0;
        const std::size_t plane = static_cast<std::size_t>(request.samplesX) *
                                  request.samplesY;
        for (; sweeps < request.maximumGrowthSweeps; ++sweeps) {
            bool changed = false;
            for (std::uint32_t k = 0; k < request.samplesZ; ++k)
                for (std::uint32_t j = 0; j < request.samplesY; ++j)
                    for (std::uint32_t i = 0; i < request.samplesX; ++i) {
                        const std::size_t at = index(
                            i,j,k,request.samplesX,request.samplesY);
                        float candidate = field.values_[at];
                        if (i) candidate = std::min(
                            candidate, field.values_[at-1] * ratio);
                        if (j) candidate = std::min(
                            candidate, field.values_[at-request.samplesX] * ratio);
                        if (k) candidate = std::min(
                            candidate, field.values_[at-plane] * ratio);
                        if (candidate < field.values_[at]) {
                            field.values_[at] = candidate;
                            changed = true; ++relaxations;
                        }
                    }
            for (std::uint32_t k = request.samplesZ; k-- > 0;)
                for (std::uint32_t j = request.samplesY; j-- > 0;)
                    for (std::uint32_t i = request.samplesX; i-- > 0;) {
                        const std::size_t at = index(
                            i,j,k,request.samplesX,request.samplesY);
                        float candidate = field.values_[at];
                        if (i+1 < request.samplesX) candidate = std::min(
                            candidate, field.values_[at+1] * ratio);
                        if (j+1 < request.samplesY) candidate = std::min(
                            candidate, field.values_[at+request.samplesX] * ratio);
                        if (k+1 < request.samplesZ) candidate = std::min(
                            candidate, field.values_[at+plane] * ratio);
                        if (candidate < field.values_[at]) {
                            field.values_[at] = candidate;
                            changed = true; ++relaxations;
                        }
                    }
            if (!changed) { ++sweeps; break; }
        }

        float minimum = (std::numeric_limits<float>::infinity)();
        float maximum = 0.0F;
        double measured = 1.0;
        for (std::uint32_t k = 0; k < request.samplesZ; ++k)
            for (std::uint32_t j = 0; j < request.samplesY; ++j)
                for (std::uint32_t i = 0; i < request.samplesX; ++i) {
                    const std::size_t at = index(
                        i,j,k,request.samplesX,request.samplesY);
                    const float value = field.values_[at];
                    if (!std::isfinite(value) || !(value > 0.0F)) {
                        out.reason = SpatialFieldCompileReason::NumericFailure;
                        return out;
                    }
                    minimum = std::min(minimum, value);
                    maximum = std::max(maximum, value);
                    const auto measure = [&](std::size_t other) {
                        const double a = field.values_[at];
                        const double b = field.values_[other];
                        measured = std::max(measured,
                            std::max(a/b, b/a));
                    };
                    if (i+1 < request.samplesX) measure(at+1);
                    if (j+1 < request.samplesY) measure(at+request.samplesX);
                    if (k+1 < request.samplesZ) measure(at+plane);
                }
        if (measured > request.maximumAdjacentGrowth * (1.0 + 2.0e-6)) {
            out.reason = SpatialFieldCompileReason::NumericFailure;
            return out;
        }
        field.hash_ = fieldHash(field);
        out.metrics.samples = sampleCount;
        out.metrics.controlEvaluations = sampleCount * request.controls.size();
        out.metrics.growthRelaxations = relaxations;
        out.metrics.growthSweeps = sweeps;
        out.metrics.logicalBytes = sampleCount * sizeof(float);
        out.metrics.minimumCompiledSize = minimum;
        out.metrics.maximumCompiledSize = maximum;
        out.metrics.measuredMaximumAdjacentGrowth = measured;
        out.compiled = true;
        out.reason = SpatialFieldCompileReason::None;
    } catch (const std::bad_alloc&) {
        out.reason = SpatialFieldCompileReason::AllocationFailed;
    } catch (...) {
        out.reason = SpatialFieldCompileReason::InternalError;
    }
    return out;
}

}  // namespace paacs::sizing
