// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// CK-2 constrained-v2 sizing boundary.
//
// The legacy meshers retain their historical unchecked MeshSpec behavior.
// meshModelConstrained()/meshFacesConstrained() replace the caller field with
// CheckedSizeField, and the three constrained construction phases route every
// actual sizing evaluation through the helpers below.  SizingEvaluationError
// is deliberately a small allocation-free exception token so an invalid user
// result cannot be confused with std::bad_alloc.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <memory>
#include <new>
#include <utility>

#include "paacs/sizing/SizeField.hpp"

namespace paacs::algo::detail {

struct SizingEvaluationError final {};

inline bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) &&
           std::isfinite(p.z);
}

inline bool positiveFinite(double value) noexcept {
    return std::isfinite(value) && value > 0.0;
}

class CheckedSizeField final : public sizing::SizeField {
public:
    explicit CheckedSizeField(std::shared_ptr<sizing::SizeField> source)
        noexcept
        : source_(std::move(source)) {}

    double eval(const Vec3& point) const override {
        if (!finite(point)) fail();
        double value = 1.0;
        try {
            if (source_) value = source_->eval(point);
        } catch (const std::bad_alloc&) {
            throw;
        } catch (const SizingEvaluationError&) {
            throw;
        } catch (...) {
            fail();
        }
        if (!positiveFinite(value)) fail();
        return value;
    }

    bool failed() const noexcept {
        return failed_.load(std::memory_order_relaxed);
    }

    [[noreturn]] void fail() const {
        failed_.store(true, std::memory_order_relaxed);
        throw SizingEvaluationError{};
    }

private:
    std::shared_ptr<sizing::SizeField> source_;
    mutable std::atomic<bool> failed_{false};
};

inline bool usesCheckedSizing(const sizing::MeshSpec& spec) noexcept {
    return dynamic_cast<const CheckedSizeField*>(spec.field.get()) != nullptr;
}

inline const CheckedSizeField* checkedField(
    const sizing::MeshSpec& spec) noexcept {
    return dynamic_cast<const CheckedSizeField*>(spec.field.get());
}

inline bool sizingEvaluationFailed(const sizing::MeshSpec& spec) noexcept {
    const auto* field = checkedField(spec);
    return field != nullptr && field->failed();
}

[[noreturn]] inline void failSizing(const sizing::MeshSpec& spec) {
    if (const auto* field = checkedField(spec)) field->fail();
    throw SizingEvaluationError{};
}

inline sizing::MeshSpec makeCheckedSpec(const sizing::MeshSpec& source) {
    if (usesCheckedSizing(source)) return source;
    sizing::MeshSpec checked = source;
    checked.field = std::make_shared<CheckedSizeField>(source.field);
    return checked;
}

// Preserve the exact legacy expression whenever a normal MeshSpec reaches a
// legacy entry.  The checked branch validates both the base value (inside the
// wrapper) and the order-adjusted value returned by MeshSpec::h().
inline double evaluatedSize(const sizing::MeshSpec& spec,
                            const Vec3& point) {
    if (!usesCheckedSizing(spec)) return spec.h(point);
    if (!finite(point)) failSizing(spec);
    const double value = spec.h(point);
    if (!positiveFinite(value)) failSizing(spec);
    return value;
}

// Compute the curvature cap in the same operation order as MeshSpec::hGeom,
// but validate every intermediate before it can feed a division.
inline double evaluatedGeometrySize(const sizing::MeshSpec& spec,
                                    const Vec3& point, double curvature) {
    if (!usesCheckedSizing(spec)) return spec.hGeom(point, curvature);
    if (!finite(point) || !std::isfinite(curvature)) failSizing(spec);
    double value = evaluatedSize(spec, point);
    if (spec.elementsPerCircle > 0 && curvature > 1.0e-12) {
        double count = static_cast<double>(spec.elementsPerCircle);
        if (spec.order >= 2) count /= spec.orderFactor;
        if (!positiveFinite(count)) failSizing(spec);
        const double denominator = count * curvature;
        if (!positiveFinite(denominator)) failSizing(spec);
        const double cap = 2.0 * 3.14159265358979323846 / denominator;
        if (!positiveFinite(cap)) failSizing(spec);
        value = std::min(value, cap);
    }
    if (!positiveFinite(value)) failSizing(spec);
    return value;
}

inline double checkedPositiveProduct(double a, double b) {
    if (!positiveFinite(a) || !positiveFinite(b))
        throw SizingEvaluationError{};
    const double value = a * b;
    if (!positiveFinite(value)) throw SizingEvaluationError{};
    return value;
}

inline double checkedPositiveRatio(double numerator, double denominator) {
    if (!positiveFinite(numerator) || !positiveFinite(denominator))
        throw SizingEvaluationError{};
    const double value = numerator / denominator;
    if (!positiveFinite(value)) throw SizingEvaluationError{};
    return value;
}

inline void requireProgress(double current, double increment) {
    if (!std::isfinite(current) || !positiveFinite(increment))
        throw SizingEvaluationError{};
    const double next = current + increment;
    if (!std::isfinite(next) || !(next > current))
        throw SizingEvaluationError{};
}

inline void requireProgress(const sizing::MeshSpec& spec, double current,
                            double increment) {
    if (!std::isfinite(current) || !positiveFinite(increment))
        failSizing(spec);
    const double next = current + increment;
    if (!std::isfinite(next) || !(next > current)) failSizing(spec);
}

}  // namespace paacs::algo::detail
