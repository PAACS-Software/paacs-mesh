// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/core/Types.hpp"
#include "paacs/sizing/SizeField.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace paacs::sizing {

// Symmetric positive-definite 3D metric.  Desired metric edge length is
// sqrt(dx^T M dx).  The public representation is compact and independent of
// any geometry kernel or eigensolver implementation.
struct MetricTensor3 {
    double xx = 1.0;
    double yy = 1.0;
    double zz = 1.0;
    double xy = 0.0;
    double xz = 0.0;
    double yz = 0.0;

    static MetricTensor3 isotropic(double targetSize) noexcept;
    bool finite() const noexcept;
    bool positiveDefinite() const noexcept;
    double quadraticForm(const Vec3& direction) const noexcept;
    double metricLength(const Vec3& delta) const noexcept;
};

struct MetricEigenBounds {
    bool valid = false;
    double minimum = 0.0;
    double maximum = 0.0;
};

MetricEigenBounds metricEigenBounds(const MetricTensor3& value) noexcept;

class MetricField {
public:
    virtual ~MetricField() = default;
    virtual MetricTensor3 evalMetric(const Vec3& point) const = 0;
};

class ConstantMetricField final : public MetricField {
public:
    explicit ConstantMetricField(MetricTensor3 value) : value_(value) {}
    MetricTensor3 evalMetric(const Vec3&) const override { return value_; }
private:
    MetricTensor3 value_{};
};

// Deterministic conservative composition.  SPD tensors are added, so each
// contributor can only make the metric stricter.  This is intentionally
// simple and robust; advanced metric intersections can be added after real
// CAD/solution fields exist without changing the consumer interface.
class AdditiveMetricField final : public MetricField {
public:
    void add(std::shared_ptr<MetricField> field);
    MetricTensor3 evalMetric(const Vec3& point) const override;
    bool empty() const noexcept { return fields_.empty(); }
private:
    std::vector<std::shared_ptr<MetricField>> fields_;
};

// Adapter allowing all existing scalar meshers to consume an SPD metric
// conservatively.  It returns the smallest desired principal size,
// 1/sqrt(lambda_max), clamped to the requested scalar envelope.
class MetricSizeFieldAdapter final : public SizeField {
public:
    MetricSizeFieldAdapter(std::shared_ptr<MetricField> metric,
                           double minimumSize,
                           double maximumSize) noexcept;
    double eval(const Vec3& point) const override;
    bool valid() const noexcept { return valid_; }
private:
    std::shared_ptr<MetricField> metric_;
    double minimumSize_ = 0.0;
    double maximumSize_ = 0.0;
    bool valid_ = false;
};

struct DirectionalMetricControl {
    Vec3 axis{1.0, 0.0, 0.0};
    double alongSize = 1.0;
    double transverseSize = 1.0;
};

// Creates an SPD metric whose principal axis follows control.axis.  The two
// perpendicular directions use transverseSize.
MetricTensor3 directionalMetric(
    const DirectionalMetricControl& control) noexcept;

}  // namespace paacs::sizing
