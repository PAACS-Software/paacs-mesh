// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sizing/MetricField.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace paacs::sizing {
namespace {

bool positiveFinite(double value) noexcept {
    return std::isfinite(value) && value > 0.0;
}

}  // namespace

MetricTensor3 MetricTensor3::isotropic(double targetSize) noexcept {
    if (!positiveFinite(targetSize)) return {};
    const double inverseSquared = 1.0 / (targetSize * targetSize);
    return {inverseSquared, inverseSquared, inverseSquared, 0.0, 0.0, 0.0};
}

bool MetricTensor3::finite() const noexcept {
    return std::isfinite(xx) && std::isfinite(yy) && std::isfinite(zz) &&
           std::isfinite(xy) && std::isfinite(xz) && std::isfinite(yz);
}

bool MetricTensor3::positiveDefinite() const noexcept {
    if (!finite() || !(xx > 0.0)) return false;
    const double l10 = xy / std::sqrt(xx);
    const double l20 = xz / std::sqrt(xx);
    const double d11 = yy - l10 * l10;
    if (!(d11 > 0.0) || !std::isfinite(d11)) return false;
    const double l21 = (yz - l20 * l10) / std::sqrt(d11);
    const double d22 = zz - l20 * l20 - l21 * l21;
    return d22 > 0.0 && std::isfinite(d22);
}

double MetricTensor3::quadraticForm(const Vec3& v) const noexcept {
    return xx * v.x * v.x + yy * v.y * v.y + zz * v.z * v.z +
           2.0 * xy * v.x * v.y + 2.0 * xz * v.x * v.z +
           2.0 * yz * v.y * v.z;
}

double MetricTensor3::metricLength(const Vec3& delta) const noexcept {
    const double squared = quadraticForm(delta);
    return squared >= 0.0 && std::isfinite(squared)
               ? std::sqrt(squared)
               : (std::numeric_limits<double>::quiet_NaN)();
}

MetricEigenBounds metricEigenBounds(const MetricTensor3& value) noexcept {
    MetricEigenBounds out;
    if (!value.positiveDefinite()) return out;
    double a[3][3] = {
        {value.xx, value.xy, value.xz},
        {value.xy, value.yy, value.yz},
        {value.xz, value.yz, value.zz}};

    // Fixed-iteration Jacobi diagonalization.  The loop count and pivot
    // ordering are deterministic and sufficient for metric validation and
    // conservative scalar adaptation; topology does not depend on iteration
    // convergence details.
    for (int sweep = 0; sweep < 12; ++sweep) {
        for (const auto pair : std::array<std::array<int, 2>, 3>{
                 std::array<int, 2>{0, 1},
                 std::array<int, 2>{0, 2},
                 std::array<int, 2>{1, 2}}) {
            const int p = pair[0];
            const int q = pair[1];
            const double apq = a[p][q];
            if (apq == 0.0) continue;
            const double tau = (a[q][q] - a[p][p]) / (2.0 * apq);
            const double t = (tau >= 0.0 ? 1.0 : -1.0) /
                (std::abs(tau) + std::sqrt(1.0 + tau * tau));
            const double c = 1.0 / std::sqrt(1.0 + t * t);
            const double s = t * c;
            const double app = a[p][p];
            const double aqq = a[q][q];
            a[p][p] = app - t * apq;
            a[q][q] = aqq + t * apq;
            a[p][q] = 0.0;
            a[q][p] = 0.0;
            for (int r = 0; r < 3; ++r) {
                if (r == p || r == q) continue;
                const double arp = a[r][p];
                const double arq = a[r][q];
                a[r][p] = c * arp - s * arq;
                a[p][r] = a[r][p];
                a[r][q] = s * arp + c * arq;
                a[q][r] = a[r][q];
            }
        }
    }
    out.minimum = std::min({a[0][0], a[1][1], a[2][2]});
    out.maximum = std::max({a[0][0], a[1][1], a[2][2]});
    out.valid = positiveFinite(out.minimum) &&
                positiveFinite(out.maximum) &&
                out.maximum >= out.minimum;
    return out;
}

void AdditiveMetricField::add(std::shared_ptr<MetricField> field) {
    if (field) fields_.push_back(std::move(field));
}

MetricTensor3 AdditiveMetricField::evalMetric(const Vec3& point) const {
    MetricTensor3 out{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (const auto& field : fields_) {
        const MetricTensor3 value = field->evalMetric(point);
        if (!value.positiveDefinite())
            return {-(std::numeric_limits<double>::infinity)(), 0.0, 0.0,
                    0.0, 0.0, 0.0};
        out.xx += value.xx;
        out.yy += value.yy;
        out.zz += value.zz;
        out.xy += value.xy;
        out.xz += value.xz;
        out.yz += value.yz;
    }
    return out;
}

MetricSizeFieldAdapter::MetricSizeFieldAdapter(
        std::shared_ptr<MetricField> metric,
        double minimumSize,
        double maximumSize) noexcept
    : metric_(std::move(metric)),
      minimumSize_(minimumSize),
      maximumSize_(maximumSize),
      valid_(metric_ && positiveFinite(minimumSize_) &&
             positiveFinite(maximumSize_) && minimumSize_ <= maximumSize_) {}

double MetricSizeFieldAdapter::eval(const Vec3& point) const {
    if (!valid_) return (std::numeric_limits<double>::quiet_NaN)();
    const auto bounds = metricEigenBounds(metric_->evalMetric(point));
    if (!bounds.valid) return (std::numeric_limits<double>::quiet_NaN)();
    const double size = 1.0 / std::sqrt(bounds.maximum);
    return std::clamp(size, minimumSize_, maximumSize_);
}

MetricTensor3 directionalMetric(
        const DirectionalMetricControl& control) noexcept {
    if (!positiveFinite(control.alongSize) ||
        !positiveFinite(control.transverseSize))
        return {};
    Vec3 axis = control.axis;
    const double length = axis.norm();
    if (!positiveFinite(length)) return {};
    axis = axis / length;
    const double transverse =
        1.0 / (control.transverseSize * control.transverseSize);
    const double along = 1.0 / (control.alongSize * control.alongSize);
    const double delta = along - transverse;
    return {
        transverse + delta * axis.x * axis.x,
        transverse + delta * axis.y * axis.y,
        transverse + delta * axis.z * axis.z,
        delta * axis.x * axis.y,
        delta * axis.x * axis.z,
        delta * axis.y * axis.z};
}

}  // namespace paacs::sizing
