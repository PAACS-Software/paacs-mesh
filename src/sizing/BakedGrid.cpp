// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sizing/SizeField.hpp"
#include "paacs/core/Parallel.hpp"

#include <algorithm>
#include <cmath>

namespace paacs::sizing {

BakedGridField::BakedGridField(const SizeField& src, const Vec3& lo,
                               const Vec3& hi, int nx, int ny, int nz)
    : lo_(lo), hi_(hi),
      nx_(std::max(2, nx)), ny_(std::max(2, ny)), nz_(std::max(2, nz)) {
    data_.resize(static_cast<size_t>(nx_) * ny_ * nz_);
    const double dx = (hi_.x - lo_.x) / (nx_ - 1);
    const double dy = (hi_.y - lo_.y) / (ny_ - 1);
    const double dz = (hi_.z - lo_.z) / (nz_ - 1);

    // Embarrassingly parallel bake over grid nodes.
    core::parallelFor(0, data_.size(), [&](size_t idx) {
        const int i = static_cast<int>(idx % static_cast<size_t>(nx_));
        const int j = static_cast<int>((idx / static_cast<size_t>(nx_)) %
                                       static_cast<size_t>(ny_));
        const int k = static_cast<int>(idx / (static_cast<size_t>(nx_) * ny_));
        data_[idx] = static_cast<float>(src.eval(
            {lo_.x + i * dx, lo_.y + j * dy, lo_.z + k * dz}));
    });
}

double BakedGridField::eval(const Vec3& p) const {
    // Trilinear interpolation, clamped to the grid.
    auto cell = [](double x, double lo, double hi, int n, int& i0, double& f) {
        const double t = (x - lo) / (hi - lo) * (n - 1);
        const double tc = std::clamp(t, 0.0, static_cast<double>(n - 1));
        i0 = std::min(n - 2, static_cast<int>(tc));
        f = tc - i0;
    };
    int i, j, k;
    double fx, fy, fz;
    cell(p.x, lo_.x, hi_.x, nx_, i, fx);
    cell(p.y, lo_.y, hi_.y, ny_, j, fy);
    cell(p.z, lo_.z, hi_.z, nz_, k, fz);

    auto at = [&](int a, int b, int c) {
        return static_cast<double>(
            data_[static_cast<size_t>(a) +
                  static_cast<size_t>(nx_) *
                      (static_cast<size_t>(b) +
                       static_cast<size_t>(ny_) * static_cast<size_t>(c))]);
    };
    const double c00 = at(i, j, k)     * (1 - fx) + at(i + 1, j, k)     * fx;
    const double c10 = at(i, j + 1, k) * (1 - fx) + at(i + 1, j + 1, k) * fx;
    const double c01 = at(i, j, k + 1)     * (1 - fx) + at(i + 1, j, k + 1)     * fx;
    const double c11 = at(i, j + 1, k + 1) * (1 - fx) + at(i + 1, j + 1, k + 1) * fx;
    const double c0 = c00 * (1 - fy) + c10 * fy;
    const double c1 = c01 * (1 - fy) + c11 * fy;
    return c0 * (1 - fz) + c1 * fz;
}

} // namespace paacs::sizing
