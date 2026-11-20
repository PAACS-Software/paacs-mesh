// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/mesh/Hex8Jacobian.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace paacs::mesh {
// det(J) of a trilinear Hex8 is degree <=2 in each coordinate. Convert
// evaluations on {0,.5,1}^3 to tensor Bernstein coefficients. Their convex
// hull bounds the complete polynomial. A deliberately conservative error
// guard covers subtract/interpolate/determinant/conversion roundoff; there
// is no subdivision, so a failed sufficient test means "unresolved".
bool positiveHex8Jacobian(const std::array<Vec3, 8> &p) noexcept {
    constexpr std::array<unsigned, 8> order{0, 1, 3, 2, 4, 5, 7, 6};
    std::array<Vec3, 8> q{};
    double scale = 0;
    for (unsigned i = 0; i < 8; ++i) {
        q[i] = p[order[i]] - p[0];
        scale = std::max({scale, std::abs(q[i].x), std::abs(q[i].y), std::abs(q[i].z)});
    }
    if (!(scale > 0) || !std::isfinite(scale))
        return false;
    // Normalize to keep the determinant test away from overflow/underflow.
    for (auto &v : q)
        v = v / scale;
    std::array<double, 27> b{};
    auto index = [](unsigned i, unsigned j, unsigned k) { return i + 3 * j + 9 * k; };
    for (unsigned k = 0; k < 3; ++k)
        for (unsigned j = 0; j < 3; ++j)
            for (unsigned i = 0; i < 3; ++i) {
                const double u = .5 * i, v = .5 * j, w = .5 * k;
                Vec3 du{}, dv{}, dw{};
                for (unsigned z = 0; z < 2; ++z)
                    for (unsigned y = 0; y < 2; ++y)
                        du = du + (q[1 + 2 * y + 4 * z] - q[2 * y + 4 * z]) *
                                      ((y ? v : 1 - v) * (z ? w : 1 - w));
                for (unsigned z = 0; z < 2; ++z)
                    for (unsigned x = 0; x < 2; ++x)
                        dv = dv + (q[x + 2 + 4 * z] - q[x + 4 * z]) *
                                      ((x ? u : 1 - u) * (z ? w : 1 - w));
                for (unsigned y = 0; y < 2; ++y)
                    for (unsigned x = 0; x < 2; ++x)
                        dw = dw + (q[x + 2 * y + 4] - q[x + 2 * y]) *
                                      ((x ? u : 1 - u) * (y ? v : 1 - v));
                b[index(i, j, k)] = du.cross(dv).dot(dw);
            }
    for (unsigned axis = 0; axis < 3; ++axis) {
        const unsigned stride = axis == 0 ? 1 : axis == 1 ? 3 : 9;
        for (unsigned k = 0; k < 3; ++k)
            for (unsigned j = 0; j < 3; ++j)
                for (unsigned i = 0; i < 3; ++i)
                    if ((axis == 0 ? i : axis == 1 ? j : k) == 0) {
                        const auto o = index(i, j, k);
                        b[o + stride] = 2 * b[o + stride] - .5 * (b[o] + b[o + 2 * stride]);
                    }
    }
    constexpr double guard = 262144 * std::numeric_limits<double>::epsilon();
    return std::all_of(b.begin(), b.end(), [](double v) { return std::isfinite(v) && v > guard; });
}
} // namespace paacs::mesh
