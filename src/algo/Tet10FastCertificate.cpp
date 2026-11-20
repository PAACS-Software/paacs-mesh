// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Tet10FastCertificate.hpp"
#include "paacs/topo/ElementTraits.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
static_assert(std::numeric_limits<double>::is_iec559,
              "Tet10 interval certification requires IEEE binary64");
namespace paacs::algo {
namespace {
struct I {
    double lo, hi;
};
constexpr double inf = std::numeric_limits<double>::infinity();
I add(I a, I b) {
    return {std::nextafter(a.lo + b.lo, -inf), std::nextafter(a.hi + b.hi, inf)};
}
I sub(I a, I b) {
    return add(a, {-b.hi, -b.lo});
}
I mul(I a, I b) {
    const std::array<double, 4> v{a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi};
    return {std::nextafter(*std::min_element(v.begin(), v.end()), -inf),
            std::nextafter(*std::max_element(v.begin(), v.end()), inf)};
}
I scale(I a, double s) {
    return mul(a, {s, s});
}
} // namespace
bool positiveTet10JacobianFast(const std::array<Vec3, 10> &nodes) noexcept {
    double largest = 0;
    for (const auto &p : nodes) {
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
            return false;
        largest = std::max({largest, std::abs(p.x - nodes[0].x), std::abs(p.y - nodes[0].y),
                            std::abs(p.z - nodes[0].z)});
    }
    if (!(largest > 0) || !std::isfinite(largest))
        return false;
    const I inverse{std::nextafter(1 / largest, -inf), std::nextafter(1 / largest, inf)};
    if (!std::isfinite(inverse.hi))
        return false;
    std::array<std::array<I, 3>, 10> p{};
    for (unsigned i = 1; i < 10; ++i) {
        const double a[3]{nodes[i].x, nodes[i].y, nodes[i].z};
        const double b[3]{nodes[0].x, nodes[0].y, nodes[0].z};
        for (unsigned d = 0; d < 3; ++d)
            p[i][d] = mul(sub({a[d], a[d]}, {b[d], b[d]}), inverse);
    }
    const auto *edges = topo::traits(topo::ElementType::Tet10).edges;
    constexpr int dl[4][3]{{-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    I j[3][3];
    for (unsigned a = 0; a < 3; ++a)
        for (unsigned d = 0; d < 3; ++d)
            j[a][d] = {inf, -inf};
    std::array<std::array<I, 3>, 6> delta{};
    for (unsigned e = 0; e < 6; ++e)
        for (unsigned a = 0; a < 3; ++a)
            delta[e][a] = sub(p[4 + e][a], scale(add(p[edges[e].a][a], p[edges[e].b][a]), .5));
    for (unsigned v = 0; v < 4; ++v)
        for (unsigned a = 0; a < 3; ++a)
            for (unsigned d = 0; d < 3; ++d) {
                I value = p[d + 1][a];
                for (unsigned e = 0; e < 6; ++e) {
                    const int derivative = (edges[e].a == v ? dl[edges[e].b][d] : 0) +
                                           (edges[e].b == v ? dl[edges[e].a][d] : 0);
                    if (derivative)
                        value = add(value, scale(delta[e][a], 4 * derivative));
                }
                j[a][d].lo = std::min(j[a][d].lo, value.lo);
                j[a][d].hi = std::max(j[a][d].hi, value.hi);
            }
    const I det = add(sub(mul(j[0][0], sub(mul(j[1][1], j[2][2]), mul(j[1][2], j[2][1]))),
                          mul(j[0][1], sub(mul(j[1][0], j[2][2]), mul(j[1][2], j[2][0])))),
                      mul(j[0][2], sub(mul(j[1][0], j[2][1]), mul(j[1][1], j[2][0]))));
    return std::isfinite(det.lo) && std::isfinite(det.hi) && det.lo > 0;
}
} // namespace paacs::algo
