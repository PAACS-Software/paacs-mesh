// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "primitive_checks.hpp"
namespace primitive_test {
struct SizeStats {
    double cv, maximumRatio, minimumQuality;
};
// Independent long-double determinants, rather than using production metrics.
inline SizeStats sizes(const mesh::Mesh &m) {
    long double sum = 0, square = 0, low = std::numeric_limits<long double>::infinity(), high = 0;
    double minq = 1;
    std::uint64_t count = 0;
    for (Index i = 0; i < m.cellCount(); ++i) {
        if (m.cell(i).onEntity.dim != 3)
            continue;
        const auto *v = m.cellNodes(i);
        long double p[4][3]{};
        for (unsigned j = 0; j < 4; ++j) {
            const auto a = m.node(v[j]), b = m.node(v[0]);
            p[j][0] = static_cast<long double>(a.x) - b.x;
            p[j][1] = static_cast<long double>(a.y) - b.y;
            p[j][2] = static_cast<long double>(a.z) - b.z;
        }
        const auto volume = (p[1][0] * (p[2][1] * p[3][2] - p[2][2] * p[3][1]) -
                             p[1][1] * (p[2][0] * p[3][2] - p[2][2] * p[3][0]) +
                             p[1][2] * (p[2][0] * p[3][1] - p[2][1] * p[3][0])) /
                            6;
        require(volume > 0, "independent positive determinant");
        long double edgeSum = 0;
        for (unsigned a = 0; a < 4; ++a)
            for (unsigned b = 0; b < a; ++b)
                for (unsigned d = 0; d < 3; ++d)
                    edgeSum += (p[a][d] - p[b][d]) * (p[a][d] - p[b][d]);
        minq = std::min(minq, static_cast<double>(12 * std::pow(3 * volume, 2.L / 3) / edgeSum));
        sum += volume;
        square += volume * volume;
        low = std::min(low, volume);
        high = std::max(high, volume);
        ++count;
    }
    const auto mean = sum / count;
    return {static_cast<double>(std::sqrt(std::max(0.L, square / count - mean * mean)) / mean),
            static_cast<double>(high / low), minq};
}
} // namespace primitive_test
