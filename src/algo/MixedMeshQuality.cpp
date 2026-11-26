// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MixedMeshQuality.hpp"

#include "paacs/mesh/Quality.hpp"
#include "paacs/predicates/Predicates.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <new>
#include <numeric>
#include <vector>

namespace paacs::algo {
namespace {
constexpr std::array<std::array<std::uint8_t, 4>, 2> kPyramid{{{0, 1, 2, 4}, {0, 2, 3, 4}}};
constexpr std::array<std::array<std::uint8_t, 4>, 3> kWedge{
    {{0, 1, 2, 3}, {1, 4, 2, 3}, {2, 4, 5, 3}}};
constexpr std::array<std::array<std::uint8_t, 4>, 6> kHex{
    {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}}};
struct CellQuality {
    bool supported = false, valid = false;
    double quality = 0, volume = 0;
};
bool volumeType(mesh::CellType t) noexcept {
    return t == mesh::CellType::Tet4 || t == mesh::CellType::Pyramid5 ||
           t == mesh::CellType::Prism6 || t == mesh::CellType::Hex8;
}
double signedTet(const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &d) noexcept {
    return mesh::tetVolume(a, b, c, d);
}
double tetMeanRatio(const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &d,
                    double volume) noexcept {
    const double sum = (b - a).norm2() + (c - a).norm2() + (d - a).norm2() + (c - b).norm2() +
                       (d - b).norm2() + (d - c).norm2();
    if (!(sum > 0) || !(volume > 0))
        return 0.0;
    return 12.0 * std::pow(3.0 * volume, 2.0 / 3.0) / sum;
}
template <std::size_t N, std::size_t M>
CellQuality decompose(const std::array<Vec3, N> &p,
                      const std::array<std::array<std::uint8_t, 4>, M> &tets) noexcept {
    CellQuality out;
    out.supported = true;
    out.valid = true;
    out.quality = (std::numeric_limits<double>::infinity)();
    for (const auto &t : tets) {
        const auto sign = predicates::orient3d(p[t[0]], p[t[1]], p[t[2]], p[t[3]]);
        double v = signedTet(p[t[0]], p[t[1]], p[t[2]], p[t[3]]);
        // A positive exact determinant can cancel to zero or a negative value
        // in the ordinary metric calculation. Recover its magnitude without
        // changing the exact orientation gate or accepting flat/inverted cells.
        if (sign == predicates::Sign::Positive && std::isfinite(v) && !(v > 0))
            v = predicates::orient3dDeterminantEstimate(
                    p[t[0]], p[t[1]], p[t[2]], p[t[3]]) / 6.0;
        if (sign != predicates::Sign::Positive || !(v > 0) || !std::isfinite(v)) {
            out.valid = false;
            out.quality = 0;
            return out;
        }
        out.volume += v;
        out.quality = std::min(out.quality, tetMeanRatio(p[t[0]], p[t[1]], p[t[2]], p[t[3]], v));
    }
    return out;
}
CellQuality quality(const mesh::Mesh &mesh, Index ci) noexcept {
    const auto &c = mesh.cell(ci);
    const Index *n = mesh.cellNodes(ci);
    const int count = mesh::nodesPerCell(c.type);
    if (!volumeType(c.type))
        return {};
    for (int i = 0; i < count; ++i)
        if (n[i] < 0 || n[i] >= mesh.nodeCount())
            return {true, false, 0, 0};
    if (c.type == mesh::CellType::Tet4) {
        std::array<Vec3, 4> p{};
        for (int i = 0; i < 4; ++i)
            p[i] = mesh.node(n[i]);
        constexpr std::array<std::array<std::uint8_t, 4>, 1> one{{{0, 1, 2, 3}}};
        return decompose(p, one);
    }
    if (c.type == mesh::CellType::Pyramid5) {
        std::array<Vec3, 5> p{};
        for (int i = 0; i < 5; ++i)
            p[i] = mesh.node(n[i]);
        return decompose(p, kPyramid);
    }
    if (c.type == mesh::CellType::Prism6) {
        std::array<Vec3, 6> p{};
        for (int i = 0; i < 6; ++i)
            p[i] = mesh.node(n[i]);
        return decompose(p, kWedge);
    }
    std::array<Vec3, 8> p{};
    for (int i = 0; i < 8; ++i)
        p[i] = mesh.node(n[i]);
    return decompose(p, kHex);
}
double percentile(const std::vector<double> &sorted, double fraction) noexcept {
    if (sorted.empty())
        return 0;
    const double pos = fraction * (sorted.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(pos), hi = std::min(lo + 1, sorted.size() - 1);
    const double t = pos - lo;
    return sorted[lo] * (1 - t) + sorted[hi] * t;
}
} // namespace

MixedMeshQualityReport analyzeMixedMeshQuality(const mesh::Mesh &value,
                                               MixedQualityMetric) noexcept {
    MixedMeshQualityReport out;
    try {
        std::vector<double> qualities;
        qualities.reserve(static_cast<std::size_t>(value.cellCount()));
        long double sum = 0, totalVolume = 0;
        double minVolume = (std::numeric_limits<double>::infinity)();
        for (Index ci = 0; ci < value.cellCount(); ++ci) {
            const auto &t = value.cell(ci).type;
            if (!volumeType(t))
                continue;
            ++out.volumeCells;
            if (t == mesh::CellType::Tet4)
                ++out.tet4;
            else if (t == mesh::CellType::Pyramid5)
                ++out.pyramid5;
            else if (t == mesh::CellType::Prism6)
                ++out.prism6;
            else if (t == mesh::CellType::Hex8)
                ++out.hex8;
            const auto q = quality(value, ci);
            if (!q.supported) {
                ++out.unsupportedCells;
                continue;
            }
            if (!q.valid) {
                ++out.invertedCells;
                continue;
            }
            if (!std::isfinite(q.quality)) {
                ++out.nonfiniteQuality;
                continue;
            }
            qualities.push_back(q.quality);
            sum += q.quality;
            totalVolume += q.volume;
            minVolume = std::min(minVolume, q.volume);
        }
        std::sort(qualities.begin(), qualities.end());
        if (!qualities.empty()) {
            out.minimumQuality = qualities.front();
            out.percentile01 = percentile(qualities, 0.001);
            out.percentile1 = percentile(qualities, 0.01);
            out.percentile5 = percentile(qualities, 0.05);
            out.medianQuality = percentile(qualities, 0.5);
            out.meanQuality = static_cast<double>(sum / qualities.size());
            out.minimumVolume = minVolume;
            out.totalVolume = static_cast<double>(totalVolume);
        }
        out.valid = out.volumeCells > 0 && out.unsupportedCells == 0 &&
                    out.invalidConnectivity == 0 && out.invertedCells == 0 &&
                    out.nonfiniteQuality == 0 && qualities.size() == out.volumeCells;
    } catch (const std::bad_alloc &) {
        out.complete = false;
        out.allocationFailed = true;
    } catch (...) {
        out.complete = false;
    }
    return out;
}

MixedMeshOptimizationOutcome
optimizeMixedMeshQuality(mesh::Mesh &value, const MixedMeshOptimizationRequest &request) noexcept {
    MixedMeshOptimizationOutcome out;
    out.before = analyzeMixedMeshQuality(value);
    if (out.before.allocationFailed) {
        out.reason = MixedMeshOptimizationReason::AllocationFailed;
        return out;
    }
    if (!out.before.complete || !out.before.valid) {
        out.reason = MixedMeshOptimizationReason::InvalidInput;
        return out;
    }
    if (request.maximumIterations == 0 || !std::isfinite(request.relaxation) ||
        request.relaxation <= 0 || request.relaxation > 1 ||
        !std::isfinite(request.minimumAcceptedImprovement) ||
        request.minimumAcceptedImprovement < 0) {
        out.reason = MixedMeshOptimizationReason::InvalidRequest;
        return out;
    }
    try {
        std::vector<std::vector<Index>> incident(static_cast<std::size_t>(value.nodeCount())),
            neighbors(static_cast<std::size_t>(value.nodeCount()));
        for (Index ci = 0; ci < value.cellCount(); ++ci) {
            const auto &t = value.cell(ci).type;
            if (!volumeType(t))
                continue;
            const int count = mesh::nodesPerCell(t);
            const Index *n = value.cellNodes(ci);
            for (int i = 0; i < count; ++i) {
                incident[static_cast<std::size_t>(n[i])].push_back(ci);
                for (int j = 0; j < count; ++j)
                    if (i != j)
                        neighbors[static_cast<std::size_t>(n[i])].push_back(n[j]);
            }
        }
        for (auto &row : neighbors) {
            std::sort(row.begin(), row.end());
            row.erase(std::unique(row.begin(), row.end()), row.end());
        }
        for (std::uint32_t iteration = 0; iteration < request.maximumIterations; ++iteration) {
            std::uint64_t acceptedThis = 0;
            for (Index node = 0; node < value.nodeCount(); ++node) {
                if (value.nodeClass(node).dim != 3 ||
                    neighbors[static_cast<std::size_t>(node)].empty())
                    continue;
                if (++out.nodeTrials > request.maximumNodeTrials) {
                    out.reason = MixedMeshOptimizationReason::BudgetExceeded;
                    out.after = analyzeMixedMeshQuality(value);
                    if (out.after.allocationFailed)
                        out.reason = MixedMeshOptimizationReason::AllocationFailed;
                    return out;
                }
                double before = (std::numeric_limits<double>::infinity)();
                for (Index ci : incident[static_cast<std::size_t>(node)])
                    before = std::min(before, quality(value, ci).quality);
                Vec3 centroid{};
                for (Index other : neighbors[static_cast<std::size_t>(node)])
                    centroid = centroid + value.node(other);
                centroid = centroid /
                           static_cast<double>(neighbors[static_cast<std::size_t>(node)].size());
                const Vec3 original = value.node(node),
                           candidate = original + (centroid - original) * request.relaxation;
                value.setNode(node, candidate);
                double after = (std::numeric_limits<double>::infinity)();
                bool valid = true;
                for (Index ci : incident[static_cast<std::size_t>(node)]) {
                    const auto q = quality(value, ci);
                    valid = valid && q.valid;
                    after = std::min(after, q.quality);
                }
                if (!valid || after < before + request.minimumAcceptedImprovement)
                    value.setNode(node, original);
                else {
                    ++acceptedThis;
                    ++out.acceptedMoves;
                }
            }
            ++out.iterations;
            if (acceptedThis == 0)
                break;
        }
        out.after = analyzeMixedMeshQuality(value);
        if (out.after.allocationFailed) {
            out.reason = MixedMeshOptimizationReason::AllocationFailed;
            return out;
        }
        out.completed = out.after.valid;
        out.reason = out.completed ? MixedMeshOptimizationReason::None
                                   : MixedMeshOptimizationReason::InvalidInput;
        return out;
    } catch (const std::bad_alloc &) {
        out.reason = MixedMeshOptimizationReason::AllocationFailed;
    } catch (...) {
        out.reason = MixedMeshOptimizationReason::InternalError;
    }
    return out;
}

} // namespace paacs::algo
