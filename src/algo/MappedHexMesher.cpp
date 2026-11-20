// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MappedHexMesher.hpp"

#include "paacs/topo/ElementTraits.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/CertifiedMeshStorage.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <set>
#include <vector>

namespace paacs::algo {
namespace {

using FaceKey = std::array<Index, 4>;

bool finite(const Vec3& value) noexcept {
    return std::isfinite(value.x) && std::isfinite(value.y) &&
           std::isfinite(value.z);
}

// A positive ideal affine frame does not certify the rounded stored nodes.
// If coordinate error can approach a cell altitude, inspect emitted cells.
bool resolvedCoordinates(const StructuredBlockGeometry &g, double determinant,
                         double du, double dv, double dw) noexcept {
    const auto magnitude=[](Vec3 p){return std::max({std::abs(p.x),std::abs(p.y),std::abs(p.z)});};
    const double scale=magnitude(g.origin)+magnitude(g.axisU)+magnitude(g.axisV)+magnitude(g.axisW);
    const double error=128*std::numeric_limits<double>::epsilon()*scale;
    const auto length=[](Vec3 p){return std::hypot(p.x,p.y,p.z);};
    return std::isfinite(error) &&
        du*(determinant/length(g.axisV.cross(g.axisW)))>error &&
        dv*(determinant/length(g.axisW.cross(g.axisU)))>error &&
        dw*(determinant/length(g.axisU.cross(g.axisV)))>error;
}

bool validClass(const mesh::Classification& value, int dimension) noexcept {
    return value.dim == dimension && value.tag >= 0;
}

bool checkedAdd(std::uint64_t a, std::uint64_t b,
                std::uint64_t& out) noexcept {
    if (b > (std::numeric_limits<std::uint64_t>::max)() - a) return false;
    out = a + b;
    return true;
}

bool checkedMul(std::uint64_t a, std::uint64_t b,
                std::uint64_t& out) noexcept {
    if (a != 0 && b > (std::numeric_limits<std::uint64_t>::max)() / a)
        return false;
    out = a * b;
    return true;
}

struct BlockCounts {
    std::uint64_t nodes = 0;
    std::uint64_t cells = 0;
    std::uint64_t boundaryFaces = 0;
    std::uint64_t connectivity = 0;
};

bool blockCounts(const StructuredBlockGeometry& g, bool boundary,
                 BlockCounts& out) noexcept {
    const std::uint64_t u = g.cellsU;
    const std::uint64_t v = g.cellsV;
    const std::uint64_t w = g.cellsW;
    std::uint64_t uv = 0, uw = 0, vw = 0;
    std::uint64_t uvw = 0, nodeUv = 0;
    if (!checkedMul(u, v, uv) || !checkedMul(u, w, uw) ||
        !checkedMul(v, w, vw) || !checkedMul(uv, w, uvw) ||
        !checkedMul(u + 1, v + 1, nodeUv) ||
        !checkedMul(nodeUv, w + 1, out.nodes))
        return false;
    out.cells = uvw;
    if (boundary) {
        std::uint64_t sum = 0;
        if (!checkedAdd(uv, uw, sum) || !checkedAdd(sum, vw, sum) ||
            !checkedMul(sum, 2, out.boundaryFaces))
            return false;
    }
    std::uint64_t hexConn = 0, quadConn = 0;
    if (!checkedMul(out.cells, 8, hexConn) ||
        !checkedMul(out.boundaryFaces, 4, quadConn) ||
        !checkedAdd(hexConn, quadConn, out.connectivity))
        return false;
    return true;
}

Index nodeIndex(std::uint64_t i, std::uint64_t j, std::uint64_t k,
                std::uint64_t u, std::uint64_t v) noexcept {
    return static_cast<Index>((k * (v + 1) + j) * (u + 1) + i);
}

FaceKey keyOf(const std::array<Index, 4>& values) noexcept {
    FaceKey key = values;
    std::sort(key.begin(), key.end());
    return key;
}

struct JacobianSample {
    double determinant = 0.0;
    double scaled = 0.0;
};

JacobianSample hexJacobian(const std::array<Vec3, 8>& p,
                           double xi, double eta, double zeta) noexcept {
    constexpr double sx[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
    constexpr double sy[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
    constexpr double sz[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
    Vec3 dx{}, dy{}, dz{};
    for (std::size_t i = 0; i < p.size(); ++i) {
        const double dxi = 0.125 * sx[i] * (1.0 + sy[i] * eta) *
                           (1.0 + sz[i] * zeta);
        const double deta = 0.125 * sy[i] * (1.0 + sx[i] * xi) *
                            (1.0 + sz[i] * zeta);
        const double dzeta = 0.125 * sz[i] * (1.0 + sx[i] * xi) *
                             (1.0 + sy[i] * eta);
        dx = dx + p[i] * dxi;
        dy = dy + p[i] * deta;
        dz = dz + p[i] * dzeta;
    }
    JacobianSample out;
    out.determinant = dx.dot(dy.cross(dz));
    const double denom = dx.norm() * dy.norm() * dz.norm();
    out.scaled = denom > 0.0 ? out.determinant / denom : 0.0;
    return out;
}

Hex8AuditReport certifiedAudit(const StructuredBlockGeometry& g,
                               const BlockCounts& counts,
                               double determinant,
                               bool boundaryQuads,
                               const StructuredAxisDistribution& uGrowth,
                               const StructuredAxisDistribution& vGrowth,
                               const StructuredAxisDistribution& wGrowth) noexcept {
    Hex8AuditReport out;
    out.nodes = counts.nodes;
    out.hexes = counts.cells;
    out.boundaryQuads = boundaryQuads ? counts.boundaryFaces : 0;
    out.boundaryFaces = counts.boundaryFaces;
    const std::uint64_t faceIncidences = counts.cells * 6U;
    out.interiorFaces =
        (faceIncidences - counts.boundaryFaces) / 2U;
    const double minimumCellScale = uGrowth.minimumCellFraction *
        vGrowth.minimumCellFraction * wGrowth.minimumCellFraction;
    const double maximumCellScale = uGrowth.maximumCellFraction *
        vGrowth.maximumCellFraction * wGrowth.maximumCellFraction;
    const double minimumCellJacobian =
        determinant * minimumCellScale / 8.0;
    const double maximumCellJacobian =
        determinant * maximumCellScale / 8.0;
    out.minimumJacobian = minimumCellJacobian;
    out.maximumJacobian = maximumCellJacobian;
    const double denominator =
        g.axisU.norm() * g.axisV.norm() * g.axisW.norm();
    out.minimumScaledJacobian =
        denominator > 0.0 ? determinant / denominator : 0.0;
    out.totalVolume = determinant;
    out.valid = std::isfinite(minimumCellJacobian) &&
        minimumCellJacobian > 0.0 &&
        std::isfinite(maximumCellJacobian) &&
        maximumCellJacobian > 0.0 &&
        std::isfinite(out.minimumScaledJacobian) &&
        out.minimumScaledJacobian > 0.0;
    return out;
}

MappedHexOutcome rejected(MappedHexReason reason) noexcept {
    MappedHexOutcome out;
    out.status = MappedHexStatus::Rejected;
    out.reason = reason;
    out.phase = MappedHexPhase::Preflight;
    return out;
}

}  // namespace

Hex8AuditReport auditHex8Mesh(const mesh::Mesh& value,
                              bool requireBoundaryQuads) noexcept {
    Hex8AuditReport out;
    out.nodes = static_cast<std::uint64_t>(value.nodeCount());
    try {
        for (Index i = 0; i < value.nodeCount(); ++i) {
            if (!finite(value.node(i))) ++out.nonfiniteNodes;
            const auto cls = value.nodeClass(i);
            if (cls.dim < 0 || cls.dim > 3 || cls.tag < 0)
                ++out.invalidClassifications;
        }

        std::map<FaceKey, std::uint8_t> faceUses;
        std::map<FaceKey, std::uint32_t> quadUses;
        double minJac = (std::numeric_limits<double>::infinity)();
        double maxJac = -(std::numeric_limits<double>::infinity)();
        double minScaled = (std::numeric_limits<double>::infinity)();
        long double totalVolume = 0.0L;
        constexpr double corner[2] = {-1.0, 1.0};
        constexpr double gauss = 0.577350269189625764509148780501957;
        constexpr double gp[2] = {-gauss, gauss};
        const auto& hexTraits = topo::traits(topo::ElementType::Hex8);

        for (Index ci = 0; ci < value.cellCount(); ++ci) {
            const auto& cell = value.cell(ci);
            if (cell.type == mesh::CellType::Quad4) {
                ++out.boundaryQuads;
                if (!validClass(cell.onEntity, 2))
                    ++out.invalidClassifications;
                std::array<Index, 4> nodes{};
                const Index* src = value.cellNodes(ci);
                bool bad = false;
                for (std::size_t n = 0; n < nodes.size(); ++n) {
                    nodes[n] = src[n];
                    if (nodes[n] < 0 || nodes[n] >= value.nodeCount())
                        bad = true;
                }
                auto unique = nodes;
                std::sort(unique.begin(), unique.end());
                if (std::adjacent_find(unique.begin(), unique.end()) !=
                    unique.end())
                    bad = true;
                if (bad) {
                    ++out.invalidConnectivity;
                } else {
                    ++quadUses[keyOf(nodes)];
                }
                continue;
            }
            if (cell.type != mesh::CellType::Hex8) {
                ++out.unsupportedCells;
                continue;
            }
            ++out.hexes;
            if (!validClass(cell.onEntity, 3))
                ++out.invalidClassifications;
            std::array<Index, 8> nodes{};
            std::array<Vec3, 8> points{};
            const Index* src = value.cellNodes(ci);
            bool bad = false;
            for (std::size_t n = 0; n < nodes.size(); ++n) {
                nodes[n] = src[n];
                if (nodes[n] < 0 || nodes[n] >= value.nodeCount()) {
                    bad = true;
                } else {
                    points[n] = value.node(nodes[n]);
                }
            }
            auto unique = nodes;
            std::sort(unique.begin(), unique.end());
            if (std::adjacent_find(unique.begin(), unique.end()) !=
                unique.end()) {
                ++out.duplicateCellNodes;
                bad = true;
            }
            if (bad) {
                ++out.invalidConnectivity;
                continue;
            }

            for (double xi : corner)
                for (double eta : corner)
                    for (double zeta : corner) {
                        const auto sample = hexJacobian(
                            points, xi, eta, zeta);
                        if (!std::isfinite(sample.determinant) ||
                            !std::isfinite(sample.scaled) ||
                            sample.determinant <= 0.0)
                            ++out.nonPositiveJacobians;
                        minJac = std::min(minJac, sample.determinant);
                        maxJac = std::max(maxJac, sample.determinant);
                        minScaled = std::min(minScaled, sample.scaled);
                    }
            long double volume = 0.0L;
            for (double xi : gp)
                for (double eta : gp)
                    for (double zeta : gp) {
                        const auto sample = hexJacobian(
                            points, xi, eta, zeta);
                        if (!std::isfinite(sample.determinant) ||
                            sample.determinant <= 0.0)
                            ++out.nonPositiveJacobians;
                        volume += static_cast<long double>(
                            sample.determinant);
                    }
            if (!std::isfinite(static_cast<double>(volume)) ||
                volume <= 0.0L)
                ++out.nonPositiveJacobians;
            else
                totalVolume += volume;

            for (std::size_t fi = 0; fi < hexTraits.faceCount; ++fi) {
                const auto& local = hexTraits.faces[fi];
                std::array<Index, 4> face{};
                for (std::size_t n = 0; n < 4; ++n)
                    face[n] = nodes[local.v[n]];
                auto& use = faceUses[keyOf(face)];
                if (use != (std::numeric_limits<std::uint8_t>::max)())
                    ++use;
            }
        }

        for (const auto& [key, count] : faceUses) {
            if (count == 1) {
                ++out.boundaryFaces;
                const auto it = quadUses.find(key);
                const std::uint32_t quads =
                    it == quadUses.end() ? 0U : it->second;
                if (requireBoundaryQuads && quads == 0)
                    ++out.missingBoundaryQuads;
                if (quads > 1)
                    out.duplicateBoundaryQuads += quads - 1U;
            } else if (count == 2) {
                ++out.interiorFaces;
                const auto it = quadUses.find(key);
                if (it != quadUses.end())
                    out.unexpectedBoundaryQuads += it->second;
            } else {
                ++out.nonmanifoldFaces;
            }
        }
        for (const auto& [key, count] : quadUses)
            if (faceUses.find(key) == faceUses.end())
                out.unexpectedBoundaryQuads += count;

        if (out.hexes == 0) {
            minJac = 0.0;
            maxJac = 0.0;
            minScaled = 0.0;
        }
        out.minimumJacobian = minJac;
        out.maximumJacobian = maxJac;
        out.minimumScaledJacobian = minScaled;
        out.totalVolume = static_cast<double>(totalVolume);
        out.valid = out.hexes != 0 && out.unsupportedCells == 0 &&
            out.nonfiniteNodes == 0 && out.invalidConnectivity == 0 &&
            out.duplicateCellNodes == 0 &&
            out.invalidClassifications == 0 &&
            out.nonPositiveJacobians == 0 && out.nonmanifoldFaces == 0 &&
            out.missingBoundaryQuads == 0 &&
            out.duplicateBoundaryQuads == 0 &&
            out.unexpectedBoundaryQuads == 0;
    } catch (const std::bad_alloc&) {
        out.complete = false;
        out.valid = false;
        out.allocationFailed = true;
    } catch (...) {
        out.complete = false;
        out.valid = false;
    }
    return out;
}

MappedHexOutcome generateMappedHexBlock(
        const MappedHexBlockRequest& request,
        mesh::Mesh& destination) noexcept {
    MappedHexOutcome out;
    const auto& g = request.geometry;
    if (g.cellsU == 0 || g.cellsV == 0 || g.cellsW == 0)
        return rejected(MappedHexReason::InvalidCellCounts);
    if (!validClass(g.volumeClassification, 3))
        return rejected(MappedHexReason::InvalidClassification);
    for (const auto cls : g.vertexClassifications)
        if (!validClass(cls, 0))
            return rejected(MappedHexReason::InvalidClassification);
    for (const auto cls : g.edgeClassifications)
        if (!validClass(cls, 1))
            return rejected(MappedHexReason::InvalidClassification);
    if (request.emitBoundaryQuads)
        for (const auto cls : g.boundaryClassifications)
            if (!validClass(cls, 2))
                return rejected(MappedHexReason::InvalidClassification);
    if (!finite(g.origin) || !finite(g.axisU) || !finite(g.axisV) ||
        !finite(g.axisW))
        return rejected(MappedHexReason::NonfiniteFrame);
    const double determinant = g.axisU.dot(g.axisV.cross(g.axisW));
    if (!std::isfinite(determinant) || determinant <= 0.0)
        return rejected(MappedHexReason::NonPositiveFrame);

    BlockCounts counts;
    if (!blockCounts(g, request.emitBoundaryQuads, counts))
        return rejected(MappedHexReason::CountOverflow);
    if (counts.nodes > g.maxNodes || counts.cells > g.maxCells ||
        counts.boundaryFaces > g.maxFaces)
        return rejected(MappedHexReason::BudgetExceeded);
    const auto indexMax = static_cast<std::uint64_t>(
        (std::numeric_limits<Index>::max)());
    if (counts.nodes > indexMax || counts.cells > indexMax ||
        counts.boundaryFaces > indexMax ||
        counts.nodes > (std::numeric_limits<std::size_t>::max)() ||
        counts.cells + counts.boundaryFaces >
            (std::numeric_limits<std::size_t>::max)() ||
        counts.connectivity >
            (std::numeric_limits<std::size_t>::max)())
        return rejected(MappedHexReason::CountOverflow);

    const auto uGrowth = buildStructuredAxisDistribution(
        g.cellsU, request.growth.u);
    const auto vGrowth = buildStructuredAxisDistribution(
        g.cellsV, request.growth.v);
    const auto wGrowth = buildStructuredAxisDistribution(
        g.cellsW, request.growth.w);
    if (!uGrowth.valid || !vGrowth.valid || !wGrowth.valid) {
        const bool allocationFailed =
            uGrowth.reason == StructuredGrowthReason::AllocationFailed ||
            vGrowth.reason == StructuredGrowthReason::AllocationFailed ||
            wGrowth.reason == StructuredGrowthReason::AllocationFailed;
        if (allocationFailed) {
            MappedHexOutcome failed;
            failed.status = MappedHexStatus::Failed;
            failed.reason = MappedHexReason::AllocationFailed;
            failed.phase = MappedHexPhase::Preflight;
            return failed;
        }
        return rejected(MappedHexReason::InvalidGrowth);
    }

    out.phase = MappedHexPhase::Construction;
    out.metrics.nodes = counts.nodes;
    out.metrics.hexes = counts.cells;
    out.metrics.boundaryQuads = counts.boundaryFaces;
    out.metrics.connectivityEntries = counts.connectivity;
    out.metrics.blockDeterminant = determinant;
    out.metrics.growthPreparationWorkUnits =
        static_cast<std::uint64_t>(uGrowth.boundaries.size() +
                                   vGrowth.boundaries.size() +
                                   wGrowth.boundaries.size());
    out.metrics.minimumCellFraction =
        uGrowth.minimumCellFraction * vGrowth.minimumCellFraction *
        wGrowth.minimumCellFraction;
    out.metrics.maximumCellFraction =
        uGrowth.maximumCellFraction * vGrowth.maximumCellFraction *
        wGrowth.maximumCellFraction;
    out.metrics.maximumAdjacentGrowth = std::max({
        uGrowth.maximumAdjacentRatio,
        vGrowth.maximumAdjacentRatio,
        wGrowth.maximumAdjacentRatio});
    out.metrics.logicalBytes = counts.nodes *
            (sizeof(Vec3) + sizeof(mesh::Classification)) +
        (counts.cells + counts.boundaryFaces) * sizeof(mesh::Cell) +
        counts.connectivity * sizeof(Index);

    try {
        mesh::detail::CertifiedMeshStorage storage(
            static_cast<std::size_t>(counts.nodes),
            static_cast<std::size_t>(counts.cells + counts.boundaryFaces),
            static_cast<std::size_t>(counts.connectivity));

        const std::uint64_t nodeU = g.cellsU + 1U;
        const std::uint64_t nodeV = g.cellsV + 1U;
        const std::uint64_t nodePlane = nodeU * nodeV;
        core::parallelForWithThreads(
            0, static_cast<std::size_t>(counts.nodes),
            request.executionThreads, [&](std::size_t ordinal) {
                const std::uint64_t value = ordinal;
                const std::uint64_t k = value / nodePlane;
                const std::uint64_t remainder = value % nodePlane;
                const std::uint64_t j = remainder / nodeU;
                const std::uint64_t i = remainder % nodeU;
                const double ui = uGrowth.coordinate(
                    static_cast<std::uint32_t>(i));
                const double vj = vGrowth.coordinate(
                    static_cast<std::uint32_t>(j));
                const double wk = wGrowth.coordinate(
                    static_cast<std::uint32_t>(k));
                storage.setNode(
                    ordinal,
                    g.origin + g.axisU * ui + g.axisV * vj + g.axisW * wk,
                    structuredBlockNodeClassification(g, i, j, k));
            });

        const std::uint64_t cellPlane = g.cellsU * g.cellsV;
        core::parallelForWithThreads(
            0, static_cast<std::size_t>(counts.cells),
            request.executionThreads, [&](std::size_t ordinal) {
                const std::uint64_t value = ordinal;
                const std::uint64_t k = value / cellPlane;
                const std::uint64_t remainder = value % cellPlane;
                const std::uint64_t j = remainder / g.cellsU;
                const std::uint64_t i = remainder % g.cellsU;
                const std::array<Index, 8> conn{
                    nodeIndex(i, j, k, g.cellsU, g.cellsV),
                    nodeIndex(i + 1, j, k, g.cellsU, g.cellsV),
                    nodeIndex(i + 1, j + 1, k, g.cellsU, g.cellsV),
                    nodeIndex(i, j + 1, k, g.cellsU, g.cellsV),
                    nodeIndex(i, j, k + 1, g.cellsU, g.cellsV),
                    nodeIndex(i + 1, j, k + 1, g.cellsU, g.cellsV),
                    nodeIndex(i + 1, j + 1, k + 1,
                              g.cellsU, g.cellsV),
                    nodeIndex(i, j + 1, k + 1, g.cellsU, g.cellsV)};
                storage.setCell(ordinal, ordinal * 8U,
                                mesh::CellType::Hex8, conn,
                                g.volumeClassification);
            });

        if (request.emitBoundaryQuads) {
            const std::uint64_t volumeConnectivity = counts.cells * 8U;
            const std::uint64_t uv = g.cellsU * g.cellsV;
            const std::uint64_t uw = g.cellsU * g.cellsW;
            const std::uint64_t vw = g.cellsV * g.cellsW;
            const auto setQuad = [&](std::uint64_t pairBase,
                                     std::uint64_t pairOrdinal,
                                     bool second,
                                     const std::array<Index,4>& nodes,
                                     StructuredBlockFace face) {
                const std::uint64_t local = pairBase + 2U * pairOrdinal +
                                            (second ? 1U : 0U);
                storage.setCell(
                    static_cast<std::size_t>(counts.cells + local),
                    static_cast<std::size_t>(volumeConnectivity + local * 4U),
                    mesh::CellType::Quad4, nodes,
                    g.boundaryClassifications[
                        static_cast<std::size_t>(face)]);
            };

            core::parallelForWithThreads(
                0, static_cast<std::size_t>(uv), request.executionThreads,
                [&](std::size_t ordinal) {
                    const std::uint64_t j = ordinal / g.cellsU;
                    const std::uint64_t i = ordinal % g.cellsU;
                    setQuad(0, ordinal, false,
                        {nodeIndex(i, j, 0, g.cellsU, g.cellsV),
                         nodeIndex(i, j + 1, 0, g.cellsU, g.cellsV),
                         nodeIndex(i + 1, j + 1, 0, g.cellsU, g.cellsV),
                         nodeIndex(i + 1, j, 0, g.cellsU, g.cellsV)},
                        StructuredBlockFace::WMin);
                    setQuad(0, ordinal, true,
                        {nodeIndex(i, j, g.cellsW, g.cellsU, g.cellsV),
                         nodeIndex(i + 1, j, g.cellsW, g.cellsU, g.cellsV),
                         nodeIndex(i + 1, j + 1, g.cellsW,
                                   g.cellsU, g.cellsV),
                         nodeIndex(i, j + 1, g.cellsW, g.cellsU, g.cellsV)},
                        StructuredBlockFace::WMax);
                });

            core::parallelForWithThreads(
                0, static_cast<std::size_t>(uw), request.executionThreads,
                [&](std::size_t ordinal) {
                    const std::uint64_t k = ordinal / g.cellsU;
                    const std::uint64_t i = ordinal % g.cellsU;
                    setQuad(2U * uv, ordinal, false,
                        {nodeIndex(i, 0, k, g.cellsU, g.cellsV),
                         nodeIndex(i + 1, 0, k, g.cellsU, g.cellsV),
                         nodeIndex(i + 1, 0, k + 1, g.cellsU, g.cellsV),
                         nodeIndex(i, 0, k + 1, g.cellsU, g.cellsV)},
                        StructuredBlockFace::VMin);
                    setQuad(2U * uv, ordinal, true,
                        {nodeIndex(i + 1, g.cellsV, k, g.cellsU, g.cellsV),
                         nodeIndex(i, g.cellsV, k, g.cellsU, g.cellsV),
                         nodeIndex(i, g.cellsV, k + 1, g.cellsU, g.cellsV),
                         nodeIndex(i + 1, g.cellsV, k + 1,
                                   g.cellsU, g.cellsV)},
                        StructuredBlockFace::VMax);
                });

            core::parallelForWithThreads(
                0, static_cast<std::size_t>(vw), request.executionThreads,
                [&](std::size_t ordinal) {
                    const std::uint64_t k = ordinal / g.cellsV;
                    const std::uint64_t j = ordinal % g.cellsV;
                    setQuad(2U * uv + 2U * uw, ordinal, false,
                        {nodeIndex(0, j + 1, k, g.cellsU, g.cellsV),
                         nodeIndex(0, j, k, g.cellsU, g.cellsV),
                         nodeIndex(0, j, k + 1, g.cellsU, g.cellsV),
                         nodeIndex(0, j + 1, k + 1, g.cellsU, g.cellsV)},
                        StructuredBlockFace::UMin);
                    setQuad(2U * uv + 2U * uw, ordinal, true,
                        {nodeIndex(g.cellsU, j, k, g.cellsU, g.cellsV),
                         nodeIndex(g.cellsU, j + 1, k, g.cellsU, g.cellsV),
                         nodeIndex(g.cellsU, j + 1, k + 1,
                                   g.cellsU, g.cellsV),
                         nodeIndex(g.cellsU, j, k + 1,
                                   g.cellsU, g.cellsV)},
                        StructuredBlockFace::UMax);
                });
        }

        mesh::Mesh local;
        storage.publishTo(local);

        out.metrics.constructionWorkUnits = counts.nodes + counts.cells +
            counts.boundaryFaces + counts.connectivity +
            out.metrics.growthPreparationWorkUnits;
        out.phase = MappedHexPhase::Audit;
        if (static_cast<std::uint64_t>(local.nodeCount()) != counts.nodes ||
            static_cast<std::uint64_t>(local.cellCount()) !=
                counts.cells + counts.boundaryFaces ||
            static_cast<std::uint64_t>(local.connectivityCount()) !=
                counts.connectivity) {
            out.status = MappedHexStatus::Failed;
            out.reason = MappedHexReason::AuditFailed;
            return out;
        }
        out.metrics.usedFullGeometricAudit = request.validation ==
                MappedHexValidationMode::FullTopologyAudit ||
            !resolvedCoordinates(g, determinant, uGrowth.minimumCellFraction,
                                 vGrowth.minimumCellFraction, wGrowth.minimumCellFraction);
        out.metrics.audit = out.metrics.usedFullGeometricAudit
            ? auditHex8Mesh(local, request.emitBoundaryQuads)
            : certifiedAudit(g, counts, determinant,
                             request.emitBoundaryQuads,
                             uGrowth, vGrowth, wGrowth);
        if (!out.metrics.audit.complete) {
            out.status = MappedHexStatus::Failed;
            out.reason = out.metrics.audit.allocationFailed
                ? MappedHexReason::AllocationFailed
                : MappedHexReason::AuditFailed;
            return out;
        }
        if (!out.metrics.audit.valid ||
            std::abs(out.metrics.audit.totalVolume - determinant) >
                128.0 * (std::numeric_limits<double>::epsilon)() *
                std::max(1.0, std::abs(determinant))) {
            out.status = MappedHexStatus::Failed;
            out.reason = MappedHexReason::AuditFailed;
            return out;
        }
        out.phase = MappedHexPhase::Publication;
        destination.swap(local);
        out.status = MappedHexStatus::Published;
        out.reason = MappedHexReason::None;
        out.destinationUnchanged = false;
        return out;
    } catch (const std::bad_alloc&) {
        out.status = MappedHexStatus::Failed;
        out.reason = MappedHexReason::AllocationFailed;
    } catch (...) {
        out.status = MappedHexStatus::Failed;
        out.reason = MappedHexReason::InternalError;
    }
    return out;
}

}  // namespace paacs::algo
