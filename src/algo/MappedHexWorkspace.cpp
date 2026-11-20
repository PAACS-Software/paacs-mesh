// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MappedHexMesher.hpp"

#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/ReusableCertifiedMeshStorage.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <vector>

namespace paacs::algo {
namespace {

bool finiteWorkspace(const Vec3& value) noexcept {
    return std::isfinite(value.x) && std::isfinite(value.y) &&
           std::isfinite(value.z);
}

bool validClassWorkspace(const mesh::Classification& value,
                         int dimension) noexcept {
    return value.dim == dimension && value.tag >= 0;
}

bool checkedAddWorkspace(std::uint64_t a, std::uint64_t b,
                         std::uint64_t& out) noexcept {
    if (b > (std::numeric_limits<std::uint64_t>::max)() - a) return false;
    out = a + b;
    return true;
}

bool checkedMulWorkspace(std::uint64_t a, std::uint64_t b,
                         std::uint64_t& out) noexcept {
    if (a != 0 && b > (std::numeric_limits<std::uint64_t>::max)() / a)
        return false;
    out = a * b;
    return true;
}

struct WorkspaceBlockCounts {
    std::uint64_t nodes = 0;
    std::uint64_t cells = 0;
    std::uint64_t boundaryFaces = 0;
    std::uint64_t connectivity = 0;
};

bool workspaceBlockCounts(const StructuredBlockGeometry& g, bool boundary,
                          WorkspaceBlockCounts& out) noexcept {
    const std::uint64_t u = g.cellsU;
    const std::uint64_t v = g.cellsV;
    const std::uint64_t w = g.cellsW;
    std::uint64_t uv = 0, uw = 0, vw = 0;
    std::uint64_t uvw = 0, nodeUv = 0;
    if (!checkedMulWorkspace(u, v, uv) ||
        !checkedMulWorkspace(u, w, uw) ||
        !checkedMulWorkspace(v, w, vw) ||
        !checkedMulWorkspace(uv, w, uvw) ||
        !checkedMulWorkspace(u + 1, v + 1, nodeUv) ||
        !checkedMulWorkspace(nodeUv, w + 1, out.nodes))
        return false;
    out.cells = uvw;
    if (boundary) {
        std::uint64_t sum = 0;
        if (!checkedAddWorkspace(uv, uw, sum) ||
            !checkedAddWorkspace(sum, vw, sum) ||
            !checkedMulWorkspace(sum, 2, out.boundaryFaces))
            return false;
    }
    std::uint64_t hexConn = 0, quadConn = 0;
    if (!checkedMulWorkspace(out.cells, 8, hexConn) ||
        !checkedMulWorkspace(out.boundaryFaces, 4, quadConn) ||
        !checkedAddWorkspace(hexConn, quadConn, out.connectivity))
        return false;
    return true;
}

Index workspaceNodeIndex(std::uint64_t i, std::uint64_t j, std::uint64_t k,
                         std::uint64_t u, std::uint64_t v) noexcept {
    return static_cast<Index>((k * (v + 1) + j) * (u + 1) + i);
}

Hex8AuditReport workspaceCertifiedAudit(
        const StructuredBlockGeometry& g,
        const WorkspaceBlockCounts& counts,
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
    out.interiorFaces = (faceIncidences - counts.boundaryFaces) / 2U;
    const double minimumCellScale = uGrowth.minimumCellFraction *
        vGrowth.minimumCellFraction * wGrowth.minimumCellFraction;
    const double maximumCellScale = uGrowth.maximumCellFraction *
        vGrowth.maximumCellFraction * wGrowth.maximumCellFraction;
    out.minimumJacobian = determinant * minimumCellScale / 8.0;
    out.maximumJacobian = determinant * maximumCellScale / 8.0;
    const double denominator =
        g.axisU.norm() * g.axisV.norm() * g.axisW.norm();
    out.minimumScaledJacobian =
        denominator > 0.0 ? determinant / denominator : 0.0;
    out.totalVolume = determinant;
    out.valid = std::isfinite(out.minimumJacobian) &&
        out.minimumJacobian > 0.0 &&
        std::isfinite(out.maximumJacobian) &&
        out.maximumJacobian > 0.0 &&
        std::isfinite(out.minimumScaledJacobian) &&
        out.minimumScaledJacobian > 0.0;
    return out;
}

MappedHexOutcome workspaceRejected(MappedHexReason reason) noexcept {
    MappedHexOutcome out;
    out.status = MappedHexStatus::Rejected;
    out.reason = reason;
    out.phase = MappedHexPhase::Preflight;
    return out;
}

}  // namespace

struct MappedHexWorkspace::Impl {
    mesh::detail::ReusableCertifiedMeshStorage storage;
};

MappedHexWorkspace::MappedHexWorkspace()
    : impl_(std::make_unique<Impl>()) {}
MappedHexWorkspace::~MappedHexWorkspace() = default;
MappedHexWorkspace::MappedHexWorkspace(MappedHexWorkspace&&) noexcept = default;
MappedHexWorkspace& MappedHexWorkspace::operator=(
    MappedHexWorkspace&&) noexcept = default;

void MappedHexWorkspace::release() noexcept {
    if (impl_) impl_->storage.release();
}

std::uint64_t MappedHexWorkspace::retainedLogicalBytes() const noexcept {
    return impl_ ? impl_->storage.logicalBytes() : 0U;
}

MappedHexOutcome generateMappedHexBlock(
        const MappedHexBlockRequest& request,
        MappedHexWorkspace& workspace,
        mesh::Mesh& destination) noexcept {
    MappedHexOutcome out;
    if (!workspace.impl_) {
        out.status = MappedHexStatus::Failed;
        out.reason = MappedHexReason::InternalError;
        out.phase = MappedHexPhase::Preflight;
        return out;
    }
    const auto& g = request.geometry;
    if (g.cellsU == 0 || g.cellsV == 0 || g.cellsW == 0)
        return workspaceRejected(MappedHexReason::InvalidCellCounts);
    if (!validClassWorkspace(g.volumeClassification, 3))
        return workspaceRejected(MappedHexReason::InvalidClassification);
    for (const auto cls : g.vertexClassifications)
        if (!validClassWorkspace(cls, 0))
            return workspaceRejected(MappedHexReason::InvalidClassification);
    for (const auto cls : g.edgeClassifications)
        if (!validClassWorkspace(cls, 1))
            return workspaceRejected(MappedHexReason::InvalidClassification);
    if (request.emitBoundaryQuads)
        for (const auto cls : g.boundaryClassifications)
            if (!validClassWorkspace(cls, 2))
                return workspaceRejected(MappedHexReason::InvalidClassification);
    if (!finiteWorkspace(g.origin) || !finiteWorkspace(g.axisU) || !finiteWorkspace(g.axisV) ||
        !finiteWorkspace(g.axisW))
        return workspaceRejected(MappedHexReason::NonfiniteFrame);
    const double determinant = g.axisU.dot(g.axisV.cross(g.axisW));
    if (!std::isfinite(determinant) || determinant <= 0.0)
        return workspaceRejected(MappedHexReason::NonPositiveFrame);

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
        return workspaceRejected(MappedHexReason::InvalidGrowth);
    }

    WorkspaceBlockCounts counts;
    if (!workspaceBlockCounts(g, request.emitBoundaryQuads, counts))
        return workspaceRejected(MappedHexReason::CountOverflow);
    if (counts.nodes > g.maxNodes || counts.cells > g.maxCells ||
        counts.boundaryFaces > g.maxFaces)
        return workspaceRejected(MappedHexReason::BudgetExceeded);
    const auto indexMax = static_cast<std::uint64_t>(
        (std::numeric_limits<Index>::max)());
    if (counts.nodes > indexMax || counts.cells > indexMax ||
        counts.boundaryFaces > indexMax ||
        counts.nodes > (std::numeric_limits<std::size_t>::max)() ||
        counts.cells + counts.boundaryFaces >
            (std::numeric_limits<std::size_t>::max)() ||
        counts.connectivity >
            (std::numeric_limits<std::size_t>::max)())
        return workspaceRejected(MappedHexReason::CountOverflow);

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
        auto& storage = workspace.impl_->storage;
        out.metrics.reusedConstructionStorage = storage.prepare(
            static_cast<std::size_t>(counts.nodes),
            static_cast<std::size_t>(counts.cells + counts.boundaryFaces),
            static_cast<std::size_t>(counts.connectivity));

        const std::uint64_t nodeU = g.cellsU + 1U;
        const std::uint64_t nodeV = g.cellsV + 1U;
        const std::uint64_t nodePlane = nodeU * nodeV;

        if (!out.metrics.reusedConstructionStorage) {
            // Preserve the V8 one-shot path. On a cold destination the large
            // vector resize/value-initialization dominates, and the original
            // fine-grained loop remains slightly faster on the benchmark CPU.
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
                        workspaceNodeIndex(i, j, k, g.cellsU, g.cellsV),
                        workspaceNodeIndex(i + 1, j, k, g.cellsU, g.cellsV),
                        workspaceNodeIndex(i + 1, j + 1, k, g.cellsU, g.cellsV),
                        workspaceNodeIndex(i, j + 1, k, g.cellsU, g.cellsV),
                        workspaceNodeIndex(i, j, k + 1, g.cellsU, g.cellsV),
                        workspaceNodeIndex(i + 1, j, k + 1, g.cellsU, g.cellsV),
                        workspaceNodeIndex(i + 1, j + 1, k + 1,
                                  g.cellsU, g.cellsV),
                        workspaceNodeIndex(i, j + 1, k + 1, g.cellsU, g.cellsV)};
                    storage.setCell(ordinal, ordinal * 8U,
                                    mesh::CellType::Hex8, conn,
                                    g.volumeClassification);
                });
        } else {
            const std::uint64_t nodeW = g.cellsW + 1U;

            // Repeated-generation path: precompute axis contributions once,
            // then emit large contiguous ranges with no per-node/cell
            // std::function dispatch and only one div/mod decode per range.
            std::vector<Vec3> uPosition(static_cast<std::size_t>(nodeU));
            std::vector<Vec3> vPosition(static_cast<std::size_t>(nodeV));
            std::vector<Vec3> wPosition(static_cast<std::size_t>(nodeW));
            for (std::uint64_t i = 0; i < nodeU; ++i)
                uPosition[static_cast<std::size_t>(i)] =
                    g.axisU * uGrowth.coordinate(static_cast<std::uint32_t>(i));
            for (std::uint64_t j = 0; j < nodeV; ++j)
                vPosition[static_cast<std::size_t>(j)] =
                    g.axisV * vGrowth.coordinate(static_cast<std::uint32_t>(j));
            for (std::uint64_t k = 0; k < nodeW; ++k)
                wPosition[static_cast<std::size_t>(k)] =
                    g.axisW * wGrowth.coordinate(static_cast<std::uint32_t>(k));

            constexpr std::size_t kNodeChunk = 32U * 1024U;
            core::parallelChunksWithThreads(
                0, static_cast<std::size_t>(counts.nodes), kNodeChunk,
                request.executionThreads,
                [&](std::size_t, std::size_t lo, std::size_t hi) {
                    std::uint64_t ordinal = lo;
                    std::uint64_t k = ordinal / nodePlane;
                    std::uint64_t remainder = ordinal - k * nodePlane;
                    std::uint64_t j = remainder / nodeU;
                    std::uint64_t i = remainder - j * nodeU;
                    for (; ordinal < hi; ++ordinal) {
                        const Vec3& up = uPosition[static_cast<std::size_t>(i)];
                        const Vec3& vp = vPosition[static_cast<std::size_t>(j)];
                        const Vec3& wp = wPosition[static_cast<std::size_t>(k)];
                        storage.setNode(
                            static_cast<std::size_t>(ordinal),
                            {g.origin.x + up.x + vp.x + wp.x,
                             g.origin.y + up.y + vp.y + wp.y,
                             g.origin.z + up.z + vp.z + wp.z},
                            structuredBlockNodeClassification(g, i, j, k));
                        if (++i == nodeU) {
                            i = 0;
                            if (++j == nodeV) {
                                j = 0;
                                ++k;
                            }
                        }
                    }
                });

            const std::uint64_t cellPlane =
                static_cast<std::uint64_t>(g.cellsU) * g.cellsV;
            constexpr std::size_t kCellChunk = 32U * 1024U;
            core::parallelChunksWithThreads(
                0, static_cast<std::size_t>(counts.cells), kCellChunk,
                request.executionThreads,
                [&](std::size_t, std::size_t lo, std::size_t hi) {
                    std::uint64_t ordinal = lo;
                    std::uint64_t k = ordinal / cellPlane;
                    std::uint64_t remainder = ordinal - k * cellPlane;
                    std::uint64_t j = remainder / g.cellsU;
                    std::uint64_t i = remainder - j * g.cellsU;
                    Index base = static_cast<Index>(
                        (k * nodeV + j) * nodeU + i);
                    const Index row = static_cast<Index>(nodeU);
                    const Index plane = static_cast<Index>(nodePlane);
                    for (; ordinal < hi; ++ordinal) {
                        const std::array<Index, 8> conn{
                            base, base + 1, base + row + 1, base + row,
                            base + plane, base + plane + 1,
                            base + plane + row + 1, base + plane + row};
                        storage.setCell(
                            static_cast<std::size_t>(ordinal),
                            static_cast<std::size_t>(ordinal * 8U),
                            mesh::CellType::Hex8, conn,
                            g.volumeClassification);
                        ++base;
                        if (++i == g.cellsU) {
                            i = 0;
                            ++j;
                            ++base;
                            if (j == g.cellsV) {
                                j = 0;
                                ++k;
                                base += static_cast<Index>(nodeU);
                            }
                        }
                    }
                });
        }

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
                        {workspaceNodeIndex(i, j, 0, g.cellsU, g.cellsV),
                         workspaceNodeIndex(i, j + 1, 0, g.cellsU, g.cellsV),
                         workspaceNodeIndex(i + 1, j + 1, 0, g.cellsU, g.cellsV),
                         workspaceNodeIndex(i + 1, j, 0, g.cellsU, g.cellsV)},
                        StructuredBlockFace::WMin);
                    setQuad(0, ordinal, true,
                        {workspaceNodeIndex(i, j, g.cellsW, g.cellsU, g.cellsV),
                         workspaceNodeIndex(i + 1, j, g.cellsW, g.cellsU, g.cellsV),
                         workspaceNodeIndex(i + 1, j + 1, g.cellsW,
                                   g.cellsU, g.cellsV),
                         workspaceNodeIndex(i, j + 1, g.cellsW, g.cellsU, g.cellsV)},
                        StructuredBlockFace::WMax);
                });

            core::parallelForWithThreads(
                0, static_cast<std::size_t>(uw), request.executionThreads,
                [&](std::size_t ordinal) {
                    const std::uint64_t k = ordinal / g.cellsU;
                    const std::uint64_t i = ordinal % g.cellsU;
                    setQuad(2U * uv, ordinal, false,
                        {workspaceNodeIndex(i, 0, k, g.cellsU, g.cellsV),
                         workspaceNodeIndex(i + 1, 0, k, g.cellsU, g.cellsV),
                         workspaceNodeIndex(i + 1, 0, k + 1, g.cellsU, g.cellsV),
                         workspaceNodeIndex(i, 0, k + 1, g.cellsU, g.cellsV)},
                        StructuredBlockFace::VMin);
                    setQuad(2U * uv, ordinal, true,
                        {workspaceNodeIndex(i + 1, g.cellsV, k, g.cellsU, g.cellsV),
                         workspaceNodeIndex(i, g.cellsV, k, g.cellsU, g.cellsV),
                         workspaceNodeIndex(i, g.cellsV, k + 1, g.cellsU, g.cellsV),
                         workspaceNodeIndex(i + 1, g.cellsV, k + 1,
                                   g.cellsU, g.cellsV)},
                        StructuredBlockFace::VMax);
                });

            core::parallelForWithThreads(
                0, static_cast<std::size_t>(vw), request.executionThreads,
                [&](std::size_t ordinal) {
                    const std::uint64_t k = ordinal / g.cellsV;
                    const std::uint64_t j = ordinal % g.cellsV;
                    setQuad(2U * uv + 2U * uw, ordinal, false,
                        {workspaceNodeIndex(0, j + 1, k, g.cellsU, g.cellsV),
                         workspaceNodeIndex(0, j, k, g.cellsU, g.cellsV),
                         workspaceNodeIndex(0, j, k + 1, g.cellsU, g.cellsV),
                         workspaceNodeIndex(0, j + 1, k + 1, g.cellsU, g.cellsV)},
                        StructuredBlockFace::UMin);
                    setQuad(2U * uv + 2U * uw, ordinal, true,
                        {workspaceNodeIndex(g.cellsU, j, k, g.cellsU, g.cellsV),
                         workspaceNodeIndex(g.cellsU, j + 1, k, g.cellsU, g.cellsV),
                         workspaceNodeIndex(g.cellsU, j + 1, k + 1,
                                   g.cellsU, g.cellsV),
                         workspaceNodeIndex(g.cellsU, j, k + 1,
                                   g.cellsU, g.cellsV)},
                        StructuredBlockFace::UMax);
                });
        }

        out.metrics.constructionWorkUnits = counts.nodes + counts.cells +
            counts.boundaryFaces + counts.connectivity +
            out.metrics.growthPreparationWorkUnits;
        out.phase = MappedHexPhase::Audit;
        if (static_cast<std::uint64_t>(storage.nodeCount()) != counts.nodes ||
            static_cast<std::uint64_t>(storage.cellCount()) !=
                counts.cells + counts.boundaryFaces ||
            static_cast<std::uint64_t>(storage.connectivityCount()) !=
                counts.connectivity) {
            out.status = MappedHexStatus::Failed;
            out.reason = MappedHexReason::AuditFailed;
            return out;
        }

        if (request.validation ==
                MappedHexValidationMode::FullTopologyAudit) {
            mesh::Mesh local;
            storage.publishTo(local);
            out.metrics.audit = auditHex8Mesh(
                local, request.emitBoundaryQuads);
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
            storage.acquireFrom(local);
        } else {
            out.metrics.audit = workspaceCertifiedAudit(
                g, counts, determinant, request.emitBoundaryQuads,
                uGrowth, vGrowth, wGrowth);
            if (!out.metrics.audit.complete || !out.metrics.audit.valid ||
                std::abs(out.metrics.audit.totalVolume - determinant) >
                    128.0 * (std::numeric_limits<double>::epsilon)() *
                    std::max(1.0, std::abs(determinant))) {
                out.status = MappedHexStatus::Failed;
                out.reason = out.metrics.audit.allocationFailed
                    ? MappedHexReason::AllocationFailed
                    : MappedHexReason::AuditFailed;
                return out;
            }
            out.phase = MappedHexPhase::Publication;
            storage.publishTo(destination);
        }
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
