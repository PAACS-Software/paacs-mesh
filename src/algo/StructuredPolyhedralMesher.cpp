// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/core/SmallSort.hpp"
#include "paacs/algo/StructuredPolyhedralMesher.hpp"

#include "paacs/topo/ElementTraits.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <vector>

namespace paacs::algo {
namespace {

struct VariableFaceKey {
    std::array<Index, 4> nodes{kInvalidIndex, kInvalidIndex,
                               kInvalidIndex, kInvalidIndex};
    std::uint8_t count = 0;
    friend bool operator<(const VariableFaceKey& a,
                          const VariableFaceKey& b) noexcept {
        return a.count < b.count ||
            (a.count == b.count && a.nodes < b.nodes);
    }
};

bool finite(const Vec3& value) noexcept {
    return std::isfinite(value.x) && std::isfinite(value.y) &&
           std::isfinite(value.z);
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

struct Counts {
    std::uint64_t nodes = 0;
    std::uint64_t cells = 0;
    std::uint64_t uFaces = 0;
    std::uint64_t vFaces = 0;
    std::uint64_t wFaces = 0;
    std::uint64_t faces = 0;
    std::uint64_t faceNodes = 0;
    std::uint64_t cellFaces = 0;
};

bool countsFor(const StructuredBlockGeometry& g, Counts& out) noexcept {
    const std::uint64_t u = g.cellsU;
    const std::uint64_t v = g.cellsV;
    const std::uint64_t w = g.cellsW;
    std::uint64_t uv = 0, nodeUv = 0;
    if (!checkedMul(u, v, uv) || !checkedMul(uv, w, out.cells) ||
        !checkedMul(u + 1, v + 1, nodeUv) ||
        !checkedMul(nodeUv, w + 1, out.nodes) ||
        !checkedMul(u + 1, v, out.uFaces) ||
        !checkedMul(out.uFaces, w, out.uFaces) ||
        !checkedMul(u, v + 1, out.vFaces) ||
        !checkedMul(out.vFaces, w, out.vFaces) ||
        !checkedMul(u, v, out.wFaces) ||
        !checkedMul(out.wFaces, w + 1, out.wFaces))
        return false;
    std::uint64_t subtotal = 0;
    if (!checkedAdd(out.uFaces, out.vFaces, subtotal) ||
        !checkedAdd(subtotal, out.wFaces, out.faces) ||
        !checkedMul(out.faces, 4, out.faceNodes) ||
        !checkedMul(out.cells, 6, out.cellFaces))
        return false;
    return true;
}

Index nodeIndex(std::uint64_t i, std::uint64_t j, std::uint64_t k,
                std::uint64_t u, std::uint64_t v) noexcept {
    return static_cast<Index>((k * (v + 1) + j) * (u + 1) + i);
}

Index cellIndex(std::uint64_t i, std::uint64_t j, std::uint64_t k,
                std::uint64_t u, std::uint64_t v) noexcept {
    return static_cast<Index>((k * v + j) * u + i);
}

Index uFaceIndex(std::uint64_t i, std::uint64_t j, std::uint64_t k,
                 std::uint64_t u, std::uint64_t v) noexcept {
    return static_cast<Index>((k * v + j) * (u + 1) + i);
}

Index vFaceIndex(std::uint64_t i, std::uint64_t j, std::uint64_t k,
                 const Counts& counts, std::uint64_t u,
                 std::uint64_t v) noexcept {
    return static_cast<Index>(counts.uFaces +
        (k * (v + 1) + j) * u + i);
}

Index wFaceIndex(std::uint64_t i, std::uint64_t j, std::uint64_t k,
                 const Counts& counts, std::uint64_t u,
                 std::uint64_t v) noexcept {
    return static_cast<Index>(counts.uFaces + counts.vFaces +
        (k * v + j) * u + i);
}


VariableFaceKey variableKeyOf(std::span<const Index> values) {
    if (values.size() > 4) throw std::invalid_argument("face key exceeds four nodes");
    VariableFaceKey out;
    out.count = static_cast<std::uint8_t>(values.size());
    std::copy(values.begin(), values.end(), out.nodes.begin());
    (void)core::sortSmallPrefix(out.nodes, out.count);
    return out;
}

mesh::PolyhedralAuditReport certifiedAudit(
        const StructuredBlockGeometry& g, const Counts& counts,
        double determinant,
        const StructuredAxisDistribution& uGrowth,
        const StructuredAxisDistribution& vGrowth,
        const StructuredAxisDistribution& wGrowth) noexcept {
    mesh::PolyhedralAuditReport out;
    out.nodes = counts.nodes;
    out.faces = counts.faces;
    out.cells = counts.cells;
    const std::uint64_t uv =
        static_cast<std::uint64_t>(g.cellsU) * g.cellsV;
    const std::uint64_t uw =
        static_cast<std::uint64_t>(g.cellsU) * g.cellsW;
    const std::uint64_t vw =
        static_cast<std::uint64_t>(g.cellsV) * g.cellsW;
    out.boundaryFaces = 2U * (uv + uw + vw);
    out.interiorFaces = counts.faces - out.boundaryFaces;
    const double minimumCellVolume = determinant *
        uGrowth.minimumCellFraction * vGrowth.minimumCellFraction *
        wGrowth.minimumCellFraction;
    const double maximumCellVolume = determinant *
        uGrowth.maximumCellFraction * vGrowth.maximumCellFraction *
        wGrowth.maximumCellFraction;
    out.minimumCellVolume = minimumCellVolume;
    out.maximumCellVolume = maximumCellVolume;
    out.totalVolume = determinant;
    out.valid = std::isfinite(minimumCellVolume) &&
        minimumCellVolume > 0.0 &&
        std::isfinite(maximumCellVolume) && maximumCellVolume > 0.0;
    return out;
}

StructuredPolyOutcome rejected(StructuredPolyReason reason) noexcept {
    StructuredPolyOutcome out;
    out.status = StructuredPolyStatus::Rejected;
    out.reason = reason;
    out.phase = StructuredPolyPhase::Preflight;
    return out;
}

bool validGeometry(const StructuredBlockGeometry& g,
                   StructuredPolyReason& reason,
                   Counts& counts,
                   double& determinant) noexcept {
    if (g.cellsU == 0 || g.cellsV == 0 || g.cellsW == 0) {
        reason = StructuredPolyReason::InvalidCellCounts;
        return false;
    }
    if (!validClass(g.volumeClassification, 3)) {
        reason = StructuredPolyReason::InvalidClassification;
        return false;
    }
    for (const auto cls : g.vertexClassifications)
        if (!validClass(cls, 0)) {
            reason = StructuredPolyReason::InvalidClassification;
            return false;
        }
    for (const auto cls : g.edgeClassifications)
        if (!validClass(cls, 1)) {
            reason = StructuredPolyReason::InvalidClassification;
            return false;
        }
    for (const auto cls : g.boundaryClassifications)
        if (!validClass(cls, 2)) {
            reason = StructuredPolyReason::InvalidClassification;
            return false;
        }
    if (!finite(g.origin) || !finite(g.axisU) || !finite(g.axisV) ||
        !finite(g.axisW)) {
        reason = StructuredPolyReason::NonfiniteFrame;
        return false;
    }
    determinant = g.axisU.dot(g.axisV.cross(g.axisW));
    if (!std::isfinite(determinant) || determinant <= 0.0) {
        reason = StructuredPolyReason::NonPositiveFrame;
        return false;
    }
    if (!countsFor(g, counts)) {
        reason = StructuredPolyReason::CountOverflow;
        return false;
    }
    if (counts.nodes > g.maxNodes || counts.cells > g.maxCells ||
        counts.faces > g.maxFaces) {
        reason = StructuredPolyReason::BudgetExceeded;
        return false;
    }
    const auto indexMax = static_cast<std::uint64_t>(
        (std::numeric_limits<Index>::max)());
    if (counts.nodes > indexMax || counts.cells > indexMax ||
        counts.faces > indexMax ||
        counts.nodes > (std::numeric_limits<std::size_t>::max)() ||
        counts.cells > (std::numeric_limits<std::size_t>::max)() ||
        counts.faces > (std::numeric_limits<std::size_t>::max)() ||
        counts.faceNodes > (std::numeric_limits<std::size_t>::max)() ||
        counts.cellFaces > (std::numeric_limits<std::size_t>::max)()) {
        reason = StructuredPolyReason::CountOverflow;
        return false;
    }
    return true;
}

}  // namespace

StructuredPolyOutcome generateStructuredPolyhedralBlock(
        const StructuredPolyBlockRequest& request,
        mesh::PolyhedralMesh& destination) noexcept {
    const auto& g = request.geometry;
    const auto validation = request.validation;
    Counts counts;
    double determinant = 0.0;
    StructuredPolyReason preflight = StructuredPolyReason::None;
    if (!validGeometry(g, preflight, counts, determinant))
        return rejected(preflight);

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
            StructuredPolyOutcome failed;
            failed.status = StructuredPolyStatus::Failed;
            failed.reason = StructuredPolyReason::AllocationFailed;
            failed.phase = StructuredPolyPhase::Preflight;
            return failed;
        }
        return rejected(StructuredPolyReason::InvalidGrowth);
    }

    StructuredPolyOutcome out;
    out.phase = StructuredPolyPhase::Construction;
    out.metrics.nodes = counts.nodes;
    out.metrics.faces = counts.faces;
    out.metrics.cells = counts.cells;
    out.metrics.faceNodeEntries = counts.faceNodes;
    out.metrics.cellFaceEntries = counts.cellFaces;
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
        counts.faces * sizeof(mesh::PolyFace) +
        counts.faceNodes * sizeof(Index) +
        counts.cells * sizeof(mesh::PolyCell) +
        counts.cellFaces * sizeof(mesh::OrientedPolyFace);

    try {
        mesh::PolyhedralMesh local;
        local.reserve(static_cast<std::size_t>(counts.nodes),
                      static_cast<std::size_t>(counts.faces),
                      static_cast<std::size_t>(counts.faceNodes),
                      static_cast<std::size_t>(counts.cells),
                      static_cast<std::size_t>(counts.cellFaces));
        for (std::uint64_t k = 0; k <= g.cellsW; ++k) {
            const double wk = wGrowth.coordinate(
                static_cast<std::uint32_t>(k));
            for (std::uint64_t j = 0; j <= g.cellsV; ++j) {
                const double vj = vGrowth.coordinate(
                    static_cast<std::uint32_t>(j));
                for (std::uint64_t i = 0; i <= g.cellsU; ++i) {
                    const double ui = uGrowth.coordinate(
                        static_cast<std::uint32_t>(i));
                    (void)local.addNode(
                        g.origin + g.axisU * ui + g.axisV * vj +
                            g.axisW * wk,
                        structuredBlockNodeClassification(g, i, j, k));
                }
            }
        }

        std::array<Index, 4> face{};
        // U-normal faces.
        for (std::uint64_t k = 0; k < g.cellsW; ++k)
            for (std::uint64_t j = 0; j < g.cellsV; ++j)
                for (std::uint64_t i = 0; i <= g.cellsU; ++i) {
                    Index owner = 0, neighbour = kInvalidIndex;
                    mesh::Classification cls = g.volumeClassification;
                    if (i == 0) {
                        owner = cellIndex(0, j, k, g.cellsU, g.cellsV);
                        face = {nodeIndex(0, j + 1, k, g.cellsU, g.cellsV),
                                nodeIndex(0, j, k, g.cellsU, g.cellsV),
                                nodeIndex(0, j, k + 1,
                                          g.cellsU, g.cellsV),
                                nodeIndex(0, j + 1, k + 1,
                                          g.cellsU, g.cellsV)};
                        cls = g.boundaryClassifications[
                            static_cast<std::size_t>(
                                StructuredBlockFace::UMin)];
                    } else {
                        owner = cellIndex(i - 1, j, k,
                                          g.cellsU, g.cellsV);
                        if (i < g.cellsU)
                            neighbour = cellIndex(i, j, k,
                                                  g.cellsU, g.cellsV);
                        face = {nodeIndex(i, j, k, g.cellsU, g.cellsV),
                                nodeIndex(i, j + 1, k,
                                          g.cellsU, g.cellsV),
                                nodeIndex(i, j + 1, k + 1,
                                          g.cellsU, g.cellsV),
                                nodeIndex(i, j, k + 1,
                                          g.cellsU, g.cellsV)};
                        if (i == g.cellsU)
                            cls = g.boundaryClassifications[
                                static_cast<std::size_t>(
                                    StructuredBlockFace::UMax)];
                    }
                    const Index created = local.addFace(
                        face, owner, neighbour, cls);
                    if (created != uFaceIndex(
                            i, j, k, g.cellsU, g.cellsV))
                        throw std::logic_error("U face ordinal mismatch");
                }

        // V-normal faces.
        for (std::uint64_t k = 0; k < g.cellsW; ++k)
            for (std::uint64_t j = 0; j <= g.cellsV; ++j)
                for (std::uint64_t i = 0; i < g.cellsU; ++i) {
                    Index owner = 0, neighbour = kInvalidIndex;
                    mesh::Classification cls = g.volumeClassification;
                    if (j == 0) {
                        owner = cellIndex(i, 0, k, g.cellsU, g.cellsV);
                        face = {nodeIndex(i, 0, k, g.cellsU, g.cellsV),
                                nodeIndex(i + 1, 0, k,
                                          g.cellsU, g.cellsV),
                                nodeIndex(i + 1, 0, k + 1,
                                          g.cellsU, g.cellsV),
                                nodeIndex(i, 0, k + 1,
                                          g.cellsU, g.cellsV)};
                        cls = g.boundaryClassifications[
                            static_cast<std::size_t>(
                                StructuredBlockFace::VMin)];
                    } else {
                        owner = cellIndex(i, j - 1, k,
                                          g.cellsU, g.cellsV);
                        if (j < g.cellsV)
                            neighbour = cellIndex(i, j, k,
                                                  g.cellsU, g.cellsV);
                        face = {nodeIndex(i + 1, j, k,
                                          g.cellsU, g.cellsV),
                                nodeIndex(i, j, k, g.cellsU, g.cellsV),
                                nodeIndex(i, j, k + 1,
                                          g.cellsU, g.cellsV),
                                nodeIndex(i + 1, j, k + 1,
                                          g.cellsU, g.cellsV)};
                        if (j == g.cellsV)
                            cls = g.boundaryClassifications[
                                static_cast<std::size_t>(
                                    StructuredBlockFace::VMax)];
                    }
                    const Index created = local.addFace(
                        face, owner, neighbour, cls);
                    if (created != vFaceIndex(
                            i, j, k, counts, g.cellsU, g.cellsV))
                        throw std::logic_error("V face ordinal mismatch");
                }

        // W-normal faces.
        for (std::uint64_t k = 0; k <= g.cellsW; ++k)
            for (std::uint64_t j = 0; j < g.cellsV; ++j)
                for (std::uint64_t i = 0; i < g.cellsU; ++i) {
                    Index owner = 0, neighbour = kInvalidIndex;
                    mesh::Classification cls = g.volumeClassification;
                    if (k == 0) {
                        owner = cellIndex(i, j, 0, g.cellsU, g.cellsV);
                        face = {nodeIndex(i, j, 0, g.cellsU, g.cellsV),
                                nodeIndex(i, j + 1, 0,
                                          g.cellsU, g.cellsV),
                                nodeIndex(i + 1, j + 1, 0,
                                          g.cellsU, g.cellsV),
                                nodeIndex(i + 1, j, 0,
                                          g.cellsU, g.cellsV)};
                        cls = g.boundaryClassifications[
                            static_cast<std::size_t>(
                                StructuredBlockFace::WMin)];
                    } else {
                        owner = cellIndex(i, j, k - 1,
                                          g.cellsU, g.cellsV);
                        if (k < g.cellsW)
                            neighbour = cellIndex(i, j, k,
                                                  g.cellsU, g.cellsV);
                        face = {nodeIndex(i, j, k, g.cellsU, g.cellsV),
                                nodeIndex(i + 1, j, k,
                                          g.cellsU, g.cellsV),
                                nodeIndex(i + 1, j + 1, k,
                                          g.cellsU, g.cellsV),
                                nodeIndex(i, j + 1, k,
                                          g.cellsU, g.cellsV)};
                        if (k == g.cellsW)
                            cls = g.boundaryClassifications[
                                static_cast<std::size_t>(
                                    StructuredBlockFace::WMax)];
                    }
                    const Index created = local.addFace(
                        face, owner, neighbour, cls);
                    if (created != wFaceIndex(
                            i, j, k, counts, g.cellsU, g.cellsV))
                        throw std::logic_error("W face ordinal mismatch");
                }

        std::array<mesh::OrientedPolyFace, 6> refs{};
        for (std::uint64_t k = 0; k < g.cellsW; ++k)
            for (std::uint64_t j = 0; j < g.cellsV; ++j)
                for (std::uint64_t i = 0; i < g.cellsU; ++i) {
                    refs = {{
                        {wFaceIndex(i, j, k, counts,
                                    g.cellsU, g.cellsV),
                         static_cast<std::int8_t>(k == 0 ? 1 : -1)},
                        {wFaceIndex(i, j, k + 1, counts,
                                    g.cellsU, g.cellsV), 1},
                        {vFaceIndex(i, j, k, counts,
                                    g.cellsU, g.cellsV),
                         static_cast<std::int8_t>(j == 0 ? 1 : -1)},
                        {uFaceIndex(i + 1, j, k,
                                    g.cellsU, g.cellsV), 1},
                        {vFaceIndex(i, j + 1, k, counts,
                                    g.cellsU, g.cellsV), 1},
                        {uFaceIndex(i, j, k,
                                    g.cellsU, g.cellsV),
                         static_cast<std::int8_t>(i == 0 ? 1 : -1)}}};
                    const Index created = local.addCell(
                        refs, g.volumeClassification);
                    if (created != cellIndex(
                            i, j, k, g.cellsU, g.cellsV))
                        throw std::logic_error("cell ordinal mismatch");
                }

        out.metrics.constructionWorkUnits = counts.nodes + counts.faces +
            counts.cells + counts.faceNodes + counts.cellFaces +
            out.metrics.growthPreparationWorkUnits;
        out.phase = StructuredPolyPhase::Audit;
        if (static_cast<std::uint64_t>(local.nodeCount()) != counts.nodes ||
            static_cast<std::uint64_t>(local.faceCount()) != counts.faces ||
            static_cast<std::uint64_t>(local.cellCount()) != counts.cells) {
            out.status = StructuredPolyStatus::Failed;
            out.reason = StructuredPolyReason::AuditFailed;
            return out;
        }
        out.metrics.audit = validation ==
                StructuredPolyValidationMode::FullAudit
            ? mesh::auditPolyhedralMesh(local)
            : certifiedAudit(g, counts, determinant,
                             uGrowth, vGrowth, wGrowth);
        if (!out.metrics.audit.complete) {
            out.status = StructuredPolyStatus::Failed;
            out.reason = out.metrics.audit.allocationFailed
                ? StructuredPolyReason::AllocationFailed
                : StructuredPolyReason::AuditFailed;
            return out;
        }
        if (!out.metrics.audit.valid ||
            std::abs(out.metrics.audit.totalVolume - determinant) >
                128.0 * (std::numeric_limits<double>::epsilon)() *
                std::max(1.0, std::abs(determinant))) {
            out.status = StructuredPolyStatus::Failed;
            out.reason = StructuredPolyReason::AuditFailed;
            return out;
        }
        out.phase = StructuredPolyPhase::Publication;
        destination.swap(local);
        out.status = StructuredPolyStatus::Published;
        out.reason = StructuredPolyReason::None;
        out.destinationUnchanged = false;
        return out;
    } catch (const std::bad_alloc&) {
        out.status = StructuredPolyStatus::Failed;
        out.reason = StructuredPolyReason::AllocationFailed;
    } catch (...) {
        out.status = StructuredPolyStatus::Failed;
        out.reason = StructuredPolyReason::InternalError;
    }
    return out;
}

StructuredPolyOutcome generateStructuredPolyhedralBlock(
        const StructuredBlockGeometry& geometry,
        mesh::PolyhedralMesh& destination,
        StructuredPolyValidationMode validation) noexcept {
    StructuredPolyBlockRequest request;
    request.geometry = geometry;
    request.validation = validation;
    return generateStructuredPolyhedralBlock(request, destination);
}

namespace {

StructuredPolyOutcome convertVolumeMeshImpl(
        const mesh::Mesh& source,
        mesh::PolyhedralMesh& destination,
        bool allowPrisms) noexcept {
    StructuredPolyOutcome out;
    out.phase = StructuredPolyPhase::Preflight;
    try {
        std::uint64_t volumeCells = 0;
        std::uint64_t surfaceCells = 0;
        std::uint64_t faceCapacity = 0;
        std::uint64_t faceNodeCapacity = 0;
        std::uint64_t cellFaceCapacity = 0;
        for (Index ci = 0; ci < source.cellCount(); ++ci) {
            const auto type = source.cell(ci).type;
            if (type == mesh::CellType::Hex8) {
                ++volumeCells;
                faceCapacity += 6U;
                faceNodeCapacity += 24U;
                cellFaceCapacity += 6U;
            } else if (allowPrisms && type == mesh::CellType::Prism6) {
                ++volumeCells;
                faceCapacity += 5U;
                faceNodeCapacity += 18U;
                cellFaceCapacity += 5U;
            } else if (allowPrisms && type == mesh::CellType::Pyramid5) {
                ++volumeCells;
                faceCapacity += 5U;
                faceNodeCapacity += 16U;
                cellFaceCapacity += 5U;
            } else if (type == mesh::CellType::Quad4 ||
                       (allowPrisms && type == mesh::CellType::Tri3)) {
                ++surfaceCells;
            } else {
                return rejected(StructuredPolyReason::UnsupportedCellType);
            }
        }
        if (volumeCells == 0)
            return rejected(StructuredPolyReason::UnsupportedCellType);
        if (faceCapacity > (std::numeric_limits<std::size_t>::max)() ||
            faceNodeCapacity > (std::numeric_limits<std::size_t>::max)() ||
            cellFaceCapacity > (std::numeric_limits<std::size_t>::max)() ||
            volumeCells > (std::numeric_limits<std::size_t>::max)())
            return rejected(StructuredPolyReason::CountOverflow);

        out.phase = StructuredPolyPhase::Conversion;
        mesh::PolyhedralMesh local;
        local.reserve(static_cast<std::size_t>(source.nodeCount()),
                      static_cast<std::size_t>(faceCapacity),
                      static_cast<std::size_t>(faceNodeCapacity),
                      static_cast<std::size_t>(volumeCells),
                      static_cast<std::size_t>(cellFaceCapacity));
        for (Index i = 0; i < source.nodeCount(); ++i)
            (void)local.addNode(source.node(i), source.nodeClass(i));

        std::map<VariableFaceKey, mesh::Classification> boundaryClasses;
        for (Index ci = 0; ci < source.cellCount(); ++ci) {
            const auto& cell = source.cell(ci);
            if (cell.type != mesh::CellType::Quad4 &&
                cell.type != mesh::CellType::Tri3)
                continue;
            if (!validClass(cell.onEntity, 2))
                return rejected(StructuredPolyReason::InvalidClassification);
            const std::size_t n = static_cast<std::size_t>(
                mesh::nodesPerCell(cell.type));
            const Index* src = source.cellNodes(ci);
            std::array<Index, 4> nodes{};
            std::copy(src, src + n, nodes.begin());
            const auto [it, inserted] = boundaryClasses.emplace(
                variableKeyOf(std::span<const Index>(nodes.data(), n)),
                cell.onEntity);
            if (!inserted)
                return rejected(StructuredPolyReason::NonmanifoldInput);
            (void)it;
        }

        std::map<VariableFaceKey, Index> faces;
        std::uint64_t polyCell = 0;
        for (Index ci = 0; ci < source.cellCount(); ++ci) {
            const auto& cell = source.cell(ci);
            topo::ElementType element{};
            if (cell.type == mesh::CellType::Hex8)
                element = topo::ElementType::Hex8;
            else if (cell.type == mesh::CellType::Prism6 && allowPrisms)
                element = topo::ElementType::Wedge6;
            else if (cell.type == mesh::CellType::Pyramid5 && allowPrisms)
                element = topo::ElementType::Pyramid5;
            else
                continue;
            if (!validClass(cell.onEntity, 3))
                return rejected(StructuredPolyReason::InvalidClassification);
            const Index* src = source.cellNodes(ci);
            const auto& traits = topo::traits(element);
            std::vector<mesh::OrientedPolyFace> refs(
                static_cast<std::size_t>(traits.faceCount));
            for (int fi = 0; fi < traits.faceCount; ++fi) {
                const auto& localFace = traits.faces[fi];
                std::array<Index, 4> ordered{};
                for (int n = 0; n < localFace.n; ++n)
                    ordered[static_cast<std::size_t>(n)] =
                        src[localFace.v[static_cast<std::size_t>(n)]];
                const auto key = variableKeyOf(std::span<const Index>(
                    ordered.data(), static_cast<std::size_t>(localFace.n)));
                const auto found = faces.find(key);
                if (found == faces.end()) {
                    const Index face = local.addFace(
                        std::span<const Index>(ordered.data(),
                            static_cast<std::size_t>(localFace.n)),
                        static_cast<Index>(polyCell), kInvalidIndex,
                        cell.onEntity);
                    faces.emplace(key, face);
                    refs[static_cast<std::size_t>(fi)] = {face, 1};
                } else {
                    const Index face = found->second;
                    if (local.face(face).neighbour != kInvalidIndex)
                        return rejected(
                            StructuredPolyReason::NonmanifoldInput);
                    local.setFaceNeighbour(face,
                                           static_cast<Index>(polyCell));
                    refs[static_cast<std::size_t>(fi)] = {face, -1};
                }
            }
            const Index created = local.addCell(refs, cell.onEntity);
            if (created != static_cast<Index>(polyCell))
                throw std::logic_error("poly conversion ordinal mismatch");
            ++polyCell;
        }

        std::uint64_t matchedBoundaryCells = 0;
        for (const auto& [key, faceIndex] : faces) {
            const auto& face = local.face(faceIndex);
            if (face.neighbour != kInvalidIndex) {
                local.setFaceClassification(
                    faceIndex, local.cell(face.owner).onEntity);
                if (boundaryClasses.find(key) != boundaryClasses.end())
                    return rejected(StructuredPolyReason::NonmanifoldInput);
                continue;
            }
            const auto boundary = boundaryClasses.find(key);
            if (boundary != boundaryClasses.end()) {
                local.setFaceClassification(faceIndex, boundary->second);
                ++matchedBoundaryCells;
            } else {
                const auto ownerClass = local.cell(face.owner).onEntity;
                local.setFaceClassification(faceIndex, {2, ownerClass.tag});
            }
        }
        if (matchedBoundaryCells != surfaceCells)
            return rejected(StructuredPolyReason::NonmanifoldInput);

        out.metrics.nodes = static_cast<std::uint64_t>(local.nodeCount());
        out.metrics.faces = static_cast<std::uint64_t>(local.faceCount());
        out.metrics.cells = static_cast<std::uint64_t>(local.cellCount());
        for (Index fi = 0; fi < local.faceCount(); ++fi)
            out.metrics.faceNodeEntries += local.face(fi).nodeCount;
        for (Index ci = 0; ci < local.cellCount(); ++ci)
            out.metrics.cellFaceEntries += local.cell(ci).faceCount;
        out.metrics.constructionWorkUnits = out.metrics.nodes +
            out.metrics.faces + out.metrics.cells +
            out.metrics.faceNodeEntries + out.metrics.cellFaceEntries;
        out.metrics.logicalBytes = out.metrics.nodes *
                (sizeof(Vec3) + sizeof(mesh::Classification)) +
            out.metrics.faces * sizeof(mesh::PolyFace) +
            out.metrics.faceNodeEntries * sizeof(Index) +
            out.metrics.cells * sizeof(mesh::PolyCell) +
            out.metrics.cellFaceEntries * sizeof(mesh::OrientedPolyFace);

        out.phase = StructuredPolyPhase::Audit;
        out.metrics.audit = mesh::auditPolyhedralMesh(local);
        if (!out.metrics.audit.complete) {
            out.status = StructuredPolyStatus::Failed;
            out.reason = out.metrics.audit.allocationFailed
                ? StructuredPolyReason::AllocationFailed
                : StructuredPolyReason::AuditFailed;
            return out;
        }
        if (!out.metrics.audit.valid) {
            out.status = StructuredPolyStatus::Failed;
            out.reason = StructuredPolyReason::AuditFailed;
            return out;
        }
        out.phase = StructuredPolyPhase::Publication;
        destination.swap(local);
        out.status = StructuredPolyStatus::Published;
        out.reason = StructuredPolyReason::None;
        out.destinationUnchanged = false;
        return out;
    } catch (const std::bad_alloc&) {
        out.status = StructuredPolyStatus::Failed;
        out.reason = StructuredPolyReason::AllocationFailed;
    } catch (...) {
        out.status = StructuredPolyStatus::Failed;
        out.reason = StructuredPolyReason::InternalError;
    }
    return out;
}

}  // namespace

StructuredPolyOutcome convertHex8ToPolyhedral(
        const mesh::Mesh& source,
        mesh::PolyhedralMesh& destination) noexcept {
    return convertVolumeMeshImpl(source, destination, false);
}

StructuredPolyOutcome convertSupportedVolumeMeshToPolyhedral(
        const mesh::Mesh& source,
        mesh::PolyhedralMesh& destination) noexcept {
    return convertVolumeMeshImpl(source, destination, true);
}

}  // namespace paacs::algo
