// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/core/SmallSort.hpp"
#include "paacs/algo/StraightPrismLayerMesher.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/CertifiedMeshStorage.hpp"

#include "paacs/topo/ElementTraits.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <set>
#include <span>
#include <vector>

namespace paacs::algo {
namespace {

struct EdgeKey {
    Index a = kInvalidIndex;
    Index b = kInvalidIndex;
    friend bool operator<(const EdgeKey& x, const EdgeKey& y) noexcept {
        return x.a < y.a || (x.a == y.a && x.b < y.b);
    }
};

struct EdgeUse {
    std::uint32_t count = 0;
    Index directedA = kInvalidIndex;
    Index directedB = kInvalidIndex;
};

struct FaceKey {
    std::array<Index, 4> nodes{kInvalidIndex, kInvalidIndex,
                               kInvalidIndex, kInvalidIndex};
    std::uint8_t count = 0;
    friend bool operator<(const FaceKey& x, const FaceKey& y) noexcept {
        if (x.count != y.count) return x.count < y.count;
        return x.nodes < y.nodes;
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

EdgeKey edgeKey(Index a, Index b) noexcept {
    return a < b ? EdgeKey{a, b} : EdgeKey{b, a};
}

FaceKey faceKey(std::span<const Index> nodes) {
    if (nodes.size() > 4) throw std::invalid_argument("face key exceeds four nodes");
    FaceKey out;
    out.count = static_cast<std::uint8_t>(nodes.size());
    std::copy(nodes.begin(), nodes.end(), out.nodes.begin());
    (void)core::sortSmallPrefix(out.nodes, out.count);
    return out;
}

struct SourceAnalysis {
    std::vector<std::array<Index, 3>> triangles;
    std::vector<mesh::Classification> triangleClasses;
    std::vector<std::array<Index, 2>> boundaryEdges;
    std::vector<std::uint8_t> boundaryNodes;
    double totalProjectedArea = 0.0;
    double minimumProjectedArea = 0.0;
    double maximumProjectedArea = 0.0;
};

PrismLayerReason analyzeSource(const mesh::Mesh& source,
                               const Vec3& direction,
                               SourceAnalysis& out) {
    if (source.nodeCount() == 0 || source.cellCount() == 0)
        return PrismLayerReason::EmptySource;
    out.boundaryNodes.assign(static_cast<std::size_t>(source.nodeCount()), 0);
    std::map<EdgeKey, EdgeUse> edges;
    double minimumArea = (std::numeric_limits<double>::infinity)();
    for (Index ci = 0; ci < source.cellCount(); ++ci) {
        const auto& cell = source.cell(ci);
        if (cell.type != mesh::CellType::Tri3)
            return PrismLayerReason::UnsupportedSourceCell;
        if (!validClass(cell.onEntity, 2))
            return PrismLayerReason::InvalidClassification;
        const Index* conn = source.cellNodes(ci);
        std::array<Index, 3> tri{conn[0], conn[1], conn[2]};
        std::array<Index, 3> sorted = tri;
        std::sort(sorted.begin(), sorted.end());
        if (sorted[0] < 0 || sorted[2] >= source.nodeCount() ||
            std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
            return PrismLayerReason::InvalidSourceConnectivity;
        const Vec3 a = source.node(tri[0]);
        const Vec3 b = source.node(tri[1]);
        const Vec3 c = source.node(tri[2]);
        if (!finite(a) || !finite(b) || !finite(c))
            return PrismLayerReason::InvalidSourceConnectivity;
        const double twiceArea = (b - a).cross(c - a).dot(direction);
        if (!std::isfinite(twiceArea) || !(twiceArea > 0.0))
            return PrismLayerReason::InconsistentSourceOrientation;
        const double area = 0.5 * twiceArea;
        minimumArea = std::min(minimumArea, area);
        out.maximumProjectedArea = std::max(out.maximumProjectedArea, area);
        out.totalProjectedArea += area;
        out.triangles.push_back(tri);
        out.triangleClasses.push_back(cell.onEntity);

        for (std::size_t e = 0; e < 3; ++e) {
            const Index p = tri[e];
            const Index q = tri[(e + 1U) % 3U];
            auto& use = edges[edgeKey(p, q)];
            if (use.count == 0) {
                use.directedA = p;
                use.directedB = q;
            } else if (use.count == 1) {
                if (use.directedA != q || use.directedB != p)
                    return PrismLayerReason::InconsistentSourceOrientation;
            }
            ++use.count;
            if (use.count > 2)
                return PrismLayerReason::NonmanifoldSource;
        }
    }
    for (const auto& [key, use] : edges) {
        (void)key;
        if (use.count == 1) {
            out.boundaryEdges.push_back({use.directedA, use.directedB});
            out.boundaryNodes[static_cast<std::size_t>(use.directedA)] = 1;
            out.boundaryNodes[static_cast<std::size_t>(use.directedB)] = 1;
        } else if (use.count != 2) {
            return PrismLayerReason::NonmanifoldSource;
        }
    }
    out.minimumProjectedArea = minimumArea;
    return PrismLayerReason::None;
}

struct Counts {
    std::uint64_t nodes = 0;
    std::uint64_t prisms = 0;
    std::uint64_t boundaryTriangles = 0;
    std::uint64_t boundaryQuads = 0;
    std::uint64_t cells = 0;
    std::uint64_t connectivity = 0;
};

bool countsFor(std::uint64_t sourceNodes, std::uint64_t sourceTriangles,
               std::uint64_t boundaryEdges, std::uint32_t layers,
               bool emitBoundary, Counts& out) noexcept {
    if (!checkedMul(sourceNodes, static_cast<std::uint64_t>(layers) + 1U,
                    out.nodes) ||
        !checkedMul(sourceTriangles, layers, out.prisms))
        return false;
    if (emitBoundary) {
        if (!checkedMul(sourceTriangles, 2U, out.boundaryTriangles) ||
            !checkedMul(boundaryEdges, layers, out.boundaryQuads))
            return false;
    }
    std::uint64_t boundary = 0;
    if (!checkedAdd(out.boundaryTriangles, out.boundaryQuads, boundary) ||
        !checkedAdd(out.prisms, boundary, out.cells))
        return false;
    std::uint64_t prismConn = 0, triConn = 0, quadConn = 0;
    if (!checkedMul(out.prisms, 6U, prismConn) ||
        !checkedMul(out.boundaryTriangles, 3U, triConn) ||
        !checkedMul(out.boundaryQuads, 4U, quadConn) ||
        !checkedAdd(prismConn, triConn, out.connectivity) ||
        !checkedAdd(out.connectivity, quadConn, out.connectivity))
        return false;
    return true;
}

PrismLayerOutcome rejected(PrismLayerReason reason) noexcept {
    PrismLayerOutcome out;
    out.status = PrismLayerStatus::Rejected;
    out.reason = reason;
    out.phase = PrismLayerPhase::Preflight;
    return out;
}

double prismVolume(const mesh::Mesh& value, const Index* conn) noexcept {
    const Vec3& a = value.node(conn[0]);
    const Vec3& b = value.node(conn[1]);
    const Vec3& c = value.node(conn[2]);
    const Vec3 d = ((value.node(conn[3]) - a) +
                    (value.node(conn[4]) - b) +
                    (value.node(conn[5]) - c)) / 3.0;
    return 0.5 * (b - a).cross(c - a).dot(d);
}

}  // namespace

PrismLayerAuditReport auditPrismLayerMesh(
        const mesh::Mesh& value, bool requireBoundaryFaces) noexcept {
    PrismLayerAuditReport out;
    out.nodes = static_cast<std::uint64_t>(value.nodeCount());
    try {
        struct FaceUse { std::uint32_t volumes = 0; std::uint32_t surfaces = 0; };
        std::map<FaceKey, FaceUse> faces;
        double minVolume = (std::numeric_limits<double>::infinity)();
        double maxVolume = 0.0;
        long double total = 0.0L;
        const auto& wedge = topo::traits(topo::ElementType::Wedge6);
        for (Index ci = 0; ci < value.cellCount(); ++ci) {
            const auto& cell = value.cell(ci);
            const Index* conn = value.cellNodes(ci);
            if (cell.type == mesh::CellType::Prism6) {
                ++out.prisms;
                if (!validClass(cell.onEntity, 3))
                    ++out.invalidClassifications;
                std::array<Index, 6> unique{};
                std::copy(conn, conn + 6, unique.begin());
                std::sort(unique.begin(), unique.end());
                bool bad = unique.front() < 0 ||
                    unique.back() >= value.nodeCount();
                if (std::adjacent_find(unique.begin(), unique.end()) !=
                    unique.end()) {
                    ++out.duplicateCellNodes;
                    bad = true;
                }
                if (bad) {
                    ++out.invalidConnectivity;
                    continue;
                }
                const double volume = prismVolume(value, conn);
                if (!std::isfinite(volume) || !(volume > 0.0)) {
                    ++out.nonPositiveVolumes;
                } else {
                    minVolume = std::min(minVolume, volume);
                    maxVolume = std::max(maxVolume, volume);
                    total += volume;
                }
                for (int fi = 0; fi < wedge.faceCount; ++fi) {
                    const auto& local = wedge.faces[fi];
                    std::array<Index, 4> nodes{};
                    for (int n = 0; n < local.n; ++n)
                        nodes[static_cast<std::size_t>(n)] =
                            conn[local.v[static_cast<std::size_t>(n)]];
                    ++faces[faceKey(std::span<const Index>(
                        nodes.data(), static_cast<std::size_t>(local.n)))].volumes;
                }
            } else if (cell.type == mesh::CellType::Tri3 ||
                       cell.type == mesh::CellType::Quad4) {
                if (!validClass(cell.onEntity, 2))
                    ++out.invalidClassifications;
                const std::size_t n = static_cast<std::size_t>(
                    mesh::nodesPerCell(cell.type));
                std::array<Index, 4> nodes{};
                std::copy(conn, conn + n, nodes.begin());
                ++faces[faceKey(std::span<const Index>(nodes.data(), n))].surfaces;
                if (cell.type == mesh::CellType::Tri3)
                    ++out.boundaryTriangles;
                else
                    ++out.boundaryQuads;
            } else {
                ++out.unsupportedCells;
            }
        }
        for (const auto& [key, use] : faces) {
            (void)key;
            if (use.volumes == 0) {
                if (use.surfaces != 0) ++out.unexpectedBoundaryFaces;
            } else if (use.volumes == 1) {
                if (requireBoundaryFaces && use.surfaces == 0)
                    ++out.missingBoundaryFaces;
                if (use.surfaces > 1) ++out.duplicateBoundaryFaces;
            } else if (use.volumes == 2) {
                ++out.interiorFaces;
                if (use.surfaces != 0) ++out.unexpectedBoundaryFaces;
            } else {
                ++out.nonmanifoldFaces;
            }
        }
        if (out.prisms == 0) minVolume = 0.0;
        out.minimumPrismVolume = minVolume;
        out.maximumPrismVolume = maxVolume;
        out.totalVolume = static_cast<double>(total);
        out.valid = out.prisms != 0 && out.unsupportedCells == 0 &&
            out.invalidConnectivity == 0 && out.duplicateCellNodes == 0 &&
            out.invalidClassifications == 0 && out.nonPositiveVolumes == 0 &&
            out.nonmanifoldFaces == 0 && out.missingBoundaryFaces == 0 &&
            out.duplicateBoundaryFaces == 0 &&
            out.unexpectedBoundaryFaces == 0;
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

PrismLayerOutcome generateStraightPrismLayers(
        const mesh::Mesh& sourceSurface,
        const StraightPrismLayerRequest& request,
        mesh::Mesh& destination) noexcept {
    if (!finite(request.direction))
        return rejected(PrismLayerReason::NonfiniteDirection);
    const double directionNorm = request.direction.norm();
    if (!std::isfinite(directionNorm) || !(directionNorm > 0.0))
        return rejected(PrismLayerReason::ZeroDirection);
    if (request.layers == 0)
        return rejected(PrismLayerReason::InvalidLayerCount);
    if (!std::isfinite(request.firstLayerThickness) ||
        !(request.firstLayerThickness > 0.0))
        return rejected(PrismLayerReason::InvalidFirstLayerThickness);
    if (!std::isfinite(request.cellToCellGrowth) ||
        request.cellToCellGrowth < 1.0)
        return rejected(PrismLayerReason::InvalidGrowthRate);
    if (!std::isfinite(request.maximumAllowedGrowth) ||
        request.maximumAllowedGrowth < 1.0 ||
        request.cellToCellGrowth > request.maximumAllowedGrowth)
        return rejected(PrismLayerReason::GrowthLimitExceeded);
    if (!std::isfinite(request.maximumTotalThickness) ||
        !(request.maximumTotalThickness > 0.0))
        return rejected(PrismLayerReason::ThicknessBudgetExceeded);
    if (!validClass(request.volumeClassification, 3) ||
        (request.emitBoundaryFaces &&
         (!validClass(request.targetClassification, 2) ||
          !validClass(request.sideClassification, 2))))
        return rejected(PrismLayerReason::InvalidClassification);

    const Vec3 direction = request.direction / directionNorm;
    SourceAnalysis source;
    PrismLayerReason sourceReason = PrismLayerReason::None;
    try {
        sourceReason = analyzeSource(sourceSurface, direction, source);
    } catch (const std::bad_alloc&) {
        PrismLayerOutcome failed;
        failed.reason = PrismLayerReason::AllocationFailed;
        failed.phase = PrismLayerPhase::Preflight;
        return failed;
    } catch (...) {
        PrismLayerOutcome failed;
        failed.reason = PrismLayerReason::InternalError;
        failed.phase = PrismLayerPhase::Preflight;
        return failed;
    }
    if (sourceReason != PrismLayerReason::None)
        return rejected(sourceReason);

    PrismLayerOutcome out;
    out.phase = PrismLayerPhase::Schedule;
    std::vector<double> offsets;
    double totalThickness = 0.0;
    try {
        offsets.resize(static_cast<std::size_t>(request.layers) + 1U);
        long double width = request.firstLayerThickness;
        long double total = 0.0L;
        for (std::uint32_t layer = 0; layer < request.layers; ++layer) {
            total += width;
            const double value = static_cast<double>(total);
            if (!std::isfinite(value) || !(value > offsets[layer]))
                return rejected(PrismLayerReason::InvalidGrowthRate);
            offsets[static_cast<std::size_t>(layer + 1U)] = value;
            width *= static_cast<long double>(request.cellToCellGrowth);
        }
        totalThickness = offsets.back();
    } catch (const std::bad_alloc&) {
        out.reason = PrismLayerReason::AllocationFailed;
        return out;
    } catch (...) {
        out.reason = PrismLayerReason::InternalError;
        return out;
    }
    if (!std::isfinite(totalThickness) ||
        totalThickness > request.maximumTotalThickness)
        return rejected(PrismLayerReason::ThicknessBudgetExceeded);

    Counts counts;
    if (!countsFor(static_cast<std::uint64_t>(sourceSurface.nodeCount()),
                   static_cast<std::uint64_t>(source.triangles.size()),
                   static_cast<std::uint64_t>(source.boundaryEdges.size()),
                   request.layers, request.emitBoundaryFaces, counts))
        return rejected(PrismLayerReason::CountOverflow);
    if (counts.nodes > request.maxNodes || counts.prisms > request.maxCells ||
        counts.boundaryTriangles + counts.boundaryQuads >
            request.maxBoundaryFaces)
        return rejected(PrismLayerReason::BudgetExceeded);
    const auto indexMax = static_cast<std::uint64_t>(
        (std::numeric_limits<Index>::max)());
    if (counts.nodes > indexMax || counts.cells > indexMax ||
        counts.nodes > (std::numeric_limits<std::size_t>::max)() ||
        counts.cells > (std::numeric_limits<std::size_t>::max)() ||
        counts.connectivity > (std::numeric_limits<std::size_t>::max)())
        return rejected(PrismLayerReason::CountOverflow);

    out.metrics.sourceNodes = static_cast<std::uint64_t>(
        sourceSurface.nodeCount());
    out.metrics.sourceTriangles = static_cast<std::uint64_t>(
        source.triangles.size());
    out.metrics.sourceBoundaryEdges = static_cast<std::uint64_t>(
        source.boundaryEdges.size());
    out.metrics.nodes = counts.nodes;
    out.metrics.prisms = counts.prisms;
    out.metrics.boundaryTriangles = counts.boundaryTriangles;
    out.metrics.boundaryQuads = counts.boundaryQuads;
    out.metrics.connectivityEntries = counts.connectivity;
    out.metrics.totalThickness = totalThickness;
    out.metrics.maximumAdjacentGrowth = request.cellToCellGrowth;
    out.metrics.logicalBytes = counts.nodes *
            (sizeof(Vec3) + sizeof(mesh::Classification)) +
        counts.cells * sizeof(mesh::Cell) +
        counts.connectivity * sizeof(Index);

    try {
        out.phase = PrismLayerPhase::Construction;
        mesh::detail::CertifiedMeshStorage storage(
            static_cast<std::size_t>(counts.nodes),
            static_cast<std::size_t>(counts.cells),
            static_cast<std::size_t>(counts.connectivity));
        const std::uint64_t sourceNodes =
            static_cast<std::uint64_t>(sourceSurface.nodeCount());
        const std::uint64_t sourceTriangles = source.triangles.size();
        const std::uint64_t sourceBoundaryEdges = source.boundaryEdges.size();

        core::parallelForWithThreads(
            0, static_cast<std::size_t>(counts.nodes),
            request.executionThreads, [&](std::size_t ordinal) {
                const std::uint64_t layer = ordinal / sourceNodes;
                const Index node = static_cast<Index>(ordinal % sourceNodes);
                mesh::Classification cls = request.volumeClassification;
                if (layer == 0) {
                    cls = sourceSurface.nodeClass(node);
                } else if (layer == request.layers) {
                    cls = request.targetClassification;
                } else if (source.boundaryNodes[
                               static_cast<std::size_t>(node)] != 0) {
                    cls = request.sideClassification;
                }
                storage.setNode(
                    ordinal,
                    sourceSurface.node(node) +
                        direction * offsets[static_cast<std::size_t>(layer)],
                    cls);
            });

        core::parallelForWithThreads(
            0, static_cast<std::size_t>(counts.prisms),
            request.executionThreads, [&](std::size_t ordinal) {
                const std::uint64_t layer = ordinal / sourceTriangles;
                const std::uint64_t triOrdinal = ordinal % sourceTriangles;
                const Index lower = static_cast<Index>(layer * sourceNodes);
                const Index upper = static_cast<Index>(
                    (layer + 1U) * sourceNodes);
                const auto& tri = source.triangles[
                    static_cast<std::size_t>(triOrdinal)];
                const std::array<Index, 6> prism{
                    lower + tri[0], lower + tri[1], lower + tri[2],
                    upper + tri[0], upper + tri[1], upper + tri[2]};
                storage.setCell(ordinal, ordinal * 6U,
                                mesh::CellType::Prism6, prism,
                                request.volumeClassification);
            });

        if (request.emitBoundaryFaces) {
            const std::uint64_t volumeConnectivity = counts.prisms * 6U;
            const std::uint64_t triangleCells = 2U * sourceTriangles;
            const Index topOffset = static_cast<Index>(
                static_cast<std::uint64_t>(request.layers) * sourceNodes);
            core::parallelForWithThreads(
                0, static_cast<std::size_t>(sourceTriangles),
                request.executionThreads, [&](std::size_t ordinal) {
                    const auto tri = source.triangles[ordinal];
                    const std::uint64_t first = 2U * ordinal;
                    const std::array<Index,3> bottom{
                        tri[0], tri[2], tri[1]};
                    const std::array<Index,3> top{
                        topOffset + tri[0], topOffset + tri[1],
                        topOffset + tri[2]};
                    storage.setCell(
                        static_cast<std::size_t>(counts.prisms + first),
                        static_cast<std::size_t>(volumeConnectivity +
                                                 first * 3U),
                        mesh::CellType::Tri3, bottom,
                        source.triangleClasses[ordinal]);
                    storage.setCell(
                        static_cast<std::size_t>(counts.prisms + first + 1U),
                        static_cast<std::size_t>(volumeConnectivity +
                                                 (first + 1U) * 3U),
                        mesh::CellType::Tri3, top,
                        request.targetClassification);
                });

            const std::uint64_t quadCellBase = counts.prisms + triangleCells;
            const std::uint64_t quadConnectivityBase =
                volumeConnectivity + triangleCells * 3U;
            core::parallelForWithThreads(
                0, static_cast<std::size_t>(counts.boundaryQuads),
                request.executionThreads, [&](std::size_t ordinal) {
                    const std::uint64_t layer = ordinal / sourceBoundaryEdges;
                    const std::uint64_t edgeOrdinal =
                        ordinal % sourceBoundaryEdges;
                    const Index lower = static_cast<Index>(layer * sourceNodes);
                    const Index upper = static_cast<Index>(
                        (layer + 1U) * sourceNodes);
                    const auto edge = source.boundaryEdges[
                        static_cast<std::size_t>(edgeOrdinal)];
                    const std::array<Index,4> face{
                        lower + edge[0], lower + edge[1],
                        upper + edge[1], upper + edge[0]};
                    storage.setCell(
                        static_cast<std::size_t>(quadCellBase + ordinal),
                        static_cast<std::size_t>(quadConnectivityBase +
                                                 ordinal * 4U),
                        mesh::CellType::Quad4, face,
                        request.sideClassification);
                });
        }

        mesh::Mesh local;
        storage.publishTo(local);

        out.metrics.constructionWorkUnits = counts.nodes + counts.cells +
            counts.connectivity + source.triangles.size() +
            source.boundaryEdges.size();
        out.phase = PrismLayerPhase::Audit;
        if (static_cast<std::uint64_t>(local.nodeCount()) != counts.nodes ||
            static_cast<std::uint64_t>(local.cellCount()) != counts.cells ||
            static_cast<std::uint64_t>(local.connectivityCount()) !=
                counts.connectivity) {
            out.reason = PrismLayerReason::AuditFailed;
            return out;
        }
        if (request.validation == PrismLayerValidationMode::FullTopologyAudit) {
            out.metrics.audit = auditPrismLayerMesh(
                local, request.emitBoundaryFaces);
        } else {
            auto& audit = out.metrics.audit;
            audit.nodes = counts.nodes;
            audit.prisms = counts.prisms;
            audit.boundaryTriangles = counts.boundaryTriangles;
            audit.boundaryQuads = counts.boundaryQuads;
            audit.minimumPrismVolume = source.minimumProjectedArea *
                request.firstLayerThickness;
            const double lastWidth = offsets.back() -
                offsets[offsets.size() - 2U];
            audit.maximumPrismVolume = source.maximumProjectedArea *
                lastWidth;
            audit.totalVolume = source.totalProjectedArea * totalThickness;
            audit.valid = std::isfinite(audit.minimumPrismVolume) &&
                audit.minimumPrismVolume > 0.0 &&
                std::isfinite(audit.maximumPrismVolume) &&
                audit.maximumPrismVolume > 0.0 &&
                std::isfinite(audit.totalVolume) && audit.totalVolume > 0.0;
        }
        if (!out.metrics.audit.complete || !out.metrics.audit.valid) {
            out.reason = out.metrics.audit.allocationFailed
                ? PrismLayerReason::AllocationFailed
                : PrismLayerReason::AuditFailed;
            return out;
        }
        const double expected = source.totalProjectedArea * totalThickness;
        if (std::abs(out.metrics.audit.totalVolume - expected) >
            256.0 * (std::numeric_limits<double>::epsilon)() *
            std::max(1.0, std::abs(expected))) {
            out.reason = PrismLayerReason::AuditFailed;
            return out;
        }
        out.phase = PrismLayerPhase::Publication;
        destination.swap(local);
        out.status = PrismLayerStatus::Published;
        out.reason = PrismLayerReason::None;
        out.destinationUnchanged = false;
        return out;
    } catch (const std::bad_alloc&) {
        out.reason = PrismLayerReason::AllocationFailed;
    } catch (...) {
        out.reason = PrismLayerReason::InternalError;
    }
    return out;
}

bool applyCompiledPrismLayerControl(
        int faceTag,
        const sizing::CompiledMeshControls& controls,
        StraightPrismLayerRequest& request) noexcept {
    const auto* layer = controls.prismLayerForFace(faceTag);
    if (!layer) return false;
    request.layers = layer->layers;
    request.firstLayerThickness = layer->firstLayerThickness;
    request.cellToCellGrowth = layer->cellToCellGrowth;
    request.maximumAllowedGrowth = controls.maximumGrowthRate();
    request.maximumTotalThickness = layer->maximumTotalThickness;
    return true;
}

}  // namespace paacs::algo
