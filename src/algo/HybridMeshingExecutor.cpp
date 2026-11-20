// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/HybridMeshingExecutor.hpp"

#include "paacs/mesh/CertifiedMeshAppender.hpp"
#include "paacs/mesh/MeshAudit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <new>
#include <numeric>
#include <set>
#include <tuple>

namespace paacs::algo {
namespace {

bool progress(const HybridExecutionRequest &request, HybridExecutionStage stage,
              std::uint64_t completed, std::uint64_t total) {
    return !request.progress || request.progress(stage, completed, total);
}
bool finite(const Vec3 &p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
bool validClass(mesh::Classification c) noexcept {
    return c.dim >= 0 && c.dim <= 3 && c.tag >= 0;
}
mesh::Classification mergeClass(mesh::Classification a, mesh::Classification b,
                                bool &compatible) noexcept {
    if (!validClass(a))
        return b;
    if (!validClass(b))
        return a;
    if (a.dim == b.dim) {
        if (a.tag != b.tag)
            compatible = false;
        return a.tag <= b.tag ? a : b;
    }
    return a.dim < b.dim ? a : b;
}
bool validateFixed(const mesh::Mesh &mesh) noexcept {
    for (Index n = 0; n < mesh.nodeCount(); ++n)
        if (!finite(mesh.node(n)) || !validClass(mesh.nodeClass(n)))
            return false;
    for (Index ci = 0; ci < mesh.cellCount(); ++ci) {
        const auto &c = mesh.cell(ci);
        const int count = mesh::nodesPerCell(c.type);
        if (count <= 0 || !validClass(c.onEntity))
            return false;
        const Index *nodes = mesh.cellNodes(ci);
        std::vector<Index> unique;
        unique.reserve(static_cast<std::size_t>(count));
        for (int k = 0; k < count; ++k) {
            if (nodes[k] < 0 || nodes[k] >= mesh.nodeCount())
                return false;
            unique.push_back(nodes[k]);
        }
        std::sort(unique.begin(), unique.end());
        if (std::adjacent_find(unique.begin(), unique.end()) != unique.end())
            return false;
    }
    return true;
}
bool isBoundary(mesh::CellType type) noexcept {
    return type == mesh::CellType::Line2 || type == mesh::CellType::Line3 ||
           type == mesh::CellType::Tri3 || type == mesh::CellType::Tri6 ||
           type == mesh::CellType::Quad4;
}
struct BoundaryKey {
    mesh::CellType type = mesh::CellType::Tri3;
    mesh::Classification classification{};
    std::array<Index, 10> nodes{};
    std::uint8_t count = 0;
    friend bool operator<(const BoundaryKey &a, const BoundaryKey &b) noexcept {
        if (a.type != b.type)
            return a.type < b.type;
        if (a.classification.dim != b.classification.dim)
            return a.classification.dim < b.classification.dim;
        if (a.classification.tag != b.classification.tag)
            return a.classification.tag < b.classification.tag;
        if (a.count != b.count)
            return a.count < b.count;
        return a.nodes < b.nodes;
    }
};
class UnionFind {
  public:
    explicit UnionFind(std::size_t n) : parent_(n), rank_(n, 0) {
        std::iota(parent_.begin(), parent_.end(), std::size_t{0});
    }
    std::size_t find(std::size_t x) {
        while (parent_[x] != x) {
            parent_[x] = parent_[parent_[x]];
            x = parent_[x];
        }
        return x;
    }
    void unite(std::size_t a, std::size_t b) {
        a = find(a);
        b = find(b);
        if (a == b)
            return;
        if (rank_[a] < rank_[b])
            std::swap(a, b);
        parent_[b] = a;
        if (rank_[a] == rank_[b])
            ++rank_[a];
    }

  private:
    std::vector<std::size_t> parent_;
    std::vector<std::uint8_t> rank_;
};
struct GeneratedBody {
    const PlannedBody *plan = nullptr;
    mesh::Mesh fixed;
    mesh::PolyhedralMesh poly;
    bool polyhedral = false;
};
bool engineMatches(const PlannedBody &plan, const HybridBodyInput &input) noexcept {
    switch (plan.selectedEngine) {
    case sizing::VolumeEngineKind::Tet4:
        return input.sourceKind == HybridBodySourceKind::PreparedFixedMesh;
    case sizing::VolumeEngineKind::MappedHex8:
        return input.sourceKind == HybridBodySourceKind::MappedHexBlock;
    case sizing::VolumeEngineKind::StructuredPolyhedral:
        return input.sourceKind == HybridBodySourceKind::StructuredPolyhedralBlock ||
               input.sourceKind == HybridBodySourceKind::PreparedPolyhedralMesh;
    case sizing::VolumeEngineKind::LayeredPrism6:
        return input.sourceKind == HybridBodySourceKind::PreparedMixedSweep;
    case sizing::VolumeEngineKind::Auto:
        return false;
    }
    return false;
}
HybridExecutionReason generateBody(const PlannedBody &plan, const HybridBodyInput &input,
                                   GeneratedBody &output) {
    output.plan = &plan;
    if (!engineMatches(plan, input))
        return HybridExecutionReason::EngineSourceMismatch;
    switch (input.sourceKind) {
    case HybridBodySourceKind::PreparedFixedMesh:
        if (!input.preparedFixedMesh || !validateFixed(*input.preparedFixedMesh))
            return HybridExecutionReason::InvalidPreparedMesh;
        output.fixed = *input.preparedFixedMesh;
        return HybridExecutionReason::None;
    case HybridBodySourceKind::MappedHexBlock: {
        const auto result = generateMappedHexBlock(input.mappedHex, output.fixed);
        return result.published() ? HybridExecutionReason::None
                                  : HybridExecutionReason::BodyGenerationFailed;
    }
    case HybridBodySourceKind::PreparedMixedSweep: {
        if (!input.sweepSourceSurface)
            return HybridExecutionReason::InvalidPreparedMesh;
        const auto result =
            generatePreparedMixedSweep(*input.sweepSourceSurface, input.sweep, output.fixed);
        return result.published() ? HybridExecutionReason::None
                                  : HybridExecutionReason::BodyGenerationFailed;
    }
    case HybridBodySourceKind::StructuredPolyhedralBlock: {
        const auto result =
            generateStructuredPolyhedralBlock(input.structuredPolyhedral, output.poly);
        output.polyhedral = result.published();
        return result.published() ? HybridExecutionReason::None
                                  : HybridExecutionReason::BodyGenerationFailed;
    }
    case HybridBodySourceKind::PreparedPolyhedralMesh:
        if (!input.preparedPolyhedralMesh ||
            !mesh::auditPolyhedralMesh(*input.preparedPolyhedralMesh).valid)
            return HybridExecutionReason::InvalidPreparedMesh;
        output.poly = *input.preparedPolyhedralMesh;
        output.polyhedral = true;
        return HybridExecutionReason::None;
    }
    return HybridExecutionReason::InternalError;
}
} // namespace

void HybridMesh::swap(HybridMesh &other) noexcept {
    fixed.swap(other.fixed);
    polyhedral.swap(other.polyhedral);
    bodies.swap(other.bodies);
    portals.swap(other.portals);
}

HybridExecutionOutcome executeHybridMeshingPlan(const HybridExecutionRequest &request,
                                                HybridMesh &destination) noexcept {
    HybridExecutionOutcome out;
    if (!request.plan) {
        out.reason = HybridExecutionReason::MissingPlan;
        return out;
    }
    if (!std::isfinite(request.portalCoordinateTolerance) ||
        request.portalCoordinateTolerance < 0.0) {
        out.reason = HybridExecutionReason::InvalidPortalBinding;
        return out;
    }
    try {
        std::map<int, const HybridBodyInput *> inputs;
        for (const auto &input : request.bodies) {
            if (input.solidTag < 0 || !inputs.emplace(input.solidTag, &input).second) {
                out.reason = HybridExecutionReason::DuplicateBody;
                return out;
            }
        }
        if (inputs.size() != request.plan->bodies().size()) {
            out.reason = HybridExecutionReason::PlanBodyMismatch;
            return out;
        }
        if (!progress(request, HybridExecutionStage::Preflight, 0, inputs.size())) {
            out.reason = HybridExecutionReason::Cancelled;
            return out;
        }
        std::vector<GeneratedBody> generated;
        generated.reserve(inputs.size());
        for (std::size_t i = 0; i < request.plan->bodies().size(); ++i) {
            const auto &planned = request.plan->bodies()[i];
            auto it = inputs.find(planned.solidTag);
            if (it == inputs.end()) {
                out.reason = HybridExecutionReason::MissingBody;
                return out;
            }
            generated.emplace_back();
            const auto reason = generateBody(planned, *it->second, generated.back());
            if (reason != HybridExecutionReason::None) {
                out.reason = reason;
                out.stage = HybridExecutionStage::BodyGeneration;
                return out;
            }
            if (!progress(request, HybridExecutionStage::BodyGeneration, i + 1, generated.size())) {
                out.reason = HybridExecutionReason::Cancelled;
                out.stage = HybridExecutionStage::BodyGeneration;
                return out;
            }
        }

        std::vector<std::size_t> fixedBodyIndices;
        std::vector<std::size_t> nodeOffsets(generated.size(), 0);
        std::size_t totalLocalNodes = 0;
        std::uint64_t fixedCells = 0, fixedConnectivity = 0;
        for (std::size_t bi = 0; bi < generated.size(); ++bi) {
            if (generated[bi].polyhedral)
                continue;
            fixedBodyIndices.push_back(bi);
            nodeOffsets[bi] = totalLocalNodes;
            totalLocalNodes += static_cast<std::size_t>(generated[bi].fixed.nodeCount());
            fixedCells += static_cast<std::uint64_t>(generated[bi].fixed.cellCount());
            fixedConnectivity +=
                static_cast<std::uint64_t>(generated[bi].fixed.connectivityCount());
        }
        if (totalLocalNodes > request.maximumOutputNodes ||
            fixedCells > request.maximumOutputCells) {
            out.reason = HybridExecutionReason::BudgetExceeded;
            return out;
        }
        UnionFind sets(totalLocalNodes);
        std::map<int, std::size_t> bodyIndex;
        for (std::size_t i = 0; i < generated.size(); ++i)
            bodyIndex.emplace(generated[i].plan->solidTag, i);
        std::map<std::tuple<int, int, int>, std::uint64_t> bindingCounts;
        for (const auto &binding : request.portalNodeBindings) {
            auto first = bodyIndex.find(binding.firstSolidTag),
                 second = bodyIndex.find(binding.secondSolidTag);
            if (binding.faceTag < 0 || first == bodyIndex.end() || second == bodyIndex.end()) {
                out.reason = HybridExecutionReason::InvalidPortalBinding;
                return out;
            }
            if (generated[first->second].polyhedral || generated[second->second].polyhedral) {
                out.reason = HybridExecutionReason::UnsupportedPolyhedralPortal;
                return out;
            }
            const auto &a = generated[first->second].fixed;
            const auto &b = generated[second->second].fixed;
            if (binding.firstLocalNode < 0 || binding.firstLocalNode >= a.nodeCount() ||
                binding.secondLocalNode < 0 || binding.secondLocalNode >= b.nodeCount()) {
                out.reason = HybridExecutionReason::InvalidPortalBinding;
                return out;
            }
            const Vec3 delta = a.node(binding.firstLocalNode) - b.node(binding.secondLocalNode);
            const double scale = std::max({1.0, a.node(binding.firstLocalNode).norm(),
                                           b.node(binding.secondLocalNode).norm()});
            if (delta.norm() > request.portalCoordinateTolerance * scale) {
                out.reason = HybridExecutionReason::PortalCoordinateMismatch;
                return out;
            }
            sets.unite(
                nodeOffsets[first->second] + static_cast<std::size_t>(binding.firstLocalNode),
                nodeOffsets[second->second] + static_cast<std::size_t>(binding.secondLocalNode));
            ++bindingCounts[{binding.faceTag, binding.firstSolidTag, binding.secondSolidTag}];
        }
        out.metrics.portalBindings = request.portalNodeBindings.size();
        if (!progress(request, HybridExecutionStage::PortalReconciliation,
                      request.portalNodeBindings.size(), request.portalNodeBindings.size())) {
            out.reason = HybridExecutionReason::Cancelled;
            out.stage = HybridExecutionStage::PortalReconciliation;
            return out;
        }

        std::vector<mesh::Classification> groupClass(totalLocalNodes, {-1, -1});
        std::vector<Vec3> groupPoint(totalLocalNodes);
        std::vector<std::uint8_t> seen(totalLocalNodes, 0);
        bool classCompatible = true;
        for (const std::size_t bi : fixedBodyIndices) {
            const auto &mesh = generated[bi].fixed;
            for (Index local = 0; local < mesh.nodeCount(); ++local) {
                const std::size_t id = nodeOffsets[bi] + static_cast<std::size_t>(local);
                const std::size_t root = sets.find(id);
                if (!seen[root]) {
                    seen[root] = 1;
                    groupPoint[root] = mesh.node(local);
                    groupClass[root] = mesh.nodeClass(local);
                } else
                    groupClass[root] =
                        mergeClass(groupClass[root], mesh.nodeClass(local), classCompatible);
            }
        }
        if (!classCompatible) {
            out.reason = HybridExecutionReason::PortalClassificationMismatch;
            return out;
        }
        std::vector<Index> rootGlobal(totalLocalNodes, kInvalidIndex),
            localGlobal(totalLocalNodes, kInvalidIndex);
        mesh::Mesh assembled;
        mesh::detail::CertifiedMeshAppender appender(assembled, totalLocalNodes,
                                                     static_cast<std::size_t>(fixedCells),
                                                     static_cast<std::size_t>(fixedConnectivity));
        for (const std::size_t bi : fixedBodyIndices) {
            const auto &mesh = generated[bi].fixed;
            for (Index local = 0; local < mesh.nodeCount(); ++local) {
                const std::size_t id = nodeOffsets[bi] + static_cast<std::size_t>(local);
                const std::size_t root = sets.find(id);
                if (rootGlobal[root] == kInvalidIndex)
                    rootGlobal[root] = appender.addNode(groupPoint[root], groupClass[root]);
                localGlobal[id] = rootGlobal[root];
            }
        }
        out.metrics.mergedPortalNodes =
            totalLocalNodes - static_cast<std::size_t>(assembled.nodeCount());
        std::set<BoundaryKey> boundaries;
        std::vector<HybridBodyRecord> fixedRecords;
        fixedRecords.reserve(fixedBodyIndices.size());
        for (const std::size_t bi : fixedBodyIndices) {
            const auto &body = generated[bi];
            HybridBodyRecord record;
            record.solidTag = body.plan->solidTag;
            record.engine = body.plan->selectedEngine;
            record.firstCell = assembled.cellCount();
            record.generatedNodes = body.fixed.nodeCount();
            record.generatedCells = body.fixed.cellCount();
            for (Index ci = 0; ci < body.fixed.cellCount(); ++ci) {
                const auto &cell = body.fixed.cell(ci);
                const int count = mesh::nodesPerCell(cell.type);
                std::array<Index, 10> conn{};
                for (int k = 0; k < count; ++k)
                    conn[static_cast<std::size_t>(k)] =
                        localGlobal[nodeOffsets[bi] +
                                    static_cast<std::size_t>(body.fixed.cellNodes(ci)[k])];
                if (request.suppressDuplicateBoundaryFaces && isBoundary(cell.type)) {
                    BoundaryKey keyRow;
                    keyRow.type = cell.type;
                    keyRow.classification = cell.onEntity;
                    keyRow.count = static_cast<std::uint8_t>(count);
                    for (int k = 0; k < count; ++k)
                        keyRow.nodes[static_cast<std::size_t>(k)] =
                            conn[static_cast<std::size_t>(k)];
                    std::sort(keyRow.nodes.begin(), keyRow.nodes.begin() + count);
                    if (!boundaries.insert(keyRow).second) {
                        ++out.metrics.duplicateBoundaryFacesSuppressed;
                        continue;
                    }
                }
                switch (count) {
                case 2: {
                    std::array<Index, 2> a{conn[0], conn[1]};
                    (void)appender.addCell(cell.type, a, cell.onEntity);
                    break;
                }
                case 3: {
                    std::array<Index, 3> a{conn[0], conn[1], conn[2]};
                    (void)appender.addCell(cell.type, a, cell.onEntity);
                    break;
                }
                case 4: {
                    std::array<Index, 4> a{conn[0], conn[1], conn[2], conn[3]};
                    (void)appender.addCell(cell.type, a, cell.onEntity);
                    break;
                }
                case 5: {
                    std::array<Index, 5> a{conn[0], conn[1], conn[2], conn[3], conn[4]};
                    (void)appender.addCell(cell.type, a, cell.onEntity);
                    break;
                }
                case 6: {
                    std::array<Index, 6> a{conn[0], conn[1], conn[2], conn[3], conn[4], conn[5]};
                    (void)appender.addCell(cell.type, a, cell.onEntity);
                    break;
                }
                case 8: {
                    std::array<Index, 8> a{conn[0], conn[1], conn[2], conn[3],
                                           conn[4], conn[5], conn[6], conn[7]};
                    (void)appender.addCell(cell.type, a, cell.onEntity);
                    break;
                }
                case 10: {
                    std::array<Index, 10> a = conn;
                    (void)appender.addCell(cell.type, a, cell.onEntity);
                    break;
                }
                default:
                    out.reason = HybridExecutionReason::InvalidPreparedMesh;
                    return out;
                }
            }
            record.cellCount = assembled.cellCount() - record.firstCell;
            record.strictHash = mesh::meshHash(body.fixed);
            fixedRecords.push_back(record);
            out.metrics.fixedBodies++;
            out.metrics.bodiesGenerated++;
            out.metrics.outputCells += record.cellCount;
            out.metrics.outputNodes = assembled.nodeCount();
            out.stage = HybridExecutionStage::Assembly;
        }

        HybridMesh local;
        local.fixed.swap(assembled);
        local.bodies = std::move(fixedRecords);
        for (std::size_t bi = 0; bi < generated.size(); ++bi) {
            if (!generated[bi].polyhedral)
                continue;
            const auto &source = generated[bi].poly;
            const Index nodeOffset = local.polyhedral.nodeCount();
            const Index faceOffset = local.polyhedral.faceCount();
            const Index cellOffset = local.polyhedral.cellCount();
            for (Index n = 0; n < source.nodeCount(); ++n)
                (void)local.polyhedral.addNode(source.node(n), source.nodeClass(n));
            for (Index f = 0; f < source.faceCount(); ++f) {
                std::vector<Index> nodes;
                for (Index n : source.faceNodes(f))
                    nodes.push_back(n + nodeOffset);
                const auto &face = source.face(f);
                (void)local.polyhedral.addFace(
                    nodes, face.owner + cellOffset,
                    face.neighbour == kInvalidIndex ? kInvalidIndex : face.neighbour + cellOffset,
                    face.onEntity);
            }
            for (Index c = 0; c < source.cellCount(); ++c) {
                std::vector<mesh::OrientedPolyFace> faces;
                for (const auto &ref : source.cellFaces(c))
                    faces.push_back({ref.face + faceOffset, ref.sense});
                (void)local.polyhedral.addCell(faces, source.cell(c).onEntity);
            }
            local.bodies.push_back({generated[bi].plan->solidTag,
                                    generated[bi].plan->selectedEngine, true, cellOffset,
                                    source.cellCount(), source.nodeCount(), source.cellCount(),
                                    mesh::polyhedralStrictHash(source)});
            ++out.metrics.polyhedralBodies;
            ++out.metrics.bodiesGenerated;
            out.metrics.outputPolyhedralCells += source.cellCount();
        }
        for (const auto &portal : request.plan->portals()) {
            std::uint64_t count = 0;
            for (const auto &[key, value] : bindingCounts)
                if (std::get<0>(key) == portal.faceTag)
                    count += value;
            local.portals.push_back({portal.faceTag, portal.firstSolidTag, portal.secondSolidTag,
                                     portal.strategy, count, portal.identityHash});
        }
        out.stage = HybridExecutionStage::Audit;
        if (local.fixed.cellCount() > 0 && !validateFixed(local.fixed)) {
            out.reason = HybridExecutionReason::AuditFailed;
            return out;
        }
        if (local.polyhedral.cellCount() > 0 &&
            !mesh::auditPolyhedralMesh(local.polyhedral).valid) {
            out.reason = HybridExecutionReason::AuditFailed;
            return out;
        }
        if (!progress(request, HybridExecutionStage::Audit, 1, 1)) {
            out.reason = HybridExecutionReason::Cancelled;
            return out;
        }
        local.swap(destination);
        out.published = true;
        out.reason = HybridExecutionReason::None;
        out.stage = HybridExecutionStage::Publication;
        out.destinationUnchanged = false;
        return out;
    } catch (const std::bad_alloc &) {
        out.reason = HybridExecutionReason::AllocationFailed;
    } catch (...) {
        out.reason = HybridExecutionReason::InternalError;
    }
    return out;
}

} // namespace paacs::algo
