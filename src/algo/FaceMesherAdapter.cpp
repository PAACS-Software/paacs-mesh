// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#if defined(_MSC_VER)
#pragma warning(push)
// TransactionAccess's production-only assert parameter is intentionally used
// only under PAACS_TOPO_TESTING. This TU is the first production library TU to
// include that BND header; keep the additive adapter /W4-clean without editing
// the frozen topology file.
#pragma warning(disable : 4100)
#endif
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include "paacs/algo/FaceMesherAdapter.hpp"
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_set>

namespace paacs::algo {
namespace {

bool sameClass(const mesh::Classification& a,
               const mesh::Classification& b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}

bool samePoint(const Vec3& a, const Vec3& b) noexcept {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

using LocalEdge = std::array<std::size_t, 2>;
LocalEdge localEdge(std::size_t a, std::size_t b) {
    return a < b ? LocalEdge{a, b} : LocalEdge{b, a};
}

const bnd::ConstraintSegment* segmentByEntity(
    const SyntheticFaceInput& input, bnd::ConstraintEntityId id) {
    for (const auto& segment : input.constraints)
        if (segment.id == id) return &segment;
    return nullptr;
}

FaceAdapterOutcome adapterRejected(FaceAdapterReason reason,
                                   bnd::SteinerReject reject =
                                       bnd::SteinerReject::None) {
    FaceAdapterOutcome out;
    out.status = FaceAdapterStatus::Rejected;
    out.reason = reason;
    out.reject = reject;
    return out;
}

FaceAdapterOutcome adapterFailed(FaceAdapterReason reason,
                                 bnd::SteinerFail fail =
                                     bnd::SteinerFail::None) {
    FaceAdapterOutcome out;
    out.status = FaceAdapterStatus::Failed;
    out.reason = reason;
    out.fail = fail;
    return out;
}

FaceMeshAudit auditLocalResult(
    const SyntheticFaceInput& input, const mesh::Mesh& global,
    const bnd::ConstraintSet& constraints,
    const bnd::SteinerSegmentsReport& recovery,
    const std::set<bnd::ConstraintEntityId>& boundaryEntities,
    const FaceMeshResult& result) {
    FaceMeshAudit audit;
    if (!result.localWorkspace) return audit;
    const auto& ws = *result.localWorkspace;

    audit.workspaceReady =
        topo::workspaceTopologyReady(ws.nodes(), ws.cells(), ws.edges(),
                                     ws.faces()) &&
        ws.adjacencyReady();
    audit.tri3Manifold =
        bnd::auditTri3Manifold(ws) == bnd::Tri3AuditResult::Ok;

    if (result.localNodes.size() != result.localParams.size() ||
        result.localNodes.size() != result.localHandles.size())
        return audit;

    std::unordered_map<topo::NodeHandle, std::size_t> indexOf;
    bool exportedNodesMatch = true;
    for (std::size_t i = 0; i < result.localHandles.size(); ++i) {
        const auto& h = result.localHandles[i];
        if (!ws.nodes().isCurrent(h) || !indexOf.emplace(h, i).second) {
            exportedNodesMatch = false;
            continue;
        }
        const auto node = ws.nodes().get(h);
        const auto* param = constraints.paramOf(h);
        if (!param || param->chart != result.chart ||
            !samePoint(node.coord, result.localNodes[i].coord) ||
            !sameClass(node.cls, result.localNodes[i].cls) ||
            param->p.u != result.localParams[i].u ||
            param->p.v != result.localParams[i].v)
            exportedNodesMatch = false;
    }

    std::set<Index> globalIds;
    bool mappingsOk = exportedNodesMatch;
    for (const auto& [local, globalIndex] : result.existingGlobalIndex) {
        if (local >= result.localNodes.size() || globalIndex < 0 ||
            globalIndex >= global.nodeCount() ||
            !globalIds.insert(globalIndex).second) {
            mappingsOk = false;
            continue;
        }
        if (!samePoint(result.localNodes[local].coord,
                       global.node(globalIndex)) ||
            !sameClass(result.localNodes[local].cls,
                       global.nodeClass(globalIndex)))
            mappingsOk = false;
    }
    audit.mappingsInjective = mappingsOk;

    bool cellsMatch = result.localCells.size() ==
                      ws.cells().activeBySequence().size();
    predicates::Sign orientation = predicates::Sign::Zero;
    std::set<std::array<std::size_t, 3>> triangleKeys;
    std::map<LocalEdge, int> edgeUse;
    std::map<LocalEdge, std::vector<std::size_t>> edgeOwners;
    std::size_t exportedCell = 0;
    for (const auto& ch : ws.cells().activeBySequence()) {
        const auto cell = ws.cells().get(ch);
        if (cell.type != topo::ElementType::Tri3 ||
            exportedCell >= result.localCells.size()) {
            cellsMatch = false;
            break;
        }
        const auto& exported = result.localCells[exportedCell++];
        if (!sameClass(cell.cls, result.faceClass) ||
            !sameClass(exported.classification, result.faceClass))
            cellsMatch = false;
        std::array<std::size_t, 3> mapped{};
        for (std::size_t i = 0; i < 3; ++i) {
            const auto it = indexOf.find(cell.nodes[i]);
            if (it == indexOf.end()) {
                cellsMatch = false;
                mapped[i] = 0;
            } else {
                mapped[i] = it->second;
            }
        }
        if (mapped != exported.nodes) cellsMatch = false;
        auto key = mapped;
        std::sort(key.begin(), key.end());
        if (!triangleKeys.insert(key).second) cellsMatch = false;
        for (std::size_t i = 0; i < 3; ++i) {
            const auto edge = localEdge(mapped[i], mapped[(i + 1) % 3]);
            ++edgeUse[edge];
            edgeOwners[edge].push_back(exportedCell - 1);
        }
        if (mapped[0] >= result.localParams.size() ||
            mapped[1] >= result.localParams.size() ||
            mapped[2] >= result.localParams.size()) {
            cellsMatch = false;
            continue;
        }
        const auto sign = predicates::orient2d(
            result.localParams[mapped[0]], result.localParams[mapped[1]],
            result.localParams[mapped[2]]);
        if (sign == predicates::Sign::Zero)
            cellsMatch = false;
        else if (orientation == predicates::Sign::Zero)
            orientation = sign;
        else if (orientation != sign)
            cellsMatch = false;
    }
    audit.orientationConsistent =
        cellsMatch && orientation != predicates::Sign::Zero;

    std::set<LocalEdge> expectedBoundary;
    bool loopOk = result.outerLoop.size() == input.outerLoop.size() &&
                  result.outerLoop.size() >= 3;
    std::set<std::size_t> loopVertices;
    for (std::size_t i = 0; i < result.outerLoop.size(); ++i) {
        const std::size_t a = result.outerLoop[i];
        const std::size_t b =
            result.outerLoop[(i + 1) % result.outerLoop.size()];
        if (a >= result.localNodes.size() || b >= result.localNodes.size() ||
            a == b || !loopVertices.insert(a).second)
            loopOk = false;
        expectedBoundary.insert(localEdge(a, b));
    }
    std::set<LocalEdge> actualBoundary;
    for (const auto& [edge, uses] : edgeUse) {
        if (uses == 1) actualBoundary.insert(edge);
        else if (uses != 2) loopOk = false;
    }
    const std::int64_t euler =
        static_cast<std::int64_t>(result.localNodes.size()) -
        static_cast<std::int64_t>(edgeUse.size()) +
        static_cast<std::int64_t>(result.localCells.size());
    std::vector<std::vector<std::size_t>> cellNeighbors(
        result.localCells.size());
    for (const auto& [edge, owners] : edgeOwners) {
        (void)edge;
        if (owners.size() == 2) {
            cellNeighbors[owners[0]].push_back(owners[1]);
            cellNeighbors[owners[1]].push_back(owners[0]);
        }
    }
    std::set<std::size_t> connectedCells;
    std::vector<std::size_t> stack;
    if (!result.localCells.empty()) stack.push_back(0);
    while (!stack.empty()) {
        const auto cell = stack.back();
        stack.pop_back();
        if (!connectedCells.insert(cell).second) continue;
        for (const auto neighbor : cellNeighbors[cell])
            stack.push_back(neighbor);
    }
    audit.oneOuterLoop =
        loopOk && actualBoundary == expectedBoundary && euler == 1 &&
        connectedCells.size() == result.localCells.size();

    bool classificationOk = cellsMatch;
    bool recovered = recovery.processed.size() == input.constraints.size() &&
                     recovery.unattempted.empty() &&
                     result.recoveredConstraints.size() ==
                         input.constraints.size();
    bool boundaryFree = true;
    std::map<bnd::ConstraintEntityId, const bnd::SteinerEntityResult*> rows;
    for (const auto& row : recovery.processed) rows[row.entity] = &row;
    for (const auto& segment : input.constraints) {
        const auto row = rows.find(segment.id);
        if (row == rows.end()) {
            recovered = false;
            continue;
        }
        const auto& rr = row->second->result;
        if (!(rr.status == bnd::SteinerStatus::AlreadyPresent ||
              rr.status == bnd::SteinerStatus::Recovered) ||
            rr.recoveredChain.size() < 2)
            recovered = false;
        if (boundaryEntities.count(segment.id) != 0) {
            if (rr.status != bnd::SteinerStatus::AlreadyPresent ||
                rr.recoveredChain.size() != 2 ||
                !rr.createdParams.empty())
                boundaryFree = false;
        }
        for (std::size_t i = 0; i + 1 < rr.recoveredChain.size(); ++i) {
            const auto edge = ws.edges().find(
                ws.nodes(), rr.recoveredChain[i], rr.recoveredChain[i + 1]);
            if (!edge) {
                recovered = false;
                continue;
            }
            const auto data = ws.edges().get(*edge);
            if (!data.requiredConstraint ||
                !sameClass(data.cls, segment.cls)) {
                recovered = false;
                classificationOk = false;
            }
        }
    }
    for (const std::size_t local : result.outerLoop)
        if (result.existingGlobalIndex.find(local) ==
            result.existingGlobalIndex.end())
            boundaryFree = false;

    audit.classificationsPreserved = classificationOk;
    audit.constraintsRecovered = recovered;
    audit.boundarySteinerFree = boundaryFree;
    return audit;
}

}  // namespace

FaceAdapterOutcome buildSyntheticFaceMesh(
    const SyntheticFaceInput& input, const mesh::Mesh& global,
    const bnd::GeometryProvider& provider,
    const bnd::GeometryProviderStamp& expectedProvider,
    const bnd::SteinerPolicy& policy) {
    if (input.periodic) {
        FaceAdapterOutcome out;
        out.status = FaceAdapterStatus::Refused;
        out.reason = FaceAdapterReason::UnsupportedPeriodicChart;
        out.refuse = bnd::SteinerRefuse::UnsupportedChart;
        return out;
    }
    try {
        if (input.faceClass.dim != 2 || input.faceClass.tag < 0 ||
            input.nodes.size() < 3 || input.triangles.empty() ||
            input.outerLoop.size() < 3 || input.constraints.empty())
            return adapterRejected(FaceAdapterReason::InvalidInput);

        std::set<std::size_t> loopNodes;
        std::set<Index> mappedGlobals;
        std::set<bnd::ConstraintVertexId> vertexIds;
        for (const std::size_t local : input.outerLoop)
            if (local >= input.nodes.size() ||
                !loopNodes.insert(local).second ||
                !input.nodes[local].existingGlobalIndex)
                return adapterRejected(FaceAdapterReason::InvalidInput);
        for (const auto& node : input.nodes) {
            if (!std::isfinite(node.uv.u) || !std::isfinite(node.uv.v))
                return adapterRejected(FaceAdapterReason::InvalidInput);
            if (node.existingGlobalIndex) {
                const Index globalIndex = *node.existingGlobalIndex;
                if (globalIndex < 0 || globalIndex >= global.nodeCount() ||
                    !mappedGlobals.insert(globalIndex).second ||
                    !sameClass(node.classification,
                               global.nodeClass(globalIndex)))
                    return adapterRejected(FaceAdapterReason::InvalidInput);
            }
            if (node.constraintVertex &&
                !vertexIds.insert(*node.constraintVertex).second)
                return adapterRejected(FaceAdapterReason::InvalidInput);
        }

        bnd::GeometryProviderStamp providerStamp;
        try {
            providerStamp = provider.stamp();
        } catch (const std::bad_alloc&) {
            return adapterFailed(FaceAdapterReason::RecoveryStopped,
                                 bnd::SteinerFail::CandidateAllocFailed);
        } catch (...) {
            return adapterFailed(FaceAdapterReason::RecoveryStopped,
                                 bnd::SteinerFail::GeometryProviderError);
        }
        if (!(providerStamp == expectedProvider))
            return adapterRejected(FaceAdapterReason::RecoveryStopped,
                                   bnd::SteinerReject::
                                       StaleGeometryProvider);

        topo::NodePool nodes;
        std::vector<topo::NodeHandle> handles;
        handles.reserve(input.nodes.size());
        for (const auto& node : input.nodes) {
            Vec3 physical{};
            if (node.existingGlobalIndex) {
                physical = global.node(*node.existingGlobalIndex);
            } else {
                try {
                    physical = provider.pointAt(input.chart, node.uv);
                } catch (const bnd::GeometryProviderDomainError&) {
                    return adapterRejected(
                        FaceAdapterReason::RecoveryStopped,
                        bnd::SteinerReject::ProviderDomainRejected);
                } catch (const std::bad_alloc&) {
                    return adapterFailed(
                        FaceAdapterReason::RecoveryStopped,
                        bnd::SteinerFail::CandidateAllocFailed);
                } catch (...) {
                    return adapterFailed(
                        FaceAdapterReason::RecoveryStopped,
                        bnd::SteinerFail::GeometryProviderError);
                }
                if (!finite(physical))
                    return adapterRejected(
                        FaceAdapterReason::RecoveryStopped,
                        bnd::SteinerReject::NonfiniteMappedPosition);
            }
            handles.push_back(nodes.add(physical, node.classification));
        }

        topo::CellPool cells;
        for (const auto& triangle : input.triangles) {
            const auto& n = triangle.nodes;
            if (n[0] >= handles.size() || n[1] >= handles.size() ||
                n[2] >= handles.size() || n[0] == n[1] ||
                n[1] == n[2] || n[2] == n[0])
                return adapterRejected(FaceAdapterReason::InvalidInput);
            cells.add(topo::ElementType::Tri3,
                      std::vector<topo::NodeHandle>{
                          handles[n[0]], handles[n[1]], handles[n[2]]},
                      input.faceClass, nodes);
        }
        topo::EdgePool edges;
        for (const auto& ch : cells.activeBySequence())
            edges.insertCellEdges(nodes, cells, ch);

        std::set<bnd::ConstraintEntityId> boundaryEntities;
        for (std::size_t i = 0; i < input.outerLoop.size(); ++i) {
            const std::size_t la = input.outerLoop[i];
            const std::size_t lb =
                input.outerLoop[(i + 1) % input.outerLoop.size()];
            if (!input.nodes[la].constraintVertex ||
                !input.nodes[lb].constraintVertex)
                return adapterRejected(FaceAdapterReason::InvalidInput);
            const auto va = *input.nodes[la].constraintVertex;
            const auto vb = *input.nodes[lb].constraintVertex;
            const bnd::ConstraintSegment* boundary = nullptr;
            for (const auto& segment : input.constraints)
                if ((segment.a == va && segment.b == vb) ||
                    (segment.a == vb && segment.b == va)) {
                    if (boundary)
                        return adapterRejected(
                            FaceAdapterReason::InvalidInput);
                    boundary = &segment;
                }
            if (!boundary || boundary->chart != input.chart ||
                boundary->cls.dim != 1 || boundary->cls.tag < 0)
                return adapterRejected(FaceAdapterReason::InvalidInput);
            const auto edge =
                edges.find(nodes, handles[la], handles[lb]);
            if (!edge)
                return adapterRejected(FaceAdapterReason::InvalidInput);
            edges.setRequiredConstraint(*edge, true);
            edges.setClassification(*edge, boundary->cls);
            boundaryEntities.insert(boundary->id);
        }

        auto workspace = std::make_unique<topo::TopologyWorkspace>(
            std::move(nodes), std::move(cells), std::move(edges),
            topo::FacePool{});
        bnd::ConstraintSet constraintSet;
        for (std::size_t i = 0; i < input.nodes.size(); ++i) {
            constraintSet.params[handles[i]] =
                bnd::ParamPoint{input.nodes[i].uv, input.chart};
            if (input.nodes[i].constraintVertex)
                constraintSet.vertices[*input.nodes[i].constraintVertex] =
                    handles[i];
        }
        constraintSet.segments = input.constraints;
        constraintSet.stamp = bnd::OwnershipStamp::of(*workspace);

        auto recovery = bnd::recoverSteinerSegments(
            *workspace, constraintSet, provider, expectedProvider, policy);
        bool allRecovered =
            recovery.processed.size() == input.constraints.size() &&
            recovery.unattempted.empty();
        const bnd::SteinerResult* stopped = nullptr;
        for (const auto& row : recovery.processed)
            if (!(row.result.status == bnd::SteinerStatus::AlreadyPresent ||
                  row.result.status == bnd::SteinerStatus::Recovered)) {
                allRecovered = false;
                if (!stopped) stopped = &row.result;
            }
        if (!allRecovered) {
            FaceAdapterOutcome out;
            out.reason = FaceAdapterReason::RecoveryStopped;
            out.recovery = std::move(recovery);
            if (!stopped) {
                out.status = FaceAdapterStatus::Failed;
                out.fail = bnd::SteinerFail::CandidateAllocFailed;
            } else {
                out.reject = stopped->reject;
                out.refuse = stopped->refuse;
                out.fail = stopped->fail;
                if (stopped->status == bnd::SteinerStatus::Rejected)
                    out.status = FaceAdapterStatus::Rejected;
                else if (stopped->status == bnd::SteinerStatus::Refused)
                    out.status = FaceAdapterStatus::Refused;
                else
                    out.status = FaceAdapterStatus::Failed;
            }
            return out;
        }

        // Every boundary entity was pre-marked and must remain a two-node
        // AlreadyPresent chain: any created boundary node is a hard scope
        // violation, even if the topology would otherwise audit.
        for (const auto& row : recovery.processed)
            if (boundaryEntities.count(row.entity) != 0 &&
                (row.result.status != bnd::SteinerStatus::AlreadyPresent ||
                 row.result.recoveredChain.size() != 2 ||
                 !row.result.createdParams.empty())) {
                FaceAdapterOutcome out =
                    adapterRejected(
                        FaceAdapterReason::BoundarySteinerProhibited);
                out.recovery = std::move(recovery);
                return out;
            }

        workspace->rebuildAdjacency();
        auto result = std::make_unique<FaceMeshResult>();
        result->chart = input.chart;
        result->faceClass = input.faceClass;
        result->localWorkspace = std::move(workspace);

        std::unordered_map<topo::NodeHandle, Index> originalGlobal;
        for (std::size_t i = 0; i < input.nodes.size(); ++i)
            if (input.nodes[i].existingGlobalIndex)
                originalGlobal[handles[i]] =
                    *input.nodes[i].existingGlobalIndex;

        std::unordered_map<topo::NodeHandle, std::size_t> indexOf;
        for (const auto& h :
             result->localWorkspace->nodes().activeBySequence()) {
            const std::size_t local = result->localNodes.size();
            indexOf[h] = local;
            result->localHandles.push_back(h);
            result->localNodes.push_back(
                result->localWorkspace->nodes().get(h));
            const auto* parameter = constraintSet.paramOf(h);
            if (!parameter)
                return adapterFailed(FaceAdapterReason::LocalAuditFailed);
            result->localParams.push_back(parameter->p);
            const auto globalIt = originalGlobal.find(h);
            if (globalIt != originalGlobal.end())
                result->existingGlobalIndex[local] = globalIt->second;
        }
        for (const auto& ch :
             result->localWorkspace->cells().activeBySequence()) {
            const auto cell = result->localWorkspace->cells().get(ch);
            FaceMeshCell exported;
            exported.classification = cell.cls;
            for (std::size_t i = 0; i < 3; ++i) {
                const auto it = indexOf.find(cell.nodes[i]);
                if (it == indexOf.end())
                    return adapterFailed(
                        FaceAdapterReason::LocalAuditFailed);
                exported.nodes[i] = it->second;
            }
            result->localCells.push_back(exported);
        }
        for (const std::size_t inputLocal : input.outerLoop) {
            const auto it = indexOf.find(handles[inputLocal]);
            if (it == indexOf.end())
                return adapterFailed(FaceAdapterReason::LocalAuditFailed);
            result->outerLoop.push_back(it->second);
        }
        for (const auto& row : recovery.processed) {
            const auto* segment = segmentByEntity(input, row.entity);
            if (!segment)
                return adapterFailed(FaceAdapterReason::LocalAuditFailed);
            FaceRecoveredConstraint exported;
            exported.entity = row.entity;
            exported.classification = segment->cls;
            exported.status = row.result.status;
            for (const auto& h : row.result.recoveredChain) {
                const auto it = indexOf.find(h);
                if (it == indexOf.end())
                    return adapterFailed(
                        FaceAdapterReason::LocalAuditFailed);
                exported.chain.push_back(it->second);
            }
            result->recoveredConstraints.push_back(std::move(exported));
        }

        result->audit = auditLocalResult(
            input, global, constraintSet, recovery, boundaryEntities, *result);
        if (!result->audit.passed()) {
            FaceAdapterOutcome out =
                adapterFailed(FaceAdapterReason::LocalAuditFailed);
            out.recovery = std::move(recovery);
            return out;
        }
        FaceAdapterOutcome out;
        out.status = FaceAdapterStatus::Ready;
        out.reason = FaceAdapterReason::None;
        out.recovery = std::move(recovery);
        out.result = std::move(result);
        return out;
    } catch (const std::bad_alloc&) {
        return adapterFailed(FaceAdapterReason::AllocationFailed,
                             bnd::SteinerFail::CandidateAllocFailed);
    } catch (...) {
        return adapterFailed(FaceAdapterReason::InternalError);
    }
}

FacePublication publishFaceMeshResult(mesh::Mesh& global,
                                      const FaceMeshResult& result) noexcept {
    try {
        if (!result.audit.passed() || !result.localWorkspace ||
            result.localNodes.size() != result.localParams.size() ||
            result.localNodes.size() != result.localHandles.size())
            return {FacePublishStatus::Rejected,
                    FacePublishReason::InvalidResult, {}};

        const auto& ws = *result.localWorkspace;
        if (!topo::workspaceTopologyReady(
                ws.nodes(), ws.cells(), ws.edges(), ws.faces()) ||
            result.localCells.size() !=
                ws.cells().activeBySequence().size())
            return {FacePublishStatus::Rejected,
                    FacePublishReason::InvalidResult, {}};

        std::unordered_map<topo::NodeHandle, std::size_t> indexOf;
        std::set<Index> reused;
        for (std::size_t i = 0; i < result.localNodes.size(); ++i) {
            const auto& h = result.localHandles[i];
            if (!ws.nodes().isCurrent(h) || !indexOf.emplace(h, i).second)
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};
            const auto live = ws.nodes().get(h);
            if (!samePoint(live.coord, result.localNodes[i].coord) ||
                !sameClass(live.cls, result.localNodes[i].cls))
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};
            const auto existing = result.existingGlobalIndex.find(i);
            if (existing == result.existingGlobalIndex.end()) continue;
            const Index globalIndex = existing->second;
            if (globalIndex < 0 || globalIndex >= global.nodeCount() ||
                !reused.insert(globalIndex).second ||
                !samePoint(result.localNodes[i].coord,
                           global.node(globalIndex)) ||
                !sameClass(result.localNodes[i].cls,
                           global.nodeClass(globalIndex)))
                return {FacePublishStatus::Rejected,
                        FacePublishReason::ExistingIdentityMismatch, {}};
        }
        for (const auto& [local, globalIndex] :
             result.existingGlobalIndex) {
            (void)globalIndex;
            if (local >= result.localNodes.size())
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};
        }
        std::map<LocalEdge, int> edgeUse;
        predicates::Sign orientation = predicates::Sign::Zero;
        std::size_t cellIndex = 0;
        for (const auto& ch : ws.cells().activeBySequence()) {
            if (cellIndex >= result.localCells.size())
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};
            const auto live = ws.cells().get(ch);
            const auto& exported = result.localCells[cellIndex++];
            if (live.type != topo::ElementType::Tri3 ||
                !sameClass(live.cls, result.faceClass) ||
                !sameClass(exported.classification, result.faceClass))
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};
            std::array<std::size_t, 3> mapped{};
            for (std::size_t i = 0; i < 3; ++i) {
                const auto it = indexOf.find(live.nodes[i]);
                if (it == indexOf.end())
                    return {FacePublishStatus::Rejected,
                            FacePublishReason::InvalidResult, {}};
                mapped[i] = it->second;
            }
            if (mapped != exported.nodes ||
                mapped[0] >= result.localParams.size() ||
                mapped[1] >= result.localParams.size() ||
                mapped[2] >= result.localParams.size())
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};
            const auto sign = predicates::orient2d(
                result.localParams[mapped[0]], result.localParams[mapped[1]],
                result.localParams[mapped[2]]);
            if (sign == predicates::Sign::Zero ||
                (orientation != predicates::Sign::Zero &&
                 sign != orientation))
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};
            orientation = sign;
            for (std::size_t i = 0; i < 3; ++i)
                ++edgeUse[
                    localEdge(mapped[i], mapped[(i + 1) % 3])];
        }

        std::set<LocalEdge> expectedBoundary;
        for (std::size_t i = 0; i < result.outerLoop.size(); ++i) {
            const std::size_t a = result.outerLoop[i];
            const std::size_t b =
                result.outerLoop[(i + 1) % result.outerLoop.size()];
            if (a >= result.localNodes.size() ||
                b >= result.localNodes.size() || a == b)
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};
            expectedBoundary.insert(localEdge(a, b));
        }
        std::set<LocalEdge> actualBoundary;
        for (const auto& [edge, use] : edgeUse) {
            if (use == 1) actualBoundary.insert(edge);
            else if (use != 2)
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};
        }
        const std::int64_t euler =
            static_cast<std::int64_t>(result.localNodes.size()) -
            static_cast<std::int64_t>(edgeUse.size()) +
            static_cast<std::int64_t>(result.localCells.size());
        if (result.outerLoop.size() < 3 ||
            expectedBoundary != actualBoundary || euler != 1)
            return {FacePublishStatus::Rejected,
                    FacePublishReason::InvalidResult, {}};

        std::set<LocalEdge> recoveredEdges;
        for (const auto& recovered : result.recoveredConstraints) {
            if (!(recovered.status == bnd::SteinerStatus::AlreadyPresent ||
                  recovered.status == bnd::SteinerStatus::Recovered) ||
                recovered.chain.size() < 2)
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};
            for (std::size_t i = 0; i + 1 < recovered.chain.size(); ++i) {
                const std::size_t a = recovered.chain[i];
                const std::size_t b = recovered.chain[i + 1];
                if (a >= result.localHandles.size() ||
                    b >= result.localHandles.size())
                    return {FacePublishStatus::Rejected,
                            FacePublishReason::InvalidResult, {}};
                const auto edge = ws.edges().find(
                    ws.nodes(), result.localHandles[a],
                    result.localHandles[b]);
                if (!edge)
                    return {FacePublishStatus::Rejected,
                            FacePublishReason::InvalidResult, {}};
                const auto data = ws.edges().get(*edge);
                if (!data.requiredConstraint ||
                    !sameClass(data.cls, recovered.classification))
                    return {FacePublishStatus::Rejected,
                            FacePublishReason::InvalidResult, {}};
                recoveredEdges.insert(localEdge(a, b));
            }
        }
        for (const auto& edge : ws.edges().activeBySequence()) {
            const auto data = ws.edges().get(edge);
            if (!data.requiredConstraint) continue;
            const auto a = indexOf.find(data.a);
            const auto b = indexOf.find(data.b);
            if (a == indexOf.end() || b == indexOf.end() ||
                recoveredEdges.count(localEdge(a->second, b->second)) == 0)
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};
        }

        for (const std::size_t local : result.outerLoop)
            if (local >= result.localNodes.size() ||
                result.existingGlobalIndex.find(local) ==
                    result.existingGlobalIndex.end())
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};

        const std::size_t newNodeCount =
            result.localNodes.size() - result.existingGlobalIndex.size();
        if (newNodeCount >
            static_cast<std::size_t>(
                std::numeric_limits<Index>::max() - global.nodeCount()))
            return {FacePublishStatus::Rejected,
                    FacePublishReason::InvalidResult, {}};

        std::vector<Index> localToGlobal(
            result.localNodes.size(), kInvalidIndex);
        Index next = global.nodeCount();
        for (std::size_t i = 0; i < result.localNodes.size(); ++i) {
            const auto existing = result.existingGlobalIndex.find(i);
            if (existing != result.existingGlobalIndex.end())
                localToGlobal[i] = existing->second;
            else
                localToGlobal[i] = next++;
        }
        std::vector<std::vector<Index>> connectivity;
        connectivity.reserve(result.localCells.size());
        for (const auto& cell : result.localCells) {
            if (!sameClass(cell.classification, result.faceClass) ||
                cell.nodes[0] >= localToGlobal.size() ||
                cell.nodes[1] >= localToGlobal.size() ||
                cell.nodes[2] >= localToGlobal.size() ||
                cell.nodes[0] == cell.nodes[1] ||
                cell.nodes[1] == cell.nodes[2] ||
                cell.nodes[2] == cell.nodes[0])
                return {FacePublishStatus::Rejected,
                        FacePublishReason::InvalidResult, {}};
            connectivity.push_back(
                {localToGlobal[cell.nodes[0]], localToGlobal[cell.nodes[1]],
                 localToGlobal[cell.nodes[2]]});
        }
        if (result.localCells.size() >
            std::numeric_limits<std::size_t>::max() / 3)
            return {FacePublishStatus::Rejected,
                    FacePublishReason::InvalidResult, {}};

        // The only mutation before commit is hidden vector-capacity growth.
        global.reserveAppend(newNodeCount, result.localCells.size(),
                             result.localCells.size() * 3);

        // All buffers and mappings are prepared. With trivially-copyable
        // payloads and sufficient capacity, these required addNode/addCell
        // calls perform no allocation and cannot fail. Keep the bounded tail
        // noexcept so a hypothetical contract violation cannot be mistaken
        // for a recoverable partial publication.
        auto commit = [&]() noexcept {
            for (std::size_t i = 0; i < result.localNodes.size(); ++i)
                if (result.existingGlobalIndex.find(i) ==
                    result.existingGlobalIndex.end()) {
                    const Index added = global.addNode(
                        result.localNodes[i].coord,
                        result.localNodes[i].cls);
                    if (added != localToGlobal[i]) std::terminate();
                }
            for (std::size_t i = 0; i < result.localCells.size(); ++i)
                global.addCell(mesh::CellType::Tri3, connectivity[i],
                               result.localCells[i].classification);
        };
        commit();
        FacePublication publication;
        publication.status = FacePublishStatus::Published;
        publication.reason = FacePublishReason::None;
        publication.localToGlobal = std::move(localToGlobal);
        return publication;
    } catch (const std::bad_alloc&) {
        return {FacePublishStatus::Failed,
                FacePublishReason::AllocationFailed, {}};
    } catch (...) {
        return {FacePublishStatus::Failed,
                FacePublishReason::InternalError, {}};
    }
}

}  // namespace paacs::algo
