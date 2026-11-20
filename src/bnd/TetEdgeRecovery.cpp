// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/bnd/TetFacetSetRecovery.hpp"

#include "paacs/bnd/TetRecoveryIndex.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/topo/TetRefinement.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <map>
#include <new>
#include <queue>
#include <set>
#include "TetEdgeGeneralCavityCadPolicy.hpp"
#include "TetEdgeManifoldClosure.hpp"

namespace paacs::bnd {

#if defined(PAACS_TOPO_TESTING)
namespace {
TetEdgeCavityDiagnostics gTetEdgeCavityDiagnostics{};
TetEdgeTransactionDiagnostics gTetEdgeTransactionDiagnostics{};
std::atomic_bool gTetEdgeCavityDiagnosticsEnabled{false};
std::atomic_bool gTetEdgeTransactionDiagnosticsEnabled{false};
std::atomic_flag gTetEdgeCavityDiagnosticsLock = ATOMIC_FLAG_INIT;
struct TetEdgeCavityDiagnosticsLock {
    TetEdgeCavityDiagnosticsLock() noexcept {
        while(gTetEdgeCavityDiagnosticsLock.test_and_set(
                  std::memory_order_acquire)) {}
    }
    ~TetEdgeCavityDiagnosticsLock() {
        gTetEdgeCavityDiagnosticsLock.clear(std::memory_order_release);
    }
};
}

void resetTetEdgeCavityDiagnosticsForTesting(
    bool enabled, VolumeArtifactEdgeKey target) noexcept {
    TetEdgeCavityDiagnosticsLock lock;
    TetEdgeCavityDiagnostics fresh;
    fresh.enabled = enabled;
    fresh.targetGlobal = target;
    gTetEdgeCavityDiagnostics = std::move(fresh);
    gTetEdgeCavityDiagnosticsEnabled.store(enabled,std::memory_order_release);
}

const TetEdgeCavityDiagnostics&
tetEdgeCavityDiagnosticsForTesting() noexcept {
    return gTetEdgeCavityDiagnostics;
}

void resetTetEdgeTransactionDiagnosticsForTesting(
    bool enabled,VolumeArtifactEdgeKey target) noexcept {
    TetEdgeCavityDiagnosticsLock lock;
    gTetEdgeTransactionDiagnostics={};
    gTetEdgeTransactionDiagnostics.enabled=enabled;
    gTetEdgeTransactionDiagnostics.targetGlobal=target;
    gTetEdgeTransactionDiagnosticsEnabled.store(
        enabled,std::memory_order_release);
}

TetEdgeTransactionDiagnostics
tetEdgeTransactionDiagnosticsForTesting() noexcept {
    TetEdgeCavityDiagnosticsLock lock;
    return gTetEdgeTransactionDiagnostics;
}
#endif

namespace {

using predicates::Sign;

bool sameClass(const mesh::Classification& a,
               const mesh::Classification& b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}

TetEdgeRecoveryOutcome fail(TetEdgeRecoveryStatus status,
                            TetEdgeRecoveryReason reason,
                            TetEdgeRecoveryPhase phase) noexcept {
    TetEdgeRecoveryOutcome out;
    out.status = status;
    out.reason = reason;
    out.phase = phase;
    return out;
}

constexpr std::array<std::array<int,2>,6> kEdges{{
    {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
constexpr std::array<std::array<int,3>,4> kFaces{{
    {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};

std::array<Vec3,4> points(const TetRecoveryContext& context,
                          const CurrentTetRecord& cell,
                          const topo::TopologyWorkspace& workspace) {
    std::array<Vec3,4> result{};
    for (std::size_t i = 0; i < 4; ++i)
        result[i] = workspace.nodes().get(
            context.denseToNode[cell.nodes[i]]).coord;
    return result;
}

bool endpointWedge(const Vec3& endpoint, const Vec3& other,
                   const std::array<Vec3,4>& tet,
                   const std::array<DenseTetNodeId,4>& nodes,
                   DenseTetNodeId endpointId) {
    if (std::find(nodes.begin(), nodes.end(), endpointId) == nodes.end())
        return false;
    bool constrained = false;
    for (const auto& face : kFaces) {
        const Sign atEndpoint = predicates::orient3d(
            tet[static_cast<std::size_t>(face[0])],
            tet[static_cast<std::size_t>(face[1])],
            tet[static_cast<std::size_t>(face[2])], endpoint);
        if (atEndpoint != Sign::Zero) continue;
        constrained = true;
        if (predicates::orient3d(
                tet[static_cast<std::size_t>(face[0])],
                tet[static_cast<std::size_t>(face[1])],
                tet[static_cast<std::size_t>(face[2])], other) ==
            Sign::Positive)
            return false;
    }
    return constrained;
}

bool segmentIntersectsCellOpen(
    const Vec3& a, const Vec3& b, DenseTetNodeId aId, DenseTetNodeId bId,
    const TetRecoveryContext& context, const CurrentTetRecord& cell,
    const topo::TopologyWorkspace& workspace) {
    const auto tet = points(context, cell, workspace);
    int crossings = 0;
    for (const auto& face : kFaces)
        crossings += tet_facet_detail::segmentCrossesOpenTriangle(
            a, b, tet[static_cast<std::size_t>(face[0])],
            tet[static_cast<std::size_t>(face[1])],
            tet[static_cast<std::size_t>(face[2])]) ? 1 : 0;
    if (crossings >= 2) return true;
    return endpointWedge(a,b,tet,cell.nodes,aId) ||
           endpointWedge(b,a,tet,cell.nodes,bId);
}

std::array<topo::NodeHandle,4> positive(
    const topo::TopologyWorkspace& workspace,
    std::array<topo::NodeHandle,4> cell) {
    if (predicates::orient3d(workspace.nodes().get(cell[0]).coord,
                             workspace.nodes().get(cell[1]).coord,
                             workspace.nodes().get(cell[2]).coord,
                             workspace.nodes().get(cell[3]).coord) ==
        Sign::Negative)
        std::swap(cell[0], cell[1]);
    return cell;
}

bool addBindingWork(std::uint64_t& target,std::uint64_t amount,
                    TetRecoveryBindingWorkCounters* work) noexcept {
    if(!work)return true;
    if(amount>(std::numeric_limits<std::uint64_t>::max)()-target) {
        work->complete=false;
        return false;
    }
    target+=amount;
    return true;
}

void recordStateCopy(const TetFacetSetState& state,
                     TetRecoveryBindingWorkCounters* work) noexcept {
    if(!work)return;
    (void)addBindingWork(work->stateCopies,1,work);
    const std::array<std::size_t,4> rows{
        state.edgeBindings.size(),state.bindings.size(),
        state.canonicalEdgeSchedule.size(),state.canonicalSchedule.size()};
    for(const auto count:rows) {
        if(count>(std::numeric_limits<std::uint64_t>::max)()) {
            work->complete=false;
            continue;
        }
        (void)addBindingWork(work->stateCopyRows,
            static_cast<std::uint64_t>(count),work);
    }
}

bool validState(const TetFacetSetState& state,
                const topo::TopologyWorkspace& workspace,
                TetRecoveryBindingWorkCounters* work) {
    if(work)(void)addBindingWork(work->edgeBindingAuditCalls,1,work);
    std::set<topo::EdgeHandle> handles;
    for (const auto& [key,binding] : state.edgeBindings) {
        if(work)(void)addBindingWork(work->edgeBindingRowsVisited,1,work);
        if (!(key == binding.key) ||
            !workspace.edges().isCurrent(binding.edge) ||
            !handles.insert(binding.edge).second)
            return false;
        const auto data = workspace.edges().get(binding.edge);
        if(work)(void)addBindingWork(work->edgeDataCopyRows,1,work);
        if (!data.requiredConstraint ||
            !sameClass(data.cls,binding.classification))
            return false;
    }
    return true;
}

void retainProtectedEdge(TetRecoveryContext& context,
                         const DenseTetEdgeKey& key) {
    const auto value = canonicalTetFacetEdge(key.nodes[0],key.nodes[1]);
    if (std::find(context.protectedEdges.begin(),context.protectedEdges.end(),
                  value) == context.protectedEdges.end()) {
        context.protectedEdges.push_back(value);
        std::sort(context.protectedEdges.begin(),context.protectedEdges.end());
    }
}

// Build one immutable local identity map. Coordinates are never used to merge
// nodes. The retained solver certifies the complete boundary, protected
// features, positive cells, cyclic target star and exact pairwise contacts.
TetEdgeRecoveryOutcome buildGeneralEdgeCandidate(
    const std::set<topo::CellHandle>& selected,
    const DenseTetEdgeKey& target, const TetRecoveryContext& context,
    const TetEdgeGeneralCavityLimits& limits, topo::TetEdgeCavitySpec& spec,
    TetEdgeGeneralRecoveryDiagnostics* diagnostics,
    topo::TetRefinementCavitySpec* staged = nullptr) {
    if (diagnostics) *diagnostics = {};
    if (staged) *staged = {};
    if (selected.size() > limits.maxCavityCells)
        return fail(TetEdgeRecoveryStatus::Refused,
                    TetEdgeRecoveryReason::ResourceBudgetReached,
                    TetEdgeRecoveryPhase::Candidate);
    std::set<DenseTetNodeId> denseSet;
    for (const auto cell : selected) {
        const auto& record = context.cells.at(cell);
        denseSet.insert(record.nodes.begin(),record.nodes.end());
        if (denseSet.size() > limits.maxCavityVertices)
            return fail(TetEdgeRecoveryStatus::Refused,
                        TetEdgeRecoveryReason::ResourceBudgetReached,
                        TetEdgeRecoveryPhase::Candidate);
    }
    std::vector<DenseTetNodeId> dense(denseSet.begin(),denseSet.end());
    std::sort(dense.begin(),dense.end(),[&](const auto a,const auto b) {
        return std::tie(context.plainVertices[a].sequence,a) <
               std::tie(context.plainVertices[b].sequence,b);
    });
    std::map<DenseTetNodeId,TetFacetVertexId> local;
    std::vector<TetFacetVertex> vertices;
    vertices.reserve(dense.size());
    for (const auto node : dense) {
        local.emplace(node,static_cast<TetFacetVertexId>(vertices.size()));
        vertices.push_back(context.plainVertices[node]);
    }
    if (!local.contains(target.nodes[0]) || !local.contains(target.nodes[1]))
        return fail(TetEdgeRecoveryStatus::Refused,
                    TetEdgeRecoveryReason::NonCycleLink,
                    TetEdgeRecoveryPhase::Candidate);
    std::vector<TetFacetCell> removed;
    removed.reserve(selected.size());
    for (const auto cell : selected) {
        const auto& record = context.cells.at(cell);
        std::array<TetFacetVertexId,4> nodes{};
        for (std::size_t i = 0; i < 4; ++i) nodes[i] = local.at(record.nodes[i]);
        removed.push_back({nodes,record.sequence,record.classification});
    }
    std::sort(removed.begin(),removed.end(),[](const auto& a,const auto& b) {
        return a.sequence < b.sequence;
    });
    // A local vertex set can also span constraints owned entirely by untouched
    // cells. Preserve those in place; do not demand they appear in the refill.
    std::set<TetFacetEdgeKey> replacedEdges;
    std::set<TetFacetFaceKey> replacedFaces;
    for (const auto cell : selected) {
        const auto& row = context.cells.at(cell);
        for (const auto& edge : kEdges)
            replacedEdges.insert(canonicalTetFacetEdge(row.nodes[edge[0]],row.nodes[edge[1]]));
        for (const auto& face : kFaces)
            replacedFaces.insert(canonicalTetFacetFace({row.nodes[face[0]],row.nodes[face[1]],row.nodes[face[2]]}));
    }
    std::vector<TetFacetEdgeKey> protectedEdges;
    for (const auto edge : context.protectedEdges)
        if (replacedEdges.contains(edge))
            protectedEdges.push_back(canonicalTetFacetEdge(
                local.at(edge.lo),local.at(edge.hi)));
    std::vector<TetFacetProtectedFace> protectedFaces;
    for (const auto& face : context.protectedFaces) {
        if (!replacedFaces.contains(face.key)) continue;
        std::array<TetFacetVertexId,3> nodes{};
        bool complete = true;
        for (std::size_t i = 0; i < 3; ++i) {
            const auto found = local.find(face.key.vertices[i]);
            if (found == local.end()) { complete = false; break; }
            nodes[i] = found->second;
        }
        if (complete) protectedFaces.push_back(
            {canonicalTetFacetFace(nodes),face.identity,face.classification});
    }
    std::sort(protectedEdges.begin(),protectedEdges.end());
    std::sort(protectedFaces.begin(),protectedFaces.end(),[](const auto& a,const auto& b) {
        return a.key < b.key;
    });
    const TetEdgeGeneralCavityInput input{
        vertices,removed,canonicalTetFacetEdge(local.at(target.nodes[0]),
                                             local.at(target.nodes[1])),
        protectedEdges,protectedFaces,spec.region};
    auto solved = solveTetEdgeGeneralCavity(input,limits);
    if (staged && limits.maximumSteinerCandidates &&
        solved.status == TetEdgeGeneralCavityStatus::Unavailable &&
        solved.metrics.runComplete && solved.metrics.workUnits < limits.maxWorkUnits) {
        const auto existingWork = solved.metrics.workUnits;
        auto bounded = limits;
        bounded.maxWorkUnits -= existingWork;
        auto added = solveTetEdgeGeneralCavityOneSteiner(input,bounded,
            limits.maximumSteinerCandidates,bounded.maxWorkUnits);
        if (diagnostics) diagnostics->steinerCandidateAttempts = added.candidateAttempts;
        solved = std::move(added.cavity);
        solved.status = added.status;
        solved.geometryReason = added.geometryReason;
        solved.metrics.workUnits = existingWork + added.cumulativeWorkUnits;
        if (added.solved()) {
            staged->stamp = spec.stamp;
            staged->stagedNode = {0,added.steinerPosition,spec.region};
            staged->removed = spec.removed;
            staged->recoverTargetEdge = true;
            staged->targetEdge = spec.target;
            staged->targetClassification = spec.targetClassification;
            for (const auto& candidate : solved.cells) {
                topo::TetRefinementCellSpec cell;
                cell.classification = spec.region;
                for (std::size_t i=0;i<4;++i)
                    cell.vertices[i] = candidate.vertices[i] == vertices.size()
                        ? topo::TetRefinementNodeRef::stagedNode()
                        : topo::TetRefinementNodeRef::existingNode(
                            context.denseToNode[dense.at(candidate.vertices[i])]);
                staged->added.push_back(cell);
            }
            if (diagnostics) diagnostics->stagedCandidateAvailable = true;
        }
    }
    if (diagnostics) {
        diagnostics->attempted = true;
        diagnostics->status = solved.status;
        diagnostics->geometryReason = solved.geometryReason;
        diagnostics->metrics = solved.metrics;
    }
    if (!solved.solved()) {
        auto reason = TetEdgeRecoveryReason::NonCycleLink;
        auto status = TetEdgeRecoveryStatus::Refused;
        if (solved.status == TetEdgeGeneralCavityStatus::BudgetReached)
            reason = TetEdgeRecoveryReason::ResourceBudgetReached;
        else if (solved.status == TetEdgeGeneralCavityStatus::AllocationFailed) {
            status = TetEdgeRecoveryStatus::Failed;
            reason = TetEdgeRecoveryReason::AllocationFailed;
        } else if (solved.status == TetEdgeGeneralCavityStatus::InvalidInput ||
                   solved.status == TetEdgeGeneralCavityStatus::InvalidCandidate)
            reason = TetEdgeRecoveryReason::InvalidCandidate;
        return fail(status,reason,TetEdgeRecoveryPhase::Candidate);
    }
    if (staged && staged->recoverTargetEdge) {
        TetEdgeRecoveryOutcome out;
        out.status = TetEdgeRecoveryStatus::Recovered;
        return out;
    }
    spec.added.reserve(solved.cells.size());
    for (const auto& candidate : solved.cells) {
        topo::TetCavityCell cell;
        for (std::size_t i = 0; i < 4; ++i)
            cell.nodes[i] = context.denseToNode[dense.at(candidate.vertices[i])];
        cell.cls = spec.region;
        spec.added.push_back(cell);
    }
    TetEdgeRecoveryOutcome out;
    out.status = TetEdgeRecoveryStatus::Recovered;
    return out;
}

TetEdgeRecoveryOutcome buildBoundedGeneralEdgeCandidate(
    std::set<topo::CellHandle>& selected, const DenseTetEdgeKey& target,
    const TetRecoveryContext& context, const TetEdgeGeneralCavityLimits& limits,
    std::uint64_t maximumAdjacentAttempts, topo::TetEdgeCavitySpec& spec,
    TetEdgeGeneralRecoveryDiagnostics* diagnostics,
    topo::TetRefinementCavitySpec* staged = nullptr) {
    TetEdgeGeneralRecoveryDiagnostics evidence;
    auto result = buildGeneralEdgeCandidate(selected,target,context,limits,spec,&evidence,staged);
    if (diagnostics) *diagnostics = evidence;
    if (detail::cadGeneralPruningEnabled() && maximumAdjacentAttempts &&
        evidence.attempted && evidence.status == TetEdgeGeneralCavityStatus::InvalidInput &&
        (evidence.geometryReason == TetFacetGeometryReason::BoundaryNotClosed ||
         evidence.geometryReason == TetFacetGeometryReason::InvalidVertexLink)) {
        const auto initialWork = evidence.metrics.workUnits;
        const auto closure = detail::prepareTetEdgeManifoldClosure(selected,context,spec.region,
            limits.maxCavityCells,limits.maxCavityVertices,limits.maxWorkUnits-initialWork,
            static_cast<unsigned>((std::min)(maximumAdjacentAttempts,std::uint64_t{32})));
        evidence.metrics.workUnits += closure.work;
        evidence.manifoldClosurePasses = closure.passes;
        if (closure.status == detail::ManifoldClosureStatus::Closed) {
            auto proposal = spec;
            proposal.removed.assign(closure.cells.begin(),closure.cells.end());
            auto bounded = limits;
            bounded.maxWorkUnits -= evidence.metrics.workUnits;
            if (bounded.maxWorkUnits) {
                TetEdgeGeneralRecoveryDiagnostics repaired;
                result = buildGeneralEdgeCandidate(closure.cells,target,context,bounded,proposal,&repaired,staged);
                repaired.metrics.workUnits += evidence.metrics.workUnits;
                repaired.manifoldClosurePasses = closure.passes;
                repaired.manifoldCellsAdded = closure.cells.size()-selected.size();
                if (diagnostics) *diagnostics = repaired;
                if (result.succeeded()) { selected = closure.cells; spec = std::move(proposal); }
                return result;
            }
        }
        if (closure.status == detail::ManifoldClosureStatus::BudgetReached ||
            evidence.metrics.workUnits == limits.maxWorkUnits) {
            evidence.status = TetEdgeGeneralCavityStatus::BudgetReached;
            evidence.metrics.runComplete = false;
            result = fail(TetEdgeRecoveryStatus::Refused,TetEdgeRecoveryReason::ResourceBudgetReached,
                          TetEdgeRecoveryPhase::Candidate);
        }
        if (diagnostics) *diagnostics = evidence;
        return result;
    }
    if (result.succeeded() || !maximumAdjacentAttempts || !evidence.attempted ||
        evidence.status != TetEdgeGeneralCavityStatus::Unavailable ||
        !evidence.metrics.runComplete) return result;

    // Only extend a certified no-solution cavity. Every proposal uses the
    // immutable original one-ring, never a chain of speculative expansions.
    // All exact solver calls and neighbour-face visits share one work budget.
    std::uint64_t remaining = limits.maxWorkUnits - evidence.metrics.workUnits;
    std::uint64_t attempts = 0;
    bool sawBudget = false;
    const auto finish = [&](TetEdgeRecoveryOutcome outcome,std::uint64_t added) {
        evidence.metrics.workUnits = limits.maxWorkUnits - remaining;
        evidence.adjacentExpansionAttempts = attempts;
        evidence.adjacentCellsAdded = added;
        if (diagnostics) *diagnostics = evidence;
        return outcome;
    };
    const auto budget = [&]() {
        evidence.status = TetEdgeGeneralCavityStatus::BudgetReached;
        evidence.metrics.runComplete = false;
        return finish(fail(TetEdgeRecoveryStatus::Refused,
                           TetEdgeRecoveryReason::ResourceBudgetReached,
                           TetEdgeRecoveryPhase::Candidate),0);
    };
    std::set<topo::CellHandle> neighbourSet;
    for (const auto cell : selected) {
        const auto& record = context.cells.at(cell);
        for (const auto& face : kFaces) {
            if (!remaining) return budget();
            --remaining;
            const auto key = DenseTetFaceKey::of({
                record.nodes[static_cast<std::size_t>(face[0])],
                record.nodes[static_cast<std::size_t>(face[1])],
                record.nodes[static_cast<std::size_t>(face[2])]});
            const auto found = context.faces.find(key);
            if (found == context.faces.end() || found->second.owners.size() != 2) continue;
            const auto plain = canonicalTetFacetFace(key.nodes);
            bool protectedFace = false;
            for (const auto& protectedRow : context.protectedFaces) {
                if (!remaining) return budget();
                --remaining;
                if (protectedRow.key == plain) { protectedFace = true; break; }
            }
            if (protectedFace) continue;
            for (const auto& owner : found->second.owners)
                if (!selected.contains(owner.cell) &&
                    sameClass(context.cells.at(owner.cell).classification,spec.region))
                    neighbourSet.insert(owner.cell);
        }
    }
    std::vector<topo::CellHandle> neighbours(neighbourSet.begin(),neighbourSet.end());
    std::sort(neighbours.begin(),neighbours.end(),[&](const auto a,const auto b) {
        return std::tie(context.cells.at(a).sequence,a) <
               std::tie(context.cells.at(b).sequence,b);
    });
    const auto attempt = [&](topo::CellHandle first,topo::CellHandle second,
                             bool pair) {
        ++attempts;
        auto expanded = selected;
        expanded.insert(first);
        if (pair) expanded.insert(second);
        topo::TetEdgeCavitySpec proposal = spec;
        proposal.removed.assign(expanded.begin(),expanded.end());
        auto bounded = limits;
        bounded.maxWorkUnits = remaining;
        result = buildGeneralEdgeCandidate(expanded,target,context,bounded,proposal,&evidence,staged);
        remaining -= (std::min)(remaining,evidence.metrics.workUnits);
        sawBudget = sawBudget || result.reason == TetEdgeRecoveryReason::ResourceBudgetReached;
        if (result.succeeded()) {
            selected.swap(expanded);
            spec = std::move(proposal);
        }
    };
    for (std::size_t i = 0; i < neighbours.size(); ++i) {
        if (attempts == maximumAdjacentAttempts || !remaining) return budget();
        attempt(neighbours[i],{},false);
        if (result.succeeded()) return finish(result,1);
        if (result.status == TetEdgeRecoveryStatus::Failed) return finish(result,0);
    }
    for (std::size_t i = 0; i < neighbours.size(); ++i)
        for (std::size_t j = i+1; j < neighbours.size(); ++j) {
            if (attempts == maximumAdjacentAttempts || !remaining) return budget();
            attempt(neighbours[i],neighbours[j],true);
            if (result.succeeded()) return finish(result,2);
            if (result.status == TetEdgeRecoveryStatus::Failed) return finish(result,0);
        }
    if (sawBudget) return budget();
    return finish(fail(TetEdgeRecoveryStatus::Refused,TetEdgeRecoveryReason::NonCycleLink,
                       TetEdgeRecoveryPhase::Candidate),0);
}

// Preserve the old staged-node transaction, with current state publication
// and context failure semantics. All throwing preparation remains local.
TetEdgeRecoveryOutcome publishSteinerEdgeCandidate(
    const TetEdgeRecoveryRequest& request,
    const topo::TetRefinementCavitySpec& spec,
    TetFacetSetState& state, topo::TopologyWorkspace& workspace,
    TetRecoveryContext& context, TetRecoveryBindingWorkCounters* bindingWork,
    std::vector<topo::CellHandle>* createdWitness) noexcept {
    bool topologyCommitted = false;
    try {
        auto transaction = workspace.beginTetRefinementTransaction();
        (void)transaction.stage(spec);
        const auto& prediction = transaction.prepare();
        TetFacetSetState scratch = state;
        recordStateCopy(state,bindingWork);
        scratch.edgeBindings[request.key] =
            {request.key,request.classification,true,prediction.targetEdge};

        // The new identity and index storage are prepared in a private copy.
        // No partially extended live context is visible on preparation failure.
        auto next = std::make_unique<TetRecoveryContext>(context);
        if (next->denseToNode.size() >= kInvalidTetFacetVertex)
            return fail(TetEdgeRecoveryStatus::Refused,
                TetEdgeRecoveryReason::ResourceBudgetReached,TetEdgeRecoveryPhase::Preparation);
        const auto dense = static_cast<DenseTetNodeId>(next->denseToNode.size());
        next->denseToNode.push_back(prediction.stagedNode);
        next->plainVertices.push_back({spec.stagedNode.coordinate,0});
        next->nodeToDense.emplace(prediction.stagedNode,dense);
        next->nodeUseCount.push_back(0);
        next->nodeStars.emplace_back();
        next->metrics.nodes = next->denseToNode.size();
        std::vector<topo::CellHandle> created;
        created.reserve(spec.added.size());
        for (const auto& cell : prediction.cells)
            if (cell.transition == topo::TetRefinementTransition::Created)
                created.push_back(cell.predicted);

        const auto committed = transaction.commit();
        topologyCommitted = true;
        if (committed.actualNode != committed.predictedNode ||
            committed.actualTargetEdge != committed.predictedTargetEdge) std::abort();
        next->plainVertices.back().sequence = workspace.nodes().seqOf(committed.actualNode);
        const auto updated = updateTetRecoveryContext(*next,workspace,spec.removed,created);
        if (!updated.updated()) {
            context.poisoned = true;
            return fail(TetEdgeRecoveryStatus::Failed,
                updated.reason == TetRecoveryContextReason::AllocationFailed
                    ? TetEdgeRecoveryReason::AllocationFailed
                    : TetEdgeRecoveryReason::StaleRecoveryContext,
                TetEdgeRecoveryPhase::StateCopy);
        }
        retainProtectedEdge(*next,DenseTetEdgeKey::of(
            next->nodeToDense.at(request.vertices[0]),next->nodeToDense.at(request.vertices[1])));
        if (createdWitness) *createdWitness = created;
        scratch.workspaceStamp = topo::WorkspaceStamp::of(workspace);
        static_assert(std::is_nothrow_move_assignable_v<TetRecoveryContext>);
        context = std::move(*next);
        state.swap(scratch);
        TetEdgeRecoveryOutcome out;
        out.status = TetEdgeRecoveryStatus::Recovered;
        out.edge = committed.actualTargetEdge;
        out.stagedNodeCreated = true;
        return out;
    } catch (const std::bad_alloc&) {
        if (topologyCommitted) context.poisoned = true;
        if (bindingWork) bindingWork->complete = false;
        return fail(TetEdgeRecoveryStatus::Failed,TetEdgeRecoveryReason::AllocationFailed,
                    TetEdgeRecoveryPhase::Preparation);
    } catch (...) {
        if (topologyCommitted) context.poisoned = true;
        return fail(TetEdgeRecoveryStatus::Failed,TetEdgeRecoveryReason::PreparationConflict,
                    TetEdgeRecoveryPhase::Preparation);
    }
}

#if defined(PAACS_TOPO_TESTING)
std::array<std::uint32_t,2> diagnosticEdge(
    const DenseTetEdgeKey& edge) noexcept {
    return {edge.nodes[0],edge.nodes[1]};
}

std::array<std::uint32_t,3> diagnosticFace(
    const DenseTetFaceKey& face) noexcept {
    return {face.nodes[0],face.nodes[1],face.nodes[2]};
}

TetFacetFaceKey plainFace(const DenseTetFaceKey& face) noexcept {
    return canonicalTetFacetFace(
        {face.nodes[0],face.nodes[1],face.nodes[2]});
}

bool diagnosticProtectedFace(const TetRecoveryContext& context,
                             const DenseTetFaceKey& face) noexcept {
    const auto key = plainFace(face);
    return std::find_if(
        context.protectedFaces.begin(),context.protectedFaces.end(),
        [&](const TetFacetProtectedFace& row) { return row.key == key; }) !=
        context.protectedFaces.end();
}

void sortDiagnosticCells(std::vector<topo::CellHandle>& cells,
                         const TetRecoveryContext& context) {
    std::sort(cells.begin(),cells.end(),[&](const auto a,const auto b) {
        const auto& ca=context.cells.at(a);
        const auto& cb=context.cells.at(b);
        return std::tie(ca.sequence,a)<std::tie(cb.sequence,b);
    });
}

void captureTetEdgeCavityDiagnostics(
    const TetEdgeRecoveryRequest& request, const DenseTetEdgeKey& target,
    std::span<const topo::CellHandle> initiallySelected,
    const std::set<topo::CellHandle>& selected,
    const TetRecoveryContext& context,
    const topo::TopologyWorkspace& workspace) noexcept {
    if (!gTetEdgeCavityDiagnosticsEnabled.load(std::memory_order_acquire))
        return;
    {
        TetEdgeCavityDiagnosticsLock lock;
        const auto requested=gTetEdgeCavityDiagnostics.targetGlobal;
        if(requested.canonicalGlobalNodes[0]!=kInvalidIndex&&
           !(requested==request.key))return;
    }
    try {
        TetEdgeCavityDiagnostics next;
        next.enabled=true;
        next.targetGlobal=request.key;
        next.targetDense=diagnosticEdge(target);
        next.initiallySelected.assign(
            initiallySelected.begin(),initiallySelected.end());
        sortDiagnosticCells(next.initiallySelected,context);
        next.selected.assign(selected.begin(),selected.end());
        sortDiagnosticCells(next.selected,context);
        std::set<topo::CellHandle> initialSet(
            initiallySelected.begin(),initiallySelected.end());
        for(const auto cell:next.selected)
            if(initialSet.count(cell)==0)next.closureAdded.push_back(cell);

        std::set<DenseTetNodeId> initialVertexIds;
        for(const auto cell:next.initiallySelected) {
            const auto& record=context.cells.at(cell);
            next.initialCells.push_back(
                {cell,record.sequence,record.nodes});
            initialVertexIds.insert(record.nodes.begin(),record.nodes.end());
        }
        for(const auto dense:initialVertexIds) {
            const auto node=context.denseToNode[dense];
            next.initialVertices.push_back({
                dense,workspace.nodes().seqOf(node),
                workspace.nodes().get(node).coord});
        }

        std::map<DenseTetEdgeKey,std::size_t> ringMultiplicities;
        for(const auto cell:next.selected) {
            const auto& record=context.cells.at(cell);
            for(const auto& localFace:kFaces) {
                const auto face=DenseTetFaceKey::of({
                    record.nodes[static_cast<std::size_t>(localFace[0])],
                    record.nodes[static_cast<std::size_t>(localFace[1])],
                    record.nodes[static_cast<std::size_t>(localFace[2])]});
                const auto resolver=context.faces.find(face);
                if(resolver==context.faces.end())continue;
                std::size_t selectedOwners=0;
                for(const auto& owner:resolver->second.owners)
                    selectedOwners+=selected.count(owner.cell)!=0?1U:0U;
                if(selectedOwners!=1)continue;
                std::array<DenseTetNodeId,2> ring{};
                std::size_t ringSize=0;
                std::size_t targets=0;
                for(const auto node:face.nodes) {
                    if(node==target.nodes[0]||node==target.nodes[1])++targets;
                    else if(ringSize<ring.size())ring[ringSize++]=node;
                }
                if(targets<next.boundaryFacesByTargetEndpoints.size())
                    ++next.boundaryFacesByTargetEndpoints[targets];
                TetEdgeBoundaryFaceDiagnostic row;
                row.denseNodes=diagnosticFace(face);
                row.directedDenseNodes={
                    record.nodes[static_cast<std::size_t>(localFace[0])],
                    record.nodes[static_cast<std::size_t>(localFace[1])],
                    record.nodes[static_cast<std::size_t>(localFace[2])]};
                row.targetEndpointCount=static_cast<std::uint8_t>(targets);
                row.currentOwnerCount=static_cast<std::uint8_t>(
                    (std::min)(resolver->second.owners.size(),
                               static_cast<std::size_t>(255)));
                row.protectedFace=diagnosticProtectedFace(context,face);
                row.regionBoundary=resolver->second.owners.size()==1;
                next.boundaryFaces.push_back(row);
                if(row.protectedFace)
                    next.protectedBoundaryFaces.push_back(row.denseNodes);
                if(row.regionBoundary)
                    next.regionBoundaryFaces.push_back(row.denseNodes);
                if(targets==1&&ringSize==2)
                    ++ringMultiplicities[
                        DenseTetEdgeKey::of(ring[0],ring[1])];
            }
        }

        std::set<DenseTetEdgeKey> scannedEdges;
        const auto p0=workspace.nodes().get(
            context.denseToNode[target.nodes[0]]).coord;
        const auto p1=workspace.nodes().get(
            context.denseToNode[target.nodes[1]]).coord;
        for(const auto cell:next.selected) {
            const auto& record=context.cells.at(cell);
            for(const auto& localEdge:kEdges) {
                const auto edge=DenseTetEdgeKey::of(
                    record.nodes[static_cast<std::size_t>(localEdge[0])],
                    record.nodes[static_cast<std::size_t>(localEdge[1])]);
                if(edge==target||!scannedEdges.insert(edge).second)continue;
                const auto a=workspace.nodes().get(
                    context.denseToNode[edge.nodes[0]]).coord;
                const auto b=workspace.nodes().get(
                    context.denseToNode[edge.nodes[1]]).coord;
                const auto contact=tet_facet_detail::segmentContact3(p0,p1,a,b);
                if(contact==tet_facet_detail::Segment2Contact::None)continue;
                const bool sharesEndpoint=edge.nodes[0]==target.nodes[0]||
                    edge.nodes[1]==target.nodes[0]||
                    edge.nodes[0]==target.nodes[1]||
                    edge.nodes[1]==target.nodes[1];
                if(contact==tet_facet_detail::Segment2Contact::Touch&&
                   sharesEndpoint)continue;
                const auto resolver=context.edges.find(edge);
                if(resolver!=context.edges.end()&&resolver->second.required)
                    next.protectedContactEdges.push_back(diagnosticEdge(edge));
            }
        }

        std::map<DenseTetNodeId,std::set<DenseTetNodeId>> link;
        bool ringMultiplicityValid=true;
        for(const auto& [edge,count]:ringMultiplicities) {
            next.ringEdges.push_back({diagnosticEdge(edge),count});
            ringMultiplicityValid=ringMultiplicityValid&&count==2;
            link[edge.nodes[0]].insert(edge.nodes[1]);
            link[edge.nodes[1]].insert(edge.nodes[0]);
        }
        std::map<DenseTetNodeId,std::size_t> componentByVertex;
        std::set<DenseTetNodeId> unseen;
        for(const auto& [vertex,neighbours]:link) {
            (void)neighbours;
            unseen.insert(vertex);
        }
        std::size_t component=0;
        while(!unseen.empty()) {
            std::queue<DenseTetNodeId> pending;
            pending.push(*unseen.begin());
            unseen.erase(unseen.begin());
            while(!pending.empty()) {
                const auto vertex=pending.front();
                pending.pop();
                componentByVertex[vertex]=component;
                for(const auto neighbour:link[vertex]) {
                    const auto it=unseen.find(neighbour);
                    if(it!=unseen.end()) {
                        pending.push(neighbour);
                        unseen.erase(it);
                    }
                }
            }
            ++component;
        }
        bool degreeValid=true;
        for(const auto& [vertex,neighbours]:link) {
            next.linkVertices.push_back(
                {vertex,neighbours.size(),componentByVertex.at(vertex)});
            degreeValid=degreeValid&&neighbours.size()==2;
        }
        next.linkComponents=component;
        next.linkIsSingleCycle=component==1&&link.size()>=3&&degreeValid&&
            ringMultiplicityValid&&next.boundaryFacesByTargetEndpoints[0]==0&&
            next.boundaryFacesByTargetEndpoints[2]==0;

        // T2-B.1 feasibility simulation: deterministically close only the
        // zero-target cap faces over existing same-region Tets. This is
        // diagnostic and never changes the workspace or accepted route.
        std::set<topo::CellHandle> completed=selected;
        // The current cell count is a finite exhaustive bound. Each pass
        // must add at least one previously unselected current cell.
        while(next.completionPasses<context.cells.size()&&
              completed.size()<context.cells.size()) {
            std::set<topo::CellHandle> additions;
            std::set<DenseTetFaceKey> visitedBoundaryFaces;
            for(const auto cell:completed) {
                const auto& record=context.cells.at(cell);
                for(const auto& localFace:kFaces) {
                    ++next.completionWorkUnits;
                    const auto face=DenseTetFaceKey::of({
                        record.nodes[static_cast<std::size_t>(localFace[0])],
                        record.nodes[static_cast<std::size_t>(localFace[1])],
                        record.nodes[static_cast<std::size_t>(localFace[2])]});
                    if(!visitedBoundaryFaces.insert(face).second)continue;
                    const auto resolver=context.faces.find(face);
                    if(resolver==context.faces.end()) {
                        next.completionBlockedByBoundary=true;
                        continue;
                    }
                    std::size_t selectedOwners=0;
                    topo::CellHandle outside{};
                    for(const auto& owner:resolver->second.owners) {
                        if(completed.count(owner.cell)!=0)++selectedOwners;
                        else outside=owner.cell;
                    }
                    if(selectedOwners!=1)continue;
                    std::size_t targets=0;
                    for(const auto node:face.nodes)
                        targets+=node==target.nodes[0]||
                                 node==target.nodes[1]?1U:0U;
                    if(targets!=0)continue;
                    if(diagnosticProtectedFace(context,face)) {
                        next.completionBlockedByProtectedFace=true;
                        continue;
                    }
                    if(resolver->second.owners.size()!=2||
                       !outside.isSpecified()) {
                        next.completionBlockedByBoundary=true;
                        continue;
                    }
                    const auto outsideCell=context.cells.find(outside);
                    if(outsideCell==context.cells.end()||
                       !sameClass(outsideCell->second.classification,
                                  context.region)) {
                        next.completionBlockedByRegion=true;
                        continue;
                    }
                    additions.insert(outside);
                }
            }
            if(additions.empty())break;
            ++next.completionPasses;
            for(const auto cell:additions)
                if(completed.insert(cell).second)
                    next.completionAdded.push_back(cell);
        }
        next.completionSelected.assign(completed.begin(),completed.end());
        sortDiagnosticCells(next.completionAdded,context);
        sortDiagnosticCells(next.completionSelected,context);

        std::map<DenseTetEdgeKey,std::size_t> completedRing;
        std::set<DenseTetFaceKey> completedBoundaryFaces;
        for(const auto cell:completed) {
            const auto& record=context.cells.at(cell);
            for(const auto& localFace:kFaces) {
                const auto face=DenseTetFaceKey::of({
                    record.nodes[static_cast<std::size_t>(localFace[0])],
                    record.nodes[static_cast<std::size_t>(localFace[1])],
                    record.nodes[static_cast<std::size_t>(localFace[2])]});
                const auto resolver=context.faces.find(face);
                if(resolver==context.faces.end())continue;
                std::size_t selectedOwners=0;
                for(const auto& owner:resolver->second.owners)
                    selectedOwners+=completed.count(owner.cell)!=0?1U:0U;
                if(selectedOwners!=1||
                   !completedBoundaryFaces.insert(face).second)continue;
                std::array<DenseTetNodeId,2> ring{};
                std::size_t ringSize=0;
                std::size_t targets=0;
                for(const auto node:face.nodes) {
                    if(node==target.nodes[0]||node==target.nodes[1])++targets;
                    else if(ringSize<ring.size())ring[ringSize++]=node;
                }
                if(targets<
                   next.completionBoundaryFacesByTargetEndpoints.size())
                    ++next.completionBoundaryFacesByTargetEndpoints[targets];
                if(targets==1&&ringSize==2)
                    ++completedRing[DenseTetEdgeKey::of(ring[0],ring[1])];
            }
        }
        std::map<DenseTetNodeId,std::set<DenseTetNodeId>> completedLink;
        bool completedMultiplicity=true;
        for(const auto& [edge,count]:completedRing) {
            completedMultiplicity=completedMultiplicity&&count==2;
            completedLink[edge.nodes[0]].insert(edge.nodes[1]);
            completedLink[edge.nodes[1]].insert(edge.nodes[0]);
        }
        bool completedDegrees=true;
        std::set<DenseTetNodeId> completedUnseen;
        for(const auto& [vertex,neighbours]:completedLink) {
            completedDegrees=completedDegrees&&neighbours.size()==2;
            completedUnseen.insert(vertex);
        }
        std::size_t completedComponents=0;
        while(!completedUnseen.empty()) {
            std::queue<DenseTetNodeId> pending;
            pending.push(*completedUnseen.begin());
            completedUnseen.erase(completedUnseen.begin());
            while(!pending.empty()) {
                const auto vertex=pending.front();
                pending.pop();
                for(const auto neighbour:completedLink[vertex]) {
                    const auto it=completedUnseen.find(neighbour);
                    if(it!=completedUnseen.end()) {
                        pending.push(neighbour);
                        completedUnseen.erase(it);
                    }
                }
            }
            ++completedComponents;
        }
        next.completionLinkComponents=completedComponents;
        next.completionLinkIsSingleCycle=completedComponents==1&&
            completedLink.size()>=3&&completedDegrees&&completedMultiplicity&&
            next.completionBoundaryFacesByTargetEndpoints[0]==0&&
            next.completionBoundaryFacesByTargetEndpoints[2]==0;

        std::vector<topo::CellHandle> allHandles;
        allHandles.reserve(context.cells.size());
        for(const auto& [handle,record]:context.cells) {
            (void)record;
            allHandles.push_back(handle);
        }
        sortDiagnosticCells(allHandles,context);
        std::vector<TetFacetCell> plainCells;
        plainCells.reserve(allHandles.size());
        std::vector<std::size_t> selectedIndices;
        selectedIndices.reserve(next.selected.size());
        for(std::size_t i=0;i<allHandles.size();++i) {
            const auto& record=context.cells.at(allHandles[i]);
            plainCells.push_back({record.nodes,record.sequence,
                                  record.classification});
            if(selected.count(allHandles[i])!=0)selectedIndices.push_back(i);
        }
        const TetFacetComplexView view{
            context.plainVertices,plainCells,{},context.protectedEdges,
            context.protectedFaces,context.region};
        const auto ball=auditTetFacetBall(view,selectedIndices);
        next.ballValid=ball.valid;
        next.ballReason=ball.reason;
        std::vector<std::size_t> completionIndices;
        completionIndices.reserve(completed.size());
        for(std::size_t i=0;i<allHandles.size();++i)
            if(completed.count(allHandles[i])!=0)completionIndices.push_back(i);
        const auto completionBall=auditTetFacetBall(view,completionIndices);
        next.completionBallValid=completionBall.valid;
        next.completionBallReason=completionBall.reason;
        next.completionAvailable=next.completionLinkIsSingleCycle&&
            next.completionBallValid&&!next.completionBlockedByBoundary&&
            !next.completionBlockedByProtectedFace&&
            !next.completionBlockedByRegion;
        next.complete=true;
        TetEdgeCavityDiagnosticsLock lock;
        gTetEdgeCavityDiagnostics=std::move(next);
    } catch(const std::bad_alloc&) {
        TetEdgeCavityDiagnostics failed;
        failed.enabled=true;
        failed.allocationFailed=true;
        TetEdgeCavityDiagnosticsLock lock;
        gTetEdgeCavityDiagnostics=std::move(failed);
    } catch(...) {
        TetEdgeCavityDiagnostics failed;
        failed.enabled=true;
        TetEdgeCavityDiagnosticsLock lock;
        gTetEdgeCavityDiagnostics=std::move(failed);
    }
}

bool transactionDiagnosticMatches(
    const TetEdgeRecoveryRequest& request) noexcept {
    if (!gTetEdgeTransactionDiagnosticsEnabled.load(
            std::memory_order_acquire))
        return false;
    TetEdgeCavityDiagnosticsLock lock;
    const auto requested=gTetEdgeTransactionDiagnostics.targetGlobal;
    return requested.canonicalGlobalNodes[0]==kInvalidIndex||
           requested==request.key;
}

void captureTetEdgeTransactionInput(
    const TetEdgeRecoveryRequest& request,
    const topo::TetEdgeCavitySpec& spec,
    const TetRecoveryContext& context,
    const topo::TopologyWorkspace& workspace) noexcept {
    if (!transactionDiagnosticMatches(request)) return;
    auto mix=[](std::uint64_t& digest,std::uint64_t value) noexcept {
        digest^=value+UINT64_C(0x9e3779b97f4a7c15)+(digest<<6)+(digest>>2);
    };
    std::uint64_t removedDigest=UINT64_C(0xcbf29ce484222325);
    for(const auto cell:spec.removed) {
        mix(removedDigest,cell.slot); mix(removedDigest,cell.gen);
        mix(removedDigest,cell.pool);
    }
    std::uint64_t candidateDigest=UINT64_C(0xcbf29ce484222325);
    for(const auto& cell:spec.added) {
        for(const auto node:cell.nodes) {
            const auto found=context.nodeToDense.find(node);
            mix(candidateDigest,found==context.nodeToDense.end()
                ? UINT64_MAX:found->second);
        }
        mix(candidateDigest,static_cast<std::uint64_t>(cell.cls.dim));
        mix(candidateDigest,static_cast<std::uint64_t>(cell.cls.tag));
    }
    TetEdgeCavityDiagnosticsLock lock;
    auto& diagnostic=gTetEdgeTransactionDiagnostics;
    diagnostic.transactionRequestStamp=request.stamp;
    diagnostic.transactionCurrentStamp=topo::WorkspaceStamp::of(workspace);
    diagnostic.transactionClassification=request.classification;
    diagnostic.transactionRegion=request.region;
    diagnostic.transactionRequiredConstraint=spec.requiredConstraint;
    diagnostic.transactionRemovedCount=spec.removed.size();
    diagnostic.transactionCandidateCount=spec.added.size();
    diagnostic.transactionRemovedDigest=removedDigest;
    diagnostic.transactionCandidateDigest=candidateDigest;
    diagnostic.transactionCheckpoint=
        TetEdgeTransactionDiagnosticCheckpoint::StageEntered;
}

#endif

}  // namespace

static TetEdgeRecoveryOutcome recoverTetEdgeIndexedImpl(
    const TetEdgeRecoveryRequest& request, TetFacetSetState& state,
    topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    bool assessmentOnly,
    std::vector<topo::CellHandle>* retiredWitness = nullptr,
    std::vector<topo::CellHandle>* createdWitness = nullptr,
    std::span<topo::CellHandle> frontierWitness = {},
    std::size_t* frontierWitnessUsed = nullptr,
    bool* frontierWitnessComplete = nullptr,
    TetRecoveryBindingWorkCounters* bindingWork = nullptr
#if defined(PAACS_TOPO_TESTING)
    ,TetEdgePreparedGuardForTesting preparedGuard = nullptr,
    void* preparedGuardState = nullptr
#endif
    ,const TetEdgeGeneralCavityLimits* generalLimits = nullptr,
    TetEdgeGeneralRecoveryDiagnostics* generalDiagnostics = nullptr,
    std::uint64_t maximumAdjacentExpansionAttempts = 0
    ) noexcept {
    bool topologyCommitted = false;
    try {
        if (generalDiagnostics) *generalDiagnostics = {};
        if (generalLimits && !generalLimits->positive())
            return fail(TetEdgeRecoveryStatus::Rejected,
                        TetEdgeRecoveryReason::InvalidIdentity,
                        TetEdgeRecoveryPhase::Preflight);
        const auto current = topo::WorkspaceStamp::of(workspace);
        if (!(request.stamp == current) || !(state.workspaceStamp == current))
            return fail(TetEdgeRecoveryStatus::Rejected,
                        TetEdgeRecoveryReason::StaleWorkspaceStamp,
                        TetEdgeRecoveryPhase::Preflight);
        if (!context.currentFor(workspace))
            return fail(TetEdgeRecoveryStatus::Failed,
                        TetEdgeRecoveryReason::StaleRecoveryContext,
                        TetEdgeRecoveryPhase::Preflight);
        if (!(request.artifactStamp == state.artifactStamp) ||
            request.artifactStamp.persistentId == 0)
            return fail(TetEdgeRecoveryStatus::Rejected,
                        TetEdgeRecoveryReason::StateKeyMismatch,
                        TetEdgeRecoveryPhase::Preflight);
        if (!validState(state, workspace,bindingWork))
            return fail(TetEdgeRecoveryStatus::Rejected,
                        TetEdgeRecoveryReason::StaleBinding,
                        TetEdgeRecoveryPhase::Preflight);
        if ((request.classification.dim != 1 &&
             request.classification.dim != 2) ||
            request.classification.tag < 0 || request.region.dim != 3 ||
            request.region.tag < 0)
            return fail(TetEdgeRecoveryStatus::Rejected,
                        TetEdgeRecoveryReason::InvalidClassification,
                        TetEdgeRecoveryPhase::Preflight);
        if (!workspace.nodes().isCurrent(request.vertices[0]) ||
            !workspace.nodes().isCurrent(request.vertices[1]) ||
            request.vertices[0] == request.vertices[1])
            return fail(TetEdgeRecoveryStatus::Rejected,
                        TetEdgeRecoveryReason::InvalidIdentity,
                        TetEdgeRecoveryPhase::Preflight);
        const auto p0 = workspace.nodes().get(request.vertices[0]).coord;
        const auto p1 = workspace.nodes().get(request.vertices[1]).coord;
        if (!std::isfinite(p0.x) || !std::isfinite(p0.y) ||
            !std::isfinite(p0.z) || !std::isfinite(p1.x) ||
            !std::isfinite(p1.y) || !std::isfinite(p1.z))
            return fail(TetEdgeRecoveryStatus::Rejected,
                        TetEdgeRecoveryReason::NonfiniteCoordinate,
                        TetEdgeRecoveryPhase::Preflight);
        if (p0.x == p1.x && p0.y == p1.y && p0.z == p1.z)
            return fail(TetEdgeRecoveryStatus::Rejected,
                        TetEdgeRecoveryReason::DegenerateTarget,
                        TetEdgeRecoveryPhase::Preflight);
        const auto d0 = context.nodeToDense.find(request.vertices[0]);
        const auto d1 = context.nodeToDense.find(request.vertices[1]);
        if (d0 == context.nodeToDense.end() || d1 == context.nodeToDense.end())
            return fail(TetEdgeRecoveryStatus::Rejected,
                        TetEdgeRecoveryReason::InvalidIdentity,
                        TetEdgeRecoveryPhase::Preflight);
        const auto denseKey = DenseTetEdgeKey::of(d0->second,d1->second);

        const auto oldBinding = state.edgeBindings.find(request.key);
        if (oldBinding != state.edgeBindings.end()) {
            if (!sameClass(oldBinding->second.classification,
                           request.classification))
                return fail(TetEdgeRecoveryStatus::Rejected,
                            TetEdgeRecoveryReason::WrongBindingClassification,
                            TetEdgeRecoveryPhase::Preflight);
            const auto actual = workspace.edges().find(
                workspace.nodes(), request.vertices[0], request.vertices[1]);
            if (!actual || *actual != oldBinding->second.edge)
                return fail(TetEdgeRecoveryStatus::Rejected,
                            TetEdgeRecoveryReason::WrongBindingEdge,
                            TetEdgeRecoveryPhase::Preflight);
        }
        const auto existing = workspace.edges().find(
            workspace.nodes(), request.vertices[0], request.vertices[1]);
        if (existing) {
            for (const auto& [key,binding] : state.edgeBindings) {
                if(bindingWork)(void)addBindingWork(
                    bindingWork->edgeAliasRowsVisited,1,bindingWork);
                if (!(key == request.key) && binding.edge == *existing)
                    return fail(TetEdgeRecoveryStatus::Rejected,
                                TetEdgeRecoveryReason::EdgeAlias,
                                TetEdgeRecoveryPhase::Preflight);
            }
            const auto data = workspace.edges().get(*existing);
            if(bindingWork)(void)addBindingWork(
                bindingWork->edgeDataCopyRows,1,bindingWork);
            const bool exact = sameClass(data.cls, request.classification);
            if (data.requiredConstraint && !exact)
                return fail(TetEdgeRecoveryStatus::Refused,
                            TetEdgeRecoveryReason::ProtectedConstraintConflict,
                            TetEdgeRecoveryPhase::Routing);
            if ((!exact && (data.cls.dim != -1 || data.cls.tag != -1)) ||
                (data.requiredConstraint && !exact))
                return fail(TetEdgeRecoveryStatus::Refused,
                            TetEdgeRecoveryReason::ProtectedConstraintConflict,
                            TetEdgeRecoveryPhase::Routing);
            if (assessmentOnly) {
                TetEdgeRecoveryOutcome out;
                out.status = exact && data.requiredConstraint
                    ? TetEdgeRecoveryStatus::AlreadyPresent
                    : TetEdgeRecoveryStatus::Marked;
                out.edge = *existing;
                return out;
            }
            TetFacetSetState scratch = state;
            recordStateCopy(state,bindingWork);
            scratch.edgeBindings[request.key] =
                {request.key, request.classification, true, *existing};
            if (exact && data.requiredConstraint) {
                scratch.workspaceStamp = current;
                state.swap(scratch);
                TetEdgeRecoveryOutcome out;
                out.status = TetEdgeRecoveryStatus::AlreadyPresent;
                out.edge = *existing;
                return out;
            }
            // Classification and the required bit are one metadata state.
            // The accepted target-edge mode publishes both in one noexcept
            // Phase 3 without changing cell/edge identity.
            topo::TetEdgeMetadataSpec metadata;
            metadata.stamp = current;
            metadata.target = *existing;
            metadata.targetClassification = request.classification;
            metadata.requiredConstraint = true;
            auto markTx = workspace.beginTetCavityTransaction();
            (void)markTx.stage(metadata);
            markTx.prepare();
            const auto predicted = markTx.predictedTargetEdge();
            const auto committed = markTx.commitEdge();
            topologyCommitted = true;
            if (predicted != *existing || committed.targetEdge != predicted)
                std::abort();
            const auto after = topo::WorkspaceStamp::of(workspace);
            scratch.workspaceStamp = after;
            state.swap(scratch);
            context.stamp = after;
            context.contextRevision = after.revision;
            auto row = context.edges.find(denseKey);
            if (row == context.edges.end())
                return fail(TetEdgeRecoveryStatus::Failed,
                            TetEdgeRecoveryReason::StaleRecoveryContext,
                            TetEdgeRecoveryPhase::StateCopy);
            row->second.classification = request.classification;
            row->second.required = true;
            retainProtectedEdge(context,denseKey);
            TetEdgeRecoveryOutcome out;
            out.status = TetEdgeRecoveryStatus::Marked;
            out.edge = *existing;
            return out;
        }

        std::set<topo::CellHandle> selected;
        for (const auto cell : context.querySegment(p0,p1)) {
            ++context.metrics.auditVisits;
            const auto found = context.cells.find(cell);
            if (found != context.cells.end() &&
                segmentIntersectsCellOpen(p0,p1,d0->second,d1->second,
                                          context,found->second,workspace))
                selected.insert(cell);
        }
        if (selected.empty())
            return fail(TetEdgeRecoveryStatus::Refused,
                        TetEdgeRecoveryReason::BoundaryStar,
                        TetEdgeRecoveryPhase::Discovery);

#if defined(PAACS_TOPO_TESTING)
        std::vector<topo::CellHandle> initiallySelectedForDiagnostics;
        if (gTetEdgeCavityDiagnosticsEnabled.load(std::memory_order_acquire))
            initiallySelectedForDiagnostics.assign(
                selected.begin(),selected.end());
#endif

        std::queue<topo::CellHandle> queue;
        for (const auto cell : selected) queue.push(cell);
        std::set<DenseTetEdgeKey> checkedEdges;
        while (!queue.empty()) {
            const auto cell = queue.front();
            queue.pop();
            const auto& record = context.cells.at(cell);
            for (const auto& edge : kEdges) {
                const auto key = DenseTetEdgeKey::of(
                    record.nodes[static_cast<std::size_t>(edge[0])],
                    record.nodes[static_cast<std::size_t>(edge[1])]);
                if (key == denseKey || !checkedEdges.insert(key).second) continue;
                const Vec3 a = workspace.nodes().get(
                    context.denseToNode[key.nodes[0]]).coord;
                const Vec3 b = workspace.nodes().get(
                    context.denseToNode[key.nodes[1]]).coord;
                const auto contact = tet_facet_detail::segmentContact3(p0,p1,a,b);
                if (contact == tet_facet_detail::Segment2Contact::None) continue;
                const bool sharesEndpoint = key.nodes[0] == d0->second ||
                    key.nodes[1] == d0->second || key.nodes[0] == d1->second ||
                    key.nodes[1] == d1->second;
                if (contact == tet_facet_detail::Segment2Contact::Touch &&
                    sharesEndpoint)
                    continue;
                const auto resolver = context.edges.find(key);
                if (resolver == context.edges.end())
                    return fail(TetEdgeRecoveryStatus::Refused,
                                TetEdgeRecoveryReason::NonManifoldStar,
                                TetEdgeRecoveryPhase::Discovery);
                if (resolver->second.required)
                    return fail(TetEdgeRecoveryStatus::Refused,
                                TetEdgeRecoveryReason::ProtectedConstraintConflict,
                                TetEdgeRecoveryPhase::Discovery);
                if (contact != tet_facet_detail::Segment2Contact::Proper)
                    return fail(TetEdgeRecoveryStatus::Refused,
                                TetEdgeRecoveryReason::TangentialTouch,
                                TetEdgeRecoveryPhase::Discovery);
                for (const auto starCell : resolver->second.completeStar)
                    if (selected.insert(starCell).second) queue.push(starCell);
            }
        }

#if defined(PAACS_TOPO_TESTING)
        captureTetEdgeCavityDiagnostics(
            request,denseKey,initiallySelectedForDiagnostics,selected,
            context,workspace);
#endif

        if(frontierWitnessUsed&&frontierWitnessComplete) {
            *frontierWitnessUsed=0;
            *frontierWitnessComplete=selected.size()<=frontierWitness.size();
            if(*frontierWitnessComplete)
                for(const auto cell:selected)
                    frontierWitness[(*frontierWitnessUsed)++]=cell;
        }

        std::map<DenseTetEdgeKey,std::size_t> ringEdges;
        bool starCycle = true;
        for (const auto cell : selected) {
            const auto& record = context.cells.at(cell);
            for (const auto& face : kFaces) {
                const DenseTetFaceKey faceKey = DenseTetFaceKey::of({
                    record.nodes[static_cast<std::size_t>(face[0])],
                    record.nodes[static_cast<std::size_t>(face[1])],
                    record.nodes[static_cast<std::size_t>(face[2])]});
                const auto resolver = context.faces.find(faceKey);
                if (resolver == context.faces.end())
                    return fail(TetEdgeRecoveryStatus::Refused,
                                TetEdgeRecoveryReason::NonBallClosure,
                                TetEdgeRecoveryPhase::Discovery);
                std::size_t insideOwners = 0;
                for (const auto& owner : resolver->second.owners)
                    insideOwners += selected.count(owner.cell) ? 1U : 0U;
                if (insideOwners != 1) continue;
                std::vector<DenseTetNodeId> ring;
                std::size_t targets = 0;
                for (const auto node : faceKey.nodes) {
                    if (node == d0->second || node == d1->second) ++targets;
                    else ring.push_back(node);
                }
                if (targets != 1 || ring.size() != 2) {
                    if (generalLimits) { starCycle = false; continue; }
                    return fail(TetEdgeRecoveryStatus::Refused,
                                TetEdgeRecoveryReason::NonCycleLink,
                                TetEdgeRecoveryPhase::Discovery);
                }
                ++ringEdges[DenseTetEdgeKey::of(ring[0],ring[1])];
            }
        }
        std::map<DenseTetNodeId,std::set<DenseTetNodeId>> link;
        for (const auto& [edge,count] : ringEdges) {
            if (count != 2) {
                if (generalLimits) starCycle = false;
                else
                return fail(TetEdgeRecoveryStatus::Refused,
                            TetEdgeRecoveryReason::NonCycleLink,
                            TetEdgeRecoveryPhase::Discovery);
            }
            link[edge.nodes[0]].insert(edge.nodes[1]);
            link[edge.nodes[1]].insert(edge.nodes[0]);
        }
        for (const auto& [node,neighbours] : link) {
            (void)node;
            if (neighbours.size() != 2) {
                if (generalLimits) starCycle = false;
                else
                return fail(TetEdgeRecoveryStatus::Refused,
                            TetEdgeRecoveryReason::NonCycleLink,
                            TetEdgeRecoveryPhase::Discovery);
            }
        }
        if (link.size() < 3) {
            if (generalLimits) starCycle = false;
            else
            return fail(TetEdgeRecoveryStatus::Refused,
                        TetEdgeRecoveryReason::NonCycleLink,
                        TetEdgeRecoveryPhase::Discovery);
        }

        topo::TetEdgeCavitySpec spec;
        topo::TetRefinementCavitySpec staged;
        spec.stamp = current;
        spec.target = request.vertices;
        spec.targetClassification = request.classification;
        spec.region = request.region;
        spec.removed.assign(selected.begin(),selected.end());
        if (!starCycle) {
            const auto candidate = buildBoundedGeneralEdgeCandidate(
                selected,denseKey,context,*generalLimits,
                maximumAdjacentExpansionAttempts,spec,generalDiagnostics,&staged);
            if (!candidate.succeeded()) return candidate;
        } else for (const auto& [edge,count] : ringEdges) {
            (void)count;
            topo::TetCavityCell cell;
            cell.nodes = positive(workspace,
                {request.vertices[0],request.vertices[1],
                 context.denseToNode[edge.nodes[0]],
                 context.denseToNode[edge.nodes[1]]});
            if (predicates::orient3d(
                    workspace.nodes().get(cell.nodes[0]).coord,
                    workspace.nodes().get(cell.nodes[1]).coord,
                    workspace.nodes().get(cell.nodes[2]).coord,
                    workspace.nodes().get(cell.nodes[3]).coord) !=
                Sign::Positive)
                return fail(TetEdgeRecoveryStatus::Refused,
                            TetEdgeRecoveryReason::InvalidCandidate,
                            TetEdgeRecoveryPhase::Candidate);
            cell.cls = request.region;
        spec.added.push_back(cell);
        }

        if (retiredWitness) *retiredWitness = spec.removed;

        if(frontierWitnessUsed&&frontierWitnessComplete) {
            *frontierWitnessUsed=0;
            *frontierWitnessComplete=spec.removed.size()<=frontierWitness.size();
            if(*frontierWitnessComplete)
                for(const auto cell:spec.removed)
                    frontierWitness[(*frontierWitnessUsed)++]=cell;
        }

        if (assessmentOnly) {
            TetEdgeRecoveryOutcome out;
            out.status = TetEdgeRecoveryStatus::Recovered;
            return out;
        }

        if (staged.recoverTargetEdge) {
#if defined(PAACS_TOPO_TESTING)
            // The old prepared guard observes a different preview type. A
            // caller requiring that guard must not silently bypass it here.
            if (preparedGuard)
                return fail(TetEdgeRecoveryStatus::Refused,
                    TetEdgeRecoveryReason::PendingConstraint,TetEdgeRecoveryPhase::Preparation);
#endif
            return publishSteinerEdgeCandidate(request,staged,state,workspace,
                context,bindingWork,createdWitness);
        }

#if defined(PAACS_TOPO_TESTING)
        captureTetEdgeTransactionInput(request,spec,context,workspace);
#endif
        auto transaction = workspace.beginTetCavityTransaction();
#if defined(PAACS_TOPO_TESTING)
        const topo::TetEdgeCavityPreview* preparedPreview=nullptr;
#endif
        try {
#if defined(PAACS_TOPO_TESTING)
            preparedPreview=&transaction.stage(spec);
#else
            (void)transaction.stage(spec);
#endif
            transaction.prepare();
        } catch (const std::bad_alloc&) {
            return fail(TetEdgeRecoveryStatus::Failed,
                        TetEdgeRecoveryReason::AllocationFailed,
                        TetEdgeRecoveryPhase::Preparation);
        } catch (...) {
            return fail(TetEdgeRecoveryStatus::Failed,
                        TetEdgeRecoveryReason::PreparationConflict,
                        TetEdgeRecoveryPhase::Preparation);
        }
        const auto predicted = transaction.predictedTargetEdge();
#if defined(PAACS_TOPO_TESTING)
        if(preparedGuard&&preparedGuard(
                request,current,predicted,*preparedPreview,context,
                preparedGuardState)==
                    TetEdgePreparedGuardDecisionForTesting::Defer) {
            transaction.abandon();
            auto deferred=fail(
                TetEdgeRecoveryStatus::Refused,
                TetEdgeRecoveryReason::PendingConstraint,
                TetEdgeRecoveryPhase::Preparation);
            deferred.capacityResiduePermitted=true;
            return deferred;
        }
#endif
        TetFacetSetState scratch;
        try {
            scratch = state;
            recordStateCopy(state,bindingWork);
            scratch.workspaceStamp.revision = current.revision + 1;
            scratch.edgeBindings[request.key] =
                {request.key,request.classification,true,predicted};
        } catch (const std::bad_alloc&) {
            if(bindingWork)bindingWork->complete=false;
            return fail(TetEdgeRecoveryStatus::Failed,
                        TetEdgeRecoveryReason::AllocationFailed,
                        TetEdgeRecoveryPhase::StateCopy);
        }
        const auto committed = transaction.commitEdge();
        topologyCommitted = true;
        if (createdWitness) *createdWitness = committed.resultCells;
        const auto updated = updateTetRecoveryContext(
            context,workspace,spec.removed,committed.resultCells);
        if (!updated.updated()) {
            // The topology commit is already visible.  A failed incremental
            // context update cannot leave a seemingly-current recovery
            // index available to the caller; preserve the established
            // post-commit fault tier by poisoning it explicitly.
            context.poisoned = true;
            return fail(TetEdgeRecoveryStatus::Failed,
                        updated.reason == TetRecoveryContextReason::AllocationFailed
                            ? TetEdgeRecoveryReason::AllocationFailed
                            : TetEdgeRecoveryReason::StaleRecoveryContext,
                        TetEdgeRecoveryPhase::StateCopy);
        }
        retainProtectedEdge(context,denseKey);
        scratch.workspaceStamp = topo::WorkspaceStamp::of(workspace);
        state.swap(scratch);
        TetEdgeRecoveryOutcome out;
        out.status = TetEdgeRecoveryStatus::Recovered;
        out.edge = committed.targetEdge;
        return out;
    } catch (const std::bad_alloc&) {
        if(bindingWork)bindingWork->complete=false;
        if (topologyCommitted || !context.currentFor(workspace))
            context.poisoned = true;
        return fail(TetEdgeRecoveryStatus::Failed,
                    TetEdgeRecoveryReason::AllocationFailed,
                    TetEdgeRecoveryPhase::Candidate);
    } catch (...) {
        if (topologyCommitted || !context.currentFor(workspace))
            context.poisoned = true;
        return fail(TetEdgeRecoveryStatus::Failed,
                    TetEdgeRecoveryReason::InvalidCandidate,
                    TetEdgeRecoveryPhase::Candidate);
    }
}

TetEdgeRecoveryOutcome recoverTetEdgeIndexed(
    const TetEdgeRecoveryRequest& request, TetFacetSetState& state,
    topo::TopologyWorkspace& workspace, TetRecoveryContext& context) noexcept {
    return recoverTetEdgeIndexedImpl(
        request,state,workspace,context,false);
}

TetEdgeRecoveryOutcome assessTetEdgeIndexed(
    const TetEdgeRecoveryRequest& request, const TetFacetSetState& state,
    const topo::TopologyWorkspace& workspace,
    TetRecoveryContext& context,
    TetRecoveryBindingWorkCounters* bindingWork) noexcept {
    // The assessment-only branch is mechanically mutation-free; const_cast
    // keeps the implementation and exact taxonomy single-sourced.
    return recoverTetEdgeIndexedImpl(
        request,const_cast<TetFacetSetState&>(state),
        const_cast<topo::TopologyWorkspace&>(workspace),context,true,
        nullptr,nullptr,{},nullptr,nullptr,bindingWork);
}

TetEdgeIndexedOutcome assessTetEdgeIndexedDetailed(
    const TetEdgeRecoveryRequest& request, const TetFacetSetState& state,
    const topo::TopologyWorkspace& workspace,
    TetRecoveryContext& context,
    TetRecoveryBindingWorkCounters* bindingWork) noexcept {
    TetEdgeIndexedOutcome out;
    out.recovery = recoverTetEdgeIndexedImpl(
        request,const_cast<TetFacetSetState&>(state),
        const_cast<topo::TopologyWorkspace&>(workspace),context,true,
        &out.retiredCells,nullptr,{},nullptr,nullptr,bindingWork);
    return out;
}

TetEdgeRecoveryOutcome assessTetEdgeIndexedFrontier(
    const TetEdgeRecoveryRequest& request,const TetFacetSetState& state,
    const topo::TopologyWorkspace& workspace,TetRecoveryContext& context,
    std::span<topo::CellHandle> selectedCells,
    std::size_t& selectedCellsUsed,bool& witnessComplete,
    TetRecoveryBindingWorkCounters* bindingWork) noexcept {
    return recoverTetEdgeIndexedImpl(
        request,const_cast<TetFacetSetState&>(state),
        const_cast<topo::TopologyWorkspace&>(workspace),context,true,
        nullptr,nullptr,selectedCells,&selectedCellsUsed,&witnessComplete,
        bindingWork);
}

TetEdgeIndexedOutcome recoverTetEdgeIndexedDetailed(
    const TetEdgeRecoveryRequest& request, TetFacetSetState& state,
    topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    TetRecoveryBindingWorkCounters* bindingWork) noexcept {
    TetEdgeIndexedOutcome out;
    out.recovery = recoverTetEdgeIndexedImpl(
        request,state,workspace,context,false,&out.retiredCells,
        &out.createdCells,{},nullptr,nullptr,bindingWork);
    return out;
}

TetEdgeIndexedOutcome assessTetEdgeIndexedDetailed(
    const TetEdgeRecoveryRequest& request, const TetFacetSetState& state,
    const topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    const TetEdgeGeneralCavityLimits& limits,
    TetEdgeGeneralRecoveryDiagnostics* diagnostics,
    TetRecoveryBindingWorkCounters* bindingWork,
    std::uint64_t maximumAdjacentExpansionAttempts) noexcept {
    TetEdgeIndexedOutcome out;
    out.recovery = recoverTetEdgeIndexedImpl(
        request,const_cast<TetFacetSetState&>(state),
        const_cast<topo::TopologyWorkspace&>(workspace),context,true,
        &out.retiredCells,nullptr,{},nullptr,nullptr,bindingWork
#if defined(PAACS_TOPO_TESTING)
        ,nullptr,nullptr
#endif
        ,&limits,diagnostics,maximumAdjacentExpansionAttempts);
    return out;
}

TetEdgeIndexedOutcome recoverTetEdgeIndexedDetailed(
    const TetEdgeRecoveryRequest& request, TetFacetSetState& state,
    topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    const TetEdgeGeneralCavityLimits& limits,
    TetEdgeGeneralRecoveryDiagnostics* diagnostics,
    TetRecoveryBindingWorkCounters* bindingWork,
    std::uint64_t maximumAdjacentExpansionAttempts) noexcept {
    TetEdgeIndexedOutcome out;
    out.recovery = recoverTetEdgeIndexedImpl(
        request,state,workspace,context,false,&out.retiredCells,
        &out.createdCells,{},nullptr,nullptr,bindingWork
#if defined(PAACS_TOPO_TESTING)
        ,nullptr,nullptr
#endif
        ,&limits,diagnostics,maximumAdjacentExpansionAttempts);
    return out;
}

TetEdgeRecoveryOutcome assessTetEdgeIndexedFrontier(
    const TetEdgeRecoveryRequest& request,const TetFacetSetState& state,
    const topo::TopologyWorkspace& workspace,TetRecoveryContext& context,
    std::span<topo::CellHandle> selectedCells,
    std::size_t& selectedCellsUsed,bool& witnessComplete,
    const TetEdgeGeneralCavityLimits& limits,
    TetEdgeGeneralRecoveryDiagnostics* diagnostics,
    TetRecoveryBindingWorkCounters* bindingWork,
    std::uint64_t maximumAdjacentExpansionAttempts) noexcept {
    return recoverTetEdgeIndexedImpl(
        request,const_cast<TetFacetSetState&>(state),
        const_cast<topo::TopologyWorkspace&>(workspace),context,true,
        nullptr,nullptr,selectedCells,&selectedCellsUsed,&witnessComplete,
        bindingWork
#if defined(PAACS_TOPO_TESTING)
        ,nullptr,nullptr
#endif
        ,&limits,diagnostics,maximumAdjacentExpansionAttempts);
}

#if defined(PAACS_TOPO_TESTING)
TetEdgeIndexedOutcome recoverTetEdgeIndexedPreparedGuardedForTesting(
    const TetEdgeRecoveryRequest& request,TetFacetSetState& state,
    topo::TopologyWorkspace& workspace,TetRecoveryContext& context,
    TetEdgePreparedGuardForTesting guard,void* guardState,
    TetRecoveryBindingWorkCounters* bindingWork) noexcept {
    TetEdgeIndexedOutcome out;
    out.recovery=recoverTetEdgeIndexedImpl(
        request,state,workspace,context,false,&out.retiredCells,
        &out.createdCells,{},nullptr,nullptr,bindingWork,guard,guardState);
    return out;
}

TetEdgeIndexedOutcome recoverTetEdgeIndexedPreparedGuardedForTesting(
    const TetEdgeRecoveryRequest& request,TetFacetSetState& state,
    topo::TopologyWorkspace& workspace,TetRecoveryContext& context,
    TetEdgePreparedGuardForTesting guard,void* guardState,
    const TetEdgeGeneralCavityLimits& limits,
    TetEdgeGeneralRecoveryDiagnostics* diagnostics,
    TetRecoveryBindingWorkCounters* bindingWork,
    std::uint64_t maximumAdjacentExpansionAttempts) noexcept {
    TetEdgeIndexedOutcome out;
    out.recovery=recoverTetEdgeIndexedImpl(
        request,state,workspace,context,false,&out.retiredCells,
        &out.createdCells,{},nullptr,nullptr,bindingWork,guard,guardState,
        &limits,diagnostics,maximumAdjacentExpansionAttempts);
    return out;
}
#endif

}  // namespace paacs::bnd
