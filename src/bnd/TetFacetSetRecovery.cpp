// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/bnd/TetFacetSetRecovery.hpp"

#include "paacs/predicates/Predicates.hpp"

#include <algorithm>
#include <bit>
#include <cstdlib>
#include <functional>
#include <map>
#include <new>
#include <queue>
#include <set>

namespace paacs::bnd {
namespace {

#if defined(PAACS_TOPO_TESTING)
thread_local bool gMixedPresentCommitZoneActive = false;

struct MixedPresentCommitZoneScope {
    MixedPresentCommitZoneScope() noexcept {
        if (gMixedPresentCommitZoneActive)
            std::abort();
        gMixedPresentCommitZoneActive = true;
    }
    ~MixedPresentCommitZoneScope() noexcept { gMixedPresentCommitZoneActive = false; }
};
#endif

bool sameClass(const mesh::Classification &a, const mesh::Classification &b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}

bool addCounter(std::uint64_t &target, std::uint64_t value) noexcept {
    if (value > (std::numeric_limits<std::uint64_t>::max)() - target)
        return false;
    target += value;
    return true;
}

bool addSize(std::size_t a, std::size_t b, std::size_t &out) noexcept {
    if (b > (std::numeric_limits<std::size_t>::max)() - a)
        return false;
    out = a + b;
    return true;
}

bool sizeAsCounter(std::size_t value, std::uint64_t &out) noexcept {
    if (value > (std::numeric_limits<std::uint64_t>::max)())
        return false;
    out = static_cast<std::uint64_t>(value);
    return true;
}

bool addMixedMetric(std::uint64_t &target, std::uint64_t value,
                    TetFacetSetMixedPresentBatchSink *sink) noexcept {
    if (!sink)
        return true;
    if (!addCounter(target, value)) {
        sink->metrics.complete = false;
        return false;
    }
    return true;
}

bool addMixedBytes(std::uint64_t &target, std::size_t count, std::size_t width,
                   TetFacetSetMixedPresentBatchSink *sink) noexcept {
    if (!sink)
        return true;
    if (count != 0 && width > (std::numeric_limits<std::uint64_t>::max)() / count) {
        sink->metrics.complete = false;
        return false;
    }
    return addMixedMetric(
        target, static_cast<std::uint64_t>(count) * static_cast<std::uint64_t>(width), sink);
}

template <class Status>
bool countMixedStatus(std::array<std::uint64_t, 6> &counts, Status status,
                      TetFacetSetMixedPresentBatchSink *sink) noexcept {
    if (!sink)
        return true;
    const auto index = static_cast<std::size_t>(status);
    if (index >= counts.size()) {
        sink->metrics.complete = false;
        return false;
    }
    return addMixedMetric(counts[index], 1, sink);
}

bool addAssessmentCounters(TetFacetGeometryCounters &target,
                           const TetFacetGeometryCounters &value) noexcept {
    return addCounter(target.cellsScanned, value.cellsScanned) &&
           addCounter(target.crossedEdgeStars, value.crossedEdgeStars) &&
           addCounter(target.closureInsertions, value.closureInsertions) &&
           addCounter(target.fixedPointPasses, value.fixedPointPasses) &&
           addCounter(target.apexAttempts, value.apexAttempts) &&
           addCounter(target.tetTetIntersectionTests, value.tetTetIntersectionTests);
}

bool addMixedGeometryCounters(TetFacetGeometryCounters &target,
                              const TetFacetGeometryCounters &value,
                              TetFacetSetMixedPresentBatchSink *sink) noexcept {
    if (!sink)
        return true;
    if (addAssessmentCounters(target, value))
        return true;
    sink->metrics.complete = false;
    return false;
}

TetFacetRecoveryOutcome
facetFailure(TetFacetRecoveryStatus status, TetFacetRecoveryReason reason,
             TetFacetRecoveryPhase phase,
             TetFacetGeometryReason geometry = TetFacetGeometryReason::None) {
    TetFacetRecoveryOutcome out;
    out.status = status;
    out.reason = reason;
    out.phase = phase;
    out.geometryReason = geometry;
    return out;
}

TetEdgeRecoveryOutcome edgeFailure(TetEdgeRecoveryStatus status, TetEdgeRecoveryReason reason,
                                   TetEdgeRecoveryPhase phase) {
    TetEdgeRecoveryOutcome out;
    out.status = status;
    out.reason = reason;
    out.phase = phase;
    return out;
}

constexpr std::array<std::array<int, 2>, 6> kEdges{
    {{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
constexpr std::array<std::array<int, 3>, 4> kFaces{
    {{{1, 2, 3}}, {{0, 3, 2}}, {{0, 1, 3}}, {{0, 2, 1}}}};

std::array<Vec3, 3> targetPoints(const TetFacetIndexedRequest &request,
                                 const topo::TopologyWorkspace &workspace) {
    return {workspace.nodes().get(request.vertices[0]).coord,
            workspace.nodes().get(request.vertices[1]).coord,
            workspace.nodes().get(request.vertices[2]).coord};
}

bool faceStraddles(const std::array<Vec3, 3> &target, const std::array<DenseTetNodeId, 3> &face,
                   const TetRecoveryContext &context, const topo::TopologyWorkspace &workspace) {
    bool positive = false, negative = false;
    for (const auto id : face) {
        const auto sign = predicates::orient3d(
            target[0], target[1], target[2], workspace.nodes().get(context.denseToNode[id]).coord);
        positive = positive || sign == predicates::Sign::Positive;
        negative = negative || sign == predicates::Sign::Negative;
    }
    return positive && negative;
}

TetFacetCell plainCell(const CurrentTetRecord &source) {
    TetFacetCell cell;
    cell.vertices = source.nodes;
    cell.sequence = source.sequence;
    cell.classification = source.classification;
    return cell;
}

struct PlainSubset {
    std::vector<TetFacetCell> cells;
    std::vector<topo::CellHandle> handles;
    std::map<topo::CellHandle, std::size_t> indexByHandle;
};

PlainSubset makeSubset(const std::set<topo::CellHandle> &handles,
                       const TetRecoveryContext &context) {
    PlainSubset result;
    result.handles.assign(handles.begin(), handles.end());
    std::sort(result.handles.begin(), result.handles.end(), [&](const auto &a, const auto &b) {
        return context.cells.at(a).sequence < context.cells.at(b).sequence;
    });
    result.cells.reserve(result.handles.size());
    for (std::size_t i = 0; i < result.handles.size(); ++i) {
        result.indexByHandle.emplace(result.handles[i], i);
        result.cells.push_back(plainCell(context.cells.at(result.handles[i])));
    }
    return result;
}

bool addBindingWork(std::uint64_t &target, std::uint64_t amount,
                    TetRecoveryBindingWorkCounters *work) noexcept {
    if (!work)
        return true;
    if (amount > (std::numeric_limits<std::uint64_t>::max)() - target) {
        work->complete = false;
        return false;
    }
    target += amount;
    return true;
}

void recordFaceDataCopy(const topo::FacePool::FaceData &data,
                        TetRecoveryBindingWorkCounters *work) noexcept {
    if (!work)
        return;
    (void)addBindingWork(work->faceDataCopyRows, 1, work);
    if (data.incident.size() > (std::numeric_limits<std::uint64_t>::max)()) {
        work->complete = false;
        return;
    }
    (void)addBindingWork(work->faceIncidenceRowsCopied,
                         static_cast<std::uint64_t>(data.incident.size()), work);
}

void recordStateCopy(const TetFacetSetState &state, TetRecoveryBindingWorkCounters *work) noexcept {
    if (!work)
        return;
    (void)addBindingWork(work->stateCopies, 1, work);
    const std::array<std::size_t, 4> rows{state.edgeBindings.size(), state.bindings.size(),
                                          state.canonicalEdgeSchedule.size(),
                                          state.canonicalSchedule.size()};
    for (const auto count : rows) {
        if (count > (std::numeric_limits<std::uint64_t>::max)()) {
            work->complete = false;
            continue;
        }
        (void)addBindingWork(work->stateCopyRows, static_cast<std::uint64_t>(count), work);
    }
}

std::array<topo::NodeHandle, 3> canonicalHandles(const topo::TopologyWorkspace &workspace,
                                                 std::array<topo::NodeHandle, 3> nodes) {
    std::sort(nodes.begin(), nodes.end(), [&](const auto a, const auto b) {
        return workspace.nodes().seqOf(a) < workspace.nodes().seqOf(b);
    });
    return nodes;
}

struct FacetBindingReceiptRow {
    std::array<topo::NodeHandle, 3> canonicalVertices{};
    VolumeFacetKey key{};
    topo::FaceHandle recoveredFace{};
};

enum class FacetBindingReceiptVerdict : std::uint8_t { Unbuilt, Valid, Invalid };

// Scheduler-private receipt. It never escapes one recoverTetFacetSetIndexed
// call and is explicitly invalidated whenever that scheduler publishes a
// facet-binding state change. Workspace/artifact stamps prevent reuse across
// an indexed topology delta.
struct TetFacetBindingValidationReceipt {
    const TetFacetSetState *state = nullptr;
    topo::WorkspaceStamp workspaceStamp{};
    algo::ArtifactStamp artifactStamp{};
    std::size_t bindingCount = 0;
    std::uint64_t bindingStateRevision = 0;
    FacetBindingReceiptVerdict verdict = FacetBindingReceiptVerdict::Unbuilt;
    std::vector<FacetBindingReceiptRow> rowsByTriple;

    void invalidate() noexcept {
        state = nullptr;
        workspaceStamp = {};
        artifactStamp = {};
        bindingCount = 0;
        bindingStateRevision = 0;
        verdict = FacetBindingReceiptVerdict::Unbuilt;
        rowsByTriple.clear();
    }
};

bool validateFacetBindingRow(const VolumeFacetKey &key, const PersistentTetFacetBinding &binding,
                             const topo::TopologyWorkspace &workspace, FacetBindingReceiptRow &row,
                             TetRecoveryBindingWorkCounters *work) {
    if (work)
        (void)addBindingWork(work->facetBindingRowsVisited, 1, work);
    if (!(key == binding.key))
        return false;
    for (const auto vertex : binding.canonicalVertices)
        if (!workspace.nodes().isCurrent(vertex))
            return false;
    if (canonicalHandles(workspace, binding.canonicalVertices) != binding.canonicalVertices)
        return false;
    if (!workspace.faces().isCurrent(binding.recoveredFace))
        return false;
    const auto data = workspace.faces().get(binding.recoveredFace);
    recordFaceDataCopy(data, work);
    if (data.cornerCount != 3)
        return false;
    std::array<topo::NodeHandle, 3> actual{data.corners[0], data.corners[1], data.corners[2]};
    for (const auto vertex : actual)
        if (!workspace.nodes().isCurrent(vertex))
            return false;
    if (canonicalHandles(workspace, actual) != binding.canonicalVertices)
        return false;
    if (!sameClass(data.cls, binding.classification))
        return false;
    row = {binding.canonicalVertices, key, binding.recoveredFace};
    return true;
}

bool receiptMatches(const TetFacetBindingValidationReceipt &receipt, const TetFacetSetState &state,
                    const topo::TopologyWorkspace &workspace,
                    std::uint64_t bindingStateRevision) noexcept {
    return receipt.state == &state &&
           receipt.workspaceStamp == topo::WorkspaceStamp::of(workspace) &&
           receipt.artifactStamp == state.artifactStamp &&
           receipt.bindingCount == state.bindings.size() &&
           receipt.bindingStateRevision == bindingStateRevision &&
           receipt.verdict != FacetBindingReceiptVerdict::Unbuilt;
}

bool validBindings(const TetFacetSetState &state, const topo::TopologyWorkspace &workspace,
                   TetRecoveryBindingWorkCounters *work,
                   TetFacetBindingValidationReceipt *receipt = nullptr,
                   std::uint64_t bindingStateRevision = 0) {
    if (work)
        (void)addBindingWork(work->facetBindingAuditCalls, 1, work);
    if (receipt && receiptMatches(*receipt, state, workspace, bindingStateRevision)) {
        if (work)
            (void)addBindingWork(work->facetBindingReceiptReuses, 1, work);
        return receipt->verdict == FacetBindingReceiptVerdict::Valid;
    }

    if (!receipt) {
        std::set<topo::FaceHandle> handles;
        for (auto current = state.bindings.begin(); current != state.bindings.end(); ++current) {
            const auto &key = current->first;
            const auto &binding = current->second;
            // Preserve the default path's single-set allocation profile, but
            // reject corrupt row identity/currentness/payload before the set
            // can allocate for this row.  Triple injectivity remains the
            // allocation-free nested comparison below.
            FacetBindingReceiptRow row;
            if (!validateFacetBindingRow(key, binding, workspace, row, work))
                return false;
            if (!handles.insert(binding.recoveredFace).second)
                return false;
            for (auto prior = state.bindings.begin(); prior != current; ++prior) {
                if (work)
                    (void)addBindingWork(work->facetBindingDuplicateRowsCompared, 1, work);
                if (prior->second.canonicalVertices == binding.canonicalVertices)
                    return false;
            }
        }
        return true;
    }

    receipt->invalidate();
    receipt->state = &state;
    receipt->workspaceStamp = topo::WorkspaceStamp::of(workspace);
    receipt->artifactStamp = state.artifactStamp;
    receipt->bindingCount = state.bindings.size();
    receipt->bindingStateRevision = bindingStateRevision;
    if (work)
        (void)addBindingWork(work->facetBindingReceiptBuilds, 1, work);
    if (receipt->rowsByTriple.capacity() < state.bindings.size()) {
        if (work)
            (void)addBindingWork(work->facetBindingReceiptReserveCalls, 1, work);
        receipt->rowsByTriple.reserve(state.bindings.size());
    }
    for (const auto &[key, binding] : state.bindings) {
        receipt->rowsByTriple.push_back({binding.canonicalVertices, key, binding.recoveredFace});
        if (work)
            (void)addBindingWork(work->facetBindingReceiptRowsCopied, 1, work);
    }
    std::sort(receipt->rowsByTriple.begin(), receipt->rowsByTriple.end(),
              [](const auto &a, const auto &b) {
                  return std::tie(a.recoveredFace, a.key) < std::tie(b.recoveredFace, b.key);
              });
    if (std::adjacent_find(receipt->rowsByTriple.begin(), receipt->rowsByTriple.end(),
                           [](const auto &a, const auto &b) {
                               return a.recoveredFace == b.recoveredFace;
                           }) != receipt->rowsByTriple.end()) {
        receipt->verdict = FacetBindingReceiptVerdict::Invalid;
        return false;
    }
    std::sort(receipt->rowsByTriple.begin(), receipt->rowsByTriple.end(),
              [](const auto &a, const auto &b) {
                  return std::tie(a.canonicalVertices, a.key) <
                         std::tie(b.canonicalVertices, b.key);
              });
    if (std::adjacent_find(receipt->rowsByTriple.begin(), receipt->rowsByTriple.end(),
                           [](const auto &a, const auto &b) {
                               return a.canonicalVertices == b.canonicalVertices;
                           }) != receipt->rowsByTriple.end()) {
        receipt->verdict = FacetBindingReceiptVerdict::Invalid;
        return false;
    }
    for (const auto &[key, binding] : state.bindings) {
        FacetBindingReceiptRow validated;
        if (!validateFacetBindingRow(key, binding, workspace, validated, work)) {
            receipt->verdict = FacetBindingReceiptVerdict::Invalid;
            return false;
        }
    }
    receipt->verdict = FacetBindingReceiptVerdict::Valid;
    return true;
}

bool validBindingRequest(const TetFacetSetState &state, const TetFacetIndexedRequest &request,
                         const std::array<topo::NodeHandle, 3> &requestTriple,
                         const topo::TopologyWorkspace &workspace,
                         TetRecoveryBindingWorkCounters *work,
                         const TetFacetBindingValidationReceipt *receipt = nullptr,
                         std::uint64_t bindingStateRevision = 0) noexcept {
    const auto exact = state.bindings.find(request.key);
    if (exact != state.bindings.end() &&
        (exact->second.canonicalVertices != requestTriple ||
         !sameClass(exact->second.classification, request.classification)))
        return false;
    if (receipt && receiptMatches(*receipt, state, workspace, bindingStateRevision) &&
        receipt->verdict == FacetBindingReceiptVerdict::Valid) {
        if (work)
            (void)addBindingWork(work->facetBindingReceiptLookups, 1, work);
        const auto found = std::lower_bound(
            receipt->rowsByTriple.begin(), receipt->rowsByTriple.end(), requestTriple,
            [](const auto &row, const auto &triple) { return row.canonicalVertices < triple; });
        return found == receipt->rowsByTriple.end() || found->canonicalVertices != requestTriple ||
               found->key == request.key;
    }
    for (const auto &[key, binding] : state.bindings) {
        if (work)
            (void)addBindingWork(work->facetBindingRequestRowsVisited, 1, work);
        if (!(key == request.key) && binding.canonicalVertices == requestTriple)
            return false;
    }
    return true;
}

void retainProtectedFace(TetRecoveryContext &context, const DenseTetFaceKey &key,
                         topo::FaceHandle handle, mesh::Classification classification) {
    const TetFacetFaceKey faceKey{key.nodes};
    const auto found = std::find_if(context.protectedFaces.begin(), context.protectedFaces.end(),
                                    [&](const auto &value) { return value.key == faceKey; });
    const TetFacetProtectedFace value{faceKey, static_cast<std::uint64_t>(handle.slot) + 1U,
                                      classification};
    if (found == context.protectedFaces.end())
        context.protectedFaces.push_back(value);
    else
        *found = value;
    std::sort(context.protectedFaces.begin(), context.protectedFaces.end(),
              [](const auto &a, const auto &b) { return a.key < b.key; });
}

TetAabb candidateBounds(const TetFacetCell &cell, const std::vector<TetFacetVertex> &vertices) {
    TetAabb bounds{vertices[cell.vertices[0]].position, vertices[cell.vertices[0]].position};
    for (std::size_t i = 1; i < 4; ++i) {
        const Vec3 p = vertices[cell.vertices[i]].position;
        bounds.lo.x = std::min(bounds.lo.x, p.x);
        bounds.lo.y = std::min(bounds.lo.y, p.y);
        bounds.lo.z = std::min(bounds.lo.z, p.z);
        bounds.hi.x = std::max(bounds.hi.x, p.x);
        bounds.hi.y = std::max(bounds.hi.y, p.y);
        bounds.hi.z = std::max(bounds.hi.z, p.z);
    }
    return bounds;
}

bool facetReady(const TetFacetRecoveryOutcome &outcome) noexcept {
    return outcome.status == TetFacetRecoveryStatus::AlreadyPresent ||
           outcome.status == TetFacetRecoveryStatus::Marked ||
           outcome.status == TetFacetRecoveryStatus::Recovered;
}

bool mapFacetRequest(const VolumeFacet &facet,
                     const std::map<Index, topo::NodeHandle> &localByGlobal,
                     const TetFacetSetState &state, const topo::TopologyWorkspace &workspace,
                     mesh::Classification temporaryRegion, TetFacetIndexedRequest &request) {
    request.stamp = topo::WorkspaceStamp::of(workspace);
    request.artifactStamp = state.artifactStamp;
    request.key = facet.key;
    request.classification = facet.classification;
    request.region = temporaryRegion;
    for (std::size_t i = 0; i < 3; ++i) {
        const auto found = localByGlobal.find(facet.directedGlobalNodes[i]);
        if (found == localByGlobal.end())
            return false;
        request.vertices[i] = found->second;
    }
    return true;
}

TetAabb facetBounds(const TetFacetIndexedRequest &request,
                    const topo::TopologyWorkspace &workspace) {
    const Vec3 first = workspace.nodes().get(request.vertices[0]).coord;
    TetAabb bounds{first, first};
    for (std::size_t i = 1; i < 3; ++i) {
        const Vec3 point = workspace.nodes().get(request.vertices[i]).coord;
        bounds.lo.x = std::min(bounds.lo.x, point.x);
        bounds.lo.y = std::min(bounds.lo.y, point.y);
        bounds.lo.z = std::min(bounds.lo.z, point.z);
        bounds.hi.x = std::max(bounds.hi.x, point.x);
        bounds.hi.y = std::max(bounds.hi.y, point.y);
        bounds.hi.z = std::max(bounds.hi.z, point.z);
    }
    return bounds;
}

bool selectedCavityTouchesTarget(std::span<const topo::CellHandle> selected,
                                 const TetFacetIndexedRequest &target,
                                 const TetRecoveryContext &context) {
    std::array<DenseTetNodeId, 3> dense{};
    for (std::size_t i = 0; i < 3; ++i) {
        const auto found = context.nodeToDense.find(target.vertices[i]);
        if (found == context.nodeToDense.end())
            return true;
        dense[i] = found->second;
    }
    for (const auto handle : selected) {
        const auto found = context.cells.find(handle);
        if (found == context.cells.end())
            return true;
        const TetFacetCell cell = plainCell(found->second);
        const TetFacetComplexView view{
            context.plainVertices,   std::span<const TetFacetCell>(&cell, 1),
            TetFacetTriangle{dense}, context.protectedEdges,
            context.protectedFaces,  target.region};
        const auto contact = classifyTetFacetContact(view, cell);
        if (contact != TetFacetContact::Disjoint &&
            contact != TetFacetContact::DeclaredBoundaryTouch)
            return true;
    }
    return false;
}

struct ScheduledFacetAssessment {
    TetFacetIndexedRequest request{};
    TetFacetIndexedOutcome detail{};
    TetAabb targetBounds{};
    std::vector<std::size_t> blockers;
    bool current = false;
};

bool mapEdgeRequest(const VolumeArtifactEdge &edge,
                    const std::map<Index, topo::NodeHandle> &localByGlobal,
                    const TetFacetSetState &state, const topo::TopologyWorkspace &workspace,
                    mesh::Classification temporaryRegion, TetEdgeRecoveryRequest &request) {
    const auto a = localByGlobal.find(edge.key.canonicalGlobalNodes[0]);
    const auto b = localByGlobal.find(edge.key.canonicalGlobalNodes[1]);
    if (a == localByGlobal.end() || b == localByGlobal.end())
        return false;
    request.stamp = topo::WorkspaceStamp::of(workspace);
    request.artifactStamp = state.artifactStamp;
    request.key = edge.key;
    request.vertices = {a->second, b->second};
    request.classification = edge.classification;
    request.region = temporaryRegion;
    return true;
}

bool exactPresentEdgeStillProtected(const TetEdgeRecoveryRequest &request,
                                    const TetEdgeRecoveryOutcome &outcome,
                                    const topo::TopologyWorkspace &workspace,
                                    const TetRecoveryContext &context) noexcept {
    if (outcome.status != TetEdgeRecoveryStatus::AlreadyPresent || !outcome.edge.isSpecified() ||
        !workspace.edges().isCurrent(outcome.edge))
        return false;
    const auto actual =
        workspace.edges().find(workspace.nodes(), request.vertices[0], request.vertices[1]);
    if (!actual || *actual != outcome.edge)
        return false;
    const auto data = workspace.edges().get(outcome.edge);
    if (!data.requiredConstraint || !sameClass(data.cls, request.classification))
        return false;
    const auto d0 = context.nodeToDense.find(request.vertices[0]);
    const auto d1 = context.nodeToDense.find(request.vertices[1]);
    if (d0 == context.nodeToDense.end() || d1 == context.nodeToDense.end())
        return false;
    const auto dense = DenseTetEdgeKey::of(d0->second, d1->second);
    const auto row = context.edges.find(dense);
    if (row == context.edges.end() || row->second.handle != outcome.edge || !row->second.required ||
        !sameClass(row->second.classification, request.classification))
        return false;
    const auto protectedKey = canonicalTetFacetEdge(dense.nodes[0], dense.nodes[1]);
    return std::binary_search(context.protectedEdges.begin(), context.protectedEdges.end(),
                              protectedKey);
}

bool exactPresentFacetStillProtected(const TetFacetIndexedRequest &request,
                                     const TetFacetRecoveryOutcome &outcome,
                                     const topo::TopologyWorkspace &workspace,
                                     const TetRecoveryContext &context) noexcept {
    if (outcome.status != TetFacetRecoveryStatus::AlreadyPresent || !outcome.face.isSpecified() ||
        !workspace.faces().isCurrent(outcome.face))
        return false;
    std::array<DenseTetNodeId, 3> dense{};
    for (std::size_t i = 0; i < dense.size(); ++i) {
        const auto found = context.nodeToDense.find(request.vertices[i]);
        if (found == context.nodeToDense.end())
            return false;
        dense[i] = found->second;
    }
    const auto key = DenseTetFaceKey::of(dense);
    const auto row = context.faces.find(key);
    if (row == context.faces.end() || row->second.handle != outcome.face ||
        row->second.owners.size() != 2 ||
        row->second.owners[0].sense * row->second.owners[1].sense != -1 ||
        !sameClass(row->second.classification, request.classification))
        return false;
    // `context` is stamp-bound to this workspace and is refreshed after
    // every committed recovery.  Its face row is therefore the exact,
    // allocation-free classification authority here.  FacePool::get()
    // returns FaceData by value and copies the incidence vector; using it in
    // this qualification predicate both repeats the same check and turns an
    // otherwise read-only assessment into a fallible allocation.
    for (std::size_t i = 0; i < dense.size(); ++i) {
        const auto edge =
            context.edges.find(DenseTetEdgeKey::of(dense[i], dense[(i + 1) % dense.size()]));
        if (edge == context.edges.end() || !edge->second.required)
            return false;
    }
    const TetFacetFaceKey protectedKey{key.nodes};
    const auto found =
        std::lower_bound(context.protectedFaces.begin(), context.protectedFaces.end(), protectedKey,
                         [](const TetFacetProtectedFace &value, const TetFacetFaceKey &wanted) {
                             return value.key < wanted;
                         });
    return found != context.protectedFaces.end() && found->key == protectedKey &&
           found->identity == static_cast<std::uint64_t>(outcome.face.slot) + 1U &&
           sameClass(found->classification, request.classification);
}

bool exactPresentFacetCandidate(const TetFacetIndexedRequest &request,
                                const topo::TopologyWorkspace &workspace,
                                const TetRecoveryContext &context,
                                TetFacetRecoveryOutcome &outcome) noexcept {
    std::array<DenseTetNodeId, 3> dense{};
    for (std::size_t i = 0; i < dense.size(); ++i) {
        const auto found = context.nodeToDense.find(request.vertices[i]);
        if (found == context.nodeToDense.end())
            return false;
        dense[i] = found->second;
    }
    const auto found = context.faces.find(DenseTetFaceKey::of(dense));
    if (found == context.faces.end())
        return false;
    outcome = {};
    outcome.status = TetFacetRecoveryStatus::AlreadyPresent;
    outcome.face = found->second.handle;
    return exactPresentFacetStillProtected(request, outcome, workspace, context);
}

bool segmentIntersectsTetOpen(const Vec3 &a, const Vec3 &b, const CurrentTetRecord &record,
                              const TetRecoveryContext &context) {
    std::array<Vec3, 4> tet{};
    for (std::size_t i = 0; i < 4; ++i)
        tet[i] = context.plainVertices[record.nodes[i]].position;
    if (tet_facet_detail::pointStrictlyInsideTet(a, tet) ||
        tet_facet_detail::pointStrictlyInsideTet(b, tet))
        return true;
    for (const auto &face : kFaces) {
        const Vec3 &p = tet[static_cast<std::size_t>(face[0])];
        const Vec3 &q = tet[static_cast<std::size_t>(face[1])];
        const Vec3 &r = tet[static_cast<std::size_t>(face[2])];
        if (tet_facet_detail::segmentCrossesOpenTriangle(a, b, p, q, r))
            return true;
    }
    return false;
}

bool selectedCavityTouchesEdge(std::span<const topo::CellHandle> selected,
                               const TetEdgeRecoveryRequest &target,
                               const TetRecoveryContext &context,
                               const topo::TopologyWorkspace &workspace) {
    const auto d0 = context.nodeToDense.find(target.vertices[0]);
    const auto d1 = context.nodeToDense.find(target.vertices[1]);
    if (d0 == context.nodeToDense.end() || d1 == context.nodeToDense.end())
        return true;
    const auto edge = context.edges.find(DenseTetEdgeKey::of(d0->second, d1->second));
    if (edge != context.edges.end())
        for (const auto handle : selected)
            if (std::find(edge->second.completeStar.begin(), edge->second.completeStar.end(),
                          handle) != edge->second.completeStar.end())
                return true;
    const Vec3 a = workspace.nodes().get(target.vertices[0]).coord;
    const Vec3 b = workspace.nodes().get(target.vertices[1]).coord;
    for (const auto handle : selected) {
        const auto cell = context.cells.find(handle);
        if (cell == context.cells.end() || segmentIntersectsTetOpen(a, b, cell->second, context))
            return true;
    }
    return false;
}

struct ScheduledEdgeAssessment {
    TetEdgeRecoveryRequest request{};
    TetEdgeIndexedOutcome detail{};
    std::vector<std::size_t> blockers;
    std::size_t frontierSupportOffset = 0;
    std::size_t frontierSupportCount = 0;
    bool current = false;
};

bool blockerGraphValid(const std::vector<std::vector<std::size_t>> &graph,
                       const std::vector<bool> &pending) noexcept {
    if (graph.size() != pending.size())
        return false;
    for (std::size_t i = 0; i < graph.size(); ++i) {
        if (!pending[i] && !graph[i].empty())
            return false;
        for (const auto blocker : graph[i])
            if (blocker >= pending.size() || !pending[blocker])
                return false;
    }
    return true;
}

std::vector<std::size_t> canonicalConflictCycle(const std::vector<std::vector<std::size_t>> &graph,
                                                const std::vector<bool> &pending,
                                                std::uint64_t &visits) {
    const std::size_t count = graph.size();
    std::vector<int> index(count, -1), low(count, -1);
    std::vector<std::size_t> stack;
    std::vector<bool> onStack(count, false);
    std::vector<std::vector<std::size_t>> cycles;
    int nextIndex = 0;
    std::function<void(std::size_t)> visit = [&](std::size_t v) {
        ++visits;
        index[v] = low[v] = nextIndex++;
        stack.push_back(v);
        onStack[v] = true;
        for (const auto w : graph[v]) {
            if (w >= count || !pending[w])
                continue;
            if (index[w] < 0) {
                visit(w);
                low[v] = std::min(low[v], low[w]);
            } else if (onStack[w]) {
                low[v] = std::min(low[v], index[w]);
            }
        }
        if (low[v] != index[v])
            return;
        std::vector<std::size_t> component;
        while (!stack.empty()) {
            const auto w = stack.back();
            stack.pop_back();
            onStack[w] = false;
            component.push_back(w);
            if (w == v)
                break;
        }
        std::sort(component.begin(), component.end());
        const bool self =
            component.size() == 1 &&
            std::find(graph[component.front()].begin(), graph[component.front()].end(),
                      component.front()) != graph[component.front()].end();
        if (component.size() > 1 || self)
            cycles.push_back(std::move(component));
    };
    for (std::size_t i = 0; i < count; ++i)
        if (pending[i] && index[i] < 0)
            visit(i);
    std::sort(cycles.begin(), cycles.end());
    return cycles.empty() ? std::vector<std::size_t>{} : std::move(cycles.front());
}

void beginFacetConflictCycleCapture(TetFacetConflictCycleSink &sink,
                                    topo::WorkspaceStamp workspaceStamp,
                                    algo::ArtifactStamp artifactStamp,
                                    std::size_t membersRequired) noexcept {
    sink.reset();
    sink.workspaceStamp = workspaceStamp;
    sink.artifactStamp = artifactStamp;
    sink.membersRequired = membersRequired;
    sink.captured = true;
    if (membersRequired > sink.members.size())
        sink.complete = false;
}

void appendFacetConflictCycleMember(TetFacetConflictCycleSink &sink, VolumeFacetKey key) noexcept {
    if (sink.membersUsed < sink.members.size())
        sink.members[sink.membersUsed++] = key;
    else
        sink.complete = false;
}

void appendFacetConflictCycleDependency(TetFacetConflictCycleSink &sink,
                                        std::array<VolumeFacetKey, 2> dependency) noexcept {
    if (sink.dependenciesRequired == (std::numeric_limits<std::size_t>::max)()) {
        sink.complete = false;
    } else {
        ++sink.dependenciesRequired;
    }
    if (sink.dependenciesUsed < sink.dependencies.size())
        sink.dependencies[sink.dependenciesUsed++] = dependency;
    else
        sink.complete = false;
}

void finishFacetConflictCycleCapture(TetFacetConflictCycleSink &sink) noexcept {
    if (sink.membersUsed > 1)
        std::sort(sink.members.begin(), sink.members.begin() + sink.membersUsed);
    if (sink.dependenciesUsed > 1) {
        std::sort(sink.dependencies.begin(), sink.dependencies.begin() + sink.dependenciesUsed);
        const auto uniqueEnd = std::unique(sink.dependencies.begin(),
                                           sink.dependencies.begin() + sink.dependenciesUsed);
        sink.dependenciesUsed = static_cast<std::size_t>(uniqueEnd - sink.dependencies.begin());
    }
    if (sink.complete)
        sink.dependenciesRequired = sink.dependenciesUsed;
}

} // namespace

#if defined(PAACS_TOPO_TESTING)
bool mixedPresentCommitZoneActiveForTesting() noexcept {
    return gMixedPresentCommitZoneActive;
}

namespace {
thread_local TetFacetConflictCycleDiagnostics gFacetConflictDiagnostics{};
}

void resetTetFacetConflictCycleDiagnosticsForTesting(bool enabled) noexcept {
    gFacetConflictDiagnostics = {};
    gFacetConflictDiagnostics.enabled = enabled;
}

const TetFacetConflictCycleDiagnostics &tetFacetConflictCycleDiagnosticsForTesting() noexcept {
    return gFacetConflictDiagnostics;
}

void fillTetFacetConflictCycleSinkForTesting(
    std::span<const VolumeFacetKey> members,
    std::span<const std::array<VolumeFacetKey, 2>> dependencies,
    topo::WorkspaceStamp workspaceStamp, algo::ArtifactStamp artifactStamp,
    TetFacetConflictCycleSink &sink) noexcept {
    beginFacetConflictCycleCapture(sink, workspaceStamp, artifactStamp, members.size());
    for (const auto member : members)
        appendFacetConflictCycleMember(sink, member);
    for (const auto dependency : dependencies)
        appendFacetConflictCycleDependency(sink, dependency);
    finishFacetConflictCycleCapture(sink);
}
#endif

std::vector<std::size_t>
canonicalTetRecoveryConflictCycleForTesting(const std::vector<std::vector<std::size_t>> &graph,
                                            const std::vector<bool> &pending,
                                            std::uint64_t &visits) {
    return canonicalConflictCycle(graph, pending, visits);
}

bool tetRecoveryBlockerGraphValidForTesting(const std::vector<std::vector<std::size_t>> &graph,
                                            const std::vector<bool> &pending) noexcept {
    return blockerGraphValid(graph, pending);
}

static TetFacetIndexedOutcome
recoverTetFacetIndexedImpl(const TetFacetIndexedRequest &request, TetFacetSetState &state,
                           topo::TopologyWorkspace &workspace, TetRecoveryContext &context,
                           bool assessmentOnly,
                           TetRecoveryBindingWorkCounters *bindingWork = nullptr,
                           TetFacetBindingValidationReceipt *bindingReceipt = nullptr,
                           std::uint64_t bindingStateRevision = 0) noexcept {
    TetFacetIndexedOutcome result;
    bool topologyCommitted = false;
    try {
        const auto current = topo::WorkspaceStamp::of(workspace);
        if (!(request.stamp == current) || !(state.workspaceStamp == current)) {
            result.recovery = facetFailure(TetFacetRecoveryStatus::Rejected,
                                           TetFacetRecoveryReason::StaleWorkspaceStamp,
                                           TetFacetRecoveryPhase::Preflight);
            return result;
        }
        if (!context.currentFor(workspace)) {
            result.recovery = facetFailure(TetFacetRecoveryStatus::Failed,
                                           TetFacetRecoveryReason::PreparationConflict,
                                           TetFacetRecoveryPhase::Preflight);
            return result;
        }
        if (!(request.artifactStamp == state.artifactStamp) ||
            request.artifactStamp.persistentId == 0 ||
            !validBindings(state, workspace, bindingWork, bindingReceipt, bindingStateRevision)) {
            result.recovery = facetFailure(TetFacetRecoveryStatus::Rejected,
                                           TetFacetRecoveryReason::StateKeyMismatch,
                                           TetFacetRecoveryPhase::Preflight);
            return result;
        }
        if (request.classification.dim != 2 || request.classification.tag < 0 ||
            request.region.dim != 3 || request.region.tag < 0) {
            result.recovery = facetFailure(TetFacetRecoveryStatus::Rejected,
                                           TetFacetRecoveryReason::InvalidClassification,
                                           TetFacetRecoveryPhase::Preflight);
            return result;
        }
        std::array<DenseTetNodeId, 3> dense{};
        for (std::size_t i = 0; i < 3; ++i) {
            if (!workspace.nodes().isCurrent(request.vertices[i])) {
                result.recovery = facetFailure(TetFacetRecoveryStatus::Rejected,
                                               TetFacetRecoveryReason::InvalidIdentity,
                                               TetFacetRecoveryPhase::Preflight);
                return result;
            }
            const auto found = context.nodeToDense.find(request.vertices[i]);
            if (found == context.nodeToDense.end()) {
                result.recovery = facetFailure(TetFacetRecoveryStatus::Rejected,
                                               TetFacetRecoveryReason::InvalidIdentity,
                                               TetFacetRecoveryPhase::Preflight);
                return result;
            }
            dense[i] = found->second;
        }
        auto sortedDense = dense;
        std::sort(sortedDense.begin(), sortedDense.end());
        if (std::adjacent_find(sortedDense.begin(), sortedDense.end()) != sortedDense.end()) {
            result.recovery = facetFailure(TetFacetRecoveryStatus::Rejected,
                                           TetFacetRecoveryReason::InvalidIdentity,
                                           TetFacetRecoveryPhase::Preflight);
            return result;
        }
        const auto requestTriple = canonicalHandles(workspace, request.vertices);
        if (!validBindingRequest(state, request, requestTriple, workspace, bindingWork,
                                 bindingReceipt, bindingStateRevision)) {
            result.recovery = facetFailure(TetFacetRecoveryStatus::Rejected,
                                           TetFacetRecoveryReason::StateKeyMismatch,
                                           TetFacetRecoveryPhase::Preflight);
            return result;
        }
        const auto target = targetPoints(request, workspace);
        if (predicates::orient2d({target[0].x, target[0].y}, {target[1].x, target[1].y},
                                 {target[2].x, target[2].y}) == predicates::Sign::Zero &&
            predicates::orient2d({target[0].y, target[0].z}, {target[1].y, target[1].z},
                                 {target[2].y, target[2].z}) == predicates::Sign::Zero &&
            predicates::orient2d({target[0].z, target[0].x}, {target[1].z, target[1].x},
                                 {target[2].z, target[2].x}) == predicates::Sign::Zero) {
            result.recovery = facetFailure(TetFacetRecoveryStatus::Rejected,
                                           TetFacetRecoveryReason::DegenerateTarget,
                                           TetFacetRecoveryPhase::Preflight);
            return result;
        }
        for (std::size_t i = 0; i < 3; ++i) {
            const auto edge = context.edges.find(DenseTetEdgeKey::of(dense[i], dense[(i + 1) % 3]));
            if (edge == context.edges.end() || !edge->second.required) {
                result.recovery = facetFailure(TetFacetRecoveryStatus::Rejected,
                                               TetFacetRecoveryReason::MissingProtectedTargetEdge,
                                               TetFacetRecoveryPhase::Preflight);
                return result;
            }
        }

        const DenseTetFaceKey targetKey = DenseTetFaceKey::of(dense);
        const auto currentFace = context.faces.find(targetKey);
        if (currentFace != context.faces.end()) {
            if (currentFace->second.owners.size() != 2 ||
                currentFace->second.owners[0].sense * currentFace->second.owners[1].sense != -1) {
                result.recovery =
                    facetFailure(TetFacetRecoveryStatus::Refused,
                                 currentFace->second.owners.size() == 1
                                     ? TetFacetRecoveryReason::OneSidedTargetUnsupported
                                     : TetFacetRecoveryReason::InvalidTargetIncidence,
                                 TetFacetRecoveryPhase::Routing);
                return result;
            }
            const auto data = workspace.faces().get(currentFace->second.handle);
            recordFaceDataCopy(data, bindingWork);
            if (data.cls.dim == 2 && data.cls.tag >= 0 &&
                !sameClass(data.cls, request.classification)) {
                result.recovery = facetFailure(TetFacetRecoveryStatus::Refused,
                                               TetFacetRecoveryReason::ProtectedFacetConflict,
                                               TetFacetRecoveryPhase::Routing);
                return result;
            }
            if (assessmentOnly) {
                result.recovery.status = sameClass(data.cls, request.classification)
                                             ? TetFacetRecoveryStatus::AlreadyPresent
                                             : TetFacetRecoveryStatus::Marked;
                result.recovery.face = currentFace->second.handle;
                return result;
            }
            TetFacetSetState scratch = state;
            recordStateCopy(state, bindingWork);
            const auto canonical = canonicalHandles(workspace, request.vertices);
            scratch.bindings[request.key] = {request.key, canonical, request.classification,
                                             currentFace->second.handle};
            if (sameClass(data.cls, request.classification)) {
                scratch.workspaceStamp = current;
                state.swap(scratch);
                result.recovery.status = TetFacetRecoveryStatus::AlreadyPresent;
                result.recovery.face = currentFace->second.handle;
                result.contextUpdated = true;
                return result;
            }
            auto transaction = workspace.beginTransaction();
            transaction.setFaceClassification(currentFace->second.handle, request.classification);
            (void)transaction.commit();
            topologyCommitted = true;
            const auto after = topo::WorkspaceStamp::of(workspace);
            scratch.workspaceStamp = after;
            state.swap(scratch);
            context.stamp = after;
            context.contextRevision = after.revision;
            context.faces.at(targetKey).classification = request.classification;
            retainProtectedFace(context, targetKey, currentFace->second.handle,
                                request.classification);
            result.recovery.status = TetFacetRecoveryStatus::Marked;
            result.recovery.face = currentFace->second.handle;
            result.contextUpdated = true;
            return result;
        }

        std::set<topo::CellHandle> local;
        std::queue<topo::CellHandle> queue;
        for (const auto cell : context.queryTriangle(target))
            if (local.insert(cell).second)
                queue.push(cell);
        while (!queue.empty()) {
            const auto handle = queue.front();
            queue.pop();
            const auto &cell = context.cells.at(handle);
            for (const auto &edge : kEdges) {
                const auto key = DenseTetEdgeKey::of(cell.nodes[static_cast<std::size_t>(edge[0])],
                                                     cell.nodes[static_cast<std::size_t>(edge[1])]);
                const Vec3 a = workspace.nodes().get(context.denseToNode[key.nodes[0]]).coord;
                const Vec3 b = workspace.nodes().get(context.denseToNode[key.nodes[1]]).coord;
                if (!tet_facet_detail::segmentCrossesOpenTriangle(a, b, target[0], target[1],
                                                                  target[2]))
                    continue;
                const auto row = context.edges.find(key);
                if (row == context.edges.end())
                    continue;
                for (const auto owner : row->second.completeStar)
                    if (local.insert(owner).second)
                        queue.push(owner);
            }
            for (const auto &face : kFaces) {
                const std::array<DenseTetNodeId, 3> nodes{
                    cell.nodes[static_cast<std::size_t>(face[0])],
                    cell.nodes[static_cast<std::size_t>(face[1])],
                    cell.nodes[static_cast<std::size_t>(face[2])]};
                if (!faceStraddles(target, nodes, context, workspace))
                    continue;
                const auto row = context.faces.find(DenseTetFaceKey::of(nodes));
                if (row == context.faces.end())
                    continue;
                for (const auto &owner : row->second.owners)
                    if (local.insert(owner.cell).second)
                        queue.push(owner.cell);
            }
        }
        if (local.empty()) {
            result.recovery = facetFailure(
                TetFacetRecoveryStatus::Refused, TetFacetRecoveryReason::UnsupportedDegeneracy,
                TetFacetRecoveryPhase::Discovery, TetFacetGeometryReason::NoProperCut);
            return result;
        }
        const auto &vertices = context.plainVertices;
        const auto &protectedEdges = context.protectedEdges;
        const auto &protectedFaces = context.protectedFaces;
        auto subset = makeSubset(local, context);
        const TetFacetComplexView view{vertices,       subset.cells,   TetFacetTriangle{dense},
                                       protectedEdges, protectedFaces, request.region};
        auto plan = planTetFacetCavityIndexed(view);
        if (!plan.planned()) {
            const bool allocation = plan.reason == TetFacetGeometryReason::AllocationFailed;
            result.recovery =
                facetFailure(allocation || plan.status == TetFacetPlanStatus::Failed
                                 ? TetFacetRecoveryStatus::Failed
                                 : TetFacetRecoveryStatus::Refused,
                             allocation ? TetFacetRecoveryReason::AllocationFailed
                                        : TetFacetRecoveryReason::UnsupportedDegeneracy,
                             allocation || plan.status == TetFacetPlanStatus::Failed
                                 ? TetFacetRecoveryPhase::Candidate
                                 : TetFacetRecoveryPhase::Discovery,
                             plan.reason);
            result.recovery.counters = plan.plan.counters;
            return result;
        }
        auto candidate = reconstructTetFacetCandidate(view, plan.plan);
        if (!candidate.built()) {
            const bool allocation = candidate.reason == TetFacetGeometryReason::AllocationFailed;
            result.recovery =
                facetFailure(allocation || candidate.status == TetFacetPlanStatus::Failed
                                 ? TetFacetRecoveryStatus::Failed
                                 : TetFacetRecoveryStatus::Refused,
                             allocation ? TetFacetRecoveryReason::AllocationFailed
                                        : TetFacetRecoveryReason::InvalidCandidate,
                             TetFacetRecoveryPhase::Candidate, candidate.reason);
            return result;
        }

        std::set<topo::CellHandle> validationHandles = local;
        std::vector<topo::CellHandle> removedHandles;
        for (const auto localIndex : plan.plan.selectedCells) {
            const auto handle = subset.handles[localIndex];
            removedHandles.push_back(handle);
            for (const auto node : context.cells.at(handle).nodes)
                validationHandles.insert(context.nodeStars[node].begin(),
                                         context.nodeStars[node].end());
        }
        for (const auto &cell : candidate.candidate.cells) {
            const auto bounds = candidateBounds(cell, vertices);
            const auto hits = context.queryBounds(bounds);
            validationHandles.insert(hits.begin(), hits.end());
        }
        auto validationSubset = makeSubset(validationHandles, context);
        std::vector<std::size_t> validationRemoved;
        for (const auto handle : removedHandles)
            validationRemoved.push_back(validationSubset.indexByHandle.at(handle));
        const TetFacetComplexView validationView{
            vertices,       validationSubset.cells, TetFacetTriangle{dense},
            protectedEdges, protectedFaces,         request.region};
        TetFacetGeometryCounters counters = candidate.candidate.counters;
        const auto audit = validateTetFacetCandidateIndexed(validationView, validationRemoved,
                                                            candidate.candidate.cells, counters);
        if (!audit.valid) {
            result.recovery = facetFailure(TetFacetRecoveryStatus::Failed,
                                           audit.reason == TetFacetGeometryReason::AllocationFailed
                                               ? TetFacetRecoveryReason::AllocationFailed
                                               : TetFacetRecoveryReason::InvalidCandidate,
                                           TetFacetRecoveryPhase::Candidate, audit.reason);
            result.recovery.counters = counters;
            return result;
        }

        // This exact selected-cell set is the scheduler's dependency witness.
        // It is populated before the assessment-only return and owns no
        // candidate topology across a workspace revision.
        result.retiredCells = removedHandles;
        if (assessmentOnly) {
            result.recovery.status = TetFacetRecoveryStatus::Recovered;
            result.recovery.counters = counters;
            return result;
        }

        topo::TetCavitySpec spec;
        spec.stamp = current;
        spec.removed = removedHandles;
        spec.target = request.vertices;
        spec.targetClassification = request.classification;
        spec.region = request.region;
        for (const auto &source : candidate.candidate.cells) {
            topo::TetCavityCell cell;
            for (std::size_t i = 0; i < 4; ++i)
                cell.nodes[i] = context.denseToNode[source.vertices[i]];
            cell.cls = source.classification;
            spec.added.push_back(cell);
        }
        const auto canonicalTarget = canonicalHandles(workspace, request.vertices);
        auto transaction = workspace.beginTetCavityTransaction();
        try {
            (void)transaction.stage(spec);
            transaction.prepare();
        } catch (const std::bad_alloc &) {
            result.recovery = facetFailure(TetFacetRecoveryStatus::Failed,
                                           TetFacetRecoveryReason::AllocationFailed,
                                           TetFacetRecoveryPhase::Preparation);
            return result;
        } catch (...) {
            result.recovery = facetFailure(TetFacetRecoveryStatus::Failed,
                                           TetFacetRecoveryReason::PreparationConflict,
                                           TetFacetRecoveryPhase::Preparation);
            return result;
        }
        const auto predicted = transaction.predictedTargetFace();
        TetFacetSetState scratch;
        try {
            scratch = state;
            recordStateCopy(state, bindingWork);
            scratch.workspaceStamp.revision = current.revision + 1;
            scratch.bindings[request.key] = {request.key, canonicalTarget, request.classification,
                                             predicted};
        } catch (const std::bad_alloc &) {
            if (bindingWork)
                bindingWork->complete = false;
            result.recovery = facetFailure(TetFacetRecoveryStatus::Failed,
                                           TetFacetRecoveryReason::AllocationFailed,
                                           TetFacetRecoveryPhase::StateCopy);
            return result;
        }
        const auto committed = transaction.commit();
        topologyCommitted = true;
        result.retiredCells = removedHandles;
        result.createdCells = committed.resultCells;
        const auto update =
            updateTetRecoveryContext(context, workspace, result.retiredCells, result.createdCells);
        if (!update.updated()) {
            // The cavity transaction has committed.  Do not expose a stale
            // but unpoisoned indexed context when its incremental refresh
            // cannot be retained under allocation failure.
            context.poisoned = true;
            result.recovery =
                facetFailure(TetFacetRecoveryStatus::Failed,
                             update.reason == TetRecoveryContextReason::AllocationFailed
                                 ? TetFacetRecoveryReason::AllocationFailed
                                 : TetFacetRecoveryReason::PreparationConflict,
                             TetFacetRecoveryPhase::StateCopy);
            return result;
        }
        retainProtectedFace(context, targetKey, committed.targetFace, request.classification);
        scratch.workspaceStamp = topo::WorkspaceStamp::of(workspace);
        state.swap(scratch);
        result.recovery.status = TetFacetRecoveryStatus::Recovered;
        result.recovery.face = committed.targetFace;
        result.recovery.counters = counters;
        result.contextUpdated = true;
        return result;
    } catch (const std::bad_alloc &) {
        if (bindingWork)
            bindingWork->complete = false;
        if (topologyCommitted || !context.currentFor(workspace))
            context.poisoned = true;
        result.recovery =
            facetFailure(TetFacetRecoveryStatus::Failed, TetFacetRecoveryReason::AllocationFailed,
                         TetFacetRecoveryPhase::Candidate);
        return result;
    } catch (...) {
        if (topologyCommitted || !context.currentFor(workspace))
            context.poisoned = true;
        result.recovery =
            facetFailure(TetFacetRecoveryStatus::Failed, TetFacetRecoveryReason::InvalidCandidate,
                         TetFacetRecoveryPhase::Candidate);
        return result;
    }
}

TetFacetIndexedOutcome
recoverTetFacetIndexed(const TetFacetIndexedRequest &request, TetFacetSetState &state,
                       topo::TopologyWorkspace &workspace, TetRecoveryContext &context,
                       TetRecoveryBindingWorkCounters *bindingWork) noexcept {
    return recoverTetFacetIndexedImpl(request, state, workspace, context, false, bindingWork);
}

TetFacetRecoveryOutcome
assessTetFacetIndexed(const TetFacetIndexedRequest &request, const TetFacetSetState &state,
                      const topo::TopologyWorkspace &workspace, TetRecoveryContext &context,
                      TetRecoveryBindingWorkCounters *bindingWork) noexcept {
    return recoverTetFacetIndexedImpl(request, const_cast<TetFacetSetState &>(state),
                                      const_cast<topo::TopologyWorkspace &>(workspace), context,
                                      true, bindingWork)
        .recovery;
}

TetFacetIndexedOutcome
assessTetFacetIndexedDetailed(const TetFacetIndexedRequest &request, const TetFacetSetState &state,
                              const topo::TopologyWorkspace &workspace, TetRecoveryContext &context,
                              TetRecoveryBindingWorkCounters *bindingWork) noexcept {
    return recoverTetFacetIndexedImpl(request, const_cast<TetFacetSetState &>(state),
                                      const_cast<topo::TopologyWorkspace &>(workspace), context,
                                      true, bindingWork);
}

static TetFacetSetRecoveryOutcome recoverTetFacetSetIndexedImpl(
    const ClosedTri3Artifact &artifact, const std::map<Index, topo::NodeHandle> &localByGlobal,
    mesh::Classification temporaryRegion, TetFacetSetState &state,
    topo::TopologyWorkspace &workspace, TetRecoveryContext &context,
    TetFacetSetFrontierSink *frontier, TetFacetConflictCycleSink *conflictCycle,
    TetFacetSetMixedPresentBatchSink *mixedPresent,
    const TetEdgeGeneralRecoveryPolicy* generalRecovery) noexcept {
    TetFacetSetRecoveryOutcome out;
    if (generalRecovery && !generalRecovery->valid()) {
        out.status = TetFacetSetRecoveryStatus::Rejected;
        out.reason = TetFacetSetRecoveryReason::InternalError;
        return out;
    }
    const bool useGeneral = generalRecovery && generalRecovery->enabled;
    std::uint64_t remainingGeneralWork = useGeneral
        ? generalRecovery->maximumTotalWorkUnits : 0;
    const auto generalLimits = [&]() {
        auto limits = generalRecovery->limits;
        limits.maxWorkUnits = (std::min)(limits.maxWorkUnits,remainingGeneralWork);
        return limits;
    };
    const auto accountGeneral = [&](const TetEdgeGeneralRecoveryDiagnostics& evidence) {
        remainingGeneralWork -= (std::min)(remainingGeneralWork,evidence.metrics.workUnits);
    };
    const auto budgetFailure = []() {
        return edgeFailure(TetEdgeRecoveryStatus::Refused,
                           TetEdgeRecoveryReason::ResourceBudgetReached,
                           TetEdgeRecoveryPhase::Candidate);
    };
    const auto assessEdgeDetailed = [&](
        const TetEdgeRecoveryRequest& request, const TetFacetSetState& currentState,
        const topo::TopologyWorkspace& currentWorkspace, TetRecoveryContext& currentContext,
        TetRecoveryBindingWorkCounters* bindingWork) {
        if (!useGeneral) return assessTetEdgeIndexedDetailed(
            request,currentState,currentWorkspace,currentContext,bindingWork);
        if (!remainingGeneralWork) {
            TetEdgeIndexedOutcome failed; failed.recovery = budgetFailure(); return failed;
        }
        TetEdgeGeneralRecoveryDiagnostics evidence;
        auto result = assessTetEdgeIndexedDetailed(
            request,currentState,currentWorkspace,currentContext,
            generalLimits(),&evidence,bindingWork,
            generalRecovery->maximumAdjacentExpansionAttempts);
        accountGeneral(evidence);
        return result;
    };
    const auto assessEdgeFrontier = [&](
        const TetEdgeRecoveryRequest& request, const TetFacetSetState& currentState,
        const topo::TopologyWorkspace& currentWorkspace, TetRecoveryContext& currentContext,
        std::span<topo::CellHandle> selectedCells, std::size_t& used, bool& complete,
        TetRecoveryBindingWorkCounters* bindingWork) {
        if (!useGeneral) return assessTetEdgeIndexedFrontier(
            request,currentState,currentWorkspace,currentContext,selectedCells,used,complete,bindingWork);
        if (!remainingGeneralWork) { used = 0; complete = true; return budgetFailure(); }
        TetEdgeGeneralRecoveryDiagnostics evidence;
        auto result = assessTetEdgeIndexedFrontier(
            request,currentState,currentWorkspace,currentContext,selectedCells,used,complete,
            generalLimits(),&evidence,bindingWork,
            generalRecovery->maximumAdjacentExpansionAttempts);
        accountGeneral(evidence);
        return result;
    };
    const auto recoverEdgeDetailed = [&](
        const TetEdgeRecoveryRequest& request, TetFacetSetState& currentState,
        topo::TopologyWorkspace& currentWorkspace, TetRecoveryContext& currentContext,
        TetRecoveryBindingWorkCounters* bindingWork) {
        if (!useGeneral) return recoverTetEdgeIndexedDetailed(
            request,currentState,currentWorkspace,currentContext,bindingWork);
        if (!remainingGeneralWork) {
            TetEdgeIndexedOutcome failed; failed.recovery = budgetFailure(); return failed;
        }
        TetEdgeGeneralRecoveryDiagnostics evidence;
        auto result = recoverTetEdgeIndexedDetailed(
            request,currentState,currentWorkspace,currentContext,
            generalLimits(),&evidence,bindingWork,
            generalRecovery->maximumAdjacentExpansionAttempts);
        accountGeneral(evidence);
        return result;
    };
    bool reportEnvelopePrepared = false;
    bool mixedSchedulesInstalled = false;
    bool mixedStateAdvanced = false;
    TetFacetBindingValidationReceipt mixedFacetBindingReceipt;
    std::uint64_t mixedFacetBindingStateRevision = 1;
    std::vector<VolumeArtifactEdgeKey> previousEdgeSchedule;
    std::vector<VolumeFacetKey> previousFacetSchedule;
    const auto rollbackUncommittedMixedSchedules = [&]() noexcept {
        if (!mixedPresent || !mixedSchedulesInstalled || mixedStateAdvanced)
            return;
        state.canonicalEdgeSchedule.swap(previousEdgeSchedule);
        state.canonicalSchedule.swap(previousFacetSchedule);
        mixedSchedulesInstalled = false;
    };
    const auto itemizeUnprocessed = [&]() noexcept {
        out.unattemptedEdges.clear();
        out.unattemptedFacets.clear();
        for (const auto &key : state.canonicalEdgeSchedule) {
            const bool found =
                std::any_of(out.processed.begin(), out.processed.end(), [&](const auto &row) {
                    return row.kind == TetFacetSetWorkKind::Edge && row.edgeKey == key;
                });
            if (!found)
                out.unattemptedEdges.push_back(key);
        }
        for (const auto &key : state.canonicalSchedule) {
            const bool found =
                std::any_of(out.processed.begin(), out.processed.end(), [&](const auto &row) {
                    return row.kind == TetFacetSetWorkKind::Facet && row.facetKey == key;
                });
            if (!found)
                out.unattemptedFacets.push_back(key);
        }
    };
    try {
        if (frontier)
            frontier->reset();
        if (conflictCycle)
            conflictCycle->reset();
        if (mixedPresent)
            mixedPresent->reset();
        const auto syncBindingWorkCompleteness = [&]() noexcept {
            if (mixedPresent && (!mixedPresent->metrics.assessmentBindingWork.complete ||
                                 !mixedPresent->metrics.recoveryBindingWork.complete))
                mixedPresent->metrics.complete = false;
        };
        const auto assessFacet = [&](const TetFacetIndexedRequest &requestRow) {
            if (!mixedPresent)
                return assessTetFacetIndexedDetailed(requestRow, state, workspace, context);
            auto result = recoverTetFacetIndexedImpl(requestRow, state, workspace, context, true,
                                                     &mixedPresent->metrics.assessmentBindingWork,
                                                     &mixedFacetBindingReceipt,
                                                     mixedFacetBindingStateRevision);
            syncBindingWorkCompleteness();
            return result;
        };
        if (!(state.artifactStamp == artifact.surfaceSetStamp) ||
            !(state.workspaceStamp == topo::WorkspaceStamp::of(workspace)) ||
            !context.currentFor(workspace) || temporaryRegion.dim != 3 || temporaryRegion.tag < 0) {
            out.status = TetFacetSetRecoveryStatus::Rejected;
            return out;
        }
        // The mixed/private route must validate its complete retained facet
        // binding state at the immutable entry stamp, before installing or
        // reserving either canonical schedule.  Every later assessment at
        // this unchanged binding/workspace revision reuses this receipt.
        // The accepted null/default route retains its existing allocation
        // profile and continues to validate through its ordinary assessment.
        if (mixedPresent) {
            if (!validBindings(state, workspace, &mixedPresent->metrics.assessmentBindingWork,
                               &mixedFacetBindingReceipt, mixedFacetBindingStateRevision)) {
                syncBindingWorkCompleteness();
                out.status = TetFacetSetRecoveryStatus::Rejected;
                out.reason = TetFacetSetRecoveryReason::None;
                return out;
            }
            syncBindingWorkCompleteness();
        }
#if defined(PAACS_TOPO_TESTING)
        if (mixedPresent && mixedPresent->maximumBatchEntitiesForTesting == 0) {
            out.status = TetFacetSetRecoveryStatus::Failed;
            out.reason = TetFacetSetRecoveryReason::MixedPresentCapacityReached;
            return out;
        }
#endif
        std::vector<const VolumeArtifactEdge *> edges;
        edges.reserve(artifact.edges.size());
        for (const auto &edge : artifact.edges)
            edges.push_back(&edge);
        std::sort(edges.begin(), edges.end(),
                  [](const auto *a, const auto *b) { return a->key < b->key; });
        std::vector<const VolumeFacet *> facets;
        facets.reserve(artifact.facets.size());
        for (const auto &facet : artifact.facets)
            facets.push_back(&facet);
        std::sort(facets.begin(), facets.end(), [](const auto *a, const auto *b) {
            const int amin = a->owners.empty() ? -1 : a->owners.front().regionTag;
            const int bmin = b->owners.empty() ? -1 : b->owners.front().regionTag;
            const int amax = a->owners.empty() ? -1 : a->owners.back().regionTag;
            const int bmax = b->owners.empty() ? -1 : b->owners.back().regionTag;
            const bool ai = a->owners.size() == 2;
            const bool bi = b->owners.size() == 2;
            return std::tie(amin, ai, amax, a->key, a->directedGlobalNodes) <
                   std::tie(bmin, bi, bmax, b->key, b->directedGlobalNodes);
        });
        if (facets.size() > (std::numeric_limits<std::size_t>::max)() - edges.size()) {
            out.status = TetFacetSetRecoveryStatus::Failed;
            out.reason = TetFacetSetRecoveryReason::MixedPresentCapacityReached;
            return out;
        }
        out.processed.reserve(edges.size() + facets.size());
        out.unattemptedEdges.reserve(edges.size());
        out.unattemptedFacets.reserve(facets.size());
        if (mixedPresent) {
            std::vector<VolumeArtifactEdgeKey> preparedEdges;
            std::vector<VolumeFacetKey> preparedFacets;
            preparedEdges.reserve(edges.size());
            preparedFacets.reserve(facets.size());
            for (const auto *edge : edges)
                preparedEdges.push_back(edge->key);
            for (const auto *facet : facets)
                preparedFacets.push_back(facet->key);
            state.canonicalEdgeSchedule.swap(previousEdgeSchedule);
            state.canonicalSchedule.swap(previousFacetSchedule);
            state.canonicalEdgeSchedule.swap(preparedEdges);
            state.canonicalSchedule.swap(preparedFacets);
            mixedSchedulesInstalled = true;
        } else {
            state.canonicalEdgeSchedule.clear();
            state.canonicalSchedule.clear();
            state.canonicalEdgeSchedule.reserve(edges.size());
            state.canonicalSchedule.reserve(facets.size());
            for (const auto *edge : edges)
                state.canonicalEdgeSchedule.push_back(edge->key);
            for (const auto *facet : facets)
                state.canonicalSchedule.push_back(facet->key);
        }
        reportEnvelopePrepared = true;

        // A construction whose complete prescribed artifact is already
        // present must not pay one growing persistent-state copy per entity.
        // Assess against one immutable workspace/state stamp, prepare every
        // binding in one scratch state, then publish with the existing
        // nothrow swap. Any nontrivial route falls through unchanged.
        std::vector<TetEdgeRecoveryRequest> presentEdgeRequests;
        std::vector<TetEdgeRecoveryOutcome> presentEdgeOutcomes;
        std::vector<TetFacetIndexedRequest> presentFacetRequests;
        std::vector<TetFacetRecoveryOutcome> presentFacetOutcomes;
        std::vector<TetEdgeIndexedOutcome> mixedEdgeDetails;
        std::vector<TetFacetIndexedOutcome> mixedFacetDetails;
        std::vector<bool> mixedFacetPreclassified;
        std::vector<std::size_t> mixedEdgeSupportOffsets;
        std::vector<std::size_t> mixedEdgeSupportCounts;
        presentEdgeRequests.reserve(edges.size());
        presentEdgeOutcomes.reserve(edges.size());
        presentFacetRequests.reserve(facets.size());
        presentFacetOutcomes.reserve(facets.size());
        if (mixedPresent) {
            mixedEdgeDetails.reserve(edges.size());
            mixedFacetDetails.resize(facets.size());
            mixedFacetPreclassified.resize(facets.size(), false);
            mixedEdgeSupportOffsets.resize(edges.size());
            mixedEdgeSupportCounts.resize(edges.size());
        }
        bool allAlreadyPresent = true;
        bool mixedFrontierCapacityReached = false;
        for (const auto *edge : edges) {
            TetEdgeRecoveryRequest edgeRequest;
            if (!mapEdgeRequest(*edge, localByGlobal, state, workspace, temporaryRegion,
                                edgeRequest)) {
                allAlreadyPresent = false;
                if (!mixedPresent)
                    break;
                TetEdgeIndexedOutcome detail;
                detail.recovery = edgeFailure(TetEdgeRecoveryStatus::Rejected,
                                              TetEdgeRecoveryReason::InvalidIdentity,
                                              TetEdgeRecoveryPhase::Preflight);
                mixedEdgeDetails.push_back(detail);
                presentEdgeRequests.push_back(edgeRequest);
                presentEdgeOutcomes.push_back(detail.recovery);
                ++out.schedulerMetrics.routeAssessments;
                if (!countMixedStatus(mixedPresent->metrics.edgeAssessmentEventStatusCounts,
                                      detail.recovery.status, mixedPresent) ||
                    !addMixedGeometryCounters(mixedPresent->metrics.assessmentCounters,
                                              detail.recovery.counters, mixedPresent))
                    throw std::logic_error("mixed edge status overflow");
                continue;
            }
            TetEdgeIndexedOutcome detail;
            if (mixedPresent && frontier) {
                const auto index = mixedEdgeDetails.size();
                mixedEdgeSupportOffsets[index] = frontier->supportCellsUsed;
                std::size_t used = 0;
                bool complete = true;
                detail.recovery = assessEdgeFrontier(
                    edgeRequest, state, workspace, context,
                    frontier->supportCells.subspan(frontier->supportCellsUsed), used, complete,
                    mixedPresent ? &mixedPresent->metrics.assessmentBindingWork : nullptr);
                syncBindingWorkCompleteness();
                if (!complete ||
                    used > frontier->supportCells.size() - frontier->supportCellsUsed) {
                    frontier->complete = false;
                    mixedFrontierCapacityReached = true;
                    used = 0;
                }
                mixedEdgeSupportCounts[index] = used;
                frontier->supportCellsUsed += used;
                if (!addAssessmentCounters(frontier->assessmentCounters, detail.recovery.counters))
                    frontier->complete = false;
            } else {
                detail = assessEdgeDetailed(
                    edgeRequest, state, workspace, context,
                    mixedPresent ? &mixedPresent->metrics.assessmentBindingWork : nullptr);
                syncBindingWorkCompleteness();
            }
            const auto assessed = detail.recovery;
            if (mixedPresent) {
                mixedEdgeDetails.push_back(std::move(detail));
                ++out.schedulerMetrics.routeAssessments;
                if (!countMixedStatus(mixedPresent->metrics.edgeAssessmentEventStatusCounts,
                                      assessed.status, mixedPresent) ||
                    !addMixedGeometryCounters(mixedPresent->metrics.assessmentCounters,
                                              assessed.counters, mixedPresent))
                    throw std::logic_error("mixed edge status overflow");
            } else {
                ++out.schedulerMetrics.batchPreflightAssessments;
            }
            // The batch route is an optimization, not an allocation-failure
            // retry boundary.  A detailed assessment has already translated
            // bad_alloc into its structured outcome; falling through to the
            // generic scheduler would silently consume that fault and could
            // publish a successful mesh.  Re-enter the set-level allocation
            // handler instead so the complete v3 route remains atomic.
            if (assessed.status == TetEdgeRecoveryStatus::Failed &&
                assessed.reason == TetEdgeRecoveryReason::AllocationFailed)
                throw std::bad_alloc{};
            if (assessed.status != TetEdgeRecoveryStatus::AlreadyPresent) {
                allAlreadyPresent = false;
                if (!mixedPresent)
                    break;
            }
            presentEdgeRequests.push_back(edgeRequest);
            presentEdgeOutcomes.push_back(assessed);
        }
        if (mixedFrontierCapacityReached) {
            out.status = TetFacetSetRecoveryStatus::Failed;
            out.reason = TetFacetSetRecoveryReason::FrontierSinkCapacityReached;
            itemizeUnprocessed();
            out.reportComplete = true;
            rollbackUncommittedMixedSchedules();
            return out;
        }
        if (!mixedPresent && allAlreadyPresent) {
            for (const auto *facet : facets) {
                TetFacetIndexedRequest facetRequest;
                if (!mapFacetRequest(*facet, localByGlobal, state, workspace, temporaryRegion,
                                     facetRequest)) {
                    allAlreadyPresent = false;
                    break;
                }
                auto detail = assessFacet(facetRequest);
                const auto assessed = detail.recovery;
                ++out.schedulerMetrics.batchPreflightAssessments;
                if (assessed.status == TetFacetRecoveryStatus::Failed &&
                    assessed.reason == TetFacetRecoveryReason::AllocationFailed)
                    throw std::bad_alloc{};
                if (assessed.status != TetFacetRecoveryStatus::AlreadyPresent) {
                    allAlreadyPresent = false;
                    if (!mixedPresent)
                        break;
                }
                presentFacetRequests.push_back(facetRequest);
                presentFacetOutcomes.push_back(assessed);
            }
        }
        if (mixedPresent) {
            presentFacetRequests.resize(facets.size());
            presentFacetOutcomes.resize(facets.size());
            for (std::size_t i = 0; i < facets.size(); ++i) {
                auto &requestRow = presentFacetRequests[i];
                if (!mapFacetRequest(*facets[i], localByGlobal, state, workspace, temporaryRegion,
                                     requestRow))
                    continue;
                TetFacetRecoveryOutcome candidate;
                if (!exactPresentFacetCandidate(requestRow, workspace, context, candidate))
                    continue;
                auto detail = assessFacet(requestRow);
                const auto assessed = detail.recovery;
                mixedFacetDetails[i] = std::move(detail);
                mixedFacetPreclassified[i] = true;
                presentFacetOutcomes[i] = assessed;
                ++out.schedulerMetrics.routeAssessments;
                if (!countMixedStatus(mixedPresent->metrics.facetAssessmentEventStatusCounts,
                                      assessed.status, mixedPresent) ||
                    !addMixedGeometryCounters(mixedPresent->metrics.assessmentCounters,
                                              assessed.counters, mixedPresent))
                    throw std::logic_error("mixed facet status overflow");
                if (assessed.status == TetFacetRecoveryStatus::Failed &&
                    assessed.reason == TetFacetRecoveryReason::AllocationFailed)
                    throw std::bad_alloc{};
            }
        }
        if (allAlreadyPresent && !mixedPresent) {
            TetFacetSetState scratch = state;
            for (std::size_t i = 0; i < edges.size(); ++i) {
                const auto &requestRow = presentEdgeRequests[i];
                const auto &outcomeRow = presentEdgeOutcomes[i];
                scratch.edgeBindings[requestRow.key] = {requestRow.key, requestRow.classification,
                                                        true, outcomeRow.edge};
                TetFacetSetRecoveryRow reportRow;
                reportRow.kind = TetFacetSetWorkKind::Edge;
                reportRow.edgeKey = requestRow.key;
                reportRow.edgeOutcome = outcomeRow;
                out.processed.push_back(reportRow);
            }
            for (std::size_t i = 0; i < facets.size(); ++i) {
                const auto &requestRow = presentFacetRequests[i];
                const auto &outcomeRow = presentFacetOutcomes[i];
                scratch.bindings[requestRow.key] = {
                    requestRow.key, canonicalHandles(workspace, requestRow.vertices),
                    requestRow.classification, outcomeRow.face};
                TetFacetSetRecoveryRow reportRow;
                reportRow.kind = TetFacetSetWorkKind::Facet;
                reportRow.facetKey = requestRow.key;
                reportRow.facetOutcome = outcomeRow;
                out.processed.push_back(reportRow);
            }
            scratch.workspaceStamp = topo::WorkspaceStamp::of(workspace);
            state.swap(scratch);
            out.schedulerMetrics.alreadyPresentBatchEntities = edges.size() + facets.size();
            out.status = TetFacetSetRecoveryStatus::Complete;
            out.reason = TetFacetSetRecoveryReason::None;
            out.reportComplete = true;
            return out;
        }

        // The protected-edge bootstrap uses the same pure-assessment rule as
        // facets. Missing edges whose exact cavities overlap become typed
        // dependencies rather than order-dependent protected conflicts.
        std::vector<ScheduledEdgeAssessment> edgeAssessments(edges.size());
        std::vector<bool> edgePending(edges.size(), true);
        auto assessEdge = [&](std::size_t i) {
            auto &assessment = edgeAssessments[i];
            assessment = {};
            if (!mapEdgeRequest(*edges[i], localByGlobal, state, workspace, temporaryRegion,
                                assessment.request)) {
                assessment.detail.recovery = edgeFailure(TetEdgeRecoveryStatus::Rejected,
                                                         TetEdgeRecoveryReason::InvalidIdentity,
                                                         TetEdgeRecoveryPhase::Preflight);
            } else {
                if (frontier) {
                    assessment.frontierSupportOffset = frontier->supportCellsUsed;
                    std::size_t used = 0;
                    bool complete = true;
                    assessment.detail.recovery = assessEdgeFrontier(
                        assessment.request, state, workspace, context,
                        frontier->supportCells.subspan(frontier->supportCellsUsed), used, complete,
                        mixedPresent ? &mixedPresent->metrics.assessmentBindingWork : nullptr);
                    syncBindingWorkCompleteness();
                    if (!complete) {
                        frontier->complete = false;
                        mixedFrontierCapacityReached = mixedPresent != nullptr;
                    }
                    if (used > frontier->supportCells.size() - frontier->supportCellsUsed) {
                        frontier->complete = false;
                        mixedFrontierCapacityReached = mixedPresent != nullptr;
                        used = 0;
                    }
                    assessment.frontierSupportCount = used;
                    frontier->supportCellsUsed += used;
                } else {
                    assessment.detail = assessEdgeDetailed(
                        assessment.request, state, workspace, context,
                        mixedPresent ? &mixedPresent->metrics.assessmentBindingWork : nullptr);
                    syncBindingWorkCompleteness();
                }
            }
            if (frontier && !addAssessmentCounters(frontier->assessmentCounters,
                                                   assessment.detail.recovery.counters))
                frontier->complete = false;
            if (mixedPresent &&
                !addMixedGeometryCounters(mixedPresent->metrics.assessmentCounters,
                                          assessment.detail.recovery.counters, mixedPresent))
                throw std::logic_error("mixed edge counter overflow");
            assessment.current = true;
            ++out.schedulerMetrics.routeAssessments;
        };
        const auto edgeSupport = [&](const ScheduledEdgeAssessment &assessment) {
            return frontier
                       ? std::span<const topo::CellHandle>(frontier->supportCells.subspan(
                             assessment.frontierSupportOffset, assessment.frontierSupportCount))
                       : std::span<const topo::CellHandle>(assessment.detail.retiredCells);
        };
        if (mixedPresent) {
            for (std::size_t i = 0; i < edges.size(); ++i) {
                auto &assessment = edgeAssessments[i];
                assessment.request = presentEdgeRequests[i];
                assessment.detail = std::move(mixedEdgeDetails[i]);
                assessment.current = true;
                if (frontier) {
                    assessment.frontierSupportOffset = mixedEdgeSupportOffsets[i];
                    assessment.frontierSupportCount = mixedEdgeSupportCounts[i];
                }
            }
        } else {
            for (std::size_t i = 0; i < edges.size(); ++i)
                assessEdge(i);
        }

        using EdgeReverseRow = std::pair<topo::EdgeHandle, VolumeArtifactEdgeKey>;
        using FaceReverseRow = std::pair<topo::FaceHandle, VolumeFacetKey>;
        std::vector<EdgeReverseRow> mixedEdgeReverse;
        std::vector<FaceReverseRow> mixedFaceReverse;
        std::vector<EdgeReverseRow> mixedEdgeBatch;
        std::vector<FaceReverseRow> mixedFaceBatch;
        std::vector<TetFacetSetRecoveryRow> mixedReportBatch;
        std::vector<std::size_t> mixedBatchOrdinals;
        if (mixedPresent) {
            const auto maximumRows = (std::max)(edges.size(), facets.size());
            std::size_t edgeReverseRows = 0, faceReverseRows = 0;
            if (!addSize(state.edgeBindings.size(), edges.size(), edgeReverseRows) ||
                !addSize(state.bindings.size(), facets.size(), faceReverseRows))
                throw std::logic_error("mixed reverse capacity");
            mixedEdgeReverse.reserve(edgeReverseRows);
            mixedFaceReverse.reserve(faceReverseRows);
            mixedEdgeBatch.reserve(edges.size());
            mixedFaceBatch.reserve(facets.size());
            mixedReportBatch.reserve(maximumRows);
            mixedBatchOrdinals.reserve(maximumRows);
            auto &metrics = mixedPresent->metrics;
            const auto reserveRequest = [&](std::size_t count, std::size_t width) {
                if (count == 0)
                    return true;
                return addMixedMetric(metrics.modeledBatchPreparationAllocationRequests, 1,
                                      mixedPresent) &&
                       addMixedBytes(metrics.modeledBatchPreparationLogicalBytes, count, width,
                                     mixedPresent);
            };
            if (!reserveRequest(edgeReverseRows, sizeof(EdgeReverseRow)) ||
                !reserveRequest(faceReverseRows, sizeof(FaceReverseRow)) ||
                !reserveRequest(edges.size(), sizeof(EdgeReverseRow)) ||
                !reserveRequest(facets.size(), sizeof(FaceReverseRow)) ||
                !reserveRequest(maximumRows, sizeof(TetFacetSetRecoveryRow)) ||
                !reserveRequest(maximumRows, sizeof(std::size_t)))
                throw std::logic_error("mixed preparation accounting");
            for (const auto &[key, binding] : state.edgeBindings)
                mixedEdgeReverse.push_back({binding.edge, key});
            for (const auto &[key, binding] : state.bindings)
                mixedFaceReverse.push_back({binding.recoveredFace, key});
            std::sort(mixedEdgeReverse.begin(), mixedEdgeReverse.end());
            std::sort(mixedFaceReverse.begin(), mixedFaceReverse.end());
            if (std::adjacent_find(mixedEdgeReverse.begin(), mixedEdgeReverse.end(),
                                   [](const auto &a, const auto &b) {
                                       return a.first == b.first && !(a.second == b.second);
                                   }) != mixedEdgeReverse.end() ||
                std::adjacent_find(mixedFaceReverse.begin(), mixedFaceReverse.end(),
                                   [](const auto &a, const auto &b) {
                                       return a.first == b.first && !(a.second == b.second);
                                   }) != mixedFaceReverse.end())
                throw std::logic_error("mixed reverse binding alias");
            std::size_t initialAuditRows = 0;
            std::uint64_t initialAuditCounter = 0;
            if (!addSize(state.edgeBindings.size(), state.bindings.size(), initialAuditRows) ||
                !sizeAsCounter(initialAuditRows, initialAuditCounter) ||
                !addMixedMetric(metrics.schedulerBindingAuditRows, initialAuditCounter,
                                mixedPresent))
                throw std::logic_error("mixed audit accounting");
        }

        const auto recordStateCopy = [&](std::size_t addedRows, std::size_t addedWidth) {
            if (!mixedPresent)
                return true;
            auto &metrics = mixedPresent->metrics;
            std::size_t copiedBindings = 0, copiedSchedules = 0, copied = 0;
            std::uint64_t copiedCounter = 0;
            if (!addSize(state.edgeBindings.size(), state.bindings.size(), copiedBindings) ||
                !addSize(state.canonicalEdgeSchedule.size(), state.canonicalSchedule.size(),
                         copiedSchedules) ||
                !addSize(copiedBindings, copiedSchedules, copied) ||
                !sizeAsCounter(copied, copiedCounter)) {
                metrics.complete = false;
                return false;
            }
            if (!addMixedMetric(metrics.stateCopies, 1, mixedPresent) ||
                !addMixedMetric(metrics.copiedRows, copiedCounter, mixedPresent) ||
                !addMixedBytes(metrics.copiedLogicalBytes, state.edgeBindings.size(),
                               sizeof(decltype(state.edgeBindings)::value_type), mixedPresent) ||
                !addMixedBytes(metrics.copiedLogicalBytes, state.bindings.size(),
                               sizeof(decltype(state.bindings)::value_type), mixedPresent) ||
                !addMixedBytes(metrics.copiedLogicalBytes, state.canonicalEdgeSchedule.size(),
                               sizeof(VolumeArtifactEdgeKey), mixedPresent) ||
                !addMixedBytes(metrics.copiedLogicalBytes, state.canonicalSchedule.size(),
                               sizeof(VolumeFacetKey), mixedPresent) ||
                !addMixedBytes(metrics.modeledBatchPreparationLogicalBytes,
                               state.edgeBindings.size(),
                               sizeof(decltype(state.edgeBindings)::value_type), mixedPresent) ||
                !addMixedBytes(metrics.modeledBatchPreparationLogicalBytes, state.bindings.size(),
                               sizeof(decltype(state.bindings)::value_type), mixedPresent) ||
                !addMixedBytes(metrics.modeledBatchPreparationLogicalBytes,
                               state.canonicalEdgeSchedule.size(), sizeof(VolumeArtifactEdgeKey),
                               mixedPresent) ||
                !addMixedBytes(metrics.modeledBatchPreparationLogicalBytes,
                               state.canonicalSchedule.size(), sizeof(VolumeFacetKey),
                               mixedPresent) ||
                !addMixedBytes(metrics.modeledBatchPreparationLogicalBytes, addedRows, addedWidth,
                               mixedPresent))
                return false;
            std::size_t allocationRows = 0;
            std::uint64_t allocations = 0;
            if (!addSize(copiedBindings, addedRows, allocationRows) ||
                !sizeAsCounter(allocationRows, allocations)) {
                metrics.complete = false;
                return false;
            }
            if (!state.canonicalEdgeSchedule.empty() && !addCounter(allocations, 1)) {
                metrics.complete = false;
                return false;
            }
            if (!state.canonicalSchedule.empty() && !addCounter(allocations, 1)) {
                metrics.complete = false;
                return false;
            }
            return addMixedMetric(metrics.modeledBatchPreparationAllocationRequests, allocations,
                                  mixedPresent);
        };

        const auto reverseEdgeAlias = [&](const EdgeReverseRow &proposed) {
            const auto range =
                std::equal_range(mixedEdgeReverse.begin(), mixedEdgeReverse.end(), proposed,
                                 [](const auto &a, const auto &b) { return a.first < b.first; });
            for (auto it = range.first; it != range.second; ++it)
                if (!(it->second == proposed.second))
                    return true;
            return false;
        };
        const auto reverseFaceAlias = [&](const FaceReverseRow &proposed) {
            const auto range =
                std::equal_range(mixedFaceReverse.begin(), mixedFaceReverse.end(), proposed,
                                 [](const auto &a, const auto &b) { return a.first < b.first; });
            for (auto it = range.first; it != range.second; ++it)
                if (!(it->second == proposed.second))
                    return true;
            return false;
        };

        std::vector<bool> edgePresent;
        std::vector<bool> edgeDone;
        if (mixedPresent) {
            edgePresent.resize(edges.size(), false);
            edgeDone.resize(edges.size(), false);
            for (std::size_t i = 0; i < edges.size(); ++i) {
                edgePresent[i] = edgeAssessments[i].detail.recovery.status ==
                                 TetEdgeRecoveryStatus::AlreadyPresent;
                edgePending[i] = !edgePresent[i];
            }
        }

        const auto bindPresentEdgePrefix = [&](std::size_t limit) {
            if (!mixedPresent)
                return std::size_t{0};
            mixedEdgeBatch.clear();
            mixedReportBatch.clear();
            mixedBatchOrdinals.clear();
            for (std::size_t i = 0; i < limit && i < edges.size(); ++i) {
                if (!edgePresent[i] || edgeDone[i])
                    continue;
#if defined(PAACS_TOPO_TESTING)
                if (mixedEdgeBatch.size() >= mixedPresent->maximumBatchEntitiesForTesting)
                    break;
#endif
                const auto &requestRow = edgeAssessments[i].request;
                const auto &outcomeRow = edgeAssessments[i].detail.recovery;
                if (!exactPresentEdgeStillProtected(requestRow, outcomeRow, workspace, context))
                    throw std::logic_error("mixed exact-present edge lost protection");
                const EdgeReverseRow reverseRow{outcomeRow.edge, requestRow.key};
                if (reverseEdgeAlias(reverseRow))
                    throw std::logic_error("mixed edge binding alias");
                mixedEdgeBatch.push_back(reverseRow);
                mixedBatchOrdinals.push_back(i);
                TetFacetSetRecoveryRow reportRow;
                reportRow.kind = TetFacetSetWorkKind::Edge;
                reportRow.edgeKey = requestRow.key;
                reportRow.edgeOutcome = outcomeRow;
                mixedReportBatch.push_back(reportRow);
            }
            if (mixedEdgeBatch.empty())
                return std::size_t{0};
            std::sort(mixedEdgeBatch.begin(), mixedEdgeBatch.end());
            if (std::adjacent_find(mixedEdgeBatch.begin(), mixedEdgeBatch.end(),
                                   [](const auto &a, const auto &b) {
                                       return a.first == b.first && !(a.second == b.second);
                                   }) != mixedEdgeBatch.end())
                throw std::logic_error("mixed batch edge alias");
            std::size_t inserted = 0;
            for (const auto &row : mixedEdgeBatch)
                if (state.edgeBindings.find(row.second) == state.edgeBindings.end())
                    ++inserted;
            if (!recordStateCopy(inserted, sizeof(decltype(state.edgeBindings)::value_type)) ||
                !addMixedMetric(mixedPresent->metrics.schedulerBindingAuditRows,
                                static_cast<std::uint64_t>(mixedEdgeBatch.size()), mixedPresent))
                throw std::logic_error("mixed edge batch accounting");
            TetFacetSetState scratch = state;
            for (const auto i : mixedBatchOrdinals) {
                const auto &requestRow = edgeAssessments[i].request;
                const auto &outcomeRow = edgeAssessments[i].detail.recovery;
                scratch.edgeBindings[requestRow.key] = {requestRow.key, requestRow.classification,
                                                        true, outcomeRow.edge};
            }
            scratch.workspaceStamp = topo::WorkspaceStamp::of(workspace);
            if (out.processed.capacity() - out.processed.size() < mixedReportBatch.size())
                throw std::logic_error("mixed report capacity");
            const auto committed = mixedEdgeBatch.size();
            auto &metrics = mixedPresent->metrics;
            if (metrics.presentBindingBatches == (std::numeric_limits<std::uint64_t>::max)() ||
                committed >
                    (std::numeric_limits<std::uint64_t>::max)() - metrics.presentBindingEntities)
                throw std::logic_error("mixed edge batch accounting");
            const auto nextBatches = metrics.presentBindingBatches + 1;
            const auto nextEntities =
                metrics.presentBindingEntities + static_cast<std::uint64_t>(committed);
#if defined(PAACS_TOPO_TESTING)
            const MixedPresentCommitZoneScope commitZone;
#endif
            out.processed.insert(out.processed.end(), mixedReportBatch.begin(),
                                 mixedReportBatch.end());
            mixedStateAdvanced = true;
            state.swap(scratch);
            for (const auto &row : mixedEdgeBatch) {
                const auto existing =
                    std::find(mixedEdgeReverse.begin(), mixedEdgeReverse.end(), row);
                if (existing == mixedEdgeReverse.end())
                    mixedEdgeReverse.push_back(row);
            }
            std::sort(mixedEdgeReverse.begin(), mixedEdgeReverse.end());
            for (const auto i : mixedBatchOrdinals)
                edgeDone[i] = true;
            metrics.presentBindingBatches = nextBatches;
            metrics.presentBindingEntities = nextEntities;
            return committed;
        };

        auto rebuildEdgeBlockers = [&](std::size_t i) {
            auto &assessment = edgeAssessments[i];
            assessment.blockers.clear();
            if (!assessment.detail.recovery.succeeded() ||
                assessment.detail.recovery.status != TetEdgeRecoveryStatus::Recovered)
                return;
            for (std::size_t j = 0; j < edges.size(); ++j) {
                if (i == j)
                    continue;
                if (mixedPresent && edgePresent[j]) {
                    if (edgePending[j] &&
                        !addMixedMetric(mixedPresent->metrics.exactPresentDependencyTests, 1,
                                        mixedPresent))
                        throw std::logic_error("mixed present dependency accounting");
                    continue;
                }
                if (!edgePending[j])
                    continue;
                if (mixedPresent &&
                    !addMixedMetric(mixedPresent->metrics.dependencyPairsTested, 1, mixedPresent))
                    throw std::logic_error("mixed dependency accounting");
                if (selectedCavityTouchesEdge(edgeSupport(assessment), edgeAssessments[j].request,
                                              context, workspace)) {
                    assessment.blockers.push_back(j);
                    ++out.schedulerMetrics.dependencyEdges;
                }
            }
        };

        auto appendEdgeStop = [&](std::size_t stopped, TetEdgeRecoveryOutcome outcome,
                                  TetFacetSetRecoveryStatus status,
                                  TetFacetSetRecoveryReason reason) {
            TetFacetSetRecoveryRow row;
            row.kind = TetFacetSetWorkKind::Edge;
            row.edgeKey = edges[stopped]->key;
            row.edgeOutcome = outcome;
            out.processed.push_back(row);
            for (std::size_t j = 0; j < edges.size(); ++j)
                if ((mixedPresent ? !edgeDone[j] : edgePending[j]) && j != stopped)
                    out.unattemptedEdges.push_back(edges[j]->key);
            for (const auto *facet : facets)
                out.unattemptedFacets.push_back(facet->key);
            out.status = status;
            out.reason = reason;
            out.reportComplete = true;
        };

        std::size_t remainingEdges =
            mixedPresent
                ? static_cast<std::size_t>(std::count(edgePending.begin(), edgePending.end(), true))
                : edges.size();
        std::size_t remainingPresentEdges =
            mixedPresent
                ? static_cast<std::size_t>(std::count(edgePresent.begin(), edgePresent.end(), true))
                : 0U;
        while (remainingEdges != 0 || remainingPresentEdges != 0) {
            for (std::size_t i = 0; i < edges.size(); ++i) {
                if (!edgePending[i])
                    continue;
                if (mixedPresent &&
                    !addMixedMetric(mixedPresent->metrics.pendingRowsVisited, 1, mixedPresent))
                    throw std::logic_error("mixed pending accounting");
                if (!edgeAssessments[i].current) {
                    assessEdge(i);
                    ++out.schedulerMetrics.dependencyInvalidations;
                    if (mixedFrontierCapacityReached) {
                        out.status = TetFacetSetRecoveryStatus::Failed;
                        out.reason = TetFacetSetRecoveryReason::FrontierSinkCapacityReached;
                        itemizeUnprocessed();
                        out.reportComplete = true;
                        return out;
                    }
                }
                rebuildEdgeBlockers(i);
            }
            std::size_t permanent = edges.size();
            for (std::size_t i = 0; i < edges.size(); ++i)
                if (edgePending[i] && !edgeAssessments[i].detail.recovery.succeeded()) {
                    permanent = i;
                    break;
                }
            if (permanent != edges.size()) {
                if (frontier) {
                    for (std::size_t i = 0; i < edges.size(); ++i) {
                        if (!edgePending[i] || edgeAssessments[i].detail.recovery.succeeded())
                            continue;
                        if (frontier->rowsUsed >= frontier->rows.size()) {
                            frontier->complete = false;
                            break;
                        }
                        auto &frontierRow = frontier->rows[frontier->rowsUsed++];
                        frontierRow = {};
                        frontierRow.kind = TetFacetSetWorkKind::Edge;
                        frontierRow.edgeKey = edges[i]->key;
                        frontierRow.workspaceStamp = topo::WorkspaceStamp::of(workspace);
                        frontierRow.artifactStamp = artifact.surfaceSetStamp;
                        frontierRow.edgeOutcome = edgeAssessments[i].detail.recovery;
                        frontierRow.supportOffset = edgeAssessments[i].frontierSupportOffset;
                        frontierRow.supportCount = edgeAssessments[i].frontierSupportCount;
                        frontierRow.dependencyOffset = frontier->dependenciesUsed;
                        for (const auto dependency : edgeAssessments[i].blockers) {
                            if (frontier->dependenciesUsed >= frontier->dependencies.size()) {
                                frontier->complete = false;
                                break;
                            }
                            frontier->dependencies[frontier->dependenciesUsed++] = dependency;
                            ++frontierRow.dependencyCount;
                        }
                    }
                }
                const auto nested = edgeAssessments[permanent].detail.recovery;
                if (nested.reason == TetEdgeRecoveryReason::ProtectedConstraintConflict)
                    ++out.schedulerMetrics.protectedConflicts;
                appendEdgeStop(permanent, nested,
                               nested.status == TetEdgeRecoveryStatus::Rejected
                                   ? TetFacetSetRecoveryStatus::Rejected
                               : nested.status == TetEdgeRecoveryStatus::Refused
                                   ? TetFacetSetRecoveryStatus::Refused
                                   : TetFacetSetRecoveryStatus::Failed,
                               TetFacetSetRecoveryReason::None);
                return out;
            }
            std::size_t selected = edges.size();
            for (std::size_t i = 0; i < edges.size(); ++i)
                if (edgePending[i] && edgeAssessments[i].blockers.empty()) {
                    selected = i;
                    break;
                }
            if (mixedPresent) {
                const auto bound =
                    bindPresentEdgePrefix(selected == edges.size() ? edges.size() : selected);
                if (bound > remainingPresentEdges)
                    throw std::logic_error("mixed edge prefix count");
                remainingPresentEdges -= bound;
                if (bound != 0)
                    continue;
            }
            if (selected == edges.size()) {
                std::vector<std::vector<std::size_t>> graph(edges.size());
                for (std::size_t i = 0; i < edges.size(); ++i)
                    graph[i] = edgeAssessments[i].blockers;
                if (!blockerGraphValid(graph, edgePending)) {
                    std::size_t first = 0;
                    while (first < edges.size() && !edgePending[first])
                        ++first;
                    appendEdgeStop(first,
                                   edgeFailure(TetEdgeRecoveryStatus::Failed,
                                               TetEdgeRecoveryReason::PreparationConflict,
                                               TetEdgeRecoveryPhase::Routing),
                                   TetFacetSetRecoveryStatus::Failed,
                                   TetFacetSetRecoveryReason::InternalError);
                    return out;
                }
                const auto cycle = canonicalConflictCycle(graph, edgePending,
                                                          out.schedulerMetrics.conflictSccVisits);
                if (cycle.empty()) {
                    std::size_t first = 0;
                    while (first < edges.size() && !edgePending[first])
                        ++first;
                    appendEdgeStop(first,
                                   edgeFailure(TetEdgeRecoveryStatus::Failed,
                                               TetEdgeRecoveryReason::PreparationConflict,
                                               TetEdgeRecoveryPhase::Routing),
                                   TetFacetSetRecoveryStatus::Failed,
                                   TetFacetSetRecoveryReason::InternalError);
                    return out;
                }
                for (const auto member : cycle)
                    out.conflictEdgeVertices.push_back(edges[member]->key);
                for (const auto member : cycle)
                    for (const auto blocker : edgeAssessments[member].blockers)
                        if (std::binary_search(cycle.begin(), cycle.end(), blocker))
                            out.conflictEdgeDependencies.push_back(
                                {edges[member]->key, edges[blocker]->key});
                std::sort(out.conflictEdgeDependencies.begin(), out.conflictEdgeDependencies.end());
                out.conflictEdgeDependencies.erase(std::unique(out.conflictEdgeDependencies.begin(),
                                                               out.conflictEdgeDependencies.end()),
                                                   out.conflictEdgeDependencies.end());
                const auto cycleFailure = edgeFailure(TetEdgeRecoveryStatus::Refused,
                                                      TetEdgeRecoveryReason::RecoveryConflictCycle,
                                                      TetEdgeRecoveryPhase::Routing);
                if (frontier) {
                    if (frontier->rowsUsed >= frontier->rows.size())
                        frontier->complete = false;
                    else {
                        auto &row = frontier->rows[frontier->rowsUsed++];
                        row = {};
                        row.kind = TetFacetSetWorkKind::Edge;
                        row.edgeKey = edges[cycle.front()]->key;
                        row.workspaceStamp = topo::WorkspaceStamp::of(workspace);
                        row.artifactStamp = artifact.surfaceSetStamp;
                        row.edgeOutcome = cycleFailure;
                        row.supportOffset = edgeAssessments[cycle.front()].frontierSupportOffset;
                        row.supportCount = edgeAssessments[cycle.front()].frontierSupportCount;
                        row.dependencyOffset = frontier->dependenciesUsed;
                        for (const auto blocker : edgeAssessments[cycle.front()].blockers) {
                            if (frontier->dependenciesUsed >= frontier->dependencies.size()) {
                                frontier->complete = false;
                                break;
                            }
                            frontier->dependencies[frontier->dependenciesUsed++] = blocker;
                            ++row.dependencyCount;
                        }
                    }
                }
                appendEdgeStop(cycle.front(), cycleFailure, TetFacetSetRecoveryStatus::Refused,
                               TetFacetSetRecoveryReason::RecoveryConflictCycle);
                return out;
            }

            if (mixedPresent && edgePresent[selected]) {
                if (!addMixedMetric(mixedPresent->metrics.exactPresentSelectedRebuilds, 1,
                                    mixedPresent))
                    throw std::logic_error("mixed present rebuild accounting");
                throw std::logic_error("mixed present edge selected");
            }
            ++out.schedulerMetrics.selectedRebuilds;
            assessEdge(selected);
            if (mixedFrontierCapacityReached) {
                out.status = TetFacetSetRecoveryStatus::Failed;
                out.reason = TetFacetSetRecoveryReason::FrontierSinkCapacityReached;
                itemizeUnprocessed();
                out.reportComplete = true;
                return out;
            }
            rebuildEdgeBlockers(selected);
            if (!edgeAssessments[selected].detail.recovery.succeeded() ||
                !edgeAssessments[selected].blockers.empty())
                continue;
            std::vector<bool> invalidated(edges.size(), false);
            for (std::size_t i = 0; i < edges.size(); ++i)
                if (edgePending[i] && i != selected)
                    invalidated[i] =
                        selectedCavityTouchesEdge(edgeSupport(edgeAssessments[selected]),
                                                  edgeAssessments[i].request, context, workspace);
            std::uint64_t preparedNontrivialRecoveries = 0;
            if (mixedPresent) {
                preparedNontrivialRecoveries = mixedPresent->metrics.nontrivialRecoveries;
                if (!addCounter(preparedNontrivialRecoveries, 1)) {
                    mixedPresent->metrics.complete = false;
                    throw std::logic_error("mixed recovery accounting");
                }
            }
            if (mixedPresent)
                mixedStateAdvanced = true;
            auto recovered = recoverEdgeDetailed(
                edgeAssessments[selected].request, state, workspace, context,
                mixedPresent ? &mixedPresent->metrics.recoveryBindingWork : nullptr);
            syncBindingWorkCompleteness();
            if (mixedPresent && recovered.recovery.status == TetEdgeRecoveryStatus::Recovered)
                mixedPresent->metrics.nontrivialRecoveries = preparedNontrivialRecoveries;
            if (mixedPresent &&
                !addMixedGeometryCounters(mixedPresent->metrics.recoveryCounters,
                                          recovered.recovery.counters, mixedPresent))
                mixedPresent->metrics.complete = false;
            TetFacetSetRecoveryRow row;
            row.kind = TetFacetSetWorkKind::Edge;
            row.edgeKey = edges[selected]->key;
            row.edgeOutcome = recovered.recovery;
            out.processed.push_back(row);
            if (!row.edgeOutcome.succeeded()) {
                edgePending[selected] = false;
                if (mixedPresent)
                    edgeDone[selected] = true;
                --remainingEdges;
                for (std::size_t i = 0; i < edges.size(); ++i)
                    if (mixedPresent ? !edgeDone[i] : edgePending[i])
                        out.unattemptedEdges.push_back(edges[i]->key);
                for (const auto *facet : facets)
                    out.unattemptedFacets.push_back(facet->key);
                out.status = row.edgeOutcome.status == TetEdgeRecoveryStatus::Rejected
                                 ? TetFacetSetRecoveryStatus::Rejected
                             : row.edgeOutcome.status == TetEdgeRecoveryStatus::Refused
                                 ? TetFacetSetRecoveryStatus::Refused
                                 : TetFacetSetRecoveryStatus::Failed;
                out.reportComplete = true;
                return out;
            }
            edgePending[selected] = false;
            if (mixedPresent) {
                edgeDone[selected] = true;
                const auto binding = state.edgeBindings.find(edgeAssessments[selected].request.key);
                if (binding == state.edgeBindings.end())
                    throw std::logic_error("mixed recovered edge binding");
                mixedEdgeReverse.push_back({binding->second.edge, binding->first});
                std::sort(mixedEdgeReverse.begin(), mixedEdgeReverse.end());
            }
            --remainingEdges;
            for (std::size_t i = 0; i < edges.size(); ++i) {
                if (!edgePending[i])
                    continue;
                if (!invalidated[i])
                    invalidated[i] = selectedCavityTouchesEdge(
                        recovered.createdCells, edgeAssessments[i].request, context, workspace);
                if (!invalidated[i] && std::find(edgeAssessments[i].blockers.begin(),
                                                 edgeAssessments[i].blockers.end(),
                                                 selected) != edgeAssessments[i].blockers.end())
                    invalidated[i] = true;
                if (invalidated[i])
                    edgeAssessments[i].current = false;
            }
        }

        // Facet scheduling retains only plain dependency witnesses. Candidate
        // topology is rebuilt for the selected entity at the exact current
        // stamp, while disjoint cached assessments survive unrelated deltas.
        std::vector<ScheduledFacetAssessment> assessments(facets.size());
        std::vector<bool> pending(facets.size(), true);
        auto assessOne = [&](std::size_t i) {
            auto &assessment = assessments[i];
            assessment = {};
            if (!mapFacetRequest(*facets[i], localByGlobal, state, workspace, temporaryRegion,
                                 assessment.request)) {
                assessment.detail.recovery = facetFailure(TetFacetRecoveryStatus::Rejected,
                                                          TetFacetRecoveryReason::InvalidIdentity,
                                                          TetFacetRecoveryPhase::Preflight);
            } else {
                assessment.targetBounds = facetBounds(assessment.request, workspace);
                assessment.detail = assessFacet(assessment.request);
            }
            if (frontier && !addAssessmentCounters(frontier->assessmentCounters,
                                                   assessment.detail.recovery.counters))
                frontier->complete = false;
            if (mixedPresent &&
                !addMixedGeometryCounters(mixedPresent->metrics.assessmentCounters,
                                          assessment.detail.recovery.counters, mixedPresent))
                throw std::logic_error("mixed facet counter overflow");
            assessment.current = true;
            ++out.schedulerMetrics.routeAssessments;
        };
        if (mixedPresent) {
            for (std::size_t i = 0; i < facets.size(); ++i) {
                if (mixedFacetPreclassified[i]) {
                    assessments[i].request = presentFacetRequests[i];
                    assessments[i].detail = std::move(mixedFacetDetails[i]);
                    assessments[i].targetBounds = facetBounds(assessments[i].request, workspace);
                    assessments[i].current = true;
                    if (frontier && !addAssessmentCounters(frontier->assessmentCounters,
                                                           assessments[i].detail.recovery.counters))
                        frontier->complete = false;
                } else {
                    assessOne(i);
                    if (!countMixedStatus(mixedPresent->metrics.facetAssessmentEventStatusCounts,
                                          assessments[i].detail.recovery.status, mixedPresent))
                        throw std::logic_error("mixed facet status overflow");
                }
            }
        } else {
            for (std::size_t i = 0; i < facets.size(); ++i) {
                assessOne(i);
            }
        }

        std::vector<bool> facetPresent;
        std::vector<bool> facetDone;
        if (mixedPresent) {
            facetPresent.resize(facets.size(), false);
            facetDone.resize(facets.size(), false);
            for (std::size_t i = 0; i < facets.size(); ++i) {
                facetPresent[i] =
                    assessments[i].detail.recovery.status == TetFacetRecoveryStatus::AlreadyPresent;
                pending[i] = !facetPresent[i];
            }
        }

        const auto bindPresentFacetPrefix = [&](std::size_t limit) {
            if (!mixedPresent)
                return std::size_t{0};
            mixedFaceBatch.clear();
            mixedReportBatch.clear();
            mixedBatchOrdinals.clear();
            for (std::size_t i = 0; i < limit && i < facets.size(); ++i) {
                if (!facetPresent[i] || facetDone[i])
                    continue;
#if defined(PAACS_TOPO_TESTING)
                if (mixedFaceBatch.size() >= mixedPresent->maximumBatchEntitiesForTesting)
                    break;
#endif
                const auto &requestRow = assessments[i].request;
                const auto &outcomeRow = assessments[i].detail.recovery;
                if (!exactPresentFacetStillProtected(requestRow, outcomeRow, workspace, context))
                    throw std::logic_error("mixed exact-present facet lost protection");
                const FaceReverseRow reverseRow{outcomeRow.face, requestRow.key};
                if (reverseFaceAlias(reverseRow))
                    throw std::logic_error("mixed facet binding alias");
                mixedFaceBatch.push_back(reverseRow);
                mixedBatchOrdinals.push_back(i);
                TetFacetSetRecoveryRow reportRow;
                reportRow.kind = TetFacetSetWorkKind::Facet;
                reportRow.facetKey = requestRow.key;
                reportRow.facetOutcome = outcomeRow;
                mixedReportBatch.push_back(reportRow);
            }
            if (mixedFaceBatch.empty())
                return std::size_t{0};
            std::sort(mixedFaceBatch.begin(), mixedFaceBatch.end());
            if (std::adjacent_find(mixedFaceBatch.begin(), mixedFaceBatch.end(),
                                   [](const auto &a, const auto &b) {
                                       return a.first == b.first && !(a.second == b.second);
                                   }) != mixedFaceBatch.end())
                throw std::logic_error("mixed batch facet alias");
            std::size_t inserted = 0;
            for (const auto &row : mixedFaceBatch)
                if (state.bindings.find(row.second) == state.bindings.end())
                    ++inserted;
            if (!recordStateCopy(inserted, sizeof(decltype(state.bindings)::value_type)) ||
                !addMixedMetric(mixedPresent->metrics.schedulerBindingAuditRows,
                                static_cast<std::uint64_t>(mixedFaceBatch.size()), mixedPresent))
                throw std::logic_error("mixed facet batch accounting");
            TetFacetSetState scratch = state;
            for (const auto i : mixedBatchOrdinals) {
                const auto &requestRow = assessments[i].request;
                const auto &outcomeRow = assessments[i].detail.recovery;
                scratch.bindings[requestRow.key] = {
                    requestRow.key, canonicalHandles(workspace, requestRow.vertices),
                    requestRow.classification, outcomeRow.face};
            }
            scratch.workspaceStamp = topo::WorkspaceStamp::of(workspace);
            if (mixedFacetBindingStateRevision == (std::numeric_limits<std::uint64_t>::max)()) {
                mixedPresent->metrics.complete = false;
                throw std::logic_error("mixed facet receipt revision");
            }
            const auto nextBindingStateRevision = mixedFacetBindingStateRevision + 1;
            if (!validBindings(scratch, workspace, &mixedPresent->metrics.assessmentBindingWork,
                               &mixedFacetBindingReceipt, nextBindingStateRevision))
                throw std::logic_error("mixed prepared facet bindings invalid");
            syncBindingWorkCompleteness();
            if (out.processed.capacity() - out.processed.size() < mixedReportBatch.size())
                throw std::logic_error("mixed report capacity");
            const auto committed = mixedFaceBatch.size();
            auto &metrics = mixedPresent->metrics;
            if (metrics.presentBindingBatches == (std::numeric_limits<std::uint64_t>::max)() ||
                committed >
                    (std::numeric_limits<std::uint64_t>::max)() - metrics.presentBindingEntities)
                throw std::logic_error("mixed facet batch accounting");
            const auto nextBatches = metrics.presentBindingBatches + 1;
            const auto nextEntities =
                metrics.presentBindingEntities + static_cast<std::uint64_t>(committed);
#if defined(PAACS_TOPO_TESTING)
            const MixedPresentCommitZoneScope commitZone;
#endif
            out.processed.insert(out.processed.end(), mixedReportBatch.begin(),
                                 mixedReportBatch.end());
            mixedStateAdvanced = true;
            state.swap(scratch);
            mixedFacetBindingReceipt.state = &state;
            mixedFacetBindingStateRevision = nextBindingStateRevision;
            for (const auto &row : mixedFaceBatch) {
                const auto existing =
                    std::find(mixedFaceReverse.begin(), mixedFaceReverse.end(), row);
                if (existing == mixedFaceReverse.end())
                    mixedFaceReverse.push_back(row);
            }
            std::sort(mixedFaceReverse.begin(), mixedFaceReverse.end());
            for (const auto i : mixedBatchOrdinals)
                facetDone[i] = true;
            metrics.presentBindingBatches = nextBatches;
            metrics.presentBindingEntities = nextEntities;
            return committed;
        };

        auto rebuildBlockers = [&](std::size_t i) {
            auto &assessment = assessments[i];
            assessment.blockers.clear();
            if (!facetReady(assessment.detail.recovery) ||
                assessment.detail.recovery.status != TetFacetRecoveryStatus::Recovered)
                return;
            for (std::size_t j = 0; j < facets.size(); ++j) {
                if (i == j)
                    continue;
                if (mixedPresent && facetPresent[j]) {
                    if (pending[j] &&
                        !addMixedMetric(mixedPresent->metrics.exactPresentDependencyTests, 1,
                                        mixedPresent))
                        throw std::logic_error("mixed present dependency accounting");
                    continue;
                }
                if (!pending[j])
                    continue;
                if (mixedPresent &&
                    !addMixedMetric(mixedPresent->metrics.dependencyPairsTested, 1, mixedPresent))
                    throw std::logic_error("mixed dependency accounting");
                if (selectedCavityTouchesTarget(assessment.detail.retiredCells,
                                                assessments[j].request, context)) {
                    assessment.blockers.push_back(j);
                    ++out.schedulerMetrics.dependencyEdges;
                }
            }
        };

        auto appendStop = [&](std::size_t stopped, TetFacetRecoveryOutcome outcome,
                              TetFacetSetRecoveryStatus status, TetFacetSetRecoveryReason reason) {
            TetFacetSetRecoveryRow row;
            row.kind = TetFacetSetWorkKind::Facet;
            row.facetKey = facets[stopped]->key;
            row.facetOutcome = outcome;
            out.processed.push_back(row);
            for (std::size_t j = 0; j < facets.size(); ++j)
                if ((mixedPresent ? !facetDone[j] : pending[j]) && j != stopped)
                    out.unattemptedFacets.push_back(facets[j]->key);
            out.status = status;
            out.reason = reason;
            out.reportComplete = true;
        };

        std::size_t remaining =
            mixedPresent
                ? static_cast<std::size_t>(std::count(pending.begin(), pending.end(), true))
                : facets.size();
        std::size_t remainingPresentFacets =
            mixedPresent ? static_cast<std::size_t>(
                               std::count(facetPresent.begin(), facetPresent.end(), true))
                         : 0U;
        while (remaining != 0 || remainingPresentFacets != 0) {
            for (std::size_t i = 0; i < facets.size(); ++i) {
                if (!pending[i])
                    continue;
                if (mixedPresent &&
                    !addMixedMetric(mixedPresent->metrics.pendingRowsVisited, 1, mixedPresent))
                    throw std::logic_error("mixed pending accounting");
                if (!assessments[i].current) {
                    assessOne(i);
                    ++out.schedulerMetrics.dependencyInvalidations;
                }
                rebuildBlockers(i);
            }

            // A permanent assessment stops the round before any currently
            // ready entity can mutate the private workspace.
            std::size_t permanent = facets.size();
            for (std::size_t i = 0; i < facets.size(); ++i)
                if (pending[i] && !facetReady(assessments[i].detail.recovery)) {
                    permanent = i;
                    break;
                }
            if (permanent != facets.size()) {
                if (frontier) {
                    for (std::size_t i = 0; i < facets.size(); ++i) {
                        if (!pending[i] || facetReady(assessments[i].detail.recovery))
                            continue;
                        if (frontier->rowsUsed >= frontier->rows.size()) {
                            frontier->complete = false;
                            break;
                        }
                        auto &frontierRow = frontier->rows[frontier->rowsUsed++];
                        frontierRow = {};
                        frontierRow.kind = TetFacetSetWorkKind::Facet;
                        frontierRow.facetKey = facets[i]->key;
                        frontierRow.workspaceStamp = topo::WorkspaceStamp::of(workspace);
                        frontierRow.artifactStamp = artifact.surfaceSetStamp;
                        frontierRow.facetOutcome = assessments[i].detail.recovery;
                        frontierRow.supportOffset = frontier->supportCellsUsed;
                        for (const auto cell : assessments[i].detail.retiredCells) {
                            if (frontier->supportCellsUsed >= frontier->supportCells.size()) {
                                frontier->complete = false;
                                break;
                            }
                            frontier->supportCells[frontier->supportCellsUsed++] = cell;
                            ++frontierRow.supportCount;
                        }
                        frontierRow.dependencyOffset = frontier->dependenciesUsed;
                        for (const auto dependency : assessments[i].blockers) {
                            if (frontier->dependenciesUsed >= frontier->dependencies.size()) {
                                frontier->complete = false;
                                break;
                            }
                            frontier->dependencies[frontier->dependenciesUsed++] = dependency;
                            ++frontierRow.dependencyCount;
                        }
                    }
                }
                const auto nested = assessments[permanent].detail.recovery;
                if (nested.reason == TetFacetRecoveryReason::ProtectedFacetConflict ||
                    nested.reason == TetFacetRecoveryReason::CrossedProtectedEdge ||
                    nested.reason == TetFacetRecoveryReason::CrossedProtectedFacet)
                    ++out.schedulerMetrics.protectedConflicts;
                appendStop(permanent, nested,
                           nested.status == TetFacetRecoveryStatus::Rejected
                               ? TetFacetSetRecoveryStatus::Rejected
                           : nested.status == TetFacetRecoveryStatus::Refused
                               ? TetFacetSetRecoveryStatus::Refused
                               : TetFacetSetRecoveryStatus::Failed,
                           TetFacetSetRecoveryReason::None);
                return out;
            }

            std::size_t selected = facets.size();
            for (std::size_t i = 0; i < facets.size(); ++i)
                if (pending[i] && assessments[i].blockers.empty()) {
                    selected = i;
                    break;
                }

            if (mixedPresent) {
                const auto bound =
                    bindPresentFacetPrefix(selected == facets.size() ? facets.size() : selected);
                if (bound > remainingPresentFacets)
                    throw std::logic_error("mixed facet prefix count");
                remainingPresentFacets -= bound;
                if (bound != 0)
                    continue;
            }

            if (selected == facets.size()) {
                // Every pending entity is temporary. Deterministic Tarjan
                // traversal identifies the minimum canonical cycle.
                std::vector<std::vector<std::size_t>> graph(facets.size());
                for (std::size_t i = 0; i < facets.size(); ++i)
                    graph[i] = assessments[i].blockers;
                if (!blockerGraphValid(graph, pending)) {
                    std::size_t first = 0;
                    while (first < facets.size() && !pending[first])
                        ++first;
                    appendStop(first,
                               facetFailure(TetFacetRecoveryStatus::Failed,
                                            TetFacetRecoveryReason::PreparationConflict,
                                            TetFacetRecoveryPhase::Routing),
                               TetFacetSetRecoveryStatus::Failed,
                               TetFacetSetRecoveryReason::InternalError);
                    return out;
                }
                const auto cycle =
                    canonicalConflictCycle(graph, pending, out.schedulerMetrics.conflictSccVisits);
                if (cycle.empty()) {
                    std::size_t first = 0;
                    while (first < facets.size() && !pending[first])
                        ++first;
                    appendStop(first,
                               facetFailure(TetFacetRecoveryStatus::Failed,
                                            TetFacetRecoveryReason::PreparationConflict,
                                            TetFacetRecoveryPhase::Routing),
                               TetFacetSetRecoveryStatus::Failed,
                               TetFacetSetRecoveryReason::InternalError);
                    return out;
                }
#if defined(PAACS_TOPO_TESTING)
                if (gFacetConflictDiagnostics.enabled) {
                    const bool enabled = gFacetConflictDiagnostics.enabled;
                    gFacetConflictDiagnostics = {};
                    auto &diagnostic = gFacetConflictDiagnostics;
                    diagnostic.enabled = enabled;
                    diagnostic.captured = true;
                    diagnostic.complete = true;
                    diagnostic.workspaceStamp = topo::WorkspaceStamp::of(workspace);
                    diagnostic.artifactStamp = artifact.surfaceSetStamp;
                    diagnostic.facetCount = facets.size();
                    diagnostic.pendingCount =
                        static_cast<std::size_t>(std::count(pending.begin(), pending.end(), true));
                    diagnostic.schedulerMetrics = out.schedulerMetrics;
                    for (std::size_t i = 0; i < assessments.size(); ++i) {
                        const auto &assessment = assessments[i];
                        if (!addAssessmentCounters(diagnostic.allAssessmentCounters,
                                                   assessment.detail.recovery.counters))
                            diagnostic.complete = false;
                        if (assessment.detail.retiredCells.size() >
                                (std::numeric_limits<std::uint64_t>::max)() ||
                            !addCounter(
                                diagnostic.allRetiredCellRows,
                                static_cast<std::uint64_t>(assessment.detail.retiredCells.size())))
                            diagnostic.complete = false;
                        diagnostic.maximumRetiredCells =
                            (std::max)(diagnostic.maximumRetiredCells,
                                       assessment.detail.retiredCells.size());
                        if (!pending[i])
                            continue;
                        if (!addAssessmentCounters(diagnostic.pendingAssessmentCounters,
                                                   assessment.detail.recovery.counters))
                            diagnostic.complete = false;
                        if (assessment.detail.retiredCells.size() >
                                (std::numeric_limits<std::uint64_t>::max)() ||
                            !addCounter(
                                diagnostic.pendingRetiredCellRows,
                                static_cast<std::uint64_t>(assessment.detail.retiredCells.size())))
                            diagnostic.complete = false;
                    }
                    for (const auto member : cycle) {
                        if (diagnostic.memberCount >= diagnostic.members.size()) {
                            diagnostic.complete = false;
                            continue;
                        }
                        diagnostic.members[diagnostic.memberCount++] = facets[member]->key;
                    }
                    for (const auto member : cycle)
                        for (const auto blocker : assessments[member].blockers)
                            if (std::binary_search(cycle.begin(), cycle.end(), blocker)) {
                                if (diagnostic.dependencyCount >= diagnostic.dependencies.size()) {
                                    diagnostic.complete = false;
                                    continue;
                                }
                                diagnostic.dependencies[diagnostic.dependencyCount++] = {
                                    facets[member]->key, facets[blocker]->key};
                            }
                }
#endif
                if (conflictCycle)
                    beginFacetConflictCycleCapture(*conflictCycle,
                                                   topo::WorkspaceStamp::of(workspace),
                                                   artifact.surfaceSetStamp, cycle.size());
                for (const auto member : cycle) {
                    out.conflictVertices.push_back(facets[member]->key);
                    if (conflictCycle)
                        appendFacetConflictCycleMember(*conflictCycle, facets[member]->key);
                }
                for (const auto member : cycle)
                    for (const auto blocker : assessments[member].blockers)
                        if (std::binary_search(cycle.begin(), cycle.end(), blocker)) {
                            out.conflictEdges.push_back(
                                {facets[member]->key, facets[blocker]->key});
                            if (conflictCycle)
                                appendFacetConflictCycleDependency(
                                    *conflictCycle, {facets[member]->key, facets[blocker]->key});
                        }
                std::sort(out.conflictEdges.begin(), out.conflictEdges.end());
                out.conflictEdges.erase(
                    std::unique(out.conflictEdges.begin(), out.conflictEdges.end()),
                    out.conflictEdges.end());
                if (conflictCycle)
                    finishFacetConflictCycleCapture(*conflictCycle);
                const auto cycleFailure = facetFailure(
                    TetFacetRecoveryStatus::Refused, TetFacetRecoveryReason::RecoveryConflictCycle,
                    TetFacetRecoveryPhase::Routing);
                if (frontier) {
                    if (frontier->rowsUsed >= frontier->rows.size())
                        frontier->complete = false;
                    else {
                        auto &row = frontier->rows[frontier->rowsUsed++];
                        row = {};
                        row.kind = TetFacetSetWorkKind::Facet;
                        row.facetKey = facets[cycle.front()]->key;
                        row.workspaceStamp = topo::WorkspaceStamp::of(workspace);
                        row.artifactStamp = artifact.surfaceSetStamp;
                        row.facetOutcome = cycleFailure;
                        row.supportOffset = frontier->supportCellsUsed;
                        for (const auto cell : assessments[cycle.front()].detail.retiredCells) {
                            if (frontier->supportCellsUsed >= frontier->supportCells.size()) {
                                frontier->complete = false;
                                break;
                            }
                            frontier->supportCells[frontier->supportCellsUsed++] = cell;
                            ++row.supportCount;
                        }
                        row.dependencyOffset = frontier->dependenciesUsed;
                        for (const auto blocker : assessments[cycle.front()].blockers) {
                            if (frontier->dependenciesUsed >= frontier->dependencies.size()) {
                                frontier->complete = false;
                                break;
                            }
                            frontier->dependencies[frontier->dependenciesUsed++] = blocker;
                            ++row.dependencyCount;
                        }
                    }
                }
                appendStop(cycle.front(), cycleFailure, TetFacetSetRecoveryStatus::Refused,
                           conflictCycle && !conflictCycle->complete
                               ? TetFacetSetRecoveryReason::ConflictCycleSinkCapacityReached
                               : TetFacetSetRecoveryReason::RecoveryConflictCycle);
                return out;
            }

            // Rebuild the selected route at the exact stamp and require the
            // retained ready/temporary discriminator to agree.
            if (mixedPresent && facetPresent[selected]) {
                if (!addMixedMetric(mixedPresent->metrics.exactPresentSelectedRebuilds, 1,
                                    mixedPresent))
                    throw std::logic_error("mixed present rebuild accounting");
                throw std::logic_error("mixed present facet selected");
            }
            ++out.schedulerMetrics.selectedRebuilds;
            assessOne(selected);
            rebuildBlockers(selected);
            if (!facetReady(assessments[selected].detail.recovery) ||
                !assessments[selected].blockers.empty())
                continue;

            // Determine exact dependency invalidation before the commit
            // retires any witness cells.
            std::vector<bool> invalidated(facets.size(), false);
            for (std::size_t i = 0; i < facets.size(); ++i) {
                if (!pending[i] || i == selected)
                    continue;
                invalidated[i] = selectedCavityTouchesTarget(
                    assessments[selected].detail.retiredCells, assessments[i].request, context);
            }

            std::uint64_t preparedNontrivialRecoveries = 0;
            if (mixedPresent) {
                preparedNontrivialRecoveries = mixedPresent->metrics.nontrivialRecoveries;
                if (!addCounter(preparedNontrivialRecoveries, 1)) {
                    mixedPresent->metrics.complete = false;
                    throw std::logic_error("mixed recovery accounting");
                }
            }
            if (mixedPresent)
                mixedStateAdvanced = true;
            std::uint64_t nextBindingStateRevision = mixedFacetBindingStateRevision;
            if (mixedPresent) {
                if (nextBindingStateRevision == (std::numeric_limits<std::uint64_t>::max)()) {
                    mixedPresent->metrics.complete = false;
                    throw std::logic_error("mixed facet receipt revision");
                }
                ++nextBindingStateRevision;
            }
            auto recovered = mixedPresent
                                 ? recoverTetFacetIndexedImpl(
                                       assessments[selected].request, state, workspace, context,
                                       false, &mixedPresent->metrics.recoveryBindingWork,
                                       &mixedFacetBindingReceipt, mixedFacetBindingStateRevision)
                                 : recoverTetFacetIndexed(assessments[selected].request, state,
                                                          workspace, context);
            if (mixedPresent) {
                mixedFacetBindingReceipt.invalidate();
                mixedFacetBindingStateRevision = nextBindingStateRevision;
            }
            syncBindingWorkCompleteness();
            if (mixedPresent && recovered.recovery.status == TetFacetRecoveryStatus::Recovered)
                mixedPresent->metrics.nontrivialRecoveries = preparedNontrivialRecoveries;
            if (mixedPresent &&
                !addMixedGeometryCounters(mixedPresent->metrics.recoveryCounters,
                                          recovered.recovery.counters, mixedPresent))
                mixedPresent->metrics.complete = false;
            TetFacetSetRecoveryRow row;
            row.kind = TetFacetSetWorkKind::Facet;
            row.facetKey = facets[selected]->key;
            row.facetOutcome = recovered.recovery;
            out.processed.push_back(row);
            if (!facetReady(row.facetOutcome)) {
                pending[selected] = false;
                if (mixedPresent)
                    facetDone[selected] = true;
                --remaining;
                for (std::size_t i = 0; i < facets.size(); ++i)
                    if (mixedPresent ? !facetDone[i] : pending[i])
                        out.unattemptedFacets.push_back(facets[i]->key);
                out.status = row.facetOutcome.status == TetFacetRecoveryStatus::Rejected
                                 ? TetFacetSetRecoveryStatus::Rejected
                             : row.facetOutcome.status == TetFacetRecoveryStatus::Refused
                                 ? TetFacetSetRecoveryStatus::Refused
                                 : TetFacetSetRecoveryStatus::Failed;
                out.reportComplete = true;
                return out;
            }
            pending[selected] = false;
            if (mixedPresent) {
                facetDone[selected] = true;
                const auto binding = state.bindings.find(assessments[selected].request.key);
                if (binding == state.bindings.end())
                    throw std::logic_error("mixed recovered facet binding");
                mixedFaceReverse.push_back({binding->second.recoveredFace, binding->first});
                std::sort(mixedFaceReverse.begin(), mixedFaceReverse.end());
            }
            --remaining;
            for (std::size_t i = 0; i < facets.size(); ++i) {
                if (!pending[i])
                    continue;
                if (!invalidated[i])
                    invalidated[i] = selectedCavityTouchesTarget(recovered.createdCells,
                                                                 assessments[i].request, context);
                if (!invalidated[i] &&
                    std::find(assessments[i].blockers.begin(), assessments[i].blockers.end(),
                              selected) != assessments[i].blockers.end())
                    invalidated[i] = true;
                if (invalidated[i])
                    assessments[i].current = false;
            }
        }
        out.status = TetFacetSetRecoveryStatus::Complete;
        out.reason = TetFacetSetRecoveryReason::None;
        out.reportComplete = true;
        return out;
    } catch (const std::bad_alloc &) {
        out.status = TetFacetSetRecoveryStatus::Failed;
        out.reason = TetFacetSetRecoveryReason::AllocationFailed;
        if (mixedPresent && reportEnvelopePrepared) {
            itemizeUnprocessed();
            out.reportComplete = true;
        }
        rollbackUncommittedMixedSchedules();
        return out;
    } catch (...) {
        out.status = TetFacetSetRecoveryStatus::Failed;
        out.reason = TetFacetSetRecoveryReason::InternalError;
        if (mixedPresent && reportEnvelopePrepared) {
            itemizeUnprocessed();
            out.reportComplete = true;
        }
        rollbackUncommittedMixedSchedules();
        return out;
    }
}

TetFacetSetRecoveryOutcome recoverTetFacetSetIndexed(
    const ClosedTri3Artifact& artifact, const std::map<Index,topo::NodeHandle>& localByGlobal,
    mesh::Classification temporaryRegion, TetFacetSetState& state,
    topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    TetFacetSetFrontierSink* frontier, TetFacetConflictCycleSink* conflictCycle,
    TetFacetSetMixedPresentBatchSink* mixedPresent) noexcept {
    return recoverTetFacetSetIndexedImpl(artifact,localByGlobal,temporaryRegion,state,
        workspace,context,frontier,conflictCycle,mixedPresent,nullptr);
}

TetFacetSetRecoveryOutcome recoverTetFacetSetIndexed(
    const ClosedTri3Artifact& artifact, const std::map<Index,topo::NodeHandle>& localByGlobal,
    mesh::Classification temporaryRegion, TetFacetSetState& state,
    topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    const TetEdgeGeneralRecoveryPolicy& generalRecovery,
    TetFacetSetFrontierSink* frontier, TetFacetConflictCycleSink* conflictCycle,
    TetFacetSetMixedPresentBatchSink* mixedPresent) noexcept {
    return recoverTetFacetSetIndexedImpl(artifact,localByGlobal,temporaryRegion,state,
        workspace,context,frontier,conflictCycle,mixedPresent,&generalRecovery);
}

#if defined(PAACS_TOPO_TESTING)
namespace {

thread_local T2CEdgeLineageSinkForTesting *gT2CEdgeLineageSink = nullptr;

struct LineageCellSnapshot {
    topo::CellHandle handle{};
    std::uint64_t sequence = 0;
    std::array<DenseTetNodeId, 4> nodes{};
    mesh::Classification classification{};
};

struct LineageEdgeSnapshot {
    DenseTetEdgeKey key{};
    topo::EdgeHandle handle{};
    std::vector<topo::CellHandle> completeStar;
    mesh::Classification classification{};
    bool required = false;
};

struct LineageFaceSnapshot {
    DenseTetFaceKey key{};
    topo::FaceHandle handle{};
    std::vector<FaceResolverIncidence> owners;
    mesh::Classification classification{};
};

struct LineageBlockerCertificate {
    TetEdgeRecoveryRequest request{};
    std::array<DenseTetNodeId, 2> targetDense{};
    std::array<Vec3, 2> targetPositions{};
    std::array<mesh::Classification, 2> targetClassifications{};
    std::vector<topo::CellHandle> selectedCells;
    std::vector<topo::CellHandle> supportCells;
    std::vector<LineageCellSnapshot> cells;
    std::vector<LineageEdgeSnapshot> edges;
    std::vector<LineageFaceSnapshot> faces;
    std::uint64_t digest = 0;
    std::size_t sinkRow = (std::numeric_limits<std::size_t>::max)();
};

struct LineageAssessment {
    TetEdgeRecoveryRequest request{};
    TetEdgeRecoveryOutcome outcome{};
    std::vector<topo::CellHandle> selectedCells;
    bool mapped = false;
    bool excluded = false;
};

bool lineageClassEqual(const mesh::Classification &a, const mesh::Classification &b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}

bool lineageVecEqual(const Vec3 &a, const Vec3 &b) noexcept {
    return std::bit_cast<std::uint64_t>(a.x) == std::bit_cast<std::uint64_t>(b.x) &&
           std::bit_cast<std::uint64_t>(a.y) == std::bit_cast<std::uint64_t>(b.y) &&
           std::bit_cast<std::uint64_t>(a.z) == std::bit_cast<std::uint64_t>(b.z);
}

bool lineageCellEqual(const LineageCellSnapshot &a, const LineageCellSnapshot &b) noexcept {
    return a.handle == b.handle && a.sequence == b.sequence && a.nodes == b.nodes &&
           lineageClassEqual(a.classification, b.classification);
}

bool lineageIncidenceEqual(const FaceResolverIncidence &a,
                           const FaceResolverIncidence &b) noexcept {
    return a.cell == b.cell && a.sense == b.sense;
}

bool lineageEdgeEqual(const LineageEdgeSnapshot &a, const LineageEdgeSnapshot &b) noexcept {
    return a.key == b.key && a.handle == b.handle && a.completeStar == b.completeStar &&
           lineageClassEqual(a.classification, b.classification) && a.required == b.required;
}

bool lineageFaceEqual(const LineageFaceSnapshot &a, const LineageFaceSnapshot &b) noexcept {
    if (!(a.key == b.key) || a.handle != b.handle ||
        !lineageClassEqual(a.classification, b.classification) ||
        a.owners.size() != b.owners.size())
        return false;
    for (std::size_t i = 0; i < a.owners.size(); ++i)
        if (!lineageIncidenceEqual(a.owners[i], b.owners[i]))
            return false;
    return true;
}

void lineageMix(std::uint64_t &digest, std::uint64_t value) noexcept {
    digest ^= value + UINT64_C(0x9e3779b97f4a7c15) + (digest << 6) + (digest >> 2);
}

void lineageMixHandle(std::uint64_t &digest, topo::CellHandle handle) noexcept {
    lineageMix(digest, handle.slot);
    lineageMix(digest, handle.gen);
    lineageMix(digest, handle.pool);
}

template <class Handle> void lineageMixEntityHandle(std::uint64_t &digest, Handle handle) noexcept {
    lineageMix(digest, handle.slot);
    lineageMix(digest, handle.gen);
    lineageMix(digest, handle.pool);
}

bool lineageAdd(std::uint64_t &target, std::uint64_t value) noexcept {
    if (value > (std::numeric_limits<std::uint64_t>::max)() - target)
        return false;
    target += value;
    return true;
}

bool lineageAddSize(std::size_t &target, std::size_t value) noexcept {
    if (value > (std::numeric_limits<std::size_t>::max)() - target)
        return false;
    target += value;
    return true;
}

bool lineageAddWork(T2CEdgeLineageSinkForTesting &sink, T2CEdgeLineageAttemptForTesting &attempt,
                    std::uint64_t coordinator, std::uint64_t nested) noexcept {
    if (!lineageAdd(attempt.coordinatorWork, coordinator) ||
        !lineageAdd(attempt.nestedRecoveryWork, nested) ||
        !lineageAdd(sink.coordinatorWork, coordinator) ||
        !lineageAdd(sink.nestedRecoveryWork, nested) ||
        attempt.coordinatorWork > sink.limits.maxCoordinatorWork ||
        attempt.nestedRecoveryWork > sink.limits.maxNestedRecoveryWork ||
        sink.coordinatorWork > sink.limits.maxCoordinatorWork ||
        sink.nestedRecoveryWork > sink.limits.maxNestedRecoveryWork) {
        sink.status = T2CEdgeLineageStatusForTesting::WorkBudgetReached;
        sink.complete = false;
        return false;
    }
    return true;
}

std::uint64_t lineageBindingWork(const TetRecoveryBindingWorkCounters &work,
                                 bool &complete) noexcept {
    complete = complete && work.complete;
    std::uint64_t total = 0;
    const std::array<std::uint64_t, 16> values{
        {work.edgeBindingAuditCalls, work.edgeBindingRowsVisited, work.edgeAliasRowsVisited,
         work.facetBindingAuditCalls, work.facetBindingRowsVisited,
         work.facetBindingRequestRowsVisited, work.facetBindingDuplicateRowsCompared,
         work.facetBindingReceiptBuilds, work.facetBindingReceiptReuses,
         work.facetBindingReceiptRowsCopied, work.facetBindingReceiptLookups,
         work.facetBindingReceiptReserveCalls, work.edgeDataCopyRows, work.faceDataCopyRows,
         work.faceIncidenceRowsCopied, work.stateCopyRows}};
    for (const auto value : values)
        complete = complete && lineageAdd(total, value);
    complete = complete && lineageAdd(total, work.stateCopies);
    return total;
}

std::uint64_t lineageContextDelta(const TetRecoveryContextMetrics &before,
                                  const TetRecoveryContextMetrics &after, bool &complete) noexcept {
    std::uint64_t total = 0;
    const auto addDelta = [&](std::uint64_t a, std::uint64_t b) {
        if (b < a || !lineageAdd(total, b - a))
            complete = false;
    };
    addDelta(before.auditVisits, after.auditVisits);
    addDelta(before.segmentQueries, after.segmentQueries);
    addDelta(before.facetQueries, after.facetQueries);
    addDelta(before.boundsQueries, after.boundsQueries);
    addDelta(before.indexNodesVisited, after.indexNodesVisited);
    addDelta(before.returnedCandidates, after.returnedCandidates);
    addDelta(before.indexRemovals, after.indexRemovals);
    addDelta(before.indexInsertions, after.indexInsertions);
    addDelta(before.indexRebuilds, after.indexRebuilds);
    return total;
}

template <class T>
bool lineageReserve(std::vector<T> &values, std::size_t capacity, std::size_t &accounted) {
    if (values.capacity() < capacity)
        values.reserve(capacity);
    if (values.capacity() != 0 &&
        sizeof(T) > (std::numeric_limits<std::size_t>::max)() / values.capacity())
        return false;
    return lineageAddSize(accounted, values.capacity() * sizeof(T));
}

bool prepareLineageStorage(T2CEdgeLineageSinkForTesting &sink) {
    if (!sink.limits.positive())
        return false;
    std::size_t accounted = 0;
    if (!lineageReserve(sink.attempts, sink.limits.maxAttempts, accounted) ||
        !lineageReserve(sink.blockers, sink.limits.maxBlockers, accounted) ||
        !lineageReserve(sink.recoveries, sink.limits.maxScratchRecoveries, accounted) ||
        !lineageReserve(sink.cellRows, sink.limits.maxCellRows, accounted) ||
        !lineageReserve(sink.supportCellRows, sink.limits.maxSupportCellRows, accounted) ||
        !lineageReserve(sink.vertexRows, sink.limits.maxVertexRows, accounted) ||
        !lineageReserve(sink.edgeRows, sink.limits.maxEdgeRows, accounted) ||
        !lineageReserve(sink.faceRows, sink.limits.maxFaceRows, accounted) ||
        !lineageReserve(sink.faceOwnerRows, sink.limits.maxFaceOwnerRows, accounted) ||
        !lineageReserve(sink.candidateRows, sink.limits.maxCandidateRows, accounted) ||
        !lineageReserve(sink.boundEdgeRows, sink.limits.maxBoundEdgeRows, accounted) ||
        !lineageReserve(sink.passRows, sink.limits.maxPassRows, accounted) ||
        !lineageReserve(sink.assessmentRows, sink.limits.maxAssessmentRows, accounted) ||
        !lineageReserve(sink.decisionRows, sink.limits.maxDecisionRows, accounted) ||
        !lineageReserve(sink.blockerIndexRows, sink.limits.maxBlockerIndexRows, accounted) ||
        accounted > sink.limits.maxAccountedBytes)
        return false;
    sink.accountedByteHighWater = (std::max)(sink.accountedByteHighWater, accounted);
    return true;
}

bool lineageHandleLess(topo::CellHandle a, topo::CellHandle b) noexcept {
    return std::tie(a.pool, a.slot, a.gen) < std::tie(b.pool, b.slot, b.gen);
}

bool buildLineageCertificate(const TetEdgeRecoveryRequest &request,
                             std::span<const topo::CellHandle> selected,
                             const TetRecoveryContext &context,
                             const topo::TopologyWorkspace &workspace,
                             LineageBlockerCertificate &result, std::uint64_t &work) {
    result = {};
    result.request = request;
    for (std::size_t i = 0; i < 2; ++i) {
        const auto found = context.nodeToDense.find(request.vertices[i]);
        if (found == context.nodeToDense.end() || found->second >= context.plainVertices.size())
            return false;
        result.targetDense[i] = found->second;
        result.targetPositions[i] = context.plainVertices[found->second].position;
        result.targetClassifications[i] = workspace.nodes().get(request.vertices[i]).cls;
    }
    result.selectedCells.assign(selected.begin(), selected.end());
    std::sort(result.selectedCells.begin(), result.selectedCells.end(),
              [&](const auto a, const auto b) {
                  const auto ia = context.cells.find(a), ib = context.cells.find(b);
                  if (ia == context.cells.end() || ib == context.cells.end())
                      return lineageHandleLess(a, b);
                  return std::tie(ia->second.sequence, a) < std::tie(ib->second.sequence, b);
              });
    if (std::adjacent_find(result.selectedCells.begin(), result.selectedCells.end()) !=
        result.selectedCells.end())
        return false;
    std::set<DenseTetEdgeKey> edgeKeys;
    std::set<DenseTetFaceKey> faceKeys;
    std::set<topo::CellHandle, decltype(&lineageHandleLess)> support(&lineageHandleLess);
    for (const auto handle : result.selectedCells) {
        const auto found = context.cells.find(handle);
        if (found == context.cells.end())
            return false;
        const auto &cell = found->second;
        support.insert(handle);
        for (const auto &edge : kEdges)
            edgeKeys.insert(DenseTetEdgeKey::of(cell.nodes[static_cast<std::size_t>(edge[0])],
                                                cell.nodes[static_cast<std::size_t>(edge[1])]));
        for (const auto &face : kFaces)
            faceKeys.insert(DenseTetFaceKey::of({cell.nodes[static_cast<std::size_t>(face[0])],
                                                 cell.nodes[static_cast<std::size_t>(face[1])],
                                                 cell.nodes[static_cast<std::size_t>(face[2])]}));
    }
    for (const auto &key : edgeKeys) {
        const auto found = context.edges.find(key);
        if (found == context.edges.end())
            return false;
        LineageEdgeSnapshot row;
        row.key = key;
        row.handle = found->second.handle;
        row.completeStar = found->second.completeStar;
        row.classification = found->second.classification;
        row.required = found->second.required;
        std::sort(row.completeStar.begin(), row.completeStar.end(), lineageHandleLess);
        for (const auto cell : row.completeStar)
            support.insert(cell);
        result.edges.push_back(std::move(row));
    }
    for (const auto &key : faceKeys) {
        const auto found = context.faces.find(key);
        if (found == context.faces.end())
            return false;
        LineageFaceSnapshot row;
        row.key = key;
        row.handle = found->second.handle;
        row.owners = found->second.owners;
        row.classification = found->second.classification;
        std::sort(row.owners.begin(), row.owners.end(), [](const auto &a, const auto &b) {
            return std::tie(a.cell, a.sense) < std::tie(b.cell, b.sense);
        });
        for (const auto &owner : row.owners)
            support.insert(owner.cell);
        result.faces.push_back(std::move(row));
    }
    result.supportCells.assign(support.begin(), support.end());
    // Complete contact authority: every selected cell plus every cell whose
    // edge-star or face-owner row was read is retained with its exact current
    // payload.  Selected-only payloads would miss a retired neighbour that
    // changes a frozen incidence row.
    result.cells.reserve(result.supportCells.size());
    for (const auto handle : result.supportCells) {
        const auto found = context.cells.find(handle);
        if (found == context.cells.end())
            return false;
        result.cells.push_back(
            {handle, found->second.sequence, found->second.nodes, found->second.classification});
    }
    if (!lineageAdd(work, result.cells.size()) || !lineageAdd(work, result.edges.size()) ||
        !lineageAdd(work, result.faces.size()) || !lineageAdd(work, result.supportCells.size()))
        return false;

    // Framed support certificate.  Workspace revision is deliberately not a
    // member: it advances after an independent scratch commit and is checked
    // separately through currentness.  Artifact/request identity and every
    // exact support row read by the decision are members.
    std::uint64_t digest = UINT64_C(0xcbf29ce484222325);
    lineageMix(digest, UINT64_C(0x4c494e4541474533));
    lineageMix(digest, request.artifactStamp.persistentId);
    lineageMix(digest, request.artifactStamp.revision);
    lineageMix(digest, request.key.canonicalGlobalNodes[0]);
    lineageMix(digest, request.key.canonicalGlobalNodes[1]);
    lineageMix(digest, static_cast<std::uint64_t>(request.classification.dim));
    lineageMix(digest, static_cast<std::uint64_t>(request.classification.tag));
    lineageMix(digest, static_cast<std::uint64_t>(request.region.dim));
    lineageMix(digest, static_cast<std::uint64_t>(request.region.tag));
    for (std::size_t i = 0; i < 2; ++i) {
        lineageMixEntityHandle(digest, request.vertices[i]);
        lineageMix(digest, result.targetDense[i]);
        lineageMix(digest, std::bit_cast<std::uint64_t>(result.targetPositions[i].x));
        lineageMix(digest, std::bit_cast<std::uint64_t>(result.targetPositions[i].y));
        lineageMix(digest, std::bit_cast<std::uint64_t>(result.targetPositions[i].z));
        lineageMix(digest, static_cast<std::uint64_t>(result.targetClassifications[i].dim));
        lineageMix(digest, static_cast<std::uint64_t>(result.targetClassifications[i].tag));
    }
    lineageMix(digest, UINT64_C(0x53454c4543544544));
    lineageMix(digest, result.selectedCells.size());
    for (const auto cell : result.selectedCells)
        lineageMixHandle(digest, cell);
    lineageMix(digest, UINT64_C(0x535550504f525443));
    lineageMix(digest, result.supportCells.size());
    for (const auto cell : result.supportCells)
        lineageMixHandle(digest, cell);
    lineageMix(digest, UINT64_C(0x43454c4c5041594c));
    lineageMix(digest, result.cells.size());
    for (const auto &cell : result.cells) {
        lineageMixHandle(digest, cell.handle);
        lineageMix(digest, cell.sequence);
        for (const auto node : cell.nodes)
            lineageMix(digest, node);
        lineageMix(digest, static_cast<std::uint64_t>(cell.classification.dim));
        lineageMix(digest, static_cast<std::uint64_t>(cell.classification.tag));
    }
    lineageMix(digest, UINT64_C(0x45444745524f5753));
    lineageMix(digest, result.edges.size());
    for (const auto &edge : result.edges) {
        lineageMix(digest, edge.key.nodes[0]);
        lineageMix(digest, edge.key.nodes[1]);
        lineageMixEntityHandle(digest, edge.handle);
        lineageMix(digest, edge.required ? 1U : 0U);
        lineageMix(digest, static_cast<std::uint64_t>(edge.classification.dim));
        lineageMix(digest, static_cast<std::uint64_t>(edge.classification.tag));
        lineageMix(digest, edge.completeStar.size());
        for (const auto cell : edge.completeStar)
            lineageMixHandle(digest, cell);
    }
    lineageMix(digest, UINT64_C(0x46414345524f5753));
    lineageMix(digest, result.faces.size());
    for (const auto &face : result.faces) {
        for (const auto node : face.key.nodes)
            lineageMix(digest, node);
        lineageMixEntityHandle(digest, face.handle);
        lineageMix(digest, static_cast<std::uint64_t>(face.classification.dim));
        lineageMix(digest, static_cast<std::uint64_t>(face.classification.tag));
        lineageMix(digest, face.owners.size());
        for (const auto &owner : face.owners) {
            lineageMixHandle(digest, owner.cell);
            lineageMix(digest, static_cast<std::uint64_t>(static_cast<std::int64_t>(owner.sense)));
        }
    }
    result.digest = digest;
    return true;
}

bool lineageCertificateEqual(const LineageBlockerCertificate &a,
                             const LineageBlockerCertificate &b) noexcept {
    if (!(a.request.key == b.request.key) || a.request.vertices != b.request.vertices ||
        !lineageClassEqual(a.request.classification, b.request.classification) ||
        !lineageClassEqual(a.request.region, b.request.region) || a.targetDense != b.targetDense ||
        !lineageVecEqual(a.targetPositions[0], b.targetPositions[0]) ||
        !lineageVecEqual(a.targetPositions[1], b.targetPositions[1]) ||
        !lineageClassEqual(a.targetClassifications[0], b.targetClassifications[0]) ||
        !lineageClassEqual(a.targetClassifications[1], b.targetClassifications[1]) ||
        a.selectedCells != b.selectedCells || a.supportCells != b.supportCells ||
        a.cells.size() != b.cells.size() || a.edges.size() != b.edges.size() ||
        a.faces.size() != b.faces.size() || a.digest != b.digest)
        return false;
    for (std::size_t i = 0; i < a.cells.size(); ++i)
        if (!lineageCellEqual(a.cells[i], b.cells[i]))
            return false;
    for (std::size_t i = 0; i < a.edges.size(); ++i)
        if (!lineageEdgeEqual(a.edges[i], b.edges[i]))
            return false;
    for (std::size_t i = 0; i < a.faces.size(); ++i)
        if (!lineageFaceEqual(a.faces[i], b.faces[i]))
            return false;
    return true;
}

bool lineageCertificateCurrent(const LineageBlockerCertificate &certificate,
                               const topo::TopologyWorkspace &workspace,
                               const TetRecoveryContext &context) noexcept {
    if (!context.currentFor(workspace))
        return false;
    for (const auto node : certificate.request.vertices)
        if (!workspace.nodes().isCurrent(node))
            return false;
    for (const auto cell : certificate.supportCells)
        if (!workspace.cells().isCurrent(cell) || context.cells.find(cell) == context.cells.end())
            return false;
    for (const auto &edge : certificate.edges)
        if (!workspace.edges().isCurrent(edge.handle))
            return false;
    for (const auto &face : certificate.faces)
        if (!workspace.faces().isCurrent(face.handle))
            return false;
    return true;
}

template <class T, class Compare>
bool lineageIntersects(const std::vector<T> &a, const std::vector<T> &b, Compare compare) {
    auto ia = a.begin(), ib = b.begin();
    while (ia != a.end() && ib != b.end()) {
        if (compare(*ia, *ib))
            ++ia;
        else if (compare(*ib, *ia))
            ++ib;
        else
            return true;
    }
    return false;
}

bool lineageEndpointWedge(const Vec3 &endpoint, const Vec3 &other, const std::array<Vec3, 4> &tet,
                          const std::array<DenseTetNodeId, 4> &nodes, DenseTetNodeId endpointId) {
    if (std::find(nodes.begin(), nodes.end(), endpointId) == nodes.end())
        return false;
    bool constrained = false;
    for (const auto &face : kFaces) {
        const auto atEndpoint = predicates::orient3d(
            tet[static_cast<std::size_t>(face[0])], tet[static_cast<std::size_t>(face[1])],
            tet[static_cast<std::size_t>(face[2])], endpoint);
        if (atEndpoint != predicates::Sign::Zero)
            continue;
        constrained = true;
        if (predicates::orient3d(
                tet[static_cast<std::size_t>(face[0])], tet[static_cast<std::size_t>(face[1])],
                tet[static_cast<std::size_t>(face[2])], other) == predicates::Sign::Positive)
            return false;
    }
    return constrained;
}

// Exact plain-value equivalent of TetEdgeRecovery.cpp's
// segmentIntersectsCellOpen(), including its endpoint-wedge branch.  This is
// intentionally not the facet scheduler's weaker open-Tet contact predicate.
bool lineageCellTouchesTarget(const LineageCellSnapshot &cell,
                              const LineageBlockerCertificate &target,
                              const TetRecoveryContext &context) {
    std::array<Vec3, 4> tet{};
    for (std::size_t i = 0; i < cell.nodes.size(); ++i) {
        if (cell.nodes[i] >= context.plainVertices.size())
            return true;
        tet[i] = context.plainVertices[cell.nodes[i]].position;
    }
    int crossings = 0;
    for (const auto &face : kFaces)
        crossings +=
            tet_facet_detail::segmentCrossesOpenTriangle(
                target.targetPositions[0], target.targetPositions[1],
                tet[static_cast<std::size_t>(face[0])], tet[static_cast<std::size_t>(face[1])],
                tet[static_cast<std::size_t>(face[2])])
                ? 1
                : 0;
    if (crossings >= 2)
        return true;
    return lineageEndpointWedge(target.targetPositions[0], target.targetPositions[1], tet,
                                cell.nodes, target.targetDense[0]) ||
           lineageEndpointWedge(target.targetPositions[1], target.targetPositions[0], tet,
                                cell.nodes, target.targetDense[1]);
}

bool lineageCellTouchesSupport(const LineageCellSnapshot &cell,
                               const LineageBlockerCertificate &blocker,
                               const TetRecoveryContext &context) {
    if (lineageCellTouchesTarget(cell, blocker, context))
        return true;
    for (const auto &edge : kEdges) {
        const auto key = DenseTetEdgeKey::of(cell.nodes[static_cast<std::size_t>(edge[0])],
                                             cell.nodes[static_cast<std::size_t>(edge[1])]);
        const auto found =
            std::lower_bound(blocker.edges.begin(), blocker.edges.end(), key,
                             [](const auto &row, const auto &value) { return row.key < value; });
        if (found != blocker.edges.end() && found->key == key)
            return true;
    }
    for (const auto &face : kFaces) {
        const auto key = DenseTetFaceKey::of({cell.nodes[static_cast<std::size_t>(face[0])],
                                              cell.nodes[static_cast<std::size_t>(face[1])],
                                              cell.nodes[static_cast<std::size_t>(face[2])]});
        const auto found =
            std::lower_bound(blocker.faces.begin(), blocker.faces.end(), key,
                             [](const auto &row, const auto &value) { return row.key < value; });
        if (found != blocker.faces.end() && found->key == key)
            return true;
    }
    return false;
}

bool lineageBlockersIndependent(const LineageBlockerCertificate &a,
                                const LineageBlockerCertificate &b,
                                const TetRecoveryContext &context, std::uint64_t &work) {
    if (!lineageAdd(work, 1))
        return false;
    if (lineageIntersects(a.supportCells, b.supportCells, lineageHandleLess) ||
        lineageIntersects(a.edges, b.edges,
                          [](const auto &x, const auto &y) { return x.key < y.key; }) ||
        lineageIntersects(a.faces, b.faces,
                          [](const auto &x, const auto &y) { return x.key < y.key; }))
        return false;
    for (const auto &cell : a.cells) {
        if (!lineageAdd(work, 1) || lineageCellTouchesSupport(cell, b, context))
            return false;
    }
    for (const auto &cell : b.cells) {
        if (!lineageAdd(work, 1) || lineageCellTouchesSupport(cell, a, context))
            return false;
    }
    return true;
}

bool appendLineageBlocker(T2CEdgeLineageSinkForTesting &sink,
                          LineageBlockerCertificate &certificate, topo::WorkspaceStamp stamp,
                          TetEdgeRecoveryOutcome outcome, bool deferred) {
    const auto &limits = sink.limits;
    std::size_t edgeStarRows = 0, faceOwnerRows = 0;
    for (const auto &edge : certificate.edges)
        if (!lineageAddSize(edgeStarRows, edge.completeStar.size()))
            return false;
    for (const auto &face : certificate.faces)
        if (!lineageAddSize(faceOwnerRows, face.owners.size()))
            return false;
    if (sink.blockers.size() >= limits.maxBlockers ||
        certificate.selectedCells.size() > limits.maxCellRows - sink.cellRows.size() ||
        certificate.supportCells.size() >
            limits.maxCellRows - sink.cellRows.size() - certificate.selectedCells.size() ||
        edgeStarRows > limits.maxCellRows - sink.cellRows.size() -
                           certificate.selectedCells.size() - certificate.supportCells.size() ||
        certificate.cells.size() > limits.maxSupportCellRows - sink.supportCellRows.size() ||
        certificate.edges.size() > limits.maxEdgeRows - sink.edgeRows.size() ||
        certificate.faces.size() > limits.maxFaceRows - sink.faceRows.size() ||
        faceOwnerRows > limits.maxFaceOwnerRows - sink.faceOwnerRows.size())
        return false;
    T2CEdgeLineageBlockerForTesting row;
    row.key = certificate.request.key;
    row.observedStamp = stamp;
    row.finalValidatedStamp = stamp;
    row.artifactStamp = certificate.request.artifactStamp;
    row.classification = certificate.request.classification;
    row.region = certificate.request.region;
    row.targetDense = certificate.targetDense;
    row.outcome = outcome;
    row.selectedCellOffset = sink.cellRows.size();
    row.selectedCellCount = certificate.selectedCells.size();
    sink.cellRows.insert(sink.cellRows.end(), certificate.selectedCells.begin(),
                         certificate.selectedCells.end());
    row.supportCellOffset = sink.cellRows.size();
    row.supportCellCount = certificate.supportCells.size();
    sink.cellRows.insert(sink.cellRows.end(), certificate.supportCells.begin(),
                         certificate.supportCells.end());
    row.supportPayloadOffset = sink.supportCellRows.size();
    row.supportPayloadCount = certificate.cells.size();
    for (const auto &cell : certificate.cells)
        sink.supportCellRows.push_back(
            {cell.handle, cell.sequence, cell.nodes, cell.classification});
    row.supportEdgeOffset = sink.edgeRows.size();
    row.supportEdgeCount = certificate.edges.size();
    for (const auto &edge : certificate.edges) {
        T2CEdgeLineageSupportEdgeForTesting retained;
        retained.key = edge.key;
        retained.handle = edge.handle;
        retained.classification = edge.classification;
        retained.required = edge.required;
        retained.completeStarOffset = sink.cellRows.size();
        retained.completeStarCount = edge.completeStar.size();
        sink.cellRows.insert(sink.cellRows.end(), edge.completeStar.begin(),
                             edge.completeStar.end());
        sink.edgeRows.push_back(retained);
    }
    row.supportFaceOffset = sink.faceRows.size();
    row.supportFaceCount = certificate.faces.size();
    for (const auto &face : certificate.faces) {
        T2CEdgeLineageSupportFaceForTesting retained;
        retained.key = face.key;
        retained.handle = face.handle;
        retained.classification = face.classification;
        retained.ownerOffset = sink.faceOwnerRows.size();
        retained.ownerCount = face.owners.size();
        for (const auto &owner : face.owners)
            sink.faceOwnerRows.push_back({owner.cell, owner.sense});
        sink.faceRows.push_back(retained);
    }
    row.supportDigest = certificate.digest;
    row.deferred = deferred;
    certificate.sinkRow = sink.blockers.size();
    sink.blockers.push_back(row);
    return true;
}

bool lineagePreviewTouchesBlocker(const topo::TetEdgeCavityPreview &preview,
                                  const LineageBlockerCertificate &blocker,
                                  const TetRecoveryContext &context,
                                  std::vector<T2CEdgeLineageCandidateCellForTesting> &candidateRows,
                                  std::size_t candidateOffset, std::uint64_t &work,
                                  bool &retiredTouches, bool &candidateTouches) noexcept {
    retiredTouches = false;
    candidateTouches = false;
    for (const auto retired : preview.removedCells) {
        if (!lineageAdd(work, 1)) {
            retiredTouches = true;
            continue;
        }
        if (std::binary_search(blocker.supportCells.begin(), blocker.supportCells.end(), retired,
                               lineageHandleLess))
            retiredTouches = true;
        const auto found = context.cells.find(retired);
        if (found == context.cells.end()) {
            retiredTouches = true;
            continue;
        }
        const LineageCellSnapshot retiredRow{retired, found->second.sequence, found->second.nodes,
                                             found->second.classification};
        if (lineageCellTouchesSupport(retiredRow, blocker, context))
            retiredTouches = true;
    }
    for (std::size_t i = 0; i < preview.addedConnectivity.size(); ++i) {
        if (candidateOffset + i >= candidateRows.size()) {
            candidateTouches = true;
            continue;
        }
        const auto &candidate = candidateRows[candidateOffset + i];
        if (!lineageAdd(work, 1)) {
            candidateTouches = true;
            continue;
        }
        LineageCellSnapshot plain;
        plain.nodes = candidate.denseNodes;
        plain.classification = candidate.classification;
        if (lineageCellTouchesSupport(plain, blocker, context))
            candidateTouches = true;
    }
    return retiredTouches || candidateTouches;
}

bool appendLineageDecision(T2CEdgeLineageSinkForTesting &sink,
                           T2CEdgeLineageDependencyKindForTesting kind, std::size_t sourceRow,
                           std::size_t targetRow, bool touches, bool retiredTouches = false,
                           bool candidateTouches = false);

struct LineagePreparedGuardState {
    T2CEdgeLineageSinkForTesting *sink = nullptr;
    T2CEdgeLineageAttemptForTesting *attempt = nullptr;
    const std::vector<LineageBlockerCertificate> *blockers = nullptr;
    const std::vector<std::pair<topo::EdgeHandle, VolumeArtifactEdgeKey>> *boundEdges = nullptr;
    T2CEdgeLineageRecoveryForTesting *recovery = nullptr;
    std::size_t recoveryRow = 0;
    bool invoked = false;
    bool capacityFailed = false;
    bool workFailed = false;
};

TetEdgePreparedGuardDecisionForTesting
lineagePreparedGuard(const TetEdgeRecoveryRequest &request, topo::WorkspaceStamp stamp,
                     topo::EdgeHandle predicted, const topo::TetEdgeCavityPreview &preview,
                     const TetRecoveryContext &context, void *opaque) noexcept {
    auto &guard = *static_cast<LineagePreparedGuardState *>(opaque);
    guard.invoked = true;
    auto &sink = *guard.sink;
    auto &recovery = *guard.recovery;
    recovery.beforeStamp = stamp;
    recovery.predictedEdge = predicted;
    recovery.retiredCellOffset = sink.cellRows.size();
    recovery.retiredCellCount = preview.removedCells.size();
    if (preview.removedCells.size() > sink.limits.maxCellRows - sink.cellRows.size()) {
        guard.capacityFailed = true;
        return TetEdgePreparedGuardDecisionForTesting::Defer;
    }
    sink.cellRows.insert(sink.cellRows.end(), preview.removedCells.begin(),
                         preview.removedCells.end());
    recovery.retiredPayloadOffset = sink.supportCellRows.size();
    recovery.retiredPayloadCount = preview.removedCells.size();
    if (preview.removedCells.size() >
        sink.limits.maxSupportCellRows - sink.supportCellRows.size()) {
        guard.capacityFailed = true;
        return TetEdgePreparedGuardDecisionForTesting::Defer;
    }
    for (const auto retired : preview.removedCells) {
        const auto found = context.cells.find(retired);
        if (found == context.cells.end()) {
            guard.capacityFailed = true;
            return TetEdgePreparedGuardDecisionForTesting::Defer;
        }
        sink.supportCellRows.push_back(
            {retired, found->second.sequence, found->second.nodes, found->second.classification});
    }
    recovery.candidateOffset = sink.candidateRows.size();
    recovery.candidateCount = preview.addedConnectivity.size();
    if (preview.addedConnectivity.size() >
        sink.limits.maxCandidateRows - sink.candidateRows.size()) {
        guard.capacityFailed = true;
        return TetEdgePreparedGuardDecisionForTesting::Defer;
    }
    for (std::size_t i = 0; i < preview.addedConnectivity.size(); ++i) {
        T2CEdgeLineageCandidateCellForTesting candidate;
        for (std::size_t node = 0; node < 4; ++node) {
            const auto found = context.nodeToDense.find(preview.addedConnectivity[i][node]);
            if (found == context.nodeToDense.end()) {
                guard.capacityFailed = true;
                return TetEdgePreparedGuardDecisionForTesting::Defer;
            }
            candidate.denseNodes[node] = found->second;
        }
        candidate.classification = preview.addedClassifications[i];
        sink.candidateRows.push_back(candidate);
    }
    std::uint64_t pairWork = 0;
    bool safe = true;
    for (const auto &blocker : *guard.blockers) {
        if (sink.transactionPairTests == sink.limits.maxPairTests ||
            !lineageAdd(sink.transactionPairTests, 1) ||
            !lineageAdd(recovery.blockerPairTests, 1)) {
            guard.workFailed = true;
            safe = false;
            break;
        }
        bool retiredTouches = false, candidateTouches = false;
        const bool touches = lineagePreviewTouchesBlocker(
            preview, blocker, context, sink.candidateRows, recovery.candidateOffset, pairWork,
            retiredTouches, candidateTouches);
        if (!appendLineageDecision(
                sink, T2CEdgeLineageDependencyKindForTesting::PreparedDeltaToFrozenBlocker,
                guard.recoveryRow, blocker.sinkRow, touches, retiredTouches, candidateTouches)) {
            guard.capacityFailed = true;
            safe = false;
            break;
        }
        if (touches)
            safe = false;
    }
    for (const auto &[handle, key] : *guard.boundEdges) {
        if (handle == predicted && !(key == request.key)) {
            recovery.predictedBindingAlias = true;
            safe = false;
        }
    }
    if (!lineageAddWork(sink, *guard.attempt, pairWork, 0)) {
        guard.workFailed = true;
        safe = false;
    }
    recovery.allowed = safe;
    return safe ? TetEdgePreparedGuardDecisionForTesting::Allow
                : TetEdgePreparedGuardDecisionForTesting::Defer;
}

bool appendCreatedRows(T2CEdgeLineageSinkForTesting &sink, T2CEdgeLineageRecoveryForTesting &row,
                       std::span<const topo::CellHandle> created) {
    if (created.size() > sink.limits.maxCellRows - sink.cellRows.size())
        return false;
    row.createdCellOffset = sink.cellRows.size();
    row.createdCellCount = created.size();
    sink.cellRows.insert(sink.cellRows.end(), created.begin(), created.end());
    return true;
}

bool lineageSelectedTouchesRequest(std::span<const topo::CellHandle> selected,
                                   const TetEdgeRecoveryRequest &target,
                                   const TetRecoveryContext &context) {
    LineageBlockerCertificate plainTarget;
    plainTarget.request = target;
    for (std::size_t i = 0; i < 2; ++i) {
        const auto found = context.nodeToDense.find(target.vertices[i]);
        if (found == context.nodeToDense.end() || found->second >= context.plainVertices.size())
            return true;
        plainTarget.targetDense[i] = found->second;
        plainTarget.targetPositions[i] = context.plainVertices[found->second].position;
    }
    for (const auto handle : selected) {
        const auto found = context.cells.find(handle);
        if (found == context.cells.end())
            return true;
        const LineageCellSnapshot row{handle, found->second.sequence, found->second.nodes,
                                      found->second.classification};
        if (lineageCellTouchesTarget(row, plainTarget, context))
            return true;
    }
    return false;
}

bool appendLineagePass(T2CEdgeLineageSinkForTesting &sink,
                       std::span<const VolumeArtifactEdge *const> edges,
                       const std::vector<LineageAssessment> &assessments,
                       const std::vector<LineageBlockerCertificate> &frozen,
                       const std::vector<bool> &deferredAtEntry, const TetRecoveryContext &context,
                       std::vector<std::size_t> &retainedAssessmentRows,
                       std::size_t &retainedPassRow) {
    if (sink.passRows.size() >= sink.limits.maxPassRows || assessments.size() != edges.size() ||
        deferredAtEntry.size() != edges.size() ||
        assessments.size() > sink.limits.maxAssessmentRows - sink.assessmentRows.size() ||
        frozen.size() > sink.limits.maxBlockerIndexRows - sink.blockerIndexRows.size())
        return false;
    std::size_t selectedRows = 0;
    for (const auto &assessment : assessments)
        if (!lineageAddSize(selectedRows, assessment.selectedCells.size()))
            return false;
    if (selectedRows > sink.limits.maxCellRows - sink.cellRows.size() ||
        selectedRows > sink.limits.maxSupportCellRows - sink.supportCellRows.size())
        return false;
    retainedAssessmentRows.clear();
    retainedAssessmentRows.reserve(assessments.size());
    T2CEdgeLineagePassForTesting pass;
    pass.stamp = context.stamp;
    pass.assessmentOffset = sink.assessmentRows.size();
    pass.assessmentCount = assessments.size();
    pass.decisionOffset = sink.decisionRows.size();
    pass.frozenBlockerOffset = sink.blockerIndexRows.size();
    pass.frozenBlockerCount = frozen.size();
    for (const auto &blocker : frozen) {
        if (blocker.sinkRow >= sink.blockers.size())
            return false;
        sink.blockerIndexRows.push_back(blocker.sinkRow);
    }
    for (std::size_t i = 0; i < assessments.size(); ++i) {
        const auto &assessment = assessments[i];
        T2CEdgeLineageAssessmentForTesting row;
        row.key = edges[i]->key;
        row.outcome = assessment.outcome;
        row.canonicalOrdinal = i;
        row.mapped = assessment.mapped;
        row.excludedFrozen = assessment.excluded;
        row.candidate = assessment.outcome.status == TetEdgeRecoveryStatus::Recovered;
        row.deferredAtEntry = deferredAtEntry[i];
        row.deferred = deferredAtEntry[i];
        if (assessment.mapped) {
            row.artifactStamp = assessment.request.artifactStamp;
            row.classification = assessment.request.classification;
            row.region = assessment.request.region;
            for (std::size_t endpoint = 0; endpoint < 2; ++endpoint) {
                const auto found = context.nodeToDense.find(assessment.request.vertices[endpoint]);
                if (found == context.nodeToDense.end())
                    return false;
                row.targetDense[endpoint] = found->second;
            }
        }
        row.selectedCellOffset = sink.cellRows.size();
        row.selectedCellCount = assessment.selectedCells.size();
        sink.cellRows.insert(sink.cellRows.end(), assessment.selectedCells.begin(),
                             assessment.selectedCells.end());
        row.selectedPayloadOffset = sink.supportCellRows.size();
        row.selectedPayloadCount = assessment.selectedCells.size();
        for (const auto handle : assessment.selectedCells) {
            const auto found = context.cells.find(handle);
            if (found == context.cells.end())
                return false;
            sink.supportCellRows.push_back({handle, found->second.sequence, found->second.nodes,
                                            found->second.classification});
        }
        retainedAssessmentRows.push_back(sink.assessmentRows.size());
        sink.assessmentRows.push_back(row);
    }
    retainedPassRow = sink.passRows.size();
    sink.passRows.push_back(pass);
    return true;
}

bool appendLineageDecision(T2CEdgeLineageSinkForTesting &sink,
                           T2CEdgeLineageDependencyKindForTesting kind, std::size_t sourceRow,
                           std::size_t targetRow, bool touches, bool retiredTouches,
                           bool candidateTouches) {
    if (sink.decisionRows.size() >= sink.limits.maxDecisionRows)
        return false;
    sink.decisionRows.push_back(
        {kind, sourceRow, targetRow, touches, retiredTouches, candidateTouches});
    return true;
}

bool appendFrontierFromLineage(const std::vector<LineageBlockerCertificate> &blockers,
                               topo::WorkspaceStamp finalStamp, algo::ArtifactStamp artifactStamp,
                               TetFacetSetFrontierSink &frontier) noexcept {
    std::size_t supportRequired = 0;
    for (const auto &blocker : blockers) {
        if (blocker.selectedCells.size() >
            (std::numeric_limits<std::size_t>::max)() - supportRequired)
            return false;
        supportRequired += blocker.selectedCells.size();
    }
    if (blockers.size() > frontier.rows.size() || supportRequired > frontier.supportCells.size())
        return false;
    frontier.reset();
    for (const auto &blocker : blockers) {
        auto &row = frontier.rows[frontier.rowsUsed++];
        row = {};
        row.kind = TetFacetSetWorkKind::Edge;
        row.edgeKey = blocker.request.key;
        row.workspaceStamp = finalStamp;
        row.artifactStamp = artifactStamp;
        row.edgeOutcome =
            edgeFailure(TetEdgeRecoveryStatus::Refused, TetEdgeRecoveryReason::NonCycleLink,
                        TetEdgeRecoveryPhase::Discovery);
        row.supportOffset = frontier.supportCellsUsed;
        row.supportCount = blocker.selectedCells.size();
        for (const auto cell : blocker.selectedCells)
            frontier.supportCells[frontier.supportCellsUsed++] = cell;
        row.dependencyOffset = frontier.dependenciesUsed;
    }
    frontier.complete = true;
    return true;
}

bool exactLineageBlocker(const TetEdgeRecoveryOutcome &outcome) noexcept {
    return outcome.status == TetEdgeRecoveryStatus::Refused &&
           outcome.reason == TetEdgeRecoveryReason::NonCycleLink &&
           outcome.phase == TetEdgeRecoveryPhase::Discovery;
}

bool originalHasExactLineageStop(const TetFacetSetRecoveryOutcome &original) noexcept {
    if (!original.reportComplete || original.processed.empty())
        return false;
    const auto &row = original.processed.back();
    return row.kind == TetFacetSetWorkKind::Edge && exactLineageBlocker(row.edgeOutcome);
}

} // namespace

void setT2CEdgeLineageSinkForTesting(T2CEdgeLineageSinkForTesting *sink) noexcept {
    gT2CEdgeLineageSink = sink;
}

T2CEdgeLineageSinkForTesting *t2cEdgeLineageSinkForTesting() noexcept {
    return gT2CEdgeLineageSink;
}

T2CEdgeLineageRunForTesting
traceT2CEdgeLineageForTesting(const ClosedTri3Artifact &artifact,
                              const std::map<Index, topo::NodeHandle> &localByGlobal,
                              mesh::Classification temporaryRegion, TetFacetSetState &state,
                              topo::TopologyWorkspace &workspace, TetRecoveryContext &context,
                              const TetFacetSetRecoveryOutcome &original,
                              TetFacetSetFrontierSink &replacementFrontier) noexcept {
    T2CEdgeLineageRunForTesting run;
    auto *sink = t2cEdgeLineageSinkForTesting();
    if (!sink || !originalHasExactLineageStop(original))
        return run;
    run.invoked = true;
    try {
        if (!prepareLineageStorage(*sink)) {
            sink->status = T2CEdgeLineageStatusForTesting::MemoryBudgetReached;
            sink->complete = false;
            run.complete = false;
            return run;
        }
        if (sink->attempts.size() >= sink->limits.maxAttempts) {
            sink->status = T2CEdgeLineageStatusForTesting::CapacityReached;
            sink->complete = false;
            run.complete = false;
            return run;
        }
        sink->attempts.emplace_back();
        auto &attempt = sink->attempts.back();
        attempt.initialStamp = topo::WorkspaceStamp::of(workspace);
        attempt.finalStamp = attempt.initialStamp;
        attempt.artifactStamp = artifact.surfaceSetStamp;
        attempt.blockerOffset = sink->blockers.size();
        attempt.recoveryOffset = sink->recoveries.size();
        attempt.passOffset = sink->passRows.size();
        attempt.vertexOffset = sink->vertexRows.size();
        attempt.vertexCount = context.plainVertices.size();
        if (attempt.vertexCount > sink->limits.maxVertexRows - sink->vertexRows.size()) {
            attempt.status = T2CEdgeLineageStatusForTesting::CapacityReached;
            sink->status = attempt.status;
            sink->complete = false;
            run.complete = false;
            return run;
        }
        for (std::size_t i = 0; i < context.plainVertices.size(); ++i)
            sink->vertexRows.push_back({static_cast<DenseTetNodeId>(i), context.denseToNode[i],
                                        context.plainVertices[i].position,
                                        workspace.nodes().get(context.denseToNode[i]).cls});
        attempt.accountedBytes = sink->accountedByteHighWater;

        std::vector<const VolumeArtifactEdge *> edges;
        edges.reserve(artifact.edges.size());
        for (const auto &edge : artifact.edges)
            edges.push_back(&edge);
        std::sort(edges.begin(), edges.end(),
                  [](const auto *a, const auto *b) { return a->key < b->key; });
        std::vector<LineageBlockerCertificate> frozen;
        frozen.reserve((std::min)(edges.size(), sink->limits.maxBlockers));
        std::vector<LineageAssessment> assessments(edges.size());
        std::vector<bool> preparedDeferred(edges.size(), false);
        std::vector<std::size_t> retainedAssessmentRows;
        bool stop = false;

        const auto assessAll = [&]() {
            std::vector<topo::CellHandle> witness(context.cells.size());
            for (std::size_t i = 0; i < edges.size(); ++i) {
                auto &assessment = assessments[i];
                assessment = {};
                if (std::any_of(frozen.begin(), frozen.end(),
                                [&](const auto &b) { return b.request.key == edges[i]->key; })) {
                    assessment.excluded = true;
                    continue;
                }
                if (!mapEdgeRequest(*edges[i], localByGlobal, state, workspace, temporaryRegion,
                                    assessment.request)) {
                    assessment.outcome = edgeFailure(TetEdgeRecoveryStatus::Rejected,
                                                     TetEdgeRecoveryReason::InvalidIdentity,
                                                     TetEdgeRecoveryPhase::Preflight);
                    continue;
                }
                assessment.mapped = true;
                std::size_t used = 0;
                bool complete = true;
                TetRecoveryBindingWorkCounters bindingWork;
                const auto before = context.metrics;
                assessment.outcome =
                    assessTetEdgeIndexedFrontier(assessment.request, state, workspace, context,
                                                 witness, used, complete, &bindingWork);
                const auto after = context.metrics;
                bool workComplete = complete;
                const auto nestedContext = lineageContextDelta(before, after, workComplete);
                const auto nestedBinding = lineageBindingWork(bindingWork, workComplete);
                std::uint64_t nested = nestedContext;
                workComplete = workComplete && lineageAdd(nested, nestedBinding) &&
                               lineageAddWork(*sink, attempt, 1, nested);
                if (!workComplete) {
                    sink->status = T2CEdgeLineageStatusForTesting::InternalError;
                    sink->complete = false;
                    return false;
                }
                assessment.selectedCells.assign(witness.begin(), witness.begin() + used);
            }
            return true;
        };

        while (!stop) {
            if (attempt.passes >= sink->limits.maxPasses) {
                attempt.status = T2CEdgeLineageStatusForTesting::CapacityReached;
                sink->status = attempt.status;
                sink->complete = false;
                break;
            }
            ++attempt.passes;
            if (!assessAll())
                break;
            std::size_t retainedPassRow = 0;
            if (!appendLineagePass(*sink, edges, assessments, frozen, preparedDeferred, context,
                                   retainedAssessmentRows, retainedPassRow)) {
                attempt.status = T2CEdgeLineageStatusForTesting::CapacityReached;
                sink->status = attempt.status;
                sink->complete = false;
                break;
            }
            auto &retainedPass = sink->passRows[retainedPassRow];
            retainedPass.canonicalScanComplete = true;

            std::size_t permanent = edges.size();
            for (std::size_t i = 0; i < edges.size(); ++i) {
                if (assessments[i].excluded)
                    continue;
                if (!assessments[i].mapped) {
                    permanent = i;
                    break;
                }
                const auto status = assessments[i].outcome.status;
                // Marked mutates required/classification metadata.  The
                // testing seam has no exact prepared metadata-write guard,
                // therefore it is a hard canonical-prefix boundary and no
                // later Recovered row may be committed across it.
                if (status != TetEdgeRecoveryStatus::AlreadyPresent &&
                    status != TetEdgeRecoveryStatus::Recovered) {
                    permanent = i;
                    break;
                }
            }
            if (permanent != edges.size()) {
                retainedPass.permanentAssessment = retainedAssessmentRows[permanent];
                sink->assessmentRows[retainedAssessmentRows[permanent]].prefixBoundary = true;
                auto &assessed = assessments[permanent];
                if (!exactLineageBlocker(assessed.outcome)) {
                    attempt.status = T2CEdgeLineageStatusForTesting::PermanentOutcome;
                    sink->status = attempt.status;
                    break;
                }
                LineageBlockerCertificate candidate;
                std::uint64_t supportWork = 0;
                if (!buildLineageCertificate(assessed.request, assessed.selectedCells, context,
                                             workspace, candidate, supportWork) ||
                    !lineageAddWork(*sink, attempt, supportWork, 0)) {
                    attempt.status = T2CEdgeLineageStatusForTesting::InternalError;
                    sink->status = attempt.status;
                    sink->complete = false;
                    break;
                }
                if (!appendLineageBlocker(*sink, candidate, topo::WorkspaceStamp::of(workspace),
                                          assessed.outcome, false)) {
                    attempt.status = T2CEdgeLineageStatusForTesting::CapacityReached;
                    sink->status = attempt.status;
                    sink->complete = false;
                    break;
                }
                retainedPass.discoveredBlocker = candidate.sinkRow;
                bool independent = true;
                for (const auto &blocker : frozen) {
                    if (sink->blockerPairTests >= sink->limits.maxPairTests ||
                        !lineageAdd(sink->blockerPairTests, 1)) {
                        independent = false;
                        sink->complete = false;
                        attempt.status = T2CEdgeLineageStatusForTesting::WorkBudgetReached;
                        break;
                    }
                    std::uint64_t pairWork = 0;
                    const bool pairIndependent =
                        lineageBlockersIndependent(blocker, candidate, context, pairWork);
                    independent = independent && pairIndependent;
                    if (!appendLineageDecision(
                            *sink, T2CEdgeLineageDependencyKindForTesting::BlockerToBlocker,
                            blocker.sinkRow, candidate.sinkRow, !pairIndependent)) {
                        independent = false;
                        sink->complete = false;
                        attempt.status = T2CEdgeLineageStatusForTesting::CapacityReached;
                        break;
                    }
                    if (!lineageAddWork(*sink, attempt, pairWork, 0)) {
                        independent = false;
                        break;
                    }
                }
                retainedPass.decisionCount =
                    sink->decisionRows.size() - retainedPass.decisionOffset;
                if (!sink->complete) {
                    auto &retained = sink->blockers[candidate.sinkRow];
                    retained.pairwiseIndependent = false;
                    retained.deferred = true;
                    break;
                }
                if (!independent) {
                    auto &retained = sink->blockers[candidate.sinkRow];
                    retained.pairwiseIndependent = false;
                    retained.deferred = true;
                    sink->allBlockersPairwiseIndependent = false;
                    ++attempt.deferredBlockers;
                    attempt.status = T2CEdgeLineageStatusForTesting::BlockerOverlap;
                    sink->status = attempt.status;
                    break;
                }
                frozen.push_back(std::move(candidate));
                continue;
            }

            std::vector<std::vector<std::size_t>> dependencies(edges.size());
            std::vector<bool> candidate(edges.size(), false);
            for (std::size_t i = 0; i < edges.size(); ++i) {
                const auto status = assessments[i].outcome.status;
                candidate[i] = status == TetEdgeRecoveryStatus::Recovered;
                sink->assessmentRows[retainedAssessmentRows[i]].candidate = candidate[i];
                sink->assessmentRows[retainedAssessmentRows[i]].deferred = preparedDeferred[i];
                if (!candidate[i] || preparedDeferred[i])
                    continue;
                for (const auto &blocker : frozen) {
                    if (sink->candidatePairTests >= sink->limits.maxPairTests ||
                        !lineageAdd(sink->candidatePairTests, 1) ||
                        !lineageAddWork(*sink, attempt, 1, 0)) {
                        attempt.status = T2CEdgeLineageStatusForTesting::WorkBudgetReached;
                        sink->status = attempt.status;
                        sink->complete = false;
                        stop = true;
                        break;
                    }
                    const bool touches = lineageSelectedTouchesRequest(assessments[i].selectedCells,
                                                                       blocker.request, context);
                    if (!appendLineageDecision(
                            *sink, T2CEdgeLineageDependencyKindForTesting::CandidateToFrozenBlocker,
                            retainedAssessmentRows[i], blocker.sinkRow, touches)) {
                        attempt.status = T2CEdgeLineageStatusForTesting::CapacityReached;
                        sink->status = attempt.status;
                        sink->complete = false;
                        stop = true;
                        break;
                    }
                    if (touches) {
                        if (!preparedDeferred[i])
                            ++attempt.deferredRecoveries;
                        preparedDeferred[i] = true;
                    }
                }
                sink->assessmentRows[retainedAssessmentRows[i]].deferred = preparedDeferred[i];
                if (stop)
                    break;
            }
            if (stop)
                break;
            for (std::size_t i = 0; i < edges.size(); ++i) {
                if (!candidate[i] || preparedDeferred[i])
                    continue;
                for (std::size_t j = 0; j < edges.size(); ++j) {
                    if (i == j || !candidate[j] || preparedDeferred[j])
                        continue;
                    if (sink->candidatePairTests >= sink->limits.maxPairTests ||
                        !lineageAdd(sink->candidatePairTests, 1) ||
                        !lineageAddWork(*sink, attempt, 1, 0)) {
                        attempt.status = T2CEdgeLineageStatusForTesting::WorkBudgetReached;
                        sink->status = attempt.status;
                        sink->complete = false;
                        stop = true;
                        break;
                    }
                    const bool touches = lineageSelectedTouchesRequest(
                        assessments[i].selectedCells, assessments[j].request, context);
                    if (!appendLineageDecision(
                            *sink, T2CEdgeLineageDependencyKindForTesting::CandidateToCandidate,
                            retainedAssessmentRows[i], retainedAssessmentRows[j], touches)) {
                        attempt.status = T2CEdgeLineageStatusForTesting::CapacityReached;
                        sink->status = attempt.status;
                        sink->complete = false;
                        stop = true;
                        break;
                    }
                    if (touches)
                        dependencies[i].push_back(j);
                }
                if (stop)
                    break;
            }
            if (stop)
                break;
            retainedPass.decisionCount = sink->decisionRows.size() - retainedPass.decisionOffset;
            std::size_t selected = edges.size();
            for (std::size_t i = 0; i < edges.size(); ++i)
                if (candidate[i] && !preparedDeferred[i] && dependencies[i].empty()) {
                    selected = i;
                    break;
                }
            if (selected == edges.size()) {
                attempt.edgeFixedPoint = true;
                retainedPass.edgeFixedPoint = true;
                attempt.status = frozen.empty() ? T2CEdgeLineageStatusForTesting::NoProgress
                                                : T2CEdgeLineageStatusForTesting::EdgeFixedPoint;
                sink->status = attempt.status;
                break;
            }
            retainedPass.selectedAssessment = retainedAssessmentRows[selected];

            if (sink->recoveries.size() >= sink->limits.maxScratchRecoveries) {
                attempt.status = T2CEdgeLineageStatusForTesting::CapacityReached;
                sink->status = attempt.status;
                sink->complete = false;
                break;
            }
            sink->recoveries.emplace_back();
            auto &recoveryRow = sink->recoveries.back();
            retainedPass.preparedRecovery = sink->recoveries.size() - 1;
            recoveryRow.key = assessments[selected].request.key;
            recoveryRow.beforeStamp = topo::WorkspaceStamp::of(workspace);
            std::vector<std::pair<topo::EdgeHandle, VolumeArtifactEdgeKey>> boundEdges;
            boundEdges.reserve(state.edgeBindings.size());
            for (const auto &[key, binding] : state.edgeBindings)
                boundEdges.push_back({binding.edge, key});
            recoveryRow.boundEdgeOffset = sink->boundEdgeRows.size();
            recoveryRow.boundEdgeCount = boundEdges.size();
            if (boundEdges.size() > sink->limits.maxBoundEdgeRows - sink->boundEdgeRows.size()) {
                attempt.status = T2CEdgeLineageStatusForTesting::CapacityReached;
                sink->status = attempt.status;
                sink->complete = false;
                break;
            }
            for (const auto &[handle, key] : boundEdges)
                sink->boundEdgeRows.push_back({handle, key});
            TetRecoveryBindingWorkCounters recoveryBinding;
            const auto beforeMetrics = context.metrics;
            TetEdgeIndexedOutcome recovered;
            LineagePreparedGuardState guard{sink,        &attempt,     &frozen,
                                            &boundEdges, &recoveryRow, sink->recoveries.size() - 1};
            if (assessments[selected].outcome.status == TetEdgeRecoveryStatus::Recovered) {
                recovered = recoverTetEdgeIndexedPreparedGuardedForTesting(
                    assessments[selected].request, state, workspace, context, lineagePreparedGuard,
                    &guard, &recoveryBinding);
            } else {
                recovered = recoverTetEdgeIndexedDetailed(assessments[selected].request, state,
                                                          workspace, context, &recoveryBinding);
                recoveryRow.allowed = recovered.recovery.succeeded();
                recoveryRow.samePreparedTransaction = false;
            }
            const auto afterMetrics = context.metrics;
            retainedPass.decisionCount = sink->decisionRows.size() - retainedPass.decisionOffset;
            bool nestedComplete = true;
            auto nested = lineageContextDelta(beforeMetrics, afterMetrics, nestedComplete);
            const auto bindingUnits = lineageBindingWork(recoveryBinding, nestedComplete);
            nestedComplete = nestedComplete && lineageAdd(nested, bindingUnits) &&
                             lineageAdd(nested, 1) && lineageAddWork(*sink, attempt, 0, nested);
            recoveryRow.outcome = recovered.recovery;
            if (guard.capacityFailed) {
                attempt.status = T2CEdgeLineageStatusForTesting::CapacityReached;
                sink->status = attempt.status;
                sink->complete = false;
                break;
            }
            if (guard.workFailed || !nestedComplete) {
                attempt.status = T2CEdgeLineageStatusForTesting::WorkBudgetReached;
                sink->status = attempt.status;
                sink->complete = false;
                break;
            }
            if (assessments[selected].outcome.status == TetEdgeRecoveryStatus::Recovered &&
                !guard.invoked) {
                sink->allPreparedTransactionsGuarded = false;
                attempt.status = T2CEdgeLineageStatusForTesting::InternalError;
                sink->status = attempt.status;
                sink->complete = false;
                break;
            }
            if (!recoveryRow.allowed) {
                preparedDeferred[selected] = true;
                sink->assessmentRows[retainedAssessmentRows[selected]].deferred = true;
                sink->assessmentRows[retainedAssessmentRows[selected]].preparedDeferred = true;
                ++attempt.deferredRecoveries;
                continue;
            }
            if (!recovered.recovery.succeeded() ||
                !appendCreatedRows(*sink, recoveryRow, recovered.createdCells)) {
                attempt.status =
                    recovered.recovery.reason == TetEdgeRecoveryReason::AllocationFailed
                        ? T2CEdgeLineageStatusForTesting::AllocationFailed
                        : T2CEdgeLineageStatusForTesting::PermanentOutcome;
                sink->status = attempt.status;
                sink->complete = false;
                break;
            }
            recoveryRow.afterStamp = topo::WorkspaceStamp::of(workspace);
            preparedDeferred.assign(edges.size(), false);

            // Currentness is the first post-commit read.  A stale support
            // handle may not be laundered by a later assessment/rebuild.
            for (const auto &blocker : frozen)
                if (!lineageCertificateCurrent(blocker, workspace, context)) {
                    attempt.status = T2CEdgeLineageStatusForTesting::StaleSupport;
                    sink->status = attempt.status;
                    sink->complete = false;
                    stop = true;
                    break;
                }
            if (stop)
                break;
            for (auto &blocker : frozen) {
                auto request = blocker.request;
                request.stamp = topo::WorkspaceStamp::of(workspace);
                std::vector<topo::CellHandle> witness(context.cells.size());
                std::size_t used = 0;
                bool complete = true;
                TetRecoveryBindingWorkCounters bindingWork;
                const auto before = context.metrics;
                const auto outcome = assessTetEdgeIndexedFrontier(
                    request, state, workspace, context, witness, used, complete, &bindingWork);
                const auto after = context.metrics;
                bool nestedWorkComplete = complete;
                auto revalidationNested = lineageContextDelta(before, after, nestedWorkComplete);
                const auto revalidationBinding =
                    lineageBindingWork(bindingWork, nestedWorkComplete);
                nestedWorkComplete = nestedWorkComplete &&
                                     lineageAdd(revalidationNested, revalidationBinding) &&
                                     lineageAddWork(*sink, attempt, 1, revalidationNested);
                LineageBlockerCertificate refreshed;
                std::uint64_t supportWork = 0;
                if (!nestedWorkComplete || !exactLineageBlocker(outcome) ||
                    !buildLineageCertificate(
                        request, std::span<const topo::CellHandle>(witness.data(), used), context,
                        workspace, refreshed, supportWork) ||
                    !lineageAddWork(*sink, attempt, supportWork, 0) ||
                    !lineageCertificateEqual(blocker, refreshed)) {
                    attempt.status = T2CEdgeLineageStatusForTesting::StaleSupport;
                    sink->status = attempt.status;
                    sink->complete = false;
                    stop = true;
                    break;
                }
                blocker.request.stamp = request.stamp;
                if (blocker.sinkRow >= sink->blockers.size()) {
                    stop = true;
                    sink->complete = false;
                    break;
                }
                auto &retained = sink->blockers[blocker.sinkRow];
                retained.finalValidatedStamp = request.stamp;
                if (!lineageAdd(retained.revalidations, 1)) {
                    stop = true;
                    sink->complete = false;
                    break;
                }
            }
        }

        attempt.finalStamp = topo::WorkspaceStamp::of(workspace);
        attempt.blockerCount = sink->blockers.size() - attempt.blockerOffset;
        attempt.recoveryCount = sink->recoveries.size() - attempt.recoveryOffset;
        attempt.passCount = sink->passRows.size() - attempt.passOffset;
        // A partial canonical scan (including blocker overlap) is evidence,
        // never a replacement frontier.  Only a fully completed edge-stage
        // fixed point may be exported to the existing facet-SCC path.
        const bool usable = !frozen.empty() && sink->complete && attempt.edgeFixedPoint &&
                            attempt.status == T2CEdgeLineageStatusForTesting::EdgeFixedPoint;
        if (usable) {
            if (!appendFrontierFromLineage(frozen, attempt.finalStamp, artifact.surfaceSetStamp,
                                           replacementFrontier)) {
                attempt.status = T2CEdgeLineageStatusForTesting::CapacityReached;
                sink->status = attempt.status;
                sink->complete = false;
            } else {
                attempt.frontierReplaced = true;
                run.frontierReplaced = true;
            }
        }
        attempt.complete = sink->complete;
        run.complete = sink->complete;
        run.workAccountingComplete = sink->workAccountingComplete;
        run.memoryAccountingComplete = sink->memoryAccountingComplete;
        run.coordinatorWork = attempt.coordinatorWork;
        run.nestedRecoveryWork = attempt.nestedRecoveryWork;
        return run;
    } catch (const std::bad_alloc &) {
        sink->status = T2CEdgeLineageStatusForTesting::AllocationFailed;
        sink->complete = false;
        run.complete = false;
        return run;
    } catch (...) {
        sink->status = T2CEdgeLineageStatusForTesting::InternalError;
        sink->complete = false;
        run.complete = false;
        return run;
    }
}
#endif

} // namespace paacs::bnd
