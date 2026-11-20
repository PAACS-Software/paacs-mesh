// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sdk/MeshingSession.hpp"
#include "paacs/mesh/MeshAudit.hpp"

namespace paacs::sdk {
void MeshingSession::report(SessionStage stage, std::uint64_t done, std::uint64_t total) const {
    if (options_.progress)
        options_.progress(stage, done, total);
}
SessionResult MeshingSession::executeHybrid(const algo::HybridExecutionRequest &request,
                                            algo::HybridMesh &destination) const noexcept {
    SessionResult out;
    out.manifest.requestedThreads = options_.threads;
    if (cancelled()) {
        out.cancelled = true;
        out.diagnostic = "cancelled before preflight";
        return out;
    }
    report(SessionStage::Preflight, 0, 1);
    algo::HybridExecutionRequest local = request;
    local.progress = [&](algo::HybridExecutionStage stage, std::uint64_t done,
                         std::uint64_t total) {
        SessionStage mapped = SessionStage::Generation;
        if (stage == algo::HybridExecutionStage::Preflight)
            mapped = SessionStage::Preflight;
        else if (stage == algo::HybridExecutionStage::Audit)
            mapped = SessionStage::Audit;
        else if (stage == algo::HybridExecutionStage::Publication)
            mapped = SessionStage::Publication;
        report(mapped, done, total);
        return !cancelled();
    };
    const auto result = algo::executeHybridMeshingPlan(local, destination);
    out.completed = result.published;
    out.cancelled = result.reason == algo::HybridExecutionReason::Cancelled;
    out.diagnostic = out.completed   ? "complete"
                     : out.cancelled ? "cancelled"
                                     : "hybrid meshing failed";
    out.manifest.planHash = request.plan ? request.plan->deterministicHash() : 0;
    out.manifest.meshHash = mesh::meshHash(destination.fixed);
    return out;
}
SessionResult MeshingSession::generateMappedHex(algo::MappedHexBlockRequest request,
                                                mesh::Mesh &destination) const noexcept {
    SessionResult out;
    out.manifest.requestedThreads = options_.threads;
    if (cancelled()) {
        out.cancelled = true;
        out.diagnostic = "cancelled before generation";
        return out;
    }
    report(SessionStage::Generation, 0, 1);
    if (options_.threads > 0)
        request.executionThreads = options_.threads;
    const auto result = algo::generateMappedHexBlock(request, destination);
    out.completed = result.published();
    out.diagnostic = out.completed ? "complete" : "mapped Hex8 generation failed";
    if (out.completed)
        out.manifest.meshHash = mesh::meshHash(destination);
    report(SessionStage::Publication, out.completed ? 1 : 0, 1);
    return out;
}

SessionResult MeshingSession::generateMappedHex(algo::MappedHexBlockRequest request,
                                                algo::MappedHexWorkspace &workspace,
                                                mesh::Mesh &destination) const noexcept {
    SessionResult out;
    out.manifest.requestedThreads = options_.threads;
    if (cancelled()) {
        out.cancelled = true;
        out.diagnostic = "cancelled before generation";
        return out;
    }
    report(SessionStage::Generation, 0, 1);
    if (options_.threads > 0)
        request.executionThreads = options_.threads;
    const auto result = algo::generateMappedHexBlock(request, workspace, destination);
    out.completed = result.published();
    out.diagnostic = out.completed ? "complete" : "mapped Hex8 generation failed";
    if (out.completed)
        out.manifest.meshHash = mesh::meshHash(destination);
    report(SessionStage::Publication, out.completed ? 1 : 0, 1);
    return out;
}

SessionResult
MeshingSession::meshPreparedTet(const geo::Model &model, const sizing::MeshSpec &spec,
                                const algo::PreOcctTetRequest &request, mesh::Mesh &destination,
                                const algo::VolumeBoundaryOwnershipProvider &ownershipProvider,
                                const algo::Tet10PlacementProvider *placementProvider,
                                algo::PreOcctMeshingOutcome *detailedOutcome) const noexcept {
    SessionResult out;
    out.manifest.requestedThreads = options_.threads;
    if (cancelled()) {
        out.cancelled = true;
        out.diagnostic = "cancelled before Tet preflight";
        return out;
    }
    report(SessionStage::Preflight, 0, 1);
    sizing::MeshSpec localSpec = spec;
    if (options_.threads > 0)
        localSpec.threads = options_.threads;
    report(SessionStage::Generation, 0, 1);
    auto result = algo::meshPreparedTetModel(model, localSpec, request, destination,
                                             ownershipProvider, placementProvider);
    out.completed = result.published();
    out.diagnostic = out.completed ? "complete"
                     : result.status == algo::PreOcctMeshingStatus::Incomplete
                         ? "Tet pipeline incomplete within configured budgets"
                         : "Tet production meshing failed";
    if (out.completed) {
        out.manifest.meshHash = result.strictHash;
        report(SessionStage::Audit, 1, 1);
        report(SessionStage::Publication, 1, 1);
    }
    if (detailedOutcome)
        *detailedOutcome = std::move(result);
    return out;
}
} // namespace paacs::sdk
