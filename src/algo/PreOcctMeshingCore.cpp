// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/PreOcctMeshingCore.hpp"
#include "paacs/algo/Smoothing.hpp"

#include <algorithm>
#include <limits>
#include <new>

namespace paacs::algo {
namespace {

PreOcctMeshingStatus mapStatus(BoundaryConformingV3Status value) noexcept {
    switch (value) {
        case BoundaryConformingV3Status::Published:
            return PreOcctMeshingStatus::Published;
        case BoundaryConformingV3Status::Rejected:
            return PreOcctMeshingStatus::Rejected;
        case BoundaryConformingV3Status::Unsupported:
        case BoundaryConformingV3Status::Refused:
            return PreOcctMeshingStatus::Unsupported;
        case BoundaryConformingV3Status::Failed:
            return PreOcctMeshingStatus::Failed;
    }
    return PreOcctMeshingStatus::Failed;
}

PreOcctMeshingStatus mapStatus(T3AStatus value) noexcept {
    switch (value) {
        case T3AStatus::Published:
            return PreOcctMeshingStatus::Published;
        case T3AStatus::Incomplete:
            return PreOcctMeshingStatus::Incomplete;
        case T3AStatus::Rejected:
            return PreOcctMeshingStatus::Rejected;
        case T3AStatus::Unsupported:
            return PreOcctMeshingStatus::Unsupported;
        case T3AStatus::Failed:
            return PreOcctMeshingStatus::Failed;
    }
    return PreOcctMeshingStatus::Failed;
}

PreOcctMeshingStatus mapStatus(T4AStatus value) noexcept {
    switch (value) {
        case T4AStatus::Published:
            return PreOcctMeshingStatus::Published;
        case T4AStatus::Incomplete:
            return PreOcctMeshingStatus::Incomplete;
        case T4AStatus::Rejected:
            return PreOcctMeshingStatus::Rejected;
        case T4AStatus::Unsupported:
            return PreOcctMeshingStatus::Unsupported;
        case T4AStatus::Failed:
            return PreOcctMeshingStatus::Failed;
    }
    return PreOcctMeshingStatus::Failed;
}

PreOcctMeshingStatus mapStatus(T4BStatus value) noexcept {
    switch (value) {
        case T4BStatus::Published:
            return PreOcctMeshingStatus::Published;
        case T4BStatus::Incomplete:
            return PreOcctMeshingStatus::Incomplete;
        case T4BStatus::Rejected:
            return PreOcctMeshingStatus::Rejected;
        case T4BStatus::Unsupported:
            return PreOcctMeshingStatus::Unsupported;
        case T4BStatus::Failed:
            return PreOcctMeshingStatus::Failed;
    }
    return PreOcctMeshingStatus::Failed;
}

bool auditAccepted(const mesh::AuditResult& audit,
                   bool exactBoundary) noexcept {
    const bool base = audit.passed() && audit.nInverted == 0 &&
                      audit.nInvalidTet10 == 0 && audit.nNonFinite == 0 &&
                      audit.malformedCells == 0 &&
                      audit.deferredAuditCells == 0 &&
                      audit.boundaryDuplicates == 0 &&
                      audit.nonManifoldFaces == 0;
    if (!base || !exactBoundary) return base;
    return audit.boundaryOrphans == 0 && audit.boundaryUnexpected == 0 &&
           audit.boundaryPrescribed == audit.boundaryMatched;
}

}  // namespace

PreOcctTetRequest makePreOcctTetRequest(
        PreOcctTetOutput output,
        VolumeBoundaryProviderStamp ownershipStamp,
        Tet10PlacementStamp placementStamp,
        const PreOcctResourceLimits& limits) noexcept {
    PreOcctTetRequest request;
    request.output = output;
    request.interiorSmoothingIterations =
        output == PreOcctTetOutput::QualityTet4 ? 2 : 0;
    request.pipeline.source.source.conformingRequest.expectedSurfaceSetStamp =
        {0x5041414353563031ULL, 1};
    request.pipeline.source.source.conformingRequest
        .expectedOwnershipProviderStamp = ownershipStamp;
    request.pipeline.source.source.conformingRequest.envelopeMode =
        RecoveryEnvelopeMode::DerivedFiniteBox;

    auto& refinement = request.pipeline.source.source.refinement;
    refinement.maxEdgeToTarget = 1.5;
    refinement.maxRadiusEdge = 3.0;
    refinement.maxInsertions = std::min<std::uint64_t>(64, limits.maximumNodes);
    refinement.maxWorkUnits = std::min<std::uint64_t>(
        limits.maximumWorkUnits, 500000);
    refinement.maxOutputNodes = limits.maximumNodes;
    refinement.maxOutputTets = limits.maximumTets;
    refinement.maxReportRows = std::max<std::uint64_t>(
        325, std::min<std::uint64_t>(limits.maximumReportRows, 65536));
    refinement.maxReportWitnesses = std::max<std::uint64_t>(
        limits.maximumReportRows, 500000);
    refinement.maxReportTopologyHandles =
        refinement.maxWorkUnits + refinement.maxInsertions;
    refinement.maxCandidateAttemptsPerViolation = 5;

    auto& quality = request.pipeline.source.quality;
    quality.mode = T4AQualityMode::ReportOnly;
    quality.minimumRadiusRatio = 0.15;
    quality.maxResultNodes = limits.maximumNodes;
    quality.maxResultTets = limits.maximumTets;
    quality.maxReportRows = limits.maximumReportRows;
    quality.maxWorstTailRows = std::min<std::uint64_t>(
        limits.maximumReportRows, 4096);

    const auto saturatingMul = [](std::uint64_t a, std::uint64_t b) {
        const auto max = (std::numeric_limits<std::uint64_t>::max)();
        return a != 0 && b > max / a ? max : a * b;
    };
    request.pipeline.expectedPlacementStamp = placementStamp;
    auto& tet10 = request.pipeline.policy;
    tet10.minimumCertifiedJacobianRatio = 1.0e-4;
    tet10.maxCertificateDepth = 6;
    tet10.maxCertificateCells = saturatingMul(limits.maximumTets, 16U);
    tet10.maxCoefficientTests = saturatingMul(limits.maximumWorkUnits, 8U);
    tet10.maxWorkUnits = saturatingMul(limits.maximumWorkUnits, 16U);
    tet10.maxResultNodes = saturatingMul(limits.maximumNodes, 4U);
    tet10.maxResultCells = limits.maximumTets;
    tet10.maxResultConnectivity = saturatingMul(limits.maximumTets, 10U);
    tet10.maxPlacementRows = saturatingMul(limits.maximumNodes, 4U);
    tet10.maxCertificateRows = limits.maximumTets;
    tet10.maxUnattemptedRows = limits.maximumTets;
    return request;
}

void enableBoundedTetQuality(
        PreOcctTetRequest& request,
        double minimumRadiusRatio,
        std::uint64_t maximumOperations) noexcept {
    auto& quality = request.pipeline.source.quality;
    quality.mode = T4AQualityMode::BoundedImprove;
    quality.minimumRadiusRatio = minimumRadiusRatio;
    quality.minimumImprovement = 1.0e-8;
    quality.maxOperations = maximumOperations;
    quality.maxSelectedGenerations =
        std::max<std::uint64_t>(64, maximumOperations * 64U);
    quality.maxCandidateAttemptsPerCell = 1;
    quality.maxCreatedTetsPerOperation = 128;
    quality.maxWorkUnits = std::max<std::uint64_t>(
        1000000U, maximumOperations * 250000U);
    quality.maxNewNodes = maximumOperations;
    quality.maxDependencyRowsPerAttempt = 1000;
    quality.maxProtectedWitnessesPerAttempt = 1000;
    quality.maxTopologyRowsPerOperation = 10000;

    // The accepted T4-A report shell validates that its configured row budget
    // can represent every bounded outcome.  Compute a conservative cap from
    // the public limits instead of leaving callers to reverse engineer that
    // internal accounting requirement.
    const auto saturatingAdd = [](std::uint64_t a, std::uint64_t b) {
        const auto max = (std::numeric_limits<std::uint64_t>::max)();
        return b > max - a ? max : a + b;
    };
    const auto saturatingMul = [](std::uint64_t a, std::uint64_t b) {
        const auto max = (std::numeric_limits<std::uint64_t>::max)();
        return a != 0 && b > max / a ? max : a * b;
    };
    std::uint64_t required = quality.maxResultTets;
    required = saturatingAdd(required,
        saturatingMul(quality.maxSelectedGenerations,
                      2U + quality.maxDependencyRowsPerAttempt +
                      quality.maxProtectedWitnessesPerAttempt));
    required = saturatingAdd(required,
        saturatingMul(maximumOperations,
                      1U + quality.maxTopologyRowsPerOperation +
                      quality.maxCreatedTetsPerOperation));
    required = saturatingAdd(required,
        saturatingMul(quality.maxWorstTailRows, 2U));
    quality.maxReportRows = std::max(quality.maxReportRows, required);
}

PreOcctMeshingOutcome meshPreparedTetModel(
        const geo::Model& model,
        const sizing::MeshSpec& spec,
        const PreOcctTetRequest& request,
        mesh::Mesh& destination,
        const VolumeBoundaryOwnershipProvider& ownershipProvider,
        const Tet10PlacementProvider* placementProvider) noexcept {
    PreOcctMeshingOutcome out;
    out.output = request.output;
    mesh::Mesh staged;
    try {
        switch (request.output) {
            case PreOcctTetOutput::BoundaryConformingTet4: {
                out.phase = PreOcctMeshingPhase::BoundaryConforming;
                auto result = std::make_unique<BoundaryConformingV3Outcome>(
                    meshModelBoundaryConforming(
                        model, spec,
                        request.pipeline.source.source.conformingRequest,
                        staged, ownershipProvider));
                out.status = mapStatus(result->status);
                if (!result->published()) {
                    out.reason =
                        PreOcctMeshingReason::BoundaryConformingStopped;
                    out.conforming = std::move(result);
                    return out;
                }
                out.conforming = std::move(result);
                break;
            }
            case PreOcctTetOutput::RefinedTet4: {
                out.phase = PreOcctMeshingPhase::Refinement;
                auto result = std::make_unique<T3AOutcome>(
                    meshModelRefinedTet4(
                        model, spec, request.pipeline.source.source,
                        staged, ownershipProvider));
                out.status = mapStatus(result->status);
                if (!result->published()) {
                    out.reason = PreOcctMeshingReason::RefinementStopped;
                    out.refined = std::move(result);
                    return out;
                }
                out.refined = std::move(result);
                break;
            }
            case PreOcctTetOutput::QualityTet4: {
                out.phase = PreOcctMeshingPhase::Quality;
                auto result = std::make_unique<T4AOutcome>(
                    meshModelQualityTet4(
                        model, spec, request.pipeline.source,
                        staged, ownershipProvider));
                out.status = mapStatus(result->status);
                if (!result->published()) {
                    out.reason = PreOcctMeshingReason::QualityStopped;
                    out.quality = std::move(result);
                    return out;
                }
                out.quality = std::move(result);
                break;
            }
            case PreOcctTetOutput::CertifiedTet10: {
                if (!placementProvider) {
                    out.status = PreOcctMeshingStatus::Rejected;
                    out.reason = PreOcctMeshingReason::
                        MissingTet10PlacementProvider;
                    out.phase = PreOcctMeshingPhase::Preflight;
                    return out;
                }
                out.phase = PreOcctMeshingPhase::Tet10;
                auto result = std::make_unique<T4BOutcome>(
                    meshModelCurvedTet10(
                        model, spec, request.pipeline, staged,
                        ownershipProvider, *placementProvider));
                out.status = mapStatus(result->status);
                if (!result->published()) {
                    out.reason = PreOcctMeshingReason::Tet10Stopped;
                    out.tet10 = std::move(result);
                    return out;
                }
                out.tet10 = std::move(result);
                break;
            }
        }

        out.tetQualityBeforePost = mesh::analyzeTets(staged);
        if (request.interiorSmoothingIterations > 0 &&
            request.output != PreOcctTetOutput::CertifiedTet10) {
            SmoothOptions smoothing;
            smoothing.maxIterations = request.interiorSmoothingIterations;
            (void)smoothInterior(staged, smoothing);
        }

        out.phase = PreOcctMeshingPhase::FinalAudit;
        mesh::AuditOptions auditOptions;
        auditOptions.boundaryMode = request.requireExactBoundaryAudit
            ? mesh::BoundaryComparisonMode::Required
            : mesh::BoundaryComparisonMode::WhenSurfacePresent;
        auditOptions.boundaryOrphanBudget = 0;
        auditOptions.boundaryUnexpectedBudget = 0;
        auditOptions.duplicateFacetsAreErrors = true;
        auditOptions.unusedNodesAreErrors = true;
        // Solid ownership is not mesh connectivity: two regions sharing a
        // conformal face form one face-connected component. Report the actual
        // count; authoritative domain extraction checks each region's topology.
        out.audit = mesh::auditMesh(staged, auditOptions);
        if (!auditAccepted(out.audit, request.requireExactBoundaryAudit)) {
            out.status = PreOcctMeshingStatus::Failed;
            out.reason = PreOcctMeshingReason::FinalAuditFailed;
            return out;
        }
        out.tetQuality = mesh::analyzeTets(staged);
        out.mixedQuality = analyzeMixedMeshQuality(staged);
        if (!out.mixedQuality.complete) {
            out.status = PreOcctMeshingStatus::Failed;
            out.reason = out.mixedQuality.allocationFailed
                ? PreOcctMeshingReason::AllocationFailed
                : PreOcctMeshingReason::InternalError;
            return out;
        }
        out.strictHash = mesh::meshHash(staged);
        out.semanticHash = mesh::semanticHashEx(staged);
        out.publishedNodes = static_cast<std::uint64_t>(staged.nodeCount());
        out.publishedCells = static_cast<std::uint64_t>(staged.cellCount());

        out.phase = PreOcctMeshingPhase::Publication;
        destination.swap(staged);
        out.status = PreOcctMeshingStatus::Published;
        out.reason = PreOcctMeshingReason::None;
        out.destinationUnchanged = false;
        return out;
    } catch (const std::bad_alloc&) {
        out.status = PreOcctMeshingStatus::Failed;
        out.reason = PreOcctMeshingReason::AllocationFailed;
    } catch (...) {
        out.status = PreOcctMeshingStatus::Failed;
        out.reason = PreOcctMeshingReason::InternalError;
    }
    return out;
}

}  // namespace paacs::algo
