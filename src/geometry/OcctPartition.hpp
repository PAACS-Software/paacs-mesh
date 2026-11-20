// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private optional-target CAD partitioning. Never installed as a core API.
#include "paacs/geometry/OcctModel.hpp"
#include "paacs/mesh/Mesh.hpp"
#include <array>
#include <atomic>

namespace paacs::geo::detail {
enum class OcctPartitionPolicy : std::uint8_t {
    CartesianMidplanes, CoaxialAngular, CartesianGrid, NativeFeatureAxes,
    // Native cylinder/cone axes plus exact circular-edge axes and interior
    // circle support planes. No fitted axes or independent torus-face axes.
    NativeCircularFeatureAxes,
    // Largest-area native cylinder axis, concave cylinder/cone axes and
    // planar supports normal to the principal axis. Existing policies retain
    // their plane order and options.
    NativeConcaveFeatureAxes,
    // Largest-area native cylinder axis and all interior planar supports.
    // Unlike the concave policy, no additional displaced curved-face axes
    // are proposed. Each bounded cut retains original native history.
    NativePrincipalFeatureAxes
};
enum class OcctPartitionReason : std::uint8_t {
    None, InvalidRequest, UnsupportedPolicy, BudgetExceeded, Cancelled,
    AllocationFailed, KernelFailure, InvalidResult, ProvenanceLost,
    InterfaceOwnershipInvalid, VolumeMismatch, SourceMutation
};
struct OcctPartitionPlane { Vec3 origin, normal; };
enum class OcctTrimToleranceReason : std::uint8_t {
    None, InvalidRequest, SupportMismatch, SourceCoverageLost, OriginalToleranceExceeded,
    DriftEnvelopeExceeded, NativeQueryFailed, BudgetExceeded, Cancelled, AllocationFailed,
    SourceMutation
};
// Sampled trimming witness, not a continuous global geometry certificate.
// Complete original CAD/mesh validation remains the transactional caller's gate.
struct OcctTrimToleranceReceipt {
    bool attempted = false, completed = false, admitted = false;
    OcctTrimToleranceReason reason = OcctTrimToleranceReason::InvalidRequest;
    bool supportsUnchanged = false, originalFacesRetained = false;
    bool originalEdgesRetained = false, originalVerticesRetained = false;
    bool sourceUnmodified = false;
    std::uint64_t samples = 0, work = 0, exteriorFaces = 0, sharedCutFaces = 0;
    std::uint64_t outsideOriginalTolerance = 0, sourceFaces = 0, sourceEdges = 0;
    std::uint64_t sourceVertices = 0, degeneratedSourceEdges = 0;
    double maximumPhysicalDrift = 0, maximumResidual = 0, maximumToleranceRatio = 0;
    double maximumOriginalVertexShift = 0, sourceToleranceVolumeEnvelope = 0;
    double nativeVolumeDrift = 0;
    int firstFailedParentDimension = -1, firstFailedParentTag = -1;
};
struct OcctPartitionOptions {
    OcctPartitionPolicy policy = OcctPartitionPolicy::CartesianMidplanes;
    // Empty selects every original body. Otherwise the result contains only
    // the selected bodies' pieces; caller-owned successful bodies stay separate.
    std::vector<int> sourceBodyTags;
    unsigned cartesianDivisions = 2; // CartesianGrid only, bounded 2..5.
    unsigned angularPlanes = 8; // N axial planes form up to 2*N sectors.
    unsigned featureAxisPlanes = 2; // 2..8; concave/principal policies 2..16.
    // NativeFeatureAxes uses two normal planes containing each native
    // cylindrical/conical axis. Planar support cuts open concave flat features.
    bool includeCartesianMidplanes = true;
    bool includePlanarSupportPlanes = true;
    std::uint64_t maximumPlanes = 48, maximumPieces = 512;
    std::uint64_t maximumHistoryWork = 2000000;
    std::uint64_t maximumSnapshotAuditBytes = 64ULL * 1024 * 1024;
    NativeCadImportOptions snapshotLimits{};
    double integrationRelativeError = 1e-10;
    double maximumRelativeVolumeChange = 1e-6;
    bool verifyNativeCommon = false;
    const std::atomic_bool* cancellation = nullptr;
    // Explicit additional native support cuts, applied after the policy.
    // Every attempt re-splits the original owned copy and retains history.
    // Empty preserves the existing policies and their work counts.
    std::vector<OcctPartitionPlane> extraPlanes;
    // Explicit new-policy retry only. Never modifies native source tolerance or
    // geometry. The caller must cap this finite drift by its CAD chord bound.
    bool allowSampledTrimTolerance = false;
    // A bounded split may retry native volume integration on B-spline spans
    // before applying the same strict conservation limit.
    bool retrySpanAwareVolumeReference = false;
    double maximumPhysicalTrimDrift = 0;
    // Disposable linear-mesh construction only. Allows at most 1e-4 native
    // split-volume drift; this is not a tolerance-certified CAD subdivision.
    // The consumer must validate the completed mesh against the unchanged
    // original CAD, including nodes, faces, material, embedding and volume.
    bool meshConstructionCandidate = false;
};
struct OcctPartitionInterface {
    int faceTag = -1;
    std::array<int,2> pieceTags{};
    std::array<int,2> senses{};
    int originalBodyTag = -1;
};
struct OcctPartitionOutcome {
    OcctPartitionReason reason = OcctPartitionReason::KernelFailure;
    std::unique_ptr<OcctModel> model;
    // Dense one-based derived native tags; index zero is a sentinel. Each
    // parent is established by copy/split history and exact incidence, never
    // a distance merge. New cut entities inside material map to dim3/body.
    std::vector<mesh::Classification> parentByVertex;
    std::vector<mesh::Classification> parentByEdge;
    std::vector<mesh::Classification> parentByFace;
    std::vector<int> parentBodyByPiece;
    std::vector<OcctPartitionInterface> interfaces;
    std::uint64_t planes = 0, pieces = 0, historyWork = 0;
    std::uint64_t mappedExteriorFaces = 0, sharedCutFaces = 0;
    double sourceVolume = 0, partitionVolume = 0, relativeVolumeChange = 0;
    double sourceIntegrationError = 0, maximumPieceIntegrationError = 0;
    bool nativeCommonChecked = false;
    double nativeCommonVolume = 0, nativeCommonIntegrationError = 0;
    bool sourceImmutabilityVerified = false;
    double preparationSeconds = 0, splittingSeconds = 0, validationSeconds = 0;
    double totalSeconds = 0;
    bool spanAwareVolumeReferenceUsed = false;
    double legacySourceVolume = 0, legacyPartitionVolume = 0;
    double sourceReferenceStability = 0, partitionReferenceStability = 0;
    OcctTrimToleranceReceipt trimTolerance;
    bool partitioned() const noexcept { return model != nullptr; }
};
// Cooperative checkpoints surround OCCT calls; an individual Boolean/parser
// call is not interruptible and its native allocation peak is not capped.
OcctPartitionOutcome partitionOcctModel(const OcctModel& original,
                                       const OcctPartitionOptions& options = {}) noexcept;
const char* occtPartitionReasonName(OcctPartitionReason reason) noexcept;
const char* occtPartitionPolicyName(OcctPartitionPolicy policy) noexcept;
} // namespace paacs::geo::detail
