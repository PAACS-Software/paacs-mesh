// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private optional-target assembly. No native partition model escapes.
#include "OcctPartition.hpp"
#include "OcctStarMesher.hpp"
#include "OcctNativeValidation.hpp"
#include "OcctRecoveryMesher.hpp"
#include "paacs/geometry/OcctBoundary.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"

namespace paacs::geo::detail {
enum class OcctPartitionMeshReason : std::uint8_t {
    None, InvalidRequest, PartitionFailed, SurfaceFailed, StarFailed,
    ProvenanceLost, InterfaceMismatch, GeometryError, VolumeMismatch,
    InvalidMesh, BudgetExceeded, Cancelled, AllocationFailed, InternalError
};
struct OcctPartitionMeshOptions {
    OcctPartitionOptions partition;
    OcctBoundaryOptions boundary;
    OcctStarOptions star;
    std::uint64_t maximumWork = 5'000'000;
    // Optional per-piece search admission. Failure charges consumed work to
    // maximumWork and leaves unused global work available for another policy.
    // The assembled exterior certificate uses only the global admission.
    std::uint64_t maximumPieceStarWork = UINT64_MAX;
    std::uint64_t maximumNodes = 1'000'000, maximumCells = 5'000'000;
    std::uint64_t maximumOutputBytes = 512ULL * 1024 * 1024;
    // Conservative PAACS logical working-payload admission estimate, excluding
    // allocator metadata and OCCT copy/Boolean/tessellation/projection storage.
    std::uint64_t maximumWorkingBytes = 512ULL * 1024 * 1024;
    double maximumRelativeVolumeError = .01;
    const std::atomic_bool* cancellation = nullptr;
    // Explicit opt-in: at most sixteen failed-piece transverse midpoint cuts.
    // Every attempt starts from the original CAD and spends one shared budget.
    // Zero preserves the existing single-attempt behavior and publication.
    unsigned maximumAdaptiveCuts = 0;
    // One retains the existing serial classifier and work receipts. Additional
    // workers own deep native copies and verify every assembled Tet4 centroid.
    unsigned maximumNativeValidationThreads = 1;
    // Candidate selection uses the original native body as well as the exact
    // piece PLC. Final original-body node/centroid checks remain mandatory.
    bool useNativeKernelAdmission = false;
    // Bounded source-face diagonal changes for opposite coincident facets.
    // Shared native cut triangles and every CAD curve chain remain immutable.
    bool repairSurfaceAliases = false;
    double tetQualityTargetSize = 0;
};
struct OcctPartitionAdaptiveCut {
    int failedPieceTag = -1;
    Vec3 bboxLow{}, bboxHigh{};
    OcctPartitionPlane plane;
    OcctStarReason starReason = OcctStarReason::InternalError;
    std::uint64_t attemptWork = 0, cumulativeWork = 0;
};
struct OcctPartitionMeshOutcome {
    bool published = false, destinationUnchanged = true;
    OcctPartitionMeshReason reason = OcctPartitionMeshReason::InternalError;
    const char* diagnostic = "internal-error";
    OcctPartitionReason partitionReason = OcctPartitionReason::KernelFailure;
    OcctBoundaryReason boundaryReason = OcctBoundaryReason::InternalError;
    OcctStarReason starReason = OcctStarReason::InternalError;
    OcctPartitionPolicy policy = OcctPartitionPolicy::CartesianMidplanes;
    std::uint64_t planes = 0, pieces = 0, interfaces = 0, pairedInterfaceTriangles = 0;
    std::uint64_t nodes = 0, cells = 0, surfaceTriangles = 0;
    std::uint64_t partitionWork = 0, surfaceWork = 0, starWork = 0;
    std::uint64_t assemblyWork = 0, validationWork = 0, sampleWork = 0, workUsed = 0;
    std::uint64_t kernelTrials = 0, exactSurfacePairs = 0;
    int failedPieceTag = -1;
    double cadVolume = 0, meshVolume = 0, relativeVolumeError = 0;
    double partitionVolume = 0, partitionVolumeChange = 0;
    OcctTrimToleranceReceipt trimTolerance;
    double sourceIntegrationError = 0, maximumPieceIntegrationError = 0;
    double maximumBoundaryNodeError = 0, maximumBoundarySampleError = 0;
    double surfaceChordLimit = 0;
    std::uint64_t maximumWorkingPayloadEstimate = 0;
    bool sourceImmutabilityVerified = false;
    OcctBoundaryCertificate exteriorCertificate;
    mesh::LinearVolumeAudit audit;
    algo::MixedMeshQualityReport quality;
    double partitionSeconds = 0, surfaceSeconds = 0, starSeconds = 0;
    double assemblySeconds = 0, validationSeconds = 0, totalSeconds = 0;
    bool failedPieceBoundsAvailable = false;
    bool failedPieceBoundaryCertified = false;
    Vec3 failedPieceBBoxLow{}, failedPieceBBoxHigh{};
    unsigned adaptiveAttempts = 0;
    int adaptivePrincipalFaceTag = -1;
    std::uint64_t adaptiveWork = 0;
    std::vector<OcctPartitionAdaptiveCut> adaptiveCutReports;
    // Receipt for the final attempt's optional parallel centroid validation.
    // Its work is already included in validationWork and workUsed above.
    OcctNativeValidationOutcome nativeCentroidValidation;
    // Included in starWork/workUsed across every adaptive attempt.
    std::uint64_t nativeCandidateQueries = 0;
    std::uint64_t surfaceDiagonalRepairs = 0;
};
// Partitions a copy of one selected original body, tessellates all pieces in
// one shared native session, and star-fills each piece. Actual derived native
// node keys establish joins; coordinates and mapped classes must agree bitwise.
// Only native-declared, exactly paired, oppositely wound cut triangles disappear.
// The compact original exterior (including cavities/genus) is independently
// certified and all original CAD ownership/distance/inclusion/volume gates run
// before one final nonthrowing swap. No partial piece mesh is published.
// Work bounds explicit PAACS stage visits/proposals, not OCCT kernel internals,
// import/audit internals, wall time or process RSS. Cancellation is cooperative.
OcctPartitionMeshOutcome buildOcctPartitionVolume(const OcctModel& original,
    int solidTag, const OcctPartitionMeshOptions&, mesh::Mesh& destination) noexcept;
// Additive recovery entry; the legacy function disables this route. Recovery
// counters are a separately returned subset of the existing work receipt.
OcctPartitionMeshOutcome buildOcctPartitionVolume(const OcctModel& original,
    int solidTag, const OcctPartitionMeshOptions&,
    const OcctPartitionRecoveryOptions&, mesh::Mesh& destination,
    OcctPartitionRecoveryReceipt* recoveryReceipt = nullptr) noexcept;
} // namespace paacs::geo::detail
