// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Additive private optional-target recovery. Existing native option/report layouts
// are unchanged. Disabled recovery preserves all previously accepted routes.
#include "OcctStarMesher.hpp"
namespace paacs::geo::detail {
enum class OcctPlcRecoveryReason : std::uint8_t {
    None, InvalidRequest, UnsupportedGeometry, BudgetExceeded, Cancelled,
    AllocationFailed, NativeQueryFailed, InvalidMesh, InternalError
};
struct OcctPartitionRecoveryOptions {
    bool enabled = false;
    bool onlyFinalAdaptiveAttempt = true;
    std::uint64_t maximumPieceWork = 20'000'000;
    unsigned maximumCoreAttempts = 8;
    unsigned maximumQuadFlipsPerAttempt = 16;
    std::uint64_t maximumGeneralWorkPerAttempt = 1'000'000;
    // Explicit construction aid for star-shaped PLCs. A kernel is admitted
    // only after exact PLC halfspace/fan checks and native material checks;
    // it must then classify strictly IN. It adds one volume point, changes
    // no source boundary nodes/classes and retains every final gate. A
    // geometric kernel-search refusal leaves the unseeded recovery available.
    // Resource, cancellation, allocation and native-query failures propagate.
    bool nativeKernelSeed = false;
    unsigned maximumSteinerCandidates = 0;
    // Small immutable piece PLCs can instead use the plain advancing front and
    // connected-region recovery, without constructing a topology workspace.
    bool useAdvancingFront = false;
    // Build a fixed-boundary Tet4 candidate with live constraint recovery.
    // The partition caller checks and repairs native material before joining.
    bool useLiveConstraints = false;
    // Retry a refused live cavity using retained cells and a wider local shell.
    bool useAdaptiveRefill = false;
    // Give the local front time to reduce the cavity before region search.
    // Hard work/output limits and cooperative cancellation remain in force.
    double maximumFrontSeconds = 3.;
};
struct OcctPartitionRecoveryReceipt {
    bool attempted = false;
    OcctPlcRecoveryReason reason = OcctPlcRecoveryReason::None;
    const char* diagnostic = "not-attempted";
    std::uint64_t piecesAttempted = 0, piecesRecovered = 0;
    std::uint64_t coreAttempts = 0, quadFlips = 0;
    // Core construction/import/artifact counters are collected after each
    // retained API phase. General-cavity work reserves its strict solve cap;
    // this is not a total CPU-work, in-phase cancellation or process-RSS cap.
    std::uint64_t workUsed = 0, generalWorkReserved = 0, nativeQueries = 0;
    std::uint64_t steinerNodes = 0;
    int failedPieceTag = -1;
};
struct OcctPlcRecoveryOutcome : OcctPartitionRecoveryReceipt {
    bool published = false;
};
// Updates the private boundary only after a complete recovered volume passes
// exact oriented boundary/positive-volume and BOTH native material checks.
// Native curve segments and all shared-cut triangles remain immutable. The
// caller retains the final Original source/ownership/chord/volume/audit gates.
OcctPlcRecoveryOutcome recoverOcctPlcVolume(mesh::Mesh& boundary,int pieceTag,
    const OcctModel& pieces,const OcctModel& original,int originalBody,
    const std::vector<mesh::Classification>& parentByFace,
    const OcctStarOptions&,const OcctPartitionRecoveryOptions&,
    mesh::Mesh& destination) noexcept;
}
