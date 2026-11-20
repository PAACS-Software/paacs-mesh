// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private optional adapter route. No native prism or cap snapshot escapes.
#include "OcctStarMesher.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/geometry/OcctBoundary.hpp"

namespace paacs::geo::detail {
enum class OcctNativeSweepReason : std::uint8_t {
    None,
    InvalidRequest,
    UnsupportedGeometry,
    RecognitionFailed,
    SurfaceFailed,
    ProvenanceLost,
    SweepFailed,
    GeometryError,
    VolumeMismatch,
    InvalidMesh,
    QualityRequirement,
    BudgetExceeded,
    Cancelled,
    AllocationFailed,
    SourceMutation,
    KernelFailure,
    InternalError
};
struct OcctNativeSweepOptions {
    OcctBoundaryOptions boundary;
    // Explicit uniform normal layers. Tangential cap topology is native.
    unsigned layers = 1;
    // Positive additionally derives at least ceil(height/limit) layers after
    // native recognition. Zero retains explicit layers; hard admission4096.
    double maximumLayerThickness = 0;
    std::uint64_t maximumWork = 20'000'000;
    std::uint64_t maximumNodes = 1'000'000, maximumCells = 5'000'000;
    std::uint64_t maximumOutputBytes = 512ULL * 1024 * 1024;
    std::uint64_t maximumWorkingBytes = 512ULL * 1024 * 1024;
    std::uint64_t maximumSnapshotBytes = 64ULL * 1024 * 1024;
    double maximumRelativeRecognitionError = 1e-6;
    double maximumRelativeVolumeError = .01;
    double minimumQuality = 0, warningQuality = .25;
    mesh::LinearVolumeAuditLimits auditLimits{};
    const std::atomic_bool *cancellation = nullptr;
};
struct OcctNativeSweepOutcome {
    bool published = false, destinationUnchanged = true;
    OcctNativeSweepReason reason = OcctNativeSweepReason::InternalError;
    const char *diagnostic = "internal-error";
    OcctBoundaryReason boundaryReason = OcctBoundaryReason::InternalError;
    int baseCapTag = -1, targetCapTag = -1;
    Vec3 direction{};
    double height = 0;
    unsigned layers = 0;
    std::uint64_t capNodes = 0, capTriangles = 0, nodes = 0, cells = 0;
    std::uint64_t recognitionWork = 0, surfaceWork = 0, sweepWork = 0;
    std::uint64_t validationWork = 0, workUsed = 0, nativeInteriorQueries = 0;
    std::uint64_t maximumWorkingPayloadEstimate = 0;
    double cadVolume = 0, expectedVolume = 0, commonVolume = 0, meshVolume = 0;
    double expectedRelativeError = 0, commonRelativeError = 0;
    double integrationStability = 0, relativeVolumeError = 0;
    double maximumBoundaryNodeError = 0, maximumBoundarySampleError = 0;
    bool sourceImmutabilityVerified = false;
    OcctBoundaryCertificate exteriorCertificate;
    mesh::LinearVolumeAudit audit;
    algo::MixedMeshQualityReport quality;
    double recognitionSeconds = 0, surfaceSeconds = 0, sweepSeconds = 0;
    double validationSeconds = 0, totalSeconds = 0;
};
// Two parallel planar caps, straight parallel side geometry, and a native
// trimmed-cap prism whose stable span-aware mass and non-destructive Common
// agree with the entire original solid are mandatory. Holes remain present.
// A native cap triangulation feeds the retained planar layered sweep and
// conforming pulling Tet4 conversion. Original topology supplies every tag;
// there is no coordinate identity merge. Exact exterior embedding, all original
// face coverage, original-CAD node/facet distances and every Tet4 centroid are
// checked before one nonthrowing swap. Refusal leaves destination untouched.
// Work bounds PAACS visits/native calls, not native Boolean/integration or
// tessellation internals, process memory or wall time. Cancellation is checked
// between calls. Working bytes admit logical payload, not native allocations.
OcctNativeSweepOutcome buildOcctNativeSweepVolume(const OcctModel &, int solidTag,
                                                  const OcctNativeSweepOptions &,
                                                  mesh::Mesh &destination) noexcept;
const char *occtNativeSweepReasonName(OcctNativeSweepReason) noexcept;
} // namespace paacs::geo::detail
