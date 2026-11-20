// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private optional-adapter route for a simple cylinder contacting an already
// accepted body only along a complete, natively coincident cap circle.
#include "OcctModelInternal.hpp"
#include "OcctStarMesher.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"

namespace paacs::geo::detail {
enum class OcctContactCylinderReason : std::uint8_t {
    None, InvalidRequest, UnsupportedCylinder, ContactCircleNotFound,
    AmbiguousContact, NativeProofFailed, InvalidRim, GeometryError,
    VolumeMismatch, StarFailed, InvalidMesh, BudgetExceeded, Cancelled,
    AllocationFailed, SourceChanged, InternalError
};
struct OcctContactCylinderOptions {
    OcctStarOptions star;
    double maximumChordError = .05;
    double maximumRelativeVolumeError = .01;
    std::uint64_t maximumWork = 2'000'000;
    std::uint64_t maximumNodes = 1'000'000, maximumCells = 5'000'000;
    std::uint64_t maximumOutputBytes = 512ULL*1024*1024;
    std::uint64_t maximumWorkingBytes = 512ULL*1024*1024;
    std::uint64_t maximumSnapshotBytes = 64ULL*1024*1024;
    const std::atomic<bool>* cancellation = nullptr;
};
struct OcctContactCylinderOutcome {
    bool published = false, destinationUnchanged = true;
    bool nativeCircleProved = false, nativeCylinderProved = false;
    bool sourceImmutabilityVerified = false;
    OcctContactCylinderReason reason = OcctContactCylinderReason::InternalError;
    OcctStarReason starReason = OcctStarReason::InternalError;
    int sourceEdgeTag = -1, targetEdgeTag = -1;
    int sourceFaceTag = -1, contactCapFaceTag = -1;
    std::uint64_t ringNodes = 0, workUsed = 0;
    double maximumBoundaryNodeError = 0, maximumBoundarySampleError = 0;
    double cadVolume = 0, meshVolume = 0, relativeVolumeError = 0;
    mesh::LinearVolumeAudit audit;
    algo::MixedMeshQualityReport quality;
};
// accepted must have stable original native classifications from this model.
// Source rim coordinates are copied bit-for-bit into DISTINCT target-body node
// IDs; no coordinate welding occurs. Complete edge/vertex Boolean Common and
// native incidence establish correspondence; proximity alone never does.
// The caller must retain its original whole interbody contact audit after
// appending the returned target body. This helper does not approve area contact.
// Work counts PAACS visits and explicit native calls, excluding kernel internals;
// logical payload bounds exclude OCCT allocations/process RSS. Cancellation is
// cooperative. Native operations use a deep copy and publication is one swap.
OcctContactCylinderOutcome buildOcctContactCylinder(const OcctModel&, int targetBodyTag,
    const mesh::Mesh& accepted, const OcctContactCylinderOptions&,
    mesh::Mesh& destination) noexcept;
const char* occtContactCylinderReasonName(OcctContactCylinderReason) noexcept;
} // namespace paacs::geo::detail
