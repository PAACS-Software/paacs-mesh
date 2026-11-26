// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Optional native query-only diagnostics. These receipts do not admit an interface,
// weld nodes, change an assembly policy, or certify a mesh.
#include "paacs/geometry/OcctMesher.hpp"
#include <atomic>
#include <optional>
#include <span>
#include <vector>

namespace paacs::geo {
enum class OcctBodyRelation : std::uint8_t {
    Disconnected, FaceContact, EdgeContact, PointContact, MaterialOverlap, Unresolved
};
enum class OcctBodyRelationEvidence : std::uint8_t {
    None, SeparatedBounds, NativeDistance, NativeMaterialCommon, NativeBoundaryCommon,
    NativeBoundarySection
};
enum class OcctBodyRelationIssue : std::uint8_t {
    None, NotStarted, CopyFailed, MaterialQueryFailed, InvalidMaterialResult,
    InvalidOverlapVolume, DistanceQueryFailed, BoundaryQueryFailed,
    InvalidBoundaryResult, InvalidContactMeasure, NearButUnproved,
    SectionQueryFailed, InvalidSectionResult
};
struct OcctBodyPairRequest { int firstSolidTag = -1, secondSolidTag = -1; };
struct OcctBodyPairRelation {
    int firstSolidTag = -1, secondSolidTag = -1;
    OcctBodyRelation relation = OcctBodyRelation::Unresolved;
    OcctBodyRelationEvidence evidence = OcctBodyRelationEvidence::None;
    OcctBodyRelationIssue issue = OcctBodyRelationIssue::NotStarted;
    // Millimetres, square millimetres, and cubic millimetres respectively.
    // A native contact is measured at the original model tolerances. It does
    // not assert matching surface triangles or shared mesh-node identity.
    double classificationTolerance = 0;
    std::optional<double> minimumDistance, nativeContactArea, nativeContactLength;
    std::optional<double> nativeOverlapVolume;
    bool materialQueryComplete = false, boundaryQueryComplete = false;
    std::uint64_t nativeCopies = 0, nativeBooleanCalls = 0, nativeDistanceCalls = 0;
    // Common operations are counted above; the lower-dimensional Section
    // fallback is counted separately so query totals do not double-count it.
    std::uint64_t nativeSectionCalls = 0;
    std::uint64_t nativeShapeValidations = 0, nativeVolumeIntegrations = 0, nativeBoundaryMeasures = 0;
};
struct OcctBodyRelationOptions {
    std::uint64_t maximumPairs = 256, maximumWork = 100'000;
    std::uint64_t maximumRetainedBytes = 64ULL * 1024 * 1024;
    std::uint64_t maximumSourceBytes = 16ULL * 1024 * 1024;
    // Separation/ambiguity threshold only; this does not change OCCT geometry
    // tolerances or Boolean fuzzy values. Effective threshold also includes
    // the two original body tolerances and is reported on every pair.
    double distanceTolerance = 1e-7;
    const std::atomic_bool* cancellation = nullptr;
};
struct OcctBodyRelationReport {
    bool complete = false, sourceAuditAvailable = false, sourceUnmodified = false;
    NativeCadMeshReason reason = NativeCadMeshReason::InternalError;
    const char* diagnostic = "not_started";
    std::uint64_t work = 0, maximumLogicalWorkingBytes = 0;
    std::uint64_t sourceSnapshots = 0;
    std::vector<OcctBodyPairRelation> pairs;
};
// Only the selected, distinct pairs are examined; there is no implicit all-pairs
// traversal. Results use sorted body tags and deterministic pair order. Limits
// bound retained diagnostics/source snapshots and admitted operations, not OCCT
// internal allocations or the running time of one kernel call. Cancellation is
// checked before and after each kernel call. Original BREP bytes are guarded.
OcctBodyRelationReport inspectOcctBodyPairs(const OcctModel&,
    std::span<const OcctBodyPairRequest>, const OcctBodyRelationOptions& = {}) noexcept;
const char* occtBodyRelationName(OcctBodyRelation) noexcept;
const char* occtBodyRelationEvidenceName(OcctBodyRelationEvidence) noexcept;
const char* occtBodyRelationIssueName(OcctBodyRelationIssue) noexcept;
} // namespace paacs::geo
