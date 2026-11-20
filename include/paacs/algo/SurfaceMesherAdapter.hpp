// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// BND-2A Increment 6: bounded full T1-A surface-set adapter. This remains
// separate from the normal FaceMesher::meshFaces() route.

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <type_traits>
#include <tuple>
#include <vector>

#include "paacs/bnd/SteinerRecovery.hpp"
#include "paacs/mesh/Mesh.hpp"

namespace paacs::algo {

namespace detail {
struct PreparedSurfaceCandidate;
}

struct ArtifactStamp {
    std::uint64_t persistentId = 0;
    std::uint64_t revision = 0;
    friend bool operator==(const ArtifactStamp&, const ArtifactStamp&) = default;
    friend bool operator<(const ArtifactStamp& a, const ArtifactStamp& b) {
        return a.persistentId < b.persistentId ||
               (a.persistentId == b.persistentId && a.revision < b.revision);
    }
};

enum class SurfaceNodeIdentityDomain : std::uint8_t {
    OriginalDiscretization,
    RefinedCurve,
    RefinedFaceEdge,
    RefinedFacetInterior,
    // Face-owned point placed before a parent surface triangle exists.
    // parentOrdinal/parentSubordinal are the one-based lattice coordinates.
    InterfacePrefill
};

// Exact persistent identity for a surface discretization point.  Parameter
// values are stored as their canonical binary64 bit patterns; they are never
// reduced to a hash.  Curve refinement additionally retains the exact dyadic
// path inside its original discretization segment.
struct SurfaceNodeIdentity {
    SurfaceNodeIdentityDomain domain =
        SurfaceNodeIdentityDomain::OriginalDiscretization;
    ArtifactStamp parentArtifact{};
    std::array<Index,3> parentFacetGlobalNodes{
        kInvalidIndex,kInvalidIndex,kInvalidIndex};
    std::uint64_t parentOrdinal = 0;
    std::uint64_t parentSubordinal = 0;
    std::uint64_t parameterUBits = 0;
    std::uint64_t parameterVBits = 0;
    std::uint64_t dyadicNumerator = 0;
    std::uint32_t dyadicLevel = 0;

    friend bool operator==(const SurfaceNodeIdentity&,
                           const SurfaceNodeIdentity&) = default;
    friend bool operator<(const SurfaceNodeIdentity& a,
                          const SurfaceNodeIdentity& b) noexcept {
        return std::tie(a.domain,a.parentArtifact,a.parentFacetGlobalNodes,
                        a.parentOrdinal,
                        a.parentSubordinal,a.parameterUBits,a.parameterVBits,
                        a.dyadicNumerator,a.dyadicLevel) <
               std::tie(b.domain,b.parentArtifact,b.parentFacetGlobalNodes,
                        b.parentOrdinal,
                        b.parentSubordinal,b.parameterUBits,b.parameterVBits,
                        b.dyadicNumerator,b.dyadicLevel);
    }
};

// A persistent point identity has one canonical representation for its
// domain.  Inactive fields must retain their documented zero/sentinel value;
// callers may not smuggle presentation-dependent data into ordering.
bool validSurfaceNodeIdentity(const SurfaceNodeIdentity& identity) noexcept;

struct SharedCurveNodeKey {
    std::uint64_t curvePersistentId = 0;
    std::uint64_t curveRevision = 0;
    std::uint64_t nodeOrdinal = 0;
    SurfaceNodeIdentity identity{};
    friend bool operator==(const SharedCurveNodeKey&,
                           const SharedCurveNodeKey&) = default;
    friend bool operator<(const SharedCurveNodeKey& a,
                          const SharedCurveNodeKey& b) {
        if (a.curvePersistentId != b.curvePersistentId)
            return a.curvePersistentId < b.curvePersistentId;
        if (a.curveRevision != b.curveRevision)
            return a.curveRevision < b.curveRevision;
        if (!(a.identity == b.identity)) return a.identity < b.identity;
        return a.nodeOrdinal < b.nodeOrdinal;
    }
};

struct SharedCurveNodeBinding {
    SharedCurveNodeKey key{};
    Index existingGlobalIndex = kInvalidIndex;
    mesh::Classification classification{};
};

struct CandidateProviderStamp {
    std::uint64_t persistentId = 0;
    std::uint64_t revision = 0;
    friend bool operator==(const CandidateProviderStamp&,
                           const CandidateProviderStamp&) = default;
};

using FacePersistentId = std::uint64_t;
using FaceLocalVertexId = std::uint64_t;
using SurfaceChainId = std::uint64_t;

struct SurfaceVertexUse {
    FaceLocalVertexId faceLocalVertexId = 0;
    bnd::ParamPoint parameter{};
    mesh::Classification classification{};
};

struct SurfaceChainVertex {
    FaceLocalVertexId faceLocalVertexId = 0;
    std::optional<SharedCurveNodeKey> sharedKey;
};

enum class SurfaceChainRole {
    OuterBoundary,
    HoleBoundary,
    EmbeddedConstraint
};

struct SurfaceChain {
    SurfaceChainId chainPersistentId = 0;
    SurfaceChainRole role = SurfaceChainRole::EmbeddedConstraint;
    mesh::Classification classification{};
    bool closed = false;
    std::vector<SurfaceChainVertex> uses;
    std::vector<bnd::ConstraintEntityId> segmentEntityIds;
    std::optional<ArtifactStamp> owningCurve;
};

struct SurfaceLoop {
    SurfaceChainRole role = SurfaceChainRole::OuterBoundary;
    std::vector<SurfaceChain> pieces;
};

struct SurfacePhysicalPolicy {
    double minPhysicalDoubleArea = 1.0e-14;
    double maxVertexMappingError = 0.0;
    double maxMappingDeviation = 1.0e-8;
    double minNormalAlignment = 0.99;
    double maxArtifactPositionError = 0.0;
};

struct SurfaceRefinedNodeUse {
    SurfaceNodeIdentity identity{};
    bnd::ParamPoint parameter{};
    mesh::Classification classification{};
};

struct SurfaceFaceRequest {
    ArtifactStamp faceStamp{};
    bnd::ChartId chart = 0;
    mesh::Classification faceClass{2, -1};
    int desiredNormalSign = 1;
    SurfacePhysicalPolicy physicalPolicy{};
    std::vector<SurfaceVertexUse> vertices;
    std::vector<SurfaceRefinedNodeUse> refinedNodes;
    SurfaceLoop outer{};
    std::vector<SurfaceLoop> holes;
    std::vector<SurfaceChain> embedded;
    bnd::GeometryProviderStamp expectedGeometryStamp{};
    ArtifactStamp expectedChartCertificateStamp{};
    CandidateProviderStamp expectedCandidateProviderStamp{};
};

struct SurfaceSetRequest {
    ArtifactStamp setStamp{};
    std::vector<SurfaceFaceRequest> faces;
    std::vector<SharedCurveNodeBinding> sharedBindings;
};

struct ChartCertificate {
    bnd::GeometryProviderStamp providerStamp{};
    bnd::ChartId chart = 0;
    ArtifactStamp certificateStamp{};
    std::vector<Vec2> certifiedRegionBoundary;
    int chartPositiveToNormalSign = 1;
    bool regular = true;
    bool periodic = false;
    bool seam = false;
    friend bool operator==(const ChartCertificate&, const ChartCertificate&) =
        default;
};

class SurfaceGeometryProvider : public bnd::GeometryProvider {
public:
    ~SurfaceGeometryProvider() override = default;
    virtual Vec3 normalAt(bnd::ChartId chart, const Vec2& uv) const = 0;
    virtual ChartCertificate chartCertificate(bnd::ChartId chart) const = 0;
};

struct SurfaceCandidateNode {
    topo::NodePool::NodeData node{};
    bnd::ParamPoint parameter{};
    bool frameOnly = false;
    std::optional<SurfaceNodeIdentity> persistentIdentity;
};

struct SurfaceCandidateCell {
    std::array<std::size_t, 3> nodes{};
    mesh::Classification classification{};
};

struct SurfaceCandidate {
    CandidateProviderStamp providerStamp{};
    ArtifactStamp faceStamp{};
    bnd::ChartId chart = 0;
    std::vector<SurfaceCandidateNode> nodes;
    std::vector<SurfaceCandidateCell> cells;
    std::map<FaceLocalVertexId, std::size_t> declaredVertexMap;
    std::vector<std::size_t> frameBoundary;
    std::size_t exteriorSeedCell = 0;
    // Private producer-to-adapter handoff. Third-party providers leave this
    // empty and take the unchanged fully validating construction path.
    std::shared_ptr<detail::PreparedSurfaceCandidate> preparedRecovery;
};

class SurfaceCandidateProvider {
public:
    virtual ~SurfaceCandidateProvider() = default;
    virtual CandidateProviderStamp stamp() const = 0;
    virtual SurfaceCandidate build(const SurfaceFaceRequest& request) const = 0;
};

enum class SurfaceCandidateDetail {
    None, MalformedCandidate, MissingExteriorSeed, DisconnectedFrame,
    PretrimmedDomain, ProviderStampChanged
};
enum class SurfaceDomainDetail {
    None, OneSidedDomainBoundary, FloodFailed, RetainedRebuildFailed
};
enum class SurfaceAllocationPhase {
    None, GeometryProvider, Candidate, Recovery, DomainExtraction, Audit,
    SetMapping, PublicationPreparation
};
enum class SurfaceReason {
    None,
    InvalidLoopTopology,
    InvalidConstraintEmbedding,
    InvalidSegmentIdentity,
    InvalidPhysicalPolicy,
    StaleGeometryProvider,
    StaleCandidateProvider,
    StaleFaceArtifact,
    StaleCurveArtifact,
    DuplicateFaceIdentity,
    DuplicateSharedBinding,
    MissingSharedBinding,
    UnusedSharedBinding,
    UnknownSharedGlobalIndex,
    SharedClassificationConflict,
    SharedIdentityConflict,
    SharedMappedPositionMismatch,
    InvalidClassification,
    NonfiniteParam,
    NonfiniteMappedPosition,
    NonfiniteMappedNormal,
    ProviderDomainRejected,
    UnsupportedPeriodicChart,
    UnsupportedSeamChart,
    UnsupportedSingularChart,
    BoundarySteinerProhibited,
    CurveArtifactSubdivisionRequired,
    MinimumFeatureReached,
    SplitBudgetExceeded,
    GeometryProviderError,
    CandidateProviderError,
    CandidateValidationFailed,
    RecoveryStopped,
    DomainExtractionFailed,
    PhysicalOrientationFailed,
    LocalAuditFailed,
    PublicationValidationFailed,
    AllocationFailed,
    InternalError
};

enum class SurfaceFaceStatus { Ready, Rejected, Refused, Failed };
enum class SurfaceSetStatus { Published, Rejected, Refused, Failed };

struct SurfaceFailure {
    SurfaceReason reason = SurfaceReason::None;
    SurfaceCandidateDetail candidate = SurfaceCandidateDetail::None;
    SurfaceDomainDetail domain = SurfaceDomainDetail::None;
    SurfaceAllocationPhase allocation = SurfaceAllocationPhase::None;
    bnd::SteinerReject bndReject = bnd::SteinerReject::None;
    bnd::SteinerRefuse bndRefuse = bnd::SteinerRefuse::None;
    bnd::SteinerFail bndFail = bnd::SteinerFail::None;
};

struct SurfaceSegmentOrigin {
    bnd::ConstraintEntityId bndEntityId = 0;
    SurfaceChainId chainPersistentId = 0;
    std::size_t originalSegmentOrdinal = 0;
    std::size_t normalizedSegmentOrdinal = 0;
    std::array<SurfaceChainVertex, 2> endpoints{};
    mesh::Classification classification{};
};

struct PersistentConstraintEdge {
    ArtifactStamp faceStamp{};
    bnd::ConstraintEntityId bndEntityId = 0;
    SurfaceChainId chainPersistentId = 0;
    std::size_t originalSegmentOrdinal = 0;
    std::size_t normalizedSegmentOrdinal = 0;
    SurfaceChainRole role = SurfaceChainRole::EmbeddedConstraint;
    std::array<topo::NodeHandle, 2> localEndpoints{};
    std::array<std::optional<SharedCurveNodeKey>, 2> endpointSharedKeys{};
    mesh::Classification classification{};
    bool required = true;
    std::optional<std::array<Index, 2>> publishedGlobalEndpoints;
};

struct SurfaceTopologyAudit {
    bool workspaceReady = false;
    bool connectedGenusZero = false;
    bool boundaryCyclesExact = false;
    bool nodeLinksValid = false;
    bool constraintsPersistent = false;
    bool classificationsPreserved = false;
    bool sharedMappingsPreserved = false;
    bool noFrameSurvivors = false;
    bool noDuplicateTriangles = false;
    bool physicalOrientation = false;
    std::size_t boundaryCycleCount = 0;
    std::int64_t eulerCharacteristic = 0;
    bool passed() const noexcept;
};

struct SurfaceFaceCell {
    std::array<std::size_t, 3> nodes{};
    mesh::Classification classification{};
};

struct SurfaceFaceResult {
    ArtifactStamp faceStamp{};
    bnd::ChartId chart = 0;
    mesh::Classification faceClass{};
    int desiredNormalSign = 1;
    SurfacePhysicalPolicy physicalPolicy{};
    bnd::GeometryProviderStamp expectedGeometryStamp{};
    ChartCertificate chartCertificate{};
    std::vector<topo::NodePool::NodeData> localNodes;
    std::vector<bnd::AuthoritativeParam> localParams;
    // Evaluated mapping coordinates paired one-to-one with localParams. The
    // symbolic token remains authoritative; this cache is audited against it.
    std::vector<Vec2> localMappingParams;
    // Sparse by design: accepted v3 artifacts remain byte/value compatible.
    // Only additive v4 refined nodes have entries.
    std::map<std::size_t, SurfaceNodeIdentity> refinedIdentityByLocalNode;
    std::vector<topo::NodeHandle> localHandles;
    std::vector<SurfaceFaceCell> localCells;
    std::map<std::size_t, Index> existingGlobalIndex;
    std::vector<Index> publishedGlobalByLocalNode;
    std::vector<Index> publishedGlobalCellByLocalCell;
    std::vector<std::vector<std::size_t>> boundaryLoops;
    std::vector<SurfaceChain> normalizedChains;
    std::map<bnd::ConstraintEntityId, SurfaceSegmentOrigin> segmentOrigins;
    std::vector<PersistentConstraintEdge> persistentEdges;
    bnd::SteinerSegmentsReport protectedRecovery;
    bnd::SteinerSegmentsReport localRecovery;
    // True only when an exact producer receipt was consumed. This is
    // diagnostic evidence; it never changes the published mesh contract.
    bool preparedRecoveryConsumed = false;
    SurfaceTopologyAudit audit{};
    std::unique_ptr<topo::TopologyWorkspace> localWorkspace;
};

struct SurfaceSetResult {
    ArtifactStamp setStamp{};
    std::vector<std::unique_ptr<SurfaceFaceResult>> faces;
    std::map<SurfaceNodeIdentity,Index> publishedGlobalByRefinedIdentity;
    std::map<SharedCurveNodeKey, SharedCurveNodeBinding> sharedLedger;
};

struct SurfaceFaceRow {
    FacePersistentId facePersistentId = 0;
    SurfaceFaceStatus status = SurfaceFaceStatus::Failed;
    SurfaceFailure failure{};
    std::unique_ptr<SurfaceFaceResult> readyResult;
};

struct SurfaceSetOutcome {
    SurfaceSetStatus setStatus = SurfaceSetStatus::Failed;
    SurfaceFailure setFailure{};
    bool reportComplete = false;
    std::vector<SurfaceFaceRow> processed;
    std::vector<FacePersistentId> unattempted;
    std::unique_ptr<SurfaceSetResult> retainedArtifact;
};

static_assert(std::is_nothrow_move_constructible_v<SurfaceFaceRow>);
static_assert(std::is_nothrow_move_constructible_v<std::unique_ptr<SurfaceSetResult>>);
static_assert(std::is_nothrow_move_constructible_v<std::vector<SurfaceFaceRow>>);
static_assert(std::is_nothrow_move_constructible_v<std::vector<FacePersistentId>>);
static_assert(std::is_nothrow_move_constructible_v<SurfaceSetOutcome>);
static_assert(std::is_nothrow_move_assignable_v<SurfaceSetOutcome>);

// One geometry and one candidate provider serve the complete set.
SurfaceSetOutcome buildAndPublishSurfaceSet(
    const SurfaceSetRequest& request, mesh::Mesh& destination,
    const SurfaceGeometryProvider& geometry,
    const SurfaceCandidateProvider& candidateProvider,
    const bnd::SteinerPolicy& localSteinerPolicy = {}) noexcept;

// Allocation-sweep retry seam for a retained unpublished artifact.
SurfaceFailure publishSurfaceSet(mesh::Mesh& destination,
                                 SurfaceSetResult& artifact,
                                 const SurfaceGeometryProvider& geometry,
                                 bool& published) noexcept;

#if defined(PAACS_TOPO_TESTING)
struct CandidateBoundaryRecoveryTestingOutcome {
    bool complete = false;
    bnd::SteinerReject reject = bnd::SteinerReject::None;
    bnd::SteinerRefuse refuse = bnd::SteinerRefuse::None;
    bnd::SteinerFail fail = bnd::SteinerFail::None;
};

CandidateBoundaryRecoveryTestingOutcome
recoverCandidateBoundaryForTesting(
    const SurfaceFaceRequest& request, SurfaceCandidate& candidate,
    const SurfaceGeometryProvider& geometry);

inline SurfaceDomainDetail testSurfaceBoundaryOwnerGate(
    std::size_t ownerCount, predicates::Sign first,
    predicates::Sign second) noexcept {
    return ownerCount == 2 && first != predicates::Sign::Zero &&
                   second != predicates::Sign::Zero && first != second
               ? SurfaceDomainDetail::None
               : SurfaceDomainDetail::OneSidedDomainBoundary;
}
#endif

}  // namespace paacs::algo
