// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T2-B: additive boundary-conforming-v3 model route.

#include "paacs/algo/ProductionSurfaceMesher.hpp"
#include "paacs/algo/SurfaceMesherAdapter.hpp"
#include "paacs/algo/TetCandidate.hpp"
#include "paacs/bnd/TetFacetSetRecovery.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace paacs::algo {

namespace boundary_conforming_detail {
// Allocation probes observe the final publication tail through this
// thread-local flag. It has no effect on routing or production outcomes.
inline bool& publicationTailActive() noexcept {
    static thread_local bool active = false;
    return active;
}
#if defined(PAACS_TOPO_TESTING)
// Probe-only byte marker. It is deliberately primitive TLS state so reading
// it from operator new and a terminate handler cannot allocate.
inline std::uint8_t& testingAllocationPhaseMarker() noexcept {
    static thread_local std::uint8_t phase = 0;
    return phase;
}
#endif
}

inline constexpr const char* kBoundaryConformingV3RouteId =
    "boundary-conforming-v3";
inline constexpr std::uint32_t kBoundaryConformingV3RouteVersion = 3;

struct VolumeBoundaryProviderStamp {
    std::uint64_t persistentId = 0;
    std::uint64_t revision = 0;
    friend bool operator==(const VolumeBoundaryProviderStamp&,
                           const VolumeBoundaryProviderStamp&) = default;
};

enum class VolumeProviderAvailability : std::uint8_t {
    Available,
    MissingAuthoritativeOrientation,
    RecoveryEnvelopeUnavailable
};

struct VolumeFaceOwnership {
    VolumeBoundaryProviderStamp providerStamp{};
    FacePersistentId facePersistentId = 0;
    VolumeProviderAvailability availability =
        VolumeProviderAvailability::Available;
    std::vector<bnd::VolumeRegionUse> owners;
};

struct ProviderRecoveryEnvelope {
    VolumeBoundaryProviderStamp providerStamp{};
    VolumeProviderAvailability availability =
        VolumeProviderAvailability::Available;
    std::vector<Vec3> guardCoordinates;
    std::vector<std::array<std::size_t, 3>> outwardHullFaces;
};

class VolumeBoundaryOwnershipProvider {
public:
    virtual ~VolumeBoundaryOwnershipProvider() = default;
    virtual VolumeBoundaryProviderStamp stamp() const noexcept = 0;
    virtual VolumeFaceOwnership ownership(FacePersistentId face) const = 0;
    virtual ProviderRecoveryEnvelope recoveryEnvelope(
        std::span<const Vec3> constructionCoordinates) const = 0;
};

enum class RecoveryEnvelopeMode : std::uint8_t {
    DerivedFiniteBox,
    ProviderSuppliedConvexHull
};

enum class BoundaryConformingV3Status : std::uint8_t {
    Published,
    Rejected,
    Unsupported,
    Refused,
    Failed
};

enum class BoundaryConformingV3Phase : std::uint8_t {
    None,
    Preflight,
    SurfaceConstruction,
    SurfaceMapping,
    Ownership,
    RecoveryEnvelope,
    TetConstruction,
    CandidateValidation,
    AuthoritativeImport,
    ArtifactIndexBuild,
    RecoveryContextBuild,
    IndexedPreflight,
    EdgeRouting,
    EdgeDiscovery,
    EdgeCandidate,
    EdgePreparation,
    EdgeStateCopy,
    FacetRouting,
    FacetDiscovery,
    FacetCandidate,
    FacetPreparation,
    FacetStateCopy,
    IndexUpdate,
    DomainExtraction,
    RegionAssignment,
    AuthoritativeRebuild,
    BoundaryAudit,
    MeshBuild,
    ResultAssembly,
    Publication
};

enum class VolumeProviderCheckpoint : std::uint8_t {
    None,
    Entry,
    PreOwnershipTable,
    PostOwnershipTable,
    PreEnvelope,
    PostEnvelope,
    PreEdgeRecovery,
    PreDomainExtraction,
    PreBoundaryAudit
};

enum class BoundaryConformingV3Reason : std::uint8_t {
    None,
    InvalidRequest,
    InvalidSizingPolicy,
    SizingEvaluationFailed,
    SurfaceConstructionRejected,
    SurfaceConstructionRefused,
    SurfaceConstructionFailed,
    StaleSurfaceArtifact,
    SurfaceNodeMappingIncomplete,
    SurfaceCellMappingIncomplete,
    SurfaceNodeMappingAlias,
    SurfaceCellMappingAlias,
    SurfaceNodePayloadMismatch,
    SurfaceCellPayloadMismatch,
    DuplicateFacetIdentity,
    DuplicateFacetGeometry,
    InvalidFacetClassification,
    InvalidArtifactEdge,
    ArtifactEdgeClassificationConflict,
    InvalidOwnerMultiplicity,
    InvalidOwnerSense,
    InvalidOwnerSide,
    NonClosedShellComponent,
    NonOrientableShellComponent,
    NonManifoldShellComponent,
    ArtifactContactConflict,
    StaleOwnershipProvider,
    MultipleOwnershipProviders,
    MissingOwnershipOrientation,
    ProviderEnvelopeUnavailable,
    OwnershipProviderStampChanged,
    OwnershipProviderException,
    OwnershipProviderContractViolation,
    EnvelopeRangeUnrepresentable,
    InvalidRecoveryEnvelope,
    TetConstructionFailed,
    InvalidTetCandidate,
    AuthoritativeImportFailed,
    ArtifactIndexInvalid,
    RecoveryContextInvalid,
    StaleRecoveryContext,
    MissingArtifactEdge,
    EdgeMetadataConflict,
    ProtectedConstraintConflict,
    UnsupportedEdgeDegeneracy,
    NonStarEdgeCavity,
    InvalidEdgeCandidate,
    UnsupportedFacetDegeneracy,
    InvalidFacetCandidate,
    RecoveryConflictCycle,
    PreparationConflict,
    IndexUpdateFailed,
    UnseededDomainComponent,
    RegionLabelConflict,
    DisconnectedRegionTagReuse,
    UnsupportedRegionJunction,
    UnexpectedRegionInterface,
    AuthoritativeRebuildFailed,
    BoundaryBijectionFailed,
    MeshBuildFailed,
    AllocationFailed,
    InternalError,
    ConstructionMetricsOverflow,
    ResourceBudgetReached
};

enum class VolumeWorkKind : std::uint8_t { None, ArtifactEdge, Facet };
enum class BoundaryConformingV3WorkStatus : std::uint8_t {
    AlreadyPresent,
    Marked,
    Recovered,
    Rejected,
    Unsupported,
    Refused,
    Failed
};

struct BoundaryConformingV3WorkKey {
    VolumeWorkKind kind = VolumeWorkKind::None;
    bnd::VolumeArtifactEdgeKey edgeKey{};
    bnd::VolumeFacetKey facetKey{};
};

struct BoundaryConformingV3Failure {
    BoundaryConformingV3Reason reason = BoundaryConformingV3Reason::None;
    BoundaryConformingV3Phase phase = BoundaryConformingV3Phase::None;
    VolumeProviderCheckpoint providerCheckpoint =
        VolumeProviderCheckpoint::None;
    VolumeWorkKind workKind = VolumeWorkKind::None;
    bnd::VolumeArtifactEdgeKey edgeKey{};
    bnd::VolumeFacetKey facetKey{};
    std::size_t scheduleOrdinal = 0;
    bnd::TetEdgeRecoveryOutcome edgeOutcome{};
    bnd::TetFacetRecoveryOutcome facetOutcome{};
    std::uint64_t exactWitnessDigest = 0;
};

struct BoundaryConformingV3ReportRow {
    BoundaryConformingV3WorkKey key{};
    BoundaryConformingV3WorkStatus status =
        BoundaryConformingV3WorkStatus::Failed;
    std::size_t scheduleOrdinal = 0;
    BoundaryConformingV3Failure failure{};
};

struct ArtifactPreflightMetrics {
    std::uint64_t edgeEntries = 0;
    std::uint64_t facetEntries = 0;
    std::uint64_t treeNodes = 0;
    std::uint64_t treeDepth = 0;
    std::uint64_t axisSelectionVisits = 0;
    std::uint64_t sortComparisons = 0;
    std::uint64_t nodePairVisits = 0;
    std::uint64_t aabbTests = 0;
    std::uint64_t edgeEdgePairs = 0;
    std::uint64_t edgeFacetPairs = 0;
    std::uint64_t facetFacetPairs = 0;
    std::uint64_t exactEdgeEdgeTests = 0;
    std::uint64_t exactEdgeFacetTests = 0;
    std::uint64_t exactFacetFacetTests = 0;
    std::uint64_t legalSharedContacts = 0;
    std::uint64_t rejectedContacts = 0;
    std::size_t scratchHighWaterBytes = 0;
    std::uint64_t allPairsFallbacks = 0;
    bool runComplete = false;
};

struct TetRecoveryMetrics {
    std::uint64_t contextNodes = 0, contextCells = 0, contextEdges = 0,
                  contextFaces = 0, contextAuditVisits = 0,
                  indexBuildComparisons = 0;
    std::size_t indexBytes = 0;
    std::uint64_t edgeAabbQueries = 0, facetAabbQueries = 0,
                  indexNodesVisited = 0, returnedTetCandidates = 0,
                  segmentTetTests = 0, targetTetTests = 0,
                  candidateSurvivorTests = 0, segmentWalkCells = 0,
                  linkCycleHighWater = 0, edgeCavityCommits = 0,
                  predictedEdgeMatches = 0, completeStarQueries = 0,
                  closureInsertions = 0, apexPairs = 0, facetCommits = 0,
                  indexRemovals = 0, indexInsertions = 0,
                  indexRebalances = 0, indexHighWater = 0,
                  routeAssessments = 0, dependencyEdges = 0,
                  dependencyInvalidations = 0, selectedRebuilds = 0,
                  batchPreflightAssessments = 0,
                  alreadyPresentBatchEntities = 0,
                  protectedImportFallbacks = 0,
                  protectedConflicts = 0, conflictSccVisits = 0,
                  fullStageRescans = 0,
                  productionGlobalScanFallbacks = 0;
    std::size_t scratchHighWaterBytes = 0;
    bool runComplete = false;
};

struct VolumeExtractionMetrics {
    std::uint64_t barrierFaces = 0, adjacencyFacesVisited = 0,
                  cellsVisited = 0, components = 0,
                  exteriorComponents = 0, retainedRegionComponents = 0,
                  labelSeeds = 0, labelPropagations = 0, retainedCells = 0,
                  discardedCells = 0, retainedNodes = 0,
                  discardedNodes = 0;
    std::size_t authoritativeRebuildBytes = 0;
    std::uint64_t rebuildCellInsertions = 0;
    topo::EdgePoolCompleteBuildMetrics edgeBuild{};
    topo::FacePoolCompleteBuildMetrics faceBuild{};
    std::uint64_t centroidClassifierCalls = 0, boundaryPeelingCalls = 0,
                  hullCapPeelingCalls = 0;
    bool runComplete = false;
};

struct BoundaryConformingV3SurfaceConstructionMetrics {
    std::uint64_t faceRows = 0, boundaryNodes = 0,
                  interiorCandidates = 0, spacingQueries = 0,
                  spacingBucketsVisited = 0, spacingDistanceTests = 0,
                  constraintPointTests = 0, protectedSegments = 0,
                  recoveredProtectedSegments = 0;
    Delaunay2DScalableMetrics scalable{};
    bool runComplete = false;
};

struct BoundaryConformingV3ConstructionMetrics {
    BoundaryConformingV3SurfaceConstructionMetrics surface{};
    TetCandidateMetrics tetCandidate{};
    std::uint64_t importAttempts = 0, importNodeVisits = 0,
                  importCellVisits = 0,
                  importEdgeValidationCellVisits = 0,
                  importEdgeBuildCellVisits = 0, importEdgeUses = 0,
                  importUniqueEdges = 0,
                  importFaceValidationCellVisits = 0,
                  importFaceBuildCellVisits = 0, importFaceUses = 0,
                  importUniqueFaces = 0, importEdgeAudits = 0,
                  importFaceAudits = 0, importPerCellGlobalAudits = 0,
                  importCandidateValidationCalls = 0,
                  importCandidateValidationPointRows = 0,
                  importCandidateValidationCellRows = 0,
                  importCandidateValidationFaceUses = 0,
                  importCandidateValidationFaceRows = 0,
                  importCandidateValidationBoundaryRows = 0,
                  importCandidateValidationDualNeighbourVisits = 0,
                  importCandidateValidationBoundaryEdgeUses = 0,
                  importCandidateValidationBoundaryNeighbourVisits = 0,
                  importCandidateValidationHullTests = 0,
                  importCandidateValidationSlotsInitialized = 0,
                  importSelectionRowsCopied = 0,
                  importSelectionSortComparisons = 0,
                  importSelectionDuplicateComparisons = 0,
                  importSelectionCellVisits = 0,
                  importSelectionVertexVisits = 0,
                  importSelectionPointSlotsInitialized = 0,
                  importGlobalMapRows = 0,
                  importGlobalMapComparisons = 0,
                  importGlobalMapInsertions = 0,
                  importLocalPointSlotsInitialized = 0,
                  importLocalGlobalSlotsInitialized = 0,
                  importLocalCellSlotsInitialized = 0,
                  importPersistentNodeSlotsInitialized = 0,
                  partialImportRowsCopied = 0,
                  partialImportGlobalSlotsInitialized = 0,
                  partialImportCandidateMappingVisits = 0,
                  partialImportEdgeValidationRows = 0,
                  partialImportFaceValidationRows = 0,
                  partialImportProtectedEndpointVisits = 0,
                  partialImportAuthorityPointValidations = 0,
                  partialImportDuplicateKeyComparisons = 0,
                  partialImportResolutionNodeVisits = 0,
                  partialImportEdgeResolutionLookups = 0,
                  partialImportFaceResolutionLookups = 0,
                  partialImportEdgeDataReads = 0,
                  partialImportFaceDataCopies = 0,
                  partialImportFaceIncidenceRowsCopied = 0,
                  partialImportSortComparisons = 0,
                  surfaceWorkUnits = 0, tetCandidateWorkUnits = 0,
                  importWorkUnits = 0, extractionWorkUnits = 0,
                  totalWorkUnits = 0;
    bool runComplete = false;
};

struct BoundaryConformingV3MemoryMetrics {
    std::size_t surfaceArtifactBytes = 0, recoveryEnvelopeBytes = 0,
                tetCandidateBytes = 0, authoritativeWorkspaceBytes = 0,
                recoveryIndexBytes = 0, extractedWorkspaceBytes = 0,
                flatMeshBytes = 0, peakAccountedOverlapBytes = 0,
                finalRetainedOverlapBytes = 0;
    BoundaryConformingV3Phase peakOverlapStage =
        BoundaryConformingV3Phase::None;
    bool runComplete = false;
};

struct BoundaryConformingV3Audit {
    std::uint64_t prescribedExterior = 0, matchedExterior = 0,
                  orphanExterior = 0, unexpectedExterior = 0,
                  prescribedInterface = 0, matchedInterface = 0,
                  missingInterface = 0, unexpectedInterface = 0,
                  duplicateFacets = 0, nonManifoldFacets = 0,
                  invertedTets = 0, unusedNodes = 0;
    bool passed = false;
};

struct BoundaryConformingV3Request {
    ArtifactStamp expectedSurfaceSetStamp{};
    VolumeBoundaryProviderStamp expectedOwnershipProviderStamp{};
    RecoveryEnvelopeMode envelopeMode = RecoveryEnvelopeMode::DerivedFiniteBox;
    std::size_t maxProviderEnvelopeNodes = 256;
    std::size_t maxProviderEnvelopeFaces = 508;
    // Explicit extension for planar faces with holes/concavity. It changes
    // only surface candidate recovery; all exact surface/volume audits remain.
    bool usePlanarSurfaceRecovery = false;
    bnd::TetEdgeGeneralRecoveryPolicy edgeGeneralRecovery{};
};

struct TetRecoveryEnvelope {
    std::vector<Vec3> guardCoordinates;
    std::vector<std::array<std::size_t, 3>> outwardHullFaces;
};

#if defined(PAACS_TOPO_TESTING)
namespace boundary_conforming_detail {
enum class BoundaryConformingWallPhase : std::uint8_t {
    Surface,
    Artifact,
    Envelope,
    TetConstruction,
    AuthoritativeImport,
    RecoveryContext,
    Recovery,
    Extraction,
    AuditAndPublication,
    Count
};
struct BoundaryConformingWallTrace {
    bool enabled = false;
    std::array<std::uint64_t,
        static_cast<std::size_t>(BoundaryConformingWallPhase::Count)>
        nanoseconds{};
    TetCandidateMetrics tetCandidateMetrics{};
    void reset(bool value = false) noexcept {
        enabled = value;
        nanoseconds.fill(0);
        tetCandidateMetrics = {};
    }
    void clear() noexcept {
        nanoseconds.fill(0);
        tetCandidateMetrics = {};
    }
};
BoundaryConformingWallTrace& boundaryConformingWallTraceForTesting() noexcept;
bool validateTetEnvelopeConstructorFrameForTesting(
    const TetRecoveryEnvelope& envelope,
    std::span<const Vec3> constructionCoordinates) noexcept;
int derivedTetEnvelopeRungForTesting(
    std::span<const Vec3> constructionCoordinates,
    const TetRecoveryEnvelope& envelope) noexcept;
bool auditClosedTri3ContactsForTesting(
    const bnd::ClosedTri3Artifact& artifact,const mesh::Mesh& surface,
    ArtifactPreflightMetrics& metrics) noexcept;
BoundaryConformingV3Reason auditClosedTri3RegionShellsForTesting(
    const bnd::ClosedTri3Artifact& artifact) noexcept;
}  // namespace boundary_conforming_detail
#endif

struct ClosedArtifactBuildOutcome {
    BoundaryConformingV3Status status = BoundaryConformingV3Status::Failed;
    BoundaryConformingV3Failure failure{};
    ArtifactPreflightMetrics metrics{};
    std::unique_ptr<bnd::ClosedTri3Artifact> artifact;
};

struct RecoveryEnvelopeBuildOutcome {
    BoundaryConformingV3Status status = BoundaryConformingV3Status::Failed;
    BoundaryConformingV3Failure failure{};
    std::unique_ptr<TetRecoveryEnvelope> envelope;
};

struct BoundaryConformingV3ReportEnvelope {
    BoundaryConformingV3ReportEnvelope();
    ~BoundaryConformingV3ReportEnvelope();
    BoundaryConformingV3ReportEnvelope(
        const BoundaryConformingV3ReportEnvelope&) = delete;
    BoundaryConformingV3ReportEnvelope& operator=(
        const BoundaryConformingV3ReportEnvelope&) = delete;
    BoundaryConformingV3ReportEnvelope(
        BoundaryConformingV3ReportEnvelope&&) = delete;
    BoundaryConformingV3ReportEnvelope& operator=(
        BoundaryConformingV3ReportEnvelope&&) = delete;

    bool reportComplete = false;
    std::vector<BoundaryConformingV3ReportRow> processed;
    std::vector<BoundaryConformingV3WorkKey> unattempted;
};

struct BoundaryConformingV3Result {
    BoundaryConformingV3Result();
    ~BoundaryConformingV3Result();
    BoundaryConformingV3Result(const BoundaryConformingV3Result&) = delete;
    BoundaryConformingV3Result& operator=(const BoundaryConformingV3Result&) =
        delete;
    BoundaryConformingV3Result(BoundaryConformingV3Result&&) = delete;
    BoundaryConformingV3Result& operator=(BoundaryConformingV3Result&&) =
        delete;

    const char* routeId = kBoundaryConformingV3RouteId;
    std::uint32_t routeVersion = kBoundaryConformingV3RouteVersion;
    std::unique_ptr<TetAuthoritativeResult> finalAuthoritative;
    bnd::TetFacetSetState constraintState;
    std::vector<Index> globalByAuthoritativeNode;
    std::vector<Index> globalByAuthoritativeCell;
    ArtifactPreflightMetrics artifactPreflightMetrics{};
    BoundaryConformingV3ConstructionMetrics constructionMetrics{};
    TetRecoveryMetrics recoveryMetrics{};
    VolumeExtractionMetrics extractionMetrics{};
    BoundaryConformingV3MemoryMetrics memoryMetrics{};
    BoundaryConformingV3Audit boundaryAudit{};
};

struct BoundaryConformingV3Outcome {
    BoundaryConformingV3Outcome();
    ~BoundaryConformingV3Outcome();
    BoundaryConformingV3Outcome(const BoundaryConformingV3Outcome&) = delete;
    BoundaryConformingV3Outcome& operator=(const BoundaryConformingV3Outcome&) =
        delete;
    BoundaryConformingV3Outcome(BoundaryConformingV3Outcome&&) noexcept;
    BoundaryConformingV3Outcome& operator=(
        BoundaryConformingV3Outcome&&) noexcept;

    BoundaryConformingV3Status status = BoundaryConformingV3Status::Failed;
    BoundaryConformingV3Failure failure{};
    bool destinationUnchanged = true;
    bool reportComplete = false;
    std::unique_ptr<SurfaceSetOutcome> surfaceOutcome;
    std::unique_ptr<BoundaryConformingV3ReportEnvelope> report;
    std::unique_ptr<BoundaryConformingV3Result> retainedResult;

    bool published() const noexcept {
        return status == BoundaryConformingV3Status::Published;
    }
};

static_assert(std::is_standard_layout_v<VolumeBoundaryProviderStamp>);
static_assert(std::is_trivially_copyable_v<VolumeBoundaryProviderStamp>);
static_assert(std::is_standard_layout_v<BoundaryConformingV3Request>);
static_assert(std::is_trivially_copyable_v<BoundaryConformingV3Request>);
static_assert(std::is_standard_layout_v<BoundaryConformingV3WorkKey>);
static_assert(std::is_trivially_copyable_v<BoundaryConformingV3WorkKey>);
static_assert(std::is_standard_layout_v<BoundaryConformingV3Failure>);
static_assert(std::is_trivially_copyable_v<BoundaryConformingV3Failure>);
static_assert(std::is_standard_layout_v<BoundaryConformingV3ReportRow>);
static_assert(std::is_trivially_copyable_v<BoundaryConformingV3ReportRow>);
static_assert(std::is_nothrow_copy_constructible_v<
              BoundaryConformingV3ReportRow>);
static_assert(std::is_nothrow_move_constructible_v<
              BoundaryConformingV3ReportRow>);
static_assert(std::is_standard_layout_v<ArtifactPreflightMetrics>);
static_assert(std::is_trivially_copyable_v<ArtifactPreflightMetrics>);
static_assert(std::is_standard_layout_v<TetRecoveryMetrics>);
static_assert(std::is_trivially_copyable_v<TetRecoveryMetrics>);
static_assert(std::is_standard_layout_v<VolumeExtractionMetrics>);
static_assert(std::is_trivially_copyable_v<VolumeExtractionMetrics>);
static_assert(std::is_standard_layout_v<
              BoundaryConformingV3SurfaceConstructionMetrics>);
static_assert(std::is_trivially_copyable_v<
              BoundaryConformingV3SurfaceConstructionMetrics>);
static_assert(std::is_standard_layout_v<
              BoundaryConformingV3ConstructionMetrics>);
static_assert(std::is_trivially_copyable_v<
              BoundaryConformingV3ConstructionMetrics>);
static_assert(std::is_standard_layout_v<BoundaryConformingV3MemoryMetrics>);
static_assert(std::is_trivially_copyable_v<BoundaryConformingV3MemoryMetrics>);
static_assert(std::is_standard_layout_v<BoundaryConformingV3Audit>);
static_assert(std::is_trivially_copyable_v<BoundaryConformingV3Audit>);
static_assert(std::is_nothrow_move_constructible_v<
              std::vector<BoundaryConformingV3ReportRow>>);
static_assert(std::is_nothrow_move_constructible_v<
              std::vector<BoundaryConformingV3WorkKey>>);
static_assert(std::is_nothrow_move_constructible_v<
              std::unique_ptr<SurfaceSetOutcome>>);
static_assert(std::is_nothrow_move_constructible_v<
              std::unique_ptr<BoundaryConformingV3ReportEnvelope>>);
static_assert(std::is_nothrow_move_constructible_v<
              std::unique_ptr<BoundaryConformingV3Result>>);
static_assert(std::is_nothrow_move_assignable_v<SurfaceSetOutcome>);
static_assert(std::is_nothrow_move_constructible_v<BoundaryConformingV3Outcome>);
static_assert(!std::is_move_constructible_v<BoundaryConformingV3Result>);
static_assert(!std::is_move_constructible_v<
              BoundaryConformingV3ReportEnvelope>);
static_assert(noexcept(std::declval<mesh::Mesh&>().swap(
    std::declval<mesh::Mesh&>())));

ClosedArtifactBuildOutcome buildClosedTri3Artifact(
    const SurfaceSetResult& surface, const mesh::Mesh& publishedSurface,
    ArtifactStamp expectedSurfaceSetStamp,
    VolumeBoundaryProviderStamp expectedProviderStamp,
    const VolumeBoundaryOwnershipProvider& provider) noexcept;

RecoveryEnvelopeBuildOutcome buildTetRecoveryEnvelope(
    RecoveryEnvelopeMode mode, std::span<const Vec3> constructionCoordinates,
    VolumeBoundaryProviderStamp expectedProviderStamp,
    const VolumeBoundaryOwnershipProvider& provider,
    std::size_t maxProviderNodes = 256,
    std::size_t maxProviderFaces = 508) noexcept;

BoundaryConformingV3Outcome meshModelBoundaryConforming(
    const geo::Model& model, const sizing::MeshSpec& spec,
    const BoundaryConformingV3Request& request, mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& provider) noexcept;

}  // namespace paacs::algo
