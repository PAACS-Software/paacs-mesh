// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T3-A: bounded adaptive interior Tet4 refinement over boundary-conforming-v3.

#include "paacs/algo/BoundaryConformingMesher.hpp"
#include "paacs/mesh/Quality.hpp"
#include "paacs/topo/TetRefinement.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <type_traits>
#include <vector>

namespace paacs::algo {

inline constexpr const char* kT3ARouteId = "interior-refined-tet4-t3a";
inline constexpr std::uint32_t kT3ARouteVersion = 1;

enum class TetSliverPolicy : std::uint8_t { ReportOnly, BoundedInsertion };

struct TetRefinementPolicy {
    double maxEdgeToTarget = 1.5;
    double maxRadiusEdge = 3.0;
    double minimumCreatedTetVolume = 0.0;
    TetSliverPolicy sliverPolicy = TetSliverPolicy::ReportOnly;
    double sliverRadiusRatio = 0.0;
    double minimumSliverImprovement = 0.0;
    std::uint32_t maxSliverCandidateAttempts = 0;
    std::uint64_t maxInsertions = 64;
    std::uint64_t maxWorkUnits = 100000;
    std::uint64_t maxOutputNodes = 1000000;
    std::uint64_t maxOutputTets = 5000000;
    std::uint64_t maxReportRows = 325;
    std::uint64_t maxReportWitnesses = 100000;
    std::uint64_t maxReportTopologyHandles = 100064;
    std::uint32_t maxCandidateAttemptsPerViolation = 5;
};

struct T3ARequest {
    BoundaryConformingV3Request conformingRequest{};
    TetRefinementPolicy refinement{};
};

enum class T3AStatus : std::uint8_t {
    Published,
    Incomplete,
    Rejected,
    Unsupported,
    Failed
};

enum class T3APhase : std::uint8_t {
    None,
    ReportEnvelope,
    Preflight,
    BoundaryConformingV3,
    SeedTransfer,
    ContextBuild,
    InitialAudit,
    QueueBuild,
    SizingEvaluation,
    CandidateGeneration,
    Encroachment,
    CavityDiscovery,
    CandidateValidation,
    TopologyPreparation,
    CoordinatorPreparation,
    TopologyCommit,
    CoordinatorCommit,
    LocalAudit,
    FinalQualityAudit,
    MeshPreparation,
    Publication
};

enum class T3AReason : std::uint8_t {
    None,
    InvalidRequest,
    InvalidRefinementPolicy,
    InvalidSizingEvaluation,
    SizingFieldThrew,
    SizingFieldChanged,
    StaleConformingResult,
    InvalidConformingResult,
    MultipleRegionsOrInterfaces,
    UnsupportedElementOrder,
    ProtectedIndexRangeUnrepresentable,
    InvalidProtectedFeature,
    SourceRejected,
    SourceRefused,
    SourceFailed,
    BoundaryRefinementRequired,
    InsertionBudgetReached,
    WorkBudgetReached,
    OutputBudgetReached,
    ReportBudgetReached,
    NoAdmissibleCandidate,
    Nonprogress,
    CreatedTetVolumeBlocked,
    SliverPolicyExhausted,
    NonfiniteCandidate,
    StaleRecoveryContext,
    InvalidCavity,
    InvalidCandidate,
    TopologyPreparationFailed,
    CoordinatorPreparationFailed,
    ProtectedIdentityChanged,
    BoundaryBijectionFailed,
    QualityAuditFailed,
    MeshPreparationFailed,
    AllocationFailed,
    InternalError
};

enum class RefinementVertexKind : std::uint8_t {
    ExistingGlobal,
    InteriorInsertion
};

struct RefinementVertexKey {
    RefinementVertexKind kind = RefinementVertexKind::ExistingGlobal;
    std::uint64_t value = 0;
    friend bool operator==(const RefinementVertexKey&,
                           const RefinementVertexKey&) = default;
    friend bool operator<(const RefinementVertexKey& a,
                          const RefinementVertexKey& b) noexcept {
        return a.kind < b.kind || (a.kind == b.kind && a.value < b.value);
    }
};

struct T3AEdgeKey {
    std::array<RefinementVertexKey,2> vertices{};
    friend bool operator==(const T3AEdgeKey&,const T3AEdgeKey&) = default;
    friend bool operator<(const T3AEdgeKey& a,const T3AEdgeKey& b) noexcept {
        return a.vertices < b.vertices;
    }
};
struct T3AFaceKey {
    std::array<RefinementVertexKey,3> vertices{};
    friend bool operator==(const T3AFaceKey&,const T3AFaceKey&) = default;
    friend bool operator<(const T3AFaceKey& a,const T3AFaceKey& b) noexcept {
        return a.vertices < b.vertices;
    }
};
struct T3ATetKey {
    std::array<RefinementVertexKey,4> vertices{};
    friend bool operator==(const T3ATetKey&,const T3ATetKey&) = default;
    friend bool operator<(const T3ATetKey& a,const T3ATetKey& b) noexcept {
        return a.vertices < b.vertices;
    }
};

enum class TetViolationKind : std::uint8_t { Size, RadiusEdge, Sliver };

struct TetViolationKey {
    TetViolationKind kind = TetViolationKind::Size;
    T3ATetKey tet{};
    bool hasSizeEdge = false;
    T3AEdgeKey sizeEdge{};
    friend bool operator==(const TetViolationKey&,
                           const TetViolationKey&) = default;
};

struct T3ASizingSample { Vec3 point{}; double value = 0.0; };

struct TetViolationRow {
    TetViolationKey key{};
    int regionTag = -1;
    bool hasSizeSamples = false;
    std::array<T3ASizingSample,3> sizeSamples{};
    double normalizedSeverity = 0.0;
    std::uint64_t severityBits = 0;
    std::uint64_t generation = 0;
};

enum class T3ACandidateKind : std::uint8_t {
    TetCircumcenter,
    FaceOffcenter
};
enum class T3AProtectedFeatureKind : std::uint8_t { Edge, Facet };

struct T3AProtectedFeatureWitness {
    T3AProtectedFeatureKind kind = T3AProtectedFeatureKind::Edge;
    bnd::VolumeArtifactEdgeKey edge{};
    bnd::VolumeFacetKey facet{};
};

enum class T3AAttemptDisposition : std::uint8_t {
    Committed,
    FormulaUnavailable,
    DuplicateProposal,
    Nonfinite,
    Encroached,
    OutsideDomain,
    InvalidCavity,
    InvalidCandidate,
    CreatedTetVolumeBlocked,
    NoProgress,
    SliverImprovementFailed,
    TerminalFailure
};

struct T3AAttemptRow {
    std::uint64_t attemptOrdinal = 0;
    TetViolationKey violation{};
    std::uint64_t violationGeneration = 0;
    T3ACandidateKind candidateKind = T3ACandidateKind::TetCircumcenter;
    bool hasFace = false;
    T3AFaceKey face{};
    Vec3 proposedCoordinate{};
    T3AAttemptDisposition disposition = T3AAttemptDisposition::TerminalFailure;
    T3AReason reason = T3AReason::None;
    std::uint64_t witnessOffset = 0;
    std::uint64_t witnessCount = 0;
    std::uint64_t cavityDigest = 0;
    std::uint64_t workBefore = 0;
    std::uint64_t workAfter = 0;
};

inline constexpr std::uint64_t kCurrentBlockedInsertion =
    (std::numeric_limits<std::uint64_t>::max)();

struct T3ABlockedViolationRow {
    TetViolationKey violation{};
    std::uint64_t violationGeneration = 0;
    std::uint64_t attemptOrdinal = 0;
    std::uint64_t witnessOffset = 0;
    std::uint64_t witnessCount = 0;
    std::uint64_t invalidatedByInsertionOrdinal = kCurrentBlockedInsertion;
};

enum class T3ATopologyKind : std::uint8_t { Node, Cell, Edge, Face };
enum class T3ATopologyTransition : std::uint8_t {
    Retained,
    Created,
    Retired
};

struct T3ATopologyHandleRow {
    T3ATopologyKind kind = T3ATopologyKind::Node;
    T3ATopologyTransition transition = T3ATopologyTransition::Retained;
    topo::NodeHandle predictedNode{},actualNode{};
    topo::CellHandle predictedCell{},actualCell{};
    topo::EdgeHandle predictedEdge{},actualEdge{};
    topo::FaceHandle predictedFace{},actualFace{};
};

struct T3ACommitRow {
    std::uint64_t insertionOrdinal = 0;
    TetViolationKey violation{};
    std::uint32_t stagedLocalId = 0;
    topo::NodeHandle predictedNode{},actualNode{};
    std::uint64_t topologyOffset = 0,topologyCount = 0;
    std::uint64_t postCommitRevision = 0;
};

struct T3AFailure {
    T3AReason reason = T3AReason::None;
    T3APhase phase = T3APhase::None;
    bool hasViolation = false;
    TetViolationKey violation{};
    bool hasProtectedWitness = false;
    T3AProtectedFeatureWitness protectedWitness{};
    std::uint64_t insertionOrdinal = 0;
    std::uint64_t completedWorkUnits = 0;
};

struct TetRefinementMetrics {
    std::uint64_t initialCellsScanned = 0,finalCellsAudited = 0,
        sizingEvaluations = 0,violationEvaluations = 0,queuePushes = 0,
        queuePops = 0,staleQueueRows = 0,queueHighWater = 0,
        candidateAttempts = 0,rejectedCandidates = 0,
        protectedIndexVisits = 0,protectedEdgeTests = 0,
        protectedFacetTests = 0,locationIndexVisits = 0,
        locationIndexUpdateVisits = 0,locationIndexUpdateWorkCharged = 0,
        coordinatorRowsPrepared = 0,coordinatorTailWorkCharged = 0,
        locationTetTests = 0,cavityCellsVisited = 0,cavityFaceVisits = 0,
        cavityHighWater = 0,orient3dTests = 0,insphereTests = 0,
        survivorIndexVisits = 0,candidateSurvivorTests = 0,
        preparedTransactions = 0,committedInsertions = 0,
        blockedRowsRecorded = 0,blockedRowsInvalidated = 0,
        blockedRowsCurrent = 0,blockedRowsCurrentHighWater = 0,
        contextDeltaUpdates = 0,contextRebuilds = 0,
        locallyReevaluatedCells = 0,perInsertionGlobalTetScans = 0,
        phase3Allocations = 0,publicationTailAllocations = 0,
        workUnitsConsumed = 0,transientAccountedHighWaterBytes = 0,
        authoritativeAccountedBytes = 0,publicationAccountedBytes = 0;
    bool runComplete = false,qualityComplete = false;
};

struct TetRefinementQualityReport {
    std::uint64_t nodeCount = 0,tetCount = 0,exactPositive = 0,
        inverted = 0,nonfinite = 0,inheritedBelowCreatedFloor = 0,
        createdBelowFloor = 0,sizeViolations = 0,sliversAttempted = 0,
        sliversRemaining = 0;
    double minimumSignedVolume = 0.0,maximumRadiusEdge = 0.0,
        p95RadiusEdge = 0.0,p99RadiusEdge = 0.0,
        minimumRadiusRatio = 0.0,p1RadiusRatio = 0.0,p5RadiusRatio = 0.0,
        worstNormalizedSize = 0.0;
    std::array<std::uint64_t,10> radiusEdgeHistogram{};
    std::array<std::uint64_t,10> radiusRatioHistogram{};
    TetSliverPolicy sliverPolicy = TetSliverPolicy::ReportOnly;
    bool sliverGoalComplete = false;
};

struct BoundaryConformingV3Receipt {
    BoundaryConformingV3Receipt() = default;
    BoundaryConformingV3Receipt(const BoundaryConformingV3Receipt&) = delete;
    BoundaryConformingV3Receipt& operator=(const BoundaryConformingV3Receipt&) = delete;
    BoundaryConformingV3Receipt(BoundaryConformingV3Receipt&&) = delete;
    BoundaryConformingV3Receipt& operator=(BoundaryConformingV3Receipt&&) = delete;
    BoundaryConformingV3Status status = BoundaryConformingV3Status::Failed;
    BoundaryConformingV3Failure failure{};
    bool destinationUnchanged = true,reportComplete = false;
    const char* routeId = kBoundaryConformingV3RouteId;
    std::uint32_t routeVersion = kBoundaryConformingV3RouteVersion;
    std::unique_ptr<SurfaceSetOutcome> surfaceOutcome;
    std::unique_ptr<BoundaryConformingV3ReportEnvelope> report;
    ArtifactPreflightMetrics artifactPreflightMetrics{};
    TetRecoveryMetrics recoveryMetrics{};
    VolumeExtractionMetrics extractionMetrics{};
    BoundaryConformingV3MemoryMetrics memoryMetrics{};
    BoundaryConformingV3Audit boundaryAudit{};
    bool authoritativeResultConsumed = false;
};

struct T3ABoundaryTriRecord {
    bnd::VolumeFacetKey facetKey{};
    std::array<Index,3> directedGlobalNodes{
        kInvalidIndex,kInvalidIndex,kInvalidIndex};
    mesh::Classification classification{};
    topo::FaceHandle authoritativeFace{};
};

struct T3AReport {
    T3AReport() = default;
    T3AReport(const T3AReport&) = delete;
    T3AReport& operator=(const T3AReport&) = delete;
    T3AReport(T3AReport&&) = delete;
    T3AReport& operator=(T3AReport&&) = delete;
    bool envelopePrepared = false,reportComplete = false;
    std::vector<T3AAttemptRow> attemptStorage;
    std::vector<T3AProtectedFeatureWitness> protectedWitnessStorage;
    std::vector<T3ABlockedViolationRow> blockedStorage;
    std::vector<TetViolationKey> unattemptedStorage;
    std::vector<T3ACommitRow> committedPrefixStorage;
    std::vector<T3ATopologyHandleRow> topologyHandleStorage;
    std::uint64_t attemptsUsed = 0,protectedWitnessesUsed = 0,blockedUsed = 0,
        unattemptedUsed = 0,committedPrefixUsed = 0,topologyHandlesUsed = 0;
    bool hasTerminalViolation = false;
    TetViolationKey terminalViolation{};

    std::span<const T3AAttemptRow> attempts() const noexcept {
        return {attemptStorage.data(),static_cast<std::size_t>(attemptsUsed)};
    }
    std::span<const T3ABlockedViolationRow> blocked() const noexcept {
        return {blockedStorage.data(),static_cast<std::size_t>(blockedUsed)};
    }
};

struct T3AResult {
    T3AResult() = default;
    T3AResult(const T3AResult&) = delete;
    T3AResult& operator=(const T3AResult&) = delete;
    T3AResult(T3AResult&&) = delete;
    T3AResult& operator=(T3AResult&&) = delete;
    const char* routeId = kT3ARouteId;
    std::uint32_t routeVersion = kT3ARouteVersion;
    std::unique_ptr<TetAuthoritativeResult> authoritative;
    bnd::TetFacetSetState constraintState;
    std::vector<Index> existingGlobalByAuthoritativeNodeSlot;
    std::vector<Index> publishedGlobalByAuthoritativeNodeSlot;
    std::vector<Index> publishedGlobalByAuthoritativeCellSlot;
    std::vector<T3ABoundaryTriRecord> immutableBoundary;
    TetRefinementMetrics metrics{};
    TetRefinementQualityReport quality{};
    BoundaryConformingV3Audit boundaryAudit{};
    bool publicationReady = false;
};

struct T3AOutcome {
    T3AOutcome() = default;
    T3AOutcome(const T3AOutcome&) = delete;
    T3AOutcome& operator=(const T3AOutcome&) = delete;
    T3AOutcome(T3AOutcome&&) noexcept = default;
    T3AOutcome& operator=(T3AOutcome&&) noexcept = default;
    T3AStatus status = T3AStatus::Failed;
    T3AFailure failure{};
    bool destinationUnchanged = true,reportComplete = false;
    std::unique_ptr<BoundaryConformingV3Receipt> sourceReceipt;
    std::unique_ptr<T3AReport> report;
    std::unique_ptr<T3AResult> retainedResult;
    bool published() const noexcept { return status == T3AStatus::Published; }
};

static_assert(std::is_standard_layout_v<TetRefinementPolicy>);
static_assert(std::is_trivially_copyable_v<TetRefinementPolicy>);
static_assert(std::is_standard_layout_v<RefinementVertexKey>);
static_assert(std::is_trivially_copyable_v<RefinementVertexKey>);
static_assert(std::is_standard_layout_v<T3AEdgeKey>);
static_assert(std::is_trivially_copyable_v<T3AEdgeKey>);
static_assert(std::is_standard_layout_v<T3AFaceKey>);
static_assert(std::is_trivially_copyable_v<T3AFaceKey>);
static_assert(std::is_standard_layout_v<T3ATetKey>);
static_assert(std::is_trivially_copyable_v<T3ATetKey>);
static_assert(std::is_standard_layout_v<TetViolationKey>);
static_assert(std::is_trivially_copyable_v<TetViolationKey>);
static_assert(std::is_standard_layout_v<T3ASizingSample>);
static_assert(std::is_trivially_copyable_v<T3ASizingSample>);
static_assert(std::is_standard_layout_v<TetViolationRow>);
static_assert(std::is_trivially_copyable_v<TetViolationRow>);
static_assert(std::is_standard_layout_v<T3AProtectedFeatureWitness>);
static_assert(std::is_trivially_copyable_v<T3AProtectedFeatureWitness>);
static_assert(std::is_standard_layout_v<T3AFailure>);
static_assert(std::is_trivially_copyable_v<T3AFailure>);
static_assert(std::is_standard_layout_v<T3AAttemptRow>);
static_assert(std::is_trivially_copyable_v<T3AAttemptRow>);
static_assert(std::is_standard_layout_v<T3ABlockedViolationRow>);
static_assert(std::is_trivially_copyable_v<T3ABlockedViolationRow>);
static_assert(std::is_standard_layout_v<T3ATopologyHandleRow>);
static_assert(std::is_trivially_copyable_v<T3ATopologyHandleRow>);
static_assert(std::is_standard_layout_v<T3ACommitRow>);
static_assert(std::is_trivially_copyable_v<T3ACommitRow>);
static_assert(std::is_standard_layout_v<T3ABoundaryTriRecord>);
static_assert(std::is_trivially_copyable_v<T3ABoundaryTriRecord>);
static_assert(std::is_standard_layout_v<TetRefinementMetrics>);
static_assert(std::is_trivially_copyable_v<TetRefinementMetrics>);
static_assert(std::is_standard_layout_v<TetRefinementQualityReport>);
static_assert(std::is_trivially_copyable_v<TetRefinementQualityReport>);
static_assert(std::is_nothrow_move_constructible_v<T3AOutcome>);
static_assert(std::is_nothrow_move_assignable_v<T3AOutcome>);
static_assert(!std::is_move_constructible_v<T3AResult>);
static_assert(!std::is_move_constructible_v<T3AReport>);
static_assert(!std::is_move_constructible_v<BoundaryConformingV3Receipt>);

T3AOutcome meshModelRefinedTet4(
    const geo::Model& model,const sizing::MeshSpec& spec,
    const T3ARequest& request,mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& ownershipProvider) noexcept;

#if defined(PAACS_TOPO_TESTING)
namespace t3a_detail {
inline bool& refinementSizingPhaseForTesting() noexcept {
    static thread_local bool value=false;
    return value;
}
inline bool& injectPostCommitAuditFailureForTesting() noexcept {
    static bool value=false;
    return value;
}
inline int& reportFaultCountdownForTesting() noexcept {
    static int value=0;
    return value;
}
inline int& blockerReportBudgetCountdownForTesting() noexcept {
    static thread_local int value=0;
    return value;
}
inline int& blockerWorkBudgetCountdownForTesting() noexcept {
    static thread_local int value=0;
    return value;
}
inline int& blockerInvalidationFaultCountdownForTesting() noexcept {
    static thread_local int value=0;
    return value;
}
enum class T3ABlockerDigestCheckpoint : std::uint8_t {
    BlockTransitionBefore,
    BlockTransitionTerminalAfter,
    InvalidationBeforeTransaction,
    InvalidationFailureAfterAbandon,
    Final
};
struct T3ABlockerStateDigest {
    T3ABlockerDigestCheckpoint checkpoint =
        T3ABlockerDigestCheckpoint::BlockTransitionBefore;
    std::uint64_t topologySemantic = 0,topologyStrict = 0,
        topologyRevision = 0,poolState = 0,queueStructural = 0,
        queueLogical = 0,heapStructural = 0,heapLogical = 0,
        positionStructural = 0,history = 0,reportState = 0,
        blockerMetrics = 0,structural = 0,logical = 0;
    std::uint64_t attemptsUsed = 0,protectedWitnessesUsed = 0,
        blockedUsed = 0,unattemptedUsed = 0,committedPrefixUsed = 0,
        currentBlocked = 0;
};
struct T3AObservedViolationGeneration {
    TetViolationKey key{};
    std::uint64_t generation = 0;
};
struct T3ABlockerDigestTrace {
    bool enabled = false;
    std::uint32_t used = 0;
    std::uint32_t checkpointMask = (std::numeric_limits<std::uint32_t>::max)();
    std::array<T3ABlockerStateDigest,16> rows{};
    std::uint32_t finalLiveUsed = 0;
    bool finalLiveComplete = false;
    std::array<T3AObservedViolationGeneration,4096> finalLive{};
    void reset(bool enable=true,
               std::uint32_t mask=(std::numeric_limits<std::uint32_t>::max)())
        noexcept {
        enabled=enable;
        used=0;
        checkpointMask=mask;
        rows={};
        finalLiveUsed=0;
        finalLiveComplete=false;
        finalLive={};
    }
};
inline T3ABlockerDigestTrace& blockerDigestTraceForTesting() noexcept {
    static thread_local T3ABlockerDigestTrace value;
    return value;
}
enum class T3AWallPhase : std::uint8_t {
    Source,
    Setup,
    ProtectedQuery,
    Location,
    Cavity,
    CandidateAudit,
    TopologyPreparation,
    CoordinatorPreparation,
    LocalAudit,
    FinalAudit,
    Count
};
struct T3AWallTimingTrace {
    bool enabled=false;
    std::array<std::uint64_t,
        static_cast<std::size_t>(T3AWallPhase::Count)> nanoseconds{};
    std::array<std::uint64_t,
        static_cast<std::size_t>(T3AWallPhase::Count)> calls{};
    void reset(bool enable=true) noexcept {
        enabled=enable;
        nanoseconds={};
        calls={};
    }
};
inline T3AWallTimingTrace& wallTimingTraceForTesting() noexcept {
    static thread_local T3AWallTimingTrace value;
    return value;
}
bool certifiedProtectedEdgeBoundsForTesting(
    const std::array<Vec3,2>& edge,Vec3& lo,Vec3& hi) noexcept;
bool certifiedProtectedFacetBoundsForTesting(
    const std::array<Vec3,3>& facet,Vec3& lo,Vec3& hi) noexcept;
bool evaluateTetViolationForTesting(
    const std::array<Vec3,4>& points,const std::array<Index,4>& globalNodes,
    int regionTag,const sizing::MeshSpec& spec,
    const TetRefinementPolicy& policy,TetViolationRow& row) noexcept;
bool exactLocationForTesting(const topo::TopologyWorkspace& workspace,
                             Vec3 point,topo::CellHandle& seed) noexcept;
bool validateRefinementCavityForTesting(
    const topo::TopologyWorkspace& workspace,
    std::span<const topo::CellHandle> removed,Vec3 stagedPoint) noexcept;
bool dynamicLocationScalingForTesting(
    std::size_t insertionCount,std::uint64_t& updateVisits,
    std::uint64_t& queryVisits) noexcept;
bool coordinatorDeltaScalingForTesting(
    std::size_t currentCellCount,std::uint64_t& affectedRows,
    std::uint64_t& workUnits,std::uint64_t& accountedRemovedBytes,
    std::uint64_t& accountedAddedBytes) noexcept;
bool blockedSchedulerReuseForTesting() noexcept;
}
#endif

} // namespace paacs::algo
