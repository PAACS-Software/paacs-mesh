// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T2-B artifact identities and persistent multi-facet recovery state.

#include "paacs/algo/SurfaceMesherAdapter.hpp"
#include "paacs/bnd/TetEdgeRecovery.hpp"
#include "paacs/bnd/TetFacetRecovery.hpp"
#include "paacs/bnd/TetRecoveryIndex.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace paacs::bnd {

struct VolumeFacetKey {
    algo::FacePersistentId facePersistentId = 0;
    std::size_t normalizedLocalCellOrdinal = 0;
    friend bool operator==(const VolumeFacetKey&, const VolumeFacetKey&) =
        default;
    friend bool operator<(const VolumeFacetKey& a,
                          const VolumeFacetKey& b) noexcept {
        return a.facePersistentId < b.facePersistentId ||
               (a.facePersistentId == b.facePersistentId &&
                a.normalizedLocalCellOrdinal <
                    b.normalizedLocalCellOrdinal);
    }
};

struct VolumeRegionUse {
    int regionTag = -1;
    int directedSense = 0;
    friend bool operator==(const VolumeRegionUse&,
                           const VolumeRegionUse&) = default;
    friend bool operator<(const VolumeRegionUse& a,
                          const VolumeRegionUse& b) noexcept {
        return a.regionTag < b.regionTag ||
               (a.regionTag == b.regionTag &&
                a.directedSense < b.directedSense);
    }
};

struct VolumeArtifactEdge {
    VolumeArtifactEdgeKey key{};
    mesh::Classification classification{};
    std::vector<VolumeFacetKey> incidentFacets;
};

struct VolumeFacet {
    VolumeFacetKey key{};
    std::array<Index, 3> directedGlobalNodes{
        kInvalidIndex, kInvalidIndex, kInvalidIndex};
    mesh::Classification classification{};
    std::vector<VolumeRegionUse> owners;
};

struct ClosedTri3Artifact {
    algo::ArtifactStamp surfaceSetStamp{};
    std::vector<VolumeFacet> facets;
    std::vector<VolumeArtifactEdge> edges;
    std::vector<int> regionTags;
};

struct PersistentTetEdge {
    VolumeArtifactEdgeKey key{};
    mesh::Classification classification{};
    bool required = true;
    topo::EdgeHandle edge{};
};

struct PersistentTetFacetBinding {
    VolumeFacetKey key{};
    std::array<topo::NodeHandle, 3> canonicalVertices{};
    mesh::Classification classification{};
    topo::FaceHandle recoveredFace{};
};

struct TetFacetSetState {
    topo::WorkspaceStamp workspaceStamp{};
    algo::ArtifactStamp artifactStamp{};
    std::map<VolumeArtifactEdgeKey, PersistentTetEdge> edgeBindings;
    std::map<VolumeFacetKey, PersistentTetFacetBinding> bindings;
    std::vector<VolumeArtifactEdgeKey> canonicalEdgeSchedule;
    std::vector<VolumeFacetKey> canonicalSchedule;

    void swap(TetFacetSetState& other) noexcept {
        using std::swap;
        swap(workspaceStamp, other.workspaceStamp);
        swap(artifactStamp, other.artifactStamp);
        edgeBindings.swap(other.edgeBindings);
        bindings.swap(other.bindings);
        canonicalEdgeSchedule.swap(other.canonicalEdgeSchedule);
        canonicalSchedule.swap(other.canonicalSchedule);
    }
};

struct TetFacetIndexedRequest {
    topo::WorkspaceStamp stamp{};
    algo::ArtifactStamp artifactStamp{};
    VolumeFacetKey key{};
    std::array<topo::NodeHandle, 3> vertices{};
    mesh::Classification classification{};
    mesh::Classification region{};
};

struct TetFacetIndexedOutcome {
    TetFacetRecoveryOutcome recovery{};
    std::vector<topo::CellHandle> retiredCells;
    std::vector<topo::CellHandle> createdCells;
    bool contextUpdated = false;
};

enum class TetFacetSetRecoveryStatus : std::uint8_t {
    Complete, Rejected, Refused, Failed
};

enum class TetFacetSetWorkKind : std::uint8_t { Edge, Facet };

struct TetFacetSetRecoveryRow {
    TetFacetSetWorkKind kind = TetFacetSetWorkKind::Edge;
    VolumeArtifactEdgeKey edgeKey{};
    VolumeFacetKey facetKey{};
    TetEdgeRecoveryOutcome edgeOutcome{};
    TetFacetRecoveryOutcome facetOutcome{};
};

enum class TetFacetSetRecoveryReason : std::uint8_t {
    None,
    RecoveryConflictCycle,
    InternalError,
    AllocationFailed,
    ConflictCycleSinkCapacityReached,
    MixedPresentCapacityReached,
    FrontierSinkCapacityReached
};

struct TetFacetSetSchedulerMetrics {
    std::uint64_t routeAssessments = 0;
    std::uint64_t dependencyEdges = 0;
    std::uint64_t dependencyInvalidations = 0;
    std::uint64_t selectedRebuilds = 0;
    std::uint64_t batchPreflightAssessments = 0;
    std::uint64_t alreadyPresentBatchEntities = 0;
    std::uint64_t protectedConflicts = 0;
    std::uint64_t conflictSccVisits = 0;
    std::uint64_t fullStageRescans = 0;
};

struct TetFacetSetRecoveryOutcome {
    TetFacetSetRecoveryStatus status = TetFacetSetRecoveryStatus::Failed;
    TetFacetSetRecoveryReason reason = TetFacetSetRecoveryReason::None;
    bool reportComplete = false;
    std::vector<TetFacetSetRecoveryRow> processed;
    std::vector<VolumeArtifactEdgeKey> unattemptedEdges;
    std::vector<VolumeFacetKey> unattemptedFacets;
    std::vector<VolumeArtifactEdgeKey> conflictEdgeVertices;
    std::vector<std::array<VolumeArtifactEdgeKey,2>> conflictEdgeDependencies;
    std::vector<VolumeFacetKey> conflictVertices;
    std::vector<std::array<VolumeFacetKey,2>> conflictEdges;
    TetFacetSetSchedulerMetrics schedulerMetrics{};
};

// Additive v4-only immutable-frontier observation.  Callers provide every
// byte of storage; the accepted v3 call supplies no sink and is unchanged.
struct TetFacetSetFrontierRow {
    TetFacetSetWorkKind kind = TetFacetSetWorkKind::Edge;
    VolumeArtifactEdgeKey edgeKey{};
    VolumeFacetKey facetKey{};
    topo::WorkspaceStamp workspaceStamp{};
    algo::ArtifactStamp artifactStamp{};
    TetEdgeRecoveryOutcome edgeOutcome{};
    TetFacetRecoveryOutcome facetOutcome{};
    std::size_t supportOffset = 0;
    std::size_t supportCount = 0;
    std::size_t dependencyOffset = 0;
    std::size_t dependencyCount = 0;
};

struct TetFacetSetFrontierSink {
    std::span<TetFacetSetFrontierRow> rows;
    std::span<topo::CellHandle> supportCells;
    std::span<std::size_t> dependencies;
    std::size_t rowsUsed = 0;
    std::size_t supportCellsUsed = 0;
    std::size_t dependenciesUsed = 0;
    TetFacetGeometryCounters assessmentCounters{};
    bool complete = true;

    void reset() noexcept {
        rowsUsed=0;
        supportCellsUsed=0;
        dependenciesUsed=0;
        assessmentCounters={};
        complete=true;
    }
};

// Additive v4-only transport for the exact canonical facet conflict SCC.
// Callers own all storage. The accepted v3 route passes no sink and retains
// its existing outcome and allocation behavior.
struct TetFacetConflictCycleSink {
    std::span<VolumeFacetKey> members;
    std::span<std::array<VolumeFacetKey,2>> dependencies;
    topo::WorkspaceStamp workspaceStamp{};
    algo::ArtifactStamp artifactStamp{};
    std::size_t membersUsed = 0;
    std::size_t dependenciesUsed = 0;
    std::size_t membersRequired = 0;
    std::size_t dependenciesRequired = 0;
    bool captured = false;
    bool complete = true;

    void reset() noexcept {
        workspaceStamp = {};
        artifactStamp = {};
        membersUsed = 0;
        dependenciesUsed = 0;
        membersRequired = 0;
        dependenciesRequired = 0;
        captured = false;
        complete = true;
    }
};

// Additive v4-private control and measurement for mixed exact-present
// scheduling.  Supplying this sink enables canonical-prefix batching;
// accepted v3 supplies null and executes the original scheduler unchanged.
// The byte counter is a portable logical payload/capacity count.  Exact
// allocator calls and requested bytes are qualified by the dedicated probe.
struct TetFacetSetMixedPresentMetrics {
    static constexpr std::size_t statusCount = 6;
    // These arrays count assessment events, including deterministic
    // re-assessments after an indexed delta. They are not unique entities.
    std::array<std::uint64_t,statusCount>
        edgeAssessmentEventStatusCounts{};
    std::array<std::uint64_t,statusCount>
        facetAssessmentEventStatusCounts{};
    std::uint64_t presentBindingBatches = 0;
    std::uint64_t presentBindingEntities = 0;
    std::uint64_t stateCopies = 0;
    std::uint64_t copiedRows = 0;
    std::uint64_t copiedLogicalBytes = 0;
    // Scheduler/batch reverse-binding checks only. Complete low-level
    // binding validation work is retained separately below.
    std::uint64_t schedulerBindingAuditRows = 0;
    std::uint64_t pendingRowsVisited = 0;
    std::uint64_t dependencyPairsTested = 0;
    std::uint64_t exactPresentSelectedRebuilds = 0;
    std::uint64_t exactPresentDependencyTests = 0;
    std::uint64_t nontrivialRecoveries = 0;
    // A deterministic model of the prefix-batch scratch containers only;
    // not an exhaustive allocator-traffic count. The allocation probe is
    // authoritative for complete route calls/bytes.
    std::uint64_t modeledBatchPreparationAllocationRequests = 0;
    std::uint64_t modeledBatchPreparationLogicalBytes = 0;
    TetFacetGeometryCounters assessmentCounters{};
    TetFacetGeometryCounters recoveryCounters{};
    TetRecoveryBindingWorkCounters assessmentBindingWork{};
    TetRecoveryBindingWorkCounters recoveryBindingWork{};
    bool complete = true;

    void reset() noexcept { *this = {}; }
    friend bool operator==(const TetFacetSetMixedPresentMetrics&,
                           const TetFacetSetMixedPresentMetrics&) = default;
};

struct TetFacetSetMixedPresentBatchSink {
    TetFacetSetMixedPresentMetrics metrics{};
#if defined(PAACS_TOPO_TESTING)
    // Qualification-only canonical-prefix oracle control. Production builds
    // have one fixed maximal-prefix policy.
    std::size_t maximumBatchEntitiesForTesting =
        (std::numeric_limits<std::size_t>::max)();
#endif
    void reset() noexcept { metrics.reset(); }
};

static_assert(static_cast<std::size_t>(TetEdgeRecoveryStatus::Failed)+1==
              TetFacetSetMixedPresentMetrics::statusCount);
static_assert(static_cast<std::size_t>(TetFacetRecoveryStatus::Failed)+1==
              TetFacetSetMixedPresentMetrics::statusCount);

static_assert(noexcept(std::declval<TetFacetSetState&>().swap(
    std::declval<TetFacetSetState&>())));

TetFacetIndexedOutcome recoverTetFacetIndexed(
    const TetFacetIndexedRequest& request, TetFacetSetState& state,
    topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    TetRecoveryBindingWorkCounters* bindingWork = nullptr) noexcept;

TetFacetRecoveryOutcome assessTetFacetIndexed(
    const TetFacetIndexedRequest& request, const TetFacetSetState& state,
    const topo::TopologyWorkspace& workspace,
    TetRecoveryContext& context,
    TetRecoveryBindingWorkCounters* bindingWork = nullptr) noexcept;

// The set scheduler additionally consumes the exact selected-cell witness.
// Like assessTetFacetIndexed(), this never prepares or mutates topology.
TetFacetIndexedOutcome assessTetFacetIndexedDetailed(
    const TetFacetIndexedRequest& request, const TetFacetSetState& state,
    const topo::TopologyWorkspace& workspace,
    TetRecoveryContext& context,
    TetRecoveryBindingWorkCounters* bindingWork = nullptr) noexcept;

TetFacetSetRecoveryOutcome recoverTetFacetSetIndexed(
    const ClosedTri3Artifact& artifact,
    const std::map<Index, topo::NodeHandle>& localByGlobal,
    mesh::Classification temporaryRegion, TetFacetSetState& state,
    topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    TetFacetSetFrontierSink* frontier = nullptr,
    TetFacetConflictCycleSink* conflictCycle = nullptr,
    TetFacetSetMixedPresentBatchSink* mixedPresent = nullptr) noexcept;

// The optional policy is one cumulative general-cavity work budget for this
// scheduler invocation. Disabled policy preserves the original route.
TetFacetSetRecoveryOutcome recoverTetFacetSetIndexed(
    const ClosedTri3Artifact& artifact,
    const std::map<Index, topo::NodeHandle>& localByGlobal,
    mesh::Classification temporaryRegion, TetFacetSetState& state,
    topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    const TetEdgeGeneralRecoveryPolicy& generalRecovery,
    TetFacetSetFrontierSink* frontier = nullptr,
    TetFacetConflictCycleSink* conflictCycle = nullptr,
    TetFacetSetMixedPresentBatchSink* mixedPresent = nullptr) noexcept;

#if defined(PAACS_TOPO_TESTING)
// True only across the verified allocation-free report/state/reverse-index
// publication tail of one mixed exact-present prefix batch.
bool mixedPresentCommitZoneActiveForTesting() noexcept;

// Testing-only bounded oracle for the proposed T2-C edge-stage closure.  It
// mutates only the already-disposable failed source workspace, stops at the
// edge-stage fixed point, and never enters facet routing.  The plain rows are
// retained after that workspace is destroyed so an independent test oracle
// can exhaust every blocker/blocker and blocker/transaction pairing.
enum class T2CEdgeLineageStatusForTesting : std::uint8_t {
    None,
    EdgeFixedPoint,
    NotApplicable,
    PermanentOutcome,
    BlockerOverlap,
    StaleSupport,
    NoProgress,
    CapacityReached,
    WorkBudgetReached,
    MemoryBudgetReached,
    AllocationFailed,
    InternalError
};

struct T2CEdgeLineageLimitsForTesting {
    std::size_t maxAttempts=8;
    std::size_t maxPasses=512;
    std::size_t maxBlockers=256;
    std::size_t maxScratchRecoveries=4096;
    std::size_t maxCellRows=1U<<20;
    std::size_t maxSupportCellRows=1U<<20;
    std::size_t maxVertexRows=1U<<18;
    std::size_t maxEdgeRows=1U<<20;
    std::size_t maxFaceRows=1U<<20;
    std::size_t maxFaceOwnerRows=1U<<20;
    std::size_t maxCandidateRows=1U<<18;
    std::size_t maxBoundEdgeRows=1U<<20;
    std::size_t maxPassRows=4096;
    std::size_t maxAssessmentRows=1U<<20;
    std::size_t maxDecisionRows=1U<<22;
    std::size_t maxBlockerIndexRows=1U<<20;
    std::uint64_t maxPairTests=UINT64_C(1)<<32;
    std::uint64_t maxCoordinatorWork=UINT64_C(1)<<34;
    std::uint64_t maxNestedRecoveryWork=UINT64_C(1)<<34;
    std::size_t maxAccountedBytes=std::size_t{1}<<30;

    bool positive() const noexcept {
        return maxAttempts&&maxPasses&&maxBlockers&&maxScratchRecoveries&&
            maxCellRows&&maxSupportCellRows&&maxVertexRows&&maxEdgeRows&&
            maxFaceRows&&maxFaceOwnerRows&&maxCandidateRows&&
            maxBoundEdgeRows&&maxPassRows&&
            maxAssessmentRows&&maxDecisionRows&&maxBlockerIndexRows&&
            maxPairTests&&maxCoordinatorWork&&maxNestedRecoveryWork&&
            maxAccountedBytes;
    }
};

struct T2CEdgeLineageCandidateCellForTesting {
    std::array<DenseTetNodeId,4> denseNodes{};
    mesh::Classification classification{};
    friend bool operator==(const T2CEdgeLineageCandidateCellForTesting&,
                           const T2CEdgeLineageCandidateCellForTesting&)=
        default;
};

struct T2CEdgeLineageSupportCellForTesting {
    topo::CellHandle handle{};
    std::uint64_t sequence=0;
    std::array<DenseTetNodeId,4> denseNodes{};
    mesh::Classification classification{};
    friend bool operator==(const T2CEdgeLineageSupportCellForTesting&,
                           const T2CEdgeLineageSupportCellForTesting&)=
        default;
};

struct T2CEdgeLineageVertexForTesting {
    DenseTetNodeId dense=kInvalidDenseTetNode;
    topo::NodeHandle handle{};
    Vec3 position{};
    mesh::Classification classification{};
};

// Complete, flat support rows.  Offsets point into the owning sink, never
// into the disposable recovery workspace.
struct T2CEdgeLineageSupportEdgeForTesting {
    DenseTetEdgeKey key{};
    topo::EdgeHandle handle{};
    mesh::Classification classification{};
    std::size_t completeStarOffset=0,completeStarCount=0;
    bool required=false;
};

struct T2CEdgeLineageSupportFaceOwnerForTesting {
    topo::CellHandle cell{};
    std::int8_t sense=0;
};

struct T2CEdgeLineageSupportFaceForTesting {
    DenseTetFaceKey key{};
    topo::FaceHandle handle{};
    mesh::Classification classification{};
    std::size_t ownerOffset=0,ownerCount=0;
};

enum class T2CEdgeLineageDependencyKindForTesting : std::uint8_t {
    BlockerToBlocker,
    CandidateToFrozenBlocker,
    CandidateToCandidate,
    PreparedDeltaToFrozenBlocker
};

struct T2CEdgeLineageAssessmentForTesting {
    VolumeArtifactEdgeKey key{};
    algo::ArtifactStamp artifactStamp{};
    mesh::Classification classification{};
    mesh::Classification region{};
    std::array<DenseTetNodeId,2> targetDense{
        kInvalidDenseTetNode,kInvalidDenseTetNode};
    TetEdgeRecoveryOutcome outcome{};
    std::size_t canonicalOrdinal=0;
    std::size_t selectedCellOffset=0,selectedCellCount=0;
    std::size_t selectedPayloadOffset=0,selectedPayloadCount=0;
    bool mapped=false;
    bool excludedFrozen=false;
    bool prefixBoundary=false;
    bool candidate=false;
    bool deferredAtEntry=false;
    bool deferred=false;
    bool preparedDeferred=false;
};

struct T2CEdgeLineageDependencyForTesting {
    T2CEdgeLineageDependencyKindForTesting kind=
        T2CEdgeLineageDependencyKindForTesting::CandidateToCandidate;
    // Row indices are assessment rows for candidate kinds and blocker rows
    // for BlockerToBlocker.
    std::size_t sourceRow=0;
    // Assessment-row index for CandidateToCandidate; blocker-row index for
    // the other two kinds.
    std::size_t targetRow=0;
    bool touches=false;
    bool retiredTouches=false;
    bool candidateTouches=false;
};

struct T2CEdgeLineagePassForTesting {
    topo::WorkspaceStamp stamp{};
    std::size_t assessmentOffset=0,assessmentCount=0;
    std::size_t decisionOffset=0,decisionCount=0;
    std::size_t frozenBlockerOffset=0,frozenBlockerCount=0;
    std::size_t permanentAssessment=(std::numeric_limits<std::size_t>::max)();
    std::size_t selectedAssessment=(std::numeric_limits<std::size_t>::max)();
    std::size_t discoveredBlocker=(std::numeric_limits<std::size_t>::max)();
    std::size_t preparedRecovery=(std::numeric_limits<std::size_t>::max)();
    bool canonicalScanComplete=false;
    bool edgeFixedPoint=false;
};

struct T2CEdgeLineageBlockerForTesting {
    VolumeArtifactEdgeKey key{};
    topo::WorkspaceStamp observedStamp{};
    topo::WorkspaceStamp finalValidatedStamp{};
    algo::ArtifactStamp artifactStamp{};
    mesh::Classification classification{};
    mesh::Classification region{};
    std::array<DenseTetNodeId,2> targetDense{};
    TetEdgeRecoveryOutcome outcome{};
    std::size_t selectedCellOffset=0,selectedCellCount=0;
    std::size_t supportCellOffset=0,supportCellCount=0;
    std::size_t supportPayloadOffset=0,supportPayloadCount=0;
    std::size_t supportEdgeOffset=0,supportEdgeCount=0;
    std::size_t supportFaceOffset=0,supportFaceCount=0;
    std::uint64_t supportDigest=0;
    std::uint64_t revalidations=0;
    bool pairwiseIndependent=true;
    bool deferred=false;
};

struct T2CEdgeLineageRecoveryForTesting {
    VolumeArtifactEdgeKey key{};
    topo::WorkspaceStamp beforeStamp{};
    topo::WorkspaceStamp afterStamp{};
    topo::EdgeHandle predictedEdge{};
    TetEdgeRecoveryOutcome outcome{};
    std::size_t retiredCellOffset=0,retiredCellCount=0;
    std::size_t retiredPayloadOffset=0,retiredPayloadCount=0;
    std::size_t candidateOffset=0,candidateCount=0;
    std::size_t createdCellOffset=0,createdCellCount=0;
    std::size_t boundEdgeOffset=0,boundEdgeCount=0;
    std::uint64_t blockerPairTests=0;
    bool allowed=false;
    bool samePreparedTransaction=true;
    bool predictedBindingAlias=false;
};

struct T2CEdgeLineageBoundEdgeForTesting {
    topo::EdgeHandle handle{};
    VolumeArtifactEdgeKey key{};
};

struct T2CEdgeLineageAttemptForTesting {
    topo::WorkspaceStamp initialStamp{};
    topo::WorkspaceStamp finalStamp{};
    algo::ArtifactStamp artifactStamp{};
    T2CEdgeLineageStatusForTesting status=
        T2CEdgeLineageStatusForTesting::None;
    std::size_t blockerOffset=0,blockerCount=0;
    std::size_t recoveryOffset=0,recoveryCount=0;
    std::size_t passOffset=0,passCount=0;
    std::size_t vertexOffset=0,vertexCount=0;
    std::size_t passes=0,deferredBlockers=0,deferredRecoveries=0;
    std::uint64_t coordinatorWork=0,nestedRecoveryWork=0;
    std::size_t accountedBytes=0;
    bool edgeFixedPoint=false;
    bool frontierReplaced=false;
    bool complete=false;
};

struct T2CEdgeLineageSinkForTesting {
    T2CEdgeLineageLimitsForTesting limits{};
    std::vector<T2CEdgeLineageAttemptForTesting> attempts;
    std::vector<T2CEdgeLineageBlockerForTesting> blockers;
    std::vector<T2CEdgeLineageRecoveryForTesting> recoveries;
    std::vector<topo::CellHandle> cellRows;
    std::vector<T2CEdgeLineageSupportCellForTesting> supportCellRows;
    std::vector<T2CEdgeLineageVertexForTesting> vertexRows;
    std::vector<T2CEdgeLineageSupportEdgeForTesting> edgeRows;
    std::vector<T2CEdgeLineageSupportFaceForTesting> faceRows;
    std::vector<T2CEdgeLineageSupportFaceOwnerForTesting> faceOwnerRows;
    std::vector<T2CEdgeLineageCandidateCellForTesting> candidateRows;
    std::vector<T2CEdgeLineageBoundEdgeForTesting> boundEdgeRows;
    std::vector<T2CEdgeLineagePassForTesting> passRows;
    std::vector<T2CEdgeLineageAssessmentForTesting> assessmentRows;
    std::vector<T2CEdgeLineageDependencyForTesting> decisionRows;
    std::vector<std::size_t> blockerIndexRows;
    std::uint64_t blockerPairTests=0;
    std::uint64_t candidatePairTests=0;
    std::uint64_t transactionPairTests=0;
    std::uint64_t coordinatorWork=0;
    std::uint64_t nestedRecoveryWork=0;
    std::size_t accountedByteHighWater=0;
    T2CEdgeLineageStatusForTesting status=
        T2CEdgeLineageStatusForTesting::None;
    bool complete=true;
    bool allBlockersPairwiseIndependent=true;
    bool allPreparedTransactionsGuarded=true;
    // The trace separates measured counters but intentionally does not claim
    // section-7-complete work/memory accounting for STL comparisons,
    // predicates, nested transaction scratch, or workspace/context growth.
    bool workAccountingComplete=false;
    bool memoryAccountingComplete=false;

    void reset() noexcept {
        attempts.clear();blockers.clear();recoveries.clear();cellRows.clear();
        supportCellRows.clear();vertexRows.clear();edgeRows.clear();
        faceRows.clear();faceOwnerRows.clear();candidateRows.clear();
        boundEdgeRows.clear();
        passRows.clear();assessmentRows.clear();decisionRows.clear();
        blockerIndexRows.clear();
        blockerPairTests=0;candidatePairTests=0;transactionPairTests=0;
        coordinatorWork=0;
        nestedRecoveryWork=0;accountedByteHighWater=0;
        status=T2CEdgeLineageStatusForTesting::None;
        complete=true;allBlockersPairwiseIndependent=true;
        allPreparedTransactionsGuarded=true;
        workAccountingComplete=false;memoryAccountingComplete=false;
    }
};

struct T2CEdgeLineageRunForTesting {
    bool invoked=false;
    bool frontierReplaced=false;
    bool complete=true;
    bool workAccountingComplete=false;
    bool memoryAccountingComplete=false;
    std::uint64_t coordinatorWork=0;
    std::uint64_t nestedRecoveryWork=0;
};

void setT2CEdgeLineageSinkForTesting(
    T2CEdgeLineageSinkForTesting* sink) noexcept;
T2CEdgeLineageSinkForTesting* t2cEdgeLineageSinkForTesting() noexcept;

T2CEdgeLineageRunForTesting traceT2CEdgeLineageForTesting(
    const ClosedTri3Artifact& artifact,
    const std::map<Index,topo::NodeHandle>& localByGlobal,
    mesh::Classification temporaryRegion,TetFacetSetState& state,
    topo::TopologyWorkspace& workspace,TetRecoveryContext& context,
    const TetFacetSetRecoveryOutcome& original,
    TetFacetSetFrontierSink& replacementFrontier) noexcept;

inline constexpr std::size_t kTetFacetConflictDiagnosticMemberCapacity = 4096;
inline constexpr std::size_t kTetFacetConflictDiagnosticEdgeCapacity = 16384;

// Fixed-capacity, allocation-free observation of the canonical facet SCC
// selected by the recovery scheduler.  The ordinary v3/v4 path never enables
// this seam.  Capture reuses the scheduler's already-computed assessments and
// graph; it performs no extra topology query or global scan.
struct TetFacetConflictCycleDiagnostics {
    bool enabled = false;
    bool captured = false;
    bool complete = true;
    topo::WorkspaceStamp workspaceStamp{};
    algo::ArtifactStamp artifactStamp{};
    std::size_t facetCount = 0;
    std::size_t pendingCount = 0;
    std::size_t memberCount = 0;
    std::size_t dependencyCount = 0;
    std::array<VolumeFacetKey,
               kTetFacetConflictDiagnosticMemberCapacity> members{};
    std::array<std::array<VolumeFacetKey,2>,
               kTetFacetConflictDiagnosticEdgeCapacity> dependencies{};
    TetFacetSetSchedulerMetrics schedulerMetrics{};
    TetFacetGeometryCounters allAssessmentCounters{};
    TetFacetGeometryCounters pendingAssessmentCounters{};
    std::uint64_t allRetiredCellRows = 0;
    std::uint64_t pendingRetiredCellRows = 0;
    std::size_t maximumRetiredCells = 0;
};

void resetTetFacetConflictCycleDiagnosticsForTesting(
    bool enabled) noexcept;
const TetFacetConflictCycleDiagnostics&
tetFacetConflictCycleDiagnosticsForTesting() noexcept;

bool tetRecoveryBlockerGraphValidForTesting(
    const std::vector<std::vector<std::size_t>>& graph,
    const std::vector<bool>& pending) noexcept;

std::vector<std::size_t> canonicalTetRecoveryConflictCycleForTesting(
    const std::vector<std::vector<std::size_t>>& graph,
    const std::vector<bool>& pending,std::uint64_t& visits);

void fillTetFacetConflictCycleSinkForTesting(
    std::span<const VolumeFacetKey> members,
    std::span<const std::array<VolumeFacetKey,2>> dependencies,
    topo::WorkspaceStamp workspaceStamp,algo::ArtifactStamp artifactStamp,
    TetFacetConflictCycleSink& sink) noexcept;
#endif

static_assert(std::is_standard_layout_v<VolumeFacetKey>);
static_assert(std::is_trivially_copyable_v<VolumeFacetKey>);
static_assert(std::is_standard_layout_v<VolumeRegionUse>);
static_assert(std::is_trivially_copyable_v<VolumeRegionUse>);
static_assert(std::is_standard_layout_v<VolumeArtifactEdgeKey>);
static_assert(std::is_trivially_copyable_v<VolumeArtifactEdgeKey>);
static_assert(std::is_standard_layout_v<TetFacetSetSchedulerMetrics>);
static_assert(std::is_trivially_copyable_v<TetFacetSetSchedulerMetrics>);
static_assert(std::is_standard_layout_v<TetFacetConflictCycleSink>);
static_assert(std::is_trivially_copyable_v<TetFacetConflictCycleSink>);

}  // namespace paacs::bnd
