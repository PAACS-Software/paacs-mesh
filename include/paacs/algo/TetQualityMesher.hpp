// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T4-A: bounded Tet4 quality / worst-tail pass over an accepted T3-A result.

#include "paacs/algo/TetRefinementMesher.hpp"
#include "paacs/mesh/MeshAudit.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace paacs::algo {

inline constexpr const char* kT4ARouteId = "quality-tet4-v4";
inline constexpr std::uint32_t kT4ARouteVersion = 1;
inline constexpr std::uint64_t kCurrentT4Operation =
    (std::numeric_limits<std::uint64_t>::max)();

enum class T4AQualityMode : std::uint8_t { ReportOnly, BoundedImprove };

struct T4AQualityPolicy {
    T4AQualityMode mode = T4AQualityMode::ReportOnly;
    double minimumRadiusRatio = 0.0;
    double minimumImprovement = 0.0;
    std::uint64_t maxOperations = 0;
    std::uint64_t maxSelectedGenerations = 0;
    std::uint32_t maxCandidateAttemptsPerCell = 0;
    std::uint32_t maxCreatedTetsPerOperation = 0;
    std::uint64_t maxWorkUnits = 0;
    std::uint64_t maxNewNodes = 0;
    std::uint64_t maxResultNodes = 0;
    std::uint64_t maxResultTets = 0;
    std::uint64_t maxReportRows = 0;
    std::uint64_t maxWorstTailRows = 0;
    std::uint32_t maxDependencyRowsPerAttempt = 0;
    std::uint32_t maxProtectedWitnessesPerAttempt = 0;
    std::uint32_t maxTopologyRowsPerOperation = 0;
};

struct T4ARequest {
    T3ARequest source{};
    T4AQualityPolicy quality{};
};

enum class T4AStatus : std::uint8_t {
    Published, Incomplete, Rejected, Unsupported, Failed
};

enum class T4APhase : std::uint8_t {
    ReportShell, RequestPreflight, SourceT3A, ReportEnvelope,
    ContextBuild, QueueBuild, CandidateGeneration, Encroachment,
    CavityDiscovery, CandidateValidation, TransactionPreparation,
    CoordinatorPreparation, LocalAudit, FinalAudit,
    PublicationPreparation, Publication
};

enum class T4AReason : std::uint16_t {
    None, InvalidRequest, CapacityOverflow, SourceT3AStopped,
    InvalidAuthoritativeSource, InvalidQualityValue,
    NoImprovingCavity, ProtectedFeatureBlocked,
    OperationBudgetReached, SelectedGenerationBudgetReached,
    WorkBudgetReached, NodeBudgetReached, ResultBudgetReached,
    ReportBudgetReached, ResidualQuality, InvalidCandidate, StaleContext,
    TopologyPreparationFailed, CoordinatorPreparationFailed,
    QualityAuditFailed, BoundaryBijectionFailed,
    PublicationPreparationFailed, AllocationFailed, InternalError
};

struct T4AQualityKey {
    std::uint64_t radiusRatioBits = 0;
    T3ATetKey tet{};
    std::uint64_t generation = 0;
    friend bool operator==(const T4AQualityKey&,
                           const T4AQualityKey&) = default;
    friend bool operator<(const T4AQualityKey& a,
                          const T4AQualityKey& b) noexcept {
        if (a.radiusRatioBits != b.radiusRatioBits)
            return a.radiusRatioBits < b.radiusRatioBits;
        if (a.tet < b.tet) return true;
        if (b.tet < a.tet) return false;
        return a.generation < b.generation;
    }
};

enum class T4ADependencyKind : std::uint8_t { Cell, Edge, Face };
struct T4ADependencyKey {
    T4ADependencyKind kind = T4ADependencyKind::Cell;
    std::array<RefinementVertexKey,4> vertices{};
    std::uint8_t vertexCount = 0;
    friend bool operator==(const T4ADependencyKey&,
                           const T4ADependencyKey&) = default;
    friend bool operator<(const T4ADependencyKey& a,
                          const T4ADependencyKey& b) noexcept {
        if (a.kind!=b.kind) return a.kind<b.kind;
        if (a.vertexCount!=b.vertexCount)
            return a.vertexCount<b.vertexCount;
        return a.vertices<b.vertices;
    }
};

struct T4AQualitySummary {
    std::uint64_t tetCount = 0, exactPositive = 0, nonfinite = 0;
    std::uint64_t belowThreshold = 0;
    double minimum = 0.0, p1 = 0.0, p5 = 0.0, p50 = 0.0;
    std::array<std::uint64_t,10> histogram{};
    std::uint64_t tailOffset = 0, tailCount = 0;
    bool tailTruncated = false;
};

enum class T4AAttemptDisposition : std::uint8_t {
    FormulaUnavailable, DuplicateProposal, Encroached, InvalidCandidate,
    NoQualityProgress, Committed, Blocked, TerminalFailure
};

struct T4AAttemptRow {
    std::uint64_t attemptOrdinal = 0, operationOrdinal = 0;
    T4AQualityKey selected{};
    Vec3 proposedCoordinate{};
    std::uint64_t removedCells = 0, createdCells = 0;
    std::uint64_t dependencyOffset = 0, dependencyCount = 0;
    std::uint64_t witnessOffset = 0, witnessCount = 0;
    double selectedQuality = 0.0, createdMinimumQuality = 0.0;
    T4AAttemptDisposition disposition =
        T4AAttemptDisposition::TerminalFailure;
    T4AReason reason = T4AReason::None;
    std::uint64_t workBefore = 0, workAfter = 0;
};

struct T4ABlockedRow {
    T4AQualityKey selected{};
    T4AReason reason = T4AReason::None;
    std::uint64_t attemptOrdinal = 0;
    std::uint64_t dependencyOffset = 0, dependencyCount = 0;
    std::uint64_t witnessOffset = 0, witnessCount = 0;
    std::uint64_t invalidatedByOperation = kCurrentT4Operation;
};

struct T4ATerminalRow {
    T4AQualityKey selected{};
    T4AReason reason = T4AReason::None;
    T4APhase phase = T4APhase::CandidateGeneration;
    std::uint64_t attemptOrdinal = 0;
};

struct T4ACommitRow {
    std::uint64_t operationOrdinal = 0, attemptOrdinal = 0;
    T4AQualityKey selected{};
    topo::NodeHandle predictedNode{}, actualNode{};
    std::uint64_t topologyOffset = 0, topologyCount = 0;
    double selectedQuality = 0.0, createdMinimumQuality = 0.0;
    topo::WorkspaceStamp postCommitStamp{};
};

struct T4AFailure {
    T4AReason reason = T4AReason::None;
    T4APhase phase = T4APhase::RequestPreflight;
    bool hasCell = false;
    T3ATetKey cell{};
};

struct T4AMetrics {
    std::uint64_t initialCellsScanned = 0, finalCellsAudited = 0;
    std::uint64_t qualityEvaluations = 0, queuePushes = 0, queuePops = 0;
    std::uint64_t staleRows = 0, queueHighWater = 0;
    std::uint64_t selectedGenerations = 0, candidateAttempts = 0;
    std::uint64_t formulaRejects = 0, duplicateRejects = 0;
    std::uint64_t protectedIndexVisits = 0, locationIndexVisits = 0;
    std::uint64_t locationTetTests = 0, cavityCellsVisited = 0;
    std::uint64_t cavityFaceVisits = 0, survivorIndexVisits = 0;
    std::uint64_t candidateSurvivorTests = 0;
    std::uint64_t localRowsPrepared = 0, localRowsApplied = 0;
    std::uint64_t blockedRows = 0, blockedInvalidations = 0;
    std::uint64_t preparedTransactions = 0, committedOperations = 0;
    std::uint64_t perOperationGlobalTetScans = 0, contextRebuilds = 0;
    std::uint64_t phase3Allocations = 0, publicationTailAllocations = 0;
    std::uint64_t workUnitsConsumed = 0;
    std::uint64_t transientAccountedHighWaterBytes = 0;
    bool runComplete = false, qualityComplete = false;
};

struct T4AReport {
    T4AReport() = default;
    T4AReport(const T4AReport&) = delete;
    T4AReport& operator=(const T4AReport&) = delete;
    T4AReport(T4AReport&&) = delete;
    T4AReport& operator=(T4AReport&&) = delete;
    bool shellPrepared = false, reportComplete = false;
    T4AQualitySummary before{}, after{};
    std::vector<T4AAttemptRow> attemptStorage;
    std::vector<T4ABlockedRow> blockedStorage;
    std::vector<T4ACommitRow> commitStorage;
    std::vector<T4ADependencyKey> dependencyStorage;
    std::vector<T3AProtectedFeatureWitness> witnessStorage;
    std::vector<T3ATopologyHandleRow> topologyStorage;
    std::vector<T4AQualityKey> unattemptedStorage;
    std::vector<T4AQualityKey> tailStorage;
    std::uint64_t attemptsUsed = 0, blockedUsed = 0, commitsUsed = 0;
    std::uint64_t dependenciesUsed = 0, witnessesUsed = 0;
    std::uint64_t topologyUsed = 0, unattemptedUsed = 0, tailUsed = 0;
    bool hasTerminal = false;
    T4ATerminalRow terminal{};
    std::span<const T4AAttemptRow> attempts() const noexcept {
        return {attemptStorage.data(),static_cast<std::size_t>(attemptsUsed)};
    }
    std::span<const T4ABlockedRow> blocked() const noexcept {
        return {blockedStorage.data(),static_cast<std::size_t>(blockedUsed)};
    }
};

struct T4AT3SourceReceipt {
    T4AT3SourceReceipt() = default;
    T4AT3SourceReceipt(const T4AT3SourceReceipt&) = delete;
    T4AT3SourceReceipt& operator=(const T4AT3SourceReceipt&) = delete;
    T4AT3SourceReceipt(T4AT3SourceReceipt&&) = delete;
    T4AT3SourceReceipt& operator=(T4AT3SourceReceipt&&) = delete;
    const char* routeId = kT3ARouteId;
    std::uint32_t routeVersion = kT3ARouteVersion;
    T3AStatus status = T3AStatus::Published;
    T3AFailure failure{};
    bool reportComplete = false;
    std::unique_ptr<BoundaryConformingV3Receipt> sourceReceipt;
    std::unique_ptr<T3AReport> sourceReport;
    TetRefinementMetrics metrics{};
    TetRefinementQualityReport quality{};
    BoundaryConformingV3Audit boundaryAudit{};
    std::uint64_t strictHash = 0;
    mesh::SemanticHash semanticHash{};
    std::uint64_t nodeCount = 0, cellCount = 0;
    std::uint64_t nodeMappingDigest = 0, cellMappingDigest = 0;
};

struct T4AMutableArtifact {
    T4AMutableArtifact() = default;
    T4AMutableArtifact(const T4AMutableArtifact&) = delete;
    T4AMutableArtifact& operator=(const T4AMutableArtifact&) = delete;
    T4AMutableArtifact(T4AMutableArtifact&&) = delete;
    T4AMutableArtifact& operator=(T4AMutableArtifact&&) = delete;
    std::unique_ptr<TetAuthoritativeResult> authoritative;
    bnd::TetFacetSetState constraintState;
    std::vector<Index> existingGlobalByAuthoritativeNodeSlot;
    std::vector<Index> sourceGlobalByAuthoritativeNodeSlot;
    std::vector<Index> sourceGlobalByAuthoritativeCellSlot;
    std::vector<Index> publishedGlobalByAuthoritativeNodeSlot;
    std::vector<Index> publishedGlobalByAuthoritativeCellSlot;
    std::vector<T3ABoundaryTriRecord> immutableBoundary;
    std::map<topo::NodeHandle,RefinementVertexKey> vertexKeyByHandle;
    std::uint64_t t3CommittedInsertionCount = 0;
};

struct T4AResult {
    T4AResult() = default;
    T4AResult(const T4AResult&) = delete;
    T4AResult& operator=(const T4AResult&) = delete;
    T4AResult(T4AResult&&) = delete;
    T4AResult& operator=(T4AResult&&) = delete;
    const char* routeId = kT4ARouteId;
    std::uint32_t routeVersion = kT4ARouteVersion;
    std::unique_ptr<T4AMutableArtifact> artifact;
    T4AQualitySummary before{}, after{};
    T4AMetrics metrics{};
    BoundaryConformingV3Audit boundaryAudit{};
    std::uint64_t finalStrictHash = 0;
    mesh::SemanticHash finalSemanticHash{};
    bool publicationReady = false;
};

struct T4AOutcome {
    T4AOutcome() = default;
    T4AOutcome(const T4AOutcome&) = delete;
    T4AOutcome& operator=(const T4AOutcome&) = delete;
    T4AOutcome(T4AOutcome&&) noexcept = default;
    T4AOutcome& operator=(T4AOutcome&&) noexcept = default;
    T4AStatus status = T4AStatus::Failed;
    T4AFailure failure{};
    bool destinationUnchanged = true, reportComplete = false;
    std::unique_ptr<T3AOutcome> sourceOutcome;
    std::unique_ptr<T4AT3SourceReceipt> sourceReceipt;
    std::unique_ptr<T4AReport> report;
    std::unique_ptr<T4AResult> retainedResult;
    bool published() const noexcept { return status==T4AStatus::Published; }
};

static_assert(std::is_standard_layout_v<T4AQualityPolicy>);
static_assert(std::is_trivially_copyable_v<T4AQualityPolicy>);
static_assert(std::is_standard_layout_v<T4AQualityKey>);
static_assert(std::is_trivially_copyable_v<T4AQualityKey>);
static_assert(std::is_standard_layout_v<T4ADependencyKey>);
static_assert(std::is_trivially_copyable_v<T4ADependencyKey>);
static_assert(std::is_standard_layout_v<T4AQualitySummary>);
static_assert(std::is_trivially_copyable_v<T4AQualitySummary>);
static_assert(std::is_standard_layout_v<T4AAttemptRow>);
static_assert(std::is_trivially_copyable_v<T4AAttemptRow>);
static_assert(std::is_standard_layout_v<T4ABlockedRow>);
static_assert(std::is_trivially_copyable_v<T4ABlockedRow>);
static_assert(std::is_standard_layout_v<T4ATerminalRow>);
static_assert(std::is_trivially_copyable_v<T4ATerminalRow>);
static_assert(std::is_standard_layout_v<T4ACommitRow>);
static_assert(std::is_trivially_copyable_v<T4ACommitRow>);
static_assert(std::is_standard_layout_v<T4AFailure>);
static_assert(std::is_trivially_copyable_v<T4AFailure>);
static_assert(std::is_standard_layout_v<T4AMetrics>);
static_assert(std::is_trivially_copyable_v<T4AMetrics>);
static_assert(std::is_nothrow_move_constructible_v<T4AOutcome>);
static_assert(std::is_nothrow_move_assignable_v<T4AOutcome>);
static_assert(!std::is_move_constructible_v<T4AReport>);
static_assert(!std::is_move_constructible_v<T4AT3SourceReceipt>);
static_assert(!std::is_move_constructible_v<T4AMutableArtifact>);
static_assert(!std::is_move_constructible_v<T4AResult>);
static_assert(noexcept(std::declval<std::vector<Index>&>().swap(
    std::declval<std::vector<Index>&>())));
static_assert(noexcept(std::declval<std::vector<T3ABoundaryTriRecord>&>().swap(
    std::declval<std::vector<T3ABoundaryTriRecord>&>())));
static_assert(noexcept(std::declval<std::map<topo::NodeHandle,
    RefinementVertexKey>&>().swap(std::declval<std::map<topo::NodeHandle,
    RefinementVertexKey>&>())));

T4AOutcome meshModelQualityTet4(
    const geo::Model& model,const sizing::MeshSpec& spec,
    const T4ARequest& request,mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& ownershipProvider) noexcept;

#if defined(PAACS_TOPO_TESTING)
namespace t4a_detail {
inline bool& allocationScopeForTesting() noexcept {
    static thread_local bool value=false;
    return value;
}
inline bool& forceCandidateForTesting() noexcept {
    static thread_local bool value=false;
    return value;
}
inline Vec3& forcedCandidateForTesting() noexcept {
    static thread_local Vec3 value{};
    return value;
}
inline bool& injectPostCommitAuditFailureForTesting() noexcept {
    static thread_local bool value=false;
    return value;
}
enum class T4ADigestCheckpoint : std::uint8_t {
    PreSelection,
    SelectedRowSecured,
    PrePreparation,
    Prepared
};
struct T4AStateDigest {
    T4ADigestCheckpoint checkpoint=T4ADigestCheckpoint::PreSelection;
    bool postFailure=false;
    std::uint64_t topology=0,pools=0,prepared=0,mappings=0,resolvers=0,
        location=0,queue=0,blockers=0,dependencyIndex=0,
        report=0,metrics=0,structural=0;
    std::uint64_t revision=0,attemptsUsed=0,blockedUsed=0,
        commitsUsed=0,dependenciesUsed=0,witnessesUsed=0,
        topologyUsed=0,unattemptedUsed=0;
};
struct T4ADigestTrace {
    bool enabled=false,faultArmed=false,faultPending=false;
    T4ADigestCheckpoint faultCheckpoint=T4ADigestCheckpoint::PreSelection;
    std::uint32_t used=0;
    std::array<T4AStateDigest,24> rows{};
    void reset(bool enable=true) noexcept {
        enabled=enable;
        faultArmed=false;
        faultPending=false;
        faultCheckpoint=T4ADigestCheckpoint::PreSelection;
        used=0;
        rows={};
    }
    void arm(T4ADigestCheckpoint checkpoint) noexcept {
        reset(true);
        faultArmed=true;
        faultCheckpoint=checkpoint;
    }
};
inline T4ADigestTrace& digestTraceForTesting() noexcept {
    static thread_local T4ADigestTrace value;
    return value;
}
struct T4AFinalLiveTrace {
    bool enabled=false,complete=false;
    std::uint32_t used=0;
    std::array<T4AQualityKey,8192> rows{};
    void reset(bool value=true) noexcept {
        enabled=value;
        complete=false;
        used=0;
        rows={};
    }
};
inline T4AFinalLiveTrace& finalLiveTraceForTesting() noexcept {
    static thread_local T4AFinalLiveTrace value;
    return value;
}
bool checkedTetQualityForTesting(
    std::array<Vec3,4> coordinates,
    std::array<RefinementVertexKey,4> keys,double& quality) noexcept;
}
#endif

} // namespace paacs::algo
