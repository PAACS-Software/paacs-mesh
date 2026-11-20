// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T2-C: additive shape-preserving protected-feature refinement over the
// accepted boundary-conforming-v3 route. The accepted v3 entry remains the
// default and never observes this policy.

#include "paacs/algo/BoundaryConformingMesher.hpp"
#include "paacs/sizing/MeshControls.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <type_traits>
#include <vector>

namespace paacs::algo {

inline constexpr const char* kBoundaryConformingV4RouteId =
    "boundary-conforming-v4";
inline constexpr std::uint32_t kBoundaryConformingV4RouteVersion = 4;

struct ProtectedFeatureProviderStamp {
    std::uint64_t persistentId = 0;
    std::uint64_t revision = 0;
    friend bool operator==(const ProtectedFeatureProviderStamp&,
                           const ProtectedFeatureProviderStamp&) = default;
};

enum class ProtectedFeaturePlacementStatus : std::uint8_t {
    Placed,
    Rejected,
    Unsupported,
    Failed
};

struct ProtectedFeaturePlacement {
    ProtectedFeatureProviderStamp providerStamp{};
    ArtifactStamp entityStamp{};
    SurfaceNodeIdentity requestedIdentity{};
    ProtectedFeaturePlacementStatus status =
        ProtectedFeaturePlacementStatus::Failed;
    Vec3 coordinate{};
};

enum class ProtectedFeatureEntityKind : std::uint8_t {
    Curve,
    FaceEdge,
    FacetInterior
};

struct ProtectedFeaturePlacementRequest {
    ProtectedFeatureEntityKind kind = ProtectedFeatureEntityKind::Curve;
    ArtifactStamp entityStamp{};
    SurfaceNodeIdentity requestedIdentity{};
    int entityTag = -1;
    double curveParameter = 0.0;
    Vec2 faceParameter{};
};

class ProtectedFeaturePlacementProvider {
public:
    virtual ~ProtectedFeaturePlacementProvider() = default;
    virtual ProtectedFeatureProviderStamp stamp() const noexcept = 0;
    virtual ProtectedFeaturePlacement placeOnEdge(
        const geo::GEdge& edge,
        const ProtectedFeaturePlacementRequest& request) const = 0;
    virtual ProtectedFeaturePlacement placeOnFace(
        const geo::GFace& face,
        const ProtectedFeaturePlacementRequest& request) const = 0;
};

struct ProtectedFeatureRefinementLimits {
    std::uint64_t maxEdgeSplits = 0;
    std::uint64_t maxRounds = 0;
    std::uint64_t maxAddedNodes = 0;
    std::uint64_t maxTetCells = 0;
    std::uint64_t maxAccountedBytes = 0;
    // This bounds coordinator/frontier/provider bookkeeping only.  Complete
    // nested v3 construction and recovery work is measured separately.
    std::uint64_t maxCoordinatorWorkUnits = 0;
    std::uint64_t maxBatchSize = 1;
    std::uint64_t maxFrontierSupportRows = 4096;
    std::uint64_t maxFrontierDependencyRows = 4096;
    std::uint64_t maxFrontierRows = 4096;

    bool positive() const noexcept {
        return maxEdgeSplits != 0 && maxRounds != 0 &&
               maxAddedNodes != 0 && maxTetCells != 0 &&
               maxAccountedBytes != 0 && maxCoordinatorWorkUnits != 0 &&
               maxBatchSize != 0 && maxFrontierSupportRows != 0 &&
               maxFrontierDependencyRows != 0 && maxFrontierRows != 0;
    }
};

struct SurfaceFaceCacheLimits {
    std::uint64_t maxFaces = 4096;
    std::uint64_t maxCachedNodeRows = 4'000'000;
    std::uint64_t maxCachedCellRows = 8'000'000;
    std::uint64_t maxCachedConstraintRows = 4'000'000;
    std::uint64_t maxCachedMappingRows = 8'000'000;
    std::uint64_t maxDirtyFacesPerBatch = 4096;
    std::uint64_t maxSurfaceCacheWorkUnits = 1'000'000'000;
    std::uint64_t maxSurfaceCacheAccountedBytes = 4'000'000'000ULL;

    bool positive() const noexcept {
        return maxFaces&&maxCachedNodeRows&&maxCachedCellRows&&
            maxCachedConstraintRows&&maxCachedMappingRows&&
            maxDirtyFacesPerBatch&&maxSurfaceCacheWorkUnits&&
            maxSurfaceCacheAccountedBytes;
    }
};

enum class FacetConflictRefinementMode : std::uint8_t {
    Refuse,
    BoundedCanonicalSccSubdivision
};

// V4-only bounds for the canonical facet-conflict escape.  Zero values are
// legal while mode==Refuse and are rejected when the opt-in mode is enabled.
struct FacetConflictRefinementLimits {
    std::uint64_t maxSccMembers = 0;
    std::uint64_t maxFacetInteriorPoints = 0;
    std::uint64_t maxFacetConflictRounds = 0;
    std::uint64_t maxProviderCalls = 0;
    std::uint64_t maxCoordinatorWorkUnits = 0;

    bool positive() const noexcept {
        return maxSccMembers != 0 && maxFacetInteriorPoints != 0 &&
               maxFacetConflictRounds != 0 && maxProviderCalls != 0 &&
               maxCoordinatorWorkUnits != 0;
    }
};

struct ProtectedFeatureChildIdentity {
    std::uint64_t parentEntityPersistentId = 0;
    SurfaceNodeIdentity childIdentity{};
    std::uint64_t splitOrdinal = 0;
    int edgeTag = -1;
    int faceTag = -1;
    double parentParameterLo = 0.0;
    double parentParameterHi = 0.0;
    double midpointParameter = 0.0;
    Vec2 faceParameter{};
    Vec3 coordinate{};
    mesh::Classification classification{};
    ProtectedFeatureProviderStamp providerStamp{};
};

enum class BoundaryConformingV4Status : std::uint8_t {
    Published,
    Rejected,
    Unsupported,
    Refused,
    Incomplete,
    Failed
};

enum class BoundaryConformingV4Reason : std::uint8_t {
    None,
    InvalidRequest,
    SourceRejected,
    SourceUnsupported,
    SourceRefused,
    SourceFailed,
    ProviderRejected,
    ProviderUnsupported,
    ProviderFailed,
    ProviderRevisionChanged,
    PlacementContractViolation,
    SharedChartReconciliationUnsupported,
    SubdivisionDepthUnsupported,
    FrontierCapacityReached,
    SubdivisionBudgetReached,
    NoProgressCycle,
    ChildEdgeRecoveryRefused,
    ChildFacetRecoveryRefused,
    AllocationFailed,
    BoundaryBijectionFailed,
    SurfaceCacheAuthorityChanged,
    SurfaceCacheInvalid,
    SurfaceCacheBudgetReached,
    InternalError,
    FacetConflictBudgetReached,
    FacetConflictSubdivisionUnsupported,
    InterfacePrefillBudgetReached,
    InterfacePrefillAuthorityChanged,
    InterfacePrefillUnsupported,
    SinglePassRecoveryIncomplete
};

enum class BoundaryConformingV4Phase : std::uint8_t {
    None,
    Preflight,
    SourceAttempt,
    EdgeSelection,
    EdgePlacement,
    SurfaceRebuild,
    ConstraintRecovery,
    SurfaceCachePreparation,
    SurfaceCacheAssembly,
    BoundaryAudit,
    Publication,
    InterfacePrefill
};

struct BoundaryConformingV4Metrics {
    static constexpr std::size_t phaseCount = 11;
    std::uint64_t rounds = 0;
    std::uint64_t edgeSplits = 0;
    std::uint64_t facetSplits = 0;
    std::uint64_t providerStampCalls = 0;
    std::uint64_t providerPlacementCalls = 0;
    std::uint64_t rebuiltSurfaceNodes = 0;
    std::uint64_t rebuiltSurfaceSegments = 0;
    std::uint64_t candidateTetCells = 0;
    std::uint64_t semanticStates = 0;
    std::uint64_t frontierRows = 0;
    std::uint64_t uniquePlacements = 0;
    std::uint64_t deferredFrontierRows = 0;
    std::uint64_t nestedConstructionWorkUnits = 0;
    std::uint64_t nestedRecoveryWorkUnits = 0;
    std::uint64_t totalMeasuredWorkUnits = 0;
    std::uint64_t coordinatorWorkUnits = 0;
    std::uint64_t surfaceCacheFacesBuilt = 0;
    std::uint64_t surfaceCacheFacesReused = 0;
    std::uint64_t surfaceCacheFacesReplayed = 0;
    std::uint64_t surfaceCacheFacesRejected = 0;
    std::uint64_t surfaceCacheKeyRowsVisited = 0;
    std::uint64_t surfaceCacheDirtyPropagationVisits = 0;
    std::uint64_t surfaceCacheSizingReplayCalls = 0;
    std::uint64_t surfaceCachePreparationWorkUnits = 0;
    std::uint64_t surfaceCacheReplayWorkUnits = 0;
    std::uint64_t surfaceCacheSetAuditWorkUnits = 0;
    std::uint64_t surfaceCacheCommits = 0;
    std::uint64_t surfaceCacheCleanNodeRows = 0;
    std::uint64_t surfaceCacheDirtyNodeRows = 0;
    std::uint64_t surfaceCacheCleanCellRows = 0;
    std::uint64_t surfaceCacheDirtyCellRows = 0;
    std::uint64_t surfaceCacheCleanConstraintRows = 0;
    std::uint64_t surfaceCacheDirtyConstraintRows = 0;
    std::uint64_t surfaceCacheCleanMappingRows = 0;
    std::uint64_t surfaceCacheDirtyMappingRows = 0;
    std::uint64_t surfaceCacheCopiedBytes = 0;
    std::uint64_t surfaceCacheHighWaterSlots = 0;
    std::uint64_t surfaceCacheHighWaterRows = 0;
    std::uint64_t surfaceCacheAllFaceScansPerDirtyFace = 0;
    // Faces whose producer-owned recovered workspace was validated and
    // consumed by the surface adapter instead of being reconstructed.
    std::uint64_t preparedSurfaceRecoveryFaces = 0;
    // The direct single-pass route performs no reusable-cache construction.
    // These counters make that strategy observable without changing the
    // ordinary complete per-face surface report.
    std::uint64_t directSurfaceBuildAttempts = 0;
    std::uint64_t directSurfaceBuildFaces = 0;
    std::uint64_t facetConflictRounds = 0;
    std::uint64_t facetConflictSccMembers = 0;
    std::uint64_t facetConflictDependencyRows = 0;
    std::uint64_t facetConflictProviderCalls = 0;
    std::uint64_t facetConflictCoordinatorWorkUnits = 0;
    std::uint64_t interfacePrefillCandidateFaces = 0;
    std::uint64_t interfacePrefillFaces = 0;
    std::uint64_t interfacePrefillDenseFacesSkipped = 0;
    std::uint64_t interfacePrefillSizingDrivenFaces = 0;
    // Entity controls that request a finer interface than the base size field
    // cannot rely on the base surface density test; these faces are planned
    // explicitly and counted here.
    std::uint64_t interfacePrefillControlForcedFaces = 0;
    std::uint64_t interfacePrefillMaximumPlannedU = 0;
    std::uint64_t interfacePrefillMaximumPlannedV = 0;
    std::uint64_t interfacePrefillPoints = 0;
    std::uint64_t interfacePrefillSkippedParameters = 0;
    std::uint64_t interfacePrefillOwnershipCalls = 0;
    std::uint64_t interfacePrefillProviderCalls = 0;
    std::uint64_t interfacePrefillWorkUnits = 0;
    // Cumulative scheduling work is retained above across every complete
    // source attempt.  This companion snapshot is the exact latest nested
    // source scheduler, so its nontrivial recovery count reconciles directly
    // to the retained edge/facet commit metrics.
    bnd::TetFacetSetMixedPresentMetrics mixedPresentScheduling{};
    bnd::TetFacetSetMixedPresentMetrics
        latestSourceMixedPresentScheduling{};
    std::size_t surfaceCacheAccountedBytes = 0;
    std::size_t peakAccountedBytes = 0;
    std::array<std::uint64_t,phaseCount> nestedPhaseNanoseconds{};
    bool runComplete = false;
};

enum class TetRecoverySchedulingMode : std::uint8_t {
    AcceptedScheduler,
    MixedExactPresentBatching
};

// The default preserves the exact iterative T2-C route. The opt-in fast
// strategies pre-fill coarse two-sided interfaces, make exactly one source
// attempt, and return Incomplete rather than publishing a nonconforming mesh.
// The direct variants bypass the reusable surface cache because a one-attempt
// route cannot obtain reuse from it; the published mesh contract is identical.
enum class ProtectedFeatureExecutionStrategy : std::uint8_t {
    IterativeExact,
    SinglePassInterfacePrefill,
    SinglePassInterfacePrefillDirect,
    // V4 performance route. It retains the V2/V3 one-pass interface
    // prefill and builds the same exact-Delaunay candidate, but initializes
    // constraints directly when every requested segment is already a literal
    // candidate edge. Any non-direct case falls back to full exact recovery.
    SinglePassPreparedConstraintPrefill
};

enum class InterfacePrefillDensityMode : std::uint8_t {
    FixedGrid,
    // Derive the per-face U/V point counts once from the composite geometry
    // size plus compiled face/body controls.  The derived bounded grid is then
    // generated with the same one-pass path as V4.
    SizingDriven
};

struct InterfacePrefillPolicy {
    // V2/V3/V4 aggregate order is intentionally preserved. Existing calls
    // using {interiorU, interiorV, maxFaces, maxPoints, threshold} therefore
    // retain source and behavioural compatibility.
    std::uint32_t interiorU = 0;
    std::uint32_t interiorV = 0;
    std::uint64_t maxInterfaceFaces = 0;
    std::uint64_t maxPoints = 0;
    // If both centre-sampled surface interval estimates meet this value, the
    // face is already dense enough and is left untouched.
    double minimumEstimatedIntervals = 0.0;

    InterfacePrefillDensityMode densityMode =
        InterfacePrefillDensityMode::FixedGrid;
    std::uint32_t minimumInteriorU = 1;
    std::uint32_t minimumInteriorV = 1;
    std::uint32_t maximumInteriorU = 64;
    std::uint32_t maximumInteriorV = 64;
    double sizingScale = 1.0;

    std::uint32_t maximumPlannedU() const noexcept {
        return densityMode == InterfacePrefillDensityMode::SizingDriven
            ? maximumInteriorU : interiorU;
    }
    std::uint32_t maximumPlannedV() const noexcept {
        return densityMode == InterfacePrefillDensityMode::SizingDriven
            ? maximumInteriorV : interiorV;
    }
    bool positive() const noexcept {
        const bool densityValid = densityMode ==
                InterfacePrefillDensityMode::FixedGrid
            ? interiorU != 0 && interiorV != 0 &&
              interiorU != (std::numeric_limits<std::uint32_t>::max)() &&
              interiorV != (std::numeric_limits<std::uint32_t>::max)()
            : densityMode == InterfacePrefillDensityMode::SizingDriven &&
              minimumInteriorU != 0 && minimumInteriorV != 0 &&
              minimumInteriorU <= maximumInteriorU &&
              minimumInteriorV <= maximumInteriorV &&
              maximumInteriorU !=
                  (std::numeric_limits<std::uint32_t>::max)() &&
              maximumInteriorV !=
                  (std::numeric_limits<std::uint32_t>::max)() &&
              std::isfinite(sizingScale) && sizingScale > 0.0;
        return densityValid && maxInterfaceFaces != 0 && maxPoints != 0 &&
               std::isfinite(minimumEstimatedIntervals) &&
               minimumEstimatedIntervals > 0.0;
    }
};

enum class SurfaceCacheFaceDisposition : std::uint8_t {
    CleanReplay,
    DirtyBuild,
    Stopped
};

struct SurfaceCacheFaceReportRow {
    FacePersistentId facePersistentId = 0;
    SurfaceCacheFaceDisposition disposition =
        SurfaceCacheFaceDisposition::Stopped;
    SurfaceFaceStatus surfaceStatus = SurfaceFaceStatus::Failed;
    SurfaceFailure originalSurfaceFailure{};
};

struct SurfaceCacheReport {
    bool reportComplete = false;
    bool hasFailedFace = false;
    FacePersistentId failedFacePersistentId = 0;
    std::vector<SurfaceCacheFaceReportRow> processed;
    std::vector<FacePersistentId> unattempted;
};

struct BoundaryConformingV4Request {
    BoundaryConformingV3Request sourceRequest{};
    ProtectedFeatureProviderStamp expectedPlacementProviderStamp{};
    ProtectedFeatureRefinementLimits limits{};
    // Optional v4-private cache policy.  Defaults are positive and bounded;
    // accepted v3 never observes this field.
    SurfaceFaceCacheLimits surfaceCacheLimits{};
    FacetConflictRefinementMode facetConflictMode =
        FacetConflictRefinementMode::Refuse;
    FacetConflictRefinementLimits facetConflictLimits{};
    TetRecoverySchedulingMode recoverySchedulingMode =
        TetRecoverySchedulingMode::AcceptedScheduler;
    ProtectedFeatureExecutionStrategy executionStrategy =
        ProtectedFeatureExecutionStrategy::IterativeExact;
    InterfacePrefillPolicy interfacePrefill{};
    // Optional immutable V6 control plan.  Fixed-grid V4 requests leave this
    // null and retain the exact prior path and output hashes.
    std::shared_ptr<const sizing::CompiledMeshControls> meshControls{};
};

struct BoundaryConformingV4Result {
    const char* routeId = kBoundaryConformingV4RouteId;
    std::uint32_t routeVersion = kBoundaryConformingV4RouteVersion;
    ProtectedFeatureExecutionStrategy executionStrategy =
        ProtectedFeatureExecutionStrategy::IterativeExact;
    std::vector<ProtectedFeatureChildIdentity> subdivisionLedger;
    SurfaceCacheReport surfaceCacheReport{};
    BoundaryConformingV4Metrics metrics{};
    std::unique_ptr<BoundaryConformingV3Outcome> sourceOutcome;
};

struct BoundaryConformingV4Outcome {
    BoundaryConformingV4Status status = BoundaryConformingV4Status::Failed;
    BoundaryConformingV4Reason reason = BoundaryConformingV4Reason::InternalError;
    BoundaryConformingV4Phase phase = BoundaryConformingV4Phase::None;
    bool destinationUnchanged = true;
    std::unique_ptr<BoundaryConformingV4Result> retainedResult;
    std::unique_ptr<BoundaryConformingV3Outcome> lastSourceOutcome;

    bool published() const noexcept {
        return status == BoundaryConformingV4Status::Published;
    }
};

BoundaryConformingV4Outcome meshModelBoundaryConformingV4(
    const geo::Model& model, const sizing::MeshSpec& spec,
    const BoundaryConformingV4Request& request, mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& ownershipProvider,
    const ProtectedFeaturePlacementProvider& placementProvider) noexcept;

namespace t2c_detail {

using SurfaceFaceCacheLimits = ::paacs::algo::SurfaceFaceCacheLimits;

enum class SurfaceFaceCacheStatus : std::uint8_t {
    None, Ready, Incomplete, Failed
};
enum class SurfaceFaceCacheReason : std::uint8_t {
    None, AuthorityChanged, Invalid, BudgetReached, AllocationFailed,
    FaceBuildStopped
};
enum class SurfaceFaceCachePhase : std::uint8_t {
    None, Envelope, Preparation, SizingReplay, FaceBuild, Rehydration,
    SetAudit, Commit
};

struct SurfaceFaceCacheDiagnostic {
    SurfaceFaceCacheStatus status = SurfaceFaceCacheStatus::None;
    SurfaceFaceCacheReason reason = SurfaceFaceCacheReason::None;
    SurfaceFaceCachePhase phase = SurfaceFaceCachePhase::None;
    bool reportComplete = false;
    FacePersistentId failedFace = 0;
    std::uint64_t processedFaces = 0;
    std::uint64_t unattemptedFaces = 0;
    std::uint64_t workUnits = 0;
    std::size_t accountedBytes = 0;
};

struct SurfaceFaceCacheMetrics {
    std::uint64_t keyRowsVisited = 0;
    std::uint64_t dirtyPropagationVisits = 0;
    std::uint64_t facesBuilt = 0;
    std::uint64_t facesReused = 0;
    std::uint64_t facesReplayed = 0;
    std::uint64_t facesRejected = 0;
    std::uint64_t cleanNodeRows = 0;
    std::uint64_t cleanCellRows = 0;
    std::uint64_t dirtyNodeRows = 0;
    std::uint64_t dirtyCellRows = 0;
    std::uint64_t constraintRows = 0;
    std::uint64_t cleanConstraintRows = 0;
    std::uint64_t dirtyConstraintRows = 0;
    std::uint64_t mappingRows = 0;
    std::uint64_t cleanMappingRows = 0;
    std::uint64_t dirtyMappingRows = 0;
    std::uint64_t sizingReplayCalls = 0;
    std::uint64_t preparationWorkUnits = 0;
    std::uint64_t replayWorkUnits = 0;
    std::uint64_t setAuditWorkUnits = 0;
    std::uint64_t cacheCommits = 0;
    std::uint64_t copiedBytes = 0;
    std::uint64_t highWaterSlots = 0;
    std::uint64_t highWaterRows = 0;
    std::uint64_t allFaceScansPerDirtyFace = 0;
    std::size_t accountedByteHighWater = 0;
};

// Caller-owned, fixed-capacity v4-private report channel.  The ordinary v3
// path carries a null pointer and observes no cache branch or storage.
struct SurfaceCacheDiagnosticChannel {
    std::span<SurfaceCacheFaceReportRow> processedStorage{};
    std::span<FacePersistentId> unattemptedStorage{};
    std::size_t processedUsed = 0;
    std::size_t unattemptedUsed = 0;
    bool reportComplete = false;
    bool hasFailedFace = false;
    FacePersistentId failedFacePersistentId = 0;
    SurfaceFaceCacheStatus status = SurfaceFaceCacheStatus::None;
    SurfaceFaceCacheReason reason = SurfaceFaceCacheReason::None;
    SurfaceFaceCachePhase phase = SurfaceFaceCachePhase::None;
    std::uint64_t cacheWorkUnits = 0;
    std::size_t cacheAccountedBytes = 0;

    void reset() noexcept {
        processedUsed=0;
        unattemptedUsed=0;
        reportComplete=false;
        hasFailedFace=false;
        failedFacePersistentId=0;
        status=SurfaceFaceCacheStatus::None;
        reason=SurfaceFaceCacheReason::None;
        phase=SurfaceFaceCachePhase::None;
        cacheWorkUnits=0;
        cacheAccountedBytes=0;
    }
};

class SurfaceFaceCacheSession;

struct DyadicEdgeSplit {
    int edgeTag = -1;
    std::uint64_t baseSegmentOrdinal = 0;
    double parameter = 0.0;
    Vec3 coordinate{};
    SurfaceNodeIdentity identity{};
};

struct PrivateFacePoint {
    FacePersistentId facePersistentId = 0;
    Vec2 parameter{};
    Vec3 coordinate{};
    SurfaceNodeIdentity identity{};
};

struct RefinedSurfaceSegment {
    std::array<Index,2> globalNodes{kInvalidIndex,kInvalidIndex};
    int edgeTag = -1;
    ArtifactStamp owningCurve{};
    std::uint64_t baseSegmentOrdinal = 0;
    std::uint64_t dyadicPath = 0;
    std::uint32_t dyadicLevel = 0;
    double parameterLo = 0.0;
    double parameterHi = 0.0;
};

struct PrivateCurveNodeIdentity {
    Index globalNode = kInvalidIndex;
    SurfaceNodeIdentity identity{};
};

struct PrivateBuildTrace {
    std::vector<RefinedSurfaceSegment> segments;
    std::vector<PrivateCurveNodeIdentity> curveNodeIdentities;
    std::uint64_t surfaceNodes = 0;
    std::uint64_t surfaceSegments = 0;
    std::uint64_t candidateTetCells = 0;
    std::uint64_t constructionWorkUnits = 0;
    std::uint64_t recoveryWorkUnits = 0;
    bool workComplete = true;
    bnd::TetFacetSetFrontierSink* frontierSink = nullptr;
    bnd::TetFacetConflictCycleSink* facetConflictCycleSink = nullptr;
    bnd::TetFacetSetMixedPresentBatchSink* mixedPresentBatchSink = nullptr;
    SurfaceFaceCacheSession* surfaceCacheSession = nullptr;
    SurfaceCacheDiagnosticChannel* surfaceCacheChannel = nullptr;
    std::span<const DyadicEdgeSplit> proposedEdgeSplits{};
    std::span<const PrivateFacePoint> proposedFacePoints{};
    SurfaceFaceCacheDiagnostic surfaceCacheDiagnostic{};
    SurfaceFaceCacheMetrics surfaceCacheMetrics{};
    bool preferDirectConstraintInitialization = false;
#if defined(PAACS_TOPO_TESTING)
    std::array<std::uint64_t,BoundaryConformingV4Metrics::phaseCount>
        phaseNanoseconds{};
#endif
};

void meshEdgesWithProtectedSplits(
    const geo::Model& model, const sizing::MeshSpec& spec,
    std::span<const DyadicEdgeSplit> splits, mesh::Mesh& output,
    PrivateBuildTrace& trace,
    std::span<const DyadicEdgeSplit> proposedSplits = {});

BoundaryConformingV3Outcome meshModelBoundaryConformingWithProtectedSplits(
    const geo::Model& model, const sizing::MeshSpec& spec,
    const BoundaryConformingV3Request& request, mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& provider,
    std::span<const DyadicEdgeSplit> splits,
    std::span<const PrivateFacePoint> facePoints,
    PrivateBuildTrace* trace) noexcept;

#if defined(PAACS_TOPO_TESTING)
// Fixed-value seam for the real-route NoProgress discriminator.  It does not
// mutate the committed semantic ledger and is inactive unless a focused test
// installs one exact structured identity on the current thread.
void setFacetConflictRepeatedIdentityForTesting(
    const SurfaceNodeIdentity& identity) noexcept;
void resetFacetConflictRepeatedIdentityForTesting() noexcept;
#endif

}  // namespace t2c_detail

static_assert(std::is_standard_layout_v<ProtectedFeatureProviderStamp>);
static_assert(std::is_trivially_copyable_v<ProtectedFeatureProviderStamp>);
static_assert(std::is_standard_layout_v<ProtectedFeatureRefinementLimits>);
static_assert(std::is_trivially_copyable_v<ProtectedFeatureRefinementLimits>);
static_assert(std::is_standard_layout_v<SurfaceCacheFaceReportRow>);
static_assert(std::is_trivially_copyable_v<SurfaceCacheFaceReportRow>);
static_assert(std::is_nothrow_move_constructible_v<BoundaryConformingV4Outcome>);

}  // namespace paacs::algo
