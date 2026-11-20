// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

// T2-C v4-private immutable surface-face cache.  This header is deliberately
// kept beside the implementation: no accepted v3/public surface entry sees a
// cache type.  Cache entries contain value rows only; every attempt rehydrates
// fresh handles and a fresh local workspace.

#include "paacs/algo/ProductionSurfaceMesher.hpp"
#include "paacs/algo/ProtectedFeatureRefinement.hpp"
#include "paacs/algo/SurfaceMesherAdapter.hpp"

#include <cstddef>
#include <cstdint>
#include <array>
#include <map>
#include <memory>
#include <span>
#include <vector>

namespace paacs::algo::t2c_detail {

#if defined(PAACS_TOPO_TESTING)
// Allocation qualification is deliberately scoped to the cache-specific
// production work.  The standalone probe replaces global operator new and
// consults this marker; ordinary builds compile no marker or seam.  Nested
// scopes restore their caller's zone, so (for example) local/global mapping
// allocations remain distinguishable from their enclosing rehydration.
enum class SurfaceCacheAllocationZone : std::uint8_t {
    None,
    DiagnosticEnvelope,
    AuthorityKey,
    DirtySet,
    FaceSnapshot,
    SizingReplay,
    Rehydration,
    LocalGlobalMapping,
    SetAudit,
    ReplacementPreparation,
    ProposedStateIntegration,
    FacetConflictStorage,
    FacetConflictPreparedStateStorage,
    FacetConflictPreparation,
    FacetConflictProviderBatch,
    CommitTail,
    Count
};

SurfaceCacheAllocationZone surfaceCacheAllocationZoneForTesting() noexcept;

class SurfaceCacheAllocationZoneScope {
public:
    explicit SurfaceCacheAllocationZoneScope(
        SurfaceCacheAllocationZone zone) noexcept;
    ~SurfaceCacheAllocationZoneScope() noexcept;
    SurfaceCacheAllocationZoneScope(const SurfaceCacheAllocationZoneScope&) =
        delete;
    SurfaceCacheAllocationZoneScope& operator=(
        const SurfaceCacheAllocationZoneScope&) = delete;
private:
    SurfaceCacheAllocationZone prior_ = SurfaceCacheAllocationZone::None;
};

enum class SurfaceCacheProviderCallKind : std::uint8_t {
    CandidateStamp,
    CandidateBuild
};
struct SurfaceCacheProviderCallRow {
    SurfaceCacheProviderCallKind kind =
        SurfaceCacheProviderCallKind::CandidateStamp;
    FacePersistentId face = 0;
    friend bool operator==(const SurfaceCacheProviderCallRow&,
                           const SurfaceCacheProviderCallRow&) = default;
};
void resetSurfaceCacheProviderCallsForTesting() noexcept;
void noteSurfaceCacheProviderCallForTesting(
    SurfaceCacheProviderCallKind kind,FacePersistentId face=0) noexcept;
std::span<const SurfaceCacheProviderCallRow>
surfaceCacheProviderCallsForTesting() noexcept;
bool surfaceCacheProviderCallsCompleteForTesting() noexcept;
#endif

struct SurfaceSizingSampleRow {
    FacePersistentId face = 0;
    std::uint64_t ordinal = 0;
    Vec3 point{};
    double curvature = 0.0;
    double result = 0.0;
    std::uint64_t globalOrdinal = 0;
    // Only successful checked evaluations can enter an immutable cache
    // entry.  Retaining the existing production taxonomy nevertheless makes
    // that fact explicit and prevents future result/status decoupling.
    ProductionSurfaceReason status = ProductionSurfaceReason::None;
};

struct CachedSteinerResult {
    bnd::SteinerStatus status = bnd::SteinerStatus::Failed;
    bnd::SteinerReject reject = bnd::SteinerReject::None;
    bnd::SteinerRefuse refuse = bnd::SteinerRefuse::None;
    bnd::SteinerFail fail = bnd::SteinerFail::None;
    std::vector<std::size_t> recoveredChain;
    std::vector<std::pair<std::size_t,bnd::AuthoritativeParam>> createdParams;
    std::uint64_t revision = 0;
};

struct CachedSteinerEntityResult {
    bnd::ConstraintEntityId entity = 0;
    CachedSteinerResult result{};
};

struct CachedSteinerReport {
    std::vector<CachedSteinerEntityResult> processed;
    std::vector<bnd::ConstraintEntityId> unattempted;
};

struct CachedPersistentConstraintEdge {
    ArtifactStamp faceStamp{};
    bnd::ConstraintEntityId bndEntityId = 0;
    SurfaceChainId chainPersistentId = 0;
    std::size_t originalSegmentOrdinal = 0;
    std::size_t normalizedSegmentOrdinal = 0;
    SurfaceChainRole role = SurfaceChainRole::EmbeddedConstraint;
    std::array<std::size_t,2> localEndpoints{};
    std::array<std::optional<SharedCurveNodeKey>,2> endpointSharedKeys{};
    mesh::Classification classification{};
    bool required = true;
};

struct CachedSurfaceFaceArtifact {
    FacePersistentId facePersistentId = 0;
    // Authority that is independent of the current discretization.  A
    // mismatch is never a dirty-cache event: it invalidates the complete
    // cache authority before any face can be rebuilt.
    std::vector<std::uint8_t> invariantAuthorityBytes;
    // Complete structural authority, including the current face chains,
    // refined-node rows and retained shared bindings.
    std::vector<std::uint8_t> authorityBytes;
    std::vector<SurfaceSizingSampleRow> sizingSamples;

    ArtifactStamp faceStamp{};
    bnd::ChartId chart = 0;
    mesh::Classification faceClass{};
    int desiredNormalSign = 1;
    SurfacePhysicalPolicy physicalPolicy{};
    bnd::GeometryProviderStamp expectedGeometryStamp{};
    ChartCertificate chartCertificate{};
    std::vector<topo::NodePool::NodeData> localNodes;
    std::vector<bnd::AuthoritativeParam> localParams;
    std::vector<Vec2> localMappingParams;
    std::map<std::size_t,SurfaceNodeIdentity> refinedIdentityByLocalNode;
    std::vector<SurfaceFaceCell> localCells;
    std::vector<std::pair<std::size_t,SharedCurveNodeKey>>
        sharedKeysByLocalNode;
    std::vector<std::vector<std::size_t>> boundaryLoops;
    std::vector<SurfaceChain> normalizedChains;
    std::map<bnd::ConstraintEntityId,SurfaceSegmentOrigin> segmentOrigins;
    std::vector<CachedPersistentConstraintEdge> persistentEdges;
    CachedSteinerReport protectedRecovery;
    CachedSteinerReport localRecovery;
    SurfaceTopologyAudit audit{};
    ProductionSurfaceFaceMetrics metrics{};
    bool complete = false;
    std::size_t accountedCapacityBytes = 0;
};

// V4-private, pre-placement capacity contract for a prepared facet-conflict
// batch.  The coordinator owns the value rows themselves; the cache proves
// that their bounded overlay, the fixed face report, the exact dirty-face
// set, and one complete replacement-slot vector all fit before the first
// facet placement provider call.  This is a proposal only: it never changes
// a live cache slot or any split/ledger row.
struct SurfaceFacetCacheProposalCapacities {
    std::size_t committedFacePointRows = 0;
    std::size_t committedFacePointCapacity = 0;
    std::size_t preparedFacePointRows = 0;
    std::size_t preparedFacePointCapacity = 0;
    std::size_t committedLedgerRows = 0;
    std::size_t committedLedgerCapacity = 0;
    std::size_t preparedLedgerRows = 0;
    std::size_t preparedLedgerCapacity = 0;
    std::size_t committedSemanticStateRows = 0;
    std::size_t committedSemanticStateCapacity = 0;
    std::size_t preparedSemanticStateRows = 0;
    std::size_t preparedSemanticStateCapacity = 0;
};

struct SurfaceFacetCacheProposal {
    bool ready = false;
    std::size_t plannedRows = 0;
    std::size_t dirtyFaceRows = 0;
    std::size_t replacementSlotRows = 0;
    std::size_t replacementSlotCapacity = 0;
    std::uint64_t workUnits = 0;
    std::size_t accountedBytes = 0;
};

class SurfaceCacheSizingEvaluator {
public:
    virtual ~SurfaceCacheSizingEvaluator() = default;
    virtual double evaluate(const Vec3& point,double curvature) const = 0;
};

class SurfaceCacheMetricsSource {
public:
    virtual ~SurfaceCacheMetricsSource() = default;
    virtual bool metricForFace(FacePersistentId face,
                               ProductionSurfaceFaceMetrics& out) const = 0;
};

class SurfaceFaceCacheSession {
public:
    explicit SurfaceFaceCacheSession(SurfaceFaceCacheLimits limits = {});

    const SurfaceFaceCacheLimits& limits() const noexcept { return limits_; }
    const SurfaceFaceCacheMetrics& metrics() const noexcept { return metrics_; }
    const SurfaceFaceCacheDiagnostic& diagnostic() const noexcept {
        return diagnostic_;
    }
    std::span<const ProductionSurfaceFaceMetrics> attemptMetrics() const
        noexcept { return attemptMetrics_; }
    std::uint64_t attemptWorkUnits() const noexcept;

    void attachDiagnosticChannel(SurfaceCacheDiagnosticChannel* channel)
        noexcept { diagnosticChannel_=channel; }
    SurfaceCacheDiagnosticChannel* diagnosticChannel() noexcept {
        return diagnosticChannel_;
    }

    void beginAttempt() noexcept;
    bool recordSizing(FacePersistentId face,const Vec3& point,
                      double curvature,double result) noexcept;
    std::span<const SurfaceSizingSampleRow> observedSizing(
        FacePersistentId face) const noexcept;
    std::span<const SurfaceSizingSampleRow> observedSizingSequence() const
        noexcept { return globalObservedSizing_; }
    bool verifyAndReplaySizing(const CachedSurfaceFaceArtifact& entry,
                               const SurfaceCacheSizingEvaluator& evaluator)
        noexcept;
    bool verifyFinalSizing(
        std::span<const std::shared_ptr<const CachedSurfaceFaceArtifact>>
            entries,
        const SurfaceCacheSizingEvaluator& evaluator) noexcept;

    std::shared_ptr<const CachedSurfaceFaceArtifact> find(
        FacePersistentId face,std::span<const std::uint8_t> authority)
        noexcept;
    bool configureExpectedDirty(
        std::span<const FacePersistentId> canonicalFaces,
        std::span<const FacePersistentId> expectedDirty,
        std::uint64_t dependencyVisits=0) noexcept;
    bool empty() const noexcept { return slots_.empty(); }
    bool noteLookup(FacePersistentId face,bool dirty) noexcept;
    bool prepareReplacement(
        std::vector<std::shared_ptr<const CachedSurfaceFaceArtifact>> next,
        std::span<const FacePersistentId> canonicalFaces) noexcept;

    // `plannedFacets` is the canonical SCC member vector.  All allocation
    // performed here is proposal scratch and therefore occurs
    // before a protected-feature placement call.  Success also primes the
    // exact slot-owner vector consumed by the next cached source rebuild.
    bool prepareFacetProposal(
        std::span<const FacePersistentId> canonicalFaces,
        std::span<const bnd::VolumeFacetKey> plannedFacets,
        const SurfaceFacetCacheProposalCapacities& capacities,
        SurfaceFacetCacheProposal& proposal) noexcept;
    std::span<const FacePersistentId> preparedFacetDirtyFaces() const
        noexcept { return facetProposalDirtyFaces_; }
    void abandonFacetProposal() noexcept;

    bool acquireReplacementScratch(
        std::size_t required,
        std::vector<std::shared_ptr<const CachedSurfaceFaceArtifact>>& out)
        noexcept;
    void recycleReplacementScratch(
        std::vector<std::shared_ptr<const CachedSurfaceFaceArtifact>>&& rows)
        noexcept;
    void commitPrepared() noexcept;
    void abandonPrepared() noexcept;

    void setDiagnostic(SurfaceFaceCacheStatus status,
                       SurfaceFaceCacheReason reason,
                       SurfaceFaceCachePhase phase,bool complete,
                       FacePersistentId failedFace,
                       std::uint64_t processed,
                       std::uint64_t unattempted) noexcept;
    bool noteBuilt(const CachedSurfaceFaceArtifact& entry) noexcept;
    bool noteReused(const CachedSurfaceFaceArtifact& entry) noexcept;
    bool noteSetAudit(std::uint64_t work) noexcept;

#if defined(PAACS_TOPO_TESTING)
    std::span<const std::shared_ptr<const CachedSurfaceFaceArtifact>>
    slotsForTesting() const noexcept { return slots_; }
    std::span<const FacePersistentId> canonicalFacesForTesting() const
        noexcept { return canonicalFaces_; }
#endif

private:
    bool consumePreparationWork(std::uint64_t amount,
                                SurfaceFaceCachePhase phase,
                                FacePersistentId face) noexcept;
    bool consumeReplayWork(std::uint64_t amount,
                           SurfaceFaceCachePhase phase,
                           FacePersistentId face) noexcept;
    bool consumeWork(std::uint64_t& bucket,std::uint64_t amount,
                     SurfaceFaceCachePhase phase,
                     FacePersistentId face) noexcept;
    SurfaceFaceCacheLimits limits_{};
    SurfaceFaceCacheMetrics metrics_{};
    SurfaceFaceCacheDiagnostic diagnostic_{};
    std::map<FacePersistentId,std::vector<SurfaceSizingSampleRow>> observed_;
    std::vector<SurfaceSizingSampleRow> globalObservedSizing_;
    std::vector<SurfaceSizingSampleRow> canonicalSizingSequence_;
    std::vector<SurfaceSizingSampleRow> preparedSizingSequence_;
    std::vector<ProductionSurfaceFaceMetrics> attemptMetrics_;
    std::uint64_t attemptReplayBaseline_ = 0;
    std::uint64_t attemptAuditBaseline_ = 0;
    std::uint64_t attemptPreparationBaseline_ = 0;
    std::uint64_t attemptDirtyFaces_ = 0;
    std::size_t attemptLookupOrdinal_ = 0;
    std::vector<FacePersistentId> attemptFaces_;
    std::vector<FacePersistentId> expectedDirtyFaces_;
    std::vector<FacePersistentId> canonicalFaces_;
    std::vector<std::shared_ptr<const CachedSurfaceFaceArtifact>> slots_;
    std::vector<FacePersistentId> preparedFaces_;
    std::vector<std::shared_ptr<const CachedSurfaceFaceArtifact>> prepared_;
    std::vector<FacePersistentId> facetProposalDirtyFaces_;
    std::vector<std::shared_ptr<const CachedSurfaceFaceArtifact>>
        replacementScratch_;
    std::uint64_t preparedHighWaterSlots_ = 0;
    std::uint64_t preparedHighWaterRows_ = 0;
    std::size_t preparedAccountedBytes_ = 0;
    bool preparedReady_ = false;
    SurfaceCacheDiagnosticChannel* diagnosticChannel_ = nullptr;
};

std::vector<std::uint8_t> encodeSurfaceFaceCacheAuthority(
    ArtifactStamp setStamp,const SurfaceFaceRequest& request,
    const ChartCertificate& certificate,
    const std::map<SharedCurveNodeKey,SharedCurveNodeBinding>& ledger,
    const mesh::Mesh& destination,
    const Vec2& parameterLo,const Vec2& parameterHi,
    const Vec2& metricScale);

std::shared_ptr<const CachedSurfaceFaceArtifact> snapshotSurfaceFaceForCache(
    const SurfaceFaceResult& face,std::span<const std::uint8_t> authority,
    std::span<const SurfaceSizingSampleRow> sizing,
    const ProductionSurfaceFaceMetrics& metrics,
    const std::map<SharedCurveNodeKey,SharedCurveNodeBinding>& ledger);

std::unique_ptr<SurfaceFaceResult> rehydrateSurfaceFaceFromCache(
    const CachedSurfaceFaceArtifact& entry,
    const std::map<SharedCurveNodeKey,SharedCurveNodeBinding>& ledger,
    const mesh::Mesh& destination);

bool auditRehydratedSurfaceFace(const SurfaceFaceResult& face,
                                const mesh::Mesh& destination) noexcept;

SurfaceSetOutcome buildAndPublishSurfaceSetCached(
    const SurfaceSetRequest& request,mesh::Mesh& destination,
    const SurfaceGeometryProvider& geometry,
    const SurfaceCandidateProvider& candidateProvider,
    const bnd::SteinerPolicy& localSteinerPolicy,
    const std::map<FacePersistentId,std::vector<std::uint8_t>>& authorities,
    SurfaceFaceCacheSession& cache,
    const SurfaceCacheSizingEvaluator& sizingEvaluator,
    const SurfaceCacheMetricsSource& metricsSource) noexcept;

} // namespace paacs::algo::t2c_detail
