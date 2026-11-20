// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Focused T2-C immutable face-cache replacement allocation qualification.
// This standalone probe replaces global operator new and therefore is not
// linked into the sanitizer or umbrella test binaries.
#include "../src/algo/T2CSurfaceCache.hpp"
#include "../src/algo/ProductionSurfaceMesherDetail.hpp"
#include "paacs/algo/Meshers.hpp"
#include "paacs/algo/ProtectedFeatureRefinement.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sizing/SizeField.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <new>
#include <set>
#include <span>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <variant>
#include <vector>

namespace cache=paacs::algo::t2c_detail;
namespace algo=paacs::algo;

namespace {
long gOrdinal=0;
long gFailOrdinal=0;
long gTailAllocations=0;
bool gTrackPreparation=false;
bool gTrackProduction=false;
bool gTrackTail=false;
bool gInjected=false;
cache::SurfaceCacheAllocationZone gInjectedZone=
    cache::SurfaceCacheAllocationZone::None;
constexpr std::size_t kZoneCount=static_cast<std::size_t>(
    cache::SurfaceCacheAllocationZone::Count);
std::array<long,kZoneCount> gZoneAllocations{};
std::array<bool,kZoneCount> gSelectedZones{};
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t bytes) {
    const auto zone=cache::surfaceCacheAllocationZoneForTesting();
    const auto zoneIndex=static_cast<std::size_t>(zone);
    if(gTrackTail||zone==cache::SurfaceCacheAllocationZone::CommitTail)
        ++gTailAllocations;
    const bool productionSelected=gTrackProduction&&
        zoneIndex<gSelectedZones.size()&&gSelectedZones[zoneIndex];
    if(gTrackPreparation||productionSelected) {
        ++gOrdinal;
        if(productionSelected)++gZoneAllocations[zoneIndex];
        if(gFailOrdinal>0&&gOrdinal==gFailOrdinal) {
            gInjected=true;
            gInjectedZone=zone;
            throw std::bad_alloc{};
        }
    }
    if(void* memory=std::malloc(bytes?bytes:1))return memory;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete(void* memory,std::size_t) noexcept { std::free(memory); }
void operator delete[](void* memory,std::size_t) noexcept {
    std::free(memory);
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
int failures=0;
#define REQUIRE(condition) do { \
    if(!(condition)) { \
        std::printf("FAIL %d: %s\n",__LINE__,#condition); \
        ++failures; \
    } \
} while(false)

using Entry=cache::CachedSurfaceFaceArtifact;
using Owner=std::shared_ptr<const Entry>;

// Logical payload digest for the immutable cache owners.  Pointer identity
// alone proves ownership stability, but not that every retained row is the
// same row.  This digest covers every CachedSurfaceFaceArtifact field without
// hashing padding, allocation addresses, capacities, or container internals.
class LogicalDigest {
public:
    void scalar(std::uint64_t value) noexcept {
        for(unsigned byte=0;byte<8;++byte) {
            value_ ^= (value>>(byte*8U))&UINT64_C(0xff);
            value_ *= UINT64_C(1099511628211);
        }
    }
    void signedScalar(std::int64_t value) noexcept {
        scalar(static_cast<std::uint64_t>(value));
    }
    void floating(double value) noexcept {
        scalar(std::bit_cast<std::uint64_t>(value));
    }
    void bytes(std::span<const std::uint8_t> values) noexcept {
        scalar(values.size());
        for(const auto value:values)scalar(value);
    }
    void string(std::string_view value) noexcept {
        scalar(value.size());
        for(const unsigned char byte:value)scalar(byte);
    }
    std::uint64_t value() const noexcept { return value_; }
private:
    std::uint64_t value_=UINT64_C(1469598103934665603);
};

void digestClassification(LogicalDigest& digest,
                          const paacs::mesh::Classification& cls) noexcept {
    digest.signedScalar(cls.dim);digest.signedScalar(cls.tag);
}

void digestVec(LogicalDigest& digest,const paacs::Vec2& value) noexcept {
    digest.floating(value.u);digest.floating(value.v);
}

void digestVec(LogicalDigest& digest,const paacs::Vec3& value) noexcept {
    digest.floating(value.x);digest.floating(value.y);
    digest.floating(value.z);
}

void digestStamp(LogicalDigest& digest,const algo::ArtifactStamp& stamp)
        noexcept {
    digest.scalar(stamp.persistentId);digest.scalar(stamp.revision);
}

void digestIdentity(LogicalDigest& digest,
                    const algo::SurfaceNodeIdentity& identity) noexcept {
    digest.scalar(static_cast<std::uint64_t>(identity.domain));
    digestStamp(digest,identity.parentArtifact);
    for(const auto node:identity.parentFacetGlobalNodes)digest.scalar(node);
    digest.scalar(identity.parentOrdinal);
    digest.scalar(identity.parentSubordinal);
    digest.scalar(identity.parameterUBits);
    digest.scalar(identity.parameterVBits);
    digest.scalar(identity.dyadicNumerator);
    digest.scalar(identity.dyadicLevel);
}

void digestSharedKey(LogicalDigest& digest,
                     const algo::SharedCurveNodeKey& key) noexcept {
    digest.scalar(key.curvePersistentId);digest.scalar(key.curveRevision);
    digest.scalar(key.nodeOrdinal);digestIdentity(digest,key.identity);
}

void digestOptionalSharedKey(
        LogicalDigest& digest,
        const std::optional<algo::SharedCurveNodeKey>& key) noexcept {
    digest.scalar(key.has_value()?1U:0U);
    if(key)digestSharedKey(digest,*key);
}

void digestChainVertex(LogicalDigest& digest,
                       const algo::SurfaceChainVertex& vertex) noexcept {
    digest.scalar(vertex.faceLocalVertexId);
    digestOptionalSharedKey(digest,vertex.sharedKey);
}

void digestChain(LogicalDigest& digest,
                 const algo::SurfaceChain& chain) noexcept {
    digest.scalar(chain.chainPersistentId);
    digest.scalar(static_cast<std::uint64_t>(chain.role));
    digestClassification(digest,chain.classification);
    digest.scalar(chain.closed?1U:0U);
    digest.scalar(chain.uses.size());
    for(const auto& use:chain.uses)digestChainVertex(digest,use);
    digest.scalar(chain.segmentEntityIds.size());
    for(const auto id:chain.segmentEntityIds)digest.scalar(id);
    digest.scalar(chain.owningCurve.has_value()?1U:0U);
    if(chain.owningCurve)digestStamp(digest,*chain.owningCurve);
}

void digestParam(LogicalDigest& digest,
                 const paacs::bnd::AuthoritativeParam& param) noexcept {
    digest.scalar(param.index());
    if(const auto* stored=std::get_if<paacs::bnd::StoredParam>(&param))
        digestVec(digest,stored->uv);
    else {
        const auto& affine=std::get<paacs::bnd::AffineParam>(param);
        digest.scalar(affine.lo);digest.scalar(affine.hi);
        digest.scalar(affine.num);digest.scalar(affine.level);
    }
}

void digestCachedSteinerReport(LogicalDigest& digest,
                               const cache::CachedSteinerReport& report)
        noexcept {
    digest.scalar(report.processed.size());
    for(const auto& row:report.processed) {
        digest.scalar(row.entity);
        const auto& result=row.result;
        digest.scalar(static_cast<std::uint64_t>(result.status));
        digest.scalar(static_cast<std::uint64_t>(result.reject));
        digest.scalar(static_cast<std::uint64_t>(result.refuse));
        digest.scalar(static_cast<std::uint64_t>(result.fail));
        digest.scalar(result.recoveredChain.size());
        for(const auto node:result.recoveredChain)digest.scalar(node);
        digest.scalar(result.createdParams.size());
        for(const auto& [node,param]:result.createdParams) {
            digest.scalar(node);digestParam(digest,param);
        }
        digest.scalar(result.revision);
    }
    digest.scalar(report.unattempted.size());
    for(const auto entity:report.unattempted)digest.scalar(entity);
}

std::uint64_t cachedEntryDigest(const Entry& entry) noexcept {
    LogicalDigest digest;
    digest.scalar(entry.facePersistentId);
    digest.bytes(entry.invariantAuthorityBytes);
    digest.bytes(entry.authorityBytes);
    digest.scalar(entry.sizingSamples.size());
    for(const auto& row:entry.sizingSamples) {
        digest.scalar(row.face);digest.scalar(row.ordinal);
        digestVec(digest,row.point);digest.floating(row.curvature);
        digest.floating(row.result);digest.scalar(row.globalOrdinal);
        digest.scalar(static_cast<std::uint64_t>(row.status));
    }
    digestStamp(digest,entry.faceStamp);digest.scalar(entry.chart);
    digestClassification(digest,entry.faceClass);
    digest.signedScalar(entry.desiredNormalSign);
    digest.floating(entry.physicalPolicy.minPhysicalDoubleArea);
    digest.floating(entry.physicalPolicy.maxVertexMappingError);
    digest.floating(entry.physicalPolicy.maxMappingDeviation);
    digest.floating(entry.physicalPolicy.minNormalAlignment);
    digest.floating(entry.physicalPolicy.maxArtifactPositionError);
    digest.scalar(entry.expectedGeometryStamp.identity);
    digest.scalar(entry.expectedGeometryStamp.revision);
    digest.scalar(entry.chartCertificate.providerStamp.identity);
    digest.scalar(entry.chartCertificate.providerStamp.revision);
    digest.scalar(entry.chartCertificate.chart);
    digestStamp(digest,entry.chartCertificate.certificateStamp);
    digest.scalar(entry.chartCertificate.certifiedRegionBoundary.size());
    for(const auto& uv:entry.chartCertificate.certifiedRegionBoundary)
        digestVec(digest,uv);
    digest.signedScalar(entry.chartCertificate.chartPositiveToNormalSign);
    digest.scalar(entry.chartCertificate.regular?1U:0U);
    digest.scalar(entry.chartCertificate.periodic?1U:0U);
    digest.scalar(entry.chartCertificate.seam?1U:0U);
    digest.scalar(entry.localNodes.size());
    for(const auto& node:entry.localNodes) {
        digestVec(digest,node.coord);digestClassification(digest,node.cls);
    }
    digest.scalar(entry.localParams.size());
    for(const auto& param:entry.localParams)digestParam(digest,param);
    digest.scalar(entry.localMappingParams.size());
    for(const auto& uv:entry.localMappingParams)digestVec(digest,uv);
    digest.scalar(entry.refinedIdentityByLocalNode.size());
    for(const auto& [node,identity]:entry.refinedIdentityByLocalNode) {
        digest.scalar(node);digestIdentity(digest,identity);
    }
    digest.scalar(entry.localCells.size());
    for(const auto& cell:entry.localCells) {
        for(const auto node:cell.nodes)digest.scalar(node);
        digestClassification(digest,cell.classification);
    }
    digest.scalar(entry.sharedKeysByLocalNode.size());
    for(const auto& [node,key]:entry.sharedKeysByLocalNode) {
        digest.scalar(node);digestSharedKey(digest,key);
    }
    digest.scalar(entry.boundaryLoops.size());
    for(const auto& loop:entry.boundaryLoops) {
        digest.scalar(loop.size());
        for(const auto node:loop)digest.scalar(node);
    }
    digest.scalar(entry.normalizedChains.size());
    for(const auto& chain:entry.normalizedChains)digestChain(digest,chain);
    digest.scalar(entry.segmentOrigins.size());
    for(const auto& [entity,origin]:entry.segmentOrigins) {
        digest.scalar(entity);digest.scalar(origin.bndEntityId);
        digest.scalar(origin.chainPersistentId);
        digest.scalar(origin.originalSegmentOrdinal);
        digest.scalar(origin.normalizedSegmentOrdinal);
        for(const auto& endpoint:origin.endpoints)
            digestChainVertex(digest,endpoint);
        digestClassification(digest,origin.classification);
    }
    digest.scalar(entry.persistentEdges.size());
    for(const auto& edge:entry.persistentEdges) {
        digestStamp(digest,edge.faceStamp);digest.scalar(edge.bndEntityId);
        digest.scalar(edge.chainPersistentId);
        digest.scalar(edge.originalSegmentOrdinal);
        digest.scalar(edge.normalizedSegmentOrdinal);
        digest.scalar(static_cast<std::uint64_t>(edge.role));
        for(const auto node:edge.localEndpoints)digest.scalar(node);
        for(const auto& key:edge.endpointSharedKeys)
            digestOptionalSharedKey(digest,key);
        digestClassification(digest,edge.classification);
        digest.scalar(edge.required?1U:0U);
    }
    digestCachedSteinerReport(digest,entry.protectedRecovery);
    digestCachedSteinerReport(digest,entry.localRecovery);
    const auto& audit=entry.audit;
    digest.scalar(audit.workspaceReady);digest.scalar(audit.connectedGenusZero);
    digest.scalar(audit.boundaryCyclesExact);digest.scalar(audit.nodeLinksValid);
    digest.scalar(audit.constraintsPersistent);
    digest.scalar(audit.classificationsPreserved);
    digest.scalar(audit.sharedMappingsPreserved);
    digest.scalar(audit.noFrameSurvivors);
    digest.scalar(audit.noDuplicateTriangles);
    digest.scalar(audit.physicalOrientation);
    digest.scalar(audit.boundaryCycleCount);
    digest.signedScalar(audit.eulerCharacteristic);
    const auto& metrics=entry.metrics;
    digest.scalar(metrics.facePersistentId);digest.scalar(metrics.boundaryNodes);
    digest.scalar(metrics.interiorCandidates);digest.scalar(metrics.spacingQueries);
    digest.scalar(metrics.spacingBucketsVisited);
    digest.scalar(metrics.spacingDistanceTests);
    digest.scalar(metrics.spacingIndexHighWater);
    digest.scalar(metrics.constraintPointTests);
    digest.scalar(metrics.protectedSegments);
    digest.scalar(metrics.recoveredProtectedSegments);
    const auto& construction=metrics.construction;
#define DIGEST_CONSTRUCTION_FIELD(field) digest.scalar(construction.field)
    DIGEST_CONSTRUCTION_FIELD(inputPoints);
    DIGEST_CONSTRUCTION_FIELD(coordinateSortComparisons);
    DIGEST_CONSTRUCTION_FIELD(insertionOrderComparisons);
    DIGEST_CONSTRUCTION_FIELD(bootstrapPoints);
    DIGEST_CONSTRUCTION_FIELD(hullOrientationTests);
    DIGEST_CONSTRUCTION_FIELD(bootstrapLegalizationEdgeVisits);
    DIGEST_CONSTRUCTION_FIELD(bootstrapLegalizationOrientationTests);
    DIGEST_CONSTRUCTION_FIELD(bootstrapLegalizationIncircleTests);
    DIGEST_CONSTRUCTION_FIELD(bootstrapCocircularTieTests);
    DIGEST_CONSTRUCTION_FIELD(bootstrapLegalizationFlips);
    DIGEST_CONSTRUCTION_FIELD(bootstrapWorklistInsertions);
    DIGEST_CONSTRUCTION_FIELD(bootstrapWorklistUpdates);
    DIGEST_CONSTRUCTION_FIELD(bootstrapStaleEdgeVisits);
    DIGEST_CONSTRUCTION_FIELD(bootstrapWorklistHighWater);
    DIGEST_CONSTRUCTION_FIELD(pointsInserted);
    DIGEST_CONSTRUCTION_FIELD(locationTriangleTests);
    DIGEST_CONSTRUCTION_FIELD(walkSteps);
    DIGEST_CONSTRUCTION_FIELD(cavityPredicateTests);
    DIGEST_CONSTRUCTION_FIELD(cavityMembers);
    DIGEST_CONSTRUCTION_FIELD(boundaryEdges);
    DIGEST_CONSTRUCTION_FIELD(trianglesCommitted);
    DIGEST_CONSTRUCTION_FIELD(fullScanFallbacks);
    DIGEST_CONSTRUCTION_FIELD(fullScanTriangleTests);
    DIGEST_CONSTRUCTION_FIELD(peakLiveTriangles);
    DIGEST_CONSTRUCTION_FIELD(materializedTriangleSlots);
    DIGEST_CONSTRUCTION_FIELD(peakCavitySize);
    DIGEST_CONSTRUCTION_FIELD(peakBoundarySize);
    DIGEST_CONSTRUCTION_FIELD(scratchHighWaterBytes);
    DIGEST_CONSTRUCTION_FIELD(accountedCapacityHighWaterBytes);
#undef DIGEST_CONSTRUCTION_FIELD
    digest.scalar(construction.adjacencyAuditPassed?1U:0U);
    digest.scalar(construction.runComplete?1U:0U);
    digest.scalar(entry.complete?1U:0U);
    digest.scalar(entry.accountedCapacityBytes);
    return digest.value();
}

Owner makeEntry(algo::FacePersistentId face,std::uint8_t authority) {
    auto entry=std::make_shared<Entry>();
    entry->facePersistentId=face;
    entry->faceStamp={face,1};
    entry->authorityBytes={authority};
    entry->complete=true;
    entry->audit.workspaceReady=true;
    entry->metrics.facePersistentId=face;
    entry->metrics.construction.adjacencyAuditPassed=true;
    entry->metrics.construction.runComplete=true;
    return entry;
}

struct Fixture {
    cache::SurfaceFaceCacheSession session;
    Owner old10=makeEntry(10,1);
    Owner old20=makeEntry(20,2);
    Owner next10=makeEntry(10,3);
    Owner next20=makeEntry(20,4);
    Owner next30=makeEntry(30,5);

    Fixture() {
        const std::array<algo::FacePersistentId,2> faces{10,20};
        std::vector<Owner> initial{old10,old20};
        REQUIRE(session.prepareReplacement(std::move(initial),faces));
        session.commitPrepared();
    }

    std::vector<Owner> replacement() const {
        return {next10,next20,next30};
    }
};

constexpr std::array<algo::FacePersistentId,3> kReplacementFaces{10,20,30};
constexpr std::array<std::uint8_t,1> kOld10{1};
constexpr std::array<std::uint8_t,1> kOld20{2};
constexpr std::array<std::uint8_t,1> kNext10{3};
constexpr std::array<std::uint8_t,1> kNext20{4};
constexpr std::array<std::uint8_t,1> kNext30{5};

void requireOldState(Fixture& fixture,
                     std::uint64_t expectedHighWater=2) {
    REQUIRE(fixture.session.find(10,kOld10)==fixture.old10);
    REQUIRE(fixture.session.find(20,kOld20)==fixture.old20);
    REQUIRE(!fixture.session.find(10,kNext10));
    REQUIRE(!fixture.session.find(20,kNext20));
    REQUIRE(!fixture.session.find(30,kNext30));
    REQUIRE(fixture.session.metrics().cacheCommits==1);
    REQUIRE(fixture.session.metrics().highWaterSlots==expectedHighWater);
}

void requireNewState(Fixture& fixture) {
    REQUIRE(!fixture.session.find(10,kOld10));
    REQUIRE(!fixture.session.find(20,kOld20));
    REQUIRE(fixture.session.find(10,kNext10)==fixture.next10);
    REQUIRE(fixture.session.find(20,kNext20)==fixture.next20);
    REQUIRE(fixture.session.find(30,kNext30)==fixture.next30);
    REQUIRE(fixture.session.metrics().cacheCommits==2);
    // Preparation simultaneously retains the two committed owners and the
    // three completely prepared replacements.  The allocation-free commit
    // then swaps ownership; high-water measures that real overlap window.
    REQUIRE(fixture.session.metrics().highWaterSlots==5);
}

bool prepare(Fixture& fixture,long failOrdinal) {
    auto next=fixture.replacement();
    gOrdinal=0;
    gFailOrdinal=failOrdinal;
    gInjected=false;
    gTrackPreparation=true;
    const bool result=fixture.session.prepareReplacement(
        std::move(next),kReplacementFaces);
    gTrackPreparation=false;
    gFailOrdinal=0;
    return result;
}

void commitWithoutAlloc(Fixture& fixture) {
    const long before=gTailAllocations;
    gTrackTail=true;
    fixture.session.commitPrepared();
    gTrackTail=false;
    REQUIRE(gTailAllocations==before);
}

void abandonWithoutAlloc(Fixture& fixture) {
    const long before=gTailAllocations;
    gTrackTail=true;
    fixture.session.abandonPrepared();
    gTrackTail=false;
    REQUIRE(gTailAllocations==before);
}

long measurePreparation() {
    Fixture fixture;
    requireOldState(fixture);
    REQUIRE(prepare(fixture,0));
    const long total=gOrdinal;
    REQUIRE(total>0);
    commitWithoutAlloc(fixture);
    requireNewState(fixture);
    return total;
}

void sweepPreparation(long total) {
    for(long ordinal=1;ordinal<=total;++ordinal) {
        Fixture fixture;
        requireOldState(fixture);
        const bool prepared=prepare(fixture,ordinal);
        REQUIRE(gInjected);
        REQUIRE(!prepared);
        REQUIRE(fixture.session.diagnostic().status==
                cache::SurfaceFaceCacheStatus::Failed);
        REQUIRE(fixture.session.diagnostic().reason==
                cache::SurfaceFaceCacheReason::AllocationFailed);
        REQUIRE(fixture.session.diagnostic().phase==
                cache::SurfaceFaceCachePhase::Preparation);
        REQUIRE(fixture.session.diagnostic().reportComplete);
        requireOldState(fixture);

        // A disarmed same-fixture retry must publish the exact replacement
        // owners and must use the allocation-free commit tail.
        REQUIRE(prepare(fixture,0));
        commitWithoutAlloc(fixture);
        requireNewState(fixture);
    }
}

void abandonProof() {
    Fixture fixture;
    REQUIRE(prepare(fixture,0));
    abandonWithoutAlloc(fixture);
    // Logical ownership rolled back, while the observed old/new overlap is
    // intentionally retained as a diagnostic high-water mark.
    requireOldState(fixture,5);

    // Abandonment must not poison a later clean preparation/commit.
    REQUIRE(prepare(fixture,0));
    commitWithoutAlloc(fixture);
    requireNewState(fixture);
}

bool sameCapacity(const paacs::mesh::MeshCapacitySnapshot& a,
                  const paacs::mesh::MeshCapacitySnapshot& b) noexcept {
    return a.nodeCoordinateCapacity==b.nodeCoordinateCapacity&&
        a.nodeClassificationCapacity==b.nodeClassificationCapacity&&
        a.cellCapacity==b.cellCapacity&&
        a.connectivityCapacity==b.connectivityCapacity&&
        a.dedupBucketCount==b.dedupBucketCount&&
        a.dedupEntryCount==b.dedupEntryCount&&
        a.dedupIndexCapacity==b.dedupIndexCapacity&&a.complete==b.complete;
}

struct ProductionFixture {
    std::unique_ptr<paacs::geo::Kernel> kernel;
    std::unique_ptr<paacs::geo::Model> model;
    paacs::sizing::MeshSpec spec;
    algo::ArtifactStamp setStamp{22001,1};
    cache::SurfaceFaceCacheSession session;
    std::vector<algo::SurfaceCacheFaceReportRow> processed;
    std::vector<algo::FacePersistentId> unattempted;
    cache::SurfaceCacheDiagnosticChannel channel;
    cache::DyadicEdgeSplit split{};
    std::vector<const Entry*> oldOwners;

    ProductionFixture() {
        kernel=paacs::geo::makeBuiltinKernel();
        model=kernel->makeBox({0,0,0},{1,1,1});
        spec.field=std::make_shared<paacs::sizing::ConstantField>(2.0);
        spec.threads=1;
        processed.resize(model->faces().size());
        unattempted.resize(model->faces().size());
        {
            std::size_t i=0;
            for(const auto* face:model->faces())
                unattempted[i++]=face&&face->tag()>=0
                    ? static_cast<algo::FacePersistentId>(face->tag())+1U:0;
            std::sort(unattempted.begin(),unattempted.end());
        }
        channel={processed,unattempted};
        session.attachDiagnosticChannel(&channel);
        const auto* edge=model->edges().front();
        REQUIRE(edge!=nullptr);
        if(!edge)return;
        algo::SurfaceNodeIdentity identity;
        identity.domain=algo::SurfaceNodeIdentityDomain::RefinedCurve;
        identity.parentArtifact={
            static_cast<std::uint64_t>(edge->tag())+1U,1};
        identity.parentOrdinal=0;
        identity.parameterUBits=std::bit_cast<std::uint64_t>(0.5);
        identity.dyadicNumerator=1;
        identity.dyadicLevel=1;
        split={edge->tag(),0,0.5,edge->pointAt(0.5),identity};

        paacs::mesh::Mesh base;
        algo::meshEdges(*model,spec,base);
        const auto initial=algo::detail::buildProductionSurfaceSetForV3(
            *model,spec,setStamp,base,{}, {},&session);
        REQUIRE(initial.surface.setStatus==algo::SurfaceSetStatus::Published);
        REQUIRE(initial.surface.retainedArtifact!=nullptr);
        REQUIRE(channel.reportComplete);
        for(const auto& owner:session.slotsForTesting()) {
            REQUIRE(owner!=nullptr);
            oldOwners.push_back(owner.get());
        }
    }
};

struct CacheSlotState {
    const Entry* owner=nullptr;
    algo::FacePersistentId canonicalFace=0;
    std::uint64_t payloadDigest=0;
    friend bool operator==(const CacheSlotState&,
                           const CacheSlotState&)=default;
};

struct ProductionPrivateState {
    std::vector<CacheSlotState> cacheSlots;
    std::uint64_t splitDigest=0;
    friend bool operator==(const ProductionPrivateState&,
                           const ProductionPrivateState&)=default;
};

std::uint64_t dyadicSplitDigest(const cache::DyadicEdgeSplit& split) noexcept {
    LogicalDigest digest;
    digest.signedScalar(split.edgeTag);
    digest.scalar(split.baseSegmentOrdinal);
    digest.floating(split.parameter);
    digestVec(digest,split.coordinate);
    digestIdentity(digest,split.identity);
    return digest.value();
}

ProductionPrivateState privateState(const ProductionFixture& fixture) {
    ProductionPrivateState state;
    const auto slots=fixture.session.slotsForTesting();
    const auto faces=fixture.session.canonicalFacesForTesting();
    REQUIRE(slots.size()==faces.size());
    state.cacheSlots.reserve(slots.size());
    for(std::size_t i=0;i<slots.size();++i) {
        REQUIRE(slots[i]!=nullptr);
        if(slots[i])state.cacheSlots.push_back(
            {slots[i].get(),faces[i],cachedEntryDigest(*slots[i])});
    }
    state.splitDigest=dyadicSplitDigest(fixture.split);
    return state;
}

struct ProductionRun {
    algo::detail::ProductionSurfaceSetStage outcome;
    paacs::mesh::Mesh destination;
    std::uint64_t beforeHash=0;
    paacs::mesh::MeshCapacitySnapshot beforeCapacity{};
    std::vector<cache::SurfaceCacheProviderCallRow> providerCalls;
    ProductionPrivateState beforePrivateState;
};

void selectProductionZones(bool proposedOnly=false) {
    gSelectedZones.fill(false);
    if(proposedOnly) {
        gSelectedZones[static_cast<std::size_t>(
            cache::SurfaceCacheAllocationZone::ProposedStateIntegration)]=true;
        return;
    }
    constexpr std::array<cache::SurfaceCacheAllocationZone,9> zones{
        cache::SurfaceCacheAllocationZone::DiagnosticEnvelope,
        cache::SurfaceCacheAllocationZone::AuthorityKey,
        cache::SurfaceCacheAllocationZone::DirtySet,
        cache::SurfaceCacheAllocationZone::FaceSnapshot,
        cache::SurfaceCacheAllocationZone::SizingReplay,
        cache::SurfaceCacheAllocationZone::Rehydration,
        cache::SurfaceCacheAllocationZone::LocalGlobalMapping,
        cache::SurfaceCacheAllocationZone::SetAudit,
        cache::SurfaceCacheAllocationZone::ReplacementPreparation};
    for(const auto zone:zones)
        gSelectedZones[static_cast<std::size_t>(zone)]=true;
}

ProductionRun runProduction(ProductionFixture& fixture,long failOrdinal) {
    ProductionRun run;
    run.beforePrivateState=privateState(fixture);
    cache::PrivateBuildTrace trace;
    const std::span<const cache::DyadicEdgeSplit> split(&fixture.split,1);
    algo::t2c_detail::meshEdgesWithProtectedSplits(
        *fixture.model,fixture.spec,split,run.destination,trace);
    run.beforeHash=paacs::mesh::meshHash(run.destination);
    run.beforeCapacity=run.destination.capacitySnapshot();
    gOrdinal=0;gFailOrdinal=failOrdinal;gInjected=false;
    gInjectedZone=cache::SurfaceCacheAllocationZone::None;
    gZoneAllocations.fill(0);
    struct TrackGuard {
        TrackGuard() noexcept { gTrackProduction=true; }
        ~TrackGuard() noexcept { gTrackProduction=false; }
    } track;
    run.outcome=algo::detail::buildProductionSurfaceSetForV3(
        *fixture.model,fixture.spec,fixture.setStamp,run.destination,
        trace.curveNodeIdentities,{},&fixture.session,split);
    gFailOrdinal=0;
    const auto calls=cache::surfaceCacheProviderCallsForTesting();
    run.providerCalls.assign(calls.begin(),calls.end());
    REQUIRE(cache::surfaceCacheProviderCallsCompleteForTesting());
    return run;
}

bool providerPrefix(
        std::span<const cache::SurfaceCacheProviderCallRow> prefix,
        std::span<const cache::SurfaceCacheProviderCallRow> full) {
    return prefix.size()<=full.size()&&
        std::equal(prefix.begin(),prefix.end(),full.begin());
}

void requireCompleteXor(const ProductionFixture& fixture) {
    REQUIRE(fixture.channel.reportComplete);
    REQUIRE(fixture.channel.processedUsed+
            fixture.channel.unattemptedUsed==fixture.model->faces().size());
    std::set<algo::FacePersistentId> ids;
    for(std::size_t i=0;i<fixture.channel.processedUsed;++i)
        REQUIRE(ids.insert(
            fixture.channel.processedStorage[i].facePersistentId).second);
    for(std::size_t i=0;i<fixture.channel.unattemptedUsed;++i)
        REQUIRE(ids.insert(fixture.channel.unattemptedStorage[i]).second);
    REQUIRE(ids.size()==fixture.model->faces().size());
}

void requireOldProductionOwners(const ProductionFixture& fixture) {
    const auto slots=fixture.session.slotsForTesting();
    REQUIRE(slots.size()==fixture.oldOwners.size());
    for(std::size_t i=0;i<slots.size()&&i<fixture.oldOwners.size();++i)
        REQUIRE(slots[i].get()==fixture.oldOwners[i]);
}

struct ProductionReference {
    std::uint64_t strict=0,semantic=0;
    bool semanticCanonical=false;
    std::vector<cache::SurfaceCacheProviderCallRow> providerCalls;
    std::array<long,kZoneCount> zoneCounts{};
    long total=0;
};

ProductionReference measureProductionRoute() {
    ProductionFixture fixture;
    selectProductionZones();
    auto run=runProduction(fixture,0);
    REQUIRE(run.outcome.surface.setStatus==algo::SurfaceSetStatus::Published);
    REQUIRE(run.outcome.surface.retainedArtifact!=nullptr);
    requireCompleteXor(fixture);
    ProductionReference reference;
    reference.strict=paacs::mesh::meshHash(run.destination);
    const auto semantic=paacs::mesh::semanticHashEx(run.destination);
    reference.semantic=semantic.value;
    reference.semanticCanonical=semantic.canonical;
    reference.providerCalls=std::move(run.providerCalls);
    reference.zoneCounts=gZoneAllocations;
    reference.total=gOrdinal;
    REQUIRE(reference.total>2);
    REQUIRE(gTailAllocations==0);
    return reference;
}

void sweepProductionRoute(const ProductionReference& reference) {
    for(long ordinal=1;ordinal<=reference.total;++ordinal) {
        ProductionFixture fixture;
        selectProductionZones();
        auto failed=runProduction(fixture,ordinal);
        REQUIRE(gInjected);
        REQUIRE(gInjectedZone!=cache::SurfaceCacheAllocationZone::None);
        REQUIRE(failed.outcome.surface.setStatus!=
                algo::SurfaceSetStatus::Published);
        REQUIRE(paacs::mesh::meshHash(failed.destination)==failed.beforeHash);
        // Capacity residue is permitted only after publication preparation.
        if(gInjectedZone!=cache::SurfaceCacheAllocationZone::SetAudit)
            REQUIRE(sameCapacity(failed.destination.capacitySnapshot(),
                                 failed.beforeCapacity));
        requireOldProductionOwners(fixture);
        // Complete logical payload + owner + canonical-slot ordering, and the
        // exact proposed split, are unchanged at the failed-preparation tier.
        REQUIRE(privateState(fixture)==failed.beforePrivateState);
        requireCompleteXor(fixture);
        REQUIRE(providerPrefix(failed.providerCalls,
                               reference.providerCalls));

        selectProductionZones();
        auto retry=runProduction(fixture,0);
        REQUIRE(retry.outcome.surface.setStatus==
                algo::SurfaceSetStatus::Published);
        REQUIRE(paacs::mesh::meshHash(retry.destination)==reference.strict);
        const auto semantic=paacs::mesh::semanticHashEx(retry.destination);
        REQUIRE(semantic.value==reference.semantic);
        REQUIRE(semantic.canonical==reference.semanticCanonical);
        REQUIRE(retry.providerCalls==reference.providerCalls);
        requireCompleteXor(fixture);
    }
}

class BoxOwnership final : public algo::VolumeBoundaryOwnershipProvider {
public:
    algo::VolumeBoundaryProviderStamp current{22002,1};
    std::map<algo::FacePersistentId,std::vector<paacs::bnd::VolumeRegionUse>>
        rows;
    algo::VolumeBoundaryProviderStamp stamp() const noexcept override {
        return current;
    }
    algo::VolumeFaceOwnership ownership(
            algo::FacePersistentId face) const override {
        algo::VolumeFaceOwnership result;
        result.providerStamp=current;result.facePersistentId=face;
        const auto found=rows.find(face);
        if(found!=rows.end())result.owners=found->second;
        return result;
    }
    algo::ProviderRecoveryEnvelope recoveryEnvelope(
            std::span<const paacs::Vec3>) const override {
        algo::ProviderRecoveryEnvelope result;result.providerStamp=current;
        result.availability=
            algo::VolumeProviderAvailability::RecoveryEnvelopeUnavailable;
        return result;
    }
};

class BoxPlacement final : public algo::ProtectedFeaturePlacementProvider {
public:
    algo::ProtectedFeatureProviderStamp current{22003,1};
    mutable std::uint64_t stampCalls=0;
    mutable std::uint64_t edgeCalls=0;
    mutable std::uint64_t faceCalls=0;
    algo::ProtectedFeatureProviderStamp stamp() const noexcept override {
        ++stampCalls;
        return current;
    }
    algo::ProtectedFeaturePlacement placeOnEdge(
            const paacs::geo::GEdge& edge,
            const algo::ProtectedFeaturePlacementRequest& request) const
            override {
        ++edgeCalls;
        return {current,request.entityStamp,request.requestedIdentity,
            algo::ProtectedFeaturePlacementStatus::Placed,
            edge.pointAt(request.curveParameter)};
    }
    algo::ProtectedFeaturePlacement placeOnFace(
            const paacs::geo::GFace& face,
            const algo::ProtectedFeaturePlacementRequest& request) const
            override {
        ++faceCalls;
        return {current,request.entityStamp,request.requestedIdentity,
            algo::ProtectedFeaturePlacementStatus::Placed,
            face.pointAt(request.faceParameter)};
    }
};

std::uint64_t committedV4StateDigest(
        const algo::BoundaryConformingV4Result& result,
        bool hasLastSource) noexcept {
    LogicalDigest digest;
    digest.string(result.routeId?result.routeId:"");
    digest.scalar(result.routeVersion);
    digest.scalar(result.subdivisionLedger.size());
    for(const auto& row:result.subdivisionLedger) {
        digest.scalar(row.parentEntityPersistentId);
        digestIdentity(digest,row.childIdentity);
        digest.scalar(row.splitOrdinal);digest.signedScalar(row.edgeTag);
        digest.signedScalar(row.faceTag);
        digest.floating(row.parentParameterLo);
        digest.floating(row.parentParameterHi);
        digest.floating(row.midpointParameter);
        digestVec(digest,row.faceParameter);digestVec(digest,row.coordinate);
        digestClassification(digest,row.classification);
        digest.scalar(row.providerStamp.persistentId);
        digest.scalar(row.providerStamp.revision);
    }
    // These are the public logical mirrors of committed split, identity and
    // ledger state.  Attempt/report/cache metrics are intentionally excluded.
    digest.scalar(result.metrics.edgeSplits);
    digest.scalar(result.metrics.semanticStates);
    digest.scalar(result.metrics.uniquePlacements);
    digest.scalar(result.sourceOutcome?1U:0U);
    digest.scalar(hasLastSource?1U:0U);
    return digest.value();
}

std::uint64_t emptyCommittedV4StateDigest() noexcept {
    const algo::BoundaryConformingV4Result empty;
    return committedV4StateDigest(empty,false);
}

struct ProposedFixture {
    std::unique_ptr<paacs::geo::Kernel> kernel=paacs::geo::makeBuiltinKernel();
    std::unique_ptr<paacs::geo::Model> model=kernel->makeBox({0,0,0},{1,1,1});
    BoxOwnership ownership;
    BoxPlacement placement;
    paacs::sizing::MeshSpec spec;
    algo::BoundaryConformingV4Request request;
    ProposedFixture() {
        paacs::Vec3 lo{},hi{};model->solids().front()->bbox(lo,hi);
        const auto center=(lo+hi)*0.5;
        for(const auto* face:model->faces()) {
            paacs::Vec2 uvLo{},uvHi{};face->paramBounds(uvLo,uvHi);
            const auto uv=(uvLo+uvHi)*0.5;
            const int sense=face->normalAt(uv).dot(
                face->pointAt(uv)-center)>0?1:-1;
            ownership.rows[static_cast<algo::FacePersistentId>(
                face->tag())+1U]={{0,sense}};
        }
        spec.field=std::make_shared<paacs::sizing::ConstantField>(0.6);
        spec.threads=1;
        request.sourceRequest.expectedSurfaceSetStamp={22004,1};
        request.sourceRequest.expectedOwnershipProviderStamp=ownership.current;
        request.expectedPlacementProviderStamp=placement.current;
        request.limits={16,17,16,100000,UINT64_C(1)<<32,
                        UINT64_C(1)<<40,16,1024,1024};
    }
};

algo::BoundaryConformingV4Outcome runProposed(
        ProposedFixture& fixture,paacs::mesh::Mesh& destination,
        long failOrdinal) {
    gOrdinal=0;gFailOrdinal=failOrdinal;gInjected=false;
    gInjectedZone=cache::SurfaceCacheAllocationZone::None;
    gZoneAllocations.fill(0);selectProductionZones(true);
    struct TrackGuard {
        TrackGuard() noexcept { gTrackProduction=true; }
        ~TrackGuard() noexcept { gTrackProduction=false; }
    } track;
    auto outcome=algo::meshModelBoundaryConformingV4(
        *fixture.model,fixture.spec,fixture.request,destination,
        fixture.ownership,fixture.placement);
    gFailOrdinal=0;
    return outcome;
}

long measureProposedStateIntegration() {
    ProposedFixture fixture;paacs::mesh::Mesh destination;
    auto outcome=runProposed(fixture,destination,0);
    REQUIRE(outcome.published());
    REQUIRE(outcome.retainedResult!=nullptr);
    REQUIRE(gTailAllocations==0);
    return gOrdinal;
}

void sweepProposedStateIntegration(long total) {
    for(long ordinal=1;ordinal<=total;++ordinal) {
        ProposedFixture fixture;paacs::mesh::Mesh destination;
        (void)destination.addNode({-3,-5,-7},{0,77});
        const auto before=paacs::mesh::meshHash(destination);
        const auto beforeCapacity=destination.capacitySnapshot();
        auto failed=runProposed(fixture,destination,ordinal);
        REQUIRE(gInjected);
        REQUIRE(gInjectedZone==cache::SurfaceCacheAllocationZone::
                                  ProposedStateIntegration);
        REQUIRE(!failed.published());
        REQUIRE(failed.status==algo::BoundaryConformingV4Status::Failed);
        REQUIRE(failed.reason==algo::BoundaryConformingV4Reason::
                                   AllocationFailed);
        REQUIRE(failed.destinationUnchanged);
        REQUIRE(paacs::mesh::meshHash(destination)==before);
        REQUIRE(sameCapacity(destination.capacitySnapshot(),beforeCapacity));
        REQUIRE(failed.retainedResult!=nullptr);
        REQUIRE(failed.retainedResult&&
                failed.retainedResult->subdivisionLedger.empty());
        REQUIRE(fixture.placement.stampCalls==0);
        REQUIRE(fixture.placement.edgeCalls==0);
        REQUIRE(fixture.placement.faceCalls==0);
        REQUIRE(failed.lastSourceOutcome==nullptr);
        REQUIRE(failed.retainedResult&&
                failed.retainedResult->sourceOutcome==nullptr);
        REQUIRE(failed.retainedResult&&
                failed.retainedResult->metrics.rounds==0);
        REQUIRE(failed.retainedResult&&
                committedV4StateDigest(*failed.retainedResult,
                    failed.lastSourceOutcome!=nullptr)==
                    emptyCommittedV4StateDigest());
        if(!failed.retainedResult||
           !failed.retainedResult->surfaceCacheReport.reportComplete||
           failed.retainedResult->surfaceCacheReport.processed.size()!=
               fixture.model->faces().size())
            std::printf("  proposed ordinal=%ld retained=%d complete=%d "
                        "processed=%zu phase=%d reason=%d\n",ordinal,
                failed.retainedResult?1:0,
                failed.retainedResult&&
                    failed.retainedResult->surfaceCacheReport.reportComplete,
                failed.retainedResult
                    ? failed.retainedResult->surfaceCacheReport.processed.size()
                    : 0U,static_cast<int>(failed.phase),
                static_cast<int>(failed.reason));
        REQUIRE(failed.retainedResult&&
                failed.retainedResult->surfaceCacheReport.reportComplete);
        REQUIRE(failed.retainedResult&&
                failed.retainedResult->surfaceCacheReport.processed.size()==
                    fixture.model->faces().size());
        ProposedFixture retryFixture;
        auto retry=runProposed(retryFixture,destination,0);
        REQUIRE(retry.published());
    }
}
}

int main() {
    const long total=measurePreparation();
    sweepPreparation(total);
    abandonProof();
    const auto production=measureProductionRoute();
    sweepProductionRoute(production);
    const long proposed=measureProposedStateIntegration();
    sweepProposedStateIntegration(proposed);
    std::printf(
        "probe_t2c_cache helper=%ld production=%ld proposed=%ld "
        "tail_allocations=%ld\n",total,production.total,proposed,
        gTailAllocations);
    for(std::size_t i=1;i<production.zoneCounts.size();++i)
        if(production.zoneCounts[i]!=0)
            std::printf("  cache_zone=%zu allocations=%ld\n",i,
                        production.zoneCounts[i]);
    if(failures!=0) {
        std::printf("probe_t2c_cache FAILED assertions=%d\n",failures);
        return 1;
    }
    std::puts("probe_t2c_cache PASSED");
    return 0;
}
