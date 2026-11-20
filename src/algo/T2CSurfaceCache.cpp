// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "T2CSurfaceCache.hpp"

#include "paacs/bnd/SegmentRecoverySession.hpp"

#include <algorithm>
#include <bit>
#include <cstring>
#include <functional>
#include <limits>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace paacs::algo::t2c_detail {

#if defined(PAACS_TOPO_TESTING)
namespace {
thread_local SurfaceCacheAllocationZone gSurfaceCacheAllocationZone =
    SurfaceCacheAllocationZone::None;
thread_local std::array<SurfaceCacheProviderCallRow,256>
    gSurfaceCacheProviderCalls{};
thread_local std::size_t gSurfaceCacheProviderCallsUsed = 0;
thread_local bool gSurfaceCacheProviderCallsComplete = true;
}

SurfaceCacheAllocationZone surfaceCacheAllocationZoneForTesting() noexcept {
    return gSurfaceCacheAllocationZone;
}

SurfaceCacheAllocationZoneScope::SurfaceCacheAllocationZoneScope(
    SurfaceCacheAllocationZone zone) noexcept
    : prior_(gSurfaceCacheAllocationZone) {
    gSurfaceCacheAllocationZone=zone;
}

SurfaceCacheAllocationZoneScope::~SurfaceCacheAllocationZoneScope() noexcept {
    gSurfaceCacheAllocationZone=prior_;
}

void resetSurfaceCacheProviderCallsForTesting() noexcept {
    gSurfaceCacheProviderCallsUsed=0;
    gSurfaceCacheProviderCallsComplete=true;
}

void noteSurfaceCacheProviderCallForTesting(
        SurfaceCacheProviderCallKind kind,FacePersistentId face) noexcept {
    if(gSurfaceCacheProviderCallsUsed>=gSurfaceCacheProviderCalls.size()) {
        gSurfaceCacheProviderCallsComplete=false;
        return;
    }
    gSurfaceCacheProviderCalls[gSurfaceCacheProviderCallsUsed++]={kind,face};
}

std::span<const SurfaceCacheProviderCallRow>
surfaceCacheProviderCallsForTesting() noexcept {
    return {gSurfaceCacheProviderCalls.data(),gSurfaceCacheProviderCallsUsed};
}

bool surfaceCacheProviderCallsCompleteForTesting() noexcept {
    return gSurfaceCacheProviderCallsComplete;
}
#endif

namespace {

using CacheSlotOwner=std::shared_ptr<const CachedSurfaceFaceArtifact>;
using CacheSlotVector=std::vector<CacheSlotOwner>;
using CacheFaceVector=std::vector<FacePersistentId>;
static_assert(std::is_nothrow_move_constructible_v<CacheSlotOwner>);
static_assert(std::is_nothrow_destructible_v<CacheSlotOwner>);
static_assert(noexcept(std::declval<CacheSlotVector&>().swap(
    std::declval<CacheSlotVector&>())));
static_assert(noexcept(std::declval<CacheSlotVector&>().clear()));
static_assert(noexcept(std::declval<CacheFaceVector&>().swap(
    std::declval<CacheFaceVector&>())));
static_assert(noexcept(std::declval<CacheFaceVector&>().clear()));
static_assert(std::is_trivially_copyable_v<SurfaceSizingSampleRow>);

template<class T>
bool addChecked(std::size_t& value,std::size_t count) noexcept {
    if(count>(std::numeric_limits<std::size_t>::max)()/sizeof(T))return false;
    const auto bytes=count*sizeof(T);
    if(bytes>(std::numeric_limits<std::size_t>::max)()-value)return false;
    value+=bytes;
    return true;
}

bool addRaw(std::size_t& value,std::size_t bytes) noexcept {
    if(bytes>(std::numeric_limits<std::size_t>::max)()-value)return false;
    value+=bytes;
    return true;
}

bool add64(std::uint64_t& value,std::uint64_t amount) noexcept {
    if(amount>(std::numeric_limits<std::uint64_t>::max)()-value)return false;
    value+=amount;
    return true;
}

bool sameSizingSample(const SurfaceSizingSampleRow& a,
                      const SurfaceSizingSampleRow& b) noexcept {
    const auto sameBits=[](double x,double y) noexcept {
        if(x==0.0)x=0.0;
        if(y==0.0)y=0.0;
        return std::bit_cast<std::uint64_t>(x)==
               std::bit_cast<std::uint64_t>(y);
    };
    return a.face==b.face&&a.ordinal==b.ordinal&&
           a.globalOrdinal==b.globalOrdinal&&a.status==b.status&&
           sameBits(a.point.x,b.point.x)&&sameBits(a.point.y,b.point.y)&&
           sameBits(a.point.z,b.point.z)&&
           sameBits(a.curvature,b.curvature)&&sameBits(a.result,b.result);
}

template<class K,class V,class Compare,class Allocator>
bool addMapEstimate(std::size_t& value,
                    const std::map<K,V,Compare,Allocator>& rows) noexcept {
    // Accounted capacity is deliberately conservative rather than a claim
    // about a standard-library node layout.  Four pointer widths cover the
    // three links plus colour/padding on the supported ordered-map layouts.
    constexpr std::size_t nodeWidth=sizeof(std::pair<const K,V>)+
        4*sizeof(void*);
    if(rows.size()>(std::numeric_limits<std::size_t>::max)()/nodeWidth)
        return false;
    return addRaw(value,rows.size()*nodeWidth);
}

template<class T>
void appendScalar(std::vector<std::uint8_t>& out,const T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    const auto* first=reinterpret_cast<const std::uint8_t*>(&value);
    out.insert(out.end(),first,first+sizeof(T));
}

void appendDouble(std::vector<std::uint8_t>& out,double value) {
    if(value==0.0)value=0.0; // canonicalize negative zero
    appendScalar(out,std::bit_cast<std::uint64_t>(value));
}

void appendVec2(std::vector<std::uint8_t>& out,const Vec2& value) {
    appendDouble(out,value.u); appendDouble(out,value.v);
}

void appendVec3(std::vector<std::uint8_t>& out,const Vec3& value) {
    appendDouble(out,value.x); appendDouble(out,value.y);
    appendDouble(out,value.z);
}

void appendClass(std::vector<std::uint8_t>& out,
                 mesh::Classification value) {
    appendScalar(out,value.dim); appendScalar(out,value.tag);
}

void appendStamp(std::vector<std::uint8_t>& out,ArtifactStamp value) {
    appendScalar(out,value.persistentId); appendScalar(out,value.revision);
}

std::optional<std::span<const std::uint8_t>> invariantAuthoritySpan(
    std::span<const std::uint8_t> authority) noexcept {
    constexpr std::size_t header=sizeof(std::uint32_t)+sizeof(std::uint64_t);
    if(authority.size()<header)return std::nullopt;
    std::uint32_t version=0;
    std::uint64_t count=0;
    std::memcpy(&version,authority.data(),sizeof(version));
    std::memcpy(&count,authority.data()+sizeof(version),sizeof(count));
    if(version!=3||count>authority.size()-header||
       count>(std::numeric_limits<std::size_t>::max)())return std::nullopt;
    return authority.subspan(header,static_cast<std::size_t>(count));
}

void appendIdentity(std::vector<std::uint8_t>& out,
                    const SurfaceNodeIdentity& value) {
    appendScalar(out,value.domain); appendStamp(out,value.parentArtifact);
    for(const auto node:value.parentFacetGlobalNodes)appendScalar(out,node);
    appendScalar(out,value.parentOrdinal);
    appendScalar(out,value.parentSubordinal);
    appendScalar(out,value.parameterUBits);
    appendScalar(out,value.parameterVBits);
    appendScalar(out,value.dyadicNumerator);
    appendScalar(out,value.dyadicLevel);
}

void appendSharedKey(std::vector<std::uint8_t>& out,
                     const SharedCurveNodeKey& key) {
    appendScalar(out,key.curvePersistentId);
    appendScalar(out,key.curveRevision);
    appendScalar(out,key.nodeOrdinal);
    appendIdentity(out,key.identity);
}

void appendChain(std::vector<std::uint8_t>& out,const SurfaceChain& chain) {
    appendScalar(out,chain.chainPersistentId); appendScalar(out,chain.role);
    appendClass(out,chain.classification); appendScalar(out,chain.closed);
    appendScalar(out,chain.uses.size());
    for(const auto& use:chain.uses) {
        appendScalar(out,use.faceLocalVertexId);
        appendScalar(out,use.sharedKey.has_value());
        if(use.sharedKey)appendSharedKey(out,*use.sharedKey);
    }
    appendScalar(out,chain.segmentEntityIds.size());
    for(const auto entity:chain.segmentEntityIds)appendScalar(out,entity);
    appendScalar(out,chain.owningCurve.has_value());
    if(chain.owningCurve)appendStamp(out,*chain.owningCurve);
}

std::optional<std::size_t> findLocal(
    const std::map<topo::NodeHandle,std::size_t>& local,
    topo::NodeHandle handle) {
    const auto found=local.find(handle);
    if(found==local.end())return std::nullopt;
    return found->second;
}

CachedSteinerReport cacheReport(
    const bnd::SteinerSegmentsReport& report,
    const std::map<topo::NodeHandle,std::size_t>& local,
    const std::vector<PersistentConstraintEdge>& persistentEdges) {
    CachedSteinerReport out;
    out.processed.reserve(report.processed.size());
    out.unattempted=report.unattempted;
    for(const auto& row:report.processed) {
        CachedSteinerEntityResult stored;
        stored.entity=row.entity;
        stored.result.status=row.result.status;
        stored.result.reject=row.result.reject;
        stored.result.refuse=row.result.refuse;
        stored.result.fail=row.result.fail;
        stored.result.revision=row.result.revision;
        stored.result.recoveredChain.reserve(row.result.recoveredChain.size());
        bool direct=true;
        for(const auto handle:row.result.recoveredChain)
            if(!findLocal(local,handle))direct=false;
        if(direct) for(const auto handle:row.result.recoveredChain) {
            const auto index=findLocal(local,handle);
            if(!index)throw std::logic_error("cache report has foreign handle");
            stored.result.recoveredChain.push_back(*index);
        } else {
            for(const auto& edge:persistentEdges) {
                if(edge.bndEntityId!=row.entity)continue;
                const auto a=findLocal(local,edge.localEndpoints[0]);
                const auto b=findLocal(local,edge.localEndpoints[1]);
                if(!a||!b)throw std::logic_error("cache edge has foreign handle");
                if(stored.result.recoveredChain.empty())
                    stored.result.recoveredChain.push_back(*a);
                else if(stored.result.recoveredChain.back()!=*a)
                    throw std::logic_error("cache chain is not contiguous");
                stored.result.recoveredChain.push_back(*b);
            }
            if(stored.result.recoveredChain.size()!=
               row.result.recoveredChain.size())
                throw std::logic_error("cache report chain reconstruction failed");
        }
        stored.result.createdParams.reserve(row.result.createdParams.size());
        for(const auto& created:row.result.createdParams) {
            const auto index=findLocal(local,created.first);
            if(!index)
                throw std::logic_error(
                    "cache cannot lower a foreign created parameter");
            stored.result.createdParams.emplace_back(*index,created.second);
        }
        out.processed.push_back(std::move(stored));
    }
    return out;
}

bnd::SteinerSegmentsReport restoreReport(
    const CachedSteinerReport& report,
    std::span<const topo::NodeHandle> handles) {
    bnd::SteinerSegmentsReport out;
    out.processed.reserve(report.processed.size());
    out.unattempted=report.unattempted;
    for(const auto& row:report.processed) {
        bnd::SteinerEntityResult restored;
        restored.entity=row.entity;
        restored.result.status=row.result.status;
        restored.result.reject=row.result.reject;
        restored.result.refuse=row.result.refuse;
        restored.result.fail=row.result.fail;
        restored.result.revision=row.result.revision;
        restored.result.recoveredChain.reserve(row.result.recoveredChain.size());
        for(const auto index:row.result.recoveredChain) {
            if(index>=handles.size())throw std::logic_error("invalid cached chain");
            restored.result.recoveredChain.push_back(handles[index]);
        }
        restored.result.createdParams.reserve(row.result.createdParams.size());
        for(const auto& created:row.result.createdParams) {
            if(created.first>=handles.size())
                throw std::logic_error("invalid cached created param");
            restored.result.createdParams.emplace_back(
                handles[created.first],created.second);
        }
        out.processed.push_back(std::move(restored));
    }
    return out;
}

bool sameClass(mesh::Classification a,mesh::Classification b) noexcept {
    return a.dim==b.dim&&a.tag==b.tag;
}

bool reportRows(const CachedSteinerReport& report,
                std::uint64_t& result) noexcept {
    result=0;
    const auto add=[&](std::size_t amount) noexcept {
        const auto value=static_cast<std::uint64_t>(amount);
        if(value>(std::numeric_limits<std::uint64_t>::max)()-result)
            return false;
        result+=value;
        return true;
    };
    if(!add(report.processed.size())||!add(report.unattempted.size()))
        return false;
    for(const auto& row:report.processed) {
        if(!add(row.result.recoveredChain.size())||
           !add(row.result.createdParams.size()))return false;
    }
    return true;
}

struct EntryRowCounts {
    std::uint64_t nodes = 0;
    std::uint64_t cells = 0;
    std::uint64_t constraints = 0;
    std::uint64_t mappings = 0;
    std::uint64_t total = 0;
};

bool countEntryRows(const CachedSurfaceFaceArtifact& entry,
                    EntryRowCounts& out) noexcept {
    out={};
    out.nodes=static_cast<std::uint64_t>(entry.localNodes.size());
    out.cells=static_cast<std::uint64_t>(entry.localCells.size());
    out.constraints=static_cast<std::uint64_t>(entry.persistentEdges.size());
    std::uint64_t protectedRows=0,localRows=0;
    if(!reportRows(entry.protectedRecovery,protectedRows)||
       !reportRows(entry.localRecovery,localRows))return false;
    const auto addMapping=[&](std::size_t count) noexcept {
        return add64(out.mappings,static_cast<std::uint64_t>(count));
    };
    if(!addMapping(entry.invariantAuthorityBytes.size())||
       !addMapping(entry.authorityBytes.size())||
       !addMapping(entry.sizingSamples.size())||
       !addMapping(entry.chartCertificate.certifiedRegionBoundary.size())||
       !addMapping(entry.localParams.size())||
       !addMapping(entry.localMappingParams.size())||
       !addMapping(entry.refinedIdentityByLocalNode.size())||
       !addMapping(entry.sharedKeysByLocalNode.size())||
       !addMapping(entry.boundaryLoops.size())||
       !addMapping(entry.normalizedChains.size())||
       !addMapping(entry.segmentOrigins.size())||
       !add64(out.mappings,protectedRows)||
       !add64(out.mappings,localRows))return false;
    for(const auto& loop:entry.boundaryLoops)
        if(!addMapping(loop.size()))return false;
    for(const auto& chain:entry.normalizedChains)
        if(!addMapping(chain.uses.size())||
           !addMapping(chain.segmentEntityIds.size()))return false;
    out.total=1; // the cache-entry row itself
    return add64(out.total,out.nodes)&&add64(out.total,out.cells)&&
           add64(out.total,out.constraints)&&add64(out.total,out.mappings);
}

// Exact logical payload bytes copied into a newly snapshotted cache entry.
// Allocator/node overhead and unused capacity belong to accountedEntry(), not
// to this traffic counter.
bool copiedEntryBytes(const CachedSurfaceFaceArtifact& entry,
                      std::uint64_t& bytes) noexcept {
    std::size_t value=sizeof(entry);
    if(!addChecked<std::uint8_t>(value,
            entry.invariantAuthorityBytes.size())||
       !addChecked<std::uint8_t>(value,entry.authorityBytes.size())||
       !addChecked<SurfaceSizingSampleRow>(value,entry.sizingSamples.size())||
       !addChecked<Vec2>(value,
            entry.chartCertificate.certifiedRegionBoundary.size())||
       !addChecked<topo::NodePool::NodeData>(value,entry.localNodes.size())||
       !addChecked<bnd::AuthoritativeParam>(value,entry.localParams.size())||
       !addChecked<Vec2>(value,entry.localMappingParams.size())||
       !addChecked<SurfaceFaceCell>(value,entry.localCells.size())||
       !addChecked<std::pair<const std::size_t,SurfaceNodeIdentity>>(
            value,entry.refinedIdentityByLocalNode.size())||
       !addChecked<std::pair<std::size_t,SharedCurveNodeKey>>(
            value,entry.sharedKeysByLocalNode.size())||
       !addChecked<std::vector<std::size_t>>(
            value,entry.boundaryLoops.size())||
       !addChecked<SurfaceChain>(value,entry.normalizedChains.size())||
       !addChecked<std::pair<const bnd::ConstraintEntityId,
                            SurfaceSegmentOrigin>>(
            value,entry.segmentOrigins.size())||
       !addChecked<CachedPersistentConstraintEdge>(
            value,entry.persistentEdges.size())||
       !addChecked<CachedSteinerEntityResult>(
            value,entry.protectedRecovery.processed.size())||
       !addChecked<bnd::ConstraintEntityId>(
            value,entry.protectedRecovery.unattempted.size())||
       !addChecked<CachedSteinerEntityResult>(
            value,entry.localRecovery.processed.size())||
       !addChecked<bnd::ConstraintEntityId>(
            value,entry.localRecovery.unattempted.size()))return false;
    for(const auto& loop:entry.boundaryLoops)
        if(!addChecked<std::size_t>(value,loop.size()))return false;
    for(const auto& chain:entry.normalizedChains)
        if(!addChecked<SurfaceChainVertex>(value,chain.uses.size())||
           !addChecked<bnd::ConstraintEntityId>(
                value,chain.segmentEntityIds.size()))return false;
    const auto addReport=[&](const CachedSteinerReport& report) noexcept {
        for(const auto& row:report.processed)
            if(!addChecked<std::size_t>(
                    value,row.result.recoveredChain.size())||
               !addChecked<std::pair<std::size_t,bnd::AuthoritativeParam>>(
                    value,row.result.createdParams.size()))return false;
        return true;
    };
    if(!addReport(entry.protectedRecovery)||!addReport(entry.localRecovery)||
       value>(std::numeric_limits<std::uint64_t>::max)())return false;
    bytes=static_cast<std::uint64_t>(value);
    return true;
}

bool accountedEntry(const CachedSurfaceFaceArtifact& entry,
                    std::size_t& bytes) noexcept {
    bytes=sizeof(entry);
    if(!addChecked<std::uint8_t>(bytes,
            entry.invariantAuthorityBytes.capacity())||
       !addChecked<std::uint8_t>(bytes,entry.authorityBytes.capacity())||
       !addChecked<SurfaceSizingSampleRow>(
            bytes,entry.sizingSamples.capacity())||
       !addChecked<topo::NodePool::NodeData>(
            bytes,entry.localNodes.capacity())||
       !addChecked<bnd::AuthoritativeParam>(
            bytes,entry.localParams.capacity())||
       !addChecked<Vec2>(bytes,entry.localMappingParams.capacity())||
       !addChecked<SurfaceFaceCell>(bytes,entry.localCells.capacity())||
       !addChecked<Vec2>(
            bytes,entry.chartCertificate.certifiedRegionBoundary.capacity())||
       !addMapEstimate(bytes,entry.refinedIdentityByLocalNode)||
       !addChecked<std::pair<std::size_t,SharedCurveNodeKey>>(
            bytes,entry.sharedKeysByLocalNode.capacity())||
       !addChecked<std::vector<std::size_t>>(
            bytes,entry.boundaryLoops.capacity())||
       !addChecked<SurfaceChain>(bytes,entry.normalizedChains.capacity())||
       !addMapEstimate(bytes,entry.segmentOrigins)||
       !addChecked<CachedPersistentConstraintEdge>(
            bytes,entry.persistentEdges.capacity())||
       !addChecked<CachedSteinerEntityResult>(
            bytes,entry.protectedRecovery.processed.capacity())||
       !addChecked<bnd::ConstraintEntityId>(
            bytes,entry.protectedRecovery.unattempted.capacity())||
       !addChecked<CachedSteinerEntityResult>(
            bytes,entry.localRecovery.processed.capacity())||
       !addChecked<bnd::ConstraintEntityId>(
            bytes,entry.localRecovery.unattempted.capacity()))return false;
    for(const auto& loop:entry.boundaryLoops)
        if(!addChecked<std::size_t>(bytes,loop.capacity()))return false;
    for(const auto& chain:entry.normalizedChains)
        if(!addChecked<SurfaceChainVertex>(bytes,chain.uses.capacity())||
           !addChecked<bnd::ConstraintEntityId>(
               bytes,chain.segmentEntityIds.capacity()))return false;
    const auto addReport=[&](const CachedSteinerReport& report) noexcept {
        for(const auto& row:report.processed)
            if(!addChecked<std::size_t>(
                    bytes,row.result.recoveredChain.capacity())||
               !addChecked<std::pair<std::size_t,bnd::AuthoritativeParam>>(
                    bytes,row.result.createdParams.capacity()))return false;
        return true;
    };
    return addReport(entry.protectedRecovery)&&addReport(entry.localRecovery);
}

} // namespace

SurfaceFaceCacheSession::SurfaceFaceCacheSession(SurfaceFaceCacheLimits limits)
    : limits_(limits) {}

void SurfaceFaceCacheSession::beginAttempt() noexcept {
#if defined(PAACS_TOPO_TESTING)
    resetSurfaceCacheProviderCallsForTesting();
#endif
    attemptReplayBaseline_=metrics_.replayWorkUnits;
    attemptAuditBaseline_=metrics_.setAuditWorkUnits;
    attemptPreparationBaseline_=metrics_.preparationWorkUnits;
    attemptDirtyFaces_=0;
    attemptLookupOrdinal_=0;
    attemptFaces_.clear();
    expectedDirtyFaces_.clear();
    facetProposalDirtyFaces_.clear();
    observed_.clear();
    globalObservedSizing_.clear();
    attemptMetrics_.clear();
    diagnostic_={};
    if(diagnosticChannel_)diagnosticChannel_->reset();
}

std::uint64_t SurfaceFaceCacheSession::attemptWorkUnits() const noexcept {
    if(metrics_.replayWorkUnits<attemptReplayBaseline_||
       metrics_.setAuditWorkUnits<attemptAuditBaseline_||
       metrics_.preparationWorkUnits<attemptPreparationBaseline_)return 0;
    const auto replay=metrics_.replayWorkUnits-attemptReplayBaseline_;
    const auto audit=metrics_.setAuditWorkUnits-attemptAuditBaseline_;
    const auto preparation=metrics_.preparationWorkUnits-
        attemptPreparationBaseline_;
    const auto maximum=(std::numeric_limits<std::uint64_t>::max)();
    if(audit>maximum-replay||preparation>maximum-replay-audit)return 0;
    return replay+audit+preparation;
}

bool SurfaceFaceCacheSession::consumeWork(
    std::uint64_t& bucket,std::uint64_t amount,
    SurfaceFaceCachePhase phase,FacePersistentId face) noexcept {
    const auto maximum=(std::numeric_limits<std::uint64_t>::max)();
    std::uint64_t total=0;
    if(!add64(total,metrics_.preparationWorkUnits)||
       !add64(total,metrics_.replayWorkUnits)||
       !add64(total,metrics_.setAuditWorkUnits)||
       amount>maximum-total||
       total+amount>limits_.maxSurfaceCacheWorkUnits||
       amount>maximum-bucket) {
        setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                      SurfaceFaceCacheReason::BudgetReached,phase,true,
                      face,0,0);
        return false;
    }
    bucket+=amount;
    return true;
}

bool SurfaceFaceCacheSession::consumePreparationWork(
    std::uint64_t amount,SurfaceFaceCachePhase phase,
    FacePersistentId face) noexcept {
    return consumeWork(metrics_.preparationWorkUnits,amount,phase,face);
}

bool SurfaceFaceCacheSession::consumeReplayWork(
    std::uint64_t amount,SurfaceFaceCachePhase phase,
    FacePersistentId face) noexcept {
    return consumeWork(metrics_.replayWorkUnits,amount,phase,face);
}

bool SurfaceFaceCacheSession::recordSizing(FacePersistentId face,
                                           const Vec3& point,
                                           double curvature,
                                           double result) noexcept {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::SizingReplay);
#endif
    if(!consumePreparationWork(1,SurfaceFaceCachePhase::Preparation,face))
        return false;
    try {
        auto& rows=observed_[face];
        SurfaceSizingSampleRow row{face,
            static_cast<std::uint64_t>(rows.size()),point,curvature,result,
            static_cast<std::uint64_t>(globalObservedSizing_.size())};
        if(!canonicalSizingSequence_.empty()&&
           (globalObservedSizing_.size()>=canonicalSizingSequence_.size()||
            !sameSizingSample(row,canonicalSizingSequence_[
                globalObservedSizing_.size()]))) {
            setDiagnostic(SurfaceFaceCacheStatus::Failed,
                SurfaceFaceCacheReason::AuthorityChanged,
                SurfaceFaceCachePhase::SizingReplay,true,face,0,0);
            return false;
        }
        rows.push_back(row);
        try { globalObservedSizing_.push_back(row); }
        catch(...) { rows.pop_back(); throw; }
        return true;
    } catch(...) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::AllocationFailed,
                      SurfaceFaceCachePhase::Preparation,false,face,0,0);
        return false;
    }
}

std::span<const SurfaceSizingSampleRow> SurfaceFaceCacheSession::observedSizing(
    FacePersistentId face) const noexcept {
    const auto found=observed_.find(face);
    if(found==observed_.end())return {};
    return found->second;
}

bool SurfaceFaceCacheSession::verifyAndReplaySizing(
    const CachedSurfaceFaceArtifact& entry,
    const SurfaceCacheSizingEvaluator& evaluator) noexcept {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::SizingReplay);
#endif
    try {
        auto [found,inserted]=observed_.try_emplace(entry.facePersistentId);
        (void)inserted;
        auto& observed=found->second;
        if(observed.size()>entry.sizingSamples.size()) {
            setDiagnostic(SurfaceFaceCacheStatus::Failed,
                          SurfaceFaceCacheReason::AuthorityChanged,
                          SurfaceFaceCachePhase::SizingReplay,true,
                          entry.facePersistentId,0,0);
            return false;
        }
        for(std::size_t i=0;i<observed.size();++i) {
            const auto& actual=observed[i];
            const auto& expected=entry.sizingSamples[i];
            if(actual.face!=entry.facePersistentId||
               expected.face!=entry.facePersistentId||
               actual.ordinal!=i||expected.ordinal!=i||
               !sameSizingSample(actual,expected)) {
                setDiagnostic(SurfaceFaceCacheStatus::Failed,
                              SurfaceFaceCacheReason::AuthorityChanged,
                              SurfaceFaceCachePhase::SizingReplay,true,
                              entry.facePersistentId,0,0);
                return false;
            }
        }
        // Complete all potentially allocating ledger preparation before the
        // first replayed sizing call.  Each successful call can then append
        // its exact immutable row without allocation.
        observed.reserve(entry.sizingSamples.size());
        globalObservedSizing_.reserve(globalObservedSizing_.size()+
            entry.sizingSamples.size()-observed.size());
        const auto prefix=observed.size();
        for(std::size_t i=prefix;i<entry.sizingSamples.size();++i) {
            const auto& expected=entry.sizingSamples[i];
            if(expected.face!=entry.facePersistentId||expected.ordinal!=i||
               expected.globalOrdinal!=globalObservedSizing_.size()) {
                setDiagnostic(SurfaceFaceCacheStatus::Failed,
                              SurfaceFaceCacheReason::Invalid,
                              SurfaceFaceCachePhase::SizingReplay,true,
                              entry.facePersistentId,0,0);
                return false;
            }
            if(metrics_.sizingReplayCalls==
                   (std::numeric_limits<std::uint64_t>::max)()||
               !consumeReplayWork(1,SurfaceFaceCachePhase::SizingReplay,
                                  entry.facePersistentId))
                return false;
            ++metrics_.sizingReplayCalls;
            const double actual=evaluator.evaluate(expected.point,
                                                   expected.curvature);
            SurfaceSizingSampleRow replayed=expected;
            replayed.result=actual;
            if(!sameSizingSample(replayed,expected)) {
                setDiagnostic(SurfaceFaceCacheStatus::Failed,
                              SurfaceFaceCacheReason::AuthorityChanged,
                              SurfaceFaceCachePhase::SizingReplay,true,
                              entry.facePersistentId,0,0);
                return false;
            }
            observed.push_back(expected);
            globalObservedSizing_.push_back(expected);
        }
        return true;
    } catch(const std::bad_alloc&) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::AllocationFailed,
                      SurfaceFaceCachePhase::SizingReplay,true,
                      entry.facePersistentId,0,0);
    } catch(...) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::AuthorityChanged,
                      SurfaceFaceCachePhase::SizingReplay,true,
                      entry.facePersistentId,0,0);
    }
    return false;
}

bool SurfaceFaceCacheSession::verifyFinalSizing(
    std::span<const std::shared_ptr<const CachedSurfaceFaceArtifact>> entries,
    const SurfaceCacheSizingEvaluator& evaluator) noexcept {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::SizingReplay);
#endif
    (void)evaluator; // retained for source compatibility; no provider replay
    try {
        FacePersistentId prior=0;
        bool havePrior=false;
        for(const auto& entry:entries) {
            if(!entry||(havePrior&&entry->facePersistentId<=prior)) {
                setDiagnostic(SurfaceFaceCacheStatus::Failed,
                              SurfaceFaceCacheReason::Invalid,
                              SurfaceFaceCachePhase::SizingReplay,true,
                              0,0,0);
                return false;
            }
            prior=entry->facePersistentId;
            havePrior=true;
            const auto observed=observedSizing(entry->facePersistentId);
            if(observed.size()!=entry->sizingSamples.size()) {
                setDiagnostic(SurfaceFaceCacheStatus::Failed,
                    SurfaceFaceCacheReason::AuthorityChanged,
                    SurfaceFaceCachePhase::SizingReplay,true,
                    entry->facePersistentId,0,0);
                return false;
            }
            for(std::size_t i=0;i<observed.size();++i) {
                if(!consumeReplayWork(1,
                       SurfaceFaceCachePhase::SizingReplay,
                       entry->facePersistentId))return false;
                const auto& actual=observed[i];
                const auto& expected=entry->sizingSamples[i];
                if(actual.face!=entry->facePersistentId||
                   expected.face!=entry->facePersistentId||
                   actual.ordinal!=i||expected.ordinal!=i||
                   !sameSizingSample(actual,expected)) {
                    setDiagnostic(SurfaceFaceCacheStatus::Failed,
                        SurfaceFaceCacheReason::AuthorityChanged,
                        SurfaceFaceCachePhase::SizingReplay,true,
                        entry->facePersistentId,0,0);
                    return false;
                }
            }
        }
        // No sizing row may be hidden under a face absent from the canonical
        // slot vector.  Empty expected ledgers need not create map entries.
        for(const auto& [face,rows]:observed_) {
            if(rows.empty())continue;
            const auto found=std::lower_bound(entries.begin(),entries.end(),
                face,[](const auto& entry,FacePersistentId value) {
                    return entry&&entry->facePersistentId<value;
                });
            if(found==entries.end()||!(*found)||
               (*found)->facePersistentId!=face) {
                setDiagnostic(SurfaceFaceCacheStatus::Failed,
                    SurfaceFaceCacheReason::AuthorityChanged,
                    SurfaceFaceCachePhase::SizingReplay,true,face,0,0);
                return false;
            }
        }
        if(!canonicalSizingSequence_.empty()&&
           (globalObservedSizing_.size()!=canonicalSizingSequence_.size()||
            !std::equal(globalObservedSizing_.begin(),
                        globalObservedSizing_.end(),
                        canonicalSizingSequence_.begin(),
                        sameSizingSample))) {
            setDiagnostic(SurfaceFaceCacheStatus::Failed,
                SurfaceFaceCacheReason::AuthorityChanged,
                SurfaceFaceCachePhase::SizingReplay,true,0,0,0);
            return false;
        }
        return true;
    } catch(const std::bad_alloc&) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::AllocationFailed,
                      SurfaceFaceCachePhase::SizingReplay,true,0,0,0);
    } catch(...) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::AuthorityChanged,
                      SurfaceFaceCachePhase::SizingReplay,true,0,0,0);
    }
    return false;
}

std::shared_ptr<const CachedSurfaceFaceArtifact> SurfaceFaceCacheSession::find(
    FacePersistentId face,std::span<const std::uint8_t> authority)
    noexcept {
    const auto pos=std::lower_bound(canonicalFaces_.begin(),
                                    canonicalFaces_.end(),face);
    if(pos==canonicalFaces_.end()||*pos!=face)return {};
    const auto index=static_cast<std::size_t>(pos-canonicalFaces_.begin());
    if(index>=slots_.size()||!slots_[index])return {};
    const auto incomingInvariant=invariantAuthoritySpan(authority);
    if(!slots_[index]->invariantAuthorityBytes.empty()&&
       (!incomingInvariant||
        slots_[index]->invariantAuthorityBytes.size()!=
            incomingInvariant->size()||
        !std::equal(incomingInvariant->begin(),incomingInvariant->end(),
                    slots_[index]->invariantAuthorityBytes.begin()))) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::AuthorityChanged,
                      SurfaceFaceCachePhase::Preparation,true,face,0,0);
        return {};
    }
    if(slots_[index]->authorityBytes.size()!=authority.size()||
       !std::equal(authority.begin(),authority.end(),
                   slots_[index]->authorityBytes.begin()))return {};
    return slots_[index];
}

bool SurfaceFaceCacheSession::configureExpectedDirty(
    std::span<const FacePersistentId> canonicalFaces,
    std::span<const FacePersistentId> expectedDirty,
    std::uint64_t dependencyVisits) noexcept {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::DirtySet);
#endif
    try {
        if(diagnostic_.status!=SurfaceFaceCacheStatus::None||
           canonicalFaces.size()>limits_.maxFaces||
           expectedDirty.size()>limits_.maxDirtyFacesPerBatch||
           !std::is_sorted(canonicalFaces.begin(),canonicalFaces.end())||
           std::adjacent_find(canonicalFaces.begin(),canonicalFaces.end())!=
               canonicalFaces.end()||
           !std::is_sorted(expectedDirty.begin(),expectedDirty.end())||
           std::adjacent_find(expectedDirty.begin(),expectedDirty.end())!=
               expectedDirty.end()) {
            setDiagnostic(expectedDirty.size()>limits_.maxDirtyFacesPerBatch
                    ? SurfaceFaceCacheStatus::Incomplete
                    : SurfaceFaceCacheStatus::Failed,
                expectedDirty.size()>limits_.maxDirtyFacesPerBatch
                    ? SurfaceFaceCacheReason::BudgetReached
                    : SurfaceFaceCacheReason::Invalid,
                SurfaceFaceCachePhase::Preparation,true,0,0,0);
            return false;
        }
        for(const auto face:expectedDirty)
            if(!std::binary_search(canonicalFaces.begin(),
                                   canonicalFaces.end(),face)) {
                setDiagnostic(SurfaceFaceCacheStatus::Failed,
                    SurfaceFaceCacheReason::Invalid,
                    SurfaceFaceCachePhase::Preparation,true,face,0,0);
                return false;
            }
        attemptFaces_.assign(canonicalFaces.begin(),canonicalFaces.end());
        expectedDirtyFaces_.assign(expectedDirty.begin(),expectedDirty.end());
        attemptLookupOrdinal_=0;
        attemptDirtyFaces_=0;
        if(dependencyVisits==0)
            dependencyVisits=static_cast<std::uint64_t>(canonicalFaces.size());
        if(dependencyVisits>(std::numeric_limits<std::uint64_t>::max)()-
                metrics_.dirtyPropagationVisits||
           !consumePreparationWork(dependencyVisits,
                SurfaceFaceCachePhase::Preparation,0))return false;
        metrics_.dirtyPropagationVisits+=dependencyVisits;
        return true;
    } catch(const std::bad_alloc&) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
            SurfaceFaceCacheReason::AllocationFailed,
            SurfaceFaceCachePhase::Preparation,false,0,0,0);
    } catch(...) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
            SurfaceFaceCacheReason::Invalid,
            SurfaceFaceCachePhase::Preparation,true,0,0,0);
    }
    return false;
}

bool SurfaceFaceCacheSession::noteLookup(FacePersistentId face,
                                         bool dirty) noexcept {
    if(diagnostic_.status!=SurfaceFaceCacheStatus::None||
       attemptLookupOrdinal_>=attemptFaces_.size()||
       attemptFaces_[attemptLookupOrdinal_]!=face) {
        if(diagnostic_.status==SurfaceFaceCacheStatus::None)
            setDiagnostic(SurfaceFaceCacheStatus::Failed,
                SurfaceFaceCacheReason::Invalid,
                SurfaceFaceCachePhase::Preparation,true,face,0,0);
        return false;
    }
    const bool expected=std::binary_search(expectedDirtyFaces_.begin(),
                                            expectedDirtyFaces_.end(),face);
    if(expected!=dirty) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
            SurfaceFaceCacheReason::Invalid,
            SurfaceFaceCachePhase::Preparation,true,face,0,0);
        return false;
    }
    ++attemptLookupOrdinal_;
    if(metrics_.keyRowsVisited==
       (std::numeric_limits<std::uint64_t>::max)()) {
        setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                      SurfaceFaceCacheReason::BudgetReached,
                      SurfaceFaceCachePhase::Preparation,true,face,0,0);
        return false;
    }
    if(!consumePreparationWork(1,SurfaceFaceCachePhase::Preparation,face))
        return false;
    ++metrics_.keyRowsVisited;
    if(!dirty)return true;
    if(metrics_.dirtyPropagationVisits==
           (std::numeric_limits<std::uint64_t>::max)()||
       !consumePreparationWork(1,SurfaceFaceCachePhase::Preparation,face)) {
        return false;
    }
    ++metrics_.dirtyPropagationVisits;
    if(attemptDirtyFaces_==
           (std::numeric_limits<std::uint64_t>::max)()||
       attemptDirtyFaces_+1>limits_.maxDirtyFacesPerBatch) {
        setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                      SurfaceFaceCacheReason::BudgetReached,
                      SurfaceFaceCachePhase::Preparation,true,face,0,0);
        return false;
    }
    ++attemptDirtyFaces_;
    return true;
}

bool SurfaceFaceCacheSession::prepareFacetProposal(
    std::span<const FacePersistentId> canonicalFaces,
    std::span<const bnd::VolumeFacetKey> plannedFacets,
    const SurfaceFacetCacheProposalCapacities& capacities,
    SurfaceFacetCacheProposal& proposal) noexcept {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::FacetConflictPreparation);
#endif
    proposal={};
    facetProposalDirtyFaces_.clear();
    const auto reportState=[&]() noexcept {
        const bool complete=diagnosticChannel_&&
            diagnosticChannel_->reportComplete;
        const auto processed=diagnosticChannel_
            ? static_cast<std::uint64_t>(
                diagnosticChannel_->processedUsed):0U;
        const auto unattempted=diagnosticChannel_
            ? static_cast<std::uint64_t>(
                diagnosticChannel_->unattemptedUsed):0U;
        return std::array<std::uint64_t,3>{
            complete ? 1U:0U,processed,unattempted};
    };
    const auto stop=[&](SurfaceFaceCacheStatus status,
                        SurfaceFaceCacheReason reason) noexcept {
        const auto report=reportState();
        setDiagnostic(status,reason,SurfaceFaceCachePhase::Preparation,
                      report[0]!=0,0,report[1],report[2]);
        return false;
    };
    try {
        const auto currentStatus=diagnostic_.status;
        if(!limits_.positive()||preparedReady_||
           (currentStatus!=SurfaceFaceCacheStatus::None&&
            currentStatus!=SurfaceFaceCacheStatus::Ready)||
           plannedFacets.empty()||
           canonicalFaces.size()>limits_.maxFaces||
           canonicalFaces.size()!=canonicalFaces_.size()||
           canonicalFaces.size()!=slots_.size()||
           !std::equal(canonicalFaces.begin(),canonicalFaces.end(),
                       canonicalFaces_.begin())||
           !std::is_sorted(canonicalFaces.begin(),canonicalFaces.end())||
           std::adjacent_find(canonicalFaces.begin(),canonicalFaces.end())!=
                canonicalFaces.end()||
           !std::is_sorted(plannedFacets.begin(),plannedFacets.end())||
           std::adjacent_find(plannedFacets.begin(),plannedFacets.end())!=
                plannedFacets.end())
            return stop(SurfaceFaceCacheStatus::Failed,
                        SurfaceFaceCacheReason::Invalid);

        const auto fits=[&](std::size_t rows,std::size_t capacity) noexcept {
            return rows<=capacity&&plannedFacets.size()<=capacity-rows;
        };
        if(!fits(capacities.committedFacePointRows,
                 capacities.committedFacePointCapacity)||
           !fits(capacities.preparedFacePointRows,
                 capacities.preparedFacePointCapacity)||
           !fits(capacities.committedLedgerRows,
                 capacities.committedLedgerCapacity)||
           !fits(capacities.preparedLedgerRows,
                 capacities.preparedLedgerCapacity)||
           !fits(capacities.committedSemanticStateRows,
                 capacities.committedSemanticStateCapacity)||
           !fits(capacities.preparedSemanticStateRows,
                 capacities.preparedSemanticStateCapacity)||
           !diagnosticChannel_||
           diagnosticChannel_->processedStorage.size()<
                canonicalFaces.size()||
           diagnosticChannel_->unattemptedStorage.size()<
                canonicalFaces.size())
            return stop(SurfaceFaceCacheStatus::Incomplete,
                        SurfaceFaceCacheReason::BudgetReached);

        std::uint64_t nodes=0,cells=0,constraints=0,mappings=0,liveRows=0;
        for(std::size_t i=0;i<slots_.size();++i) {
            const auto& entry=slots_[i];
            if(!entry||!entry->complete||
               entry->facePersistentId!=canonicalFaces[i])
                return stop(SurfaceFaceCacheStatus::Failed,
                            SurfaceFaceCacheReason::Invalid);
            EntryRowCounts counts;
            if(!countEntryRows(*entry,counts)||
               !add64(nodes,counts.nodes)||!add64(cells,counts.cells)||
               !add64(constraints,counts.constraints)||
               !add64(mappings,counts.mappings)||
               !add64(liveRows,counts.total))
                return stop(SurfaceFaceCacheStatus::Incomplete,
                            SurfaceFaceCacheReason::BudgetReached);
        }
        if(nodes>limits_.maxCachedNodeRows||
           cells>limits_.maxCachedCellRows||
           constraints>limits_.maxCachedConstraintRows||
           mappings>limits_.maxCachedMappingRows)
            return stop(SurfaceFaceCacheStatus::Incomplete,
                        SurfaceFaceCacheReason::BudgetReached);

        facetProposalDirtyFaces_.reserve((std::min)(
            canonicalFaces.size(),plannedFacets.size()));
        FacePersistentId prior=0;
        bool havePrior=false;
        for(const auto& facet:plannedFacets) {
            const auto face=facet.facePersistentId;
            if(face==0||!std::binary_search(canonicalFaces.begin(),
                                            canonicalFaces.end(),face))
                return stop(SurfaceFaceCacheStatus::Failed,
                            SurfaceFaceCacheReason::Invalid);
            if(!havePrior||face!=prior) {
                facetProposalDirtyFaces_.push_back(face);
                prior=face;
                havePrior=true;
            }
        }
        if(facetProposalDirtyFaces_.size()>
               limits_.maxDirtyFacesPerBatch)
            return stop(SurfaceFaceCacheStatus::Incomplete,
                        SurfaceFaceCacheReason::BudgetReached);

        replacementScratch_.reserve(canonicalFaces.size());
        if(replacementScratch_.capacity()<canonicalFaces.size()||
           facetProposalDirtyFaces_.capacity()<
                facetProposalDirtyFaces_.size())
            return stop(SurfaceFaceCacheStatus::Incomplete,
                        SurfaceFaceCacheReason::BudgetReached);

        std::uint64_t work=0;
        const auto faceRows=static_cast<std::uint64_t>(
            canonicalFaces.size());
        const auto plannedRows=static_cast<std::uint64_t>(
            plannedFacets.size());
        if(!add64(work,faceRows)||!add64(work,faceRows)||
           !add64(work,plannedRows)||!add64(work,liveRows)||
           !add64(work,liveRows)||
           !consumePreparationWork(work,
                SurfaceFaceCachePhase::Preparation,0))
            return stop(SurfaceFaceCacheStatus::Incomplete,
                        SurfaceFaceCacheReason::BudgetReached);

        std::size_t accounted=sizeof(*this);
        const auto controlBlockEstimate=4*sizeof(void*)+2*sizeof(std::size_t);
        if(!addChecked<FacePersistentId>(
                accounted,canonicalFaces_.capacity())||
           !addChecked<CacheSlotOwner>(accounted,slots_.capacity())||
           !addChecked<FacePersistentId>(
                accounted,preparedFaces_.capacity())||
           !addChecked<CacheSlotOwner>(accounted,prepared_.capacity())||
           !addChecked<FacePersistentId>(
                accounted,facetProposalDirtyFaces_.capacity())||
           !addChecked<CacheSlotOwner>(
                accounted,replacementScratch_.capacity())||
           !addChecked<FacePersistentId>(
                accounted,attemptFaces_.capacity())||
           !addChecked<FacePersistentId>(
                accounted,expectedDirtyFaces_.capacity())||
           !addChecked<ProductionSurfaceFaceMetrics>(
                accounted,attemptMetrics_.capacity())||
           !addChecked<SurfaceSizingSampleRow>(
                accounted,globalObservedSizing_.capacity())||
           !addChecked<SurfaceSizingSampleRow>(
                accounted,canonicalSizingSequence_.capacity())||
           !addChecked<SurfaceSizingSampleRow>(
                accounted,preparedSizingSequence_.capacity())||
           !addMapEstimate(accounted,observed_)||
           slots_.size()>
                (std::numeric_limits<std::size_t>::max)()/
                    controlBlockEstimate||
           !addRaw(accounted,slots_.size()*controlBlockEstimate))
            return stop(SurfaceFaceCacheStatus::Incomplete,
                        SurfaceFaceCacheReason::BudgetReached);
        for(const auto& [face,rows]:observed_) {
            (void)face;
            if(!addChecked<SurfaceSizingSampleRow>(
                    accounted,rows.capacity()))
                return stop(SurfaceFaceCacheStatus::Incomplete,
                            SurfaceFaceCacheReason::BudgetReached);
        }
        for(const auto& entry:slots_) {
            std::size_t entryBytes=0;
            if(!accountedEntry(*entry,entryBytes)||
               !addRaw(accounted,entryBytes))
                return stop(SurfaceFaceCacheStatus::Incomplete,
                            SurfaceFaceCacheReason::BudgetReached);
        }
        metrics_.accountedByteHighWater=(std::max)(
            metrics_.accountedByteHighWater,accounted);
        if(accounted>limits_.maxSurfaceCacheAccountedBytes)
            return stop(SurfaceFaceCacheStatus::Incomplete,
                        SurfaceFaceCacheReason::BudgetReached);

        proposal.ready=true;
        proposal.plannedRows=plannedFacets.size();
        proposal.dirtyFaceRows=facetProposalDirtyFaces_.size();
        proposal.replacementSlotRows=canonicalFaces.size();
        proposal.replacementSlotCapacity=replacementScratch_.capacity();
        proposal.workUnits=work;
        proposal.accountedBytes=accounted;
        return true;
    } catch(const std::bad_alloc&) {
        return stop(SurfaceFaceCacheStatus::Failed,
                    SurfaceFaceCacheReason::AllocationFailed);
    } catch(...) {
        return stop(SurfaceFaceCacheStatus::Failed,
                    SurfaceFaceCacheReason::Invalid);
    }
}

void SurfaceFaceCacheSession::abandonFacetProposal() noexcept {
    facetProposalDirtyFaces_.clear();
}

bool SurfaceFaceCacheSession::acquireReplacementScratch(
    std::size_t required,
    std::vector<std::shared_ptr<const CachedSurfaceFaceArtifact>>& out)
    noexcept {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::ReplacementPreparation);
#endif
    if(!out.empty()) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::Invalid,
                      SurfaceFaceCachePhase::Preparation,true,0,0,0);
        return false;
    }
    try {
        if(replacementScratch_.capacity()<required)
            replacementScratch_.reserve(required);
        if(replacementScratch_.capacity()<required) {
            setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                          SurfaceFaceCacheReason::BudgetReached,
                          SurfaceFaceCachePhase::Preparation,true,0,0,0);
            return false;
        }
        out.swap(replacementScratch_);
        out.clear();
        return true;
    } catch(const std::bad_alloc&) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::AllocationFailed,
                      SurfaceFaceCachePhase::Preparation,true,0,0,0);
    } catch(...) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::Invalid,
                      SurfaceFaceCachePhase::Preparation,true,0,0,0);
    }
    return false;
}

void SurfaceFaceCacheSession::recycleReplacementScratch(
    std::vector<std::shared_ptr<const CachedSurfaceFaceArtifact>>&& rows)
    noexcept {
    rows.clear();
    if(rows.capacity()>replacementScratch_.capacity())
        replacementScratch_.swap(rows);
}

bool SurfaceFaceCacheSession::prepareReplacement(
    std::vector<std::shared_ptr<const CachedSurfaceFaceArtifact>> next,
    std::span<const FacePersistentId> canonicalFaces) noexcept {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::ReplacementPreparation);
#endif
    abandonPrepared();
    try {
        if(!limits_.positive()||attemptLookupOrdinal_!=attemptFaces_.size()||
           metrics_.cacheCommits==
               (std::numeric_limits<std::uint64_t>::max)()||
           canonicalFaces.size()>limits_.maxFaces||
           next.size()!=canonicalFaces.size()||
           !std::is_sorted(canonicalFaces.begin(),canonicalFaces.end())||
           std::adjacent_find(canonicalFaces.begin(),canonicalFaces.end())!=
               canonicalFaces.end()) {
            setDiagnostic(SurfaceFaceCacheStatus::Failed,
                          SurfaceFaceCacheReason::Invalid,
                          SurfaceFaceCachePhase::Preparation,true,0,0,0);
            return false;
        }
        std::uint64_t nodes=0,cells=0,constraints=0,mappings=0,work=0;
        for(const auto& entry:next) {
            if(!entry||!entry->complete) {
                setDiagnostic(SurfaceFaceCacheStatus::Failed,
                              SurfaceFaceCacheReason::Invalid,
                              SurfaceFaceCachePhase::Preparation,true,0,0,0);
                return false;
            }
            EntryRowCounts counts;
            if(!countEntryRows(*entry,counts)||
               !add64(nodes,counts.nodes)||!add64(cells,counts.cells)||
               !add64(constraints,counts.constraints)||
               !add64(mappings,counts.mappings)||
               !add64(work,counts.total)) {
                setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                              SurfaceFaceCacheReason::BudgetReached,
                              SurfaceFaceCachePhase::Preparation,true,0,0,0);
                return false;
            }
        }
        if(nodes>limits_.maxCachedNodeRows||cells>limits_.maxCachedCellRows||
           constraints>limits_.maxCachedConstraintRows||
           mappings>limits_.maxCachedMappingRows) {
            setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                          SurfaceFaceCacheReason::BudgetReached,
                          SurfaceFaceCachePhase::Preparation,true,0,0,0);
            return false;
        }
        if(!consumePreparationWork(work,SurfaceFaceCachePhase::Preparation,0))
            return false;
        preparedSizingSequence_=globalObservedSizing_;
        preparedFaces_.assign(canonicalFaces.begin(),canonicalFaces.end());
        prepared_=std::move(next);
        std::vector<const CachedSurfaceFaceArtifact*> uniqueEntries;
        uniqueEntries.reserve(slots_.size()+prepared_.size());
        for(const auto& entry:slots_)if(entry)uniqueEntries.push_back(entry.get());
        for(const auto& entry:prepared_)
            if(entry)uniqueEntries.push_back(entry.get());
        std::sort(uniqueEntries.begin(),uniqueEntries.end(),
            std::less<const CachedSurfaceFaceArtifact*>{});
        uniqueEntries.erase(std::unique(uniqueEntries.begin(),
                                        uniqueEntries.end()),
                            uniqueEntries.end());
        std::uint64_t liveRows=0;
        for(const auto* entry:uniqueEntries) {
            EntryRowCounts counts;
            if(!countEntryRows(*entry,counts)||
               !add64(liveRows,counts.total)) {
                abandonPrepared();
                setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                              SurfaceFaceCacheReason::BudgetReached,
                              SurfaceFaceCachePhase::Preparation,true,0,0,0);
                return false;
            }
        }
        const auto slotCount=static_cast<std::uint64_t>(slots_.size());
        const auto preparedCount=static_cast<std::uint64_t>(prepared_.size());
        if(preparedCount>
           (std::numeric_limits<std::uint64_t>::max)()-slotCount) {
            abandonPrepared();
            setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                          SurfaceFaceCacheReason::BudgetReached,
                          SurfaceFaceCachePhase::Preparation,true,0,0,0);
            return false;
        }
        preparedHighWaterSlots_=slotCount+preparedCount;
        preparedHighWaterRows_=liveRows;
        metrics_.highWaterSlots=(std::max)(metrics_.highWaterSlots,
                                           preparedHighWaterSlots_);
        metrics_.highWaterRows=(std::max)(metrics_.highWaterRows,
                                          preparedHighWaterRows_);
        // Accounting is a distinct pass over the simultaneous old/new
        // immutable owners, so charge its canonical logical rows once.
        if(!consumePreparationWork(liveRows,
                                   SurfaceFaceCachePhase::Preparation,0)) {
            abandonPrepared();
            return false;
        }
        std::size_t accounted=sizeof(*this);
        const auto controlBlockEstimate=4*sizeof(void*)+2*sizeof(std::size_t);
        if(!addChecked<FacePersistentId>(
                accounted,canonicalFaces_.capacity())||
           !addChecked<CacheSlotOwner>(accounted,slots_.capacity())||
           !addChecked<FacePersistentId>(
                accounted,preparedFaces_.capacity())||
           !addChecked<CacheSlotOwner>(accounted,prepared_.capacity())||
           !addChecked<const CachedSurfaceFaceArtifact*>(
                accounted,uniqueEntries.capacity())||
           !addChecked<ProductionSurfaceFaceMetrics>(
                accounted,attemptMetrics_.capacity())||
           !addChecked<SurfaceSizingSampleRow>(
                accounted,globalObservedSizing_.capacity())||
           !addChecked<SurfaceSizingSampleRow>(
                accounted,canonicalSizingSequence_.capacity())||
           !addChecked<SurfaceSizingSampleRow>(
                accounted,preparedSizingSequence_.capacity())||
           !addMapEstimate(accounted,observed_)||
           uniqueEntries.size()>
                (std::numeric_limits<std::size_t>::max)()/
                    controlBlockEstimate||
           !addRaw(accounted,
                uniqueEntries.size()*controlBlockEstimate)) {
            abandonPrepared();
            setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                          SurfaceFaceCacheReason::BudgetReached,
                          SurfaceFaceCachePhase::Preparation,true,0,0,0);
            return false;
        }
        for(const auto& [face,rows]:observed_) {
            (void)face;
            if(!addChecked<SurfaceSizingSampleRow>(
                    accounted,rows.capacity())) {
                abandonPrepared();
                setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                              SurfaceFaceCacheReason::BudgetReached,
                              SurfaceFaceCachePhase::Preparation,true,0,0,0);
                return false;
            }
        }
        for(const auto* entry:uniqueEntries) {
            std::size_t entryBytes=0;
            if(!accountedEntry(*entry,entryBytes)||
               !addRaw(accounted,entryBytes)) {
                abandonPrepared();
                setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                              SurfaceFaceCacheReason::BudgetReached,
                              SurfaceFaceCachePhase::Preparation,true,0,0,0);
                return false;
            }
        }
        metrics_.accountedByteHighWater=(std::max)(
            metrics_.accountedByteHighWater,accounted);
        if(accounted>limits_.maxSurfaceCacheAccountedBytes) {
            abandonPrepared();
            setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                          SurfaceFaceCacheReason::BudgetReached,
                          SurfaceFaceCachePhase::Preparation,true,0,0,0);
            return false;
        }
        preparedAccountedBytes_=accounted;
        preparedReady_=true;
        diagnostic_.workUnits=attemptWorkUnits();
        diagnostic_.accountedBytes=accounted;
        return true;
    } catch(const std::bad_alloc&) {
        abandonPrepared();
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::AllocationFailed,
                      SurfaceFaceCachePhase::Preparation,true,0,0,0);
    } catch(...) {
        abandonPrepared();
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::Invalid,
                      SurfaceFaceCachePhase::Preparation,true,0,0,0);
    }
    return false;
}

void SurfaceFaceCacheSession::commitPrepared() noexcept {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::CommitTail);
#endif
    if(!preparedReady_)return;
    canonicalFaces_.swap(preparedFaces_);
    slots_.swap(prepared_);
    canonicalSizingSequence_.swap(preparedSizingSequence_);
    ++metrics_.cacheCommits;
    metrics_.highWaterSlots=(std::max)(metrics_.highWaterSlots,
        preparedHighWaterSlots_);
    metrics_.highWaterRows=(std::max)(metrics_.highWaterRows,
        preparedHighWaterRows_);
    metrics_.accountedByteHighWater=(std::max)(
        metrics_.accountedByteHighWater,preparedAccountedBytes_);
    abandonPrepared();
}

void SurfaceFaceCacheSession::abandonPrepared() noexcept {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::CommitTail);
#endif
    preparedFaces_.clear();
    recycleReplacementScratch(std::move(prepared_));
    preparedSizingSequence_.clear();
    preparedHighWaterSlots_=0;
    preparedHighWaterRows_=0;
    preparedAccountedBytes_=0;
    preparedReady_=false;
}

void SurfaceFaceCacheSession::setDiagnostic(
    SurfaceFaceCacheStatus status,SurfaceFaceCacheReason reason,
    SurfaceFaceCachePhase phase,bool complete,FacePersistentId failedFace,
    std::uint64_t processed,std::uint64_t unattempted) noexcept {
    const bool rejected=status!=SurfaceFaceCacheStatus::None&&
        status!=SurfaceFaceCacheStatus::Ready&&failedFace!=0&&
        !(diagnostic_.status!=SurfaceFaceCacheStatus::None&&
          diagnostic_.status!=SurfaceFaceCacheStatus::Ready&&
          diagnostic_.failedFace==failedFace);
    if(rejected&&metrics_.facesRejected!=
       (std::numeric_limits<std::uint64_t>::max)())
        ++metrics_.facesRejected;
    diagnostic_.status=status; diagnostic_.reason=reason;
    diagnostic_.phase=phase; diagnostic_.reportComplete=complete;
    diagnostic_.failedFace=failedFace; diagnostic_.processedFaces=processed;
    diagnostic_.unattemptedFaces=unattempted;
    diagnostic_.workUnits=attemptWorkUnits();
    diagnostic_.accountedBytes=metrics_.accountedByteHighWater;
    if(diagnosticChannel_) {
        diagnosticChannel_->status=status;
        diagnosticChannel_->reason=reason;
        diagnosticChannel_->phase=phase;
        diagnosticChannel_->cacheWorkUnits=diagnostic_.workUnits;
        diagnosticChannel_->cacheAccountedBytes=diagnostic_.accountedBytes;
    }
}

bool SurfaceFaceCacheSession::noteBuilt(
    const CachedSurfaceFaceArtifact& entry) noexcept {
    const auto maximum=(std::numeric_limits<std::uint64_t>::max)();
    EntryRowCounts counts;
    std::uint64_t copied=0;
    if(!countEntryRows(entry,counts)||!copiedEntryBytes(entry,copied)) {
        setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                      SurfaceFaceCacheReason::BudgetReached,
                      SurfaceFaceCachePhase::Preparation,true,
                      entry.facePersistentId,0,0);
        return false;
    }
    if(metrics_.facesBuilt==maximum||
       counts.nodes>maximum-metrics_.dirtyNodeRows||
       counts.cells>maximum-metrics_.dirtyCellRows||
       counts.constraints>maximum-metrics_.constraintRows||
       counts.constraints>maximum-metrics_.dirtyConstraintRows||
       counts.mappings>maximum-metrics_.mappingRows||
       counts.mappings>maximum-metrics_.dirtyMappingRows||
       copied>maximum-metrics_.copiedBytes) {
        setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                      SurfaceFaceCacheReason::BudgetReached,
                      SurfaceFaceCachePhase::Preparation,true,
                      entry.facePersistentId,0,0);
        return false;
    }
    // Snapshotting a dirty face is preparation traffic, never replay work.
    if(!consumePreparationWork(counts.total,
                               SurfaceFaceCachePhase::Preparation,
                               entry.facePersistentId)) {
        return false;
    }
    metrics_.dirtyNodeRows+=counts.nodes;
    metrics_.dirtyCellRows+=counts.cells;
    metrics_.constraintRows+=counts.constraints;
    metrics_.dirtyConstraintRows+=counts.constraints;
    metrics_.mappingRows+=counts.mappings;
    metrics_.dirtyMappingRows+=counts.mappings;
    metrics_.copiedBytes+=copied;
    try { attemptMetrics_.push_back(entry.metrics); }
    catch(...) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::AllocationFailed,
                      SurfaceFaceCachePhase::Preparation,true,
                      entry.facePersistentId,0,0);
        return false;
    }
    ++metrics_.facesBuilt;
    return true;
}

bool SurfaceFaceCacheSession::noteReused(
    const CachedSurfaceFaceArtifact& entry) noexcept {
    const auto maximum=(std::numeric_limits<std::uint64_t>::max)();
    EntryRowCounts counts;
    if(!countEntryRows(entry,counts)) {
        setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                      SurfaceFaceCacheReason::BudgetReached,
                      SurfaceFaceCachePhase::Rehydration,true,
                      entry.facePersistentId,0,0);
        return false;
    }
    if(metrics_.facesReused==maximum||metrics_.facesReplayed==maximum||
       counts.nodes>maximum-metrics_.cleanNodeRows||
       counts.cells>maximum-metrics_.cleanCellRows||
       counts.constraints>maximum-metrics_.constraintRows||
       counts.constraints>maximum-metrics_.cleanConstraintRows||
       counts.mappings>maximum-metrics_.mappingRows||
       counts.mappings>maximum-metrics_.cleanMappingRows) {
        setDiagnostic(SurfaceFaceCacheStatus::Incomplete,
                      SurfaceFaceCacheReason::BudgetReached,
                      SurfaceFaceCachePhase::Rehydration,true,
                      entry.facePersistentId,0,0);
        return false;
    }
    if(!consumeReplayWork(counts.total,
                          SurfaceFaceCachePhase::Rehydration,
                          entry.facePersistentId)) {
        return false;
    }
    ++metrics_.facesReused;
    ++metrics_.facesReplayed;
    metrics_.cleanNodeRows+=counts.nodes;
    metrics_.cleanCellRows+=counts.cells;
    metrics_.constraintRows+=counts.constraints;
    metrics_.cleanConstraintRows+=counts.constraints;
    metrics_.mappingRows+=counts.mappings;
    metrics_.cleanMappingRows+=counts.mappings;
    try {
        ProductionSurfaceFaceMetrics replay;
        replay.facePersistentId=entry.facePersistentId;
        replay.construction.adjacencyAuditPassed=true;
        replay.construction.runComplete=true;
        attemptMetrics_.push_back(replay);
    }
    catch(...) {
        setDiagnostic(SurfaceFaceCacheStatus::Failed,
                      SurfaceFaceCacheReason::AllocationFailed,
                      SurfaceFaceCachePhase::Rehydration,true,
                      entry.facePersistentId,0,0);
        return false;
    }
    return true;
}

bool SurfaceFaceCacheSession::noteSetAudit(std::uint64_t work) noexcept {
    return consumeWork(metrics_.setAuditWorkUnits,work,
                       SurfaceFaceCachePhase::SetAudit,0);
}

std::vector<std::uint8_t> encodeSurfaceFaceCacheAuthority(
    ArtifactStamp setStamp,const SurfaceFaceRequest& request,
    const ChartCertificate& certificate,
    const std::map<SharedCurveNodeKey,SharedCurveNodeBinding>& ledger,
    const mesh::Mesh& destination,
    const Vec2& parameterLo,const Vec2& parameterHi,
    const Vec2& metricScale) {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::AuthorityKey);
#endif
    // Invariant authority is deliberately encoded separately from mutable
    // discretization structure.  A changed provider/certificate/chart/domain
    // can therefore never masquerade as an ordinary dirty-face rebuild.
    std::vector<std::uint8_t> invariant;
    invariant.reserve(384+certificate.certifiedRegionBoundary.size()*16);
    appendScalar(invariant,std::uint32_t{1}); appendStamp(invariant,setStamp);
    appendStamp(invariant,request.faceStamp);
    appendScalar(invariant,request.chart);
    appendClass(invariant,request.faceClass);
    appendScalar(invariant,request.desiredNormalSign);
    appendScalar(invariant,request.expectedGeometryStamp.identity);
    appendScalar(invariant,request.expectedGeometryStamp.revision);
    appendStamp(invariant,request.expectedChartCertificateStamp);
    appendScalar(invariant,request.expectedCandidateProviderStamp.persistentId);
    appendScalar(invariant,request.expectedCandidateProviderStamp.revision);
    appendScalar(invariant,certificate.providerStamp.identity);
    appendScalar(invariant,certificate.providerStamp.revision);
    appendScalar(invariant,certificate.chart);
    appendStamp(invariant,certificate.certificateStamp);
    appendScalar(invariant,certificate.chartPositiveToNormalSign);
    appendScalar(invariant,certificate.regular);
    appendScalar(invariant,certificate.periodic);
    appendScalar(invariant,certificate.seam);
    appendScalar(invariant,certificate.certifiedRegionBoundary.size());
    for(const auto p:certificate.certifiedRegionBoundary)
        appendVec2(invariant,p);
    appendVec2(invariant,parameterLo); appendVec2(invariant,parameterHi);
    appendVec2(invariant,metricScale);

    std::vector<std::uint8_t> out;
    out.reserve(invariant.size()+1024+request.vertices.size()*96+
                request.refinedNodes.size()*128);
    appendScalar(out,std::uint32_t{3});
    appendScalar(out,static_cast<std::uint64_t>(invariant.size()));
    out.insert(out.end(),invariant.begin(),invariant.end());
    // Sizing-derived physical-policy values deliberately do not participate
    // in the structural key.  Their complete evaluation sequence and exact
    // result bits are bound by sizingSamples; including the derived outputs
    // here would turn a changed/stateful sizing field into a cache miss and
    // silently rebuild it instead of rejecting the authority change.
    appendScalar(out,request.vertices.size());
    for(const auto& vertex:request.vertices) {
        appendScalar(out,vertex.faceLocalVertexId);
        appendVec2(out,vertex.parameter.p); appendScalar(out,vertex.parameter.chart);
        appendClass(out,vertex.classification);
    }
    appendScalar(out,request.refinedNodes.size());
    for(const auto& node:request.refinedNodes) {
        appendIdentity(out,node.identity); appendVec2(out,node.parameter.p);
        appendScalar(out,node.parameter.chart); appendClass(out,node.classification);
    }
    appendScalar(out,request.outer.pieces.size());
    for(const auto& chain:request.outer.pieces)appendChain(out,chain);
    appendScalar(out,request.holes.size());
    for(const auto& loop:request.holes) {
        appendScalar(out,loop.pieces.size());
        for(const auto& chain:loop.pieces)appendChain(out,chain);
    }
    appendScalar(out,request.embedded.size());
    for(const auto& chain:request.embedded)appendChain(out,chain);
    std::set<SharedCurveNodeKey> referenced;
    const auto collect=[&](const SurfaceChain& chain) {
        for(const auto& use:chain.uses)if(use.sharedKey)referenced.insert(*use.sharedKey);
    };
    for(const auto& chain:request.outer.pieces)collect(chain);
    for(const auto& loop:request.holes)for(const auto& chain:loop.pieces)collect(chain);
    for(const auto& chain:request.embedded)collect(chain);
    appendScalar(out,referenced.size());
    for(const auto& key:referenced) {
        const auto row=ledger.find(key);
        if(row==ledger.end())throw std::logic_error("missing cache binding");
        const auto global=row->second.existingGlobalIndex;
        if(global<0||global>=destination.nodeCount())
            throw std::logic_error("invalid cache binding payload");
        // The global index is attempt-local replay state.  The persistent
        // structured key and exact classification are the cache authority;
        // the current global binding is resolved and audited on rehydration.
        appendSharedKey(out,key);
        appendVec3(out,destination.node(global));
        appendClass(out,row->second.classification);
    }
    return out;
}

std::shared_ptr<const CachedSurfaceFaceArtifact> snapshotSurfaceFaceForCache(
    const SurfaceFaceResult& face,std::span<const std::uint8_t> authority,
    std::span<const SurfaceSizingSampleRow> sizing,
    const ProductionSurfaceFaceMetrics& metrics,
    const std::map<SharedCurveNodeKey,SharedCurveNodeBinding>& ledger) {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::FaceSnapshot);
#endif
    if(!face.localWorkspace||face.localHandles.size()!=face.localNodes.size())
        throw std::logic_error("cannot cache incomplete face");
    auto entry=std::make_shared<CachedSurfaceFaceArtifact>();
    entry->facePersistentId=face.faceStamp.persistentId;
    entry->authorityBytes.assign(authority.begin(),authority.end());
    const auto invariant=invariantAuthoritySpan(authority);
    if(!invariant)throw std::logic_error("invalid cache authority envelope");
    entry->invariantAuthorityBytes.assign(
        invariant->begin(),invariant->end());
    entry->sizingSamples.assign(sizing.begin(),sizing.end());
    entry->faceStamp=face.faceStamp; entry->chart=face.chart;
    entry->faceClass=face.faceClass; entry->desiredNormalSign=face.desiredNormalSign;
    entry->physicalPolicy=face.physicalPolicy;
    entry->expectedGeometryStamp=face.expectedGeometryStamp;
    entry->chartCertificate=face.chartCertificate;
    entry->localNodes=face.localNodes; entry->localParams=face.localParams;
    entry->localMappingParams=face.localMappingParams;
    entry->refinedIdentityByLocalNode=face.refinedIdentityByLocalNode;
    entry->localCells=face.localCells; entry->boundaryLoops=face.boundaryLoops;
    entry->normalizedChains=face.normalizedChains;
    entry->segmentOrigins=face.segmentOrigins; entry->audit=face.audit;
    entry->metrics=metrics;
    std::map<topo::NodeHandle,std::size_t> local;
    for(std::size_t i=0;i<face.localHandles.size();++i)
        local.emplace(face.localHandles[i],i);
    entry->persistentEdges.reserve(face.persistentEdges.size());
    for(const auto& edge:face.persistentEdges) {
        const auto a=findLocal(local,edge.localEndpoints[0]);
        const auto b=findLocal(local,edge.localEndpoints[1]);
        if(!a||!b)throw std::logic_error("foreign cached persistent edge");
        entry->persistentEdges.push_back({edge.faceStamp,edge.bndEntityId,
            edge.chainPersistentId,edge.originalSegmentOrdinal,
            edge.normalizedSegmentOrdinal,edge.role,{*a,*b},
            edge.endpointSharedKeys,edge.classification,edge.required});
    }
    std::set<SharedCurveNodeKey> referencedSharedKeys;
    for(const auto& chain:face.normalizedChains)
        for(const auto& use:chain.uses)
            if(use.sharedKey)referencedSharedKeys.insert(*use.sharedKey);
    entry->sharedKeysByLocalNode.reserve(referencedSharedKeys.size());
    for(const auto& key:referencedSharedKeys) {
        const auto binding=ledger.find(key);
        if(binding==ledger.end())
            throw std::logic_error("cached shared key has no binding");
        std::optional<std::size_t> localIndex;
        for(const auto& [localOrdinal,global]:face.existingGlobalIndex)
            if(global==binding->second.existingGlobalIndex) {
                if(localIndex&&*localIndex!=localOrdinal)
                    throw std::logic_error(
                        "cached shared binding maps multiple local nodes");
                localIndex=localOrdinal;
            }
        if(!localIndex)
            throw std::logic_error("cached shared binding has no local node");
        entry->sharedKeysByLocalNode.emplace_back(*localIndex,key);
    }
    entry->protectedRecovery=cacheReport(
        face.protectedRecovery,local,face.persistentEdges);
    entry->localRecovery=cacheReport(
        face.localRecovery,local,face.persistentEdges);
    entry->complete=true;
    if(!accountedEntry(*entry,entry->accountedCapacityBytes))
        throw std::length_error("surface cache accounting overflow");
    return entry;
}

std::unique_ptr<SurfaceFaceResult> rehydrateSurfaceFaceFromCache(
    const CachedSurfaceFaceArtifact& entry,
    const std::map<SharedCurveNodeKey,SharedCurveNodeBinding>& ledger,
    const mesh::Mesh& destination) {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::Rehydration);
#endif
    if(!entry.complete||entry.localNodes.size()!=entry.localParams.size()||
       entry.localNodes.size()!=entry.localMappingParams.size())
        throw std::logic_error("invalid cached face row counts");
    topo::NodePool nodes;
    std::vector<topo::NodeHandle> handles;
    handles.reserve(entry.localNodes.size());
    for(const auto& node:entry.localNodes)handles.push_back(nodes.add(node.coord,node.cls));
    topo::CellPool cells;
    for(const auto& cell:entry.localCells) {
        if(cell.nodes[0]>=handles.size()||cell.nodes[1]>=handles.size()||
           cell.nodes[2]>=handles.size())throw std::logic_error("invalid cached cell");
        const std::array<topo::NodeHandle,3> cn{handles[cell.nodes[0]],
            handles[cell.nodes[1]],handles[cell.nodes[2]]};
        cells.add(topo::ElementType::Tri3,cn,cell.classification,nodes);
    }
    // Rehydration owns a complete immutable Tri3 cell set.  Use the accepted
    // deterministic bulk factory so one final edge audit replaces the
    // historical complete-edge audit after every individual cell insertion.
    auto edges=topo::EdgePool::buildComplete(nodes,cells);
    for(const auto& edge:entry.persistentEdges) {
        if(edge.localEndpoints[0]>=handles.size()||
           edge.localEndpoints[1]>=handles.size())
            throw std::logic_error("invalid cached edge");
        const auto found=edges.find(nodes,handles[edge.localEndpoints[0]],
                                    handles[edge.localEndpoints[1]]);
        if(!found)throw std::logic_error("cached required edge missing");
        edges.setRequiredConstraint(*found,edge.required);
        edges.setClassification(*found,edge.classification);
    }
    auto result=std::make_unique<SurfaceFaceResult>();
    result->faceStamp=entry.faceStamp; result->chart=entry.chart;
    result->faceClass=entry.faceClass;
    result->desiredNormalSign=entry.desiredNormalSign;
    result->physicalPolicy=entry.physicalPolicy;
    result->expectedGeometryStamp=entry.expectedGeometryStamp;
    result->chartCertificate=entry.chartCertificate;
    result->localNodes=entry.localNodes; result->localParams=entry.localParams;
    result->localMappingParams=entry.localMappingParams;
    result->refinedIdentityByLocalNode=entry.refinedIdentityByLocalNode;
    result->localHandles=handles; result->localCells=entry.localCells;
    result->boundaryLoops=entry.boundaryLoops;
    result->normalizedChains=entry.normalizedChains;
    result->segmentOrigins=entry.segmentOrigins; result->audit=entry.audit;
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope mappingAllocationZone(
        SurfaceCacheAllocationZone::LocalGlobalMapping);
#endif
    std::map<Index,std::size_t> localByGlobal;
    std::map<SharedCurveNodeKey,std::size_t> localBySharedKey;
    for(const auto& [local,key]:entry.sharedKeysByLocalNode) {
        if(local>=entry.localNodes.size())
            throw std::logic_error("invalid cached shared local ordinal");
        const auto binding=ledger.find(key);
        if(binding==ledger.end()||binding->second.existingGlobalIndex<0||
           binding->second.existingGlobalIndex>=destination.nodeCount())
            throw std::logic_error("invalid cached binding");
        const auto global=binding->second.existingGlobalIndex;
        const auto& localNode=entry.localNodes[local];
        const auto& globalPoint=destination.node(global);
        if(localNode.coord.x!=globalPoint.x||
           localNode.coord.y!=globalPoint.y||
           localNode.coord.z!=globalPoint.z||
           !sameClass(localNode.cls,destination.nodeClass(global))||
           !sameClass(localNode.cls,binding->second.classification))
            throw std::logic_error("cached binding payload mismatch");
        const auto [found,inserted]=result->existingGlobalIndex.emplace(
            local,global);
        if(!inserted&&found->second!=global)
            throw std::logic_error("conflicting cached shared global mapping");
        const auto [globalRow,globalInserted]=localByGlobal.emplace(global,local);
        if(!globalInserted&&globalRow->second!=local)
            throw std::logic_error("cached global binding aliases local nodes");
        const auto [keyRow,keyInserted]=localBySharedKey.emplace(key,local);
        if(!keyInserted&&keyRow->second!=local)
            throw std::logic_error("cached shared key aliases local nodes");
    }
    result->persistentEdges.reserve(entry.persistentEdges.size());
    for(const auto& edge:entry.persistentEdges) {
        for(std::size_t endpoint=0;endpoint<2;++endpoint) {
            if(!edge.endpointSharedKeys[endpoint])continue;
            const auto found=localBySharedKey.find(
                *edge.endpointSharedKeys[endpoint]);
            if(found==localBySharedKey.end()||
               found->second!=edge.localEndpoints[endpoint])
                throw std::logic_error(
                    "cached persistent endpoint binding mismatch");
        }
        PersistentConstraintEdge restored;
        restored.faceStamp=edge.faceStamp; restored.bndEntityId=edge.bndEntityId;
        restored.chainPersistentId=edge.chainPersistentId;
        restored.originalSegmentOrdinal=edge.originalSegmentOrdinal;
        restored.normalizedSegmentOrdinal=edge.normalizedSegmentOrdinal;
        restored.role=edge.role; restored.localEndpoints={
            handles[edge.localEndpoints[0]],handles[edge.localEndpoints[1]]};
        restored.endpointSharedKeys=edge.endpointSharedKeys;
        restored.classification=edge.classification; restored.required=edge.required;
        result->persistentEdges.push_back(std::move(restored));
    }
    result->protectedRecovery=restoreReport(entry.protectedRecovery,handles);
    result->localRecovery=restoreReport(entry.localRecovery,handles);
    result->localWorkspace=std::make_unique<topo::TopologyWorkspace>(
        std::move(nodes),std::move(cells),std::move(edges),topo::FacePool{});
    return result;
}

bool auditRehydratedSurfaceFace(const SurfaceFaceResult& face,
                                const mesh::Mesh& destination) noexcept {
#if defined(PAACS_TOPO_TESTING)
    const SurfaceCacheAllocationZoneScope allocationZone(
        SurfaceCacheAllocationZone::SetAudit);
#endif
    try {
        if(!face.audit.passed()||!face.localWorkspace||
           face.localHandles.size()!=face.localNodes.size()||
           face.localParams.size()!=face.localNodes.size()||
           face.localMappingParams.size()!=face.localNodes.size()||
           face.localCells.empty())return false;
        const auto& ws=*face.localWorkspace;
        if(!topo::workspaceTopologyReady(ws.nodes(),ws.cells(),ws.edges(),
                                         ws.faces())||!ws.adjacencyReady()||
           bnd::auditTri3Manifold(ws)!=bnd::Tri3AuditResult::Ok)return false;
        const auto activeNodes=ws.nodes().activeBySequence();
        const auto activeCells=ws.cells().activeBySequence();
        if(activeNodes.size()!=face.localNodes.size()||
           activeCells.size()!=face.localCells.size())return false;
        for(std::size_t i=0;i<activeNodes.size();++i) {
            if(!(activeNodes[i]==face.localHandles[i]))return false;
            const auto row=ws.nodes().get(activeNodes[i]);
            const auto& expected=face.localNodes[i];
            if(row.coord.x!=expected.coord.x||row.coord.y!=expected.coord.y||
               row.coord.z!=expected.coord.z||
               !sameClass(row.cls,expected.cls)||
               !std::isfinite(face.localMappingParams[i].u)||
               !std::isfinite(face.localMappingParams[i].v))return false;
            if(const auto* stored=std::get_if<bnd::StoredParam>(
                    &face.localParams[i])) {
                if(!std::isfinite(stored->uv.u)||
                   !std::isfinite(stored->uv.v)||
                   stored->uv.u!=face.localMappingParams[i].u||
                   stored->uv.v!=face.localMappingParams[i].v)return false;
            } else {
                bnd::validateAffineParam(
                    std::get<bnd::AffineParam>(face.localParams[i]));
            }
        }
        std::set<std::array<std::size_t,3>> triangles;
        struct EdgeUse { std::size_t count=0; std::array<std::size_t,2> first{}; };
        std::map<std::array<std::size_t,2>,EdgeUse> edgeUses;
        std::set<std::size_t> usedNodes;
        for(std::size_t ci=0;ci<face.localCells.size();++ci) {
            const auto& cell=face.localCells[ci];
            if(cell.nodes[0]>=face.localNodes.size()||
               cell.nodes[1]>=face.localNodes.size()||
               cell.nodes[2]>=face.localNodes.size()||
               !sameClass(cell.classification,face.faceClass))return false;
            auto key=cell.nodes; std::sort(key.begin(),key.end());
            if(!triangles.insert(key).second)return false;
            const auto stored=ws.cells().get(activeCells[ci]);
            if(stored.type!=topo::ElementType::Tri3||stored.count!=3||
               !sameClass(stored.cls,cell.classification))return false;
            for(std::size_t i=0;i<3;++i) {
                if(!(stored.nodes[i]==face.localHandles[cell.nodes[i]]))
                    return false;
                usedNodes.insert(cell.nodes[i]);
                const std::array<std::size_t,2> directed{
                    cell.nodes[i],cell.nodes[(i+1)%3]};
                auto canonical=directed;
                if(canonical[1]<canonical[0])std::swap(canonical[0],canonical[1]);
                auto& use=edgeUses[canonical];
                if(use.count==0)use.first=directed;
                else if(use.count!=1||use.first==directed)return false;
                ++use.count;
                if(use.count>2)return false;
            }
        }
        if(usedNodes.size()!=face.localNodes.size())return false;
        std::set<std::array<std::size_t,2>> directedBoundary;
        for(const auto& [key,use]:edgeUses) {
            (void)key;
            if(use.count==1&&!directedBoundary.insert(use.first).second)
                return false;
        }
        std::set<std::array<std::size_t,2>> loopBoundary;
        for(const auto& loop:face.boundaryLoops) {
            if(loop.size()<3)return false;
            for(std::size_t i=0;i<loop.size();++i) {
                if(loop[i]>=face.localNodes.size()||
                   loop[(i+1)%loop.size()]>=face.localNodes.size()||
                   !loopBoundary.insert(
                        {loop[i],loop[(i+1)%loop.size()]}).second)return false;
            }
        }
        if(loopBoundary!=directedBoundary||
           face.audit.boundaryCycleCount!=face.boundaryLoops.size())
            return false;
        const auto euler=static_cast<std::int64_t>(usedNodes.size())-
            static_cast<std::int64_t>(edgeUses.size())+
            static_cast<std::int64_t>(face.localCells.size());
        if(euler!=face.audit.eulerCharacteristic)return false;
        std::set<SurfaceNodeIdentity> refinedIdentities;
        for(const auto& [local,identity]:face.refinedIdentityByLocalNode)
            if(local>=face.localNodes.size()||
               !validSurfaceNodeIdentity(identity)||
               identity.domain==
                   SurfaceNodeIdentityDomain::OriginalDiscretization||
               !refinedIdentities.insert(identity).second)return false;
        std::map<Index,std::size_t> localByGlobal;
        for(const auto& [local,global]:face.existingGlobalIndex) {
            if(local>=face.localNodes.size()||global<0||
               global>=destination.nodeCount()||
               face.localNodes[local].coord.x!=destination.node(global).x||
               face.localNodes[local].coord.y!=destination.node(global).y||
               face.localNodes[local].coord.z!=destination.node(global).z||
               !sameClass(face.localNodes[local].cls,
                          destination.nodeClass(global)))return false;
            const auto [found,inserted]=localByGlobal.emplace(global,local);
            if(!inserted&&found->second!=local)return false;
        }
        std::map<bnd::ConstraintEntityId,const PersistentConstraintEdge*>
            persistentByEntity;
        for(const auto& edge:face.persistentEdges) {
            // Every row in the retained surface artifact is a prescribed
            // constraint.  Rehydration must not allow a corrupted cache row
            // to erase that fact merely by applying the same false bit to the
            // fresh EdgePool before this comparison.
            if(!edge.required)return false;
            if(!persistentByEntity.emplace(edge.bndEntityId,&edge).second)
                return false;
            const auto found=ws.edges().find(ws.nodes(),edge.localEndpoints[0],
                                             edge.localEndpoints[1]);
            if(!found)return false;
            const auto row=ws.edges().get(*found);
            if(row.requiredConstraint!=edge.required||
               !sameClass(row.cls,edge.classification))return false;
            const auto origin=face.segmentOrigins.find(edge.bndEntityId);
            if(origin==face.segmentOrigins.end()||
               origin->second.chainPersistentId!=edge.chainPersistentId||
               origin->second.originalSegmentOrdinal!=
                    edge.originalSegmentOrdinal||
               origin->second.normalizedSegmentOrdinal!=
                    edge.normalizedSegmentOrdinal||
               origin->second.bndEntityId!=edge.bndEntityId||
               origin->second.endpoints[0].sharedKey!=
                    edge.endpointSharedKeys[0]||
               origin->second.endpoints[1].sharedKey!=
                    edge.endpointSharedKeys[1]||
               !sameClass(origin->second.classification,
                          edge.classification))return false;
        }
        if(face.segmentOrigins.size()!=persistentByEntity.size())return false;
        std::set<bnd::ConstraintEntityId> reportedEntities;
        const auto validReport=[&](const bnd::SteinerSegmentsReport& report) {
            if(!report.unattempted.empty())return false;
            for(const auto& row:report.processed) {
                if(!reportedEntities.insert(row.entity).second||
                   (row.result.status!=bnd::SteinerStatus::AlreadyPresent&&
                    row.result.status!=bnd::SteinerStatus::Recovered))
                    return false;
                const auto edge=persistentByEntity.find(row.entity);
                if(edge==persistentByEntity.end()||
                   row.result.recoveredChain.size()<2)return false;
                for(const auto node:row.result.recoveredChain)
                    if(!ws.nodes().isCurrent(node))return false;
                const auto& endpoints=edge->second->localEndpoints;
                const auto first=row.result.recoveredChain.front();
                const auto last=row.result.recoveredChain.back();
                if(!((first==endpoints[0]&&last==endpoints[1])||
                     (first==endpoints[1]&&last==endpoints[0])))return false;
                for(const auto& created:row.result.createdParams) {
                    if(!ws.nodes().isCurrent(created.first))return false;
                    if(const auto* stored=std::get_if<bnd::StoredParam>(
                            &created.second)) {
                        if(!std::isfinite(stored->uv.u)||
                           !std::isfinite(stored->uv.v))return false;
                    } else {
                        bnd::validateAffineParam(
                            std::get<bnd::AffineParam>(created.second));
                    }
                }
            }
            return true;
        };
        if(!validReport(face.protectedRecovery)||
           !validReport(face.localRecovery)||
           reportedEntities.size()!=persistentByEntity.size())return false;
        return true;
    } catch(...) { return false; }
}

} // namespace paacs::algo::t2c_detail
