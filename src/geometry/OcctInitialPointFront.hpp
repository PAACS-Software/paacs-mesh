// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Private initial volume point construction, before constrained Delaunay.
// Visit original source facets once, fill intrinsic source-clear Delaunay voids,
// then resume bounded wall-front closure. Source geometry remains immutable.
// The temporary cones are exact noncrossing proposal occupancy only. They are
// never exported as final volume cells, nor treated as a coverage certificate.
#pragma once

#include "../bnd/TetPlcGiftWrapping.hpp"
#include "OcctInitialDelaunayScorer.hpp"
#include "PlcVoidGeometry.hpp"
#include <queue>
#include "paacs/mesh/Mesh.hpp"
#include "paacs/mesh/Quality.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace paacs::geo::detail {

struct OcctInitialPointFrontLimits {
    std::uint64_t work=4'000'000, queries=4096;
    std::uint64_t workingBytes=128ULL*1024*1024;
    std::size_t sites=512, cones=2048, candidates=16'384, frontStates=16'384;
    // Zero admits the complete connected conflict cavity within live storage/work.
    // A nonzero value is an explicit private diagnostic/test cap.
    std::size_t maximumDelaunayCavityCells=0;
    bool enableAuxiliaryCompaction=true;
    bool requireManifoldFront=false;
    const std::atomic_bool* cancellation=nullptr;
    std::size_t maximumAuxiliaryCells=65'536;
    bool balanceSourceAndBulk=false;
    bool progressiveResources=false;
    std::uint64_t reservedBulkCandidates=0,reservedBulkWork=0,reservedBulkQueries=0;
    std::size_t reservedBulkSites=0;
};
struct OcctInitialPointFrontReport {
    bool initialized=false,budgetLimited=false,frontClosed=false,cancelled=false;
    std::uint64_t workUsed=0,nativeQueries=0,candidates=0;
    std::size_t sites=0,seedCones=0,reusedApices=0,propagatedFaceAccepted=0;
    std::size_t originalFaceActivations=0,propagatedFaceActivations=0,deepSites=0;
    std::size_t spacingRejects=0,shapeRejects=0,geometryRejects=0,nativeRejects=0;
    std::size_t remainingFrontFaces=0,queueRounds=0,maximumDepth=0;
    double minimumConeMeanRatio=1.,minimumConeScaledJacobian=1.;
    bool auxiliaryInitialized=false;
    std::uint64_t auxiliaryBootstrapWork=0,auxiliaryBootstrapPeakTemporaryBytes=0;
    std::size_t auxiliaryEvaluations=0,auxiliaryUnavailable=0,auxiliaryQualityRejects=0;
    std::size_t auxiliaryInvalidRejects=0,auxiliaryCommits=0;
    std::size_t auxiliaryComparisonPlans=0,auxiliaryComparisonCommits=0,auxiliaryRankedCommits=0;
    std::uint64_t auxiliaryScorerEvaluations=0,auxiliaryFullCertifications=0,auxiliaryCertificationRefusals=0;
    std::uint64_t auxiliaryScoreStateRejects=0,auxiliaryScoreSizeRejects=0,auxiliaryScoreMappingRejects=0;
    std::uint64_t auxiliaryScoreLocateMisses=0,auxiliaryScoreCavityUnavailable=0;
    std::uint64_t auxiliaryScoreStorageRejects=0,auxiliaryScorePositiveRejects=0;
    std::uint64_t auxiliaryCommitStateRejects=0,auxiliaryCommitStorageRejects=0;
    std::uint64_t auxiliaryCommitPositiveRejects=0,auxiliaryCommitValidityRejects=0,auxiliaryCommitSurvivorRejects=0;
    std::size_t auxiliaryCavityLimit=0,auxiliaryMaximumPlannedCavity=0,auxiliaryMaximumAcceptedCavity=0;
    std::uint64_t auxiliaryCompactionAttempts=0,auxiliaryCompactions=0;
    std::uint64_t auxiliaryDeadCellsReclaimed=0,auxiliaryCompactionWork=0;
    std::uint64_t auxiliaryMaximumCompactionModeledBytes=0,auxiliaryWorkingBytes=0;
    std::size_t auxiliaryStoredPeak=0,auxiliaryStoredCells=0,auxiliaryLiveCells=0;
    std::size_t auxiliaryStorageCapacity=0;
    std::size_t auxiliaryWitnessQueries=0,auxiliaryWitnessOutside=0;
    std::uint64_t auxiliaryRelevantRemoved=0,auxiliaryRelevantAdded=0;
    std::uint64_t auxiliaryUnclassifiedRemoved=0,auxiliaryUnclassifiedAdded=0;
    std::uint64_t auxiliaryGuardRemoved=0,auxiliaryGuardAdded=0;
    std::uint64_t auxiliaryRemovedPoor10=0,auxiliaryRemovedPoor25=0;
    std::uint64_t auxiliaryAddedPoor10=0,auxiliaryAddedPoor25=0;
    double auxiliaryRemovedMinimum=1.,auxiliaryAddedMinimum=1.,auxiliaryComparedAddedMinimum=1.;
    // Original nodes use min(requested h, median incident surface edge).
    // Wall-front nodes use min(requested h, 1.15*current-face median edge).
    // Bulk nodes use their cached requested field size.
    // These construction scales allow a gradual transition away from a finer
    // immutable CAD surface; they are not a mutation of the requested field.
    // Source wave and intrinsic void construction share all limits above.
    std::size_t sourcePhaseSites=0,bulkSites=0,bulkCandidates=0;
    bool sourcePhaseCompleted=false,bulkQueueExhausted=false;
    // Conservative refused/pending potential at bulk termination; not a final
    // covering-radius certificate (later closure may repair another cavity).
    double bulkMaximumUnresolvedRatio=0.;
    std::vector<double> effectiveScales;
    const char* budgetReason="none";
    std::uint64_t candidateLimit=0,siteLimit=0,queryLimit=0,workLimit=0;
    std::uint64_t reservedBulkCandidates=0,reservedBulkSites=0,sourceCandidates=0,sourceDeferredFaces=0;
    bool sourcePhaseDeferred=false;
    std::uint64_t auxiliaryImportedCells=0,auxiliaryImportWork=0,auxiliaryImportModeledBytes=0;
    std::uint64_t auxiliaryInitialCapacity=0,auxiliaryHardCapacity=0;
    std::uint64_t auxiliaryGrowthAttempts=0,auxiliaryGrowths=0,auxiliaryGrowthWork=0;
    std::uint64_t maximumCombinedModeledBytes=0,maximumHeapCapacity=0,resourceRetries=0;
    const char* auxiliaryResourceDiagnostic="not_attempted";
};

namespace occt_initial_point_front_detail {
inline constexpr double separationFactor=.65, faceGrowth=1.15;
inline std::uint64_t saturatedAdd(std::uint64_t a,std::uint64_t b) {
    return b>UINT64_MAX-a?UINT64_MAX:a+b;
}
inline std::uint64_t saturatedMultiply(std::uint64_t a,std::uint64_t b) {
    return a&&b>UINT64_MAX/a?UINT64_MAX:a*b;
}
inline std::uint64_t initialFrontStorageBytes(std::uint64_t nodes,std::uint64_t faces,
    std::uint64_t sites,std::uint64_t cones,std::uint64_t states) {
    auto bytes=saturatedAdd(8ULL*1024*1024,saturatedMultiply(saturatedAdd(nodes,sites),2048));
    for(const auto count:std::array{faces,cones,states})bytes=saturatedAdd(bytes,saturatedMultiply(count,4096));
    if(sites>1024)bytes=saturatedAdd(bytes,saturatedMultiply(sites-1024,512));
    return bytes;
}
inline bool finite(Vec3 p) {
    return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);
}
struct Shape {double mean=0.,jacobian=0.;};
inline Shape shape(const std::array<Vec3,4>& p) {
    if(predicates::orient3dFiltered(p[0],p[1],p[2],p[3])!=predicates::Sign::Positive) return {};
    const auto determinant=predicates::orient3dDeterminantEstimate(p[0],p[1],p[2],p[3]);
    double edges=0.;for(unsigned i=0;i<4;++i)for(unsigned j=i+1;j<4;++j)edges+=(p[j]-p[i]).norm2();
    if(!(determinant>0.)||!(edges>0.))return {};
    const auto mean=12.*std::pow(determinant/2.,2./3.)/edges;
    const auto jacobian=mesh::tetScaledJacobian(p[0],p[1],p[2],p[3]);
    return {std::isfinite(mean)?mean:0.,std::isfinite(jacobian)?jacobian:0.};
}
inline double meanSize(double a,double b) {
    if(b<a)std::swap(a,b);if(a==b)return a;
    const auto d=b-a,r=d/a;
    return d/(std::isfinite(r)?std::log1p(r):std::log(b)-std::log(a));
}
}

// SizeField requires eval(Vec3) and queryWork(), with 0<h(p)<=baseTarget.
// AdmitNode must certify strict material IN. It is used for temporary cone
// centroids as well as each new point. Input is immutable on every refusal or
// external callback exception. Internal budget stops return admitted points.
// Observed cancellation throws the private stop before returning a mesh.
template<class SizeField,class AdmitNode,class Charge>
mesh::Mesh appendOcctInitialPointFrontSites(const mesh::Mesh& boundary,
    const SizeField& field,double baseTarget,mesh::Classification region,
    const OcctInitialPointFrontLimits& limits,AdmitNode admitNode,Charge charge,
    OcctInitialPointFrontReport* report=nullptr) {
    using namespace occt_initial_point_front_detail;
    using Front=bnd::detail::PlcGiftWrapping;
    using Vertex=bnd::TetFacetVertexId;
    using Face=Front::Face;
    mesh::Mesh result=boundary;
    OcctInitialPointFrontReport out;
    out.candidateLimit=limits.candidates;out.siteLimit=limits.sites;
    out.queryLimit=limits.queries;out.workLimit=limits.work;
    static_assert(std::is_nothrow_move_assignable_v<OcctInitialPointFrontReport>);
    const auto done=[&]() -> mesh::Mesh {
        // Moving a Mesh can allocate on supported standard libraries. Deliver
        // observed cancellation before returning any mesh, even if the token
        // has since been reset. Moving the diagnostic itself is non-throwing.
        const bool cancelled=out.cancelled;
        if(report)*report=std::move(out);
        if(cancelled)throw OcctInitialConstructionCancelled{};
        return std::move(result);
    };
    if(limits.cancellation&&limits.cancellation->load()){out.cancelled=true;return done();}
    if(region.dim!=3||!(baseTarget>0.)||!std::isfinite(baseTarget)||!boundary.nodeCount())return done();
    const auto siteCap=limits.sites;
    const auto coneCap=limits.cones;
    const auto ns=std::size_t(boundary.nodeCount()),nf=std::size_t(boundary.cellCount());
    const auto stateCap=std::min<std::uint64_t>(limits.frontStates,
        saturatedAdd(nf,saturatedMultiply(3,coneCap)));
    if(!siteCap||!coneCap||!limits.work||!limits.queries||nf>stateCap||
        ns>=bnd::kInvalidTetFacetVertex||siteCap>=bnd::kInvalidTetFacetVertex-ns) {
        out.budgetLimited=true;out.budgetReason="initial_limit";return done();
    }
    // Conservative admission for copied source/result, two AABB indices,
    // boundary BVH, node hash/scales (including one exclusion double per node),
    // occupied cones and all queue/front states.
    // The auxiliary's historical fixed mapping reserve covers 1,024 sites.
    // Admit the extra mapping/vertex envelope when explicit larger clouds are used.
    const auto initialSites=limits.progressiveResources?std::min<std::size_t>(siteCap,64):siteCap;
    const auto initialCones=limits.progressiveResources?std::min<std::size_t>(coneCap,64):coneCap;
    const auto initialStates=limits.progressiveResources?
        std::min<std::uint64_t>(stateCap,saturatedAdd(nf,saturatedMultiply(3,initialCones))):stateCap;
    const auto payload=initialFrontStorageBytes(ns,nf,initialSites,initialCones,initialStates);
    if(payload>limits.workingBytes){out.budgetLimited=true;out.budgetReason="modeled_memory";return done();}
    struct Limit {const char* reason="work_limit";};
    struct SourceYield {};
    struct BulkYield {};
    bool sourceBudgetActive=false,bulkBudgetActive=false;
    std::uint64_t sourceWorkStop=limits.work,bulkWorkStop=limits.work;
    std::uint64_t sourceQueryStop=limits.queries,bulkQueryStop=limits.queries;
    std::uint64_t sourceCandidateStop=limits.candidates,bulkCandidateStop=limits.candidates;
    std::size_t sourceSiteStop=siteCap,bulkSiteStop=siteCap;
    const auto account=[&](std::uint64_t n=1) {
        // A soft source refusal leaves its whole face batch uncommitted. It
        // must yield before the total guard, including one expensive query.
        if(sourceBudgetActive&&(out.workUsed>sourceWorkStop||n>sourceWorkStop-out.workUsed))throw SourceYield{};
        if(bulkBudgetActive&&(out.workUsed>bulkWorkStop||n>bulkWorkStop-out.workUsed))throw BulkYield{};
        if(out.workUsed>limits.work||n>limits.work-out.workUsed)throw Limit{};
        charge(n);out.workUsed+=n;
    };
    const auto checkCancellation=[&](){
        if(out.cancelled||(limits.cancellation&&limits.cancellation->load())){
            out.cancelled=true;throw OcctInitialConstructionCancelled{};
        }
    };
    const auto tick=[&](std::uint64_t n=1){checkCancellation();account(n);};
    const auto sizeAt=[&](Vec3 p) {
        tick(field.queryWork());const auto h=field.eval(p);
        if(!(h>0.)||!std::isfinite(h)||h>baseTarget*(1+1e-12))
            throw std::runtime_error("invalid initial-point-front size");
        return std::min(h,baseTarget);
    };
    const auto native=[&](Vec3 p) {
        tick();
        if(sourceBudgetActive&&out.nativeQueries>=sourceQueryStop)throw SourceYield{};
        if(bulkBudgetActive&&out.nativeQueries>=bulkQueryStop)throw BulkYield{};
        if(out.nativeQueries>=limits.queries)throw Limit{"query_limit"};
        ++out.nativeQueries;
        if(!admitNode(p)){++out.nativeRejects;return false;}return true;
    };
    const auto materialWitness=[&](Vec3 p) {
        tick();
        if(sourceBudgetActive&&out.nativeQueries>=sourceQueryStop)throw SourceYield{};
        if(bulkBudgetActive&&out.nativeQueries>=bulkQueryStop)throw BulkYield{};
        if(out.nativeQueries>=limits.queries)throw Limit{"query_limit"};
        ++out.nativeQueries;++out.auxiliaryWitnessQueries;
        const auto inside=admitNode(p);
        if(!inside)++out.auxiliaryWitnessOutside;
        return inside;
    };
    using Key=std::array<std::int64_t,3>;
    struct Hash {std::size_t operator()(const Key& k)const noexcept {
        std::uint64_t h=1469598103934665603ULL;
        for(auto x:k){h^=std::uint64_t(x);h*=1099511628211ULL;}return std::size_t(h);
    }};
    Front front;
    front.requireManifoldFront=limits.requireManifoldFront;
    front.useLocalCandidates=true;front.useExactContactShortcut=true;
    front.chargeWork=tick;front.interiorAdmission=native;front.centroidAdmission=native;
    front.vertices.reserve(ns+initialSites);front.cells.reserve(initialCones);front.cellBoxes.reserve(initialCones);
    std::vector<std::size_t> depths;depths.reserve(ns+initialSites);
    std::vector<double> exclusionScales;exclusionScales.reserve(ns+initialSites);
    std::unordered_map<Key,std::vector<Vertex>,Hash> buckets;
    const auto origin=boundary.node(0);
    const auto width=separationFactor*baseTarget;
    const auto keyAt=[&](Vec3 p) {
        Key key{};const auto q=(p-origin)/width;const std::array a{q.x,q.y,q.z};
        for(unsigned j=0;j<3;++j){
            if(!std::isfinite(a[j])||std::abs(a[j])>double(std::numeric_limits<std::int64_t>::max()/4))
                throw std::runtime_error("initial-point-front spatial range");
            key[j]=std::int64_t(std::floor(a[j]));
        }return key;
    };
    const auto separated=[&](Vec3 p,double birth,double requested) {
        const auto k=keyAt(p);
        for(int x=-1;x<=1;++x)for(int y=-1;y<=1;++y)for(int z=-1;z<=1;++z){
            tick();const auto it=buckets.find({k[0]+x,k[1]+y,k[2]+z});
            if(it==buckets.end())continue;
            for(auto v:it->second){tick();
                // Fixed boundary geometry keeps the original birth-scale rule.
                // Interior pairs pack against their requested field sizes, so
                // a finer geometric birth step cannot leave permanent clusters.
                const auto h=v<ns?birth:requested;
                const auto radius=separationFactor*std::min(h,exclusionScales[v]);
                if((front.vertices[v].position-p).norm()<radius)return false;}
        }return true;
    };
    const auto nearby=[&](Vec3 p,double reach,Face face) {
        std::vector<std::pair<double,Vertex>> neighbors;
        const auto check=[&](Vertex v){
            tick();if(std::find(face.begin(),face.end(),v)!=face.end())return;
            const auto d=(front.vertices[v].position-p).norm2();
            if(d<=reach*reach)neighbors.emplace_back(d,v);
        };
        const auto k=keyAt(p);const auto radius=std::ceil(reach/width)+1.;
        if(radius>64.||std::pow(2*radius+1,3)>double(front.vertices.size())){
            for(Vertex v=0;v<front.vertices.size();++v)check(v);
        }else{
            const auto r=int(radius);
            for(int x=-r;x<=r;++x)for(int y=-r;y<=r;++y)for(int z=-r;z<=r;++z){
                tick();const auto it=buckets.find({k[0]+x,k[1]+y,k[2]+z});
                if(it!=buckets.end())for(auto v:it->second)check(v);
            }
        }
        std::ranges::sort(neighbors,[&](const auto& a,const auto& b){tick();return a<b;});
        if(neighbors.size()>12)neighbors.resize(12);
        return neighbors;
    };
    std::deque<Face> queue,propagationQueue;
    std::set<Face> originalFaceKeys;
    std::map<Face,std::size_t> attemptedAt;
    std::map<Face,std::size_t> compactionRetryAtSite;
    std::map<Face,std::array<std::uint64_t,3>> resourceRetryAtState;
    try {
        std::vector<std::vector<double>> incident(ns);
        for(Index i=0;i<boundary.cellCount();++i){
            tick();const auto& cell=boundary.cell(i);if(cell.type!=mesh::CellType::Tri3)return done();
            const auto* ids=boundary.cellNodes(i);Face f{};
            for(unsigned j=0;j<3;++j){
                if(ids[j]<0||std::size_t(ids[j])>=ns)return done();f[j]=Vertex(ids[j]);
                const auto length=(boundary.node(ids[j])-boundary.node(ids[(j+1)%3])).norm();
                if(!(length>0.)||!std::isfinite(length))return done();
                incident[ids[j]].push_back(length);incident[ids[(j+1)%3]].push_back(length);
            }front.boundary.push_back(f);originalFaceKeys.insert(Front::key(f));
        }
        out.effectiveScales.reserve(ns+initialSites);
        for(Index i=0;i<boundary.nodeCount();++i){
            tick();const auto p=boundary.node(i);if(!finite(p)||incident[i].empty())return done();
            auto& edges=incident[i];std::ranges::sort(edges,[&](double a,double b){tick();return a<b;});
            const auto local=edges[edges.size()/2];
            out.effectiveScales.push_back(std::min(sizeAt(p),local));
            exclusionScales.push_back(out.effectiveScales.back());
            front.vertices.push_back({p,std::uint64_t(i)+1});depths.push_back(0);
            buckets[keyAt(p)].push_back(Vertex(i));
        }
        if(!front.initialize())return done();
        // This owner keeps the initialized table fixed except for finite points
        // appended with sequence ID+1 and temporary removals of its last entry.
        front.initialOwnerCanonical=true;
        out.initialized=true;
        if(limits.workingBytes-payload<=32ULL*1024*1024){out.budgetLimited=true;out.budgetReason="auxiliary_memory";return done();}
        const auto availableAuxiliaryBytes=limits.workingBytes-payload;
        // Bounded reserves the complete lazy heap at its admitted cell ceiling.
        // Normal grows it by demand, accounting for overlap with compaction.
        constexpr std::uint64_t heapBytesPerCell=128;
        const auto auxiliaryCells=limits.progressiveResources?std::uint64_t(limits.maximumAuxiliaryCells):
            std::min<std::uint64_t>(limits.maximumAuxiliaryCells,
                (availableAuxiliaryBytes-32ULL*1024*1024)/(8192+heapBytesPerCell));
        const auto auxiliaryBytes=limits.progressiveResources?availableAuxiliaryBytes:
            availableAuxiliaryBytes-auxiliaryCells*heapBytesPerCell;
        OcctInitialDelaunayScorer auxiliary;
        bool initialized=false;
        try{
            initialized=auxiliary.initialize(boundary,region,std::size_t(auxiliaryCells),auxiliaryBytes,account,
                materialWitness,limits.maximumDelaunayCavityCells,limits.work-out.workUsed,limits.cancellation,
                limits.progressiveResources);
        }catch(const OcctInitialConstructionCancelled&){
            out.auxiliaryBootstrapWork=auxiliary.bootstrapWork;
            out.auxiliaryBootstrapPeakTemporaryBytes=auxiliary.bootstrapPeakTemporaryBytes;
            out.cancelled=true;throw;
        }
        out.auxiliaryBootstrapWork=auxiliary.bootstrapWork;
        out.auxiliaryBootstrapPeakTemporaryBytes=auxiliary.bootstrapPeakTemporaryBytes;
        out.cancelled=out.cancelled||auxiliary.cancelled;
        if(limits.progressiveResources)
            out.maximumCombinedModeledBytes=saturatedAdd(payload,
                std::min(auxiliary.maximumResourceModeledBytes,auxiliary.workingBytes()));
        if(!initialized){
            out.budgetLimited=true;out.budgetReason=limits.progressiveResources?
                auxiliary.resourceDiagnostic:"auxiliary_bootstrap";
            out.auxiliaryResourceDiagnostic=auxiliary.resourceDiagnostic;
            out.auxiliaryImportedCells=auxiliary.bootstrapImportedCells;
            out.auxiliaryImportWork=auxiliary.bootstrapImportWork;
            out.auxiliaryImportModeledBytes=auxiliary.bootstrapImportModeledBytes;
            return done();
        }
        auxiliary.compactionEnabled=limits.enableAuxiliaryCompaction;
        auxiliary.useCertifiedConeCavities=true;
        out.auxiliaryInitialized=true;
        std::uint64_t currentReservation=payload;
        const auto syncAuxiliary=[&](){
            out.budgetLimited=out.budgetLimited||auxiliary.budgetLimited;
            if(auxiliary.budgetLimited&&std::string_view(out.budgetReason)=="none")
                out.budgetReason="auxiliary_storage";
            out.cancelled=out.cancelled||auxiliary.cancelled;
            out.auxiliaryScorerEvaluations=auxiliary.evaluations;
            out.auxiliaryFullCertifications=auxiliary.fullCertifications;
            out.auxiliaryCertificationRefusals=auxiliary.certificationRefusals;
            out.auxiliaryScoreStateRejects=auxiliary.scoreStateRejects;
            out.auxiliaryScoreSizeRejects=auxiliary.scoreSizeRejects;
            out.auxiliaryScoreMappingRejects=auxiliary.scoreMappingRejects;
            out.auxiliaryScoreLocateMisses=auxiliary.scoreLocateMisses;
            out.auxiliaryScoreCavityUnavailable=auxiliary.scoreCavityUnavailable;
            out.auxiliaryScoreStorageRejects=auxiliary.scoreStorageRejects;
            out.auxiliaryScorePositiveRejects=auxiliary.scorePositiveRejects;
            out.auxiliaryCommitStateRejects=auxiliary.commitStateRejects;
            out.auxiliaryCommitStorageRejects=auxiliary.commitStorageRejects;
            out.auxiliaryCommitPositiveRejects=auxiliary.commitPositiveRejects;
            out.auxiliaryCommitValidityRejects=auxiliary.commitValidityRejects;
            out.auxiliaryCommitSurvivorRejects=auxiliary.commitSurvivorRejects;
            out.auxiliaryCavityLimit=auxiliary.cavityLimit;
            out.auxiliaryMaximumPlannedCavity=auxiliary.maximumPlannedCavity;
            out.auxiliaryMaximumAcceptedCavity=auxiliary.maximumAcceptedCavity;
            out.auxiliaryCompactionAttempts=auxiliary.compactionAttempts;
            out.auxiliaryCompactions=auxiliary.compactions;
            out.auxiliaryDeadCellsReclaimed=auxiliary.deadCellsReclaimed;
            out.auxiliaryCompactionWork=auxiliary.compactionWork;
            out.auxiliaryMaximumCompactionModeledBytes=auxiliary.maximumCompactionModeledBytes;
            out.auxiliaryWorkingBytes=auxiliary.workingBytes();
            out.auxiliaryStorageCapacity=auxiliary.storageCapacity();
            out.auxiliaryStoredPeak=auxiliary.peakStoredCells;
            out.auxiliaryStoredCells=auxiliary.storedCells();
            out.auxiliaryLiveCells=auxiliary.liveCells();
            out.auxiliaryImportedCells=auxiliary.bootstrapImportedCells;
            out.auxiliaryImportWork=auxiliary.bootstrapImportWork;
            out.auxiliaryImportModeledBytes=auxiliary.bootstrapImportModeledBytes;
            out.auxiliaryInitialCapacity=auxiliary.initialStorageCapacity;
            out.auxiliaryHardCapacity=auxiliary.hardStorageLimit();
            out.auxiliaryGrowthAttempts=auxiliary.capacityGrowthAttempts;
            out.auxiliaryGrowths=auxiliary.capacityGrowths;
            out.auxiliaryGrowthWork=auxiliary.growthWork;
            out.auxiliaryResourceDiagnostic=auxiliary.resourceDiagnostic;
            if(limits.progressiveResources)
                out.maximumCombinedModeledBytes=std::max(out.maximumCombinedModeledBytes,
                    saturatedAdd(currentReservation,std::min(auxiliary.maximumResourceModeledBytes,auxiliary.workingBytes())));
        };
        std::uint64_t heapSlots=0;
        const auto refreshWorkingBudget=[&](std::uint64_t peakHeapSlots){
            if(!limits.progressiveResources)return;
            const auto cones=std::max<std::uint64_t>(initialCones,saturatedAdd(front.cells.size(),4));
            const auto states=std::min<std::uint64_t>(stateCap,saturatedAdd(nf,saturatedMultiply(3,cones)));
            const auto current=initialFrontStorageBytes(ns,nf,
                std::max<std::uint64_t>(initialSites,saturatedAdd(out.sites,2)),cones,states);
            const auto heap=saturatedMultiply(peakHeapSlots,heapBytesPerCell);
            const auto reserved=saturatedAdd(current,heap);
            if(reserved>=limits.workingBytes||!auxiliary.setWorkingBudget(limits.workingBytes-reserved))
                throw Limit{"combined_working_memory"};
            currentReservation=reserved;
            const auto combined=saturatedAdd(reserved,auxiliary.currentModeledBytes());
            if(combined>limits.workingBytes)throw Limit{"combined_working_memory"};
            out.maximumCombinedModeledBytes=std::max(out.maximumCombinedModeledBytes,combined);
            out.maximumHeapCapacity=std::max(out.maximumHeapCapacity,peakHeapSlots);
        };
        refreshWorkingBudget(0);

        const auto runBulk=[&](){
            bulkBudgetActive=limits.balanceSourceAndBulk;
            try {
            struct Void {double ratio=0.,target=0.;Vec3 point{};std::size_t cell=0;};
            static_assert(sizeof(Void)<=128);
            const auto lower=[&](const Void& a,const Void& b){
                tick();if(a.ratio!=b.ratio)return a.ratio<b.ratio;
                return a.cell>b.cell;
            };
            using BaseHeap=std::priority_queue<Void,std::vector<Void>,decltype(lower)>;
            struct VoidHeap:BaseHeap {
                explicit VoidHeap(decltype(lower) compare,std::vector<Void>&& values):BaseHeap(compare,std::move(values)){}
                std::size_t capacity()const{return this->c.capacity();}
                void reserve(std::size_t count){this->c.reserve(count);}
            };
            std::vector<Void> storage;
            if(!limits.progressiveResources){tick(auxiliary.storageCapacity());storage.reserve(auxiliary.storageCapacity());}
            VoidHeap heap(lower,std::move(storage));
            heapSlots=heap.capacity();refreshWorkingBudget(heapSlots);
            const auto pushHeap=[&](Void item){
                if(limits.progressiveResources&&heap.size()==heap.capacity()){
                    const auto capacity=std::min<std::uint64_t>(limits.maximumAuxiliaryCells,
                        std::max<std::uint64_t>(64,saturatedMultiply(heap.capacity(),2)));
                    if(capacity<=heap.capacity())throw Limit{"auxiliary_heap_limit"};
                    refreshWorkingBudget(saturatedAdd(heap.capacity(),capacity));
                    // Charge entries that move, rather than unoccupied future slots.
                    tick(heap.size());heap.reserve(std::size_t(capacity));
                    heapSlots=heap.capacity();refreshWorkingBudget(heapSlots);
                }
                heap.push(item);
            };
            PlcVoidGeometry plc(front,tick);
            const auto push=[&](std::size_t id){
                tick();const auto& live=auxiliary.complexForDiagnostics();
                if(id>=live.cells.size()||!live.alive[id])return;
                const auto& v=live.cells[id].vertices;
                if(std::ranges::any_of(v,[&](Vertex n){return n>=ns&&n<ns+8;})){
                    return;
                }
                const auto a=live.vertices[v[0]].position;
                const auto u=live.vertices[v[1]].position-a;
                const auto v2=live.vertices[v[2]].position-a;
                const auto w=live.vertices[v[3]].position-a;
                tick(8);const auto determinant=u.dot(v2.cross(w));
                if(!(determinant>0.)||!std::isfinite(determinant)){return;}
                const auto center=a+(v2.cross(w)*u.norm2()+w.cross(u)*v2.norm2()+u.cross(v2)*w.norm2())/(2*determinant);
                if(!finite(center)){return;}
                if(center.x<=plc.bounds.lo.x||center.x>=plc.bounds.hi.x||
                   center.y<=plc.bounds.lo.y||center.y>=plc.bounds.hi.y||
                   center.z<=plc.bounds.lo.z||center.z>=plc.bounds.hi.z){return;}
                // One actual field evaluation per newly exposed cell. Cache it
                // with the immutable cell candidate; selected admission reuses it.
                const auto h=sizeAt(center);
                const auto radius=(center-a).norm();
                if(!std::isfinite(radius)){return;}
                if(!(radius>=separationFactor*h)){return;}
                const auto clipped=plc.clearance(center,radius);
                if(!(clipped>=separationFactor*h)||!std::isfinite(clipped/h)){return;}
                if(!limits.progressiveResources&&heap.size()>=auxiliary.storageCapacity())throw Limit{"auxiliary_heap_limit"};
                pushHeap({clipped/h,h,center,id});
            };
            const auto rebuildVoids=[&](){
                while(!heap.empty()){tick();heap.pop();}
                const auto& live=auxiliary.complexForDiagnostics();
                for(std::size_t id=0;id<live.cells.size();++id)push(id);
            };
            auto compactions=auxiliary.compactions;
            rebuildVoids();
            while(!heap.empty()&&out.sites<bulkSiteStop){
                // No Plan survives between batches. Compaction invalidates IDs;
                // only that storage event causes a whole-heap rebuild.
                refreshWorkingBudget(heapSlots);
                try{(void)auxiliary.prepareForBatch();}catch(...){syncAuxiliary();throw;}
                refreshWorkingBudget(heapSlots);
                syncAuxiliary();
                if(auxiliary.compactions!=compactions){
                    compactions=auxiliary.compactions;rebuildVoids();
                    if(heap.empty())break;
                }
                tick();const auto item=heap.top();heap.pop();
                const auto& live=auxiliary.complexForDiagnostics();
                if(item.cell>=live.alive.size()||!live.alive[item.cell])continue;
                if(item.ratio<separationFactor)break;
                if(bulkBudgetActive&&out.candidates>=bulkCandidateStop)throw BulkYield{};
                if(out.candidates>=limits.candidates)throw Limit{"candidate_limit"};
                ++out.candidates;++out.bulkCandidates;
                const auto reject=[&](){out.bulkMaximumUnresolvedRatio=std::max(out.bulkMaximumUnresolvedRatio,item.ratio);};
                if(!separated(item.point,item.target,item.target)){
                    ++out.spacingRejects;reject();continue;
                }
                const auto inclusion=plc.inside(item.point);
                if(inclusion!=PlcVoidGeometry::Inclusion::Inside){
                    // Exterior cells are numerical geometry, never a global
                    // mesh refusal and never a claim of material coverage.
                    continue;
                }
                if(!native(item.point)){reject();continue;}
                std::optional<OcctInitialDelaunayScorer::Plan> plan;
                const auto vertex=Vertex(front.vertices.size());
                front.vertices.push_back({item.point,front.vertices.size()+1});
                try{++out.auxiliaryEvaluations;plan=auxiliary.score(item.point,front);}
                catch(...){syncAuxiliary();front.vertices.pop_back();throw;}
                front.vertices.pop_back();syncAuxiliary();
                if(!plan||!plan->valid){
                    if(limits.progressiveResources&&auxiliary.resourceRetryRequested()){
                        ++out.resourceRetries;pushHeap(item);continue;
                    }
                    if(limits.progressiveResources&&auxiliary.budgetLimited)throw Limit{auxiliary.resourceDiagnostic};
                    ++out.auxiliaryInvalidRejects;reject();
                    if(auxiliary.compactionRequested())pushHeap(item);
                    continue;
                }
                if(!plan->available)++out.auxiliaryUnavailable;
                if(plan->comparisonAvailable)++out.auxiliaryComparisonPlans;
                out.auxiliaryRelevantRemoved+=plan->relevantRemoved;out.auxiliaryRelevantAdded+=plan->relevantAdded;
                out.auxiliaryUnclassifiedRemoved+=plan->unclassifiedRemoved;out.auxiliaryUnclassifiedAdded+=plan->unclassifiedAdded;
                out.auxiliaryGuardRemoved+=plan->guardRemoved;out.auxiliaryGuardAdded+=plan->guardAdded;
                if(!plan->qualityAdmissible){
                    ++out.auxiliaryQualityRejects;
                    reject();continue;
                }
                front.vertices.push_back({item.point,front.vertices.size()+1});
                const auto firstAdded=auxiliary.storedCells();
                bool committed=false;
                try{committed=auxiliary.commit(*plan,vertex);}
                catch(...){syncAuxiliary();front.vertices.pop_back();throw;}
                syncAuxiliary();
                if(!committed){front.vertices.pop_back();++out.auxiliaryInvalidRejects;reject();continue;}
                ++out.auxiliaryCommits;
                if(plan->available){++out.auxiliaryRankedCommits;out.auxiliaryAddedMinimum=std::min(out.auxiliaryAddedMinimum,plan->addedMinMean);}
                if(plan->comparisonAvailable){
                    ++out.auxiliaryComparisonCommits;
                    out.auxiliaryRemovedMinimum=std::min(out.auxiliaryRemovedMinimum,plan->removedMinMean);
                    out.auxiliaryComparedAddedMinimum=std::min(out.auxiliaryComparedAddedMinimum,plan->addedMinMean);
                    out.auxiliaryRemovedPoor10+=plan->removedPoor10;out.auxiliaryRemovedPoor25+=plan->removedPoor25;
                    out.auxiliaryAddedPoor10+=plan->addedPoor10;out.auxiliaryAddedPoor25+=plan->addedPoor25;
                }
                result.addNode(item.point,region);++out.sites;++out.bulkSites;
                out.effectiveScales.push_back(item.target);exclusionScales.push_back(item.target);depths.push_back(0);
                buckets[keyAt(item.point)].push_back(vertex);
                // Current cell geometry changes only in this accepted cavity.
                for(std::size_t id=firstAdded;id<auxiliary.storedCells();++id)push(id);
            }
            while(!heap.empty()){
                tick();const auto item=heap.top();
                const auto& live=auxiliary.complexForDiagnostics();
                if(item.cell<live.alive.size()&&live.alive[item.cell])break;
                heap.pop();
            }
            // Conservative diagnostic at bulk termination: includes historical
            // refused potentials, possibly repaired by another later cavity.
            // It is deliberately not a formal covering-radius certificate.
            out.bulkQueueExhausted=heap.empty();
            if(!heap.empty())out.bulkMaximumUnresolvedRatio=std::max(out.bulkMaximumUnresolvedRatio,heap.top().ratio);
            if(!heap.empty()&&out.sites>=bulkSiteStop&&bulkBudgetActive){
                out.budgetLimited=true;out.budgetReason="bulk_phase_sites";
            }
            }catch(const BulkYield&){
                // A partially scanned/rebuilt heap cannot bound all remaining
                // voids. Infinity is serialized as unavailable, never as zero.
                out.bulkQueueExhausted=false;
                out.bulkMaximumUnresolvedRatio=std::numeric_limits<double>::infinity();
                out.budgetLimited=true;out.budgetReason="bulk_phase_allowance";
            }
            bulkBudgetActive=false;
            // The points were individually PLC- and native-certified. Treat
            // them as fixed NON-SOURCE occupancy vertices so later wall cones
            // cannot hide them; the prescribed Tri3 set stays exactly unchanged.
            front.fixedVertexCount=Vertex(front.vertices.size());
            // The heap is destroyed on return; release its modeled allowance
            // before resuming the source/closure frontier.
            if(limits.progressiveResources)heapSlots=0;
        };

        // Visit an initial snapshot of original source keys once. Newly
        // exposed faces are parked until intrinsic bulk placement has run.
        bool sourcePhase=true,bulkRan=false;
        bool sourceDeferred=false;
        const auto rebuild=[&](){
            ++out.queueRounds;
            for(const auto& [key,row]:front.front){(void)row;tick();
                const auto it=attemptedAt.find(key);
                if(it==attemptedAt.end()||it->second!=front.cells.size()){
                    if(originalFaceKeys.contains(key))queue.push_back(key);
                    else propagationQueue.push_back(key);
                }
            }
        };
        rebuild();auto roundStart=front.cells.size();
        if(limits.balanceSourceAndBulk){
            // Bootstrap and queue preparation are shared prerequisites. Split
            // only the remaining actual allowance, with separate closure room.
            const auto reserve=[&](std::uint64_t total,std::uint64_t requested){
                return std::min(total-total/8,requested?requested:total/3+(total%3!=0));
            };
            const auto workRoom=limits.work-out.workUsed,queryRoom=limits.queries-out.nativeQueries;
            const auto bulkWork=reserve(workRoom,limits.reservedBulkWork);
            const auto bulkQueries=reserve(queryRoom,limits.reservedBulkQueries);
            out.reservedBulkCandidates=reserve(limits.candidates,limits.reservedBulkCandidates);
            out.reservedBulkSites=reserve(siteCap,limits.reservedBulkSites);
            bulkWorkStop=limits.work-workRoom/8;sourceWorkStop=bulkWorkStop-bulkWork;
            bulkQueryStop=limits.queries-queryRoom/8;sourceQueryStop=bulkQueryStop-bulkQueries;
            bulkCandidateStop=limits.candidates-limits.candidates/8;
            sourceCandidateStop=bulkCandidateStop-out.reservedBulkCandidates;
            bulkSiteStop=siteCap-siteCap/8;
            sourceSiteStop=bulkSiteStop-std::size_t(out.reservedBulkSites);
            sourceBudgetActive=true;
        }
        while(!bulkRan||!propagationQueue.empty()||!queue.empty()||front.cells.size()!=roundStart){
            if(sourcePhase&&limits.balanceSourceAndBulk&&out.sites>=sourceSiteStop)sourceDeferred=true;
            if(sourcePhase&&(queue.empty()||sourceDeferred)){
                sourceBudgetActive=false;
                out.sourcePhaseCompleted=queue.empty();
                if(sourceDeferred&&!queue.empty()){
                    out.sourcePhaseDeferred=true;
                    // Unpolled original source-face queue entries, not a
                    // statement about missing CAD faces or material coverage.
                    out.sourceDeferredFaces=queue.size();
                }
                const auto beforeBulk=out.sites;
                runBulk();bulkRan=true;sourcePhase=false;
                // Bulk changes the point set without adding occupancy cones.
                // Failed source faces must see the new apices on closure.
                if(out.sites!=beforeBulk)attemptedAt.clear();
                if(front.front.empty())break;
            }
            std::optional<Face> interruptedSource;
            try {
            // Prior candidate vectors have been destroyed. Compact only here,
            // before any point/Plan for the next face is staged.
            refreshWorkingBudget(heapSlots);
            try{(void)auxiliary.prepareForBatch();}
            catch(...){syncAuxiliary();throw;}
            refreshWorkingBudget(heapSlots);
            syncAuxiliary();
            if(propagationQueue.empty()&&queue.empty()){
                roundStart=front.cells.size();rebuild();
                if(propagationQueue.empty()&&queue.empty())break;
            }
            // Park propagated faces during the source wave. After bulk,
            // consume that frontier before retrying remaining source faces.
            tick();const bool propagating=!sourcePhase&&!propagationQueue.empty();
            const auto key=propagating?propagationQueue.front():queue.front();
            if(propagating)propagationQueue.pop_front();else queue.pop_front();
            if(sourcePhase)interruptedSource=key;
            const auto found=front.front.find(key);
            if(found==front.front.end())continue;
            if(const auto at=attemptedAt.find(key);at!=attemptedAt.end()&&at->second==front.cells.size())continue;
            if(attemptedAt.size()>=stateCap&&!attemptedAt.contains(key)){
                if(sourceBudgetActive)throw SourceYield{};
                throw Limit{"front_state_limit"};
            }
            attemptedAt[key]=front.cells.size();const auto f=found->second.directed;
            const auto depth=1+std::max({depths[f[0]],depths[f[1]],depths[f[2]]});
            if(depth==1)++out.originalFaceActivations;else ++out.propagatedFaceActivations;
            const auto a=front.vertices[f[0]].position,b=front.vertices[f[1]].position,c=front.vertices[f[2]].position;
            const auto u=b-a,v=c-a,n=u.cross(v);const auto n2=n.norm2();if(!(n2>0.))continue;
            const auto middle=(a+b+c)/3.;auto center=middle;
            if(u.dot(v)>0.&&(a-b).dot(c-b)>0.&&(a-c).dot(b-c)>0.)
                center=a+(v.cross(n)*u.norm2()+n.cross(u)*v.norm2())/(2*n2);
            std::array lengths{(a-b).norm(),(b-c).norm(),(c-a).norm()};std::ranges::sort(lengths);
            const auto h=std::min(sizeAt(middle),faceGrowth*lengths[1]);
            const auto r2=std::max({(center-a).norm2(),(center-b).norm2(),(center-c).norm2()});
            // A coarse/obtuse frozen face may already exceed requested h.
            // In that case this is a centroid/height proposal, not a claim that
            // all three new edges can equal h.
            const auto height=r2<h*h?std::min(std::sqrt(h*h-r2),std::sqrt(2./3.)*h):std::sqrt(2./3.)*h;
            struct Candidate {
                Vec3 point;Vertex vertex;double scale,requestedSize,score,mean,jacobian;bool fresh;
                std::optional<OcctInitialDelaunayScorer::Plan> auxiliaryPlan;
            };
            std::vector<Candidate> candidates;
            const auto offer=[&](Vec3 p,Vertex vertex,bool fresh){
                if(sourceBudgetActive&&out.candidates>=sourceCandidateStop)throw SourceYield{};
                tick();if(out.candidates>=limits.candidates)throw Limit{"candidate_limit"};++out.candidates;
                if(sourcePhase)++out.sourceCandidates;
                if(!finite(p))return;
                // Reuse this one field sample for both geometric construction
                // and persistent interior exclusion; no extra CAD query.
                const auto requested=fresh?sizeAt(p):out.effectiveScales[vertex];
                const auto scale=fresh?std::min(requested,h):out.effectiveScales[vertex];
                if(fresh&&!separated(p,scale,requested)){++out.spacingRejects;return;}
                const std::array q{a,c,b,p};const auto measured=shape(q);
                if(measured.mean<.01||measured.jacobian<1e-6){++out.shapeRejects;return;}
                const std::array scales{out.effectiveScales[f[0]],out.effectiveScales[f[2]],out.effectiveScales[f[1]],scale};
                double error=0.;
                for(unsigned i=0;i<4;++i)for(unsigned j=i+1;j<4;++j){
                    tick();const auto target=meanSize(scales[i],scales[j]);
                    const auto ratio=(q[j]-q[i]).norm()/target;
                    if(!(ratio>0.)||!std::isfinite(ratio))return;error+=std::abs(std::log(ratio));
                }
                const auto score=std::sqrt(measured.mean*measured.jacobian)/(1+error/6.);
                std::optional<OcctInitialDelaunayScorer::Plan> plan;
                auto rankedScore=score;
                if(fresh){
                    ++out.auxiliaryEvaluations;
                    front.vertices.push_back({p,front.vertices.size()+1});
                    try{plan=auxiliary.score(p,front);}
                    catch(...){syncAuxiliary();front.vertices.pop_back();throw;}
                    front.vertices.pop_back();
                    syncAuxiliary();
                    if(!plan||!plan->valid){
                        if(limits.progressiveResources&&auxiliary.budgetLimited&&!auxiliary.resourceRetryRequested())
                            throw Limit{auxiliary.resourceDiagnostic};
                        if(!(limits.progressiveResources&&auxiliary.resourceRetryRequested()))++out.auxiliaryInvalidRejects;
                        return;
                    }
                    if(!plan->available)++out.auxiliaryUnavailable;
                    if(plan->comparisonAvailable)++out.auxiliaryComparisonPlans;
                    out.auxiliaryRelevantRemoved+=plan->relevantRemoved;
                    out.auxiliaryRelevantAdded+=plan->relevantAdded;
                    out.auxiliaryUnclassifiedRemoved+=plan->unclassifiedRemoved;
                    out.auxiliaryUnclassifiedAdded+=plan->unclassifiedAdded;
                    out.auxiliaryGuardRemoved+=plan->guardRemoved;
                    out.auxiliaryGuardAdded+=plan->guardAdded;
                    if(!plan->qualityAdmissible){++out.auxiliaryQualityRejects;return;}
                    if(plan->available)rankedScore=plan->score/(1+error/6.);
                }
                candidates.push_back({p,vertex,scale,requested,rankedScore,measured.mean,measured.jacobian,fresh,std::move(plan)});
            };
            for(const auto [distance,vertex]:nearby(center,2.5*lengths[2],f)){
                (void)distance;offer(front.vertices[vertex].position,vertex,false);
            }
            if(out.sites<siteCap)for(const auto factor:std::array{1.,.85,1.15})
                offer(center-n*(height*factor/std::sqrt(n2)),Vertex(front.vertices.size()),true);
            if(limits.progressiveResources&&auxiliary.resourceRetryRequested()){
                candidates.clear();
                const std::array<std::uint64_t,3> stamp{out.sites,auxiliary.storageCapacity(),auxiliary.planEpoch()};
                const auto retried=resourceRetryAtState.find(key);
                if(retried!=resourceRetryAtState.end()&&retried->second==stamp)
                    throw Limit{"auxiliary_resource_stalled"};
                if(resourceRetryAtState.size()>=stateCap&&!resourceRetryAtState.contains(key))throw Limit{"front_state_limit"};
                resourceRetryAtState[key]=stamp;++out.resourceRetries;
                attemptedAt.erase(key);
                if(sourcePhase)queue.push_front(key);else propagationQueue.push_front(key);
                continue;
            }
            if(auxiliary.compactionRequested()){
                // Storage history must not remove a better fresh proposal and
                // thereby bias this face toward a worse fresh/reused apex.
                // Discard the ENTIRE batch before any occupancy/native commit.
                candidates.clear();
                const auto sites=auxiliary.acceptedSites();
                const auto retried=compactionRetryAtSite.find(key);
                if(retried!=compactionRetryAtSite.end()&&retried->second==sites){
                    // At most one forced storage retry per face and accepted-
                    // point generation. Genuine live-capacity pressure remains
                    // an explicit budget limit instead of a retry loop.
                    out.budgetLimited=true;
                    if(std::string_view(out.budgetReason)=="none")out.budgetReason="auxiliary_storage";
                    continue;
                }
                if(compactionRetryAtSite.size()>=stateCap&&!compactionRetryAtSite.contains(key))throw Limit{"front_state_limit"};
                compactionRetryAtSite[key]=sites;
                attemptedAt.erase(key);propagationQueue.push_front(key);
                continue;
            }
            std::stable_sort(candidates.begin(),candidates.end(),[&](const Candidate& x,const Candidate& y){
                tick();
                // Site construction is the purpose of this stage. Existing
                // apices close collisions after every admissible fresh plan
                // has failed; unchanged DT quality must not suppress new sites.
                if(x.fresh!=y.fresh)return x.fresh;
                if(x.fresh&&x.auxiliaryPlan->available!=y.auxiliaryPlan->available)
                    return x.auxiliaryPlan->available;
                if(x.score!=y.score)return x.score>y.score;
                return x.vertex<y.vertex;
            });
            bool accepted=false;
            for(const auto& candidate:candidates){
                if(front.cells.size()>=coneCap||front.front.size()+3>stateCap){
                    if(sourceBudgetActive)throw SourceYield{};
                    throw Limit{"cone_or_front_limit"};
                }
                if(candidate.fresh)front.vertices.push_back({candidate.point,front.vertices.size()+1});
                bnd::TetFacetCell cone;
                bool valid=false;
                try{valid=front.valid(Front::Front{f},candidate.vertex,cone);}
                catch(...){if(candidate.fresh)front.vertices.pop_back();throw;}
                if(!valid){if(candidate.fresh)front.vertices.pop_back();++out.geometryRejects;continue;}
                if(candidate.fresh){
                    bool committed=false;
                    try{committed=auxiliary.commit(*candidate.auxiliaryPlan,candidate.vertex);}
                    catch(...){syncAuxiliary();front.vertices.pop_back();throw;}
                    syncAuxiliary();
                    if(!committed){front.vertices.pop_back();++out.auxiliaryInvalidRejects;continue;}
                    ++out.auxiliaryCommits;
                    if(candidate.auxiliaryPlan->available){
                        ++out.auxiliaryRankedCommits;
                        out.auxiliaryAddedMinimum=std::min(out.auxiliaryAddedMinimum,candidate.auxiliaryPlan->addedMinMean);
                    }
                    if(candidate.auxiliaryPlan->comparisonAvailable){
                        ++out.auxiliaryComparisonCommits;
                        const auto& p=*candidate.auxiliaryPlan;
                        out.auxiliaryRemovedMinimum=std::min(out.auxiliaryRemovedMinimum,p.removedMinMean);
                        out.auxiliaryComparedAddedMinimum=std::min(out.auxiliaryComparedAddedMinimum,p.addedMinMean);
                        out.auxiliaryRemovedPoor10+=p.removedPoor10;out.auxiliaryRemovedPoor25+=p.removedPoor25;
                        out.auxiliaryAddedPoor10+=p.addedPoor10;out.auxiliaryAddedPoor25+=p.addedPoor25;
                    }
                }
                front.insert(cone);accepted=true;++out.seedCones;
                out.minimumConeMeanRatio=std::min(out.minimumConeMeanRatio,candidate.mean);
                out.minimumConeScaledJacobian=std::min(out.minimumConeScaledJacobian,candidate.jacobian);
                if(depth>1)++out.propagatedFaceAccepted;
                if(candidate.fresh){
                    result.addNode(candidate.point,region);++out.sites;
                    if(sourcePhase)++out.sourcePhaseSites;
                    out.effectiveScales.push_back(candidate.scale);
                    exclusionScales.push_back(candidate.requestedSize);depths.push_back(depth);
                    buckets[keyAt(candidate.point)].push_back(candidate.vertex);
                    out.maximumDepth=std::max(out.maximumDepth,depth);if(depth>1)++out.deepSites;
                }else ++out.reusedApices;
                for(unsigned side=0;side<4;++side){
                    const auto k=Front::key(bnd::tet_facet_detail::outwardFace(cone,side));
                    if(front.front.contains(k))propagationQueue.push_back(k);
                }
                break;
            }
            (void)accepted;
            if(front.front.empty()&&!sourcePhase)break;
            }catch(const SourceYield&){
                // Every charged auxiliary commit check precedes apply; the
                // front insert has no charged interruption after that commit.
                // Existing point/Plan guards have therefore restored staging.
                if(interruptedSource&&front.front.contains(*interruptedSource)){
                    queue.push_front(*interruptedSource);attemptedAt.erase(*interruptedSource);
                }
                sourceBudgetActive=false;sourceDeferred=true;
            }
        }
        out.frontClosed=front.front.empty();
        checkCancellation();
    }catch(const OcctInitialConstructionCancelled&){out.cancelled=true;
    }catch(const Limit& limit){out.budgetLimited=true;out.budgetReason=limit.reason;}
    out.remainingFrontFaces=front.front.size();
    if(out.sites>=siteCap&&!out.frontClosed){out.budgetLimited=true;
        if(std::string_view(out.budgetReason)=="none")out.budgetReason="site_limit";}
    return done();
}

} // namespace paacs::geo::detail
