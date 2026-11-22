// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Private initial volume point construction, before constrained Delaunay.
// The temporary cones are exact noncrossing proposal occupancy only. They are
// never exported as final volume cells, nor treated as a coverage certificate.
#pragma once

#include "../bnd/TetPlcGiftWrapping.hpp"
#include "OcctInitialDelaunayScorer.hpp"
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
    // New nodes use min(requested h, 1.15*current-face median edge).
    // These construction scales allow a gradual transition away from a finer
    // immutable CAD surface; they are not a mutation of the requested field.
    std::vector<double> effectiveScales;
};

namespace occt_initial_point_front_detail {
inline constexpr double separationFactor=.65, faceGrowth=1.15;
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
    const auto siteCap=std::min<std::size_t>(limits.sites,1024);
    const auto coneCap=std::min<std::size_t>(limits.cones,4096);
    const auto ns=std::size_t(boundary.nodeCount()),nf=std::size_t(boundary.cellCount());
    const auto stateCap=std::min<std::size_t>(limits.frontStates,nf+3*coneCap);
    if(!siteCap||!coneCap||!limits.work||!limits.queries||nf>stateCap||
        ns+siteCap>=bnd::kInvalidTetFacetVertex) {
        out.budgetLimited=true;return done();
    }
    // Conservative admission for copied source/result, two AABB indices,
    // boundary BVH, node hash/scales (including one exclusion double per node),
    // occupied cones and all queue/front states.
    const auto payload=8ULL*1024*1024+std::uint64_t(ns+siteCap)*2048+
        std::uint64_t(nf)*4096+std::uint64_t(coneCap)*4096+
        std::uint64_t(stateCap)*4096;
    if(payload>limits.workingBytes){out.budgetLimited=true;return done();}
    struct Limit {};
    const auto account=[&](std::uint64_t n=1) {
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
        tick();if(out.nativeQueries>=limits.queries)throw Limit{};
        ++out.nativeQueries;
        if(!admitNode(p)){++out.nativeRejects;return false;}return true;
    };
    const auto materialWitness=[&](Vec3 p) {
        tick();if(out.nativeQueries>=limits.queries)throw Limit{};
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
    front.vertices.reserve(ns+siteCap);front.cells.reserve(coneCap);front.cellBoxes.reserve(coneCap);
    std::vector<std::size_t> depths;depths.reserve(ns+siteCap);
    std::vector<double> exclusionScales;exclusionScales.reserve(ns+siteCap);
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
        out.effectiveScales.reserve(ns+siteCap);
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
        if(limits.workingBytes<=payload+32ULL*1024*1024){out.budgetLimited=true;return done();}
        const auto auxiliaryBytes=limits.workingBytes-payload;
        const auto auxiliaryCells=std::min<std::uint64_t>(65'536,
            (auxiliaryBytes-32ULL*1024*1024)/8192);
        OcctInitialDelaunayScorer auxiliary;
        bool initialized=false;
        try{
            initialized=auxiliary.initialize(boundary,region,std::size_t(auxiliaryCells),auxiliaryBytes,account,
                materialWitness,limits.maximumDelaunayCavityCells,limits.work-out.workUsed,limits.cancellation);
        }catch(const OcctInitialConstructionCancelled&){
            out.auxiliaryBootstrapWork=auxiliary.bootstrapWork;
            out.auxiliaryBootstrapPeakTemporaryBytes=auxiliary.bootstrapPeakTemporaryBytes;
            out.cancelled=true;throw;
        }
        out.auxiliaryBootstrapWork=auxiliary.bootstrapWork;
        out.auxiliaryBootstrapPeakTemporaryBytes=auxiliary.bootstrapPeakTemporaryBytes;
        out.cancelled=out.cancelled||auxiliary.cancelled;
        if(!initialized){
            out.budgetLimited=true;return done();
        }
        auxiliary.compactionEnabled=limits.enableAuxiliaryCompaction;
        auxiliary.useCertifiedConeCavities=true;
        out.auxiliaryInitialized=true;
        const auto syncAuxiliary=[&](){
            out.budgetLimited=out.budgetLimited||auxiliary.budgetLimited;
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
        };
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
        while(!propagationQueue.empty()||!queue.empty()||front.cells.size()!=roundStart){
            // Prior candidate vectors have been destroyed. Compact only here,
            // before any point/Plan for the next face is staged.
            try{(void)auxiliary.prepareForBatch();}
            catch(...){syncAuxiliary();throw;}
            syncAuxiliary();
            if(propagationQueue.empty()&&queue.empty()){
                roundStart=front.cells.size();rebuild();
                if(propagationQueue.empty()&&queue.empty())break;
            }
            // Continue the live advancing front before activating another
            // untouched source face. Stalled propagation falls back to the
            // source queue; fixed budgets cannot be spent on a complete first
            // surface layer before any accepted point can advance inward.
            tick();const bool propagating=!propagationQueue.empty();
            const auto key=propagating?propagationQueue.front():queue.front();
            if(propagating)propagationQueue.pop_front();else queue.pop_front();
            const auto found=front.front.find(key);if(found==front.front.end())continue;
            if(const auto at=attemptedAt.find(key);at!=attemptedAt.end()&&at->second==front.cells.size())continue;
            if(attemptedAt.size()>=stateCap&&!attemptedAt.contains(key))throw Limit{};
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
                tick();if(out.candidates>=limits.candidates)throw Limit{};++out.candidates;
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
                    if(!plan||!plan->valid){++out.auxiliaryInvalidRejects;return;}
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
                    out.budgetLimited=true;continue;
                }
                if(compactionRetryAtSite.size()>=stateCap&&!compactionRetryAtSite.contains(key))throw Limit{};
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
                if(front.cells.size()>=coneCap||front.front.size()+3>stateCap)throw Limit{};
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
            if(front.front.empty())break;
        }
        out.frontClosed=front.front.empty();
        checkCancellation();
    }catch(const OcctInitialConstructionCancelled&){out.cancelled=true;
    }catch(const Limit&){out.budgetLimited=true;}
    out.remainingFrontFaces=front.front.size();
    if(out.sites>=siteCap&&!out.frontClosed)out.budgetLimited=true;
    return done();
}

} // namespace paacs::geo::detail
