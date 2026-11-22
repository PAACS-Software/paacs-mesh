// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Private persistent unprotected Delaunay geometry for initial
// point proposals. This is not a constrained material mesh or a publication.
#pragma once

#include "OcctDelaunayCavity.hpp"
#include "../bnd/TetPlcGiftWrapping.hpp"
#include "paacs/algo/TetCandidate.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/mesh/Quality.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <new>
#include <optional>
#include <set>
#include <type_traits>
#include <utility>
#include <unordered_set>
#include <vector>

namespace paacs::geo::detail {

struct OcctInitialConstructionCancelled {};

class OcctInitialDelaunayScorer {
public:
    OcctInitialDelaunayScorer()=default;
    OcctInitialDelaunayScorer(const OcctInitialDelaunayScorer&)=delete;
    OcctInitialDelaunayScorer& operator=(const OcctInitialDelaunayScorer&)=delete;
    OcctInitialDelaunayScorer(OcctInitialDelaunayScorer&&)=delete;
    OcctInitialDelaunayScorer& operator=(OcctInitialDelaunayScorer&&)=delete;
    bool budgetLimited=false;
    mutable bool cancelled=false;
    bool compactionEnabled=true;
    // Only an owner preserving validated bootstrap/commit/compaction embedding
    // may enable the local certificate. Generic direct users retain full proof.
    bool useCertifiedConeCavities=false;
    std::uint64_t coneCertificateAttempts=0,coneCertificates=0,coneCertificateFallbacks=0;
    std::uint64_t bootstrapWork=0,bootstrapPeakTemporaryBytes=0;
    std::uint64_t bootstrapWorkLimit=0,bootstrapTemporaryLimitBytes=0;
    std::uint64_t evaluations=0,fullCertifications=0,certificationRefusals=0;
    std::uint64_t scoreStateRejects=0,scoreSizeRejects=0,scoreMappingRejects=0;
    std::uint64_t scoreLocateMisses=0,scoreCavityUnavailable=0,scoreStorageRejects=0,scorePositiveRejects=0;
    std::uint64_t commitStateRejects=0,commitStorageRejects=0,commitPositiveRejects=0;
    std::uint64_t commitValidityRejects=0,commitSurvivorRejects=0;
    std::size_t cavityLimit=0,maximumPlannedCavity=0,maximumAcceptedCavity=0;
    std::uint64_t compactionAttempts=0,compactions=0,deadCellsReclaimed=0,compactionWork=0;
    std::uint64_t maximumCompactionModeledBytes=0;
    std::size_t peakStoredCells=0;
    using Complex=bnd::detail::BistellarComplex;
    using Front=bnd::detail::PlcGiftWrapping;
    using Vertex=bnd::TetFacetVertexId;
    using Cell=bnd::TetFacetCell;
    using Face=Complex::Face;
    struct Plan {
        Complex::Proposal patch;
        Vec3 point{};
        std::uint64_t generation=0;
        Vertex stagedAux=bnd::kInvalidTetFacetVertex;
        bool valid=false,certified=false,available=false,comparisonAvailable=false,qualityAdmissible=true;
        bool materialInherited=false;
        double addedMinMean=1.,addedMinJac=1.,removedMinMean=1.,removedMinJac=1.,score=0.;
        std::size_t removedPoor10=0,removedPoor25=0,addedPoor10=0,addedPoor25=0;
        std::size_t relevantRemoved=0,relevantAdded=0,guardRemoved=0,guardAdded=0;
        std::size_t unclassifiedRemoved=0,unclassifiedAdded=0;
        // New cells become known material only after the caller admits this
        // same point natively and its temporary occupancy cone geometrically.
        std::vector<unsigned char> addedRelevance;
    };

    bool initialize(const mesh::Mesh& boundary,mesh::Classification region,
        std::size_t maximumStoredCells,std::uint64_t workingBytes,
        std::function<void(std::uint64_t)> charge,
        std::function<bool(Vec3)> materialWitness={},std::size_t maximumCavityCells=0,
        std::uint64_t maximumBootstrapWork=4'000'000,
        const std::atomic_bool* cancellation=nullptr) {
        cancellation_=cancellation;checkCancellation();
        if(initialized_||region.dim!=3||boundary.nodeCount()<4||!maximumStoredCells)return false;
        // Zero means a complete connected conflict cavity, bounded by the
        // admitted numerical complex and charged work. A nonzero cap is an
        // explicit private diagnostic/test limit. Never truncate a cavity.
        cavityLimit=maximumCavityCells;
        const auto ns=std::size_t(boundary.nodeCount());
        if(ns+1032>=bnd::kInvalidTetFacetVertex)return false;
        const auto fixed=8ULL*1024*1024+std::uint64_t(ns+1032)*512;
        if(fixed>=workingBytes){budgetLimited=true;return false;}
        maximumStoredCells=std::min<std::size_t>(maximumStoredCells,
            std::size_t((workingBytes-fixed)/8192));
        if(!maximumStoredCells){budgetLimited=true;return false;}
        workingBytes_=workingBytes;fixedBytes_=fixed;
        charge_=std::move(charge);materialWitness_=std::move(materialWitness);
        auto remainingBootstrapWork=maximumBootstrapWork;
        const auto accountInitialization=[&](std::uint64_t n=1){
            if(n>remainingBootstrapWork){budgetLimited=true;return false;}
            // Accounting a completed bounded attempt must not forget its work
            // when its cancellation token has since been reset by a caller.
            if(charge_)charge_(n);
            remainingBootstrapWork-=n;return true;
        };
        const auto initializationTick=[&](std::uint64_t n=1){
            checkCancellation();return accountInitialization(n);
        };
        live_.walkConstraints=true;live_.maximumStoredCells=maximumStoredCells;
        live_.chargeWork=[this](std::uint64_t n){tick(n);};
        originalNodes_=ns;
        std::vector<Vec3> points(boundary.nodes().begin(),boundary.nodes().end());
        auto lo=points.front(),hi=lo;
        for(const auto p:points){
            if(!initializationTick())return false;
            if(!finite(p))return false;
            lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};
            hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};
        }
        const auto span=std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z});
        if(!(span>0.)||!std::isfinite(span))return false;
        for(unsigned i=0;i<8;++i){
            const Vec3 p{(i&1)?hi.x+span:lo.x-span,(i&2)?hi.y+span:lo.y-span,
                (i&4)?hi.z+span:lo.z-span};
            if(!finite(p))return false;points.push_back(p);
        }
        for(Index i=0;i<boundary.cellCount();++i){
            if(!initializationTick())return false;
            if(boundary.cell(i).type!=mesh::CellType::Tri3)return false;
            const auto* v=boundary.cellNodes(i);
            if(std::any_of(v,v+3,[&](Index id){return id<0||std::size_t(id)>=ns;}))return false;
            const Face f{Vertex(v[0]),Vertex(v[1]),Vertex(v[2])};
            if(!sourceFaces_.emplace(Complex::key(f),f).second)return false;
        }
        // Reserve exact mapping/positive-import accounting before granting the
        // remaining work to construction. The online builder covers all of its
        // retries, canonicalization and validation; never run an unbounded retry.
        const auto importReserve=std::uint64_t(points.size())+12*std::uint64_t(maximumStoredCells);
        if(importReserve>=remainingBootstrapWork){budgetLimited=true;return false;}
        algo::TetCandidateBuildLimits buildLimits;
        buildLimits.maxWork=remainingBootstrapWork-importReserve;
        buildLimits.maxCells=maximumStoredCells;
        // Caller input points and source/mapping envelopes are already retained
        // in fixedBytes_; the core's temporary model includes its own workspace.
        buildLimits.maxTemporaryBytes=workingBytes-fixed;
        buildLimits.cancellation=cancellation_;
        bootstrapWorkLimit=buildLimits.maxWork;bootstrapTemporaryLimitBytes=buildLimits.maxTemporaryBytes;
        auto built=algo::buildTetCandidate(points,buildLimits);
        const auto& m=built.metrics;
        bootstrapWork=m.boundedWork;bootstrapPeakTemporaryBytes=m.boundedPeakTemporaryBytes;
        if(built.reason==algo::TetCandidateReason::Cancelled)cancelled=true;
        if(!accountInitialization(m.boundedWork))return false;
        if(cancelled)throw OcctInitialConstructionCancelled{};
        if(built.reason==algo::TetCandidateReason::BudgetExceeded){budgetLimited=true;return false;}
        if(built.reason==algo::TetCandidateReason::AllocationFailed)throw std::bad_alloc{};
        if(!built.built()||built.candidate->points.size()!=points.size())return false;
        // Construction scratch has gone away, but the returned candidate and
        // imported numerical complex overlap until this method returns.
        const auto retainedImportBytes=4096*std::uint64_t(maximumStoredCells);
        if(m.candidateCells>maximumStoredCells||
            m.candidateAccountedBytes>buildLimits.maxTemporaryBytes-retainedImportBytes){
            budgetLimited=true;return false;
        }
        live_.vertices.reserve(ns+1032);auxToFront_.reserve(ns+1032);
        relevance_.reserve(maximumStoredCells);
        for(std::size_t i=0;i<points.size();++i){
            if(!initializationTick())return false;
            live_.vertices.push_back({points[i],std::uint64_t(i)+1});
            auxToFront_.push_back(i<ns?Vertex(i):bnd::kInvalidTetFacetVertex);
        }
        for(const auto& src:built.candidate->cells){
            Cell cell{{},0,region};
            for(unsigned j=0;j<4;++j)
                cell.vertices[j]=Vertex(built.candidate->pointToInput[src.vertices[j]]);
            if(!initializationTick(12))return false;
            const auto& v=cell.vertices;
            if(predicates::orient3dFiltered(live_.vertices[v[0]].position,live_.vertices[v[1]].position,
                live_.vertices[v[2]].position,live_.vertices[v[3]].position)!=predicates::Sign::Positive)return false;
            live_.insert(cell,false);relevance_.push_back(0);
        }
        activeCells_=live_.cells.size();peakStoredCells=live_.cells.size();
        initialized_=true;return true;
    }

    // The caller temporarily appends site to sourceFront.vertices before this
    // call. All accepted earlier points have exactly the recorded identities.
    // No auxiliary topology changes occur while candidates are compared. The
    // connected all-positive proposal is ranked here; complete exact patch
    // and survivor admission is performed once for the selected plan at commit.
    std::optional<Plan> score(Vec3 site,Front& sourceFront) {
        ++evaluations;
        if(!initialized_||!finite(site)){++scoreStateRejects;return {};}
        if(live_.vertices.size()<8||sourceFront.vertices.size()+7!=live_.vertices.size()){
            ++scoreSizeRejects;return {};
        }
        const auto staged=Vertex(live_.vertices.size());
        const auto stagedFront=Vertex(sourceFront.vertices.size()-1);
        if(stagedFront!=staged-8||!equal(site,sourceFront.vertices.back().position)){
            ++scoreMappingRejects;return {};
        }
        const auto containing=locate(site);
        if(containing==Complex::noCell){++scoreLocateMisses;return {};}
        std::size_t owner=Complex::noCell;
        // Stored cells include tombstones, so this is a conservative upper
        // bound on every complete live conflict cavity without an arbitrary
        // geometric cutoff. Shared work/storage checks still bound the attempt.
        const auto maximumCells=cavityLimit?std::min(cavityLimit,live_.cells.size()):live_.cells.size();
        auto patch=buildOcctDelaunayCavity(live_,containing,site,maximumCells,owner,true,0.);
        if(!patch){++scoreCavityUnavailable;return {};}
        maximumPlannedCavity=std::max(maximumPlannedCavity,patch->removed.size());
        if(patch->added.size()>live_.maximumStoredCells-live_.cells.size()){
            if(compactionEnabled&&live_.cells.size()>activeCells_)compactionRequested_=true;else budgetLimited=true;
            ++scoreStorageRejects;return {};
        }
        live_.vertices.push_back({site,std::uint64_t(staged)+1});
        try{
            for(const auto& cell:patch->added)if(!positive(cell)){
                ++scorePositiveRejects;live_.vertices.pop_back();return {};
            }
            Plan plan;plan.patch=std::move(*patch);plan.point=site;
            plan.generation=planEpoch_;plan.stagedAux=staged;plan.valid=true;
            for(const auto id:plan.patch.removed){
                const auto state=oldRelevance(id,sourceFront);
                if(state==3){++plan.guardRemoved;continue;}
                if(state!=1){++plan.unclassifiedRemoved;continue;}
                ++plan.relevantRemoved;const auto q=shape(live_.cells[id]);
                plan.removedMinMean=std::min(plan.removedMinMean,q.first);
                plan.removedMinJac=std::min(plan.removedMinJac,q.second);
                plan.removedPoor10+=q.first<.1;plan.removedPoor25+=q.first<.25;
            }
            // Positive cones over the complete cavity shell replace the same
            // material union. Internal original source facets must also survive;
            // otherwise retain the complete per-cell source-contact checks.
            plan.materialInherited=!plan.patch.removed.empty()&&
                plan.relevantRemoved==plan.patch.removed.size()&&!plan.guardRemoved&&
                !plan.unclassifiedRemoved&&materialInheritanceAdmissible(plan.patch);
            for(const auto& cell:plan.patch.added){
                const unsigned char state=plan.materialInherited?1:newRelevance(cell,sourceFront,staged,stagedFront);
                plan.addedRelevance.push_back(state);
                if(state==3){++plan.guardAdded;continue;}
                if(state!=1){++plan.unclassifiedAdded;continue;}
                ++plan.relevantAdded;const auto q=shape(cell);
                plan.addedMinMean=std::min(plan.addedMinMean,q.first);
                plan.addedMinJac=std::min(plan.addedMinJac,q.second);
                plan.addedPoor10+=q.first<.1;plan.addedPoor25+=q.first<.25;
            }
            // A known added material subset still supplies useful actual-DT
            // shape information when other cells cross or lie outside the PLC.
            // Comparing whole-patch poor counts requires both subsets to be
            // classified without those missing material contributions.
            plan.available=plan.relevantAdded>0;
            plan.comparisonAvailable=plan.available&&plan.relevantRemoved&&
                !plan.unclassifiedRemoved&&!plan.unclassifiedAdded;
            if(plan.available){
                plan.qualityAdmissible=plan.addedMinMean>=.01&&plan.addedMinJac>=1e-6;
                plan.score=std::sqrt(plan.addedMinMean*plan.addedMinJac);
            }
            if(plan.comparisonAvailable)plan.qualityAdmissible=plan.qualityAdmissible&&
                plan.addedMinMean+1e-12>=std::min(plan.removedMinMean,.25)&&
                plan.addedPoor10<=plan.removedPoor10&&plan.addedPoor25<=plan.removedPoor25;
            live_.vertices.pop_back();return plan;
        }catch(...){live_.vertices.pop_back();throw;}
    }

    // Call after sourceFront.valid/native admission and before sourceFront's
    // insert. Any bounded-work refusal precedes auxiliary topology mutation.
    bool commit(const Plan& plan,Vertex frontPointId) {
        if(!initialized_||!plan.valid||!finite(plan.point)||plan.generation!=planEpoch_||
            plan.stagedAux!=live_.vertices.size()||frontPointId!=live_.vertices.size()-8||
            plan.addedRelevance.size()!=plan.patch.added.size()){
            ++commitStateRejects;return false;
        }
        if(plan.patch.added.size()>live_.maximumStoredCells-live_.cells.size()){
            if(compactionEnabled&&live_.cells.size()>activeCells_)compactionRequested_=true;else budgetLimited=true;
            ++commitStorageRejects;return false;
        }
        // Plans are mutable values. Recheck the actual removed set before its
        // material classification can be inherited by the certified replacement.
        if(plan.materialInherited&&!materialInheritanceAdmissible(plan.patch)){
            ++commitStateRejects;return false;
        }
        live_.vertices.push_back({plan.point,std::uint64_t(plan.stagedAux)+1});
        bool applied=false;
        try{
            auto patch=plan.patch;
            // This planner emits cone patches, never the separate bistellar
            // bipyramid certificate that bypasses pair/survivor admission.
            patch.convexBipyramid=false;
            ++fullCertifications;
            for(const auto& cell:patch.added)if(std::ranges::any_of(cell.vertices,[&](Vertex v){
                return v>=live_.vertices.size()||!finite(live_.vertices[v].position);
            })||!positive(cell)){
                ++certificationRefusals;++commitPositiveRejects;
                live_.vertices.pop_back();return false;
            }
            if(!live_.valid(patch,std::nullopt,plan.stagedAux)){
                ++certificationRefusals;++commitValidityRejects;
                live_.vertices.pop_back();return false;
            }
            bool coneCertified=false;
            if(useCertifiedConeCavities){
                ++coneCertificateAttempts;
                coneCertified=certifiesConeCavity(patch,plan.stagedAux);
                if(coneCertified)++coneCertificates;else ++coneCertificateFallbacks;
            }
            if(!coneCertified&&!live_.preservesSurvivors(patch)){
                ++certificationRefusals;++commitSurvivorRejects;
                live_.vertices.pop_back();return false;
            }
            tick(11*(patch.removed.size()+patch.added.size()));
            auto work=std::move(live_.chargeWork);
            live_.chargeWork={};
            try{live_.apply(patch);applied=true;live_.chargeWork=std::move(work);}
            catch(...){live_.chargeWork=std::move(work);throw;}
            auxToFront_.push_back(frontPointId);
            if(plan.materialInherited)relevance_.insert(relevance_.end(),patch.added.size(),1);
            else relevance_.insert(relevance_.end(),plan.addedRelevance.begin(),plan.addedRelevance.end());
            activeCells_=activeCells_-patch.removed.size()+patch.added.size();
            peakStoredCells=std::max(peakStoredCells,live_.cells.size());
            maximumAcceptedCavity=std::max(maximumAcceptedCavity,patch.removed.size());
            ++generation_;++planEpoch_;return true;
        }catch(...){if(!applied)live_.vertices.pop_back();throw;}
    }

    bool initialized()const noexcept{return initialized_;}
    std::size_t storedCells()const noexcept{return live_.cells.size();}
    std::size_t liveCells()const noexcept{return activeCells_;}
    std::uint64_t workingBytes()const noexcept{return workingBytes_;}
    std::size_t storageCapacity()const noexcept{return live_.maximumStoredCells;}
    std::size_t acceptedSites()const noexcept{return std::size_t(generation_);}
    const std::vector<Vertex>& auxiliaryToFront()const noexcept{return auxToFront_;}
    const Complex& complexForDiagnostics()const noexcept{return live_;}
    const std::vector<unsigned char>& relevanceForDiagnostics()const noexcept{return relevance_;}
    std::uint64_t planEpoch()const noexcept{return planEpoch_;}
    bool compactionRequested()const noexcept {
        return compactionEnabled&&!compactionBlocked_&&compactionRequested_&&live_.cells.size()>activeCells_;
    }

    bool compactionPending()const noexcept {
        if(!compactionEnabled||!initialized_||compactionBlocked_||live_.cells.size()<=activeCells_)return false;
        const auto stored=live_.cells.size(),dead=stored-activeCells_;
        return compactionRequested_||
            (stored>live_.maximumStoredCells*3/4&&dead>=std::max<std::size_t>(64,stored/8));
    }

    // Caller runs this between face batches, before any candidate Plan exists.
    // Force is a private diagnostic/test control; it never raises capacity.
    bool prepareForBatch(bool force=false) {
        if(!compactionEnabled||!initialized_||compactionBlocked_||live_.cells.size()<=activeCells_||
            (!force&&!compactionPending()))return false;
        ++compactionAttempts;
        constexpr std::uint64_t perComplexCell=4096;
        using FaceValue=std::pair<const Face,std::set<std::size_t>>;
        using EdgeValue=std::pair<const Complex::Edge,std::set<std::size_t>>;
        // Modeled payload, not process RSS. This unprotected auxiliary complex
        // has at most 4 face keys, 6 edge keys and 10 owner-set entries per cell.
        // Include 2x vector capacity, 4 BVH nodes/order slots, four pointer links
        // and 64 allocator/alignment bytes per tree allocation, plus scratch.
        constexpr auto modeledCell=2*(sizeof(Cell)+sizeof(unsigned char)+sizeof(Complex::Box)+
            sizeof(std::array<std::size_t,4>))+4*(sizeof(Complex::Box)+4*sizeof(std::size_t))+
            4*(sizeof(FaceValue)+4*sizeof(void*)+64)+6*(sizeof(EdgeValue)+4*sizeof(void*)+64)+
            10*(sizeof(std::size_t)+4*sizeof(void*)+64)+512;
        static_assert(modeledCell<=perComplexCell,"auxiliary compaction payload model");
        static_assert(std::is_nothrow_move_assignable_v<Complex>,"atomic auxiliary compaction move");
        const auto stored=live_.cells.size(),dead=stored-activeCells_;
        // The admitted 8192/cell envelope reserves 4096/cell for each complex.
        // Old capacity remains reserved while the staged live-only copy exists.
        const auto peak=fixedBytes_+perComplexCell*(live_.maximumStoredCells+activeCells_);
        if(peak>workingBytes_){budgetLimited=true;compactionBlocked_=true;return false;}
        maximumCompactionModeledBytes=std::max(maximumCompactionModeledBytes,peak);
        if(!live_.protectedFaces.empty()||!live_.protectedEdges.empty()||
            !live_.requiredFaces.empty()||!live_.requiredEdges.empty())
            throw std::logic_error("auxiliary compaction requires unprotected numerical complex");
        const auto chargeCompaction=[&](std::uint64_t n=1){tick(n);compactionWork+=n;};
        Complex next;
        next.maximumStoredCells=live_.maximumStoredCells;
        next.useLookahead=live_.useLookahead;next.walkConstraints=live_.walkConstraints;next.profile=live_.profile;
        next.metrics=live_.metrics;next.edgeWalks=live_.edgeWalks;
        next.facetWalks=live_.facetWalks;next.walkFallbacks=live_.walkFallbacks;
        next.chargeWork=chargeCompaction;
        chargeCompaction(live_.vertices.size());next.vertices=live_.vertices;
        next.cells.reserve(activeCells_);next.alive.reserve(activeCells_);
        next.cellBoxes.reserve(activeCells_);next.neighbours.reserve(activeCells_);
        std::vector<unsigned char> relevance;relevance.reserve(activeCells_);
        for(std::size_t id=0;id<stored;++id){
            chargeCompaction();if(!live_.alive[id])continue;
            // insert renumbers cell sequences/IDs but keeps the exact ordered
            // vertex quartet and classification. Point IDs/sequences never move.
            next.insert(live_.cells[id]);relevance.push_back(relevance_[id]);
        }
        if(next.cells.size()!=activeCells_||relevance.size()!=activeCells_)
            throw std::logic_error("auxiliary compaction live count mismatch");
        next.index.rebuild(next.cells.size(),[&](std::size_t id){return next.cellBoxes[id];},chargeCompaction);
        chargeCompaction(11*stored); // Retire old cell/index/owner allocations.
        next.chargeWork=[this](std::uint64_t n){tick(n);}; // No staging capture may survive.
        // Every allocation and charged operation precedes this noexcept commit.
        live_=std::move(next);relevance_.swap(relevance);
        ++planEpoch_;++compactions;deadCellsReclaimed+=dead;compactionRequested_=false;
        return true;
    }

private:
    Complex live_;
    std::function<void(std::uint64_t)> charge_;
    std::function<bool(Vec3)> materialWitness_;
    std::map<Face,Face> sourceFaces_;
    std::vector<Vertex> auxToFront_;
    // 0 unqueried, 1 material, 2 outside/ambiguous, 3 guard, 4 crossing source.
    std::vector<unsigned char> relevance_;
    std::size_t originalNodes_=0;
    std::uint64_t generation_=0;
    std::uint64_t planEpoch_=0;
    std::uint64_t workingBytes_=0,fixedBytes_=0;
    const std::atomic_bool* cancellation_=nullptr;
    std::size_t activeCells_=0;
    bool compactionRequested_=false,compactionBlocked_=false;
    bool initialized_=false;
    void checkCancellation()const{
        if(cancelled||(cancellation_&&cancellation_->load())){
            cancelled=true;throw OcctInitialConstructionCancelled{};
        }
    }
    void tick(std::uint64_t n=1)const{checkCancellation();if(charge_)charge_(n);}
    static bool finite(Vec3 p){return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);}
    static bool equal(Vec3 a,Vec3 b){return a.x==b.x&&a.y==b.y&&a.z==b.z;}
    bool materialInheritanceAdmissible(const Complex::Proposal& patch)const{
        if(patch.removed.empty())return false;
        std::unordered_set<std::size_t> removed;
        tick(patch.removed.size());removed.reserve(patch.removed.size());
        for(const auto id:patch.removed){
            tick();
            if(id>=live_.alive.size()||!live_.alive[id]||id>=relevance_.size()||
                relevance_[id]!=1||!removed.insert(id).second)return false;
        }
        // Two individually admissible old cells can touch opposite sides of an
        // original source facet. Erasing that internal facet would not preserve
        // its source-contact barrier, so inheritance must fall back to contacts.
        for(const auto id:patch.removed)for(unsigned side=0;side<4;++side){
            tick();const auto face=Complex::key(bnd::tet_facet_detail::outwardFace(live_.cells[id],side));
            if(!sourceFaces_.contains(face))continue;
            const auto owner=live_.faces.find(face);if(owner==live_.faces.end())return false;
            unsigned oldOwners=0;
            for(const auto cell:owner->second){tick();oldOwners+=removed.contains(cell);}
            if(oldOwners>1)return false;
        }
        return true;
    }
    bool certifiesConeCavity(const Complex::Proposal& patch,Vertex staged)const{
        if(!live_.walkConstraints||patch.removed.empty()||patch.added.empty()||
            staged+std::size_t(1)!=live_.vertices.size()||!finite(live_.vertices[staged].position)||
            !live_.protectedFaces.empty()||!live_.protectedEdges.empty()||
            !live_.requiredFaces.empty()||!live_.requiredEdges.empty())return false;
        struct Use {std::size_t cell;unsigned side;Face directed;};
        struct Uses {std::array<Use,2> owners{};unsigned count=0;bool cone=false;};
        using Edge=Complex::Edge;
        struct EdgeUse {std::size_t face;unsigned side;bool forward;};
        struct EdgeUses {std::array<EdgeUse,2> owners{};unsigned count=0;};
        // Modeled payload, not RSS. Let L be stored cells and R unused capacity.
        // Storage admission already guarantees C removed<=L and S added<=R.
        // Bound faces by4C, old vertices by4C, shell edges/vertices by3S (even
        // before pairing), and include 64 allocator bytes per node/vector,
        // 3x incident-vector payload for growth overlap, selected patch copy,
        // and three retained fresh plans at2x vector capacity. The weighted
        // bound is modeledRemoved*L + modeledAdded*R, with L+R=capacity.
        constexpr auto modeledRemoved=
            4*(sizeof(std::pair<const Face,Uses>)+4*sizeof(void*)+64)+
            4*(sizeof(Vertex)+5*sizeof(void*)+64)+
            (sizeof(std::size_t)+6*sizeof(void*)+64)+4*sizeof(Face)+7*sizeof(std::size_t);
        constexpr auto modeledAdded=
            3*(sizeof(std::pair<const Edge,EdgeUses>)+4*sizeof(void*)+64)+
            3*(sizeof(std::pair<const Vertex,std::vector<std::size_t>>)+4*sizeof(void*)+128)+
            11*sizeof(std::size_t)+sizeof(std::array<std::size_t,3>)+sizeof(Cell)+
            6*(sizeof(Cell)+sizeof(unsigned char));
        static_assert(modeledRemoved<=4096&&modeledAdded<=4096,"local cone certificate payload envelope");
        // Full valid() scratch has returned; no compaction overlaps a plan.
        // The existing 8192/capacity admission reserves 4096 for the live
        // complex and 4096 for this scratch; fixed container heads use the
        // existing fixed reserve. No new storage limit or uncharged reserve.
        std::unordered_set<std::size_t> removed;
        std::map<Face,Uses> faces;
        std::set<Vertex> oldVertices;
        tick(patch.removed.size());removed.reserve(patch.removed.size());
        for(const auto id:patch.removed){
            tick();
            if(id>=live_.alive.size()||!live_.alive[id]||id>=live_.neighbours.size()||
                !removed.insert(id).second)return false;
            const auto& cell=live_.cells[id];
            for(const auto v:cell.vertices){tick();
                if(v>=staged||!finite(live_.vertices[v].position))return false;
                oldVertices.insert(v);
            }
            if(!positive(cell))return false;
        }
        bool pointInRemoved=false;
        const auto point=live_.vertices[staged].position;
        for(const auto id:patch.removed){
            bool pointInCell=!pointInRemoved;
            for(unsigned side=0;side<4;++side){
                tick();const auto directed=bnd::tet_facet_detail::outwardFace(live_.cells[id],side);
                auto& uses=faces[Complex::key(directed)];
                if(uses.count==2)return false;
                uses.owners[uses.count++]={id,side,directed};
                if(pointInCell){tick();
                    if(predicates::orient3dFiltered(live_.vertices[directed[0]].position,
                        live_.vertices[directed[1]].position,live_.vertices[directed[2]].position,point)
                        ==predicates::Sign::Positive)pointInCell=false;
                }
            }
            pointInRemoved=pointInRemoved||pointInCell;
        }
        if(!pointInRemoved)return false;
        std::vector<Face> shell;
        const auto shellSlots=4*patch.removed.size();
        tick(shellSlots);shell.reserve(shellSlots);
        for(const auto& [key,uses]:faces){
            tick();const auto found=live_.faces.find(key);
            if(found==live_.faces.end()||found->second.empty()||found->second.size()>2)return false;
            unsigned oldOwners=0;std::array<Face,2> directed{};unsigned at=0;
            for(const auto id:found->second){
                tick();
                if(id>=live_.alive.size()||!live_.alive[id]||id>=live_.neighbours.size())return false;
                const auto keys=Complex::faceKeys(live_.cells[id]);
                const auto face=std::ranges::find(keys,key);
                if(face==keys.end())return false;
                const auto side=unsigned(face-keys.begin());
                directed[at++]=bnd::tet_facet_detail::outwardFace(live_.cells[id],side);
                const auto other=found->second.size()==1?Complex::noCell:
                    (*found->second.begin()==id?*found->second.rbegin():*found->second.begin());
                if(live_.neighbours[id][side]!=other)return false;
                oldOwners+=removed.contains(id);
            }
            if(oldOwners!=uses.count)return false;
            if(at==2&&!sameCycle(directed[0],Face{directed[1][0],directed[1][2],directed[1][1]}))return false;
            for(unsigned i=0;i<uses.count;++i)
                if(!found->second.contains(uses.owners[i].cell))return false;
            if(uses.count==1)shell.push_back(uses.owners[0].directed);
        }
        if(shell.size()<4||shell.size()!=patch.added.size())return false;
        for(const auto& cell:patch.added){
            tick();if(cell.vertices[0]!=staged)return false;
            const Face base{cell.vertices[1],cell.vertices[2],cell.vertices[3]};
            const auto found=faces.find(Complex::key(base));
            if(found==faces.end()||found->second.count!=1||found->second.cone||
                base!=found->second.owners[0].directed||!positive(cell))return false;
            found->second.cone=true;
        }
        for(const auto& [key,uses]:faces){(void)key;tick();if(uses.count==1&&!uses.cone)return false;}
        std::map<Edge,EdgeUses> edges;
        std::map<Vertex,std::vector<std::size_t>> vertexFaces;
        for(std::size_t i=0;i<shell.size();++i)for(unsigned side=0;side<3;++side){
            tick(2);const auto a=shell[i][side],b=shell[i][(side+1)%3];
            auto& uses=edges[Complex::key(Edge{a,b})];
            if(uses.count==2)return false;
            uses.owners[uses.count++]={i,side,a<b};
            auto& incident=vertexFaces[a];
            if(incident.size()==incident.capacity()){
                const auto slots=std::min(shell.size(),std::max<std::size_t>(1,2*incident.size()));
                tick(incident.size()+slots);incident.reserve(slots);
            }
            incident.push_back(i);
        }
        if(oldVertices.size()!=vertexFaces.size())return false;
        for(const auto v:oldVertices){tick();if(!vertexFaces.contains(v))return false;}
        if(std::int64_t(vertexFaces.size())-std::int64_t(edges.size())+std::int64_t(shell.size())!=2)return false;
        tick(shell.size());std::vector<std::array<std::size_t,3>> adjacent(shell.size());
        for(const auto& [edge,uses]:edges){(void)edge;tick();
            if(uses.count!=2||uses.owners[0].forward==uses.owners[1].forward)return false;
            const auto a=uses.owners[0],b=uses.owners[1];
            adjacent[a.face][a.side]=b.face;adjacent[b.face][b.side]=a.face;
        }
        tick(2*shell.size());std::vector<std::size_t> queue;queue.reserve(shell.size());
        std::vector<std::size_t> seen(shell.size(),0);queue.push_back(0);seen[0]=1;
        for(std::size_t i=0;i<queue.size();++i)for(const auto next:adjacent[queue[i]]){
            tick();if(!seen[next]){seen[next]=1;queue.push_back(next);}
        }
        if(queue.size()!=shell.size())return false;
        // Edge incidence alone permits pinched vertices. Each vertex's incident
        // faces must form one cycle, so radial coning preserves exact simplex
        // contacts as well as positive-volume occupancy.
        tick(seen.size());std::fill(seen.begin(),seen.end(),0);std::size_t stamp=0;
        for(const auto& [vertex,incident]:vertexFaces){
            tick();++stamp;queue.clear();queue.push_back(incident.front());seen[incident.front()]=stamp;
            for(std::size_t i=0;i<queue.size();++i){
                const auto face=queue[i];const auto found=std::ranges::find(shell[face],vertex);
                if(found==shell[face].end())return false;
                const auto side=unsigned(found-shell[face].begin());
                for(const auto edgeSide:std::array{(side+2)%3,side}){
                    tick();const auto next=adjacent[face][edgeSide];
                    if(seen[next]!=stamp){seen[next]=stamp;queue.push_back(next);}
                }
            }
            if(queue.size()!=incident.size())return false;
        }
        // The trusted old complex is embedded. Its closed manifold shell is
        // strictly visible from the contained site (the exact positive cones).
        // Every ray therefore exits once: these cones form the same embedded,
        // face-to-face cavity and can touch survivors only on that old shell.
        return true;
    }
    bool guard(const Cell& cell)const{
        return std::ranges::any_of(cell.vertices,[&](Vertex v){return v>=originalNodes_&&v<originalNodes_+8;});
    }
    bool positive(const Cell& cell)const{
        tick();const auto& v=cell.vertices;
        return predicates::orient3dFiltered(live_.vertices[v[0]].position,live_.vertices[v[1]].position,
            live_.vertices[v[2]].position,live_.vertices[v[3]].position)==predicates::Sign::Positive;
    }
    std::pair<double,double> shape(const Cell& cell)const{
        tick(8);std::array<Vec3,4> p{};
        for(unsigned i=0;i<4;++i)p[i]=live_.vertices[cell.vertices[i]].position;
        const auto det=predicates::orient3dDeterminantEstimate(p[0],p[1],p[2],p[3]);
        double sum=0.;for(unsigned i=0;i<4;++i)for(unsigned j=i+1;j<4;++j)sum+=(p[i]-p[j]).norm2();
        const auto mean=det>0.&&sum>0.?12*std::pow(det/2.,2./3.)/sum:0.;
        const auto jac=mesh::tetScaledJacobian(p[0],p[1],p[2],p[3]);
        return {std::isfinite(mean)?mean:0.,std::isfinite(jac)?jac:0.};
    }
    std::size_t locate(Vec3 p){
        auto result=Complex::noCell;
        live_.visit({p,p},[&](std::size_t id){
            for(unsigned side=0;side<4;++side){
                tick();const auto f=bnd::tet_facet_detail::outwardFace(live_.cells[id],side);
                if(predicates::orient3dFiltered(live_.vertices[f[0]].position,live_.vertices[f[1]].position,
                    live_.vertices[f[2]].position,p)==predicates::Sign::Positive)return true;
            }result=id;return false;
        });return result;
    }
    std::optional<Cell> mapped(const Cell& cell,const Front& front,Vertex staged,
        Vertex stagedFront)const{
        Cell out=cell;
        for(unsigned i=0;i<4;++i){
            tick();const auto aux=cell.vertices[i];
            const auto id=aux==staged?stagedFront:(aux<auxToFront_.size()?auxToFront_[aux]:bnd::kInvalidTetFacetVertex);
            if(id==bnd::kInvalidTetFacetVertex||id>=front.vertices.size()||
                !equal(live_.vertices[aux].position,front.vertices[id].position))return {};
            out.vertices[i]=id;
        }return out;
    }
    static bool sameCycle(Face a,Face b){
        for(unsigned i=0;i<3;++i)if(a[0]==b[i]&&a[1]==b[(i+1)%3]&&a[2]==b[(i+2)%3])return true;
        return false;
    }
    unsigned char oldRelevance(std::size_t id,Front& front){
        if(relevance_[id])return relevance_[id];
        const auto& cell=live_.cells[id];
        if(guard(cell))return relevance_[id]=3;
        const auto mappedCell=mapped(cell,front,bnd::kInvalidTetFacetVertex,bnd::kInvalidTetFacetVertex);
        if(!mappedCell)return relevance_[id]=2;
        if(!front.boundaryAccepts(*mappedCell,front.box(mappedCell->vertices),0))return relevance_[id]=4;
        bool witnessed=false;
        for(auto v:mappedCell->vertices)witnessed=witnessed||v>=originalNodes_;
        for(unsigned side=0;side<4&&!witnessed;++side){
            const auto f=bnd::tet_facet_detail::outwardFace(*mappedCell,side);
            const auto it=sourceFaces_.find(Complex::key(f));
            witnessed=it!=sourceFaces_.end()&&sameCycle(f,it->second);
        }
        if(!witnessed&&materialWitness_){
            Vec3 center{};for(auto v:cell.vertices)center=center+live_.vertices[v].position*.25;
            witnessed=materialWitness_(center);
        }
        return relevance_[id]=witnessed?1:2;
    }
    unsigned char newRelevance(const Cell& cell,Front& front,Vertex staged,Vertex stagedFront){
        if(guard(cell))return 3;
        const auto mappedCell=mapped(cell,front,staged,stagedFront);
        if(!mappedCell)return 2;
        return front.boundaryAccepts(*mappedCell,front.box(mappedCell->vertices),0)?1:4;
    }
};

} // namespace paacs::geo::detail
