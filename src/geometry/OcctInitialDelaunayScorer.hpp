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
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <new>
#include <optional>
#include <utility>
#include <vector>

namespace paacs::geo::detail {

class OcctInitialDelaunayScorer {
public:
    bool budgetLimited=false;
    std::uint64_t evaluations=0,fullCertifications=0,certificationRefusals=0;
    std::uint64_t scoreStateRejects=0,scoreSizeRejects=0,scoreMappingRejects=0;
    std::uint64_t scoreLocateMisses=0,scoreCavityUnavailable=0,scoreStorageRejects=0,scorePositiveRejects=0;
    std::uint64_t commitStateRejects=0,commitStorageRejects=0,commitPositiveRejects=0;
    std::uint64_t commitValidityRejects=0,commitSurvivorRejects=0;
    std::size_t cavityLimit=0,maximumPlannedCavity=0,maximumAcceptedCavity=0;
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
        std::function<bool(Vec3)> materialWitness={},std::size_t maximumCavityCells=0) {
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
        charge_=std::move(charge);materialWitness_=std::move(materialWitness);
        live_.walkConstraints=true;live_.maximumStoredCells=maximumStoredCells;
        live_.chargeWork=charge_;
        originalNodes_=ns;
        std::vector<Vec3> points(boundary.nodes().begin(),boundary.nodes().end());
        auto lo=points.front(),hi=lo;
        for(const auto p:points){
            tick();if(!finite(p))return false;
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
            tick();if(boundary.cell(i).type!=mesh::CellType::Tri3)return false;
            const auto* v=boundary.cellNodes(i);
            if(std::any_of(v,v+3,[&](Index id){return id<0||std::size_t(id)>=ns;}))return false;
            const Face f{Vertex(v[0]),Vertex(v[1]),Vertex(v[2])};
            if(!sourceFaces_.emplace(Complex::key(f),f).second)return false;
        }
        auto built=algo::buildTetCandidate(points);
        const auto& m=built.metrics;
        tick(m.construction.conflictInsphereTests+m.construction.conflictOrientationTests+
            m.validationFaceUses+m.validationHullTests+m.candidateCells+m.discardedCandidateWork);
        if(built.reason==algo::TetCandidateReason::AllocationFailed)throw std::bad_alloc{};
        if(!built.built()||built.candidate->points.size()!=points.size())return false;
        if(m.candidateCells>maximumStoredCells||m.constructionCandidateConservativeBytes>workingBytes){
            budgetLimited=true;return false;
        }
        live_.vertices.reserve(ns+1032);auxToFront_.reserve(ns+1032);
        relevance_.reserve(maximumStoredCells);
        for(std::size_t i=0;i<points.size();++i){
            tick();live_.vertices.push_back({points[i],std::uint64_t(i)+1});
            auxToFront_.push_back(i<ns?Vertex(i):bnd::kInvalidTetFacetVertex);
        }
        for(const auto& src:built.candidate->cells){
            Cell cell{{},0,region};
            for(unsigned j=0;j<4;++j)
                cell.vertices[j]=Vertex(built.candidate->pointToInput[src.vertices[j]]);
            if(!positive(cell))return false;
            live_.insert(cell);relevance_.push_back(0);
        }
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
            budgetLimited=true;++scoreStorageRejects;return {};
        }
        live_.vertices.push_back({site,std::uint64_t(staged)+1});
        try{
            for(const auto& cell:patch->added)if(!positive(cell)){
                ++scorePositiveRejects;live_.vertices.pop_back();return {};
            }
            Plan plan;plan.patch=std::move(*patch);plan.point=site;
            plan.generation=generation_;plan.stagedAux=staged;plan.valid=true;
            for(const auto id:plan.patch.removed){
                const auto state=oldRelevance(id,sourceFront);
                if(state==3){++plan.guardRemoved;continue;}
                if(state!=1){++plan.unclassifiedRemoved;continue;}
                ++plan.relevantRemoved;const auto q=shape(live_.cells[id]);
                plan.removedMinMean=std::min(plan.removedMinMean,q.first);
                plan.removedMinJac=std::min(plan.removedMinJac,q.second);
                plan.removedPoor10+=q.first<.1;plan.removedPoor25+=q.first<.25;
            }
            for(const auto& cell:plan.patch.added){
                const auto state=newRelevance(cell,sourceFront,staged,stagedFront);
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
        if(!initialized_||!plan.valid||plan.generation!=generation_||
            plan.stagedAux!=live_.vertices.size()||frontPointId!=live_.vertices.size()-8||
            plan.addedRelevance.size()!=plan.patch.added.size()){
            ++commitStateRejects;return false;
        }
        if(plan.patch.added.size()>live_.maximumStoredCells-live_.cells.size()){
            budgetLimited=true;++commitStorageRejects;return false;
        }
        live_.vertices.push_back({plan.point,std::uint64_t(plan.stagedAux)+1});
        bool applied=false;
        try{
            auto patch=plan.patch;
            ++fullCertifications;
            for(const auto& cell:patch.added)if(!positive(cell)){
                ++certificationRefusals;++commitPositiveRejects;
                live_.vertices.pop_back();return false;
            }
            if(!live_.valid(patch,std::nullopt,plan.stagedAux)){
                ++certificationRefusals;++commitValidityRejects;
                live_.vertices.pop_back();return false;
            }
            if(!live_.preservesSurvivors(patch)){
                ++certificationRefusals;++commitSurvivorRejects;
                live_.vertices.pop_back();return false;
            }
            tick(11*(patch.removed.size()+patch.added.size()));
            auto work=std::move(live_.chargeWork);
            live_.chargeWork={};
            try{live_.apply(patch);applied=true;live_.chargeWork=std::move(work);}
            catch(...){live_.chargeWork=std::move(work);throw;}
            auxToFront_.push_back(frontPointId);
            relevance_.insert(relevance_.end(),plan.addedRelevance.begin(),plan.addedRelevance.end());
            maximumAcceptedCavity=std::max(maximumAcceptedCavity,patch.removed.size());
            ++generation_;return true;
        }catch(...){if(!applied)live_.vertices.pop_back();throw;}
    }

    bool initialized()const noexcept{return initialized_;}
    std::size_t storedCells()const noexcept{return live_.cells.size();}
    std::size_t acceptedSites()const noexcept{return std::size_t(generation_);}
    const std::vector<Vertex>& auxiliaryToFront()const noexcept{return auxToFront_;}

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
    bool initialized_=false;
    void tick(std::uint64_t n=1)const{if(charge_)charge_(n);}
    static bool finite(Vec3 p){return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);}
    static bool equal(Vec3 a,Vec3 b){return a.x==b.x&&a.y==b.y&&a.z==b.z;}
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
