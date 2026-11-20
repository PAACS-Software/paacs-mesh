// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Bounded, plain-value tetrahedralization of one closed repair region. Search
// changes four face-ledger rows per step; it never copies the surrounding mesh.
#include "TetPlcGiftWrapping.hpp"
#include "paacs/bnd/TetPlcCavityRecovery.hpp"

namespace paacs::bnd::detail {
struct PlcPatchSearchMetrics {
    std::uint64_t states=0,backtracks=0,candidates=0,pairCacheHits=0;
    bool bounded=false;
    bool stateBound=false;
};
class PlcPatchSearch {
    using Face=PlcGiftWrapping::Face;
    struct Candidate {TetFacetCell cell;PlcGiftWrapping::Box box;bool valid=false;};
    struct Use {unsigned count=0;int sense=0;};
    struct Change {Face key,directed;bool wasFront=false;};
    struct Bound {};
    PlcGiftWrapping geometry;
    const TetPlcPatchSearchOptions& limits;
    unsigned maximumCells;
    PlcPatchSearchMetrics& metrics;
    std::map<Face,int> boundary;
    std::map<Face,Face> front;
    std::map<Face,Use> uses;
    std::map<std::array<unsigned,4>,std::size_t> candidateIndex;
    std::set<std::array<unsigned,4>> rejectedGeometry;
    std::map<std::array<unsigned,4>,std::vector<std::size_t>> faceOptions;
    std::vector<Candidate> candidates;
    std::vector<std::size_t> selected;
    // Pair contacts depend only on immutable coordinates and identities.
    // Bound this cache independently of the number of search states.
    std::map<std::pair<std::size_t,std::size_t>,bool> pairCache;
    TetFacetComplexView view()const{return {geometry.vertices,{},{},{},{},{3,0}};}
    void tick(std::uint64_t n=1)const {geometry.charge(n);}
    int sense(Face f)const{return tet_facet_detail::faceSense(view(),f);}
    const std::vector<std::size_t>& options(Face directed){
        const auto k=PlcGiftWrapping::key(directed);
        const std::array<unsigned,4> optionKey{k[0],k[1],k[2],unsigned(sense(directed)>0)};
        if(!faceOptions.contains(optionKey)&&faceOptions.size()>=std::size_t(limits.maximumCandidateTets)*4)throw Bound{};
        auto [row,added]=faceOptions.try_emplace(optionKey);if(!added)return row->second;
        for(unsigned apex=0;apex<geometry.vertices.size();++apex){
            tick();if(std::find(k.begin(),k.end(),apex)!=k.end())continue;
            if(predicates::orient3d(geometry.vertices[directed[0]].position,geometry.vertices[directed[1]].position,geometry.vertices[directed[2]].position,geometry.vertices[apex].position)!=predicates::Sign::Negative)continue;
            std::array<unsigned,4> ids{k[0],k[1],k[2],apex};std::sort(ids.begin(),ids.end());
            auto known=candidateIndex.find(ids);std::size_t index;
            if(known==candidateIndex.end()){
                if(rejectedGeometry.contains(ids))continue;
                ++metrics.candidates;TetFacetCell cell;const bool valid=geometry.valid(PlcGiftWrapping::Front{directed},apex,cell);
                // Invalid proposals spend work, but cannot consume the entire
                // retained-candidate budget before a large refill can close.
                // The negative cache is independently bounded.
                if(!valid){if(rejectedGeometry.size()<limits.maximumCandidateTets)rejectedGeometry.insert(ids);continue;}
                if(candidates.size()>=limits.maximumCandidateTets)throw Bound{};
                index=candidates.size();candidateIndex.emplace(ids,index);
                candidates.push_back({cell,geometry.box(cell.vertices),true});
            }else index=known->second;
            if(candidates[index].valid)row->second.push_back(index);
        }
        return row->second;
    }
    bool ledgerAccepts(const TetFacetCell& cell)const {
        for(unsigned side=0;side<4;++side){
            tick();const auto f=tet_facet_detail::outwardFace(cell,side),k=PlcGiftWrapping::key(f);const int s=sense(f);
            const auto at=uses.find(k);const unsigned count=at==uses.end()?0:at->second.count;
            if(auto expected=boundary.find(k);expected!=boundary.end()){
                if(count||s!=expected->second)return false;
            }else if(count>=2||(count==1&&s==at->second.sense))return false;
            if(auto current=front.find(k);current!=front.end()&&s!=sense(current->second))return false;
        }
        return true;
    }
    bool pairsAccept(std::size_t index){
        const auto& a=candidates[index];
        for(auto old:selected){tick();const auto& b=candidates[old];if(!PlcGiftWrapping::overlap(a.box,b.box))continue;
            const auto pair=std::minmax(index,old);const std::pair<std::size_t,std::size_t> k{pair.first,pair.second};
            if(auto at=pairCache.find(k);at!=pairCache.end()){++metrics.pairCacheHits;if(!at->second)return false;continue;}
            const bool valid=!tet_facet_detail::tetImproperlyIntersects(view(),a.cell,b.cell);
            if(pairCache.size()<65536)pairCache.emplace(k,valid);
            if(!valid)return false;
        }
        return true;
    }
    std::array<Change,4> apply(const TetFacetCell& cell){
        std::array<Change,4> undo;
        for(unsigned side=0;side<4;++side){
            auto f=tet_facet_detail::outwardFace(cell,side);const auto k=PlcGiftWrapping::key(f);
            auto& use=uses[k];if(!use.count)use.sense=sense(f);++use.count;
            if(auto at=front.find(k);at!=front.end()){undo[side]={k,at->second,true};front.erase(at);}
            else {std::swap(f[1],f[2]);undo[side]={k,f,false};front.emplace(k,f);}
        }
        return undo;
    }
    void undo(const std::array<Change,4>& changes){
        for(const auto& change:changes){auto at=uses.find(change.key);if(!--at->second.count)uses.erase(at);
            if(change.wasFront)front.emplace(change.key,change.directed);else front.erase(change.key);}
    }
    bool search(){
        // Keep the same depth-first choices and four-row undo records on a
        // bounded heap stack. Large valid refills need not consume the native
        // thread stack or inherit the old 512-tetrahedron recursion ceiling.
        struct Frame {std::vector<std::size_t> choices;std::size_t next=0;std::array<Change,4> changes{};bool entered=false,parent=false;};
        std::vector<Frame> stack;stack.reserve(std::size_t(maximumCells)+1);stack.emplace_back();
        const auto backtrack=[&](){const auto changes=stack.back().changes;const bool parent=stack.back().parent;stack.pop_back();if(parent){selected.pop_back();undo(changes);++metrics.backtracks;}};
        while(!stack.empty()){
            auto& frame=stack.back();
            if(!frame.entered){
                tick();if(front.empty())return true;
                if(metrics.states>=limits.maximumStates){metrics.stateBound=true;throw Bound{};}++metrics.states;
                if(selected.size()>=maximumCells){metrics.bounded=true;backtrack();continue;}
                bool first=true;unsigned examined=0;
                // Bounded most-constrained-face lookahead; unchanged geometry
                // and face-ledger admission guard every candidate.
                for(const auto& [k,f]:front){(void)k;std::vector<std::size_t> eligible;
                    for(auto index:options(f))if(ledgerAccepts(candidates[index].cell))eligible.push_back(index);
                    if(first||eligible.size()<frame.choices.size()){frame.choices=std::move(eligible);first=false;}
                    if(frame.choices.empty()||++examined==8)break;
                }
                const auto covers=[&](std::size_t index){unsigned n=0;for(unsigned side=0;side<4;++side)n+=front.contains(PlcGiftWrapping::key(tet_facet_detail::outwardFace(candidates[index].cell,side)));return n;};
                std::stable_sort(frame.choices.begin(),frame.choices.end(),[&](auto a,auto b){return covers(a)>covers(b);});frame.entered=true;
            }
            if(frame.next==frame.choices.size()){backtrack();continue;}
            const auto index=frame.choices[frame.next++];tick();if(!pairsAccept(index))continue;
            const auto changes=apply(candidates[index].cell);selected.push_back(index);
            Frame child;child.changes=changes;child.parent=true;stack.push_back(std::move(child));
        }
        return false;
    }
public:
    PlcPatchSearch(std::vector<TetFacetVertex> vertices,std::vector<Face> faces,
        const TetPlcPatchSearchOptions& options,unsigned maxCells,PlcPatchSearchMetrics& counters,
        std::function<void(std::uint64_t)> charge,std::function<bool(Vec3)> admission,std::function<bool(Vec3)> cellAdmission={})
        :limits(options),maximumCells(std::min(2048U,maxCells)),metrics(counters){
        geometry.vertices=std::move(vertices);geometry.boundary=std::move(faces);
        geometry.chargeWork=std::move(charge);geometry.interiorAdmission=std::move(admission);geometry.centroidAdmission=std::move(cellAdmission);
    }
    bool run(std::vector<TetFacetCell>& result){
        try{
            tick(geometry.vertices.size()+geometry.boundary.size());if(!geometry.initialize())return false;
            {geometry.fixedVertexCount=0;for(auto face:geometry.boundary)for(auto id:face)geometry.fixedVertexCount=std::max(geometry.fixedVertexCount,TetFacetVertexId(id+1));}
            for(auto f:geometry.boundary){const auto k=PlcGiftWrapping::key(f);boundary.emplace(k,sense(f));front.emplace(k,f);}
            if(!search())return false;
            for(auto index:selected)result.push_back(candidates[index].cell);
            return true;
        }catch(const Bound&){metrics.bounded=true;return false;}
    }
};
} // namespace paacs::bnd::detail
