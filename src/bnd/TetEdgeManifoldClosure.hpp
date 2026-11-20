// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Pure cavity preparation. A closed boundary here is only a proposal: the
// exact ball, replacement, protected-entity and transaction audits still run.
#include "paacs/bnd/TetRecoveryIndex.hpp"
#include <algorithm>
#include <map>
#include <queue>
#include <set>

namespace paacs::bnd::detail {
enum class ManifoldClosureStatus { Closed, Unavailable, BudgetReached };
struct ManifoldClosureOutcome {
    ManifoldClosureStatus status=ManifoldClosureStatus::Unavailable;
    std::set<topo::CellHandle> cells;
    std::uint64_t work=0,passes=0;
};
inline ManifoldClosureOutcome prepareTetEdgeManifoldClosure(
    const std::set<topo::CellHandle>& selected,const TetRecoveryContext& context,
    mesh::Classification region,std::uint64_t maxCells,std::uint64_t maxVertices,
    std::uint64_t maxWork,unsigned maxPasses=32) {
    ManifoldClosureOutcome out;
    auto proposal=selected;
    const auto visit=[&]() {if(out.work==maxWork)return false;++out.work;return true;};
    const auto budget=[&]() {out.status=ManifoldClosureStatus::BudgetReached;return out;};
    constexpr std::array<std::array<unsigned,3>,4> faces{{{1,2,3},{0,3,2},{0,1,3},{0,2,1}}};
    std::set<TetFacetFaceKey> protectedFaces;
    for(const auto& face:context.protectedFaces){if(!visit())return budget();protectedFaces.insert(face.key);}
    for(unsigned pass=0;;++pass) {
        if(proposal.size()>maxCells)return budget();
        std::set<DenseTetNodeId> vertices;
        std::map<DenseTetEdgeKey,unsigned> edgeUses;
        std::map<DenseTetNodeId,std::vector<std::array<DenseTetNodeId,2>>> links;
        for(const auto cell:proposal) {
            if(!visit())return budget();
            const auto& row=context.cells.at(cell);
            vertices.insert(row.nodes.begin(),row.nodes.end());
            if(vertices.size()>maxVertices)return budget();
            for(const auto& f:faces) {
                if(!visit())return budget();
                const auto key=DenseTetFaceKey::of({row.nodes[f[0]],row.nodes[f[1]],row.nodes[f[2]]});
                const auto found=context.faces.find(key);if(found==context.faces.end())return out;
                unsigned owners=0;
                for(const auto& owner:found->second.owners){if(!visit())return budget();owners+=proposal.contains(owner.cell)?1U:0U;}
                if(owners!=1)continue;
                for(unsigned j=0;j<3;++j) {
                    ++edgeUses[DenseTetEdgeKey::of(key.nodes[j],key.nodes[(j+1)%3])];
                    links[key.nodes[j]].push_back({key.nodes[(j+1)%3],key.nodes[(j+2)%3]});
                }
            }
        }
        std::vector<DenseTetNodeId> entity;
        for(const auto& [edge,count]:edgeUses){if(!visit())return budget();if(count!=2){entity.assign(edge.nodes.begin(),edge.nodes.end());break;}}
        if(entity.empty())for(const auto& [vertex,rows]:links) {
            std::map<DenseTetNodeId,std::vector<DenseTetNodeId>> graph;
            for(const auto& row:rows){if(!visit())return budget();graph[row[0]].push_back(row[1]);graph[row[1]].push_back(row[0]);}
            bool bad=false;
            for(auto& [node,next]:graph){(void)node;if(!visit())return budget();std::sort(next.begin(),next.end());bad|=next.size()!=2||std::adjacent_find(next.begin(),next.end())!=next.end();}
            std::set<DenseTetNodeId> seen;
            std::vector<DenseTetNodeId> todo{graph.begin()->first};seen.insert(todo[0]);
            while(!todo.empty()) {
                const auto at=todo.back();todo.pop_back();
                for(const auto next:graph.at(at)){if(!visit())return budget();if(seen.insert(next).second)todo.push_back(next);}
            }
            if(bad||seen.size()!=graph.size()){entity.push_back(vertex);break;}
        }
        if(entity.empty()){out.status=ManifoldClosureStatus::Closed;out.cells=std::move(proposal);return out;}
        if(pass==maxPasses)return budget();
        std::set<topo::CellHandle> missing;
        const auto insertMissing=[&](const auto& star) {
            for(const auto cell:star){if(!visit())return false;if(!proposal.contains(cell))missing.insert(cell);}return true;
        };
        if(entity.size()==2) {
            const auto star=context.edges.find(DenseTetEdgeKey::of(entity[0],entity[1]));
            if(star==context.edges.end())return out;
            if(!insertMissing(star->second.completeStar))return budget();
        }else {
            if(entity[0]>=context.nodeStars.size())return out;
            if(!insertMissing(context.nodeStars[entity[0]]))return budget();
        }
        std::vector<std::vector<topo::CellHandle>> gaps;
        while(!missing.empty()) {
            std::vector<topo::CellHandle> gap;
            std::queue<topo::CellHandle> pending;pending.push(*missing.begin());missing.erase(missing.begin());
            bool admissible=true;unsigned contacts=0;
            while(!pending.empty()) {
                const auto cell=pending.front();pending.pop();gap.push_back(cell);
                const auto& row=context.cells.at(cell);
                admissible&=row.classification==region;
                for(const auto& f:faces) {
                    if(!visit())return budget();
                    const auto key=DenseTetFaceKey::of({row.nodes[f[0]],row.nodes[f[1]],row.nodes[f[2]]});
                    const auto found=context.faces.find(key);if(found==context.faces.end())return out;
                    bool touches=false;
                    for(const auto& owner:found->second.owners){if(!visit())return budget();touches|=proposal.contains(owner.cell);}
                    if(touches){++contacts;if(protectedFaces.contains(canonicalTetFacetFace(key.nodes)))admissible=false;}
                    bool contains=true;
                    for(const auto node:entity)contains&=std::find(key.nodes.begin(),key.nodes.end(),node)!=key.nodes.end();
                    if(contains)for(const auto& owner:found->second.owners){if(!visit())return budget();if(missing.erase(owner.cell))pending.push(owner.cell);}
                }
            }
            if(admissible&&contacts>=2){
                std::sort(gap.begin(),gap.end(),[&](auto a,auto b){return context.cells.at(a).sequence<context.cells.at(b).sequence;});
                gaps.push_back(std::move(gap));
            }
        }
        if(gaps.empty())return out;
        std::sort(gaps.begin(),gaps.end(),[&](const auto& a,const auto& b){return a.size()!=b.size()?a.size()<b.size():context.cells.at(a.front()).sequence<context.cells.at(b.front()).sequence;});
        proposal.insert(gaps.front().begin(),gaps.front().end());++out.passes;
    }
}
} // namespace paacs::bnd::detail
