// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private actual-use authority for closed native body boundaries.
#include "paacs/geometry/OcctModel.hpp"
#include "paacs/mesh/Mesh.hpp"
#include <algorithm>
#include <iterator>
#include <set>
#include <utility>
namespace paacs::geo::detail {
// Positive: unique actual wire. Zero: an internal face diagonal/no wire.
// Negative: ambiguous authority. A single face supplies a true wire only
// through repeated uses (a periodic seam), never through endpoint membership
// alone. charge accounts for every native use examined by this query.
template<class FaceTags,class Charge>
int actualOcctWire(const OcctModel& model,mesh::Classification a,
                  mesh::Classification b,const FaceTags& faces,Charge&& charge) {
    if(faces.empty()||a.dim>1||b.dim>1)return 0;
    const auto on=[](mesh::Classification c,const OcctEdgeUse& use){
        return c.dim==1?c.tag==use.edgeTag:
            c.dim==0&&(c.tag==use.startVertexTag||c.tag==use.endVertexTag);
    };
    std::set<int> candidates;bool first=true;
    for(const int face:faces){
        charge();std::set<int> available;
        for(const auto& wire:model.boundaryUses(face))for(const auto& use:wire.edges){
            charge();if(!use.degenerated&&on(a,use)&&on(b,use))available.insert(use.edgeTag);
        }
        if(first){candidates=std::move(available);first=false;}
        else{std::set<int> common;std::set_intersection(candidates.begin(),candidates.end(),
            available.begin(),available.end(),std::inserter(common,common.end()));candidates=std::move(common);}
    }
    if(candidates.empty())return 0;
    if(candidates.size()!=1)return -1;
    const int tag=*candidates.begin();
    if(faces.size()==1){
        unsigned uses=0;
        for(const auto& wire:model.boundaryUses(*faces.begin()))for(const auto& use:wire.edges){
            charge();if(!use.degenerated&&use.edgeTag==tag)++uses;
        }
        if(uses<2)return 0;
    }
    return tag;
}
} // namespace paacs::geo::detail
