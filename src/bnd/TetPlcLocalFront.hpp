// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "TetPlcGiftWrapping.hpp"

namespace paacs::bnd::detail {

// Retain the occupied volume and propose one Tet4 from incident open faces.
// Three faces collapse to one, or four close a pocket. A two-face advance
// must expose a second, reducing move and both are checked before insertion.
// Gift wrapping supplies the exact orientation, contact and native checks.
template<class Charge, class Accepted>
void advancePlcLocalFront(PlcGiftWrapping& mesh, unsigned maximumTets, Charge tick, Accepted accept) {
    using Face=PlcGiftWrapping::Face;
    using Edge=std::array<TetFacetVertexId,2>;
    using Tet=std::array<TetFacetVertexId,4>;
    constexpr unsigned maximumFront=4096;
    constexpr unsigned maximumEdgeUses=16;
    unsigned inserted=0;
    while(inserted<maximumTets&&!mesh.front.empty()&&mesh.front.size()<=maximumFront){
        std::map<Edge,std::vector<Face>> edges;
        for(const auto& [face,row]:mesh.front){
            (void)row;
            tick(3);
            for(unsigned i=0;i<3;++i){
                const auto [lo,hi]=std::minmax(face[i],face[(i+1)%3]);
                edges[{lo,hi}].push_back(face);
            }
        }
        std::set<Tet> unique;
        std::array<std::vector<Tet>,3> candidates;
        const auto third=[](Face f,Edge edge){
            for(auto v:f)if(v!=edge[0]&&v!=edge[1])return v;
            return kInvalidTetFacetVertex;
        };
        const auto matching=[&](const Tet& ids){
            tick(4);
            const TetFacetCell cell{ids,0,{3,0}};
            unsigned faces=0;
            for(unsigned side=0;side<4;++side)
                faces+=mesh.front.contains(PlcGiftWrapping::key(tet_facet_detail::outwardFace(cell,side)));
            return faces;
        };
        for(const auto& [edge,faces]:edges){
            tick(1);
            if(faces.size()>maximumEdgeUses)continue;
            for(std::size_t i=0;i<faces.size();++i)for(std::size_t j=0;j<i;++j){
                tick(1);
                Tet ids{edge[0],edge[1],third(faces[i],edge),third(faces[j],edge)};
                std::sort(ids.begin(),ids.end());
                const auto count=matching(ids);
                if(count>=2&&unique.insert(ids).second){
                    if(unique.size()>32768)return;
                    candidates[4-count].push_back(ids);
                }
            }
        }
        const auto before=inserted;
        const auto fromOpenFace=[&](const Tet& ids,TetFacetCell& cell){
            const TetFacetCell candidate{ids,0,{3,0}};
            for(unsigned side=0;side<4;++side){
                const auto face=PlcGiftWrapping::key(tet_facet_detail::outwardFace(candidate,side));
                if(auto row=mesh.front.find(face);row!=mesh.front.end())
                    return mesh.valid(row->second,ids[side],cell);
            }
            return false;
        };
        for(unsigned batch=0;batch<2;++batch)for(const auto& ids:candidates[batch]){
            if(inserted==maximumTets)return;
            const auto count=matching(ids);
            if(count<3)continue;
            TetFacetCell cell;
            if(fromOpenFace(ids,cell)){
                mesh.insert(cell);
                ++inserted;
                accept(count==4);
            }
        }
        if(inserted!=before)continue;
        if(maximumTets-inserted<2)break;
        const auto tryPair=[&](const Tet& ids){
            TetFacetCell first;
            if(!fromOpenFace(ids,first))return false;
            for(unsigned side=0;side<4;++side){
                auto joint=tet_facet_detail::outwardFace(first,side);
                if(mesh.front.contains(PlcGiftWrapping::key(joint)))continue;
                std::swap(joint[1],joint[2]);
                for(unsigned i=0;i<3;++i){
                    const auto [lo,hi]=std::minmax(joint[i],joint[(i+1)%3]);
                    const Edge edge{lo,hi};
                    const auto adjacent=edges.find(edge);
                    if(adjacent==edges.end()||adjacent->second.size()>maximumEdgeUses)continue;
                    for(auto face:adjacent->second){
                        tick(1);
                        const auto apex=third(face,edge);
                        if(std::ranges::find(ids,apex)!=ids.end())continue;
                        // Distinct tets share just the proposed joint. The
                        // other three faces still use the current ledger.
                        const Tet secondIds{joint[0],joint[1],joint[2],apex};
                        const unsigned count=1+matching(secondIds);
                        if(count<3)continue;
                        TetFacetCell second;
                        if(!mesh.valid(PlcGiftWrapping::Front{joint},apex,second))continue;
                        const TetFacetComplexView view{mesh.vertices,mesh.cells,{},{},{},{3,0}};
                        if(mesh.cellsConflict(view,first,second))continue;
                        second.sequence=first.sequence+1;
                        mesh.insert(first);
                        mesh.insert(second);
                        inserted+=2;
                        accept(false);
                        accept(count==4);
                        return true;
                    }
                }
            }
            return false;
        };
        unsigned trials=0;
        for(const auto& ids:candidates[2]){
            tick(1);
            if(++trials>96||tryPair(ids))break;
        }
        if(inserted==before)break;
    }
}
}
