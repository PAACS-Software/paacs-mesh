// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private bounded repair of opposite coincident triangles at a curved trim.
// Native curve segments and all unselected (in particular shared cut) faces
// stay fixed. Callers still certify embedding and original-CAD fidelity.
#include "paacs/mesh/Mesh.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <span>
#include <vector>

namespace paacs::geo::detail {
using OcctAliasTriangle=std::array<int,3>;
using OcctAliasEdge=std::array<int,2>;
template<class Position,class Charge,class Admission>
unsigned repairOcctSurfaceAliases(std::vector<OcctAliasTriangle>& triangles,
    std::span<const int> faceTags,const std::set<OcctAliasEdge>& protectedEdges,
    std::span<const int> allowedFaces,unsigned maximumFlips,
    Position position,Charge charge,Admission admission) {
    if(!maximumFlips||allowedFaces.empty())return 0;
    std::set<int> allowed;for(int tag:allowedFaces){charge(1);allowed.insert(tag);}
    unsigned flips=0;
    while(flips<maximumFlips){
        std::map<OcctAliasTriangle,std::vector<std::size_t>> duplicates;
        std::map<OcctAliasEdge,std::vector<std::size_t>> edgeUses;
        for(std::size_t i=0;i<triangles.size();++i){charge(1);auto k=triangles[i];std::sort(k.begin(),k.end());duplicates[k].push_back(i);
            for(unsigned j=0;j<3;++j){const auto a=k[j],b=k[(j+1)%3];edgeUses[{std::min(a,b),std::max(a,b)}].push_back(i);}}
        bool changed=false;
        for(const auto& [key,pair]:duplicates){(void)key;charge(1);
            if(pair.size()!=2||faceTags[pair[0]]==faceTags[pair[1]])continue;
            const auto a=triangles[pair[0]],b=triangles[pair[1]];bool opposite=false;
            for(unsigned j=0;j<3;++j)opposite|=a[0]==b[j]&&a[1]==b[(j+2)%3]&&a[2]==b[(j+1)%3];
            if(!opposite)continue;
            for(auto i:pair){if(!allowed.contains(faceTags[i]))continue;const auto t=triangles[i];
                for(unsigned side=0;side<3&&!changed;++side){charge(1);
                    const int u=t[side],v=t[(side+1)%3],w=t[(side+2)%3];const OcctAliasEdge edge{std::min(u,v),std::max(u,v)};
                    if(protectedEdges.contains(edge))continue;
                    std::optional<std::size_t> other;unsigned owners=0;
                    for(auto j:edgeUses.at(edge)){charge(1);if(faceTags[j]==faceTags[i]){++owners;if(j!=i)other=j;}}
                    if(owners!=2||!other)continue;const auto j=*other;const auto q=triangles[j];int d=-1;
                    for(unsigned z=0;z<3;++z)if(q[z]==v&&q[(z+1)%3]==u)d=q[(z+2)%3];
                    if(d<0||d==w)continue;
                    const OcctAliasEdge diagonal{std::min(w,d),std::max(w,d)};
                    if(edgeUses.contains(diagonal)||protectedEdges.contains(diagonal))continue;
                    const OcctAliasTriangle first{w,u,d},second{w,d,v};
                    const auto area=[&](auto f){return (position(f[1])-position(f[0])).cross(position(f[2])-position(f[0])).norm();};
                    const auto firstArea=area(first),secondArea=area(second);
                    if(!(firstArea>0&&secondArea>0)||!std::isfinite(firstArea)||!std::isfinite(secondArea)||!admission(faceTags[i],first,second))continue;
                    triangles[i]=first;triangles[j]=second;++flips;changed=true;break;
                }if(changed)break;
            }if(changed)break;
        }if(!changed)break;
    }
    return flips;
}
} // namespace paacs::geo::detail
