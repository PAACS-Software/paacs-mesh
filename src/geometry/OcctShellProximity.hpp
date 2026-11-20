// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "../bnd/PlcSpatialIndex.hpp"
#include "paacs/mesh/Mesh.hpp"
#include <array>
#include <cmath>
#include <span>
#include <vector>

namespace paacs::geo::detail {

// The distance filter selects CAD queries for repair proposals. Publication
// still checks every centroid of the final volume against the original CAD.
template<class Charge>
std::vector<std::size_t> occtCentroidsNearShell(const mesh::Mesh& shell,
    std::span<const Vec3> centroids,double distance,std::uint64_t maximumWork,Charge charge) {
    using SpatialIndex=bnd::detail::PlcGrowingAabbIndex;
    struct Limit {};
    const auto tick=[&](std::uint64_t work){
        if(work>maximumWork)throw Limit{};
        maximumWork-=work;
        charge(work);
    };
    const auto distanceSquared=[](Vec3 p,Vec3 a,Vec3 b,Vec3 c){
        const auto u=b-a,v=c-a,w=p-a;
        const auto normal=u.cross(v);
        const auto square=normal.norm2();
        if(square>0.){
            const auto beta=w.cross(v).dot(normal)/square;
            const auto gamma=u.cross(w).dot(normal)/square;
            if(beta>=0.&&gamma>=0.&&beta+gamma<=1.){
                const auto height=w.dot(normal);
                return height*height/square;
            }
        }
        const auto segment=[&](Vec3 from,Vec3 to){
            const auto edge=to-from;
            const auto length=edge.norm2();
            const auto fraction=length>0.?std::clamp((p-from).dot(edge)/length,0.,1.):0.;
            return (p-from-edge*fraction).norm2();
        };
        return std::min({segment(a,b),segment(b,c),segment(c,a)});
    };
    std::vector<SpatialIndex::Box> boxes;
    boxes.reserve(shell.cellCount());
    try {
        for(Index i=0;i<shell.cellCount();++i){
            tick(1);
            const auto* ids=shell.cellNodes(i);
            const auto a=shell.node(ids[0]),b=shell.node(ids[1]),c=shell.node(ids[2]);
            boxes.push_back({{std::min({a.x,b.x,c.x}),std::min({a.y,b.y,c.y}),std::min({a.z,b.z,c.z})},
                             {std::max({a.x,b.x,c.x}),std::max({a.y,b.y,c.y}),std::max({a.z,b.z,c.z})}});
        }
        const auto bounds=[&](std::size_t id){return boxes[id];};
        SpatialIndex index;
        index.rebuild(boxes.size(),bounds,tick);
        std::vector<std::size_t> near;
        const Vec3 margin{distance,distance,distance};
        for(std::size_t i=0;i<centroids.size();++i){
            tick(1);
            const auto p=centroids[i];
            const auto outside=index.all(boxes.size(),{p-margin,p+margin},bounds,tick,
                [&](std::size_t facet){
                    tick(1);
                    const auto* ids=shell.cellNodes(Index(facet));
                    const auto square=distanceSquared(p,shell.node(ids[0]),shell.node(ids[1]),shell.node(ids[2]));
                    return square>distance*distance;
                });
            if(!outside)near.push_back(i);
        }
        return near;
    }catch(const Limit&){
        return {};
    }
}

}
