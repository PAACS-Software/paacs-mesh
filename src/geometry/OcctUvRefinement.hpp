// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/Delaunay2D.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

namespace paacs::geo::detail {

struct OcctUvRefinementResult {
    unsigned added=0;
    bool completed=false;
    std::uint64_t samples=0,uniqueSamples=0,physicalQueries=0,physicalFailures=0;
    double remainingError=0;
};

template<class Position,class Evaluate,class Inside,class Append,class Charge,class Deviation=std::nullptr_t>
OcctUvRefinementResult refineOcctUvTriangles(std::vector<Vec2>& uv,std::vector<int>& global,
    std::vector<algo::Tri>& triangles,const std::set<std::array<int,2>>& protectedEdges,
    double maximumError,unsigned maximumAdded,Position position,Evaluate evaluate,Inside inside,
    Append append,Charge charge,Deviation deviation=nullptr) {
    using Edge=std::array<int,2>;
    struct Observation { Vec2 parameter{}; double error=0; };
    struct TriangleSamples {
        std::array<const Observation*,3> edges{};
        Observation center;
    };

    OcctUvRefinementResult result;
    std::map<std::pair<double,double>,std::optional<Vec3>> samples;
    std::map<std::array<double,5>,double> physicalErrors;
    std::map<Edge,Observation> edges;
    const auto key=[](Vec2 p){return std::pair{p.u,p.v};};
    const auto edgeKey=[](int a,int b){return Edge{std::min(a,b),std::max(a,b)};};
    const auto observe=[&](Vec2 p,Vec3 linear) {
        charge(1);
        auto [entry,added]=samples.try_emplace(key(p));
        if(added) {
            ++result.uniqueSamples;
            if(inside(p))entry->second=evaluate(p);
        }
        if(!entry->second)return Observation{p,0};
        auto error=(*entry->second-linear).norm();
        if constexpr(!std::is_same_v<Deviation,std::nullptr_t>) {
            if(error>maximumError) {
                charge(1);
                const std::array<double,5> query{p.u,p.v,linear.x,linear.y,linear.z};
                auto [cached,newQuery]=physicalErrors.try_emplace(query);
                if(newQuery) {
                    ++result.physicalQueries;
                    cached->second=deviation(linear);
                    result.physicalFailures+=!std::isfinite(cached->second)||cached->second<0;
                }
                if(std::isfinite(cached->second)&&cached->second>=0)
                    error=std::min(error,cached->second);
            }
        }
        return Observation{p,error};
    };
    const auto observeTriangle=[&](const algo::Tri& triangle) {
        TriangleSamples observed;
        const auto& f=triangle.v;
        for(unsigned side=0;side<3;++side) {
            charge(1);
            const auto a=f[side],b=f[(side+1)%3];
            const auto edge=edgeKey(int(a),int(b));
            if(protectedEdges.contains(edge))continue;
            auto [entry,added]=edges.try_emplace(edge);
            if(added)entry->second=observe((uv[a]+uv[b])*.5,(position(global[a])+position(global[b]))*.5);
            observed.edges[side]=&entry->second;
        }
        const Vec2 center{(uv[f[0]].u+uv[f[1]].u+uv[f[2]].u)/3.,
                          (uv[f[0]].v+uv[f[1]].v+uv[f[2]].v)/3.};
        observed.center=observe(center,(position(global[f[0]])+position(global[f[1]])+position(global[f[2]]))/3.);
        return observed;
    };

    std::set<std::pair<double,double>> points;
    for(const auto p:uv) {
        charge(1);
        points.insert(key(p));
    }
    std::vector<TriangleSamples> observations;
    observations.reserve(triangles.size());
    for(const auto& triangle:triangles)observations.push_back(observeTriangle(triangle));
    std::vector<algo::Tri> next;
    std::vector<TriangleSamples> nextObservations;
    for(;;) {
        std::optional<std::size_t> selected;
        std::optional<Edge> splitEdge;
        Vec2 proposed;
        double selectedError=maximumError;
        bool unresolved=false;
        const auto select=[&](std::size_t i,const Observation& observed,std::optional<Edge> edge) {
            charge(1);
            if(points.contains(key(observed.parameter))) {
                unresolved=true;
                return;
            }
            selected=i;
            splitEdge=edge;
            proposed=observed.parameter;
            selectedError=observed.error;
        };
        for(std::size_t i=0;i<triangles.size();++i) {
            charge(1);
            const auto& f=triangles[i].v;
            const auto& observed=observations[i];
            for(unsigned side=0;side<3;++side) {
                const auto* sample=observed.edges[side];
                if(!sample)continue;
                ++result.samples;
                if(sample->error>maximumError&&(!splitEdge||sample->error>selectedError))
                    select(i,*sample,edgeKey(int(f[side]),int(f[(side+1)%3])));
            }
            ++result.samples;
            // Resolve edges before face centers and retain the earliest witness
            // in chart order when several samples report the same error.
            if(!splitEdge&&observed.center.error>selectedError)select(i,observed.center,{});
        }
        result.remainingError=selected?selectedError:0;
        if(!selected) {
            result.completed=!unresolved;
            return result;
        }
        if(result.added==maximumAdded||uv.size()>=4096)return result;

        next.clear();
        nextObservations.clear();
        next.reserve(triangles.size()+2);
        nextObservations.reserve(triangles.size()+2);
        const auto index=Index(uv.size());
        std::vector<std::size_t> newTriangles;
        const auto appendTriangle=[&](algo::Tri triangle) {
            newTriangles.push_back(next.size());
            next.push_back(triangle);
            nextObservations.emplace_back();
        };
        unsigned edgeUses=0;
        for(std::size_t i=0;i<triangles.size();++i) {
            charge(1);
            const auto& f=triangles[i].v;
            if(splitEdge) {
                const auto previousUses=edgeUses;
                for(unsigned side=0;side<3;++side) {
                    const auto a=f[side],b=f[(side+1)%3],c=f[(side+2)%3];
                    if(edgeKey(int(a),int(b))!=*splitEdge)continue;
                    appendTriangle({{a,index,c}});
                    appendTriangle({{index,b,c}});
                    ++edgeUses;
                    break;
                }
                if(edgeUses!=previousUses)continue;
            } else if(i==*selected) {
                for(unsigned side=0;side<3;++side)appendTriangle({{f[side],f[(side+1)%3],index}});
                continue;
            }
            next.push_back(triangles[i]);
            nextObservations.push_back(observations[i]);
        }
        if(splitEdge&&edgeUses!=2)return result;
        const int id=append(proposed,*samples.at(key(proposed)));
        uv.push_back(proposed);
        global.push_back(id);
        points.insert(key(proposed));
        for(const auto i:newTriangles)nextObservations[i]=observeTriangle(next[i]);
        triangles.swap(next);
        observations.swap(nextObservations);
        ++result.added;
    }
}

} // namespace paacs::geo::detail
