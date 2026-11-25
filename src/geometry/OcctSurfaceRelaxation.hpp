// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/Delaunay2D.hpp"
#include "OcctChartTrim.hpp"
#include "OcctTrimMemo.hpp"
#include "paacs/geometry/OcctModel.hpp"
#include "paacs/predicates/Predicates.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace paacs::geo::detail {

inline double occtTriangleMeanRatio(const std::array<Vec3,3>& p) {
    const auto a=p[1]-p[0],b=p[2]-p[1],c=p[0]-p[2];
    const auto squared=a.norm2()+b.norm2()+c.norm2();
    return squared>0?2*std::sqrt(3.)*a.cross(b).norm()/squared:0.;
}

inline bool occtNearlyEquilateral(const std::array<Vec3,3>& p) {
    // Check angles of 45–75 degrees using their squared cosines.
    constexpr double minimumCosineSquared=.0669872981077807; // cos²(75°)
    for(unsigned i=0;i<3;++i){
        const auto a=p[(i+1)%3]-p[i],b=p[(i+2)%3]-p[i];
        const auto dot=a.dot(b),squared=a.norm2()*b.norm2();
        if(!(dot>0&&dot*dot>=minimumCosineSquared*squared&&dot*dot<=.5*squared))return false;
    }
    return true;
}

// Sweep interior diagonals once. Both triangles must keep their UV
// orientation, and the original CAD wires remain constraints.
template<class Position,class Charge>
void improveOcctSurfaceDiagonals(const GFace& face,const std::vector<Vec2>& uv,
    const std::vector<int>& ids,std::vector<algo::Tri>& triangles,
    const std::set<std::array<int,2>>& protectedEdges,double chord,Position position,Charge charge,
    OcctTrimMemo* stockTrim=nullptr) {
    using Edge=std::array<Index,2>;
    constexpr auto unused=std::numeric_limits<std::size_t>::max();
    const auto edge=[](Index a,Index b){return Edge{std::min(a,b),std::max(a,b)};};
    std::map<Edge,std::array<std::size_t,2>> owners;
    const auto add=[&](Edge e,std::size_t i){
        charge(1);
        auto [found,inserted]=owners.try_emplace(e,std::array{i,unused});
        if(!inserted)found->second[1]=i;
    };
    const auto sides=[&](const algo::Tri& t){
        const auto [a,b,c]=t.v;
        return std::array{edge(a,b),edge(b,c),edge(c,a)};
    };
    for(std::size_t i=0;i<triangles.size();++i)
        for(const auto e:sides(triangles[i]))add(e,i);
    const auto points=[&](const algo::Tri& t){
        return std::array{position(ids[t.v[0]]),position(ids[t.v[1]]),position(ids[t.v[2]])};
    };
    const auto admitted=[&](const algo::Tri& t){
        const auto [a,b,c]=t.v;
        if(predicates::orient2d(uv[a],uv[b],uv[c])!=predicates::Sign::Positive)return false;
        charge(1);
        const auto parameter=(uv[a]+uv[b]+uv[c])*(1./3);
        const auto p=points(t);
        return (face.pointAt(parameter)-(p[0]+p[1]+p[2])/3.).norm()<=chord&&
            (stockTrim?stockTrim->contains(parameter):face.containsParam(parameter));
    };
    std::vector<Edge> candidates;
    for(const auto& [e,uses]:owners){
        charge(1);
        if(uses[1]!=unused&&!protectedEdges.contains({int(e[0]),int(e[1])}))candidates.push_back(e);
    }
    for(const auto& e:candidates){
        charge(1);
        const auto found=owners.find(e);
        if(found==owners.end()||found->second[1]==unused)continue;
        const auto [i,j]=found->second;
        const auto a=triangles[i],b=triangles[j];
        const auto pa=points(a),pb=points(b);
        const auto qa=occtTriangleMeanRatio(pa),qb=occtTriangleMeanRatio(pb);
        if(std::min(qa,qb)>=.86)continue;
        auto [u,v]=e;
        Index w=-1,d=-1;
        for(unsigned k=0;k<3;++k)
            if(edge(a.v[k],a.v[(k+1)%3])==e){u=a.v[k];v=a.v[(k+1)%3];w=a.v[(k+2)%3];break;}
        for(unsigned k=0;k<3;++k)
            if(b.v[k]==v&&b.v[(k+1)%3]==u){d=b.v[(k+2)%3];break;}
        if(w<0||d<0||w==d||owners.contains(edge(w,d)))continue;
        std::array vertexIds{ids[u],ids[v],ids[w],ids[d]};
        std::ranges::sort(vertexIds);
        if(std::ranges::adjacent_find(vertexIds)!=vertexIds.end())continue;
        const algo::Tri first{{w,u,d}},second{{w,d,v}};
        const auto p1=points(first),p2=points(second);
        const auto q1=occtTriangleMeanRatio(p1),q2=occtTriangleMeanRatio(p2);
        if(std::min(q1,q2)<=std::min(qa,qb)+.015||q1+q2<qa+qb)continue;
        if(occtNearlyEquilateral(p1)+occtNearlyEquilateral(p2)<
           occtNearlyEquilateral(pa)+occtNearlyEquilateral(pb))continue;
        const auto mid=(uv[w]+uv[d])*.5;
        charge(1);
        const auto midpointError=(face.pointAt(mid)-(position(ids[w])+position(ids[d]))*.5).norm();
        if(!(midpointError<=chord)||
            !(stockTrim?stockTrim->contains(mid):face.containsParam(mid))||
            !admitted(first)||!admitted(second))continue;
        for(const auto owner:{i,j})for(const auto key:sides(triangles[owner])) {
            charge(1);
            const auto at=owners.find(key);
            auto& uses=at->second;
            if(uses[0]==owner)uses[0]=uses[1];
            uses[1]=unused;
            if(uses[0]==unused)owners.erase(at);
        }
        triangles[i]=first;
        triangles[j]=second;
        for(const auto owner:{i,j})for(const auto key:sides(triangles[owner]))add(key,owner);
    }
}

// One sweep of the face-interior nodes. Every accepted move keeps the
// constrained chart and has a CAD distance witness at each affected sample.
template<class Position,class Update,class Charge>
void relaxOcctSurface(const OcctModel& model,int faceTag,std::vector<Vec2>& uv,
    const std::vector<int>& ids,const std::vector<algo::Tri>& triangles,
    const std::set<std::array<int,2>>& protectedEdges,double chord,
    Position position,Update update,Charge charge,bool physicalWeights=false,OcctChartTrim* trim=nullptr,
    OcctTrimMemo* stockTrim=nullptr) {
    const auto* face=model.face(faceTag);
    std::vector<bool> fixed(uv.size());
    for(const auto& [a,b]:protectedEdges){
        charge(1);
        fixed[a]=fixed[b]=true;
    }
    std::vector<std::vector<std::size_t>> incident(uv.size());
    for(std::size_t i=0;i<triangles.size();++i){
        charge(3);
        for(const auto v:triangles[i].v)if(!fixed[v])incident[v].push_back(i);
    }
    // A chart stays on the same original trimmed face throughout this sweep.
    // Cache trim decisions by parameter and test each current chord sample.
    constexpr std::size_t maximumTrimDecisions=4096;
    std::map<std::pair<double,double>,bool> trimDecisions;
    const auto contains=[&](Vec2 parameter){
        if(trim)return trim->contains(parameter);
        if(stockTrim&&stockTrim->enabled())return stockTrim->contains(parameter);
        if(!std::isfinite(parameter.u)||!std::isfinite(parameter.v)){
            if(stockTrim)stockTrim->nonfinite();
            return false;
        }
        const auto key=std::pair{parameter.u,parameter.v};
        if(const auto found=trimDecisions.find(key);found!=trimDecisions.end()){
            if(stockTrim)stockTrim->legacyHit();
            return found->second;
        }
        const auto accepted=stockTrim?stockTrim->contains(parameter):face->containsParam(parameter);
        if(trimDecisions.size()<maximumTrimDecisions)trimDecisions.try_emplace(key,accepted);
        return accepted;
    };
    const auto witness=[&](Vec2 parameter,Vec3 sample){
        charge(1);
        return (face->pointAt(parameter)-sample).norm()<=chord&&contains(parameter);
    };
    std::vector<Index> neighbours;
    for(Index v=0;v<static_cast<Index>(uv.size());++v){
        charge(1);
        if(incident[v].empty())continue;
        neighbours.clear();
        double before=1.;
        unsigned regularBefore=0;
        for(const auto i:incident[v]){
            charge(1);
            const auto& t=triangles[i].v;
            const std::array points{position(ids[t[0]]),position(ids[t[1]]),position(ids[t[2]])};
            before=std::min(before,occtTriangleMeanRatio(points));
            if(physicalWeights&&occtNearlyEquilateral(points))++regularBefore;
            for(const auto n:t)if(n!=v)neighbours.push_back(n);
        }
        if(physicalWeights&&before>=.92)continue;
        std::sort(neighbours.begin(),neighbours.end());
        neighbours.erase(std::unique(neighbours.begin(),neighbours.end()),neighbours.end());
        const auto original=position(ids[v]);
        Vec2 displacement{};
        double weightSum=0;
        for(const auto n:neighbours){
            const auto delta=uv[n]-uv[v];
            double weight=1;
            if(physicalWeights){
                charge(1);
                const auto paramLength=delta.norm();
                if(!(paramLength>0))continue;
                // The CAD chord gives the physical scale in this UV direction
                // using the mesh points already evaluated on the face.
                weight=(position(ids[n])-original).norm()/paramLength;
                if(!(std::isfinite(weight)&&weight>0))continue;
            }
            displacement=displacement+delta*weight;
            weightSum+=weight;
        }
        if(!(std::isfinite(weightSum)&&weightSum>0))continue;
        displacement=displacement*(1./weightSum);
        const std::array fractions{1.,.5,.25};
        for(std::size_t attempt=physicalWeights?0:1;attempt<fractions.size();++attempt){
            const double fraction=fractions[attempt];
            charge(1);
            const auto trial=uv[v]+displacement*fraction;
            if((trial.u==uv[v].u&&trial.v==uv[v].v)||
               !std::isfinite(trial.u)||!std::isfinite(trial.v))continue;
            const auto point=face->pointAt(trial);
            if(point.x==original.x&&point.y==original.y&&point.z==original.z)continue;
            const auto proposedUv=[&](Index n){return n==v?trial:uv[n];};
            const auto proposedPoint=[&](Index n){return n==v?point:position(ids[n]);};
            unsigned regularAfter=0;
            const auto improves=[&](std::size_t i){
                charge(1);
                const auto& t=triangles[i].v;
                const std::array points{proposedPoint(t[0]),proposedPoint(t[1]),proposedPoint(t[2])};
                if(physicalWeights&&occtNearlyEquilateral(points))++regularAfter;
                return predicates::orient2d(proposedUv(t[0]),proposedUv(t[1]),proposedUv(t[2]))==
                        predicates::Sign::Positive&&
                    occtTriangleMeanRatio(points)>before+1e-8;
            };
            if(!std::ranges::all_of(incident[v],improves)||regularAfter<regularBefore)continue;
            if(!contains(trial))continue;
            const auto fits=[&](std::size_t i){
                charge(1);
                const auto& t=triangles[i].v;
                const std::array u{proposedUv(t[0]),proposedUv(t[1]),proposedUv(t[2])};
                const std::array p{proposedPoint(t[0]),proposedPoint(t[1]),proposedPoint(t[2])};
                const auto center=(u[0]+u[1]+u[2])*(1./3);
                if(!witness(center,(p[0]+p[1]+p[2])/3.))return false;
                for(unsigned j=0;j<3;++j){
                    const auto k=(j+1)%3;
                    if(t[j]!=v&&t[k]!=v)continue;
                    const auto mid=(u[j]+u[k])*.5;
                    const auto sample=(p[j]+p[k])*.5;
                    if(!witness(mid,sample)&&!witness(mid+(center-mid)*.125,sample))return false;
                }
                return true;
            };
            if(!std::ranges::all_of(incident[v],fits))continue;
            update(ids[v],point);
            uv[v]=trial;
            break;
        }
    }
}

} // namespace paacs::geo::detail
