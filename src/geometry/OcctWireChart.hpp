// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "OcctCurveMetric.hpp"
#include "OcctModelInternal.hpp"
#include "paacs/algo/Delaunay2DScalable.hpp"
#include "paacs/algo/PlanarConstraintRecovery.hpp"
#include "paacs/geometry/OcctModel.hpp"
#include "paacs/predicates/Predicates.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <span>
#include <vector>

namespace paacs::geo::detail {

struct OcctChartTriangle {
    std::array<int,3> nodes;
    std::array<Vec2,3> uv;
};

// Build a chart from shared CAD edge samples. Each use supplies its own UV
// corners, so a seam or pole can have multiple UV nodes for one 3D identity.
template<class Chains, class Position, class Tick, class Remaining, class Refine>
bool rebuildOcctWireChart(const OcctModel& model, int faceTag, const Chains& chains,
    std::span<const std::pair<int,Vec2>> interior, Position position, Tick tick,
    Remaining remaining, std::vector<OcctChartTriangle>& output,
    std::set<int>* crossingEdges, bool physicalMetric, Refine refine, bool fitPcurves=false) {
    struct Segment { Index a, b; bool degenerate; double tolerance; int tag; };
    constexpr std::size_t maximumChartNodes=4096;
    const auto* face=model.face(faceTag);
    std::vector<int> ids;
    std::vector<Vec2> uv;
    std::vector<std::vector<Index>> wires;
    std::vector<Segment> segments;
    std::set<int> poleIds;
    for(const auto& wire:model.boundaryUses(faceTag)){
        std::vector<Index> corners;
        std::vector<int> ends;
        const auto firstSegment=segments.size();
        for(const auto& use:wire.edges){
            tick(1);
            const auto& chain=chains.at(use.edgeTag);
            if(chain.ids.size()<2||use.sense==0)return false;
            const auto* edge=model.edge(use.edgeTag);
            double first,last;
            edge->parameterBounds(first,last);
            const auto startTolerance=model.tolerance(0,edge->startVertexTag());
            const auto endTolerance=model.tolerance(0,edge->endVertexTag());
            const auto tolerance=model.tolerance(1,use.edgeTag)+model.tolerance(2,faceTag)+
                std::max(startTolerance,endTolerance)+1e-7;
            if(use.degenerated)poleIds.insert(chain.ids.front());
            for(std::size_t step=0;step+1<chain.ids.size();++step){
                tick(1);
                const auto i=use.sense>0?step:chain.ids.size()-1-step;
                const auto next=use.sense>0?i+1:i-1;
                const auto t=i==0?first:(i+1==chain.ids.size()?last:chain.parameters[i]);
                Vec2 parameter;
                if(ids.size()>=maximumChartNodes)return false;
                if(!model.parameterOnUse(use.useId,t,parameter)){
                    if(!fitPcurves||use.seam||use.degenerated||i==0||i+1==chain.ids.size())return false;
                    // The fit evaluates 96 ternary-search samples, one residual
                    // and one final UV. Charge it before touching the CAD.
                    tick(98);
                    if(!OcctModelAccess::fitParameterOnUse(model,use.useId,chain.parameters[i-1],
                                                           t,chain.parameters[i+1],parameter))return false;
                }
                double allowed=tolerance;
                if(i==0)allowed+=startTolerance;
                else if(i+1==chain.ids.size())allowed+=endTolerance;
                const double error=(face->pointAt(parameter)-position(chain.ids[i])).norm();
                if(!std::isfinite(error)||error>allowed*1.01+1e-9)return false;
                const Index corner=static_cast<Index>(ids.size());
                ids.push_back(chain.ids[i]);
                uv.push_back(parameter);
                corners.push_back(corner);
                ends.push_back(chain.ids[next]);
                segments.push_back({corner,-1,use.degenerated,allowed,use.edgeTag});
            }
        }
        if(corners.size()<3)return false;
        for(std::size_t i=0;i<corners.size();++i){
            const auto next=corners[(i+1)%corners.size()];
            if(ids[next]!=ends[i])return false;
            segments[firstSegment+i].b=next;
        }
        wires.push_back(std::move(corners));
    }
    if(wires.empty())return false;

    // Collapse a tiny UV detour only along a native degenerate edge whose
    // physical samples all stay within the tolerance of the same CAD vertex.
    Vec2 lower,upper;
    face->paramBounds(lower,upper);
    const auto du=upper.u-lower.u,dv=upper.v-lower.v;
    std::vector<Index> parent(ids.size());
    std::iota(parent.begin(),parent.end(),Index{0});
    const auto root=[&](Index i){while(parent[i]!=i)i=parent[i];return i;};
    for(const auto& segment:segments){
        tick(1);
        if(!segment.degenerate||ids[segment.a]!=ids[segment.b])continue;
        const auto a=uv[segment.a],b=uv[segment.b];
        if(!(du>0&&dv>0&&std::hypot((a.u-b.u)/du,(a.v-b.v)/dv)<=.001))continue;
        bool within=true;
        for(unsigned i=0;i<=16&&within;++i){
            tick(1);
            const auto t=double(i)/16;
            const Vec2 p{a.u+(b.u-a.u)*t,a.v+(b.v-a.v)*t};
            within=(face->pointAt(p)-position(ids[segment.a])).norm()<=segment.tolerance;
        }
        if(within){
            const auto first=root(segment.a),second=root(segment.b);
            parent[std::max(first,second)]=std::min(first,second);
        }
    }
    std::vector<Index> remap(ids.size());
    std::size_t kept=0;
    for(Index i=0;i<static_cast<Index>(ids.size());++i){
        if(root(i)==i){remap[i]=static_cast<Index>(kept);ids[kept]=ids[i];uv[kept++]=uv[i];}
        else remap[i]=remap[root(i)];
    }
    ids.resize(kept);
    uv.resize(kept);
    for(auto& wire:wires){
        for(auto& id:wire)id=remap[id];
        wire.erase(std::unique(wire.begin(),wire.end()),wire.end());
        if(wire.size()>1&&wire.front()==wire.back())wire.pop_back();
        if(wire.size()<3)return false;
    }
    std::vector<algo::PlanarSegment> constraints;
    std::vector<int> constraintTags;
    for(const auto& s:segments){
        const auto a=remap[s.a],b=remap[s.b];
        if(a!=b){
            constraints.push_back({a,b});
            if(crossingEdges)constraintTags.push_back(s.degenerate?-1:s.tag);
        }
        else if(!s.degenerate)return false;
    }
    if(crossingEdges){
        const auto proper=[](Vec2 a,Vec2 b,Vec2 c,Vec2 d){
            const auto s=predicates::orient2d(a,b,c), t=predicates::orient2d(a,b,d);
            const auto x=predicates::orient2d(c,d,a), y=predicates::orient2d(c,d,b);
            return s!=predicates::Sign::Zero&&t!=predicates::Sign::Zero&&s!=t&&
                   x!=predicates::Sign::Zero&&y!=predicates::Sign::Zero&&x!=y;
        };
        for(std::size_t i=0;i<constraints.size();++i){
            const auto [a,b]=constraints[i];
            for(std::size_t j=i+1;j<constraints.size();++j){
                tick(1);
                const auto [c,d]=constraints[j];
                if(a==c||a==d||b==c||b==d)continue;
                if(std::max(uv[a].u,uv[b].u)<std::min(uv[c].u,uv[d].u)||
                   std::max(uv[c].u,uv[d].u)<std::min(uv[a].u,uv[b].u)||
                   std::max(uv[a].v,uv[b].v)<std::min(uv[c].v,uv[d].v)||
                   std::max(uv[c].v,uv[d].v)<std::min(uv[a].v,uv[b].v))continue;
                if(!proper(uv[a],uv[b],uv[c],uv[d]))continue;
                if(constraintTags[i]>0)crossingEdges->insert(constraintTags[i]);
                if(constraintTags[j]>0)crossingEdges->insert(constraintTags[j]);
            }
        }
        return true;
    }
    const auto inDomain=[&](Vec2 p){
        bool inside=false;
        for(const auto& wire:wires)for(std::size_t i=0,j=wire.size()-1;i<wire.size();j=i++){
            tick(1);
            const auto a=uv[wire[j]],b=uv[wire[i]];
            if((a.v>p.v)==(b.v>p.v))continue;
            const auto sign=predicates::orient2d(a,b,p);
            if((b.v>a.v&&sign==predicates::Sign::Positive)||
               (b.v<a.v&&sign==predicates::Sign::Negative))inside=!inside;
        }
        return inside;
    };
    for(const auto& [id,parameter]:interior){
        if(!inDomain(parameter))continue;
        if(ids.size()>=maximumChartNodes)return false;
        ids.push_back(id);
        uv.push_back(parameter);
    }
    tick(std::uint64_t(uv.size())*uv.size());
    std::vector<Vec2> metricCoordinates;
    if(physicalMetric&&!metricOcctChart(model,faceTag,uv,tick,metricCoordinates))return false;
    const auto& coordinates=physicalMetric?metricCoordinates:uv;
    auto triangles=algo::delaunay2DScalable(coordinates);
    const auto recovery=algo::recoverPlanarConstraints(coordinates,constraints,triangles,
        std::min<std::uint64_t>(5'000'000,remaining()));
    tick(recovery.work);
    if(!recovery.completed)return false;
    if(physicalMetric){
        std::set<std::array<Index,2>> barriers;
        std::set<std::array<int,2>> protectedEdges;
        for(const auto& [a,b]:constraints){
            tick(1);
            const auto [lo,hi]=std::minmax(a,b);
            barriers.insert({lo,hi});
            protectedEdges.insert({int(lo),int(hi)});
        }
        const auto parity=labelOcctChartRegions(triangles,barriers,tick);
        if(!parity)return false;
        std::size_t retained=0;
        for(std::size_t i=0;i<triangles.size();++i){
            tick(1);
            if((*parity)[i])triangles[retained++]=triangles[i];
        }
        triangles.resize(retained);
        if(!refine(uv,ids,triangles,protectedEdges))return false;
    }
    std::map<std::array<Index,2>,unsigned> incidence;
    for(auto triangle:triangles){
        tick(1);
        const auto v=triangle.v;
        const auto a=uv[v[0]],b=uv[v[1]],c=uv[v[2]];
        if(!physicalMetric&&!inDomain({(a.u+b.u+c.u)/3.,(a.v+b.v+c.v)/3.}))continue;
        if(predicates::orient2d(a,b,c)!=predicates::Sign::Positive)return false;
        for(unsigned i=0;i<3;++i){
            const auto [lo,hi]=std::minmax(v[i],v[(i+1)%3]);
            ++incidence[{lo,hi}];
        }
        bool collapsed=false;
        for(unsigned i=0;i<3;++i)if(ids[v[i]]==ids[v[(i+1)%3]]){
            if(!poleIds.contains(ids[v[i]]))return false;
            collapsed=true;
        }
        if(!collapsed)output.push_back({{ids[v[0]],ids[v[1]],ids[v[2]]},{a,b,c}});
    }
    for(const auto& segment:constraints){
        const auto [lo,hi]=std::minmax(segment.first,segment.second);
        if(incidence[{lo,hi}]!=1)return false;
    }
    return !output.empty();
}

} // namespace paacs::geo::detail
