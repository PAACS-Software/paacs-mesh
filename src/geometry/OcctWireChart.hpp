// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "OcctBoundaryInternal.hpp"
#include "OcctCurveMetric.hpp"
#include "OcctModelInternal.hpp"
#include "paacs/algo/Delaunay2DScalable.hpp"
#include "paacs/algo/PlanarConstraintRecovery.hpp"
#include "paacs/geometry/OcctModel.hpp"
#include "paacs/predicates/Predicates.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <iterator>
#include <limits>
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
struct OcctChartSplitInterval {
    int edgeTag=-1;
    double first=0,last=0;
    friend auto operator<=>(const OcctChartSplitInterval&,const OcctChartSplitInterval&)=default;
};

// Refine only canonical intervals implicated by proper chart crossings. The
// source curves and existing node identities stay fixed. Grouping by edge
// makes one linear chain pass instead of one vector insertion per split.
template<class Chains,class Append,class Tick,class Admit>
bool bisectOcctChartIntervals(Chains& chains,
    const std::set<OcctChartSplitInterval>& intervals,Append append,Tick tick,Admit admit,
    OcctBoundaryConstructionWitness* witness=nullptr) {
    if(witness)*witness={};
    for(auto begin=intervals.begin();begin!=intervals.end();){
        const int tag=begin->edgeTag;
        auto end=begin;
        while(end!=intervals.end()&&end->edgeTag==tag){tick(1);++end;}
        const auto count=std::size_t(std::distance(begin,end));
        auto& chain=chains.at(tag);
        const auto refuse=[&](const char* detail,bool geometric){
            if(witness)*witness={"shared_chart_sampling",detail,tag,-1,std::nullopt,geometric};
            return false;
        };
        if(chain.ids.size()!=chain.parameters.size()||chain.ids.size()<2)
            return refuse("canonical_chain_invalid",false);
        if(count>(std::numeric_limits<std::size_t>::max)()-chain.ids.size())
            return refuse("canonical_chain_storage",false);
        admit(chain.ids.size()+count);
        std::vector<int> ids;
        std::vector<double> parameters;
        ids.reserve(chain.ids.size()+count);parameters.reserve(chain.parameters.size()+count);
        ids.push_back(chain.ids.front());parameters.push_back(chain.parameters.front());
        auto selected=begin;
        std::uint64_t added=0;
        for(std::size_t i=1;i<chain.ids.size();++i){
            tick(1);
            const double first=chain.parameters[i-1],last=chain.parameters[i];
            if(!(std::isfinite(first)&&std::isfinite(last)&&first<last))
                return refuse("canonical_parameter_order",false);
            if(selected!=end&&selected->first==first){
                if(selected->last!=last)return refuse("canonical_interval_missing",false);
                const double at=std::midpoint(first,last);
                if(!(first<at&&at<last))return refuse("chart_parameter_resolution",true);
                tick(1);
                ids.push_back(append(tag,at,std::uint64_t(chain.ids.size())+added));
                parameters.push_back(at);
                ++selected;++added;
            }
            ids.push_back(chain.ids[i]);parameters.push_back(last);
        }
        if(selected!=end)return refuse("canonical_interval_missing",false);
        chain.ids.swap(ids);chain.parameters.swap(parameters);
        begin=end;
    }
    return true;
}

// Build a chart from shared CAD edge samples. Each use supplies its own UV
// corners, so a seam or pole can have multiple UV nodes for one 3D identity.
template<class Chains, class Position, class Tick, class Remaining, class Refine>
bool rebuildOcctWireChart(const OcctModel& model, int faceTag, const Chains& chains,
    std::span<const std::pair<int,Vec2>> interior, Position position, Tick tick,
    Remaining remaining, std::vector<OcctChartTriangle>& output,
    std::set<int>* crossingEdges, bool physicalMetric, Refine refine, bool fitPcurves=false,
    std::size_t maximumChartNodes=4096,OcctBoundaryConstructionWitness* witness=nullptr,
    std::set<OcctChartSplitInterval>* splitIntervals=nullptr) {
    if(witness)*witness={};
    const auto refuse=[&](const char* detail,bool geometric=false,int edge=-1,
                         std::optional<std::uint64_t> use=std::nullopt){
        if(witness)*witness={"wire_chart",detail,edge,faceTag,use,geometric};
        return false;
    };
    struct Segment {
        Index a,b;bool degenerate;double tolerance;int tag;
        double first,last;
    };
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
            if(chain.ids.size()<2)return refuse("curve_samples_missing",false,use.edgeTag,use.useId);
            if(use.sense==0)return refuse("edge_use_orientation",false,use.edgeTag,use.useId);
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
                if(ids.size()>=maximumChartNodes)return refuse("chart_node_limit",false,use.edgeTag,use.useId);
                if(!model.parameterOnUse(use.useId,t,parameter)){
                    if(!fitPcurves||use.seam||use.degenerated||i==0||i+1==chain.ids.size())
                        return refuse("pcurve_unavailable",true,use.edgeTag,use.useId);
                    // The fit evaluates 96 ternary-search samples, one residual
                    // and one final UV. Charge it before touching the CAD.
                    tick(98);
                    if(!OcctModelAccess::fitParameterOnUse(model,use.useId,chain.parameters[i-1],
                                                           t,chain.parameters[i+1],parameter))
                        return refuse("pcurve_fit_failed",true,use.edgeTag,use.useId);
                }
                double allowed=tolerance;
                if(i==0)allowed+=startTolerance;
                else if(i+1==chain.ids.size())allowed+=endTolerance;
                const double error=(face->pointAt(parameter)-position(chain.ids[i])).norm();
                if(!std::isfinite(error))return refuse("pcurve_nonfinite_position",false,use.edgeTag,use.useId);
                if(error>allowed*1.01+1e-9)return refuse("pcurve_position_mismatch",true,use.edgeTag,use.useId);
                const Index corner=static_cast<Index>(ids.size());
                ids.push_back(chain.ids[i]);
                uv.push_back(parameter);
                corners.push_back(corner);
                ends.push_back(chain.ids[next]);
                const auto end=next==0?first:(next+1==chain.ids.size()?last:chain.parameters[next]);
                segments.push_back({corner,-1,use.degenerated,allowed,use.edgeTag,std::min(t,end),std::max(t,end)});
            }
        }
        if(corners.size()<3)return refuse("wire_too_short");
        for(std::size_t i=0;i<corners.size();++i){
            const auto next=corners[(i+1)%corners.size()];
            if(ids[next]!=ends[i])return refuse("wire_disconnected");
            segments[firstSegment+i].b=next;
        }
        wires.push_back(std::move(corners));
    }
    if(wires.empty())return refuse("wires_missing");

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
        if(wire.size()<3)return refuse("collapsed_wire");
    }
    std::vector<algo::PlanarSegment> constraints;
    std::vector<int> constraintTags;
    std::vector<OcctChartSplitInterval> constraintIntervals;
    for(const auto& s:segments){
        const auto a=remap[s.a],b=remap[s.b];
        if(a!=b){
            constraints.push_back({a,b});
            if(crossingEdges||splitIntervals){
                constraintTags.push_back(s.degenerate?-1:s.tag);
                if(splitIntervals)constraintIntervals.push_back({s.tag,s.first,s.last});
            }
        }
        else if(!s.degenerate)return refuse("nondegenerate_edge_collapsed",false,s.tag);
    }
    if(crossingEdges||splitIntervals){
        const auto proper=[](Vec2 a,Vec2 b,Vec2 c,Vec2 d){
            const auto s=predicates::orient2d(a,b,c), t=predicates::orient2d(a,b,d);
            const auto x=predicates::orient2d(c,d,a), y=predicates::orient2d(c,d,b);
            return s!=predicates::Sign::Zero&&t!=predicates::Sign::Zero&&s!=t&&
                   x!=predicates::Sign::Zero&&y!=predicates::Sign::Zero&&x!=y;
        };
        // Sweep the segment boxes before invoking exact predicates. Charge
        // sorting and every active scan as well as real candidates; a
        // pathological overlapping-box arrangement remains budget bounded.
        struct Bounds {double minU,maxU,minV,maxV;};
        std::vector<Bounds> bounds;
        std::vector<std::size_t> order(constraints.size()),active;
        bounds.reserve(constraints.size());active.reserve(constraints.size());
        std::iota(order.begin(),order.end(),std::size_t{0});
        for(const auto [a,b]:constraints){
            tick(1);
            bounds.push_back({std::min(uv[a].u,uv[b].u),std::max(uv[a].u,uv[b].u),
                              std::min(uv[a].v,uv[b].v),std::max(uv[a].v,uv[b].v)});
        }
        std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){
            tick(1);
            return bounds[a].minU<bounds[b].minU||
                   (bounds[a].minU==bounds[b].minU&&a<b);
        });
        for(const auto i:order){
            const auto [a,b]=constraints[i];std::size_t keptActive=0;
            for(const auto j:active){
                tick(1);
                if(bounds[j].maxU<bounds[i].minU)continue;
                active[keptActive++]=j;
                const auto [c,d]=constraints[j];
                if(a==c||a==d||b==c||b==d)continue;
                if(bounds[i].maxV<bounds[j].minV||bounds[j].maxV<bounds[i].minV)continue;
                if(!proper(uv[a],uv[b],uv[c],uv[d]))continue;
                if(constraintTags[i]>0){
                    if(crossingEdges)crossingEdges->insert(constraintTags[i]);
                    if(splitIntervals)splitIntervals->insert(constraintIntervals[i]);
                }
                if(constraintTags[j]>0){
                    if(crossingEdges)crossingEdges->insert(constraintTags[j]);
                    if(splitIntervals)splitIntervals->insert(constraintIntervals[j]);
                }
            }
            active.resize(keptActive);active.push_back(i);
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
        if(ids.size()>=maximumChartNodes)return refuse("chart_node_limit");
        ids.push_back(id);
        uv.push_back(parameter);
    }
    tick(std::uint64_t(uv.size())*uv.size());
    std::vector<Vec2> metricCoordinates;
    if(physicalMetric&&!metricOcctChart(model,faceTag,uv,tick,metricCoordinates))
        return refuse("chart_metric_singular",true);
    const auto& coordinates=physicalMetric?metricCoordinates:uv;
    auto triangles=algo::delaunay2DScalable(coordinates);
    const auto recovery=algo::recoverPlanarConstraints(coordinates,constraints,triangles,
        std::min<std::uint64_t>(5'000'000,remaining()));
    tick(recovery.work);
    if(!recovery.completed){
        switch(recovery.reason){
        case algo::PlanarRecoveryReason::CrossingConstraints:return refuse("planar_crossing_constraints",true);
        case algo::PlanarRecoveryReason::UnsplitSegment:return refuse("planar_unsplit_segment",true);
        case algo::PlanarRecoveryReason::NotRecoverable:return refuse("planar_constraints_not_recoverable",true);
        case algo::PlanarRecoveryReason::WorkBudget:return refuse("planar_recovery_work_budget");
        case algo::PlanarRecoveryReason::AllocationFailed:return refuse("planar_recovery_allocation_failed");
        case algo::PlanarRecoveryReason::InvalidInput:return refuse("planar_recovery_invalid_input");
        case algo::PlanarRecoveryReason::InvalidTriangulation:return refuse("planar_recovery_invalid_triangulation");
        default:return refuse("planar_recovery_internal_error");
        }
    }
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
        if(!parity)return refuse("region_parity_refused");
        std::size_t retained=0;
        for(std::size_t i=0;i<triangles.size();++i){
            tick(1);
            if((*parity)[i])triangles[retained++]=triangles[i];
        }
        triangles.resize(retained);
        if(!refine(uv,ids,triangles,protectedEdges))return refuse("chart_refinement_refused");
    }
    std::map<std::array<Index,2>,unsigned> incidence;
    for(auto triangle:triangles){
        tick(1);
        const auto v=triangle.v;
        const auto a=uv[v[0]],b=uv[v[1]],c=uv[v[2]];
        if(!physicalMetric&&!inDomain({(a.u+b.u+c.u)/3.,(a.v+b.v+c.v)/3.}))continue;
        if(predicates::orient2d(a,b,c)!=predicates::Sign::Positive)return refuse("triangle_orientation");
        for(unsigned i=0;i<3;++i){
            const auto [lo,hi]=std::minmax(v[i],v[(i+1)%3]);
            ++incidence[{lo,hi}];
        }
        bool collapsed=false;
        for(unsigned i=0;i<3;++i)if(ids[v[i]]==ids[v[(i+1)%3]]){
            if(!poleIds.contains(ids[v[i]]))return refuse("unexpected_physical_collapse");
            collapsed=true;
        }
        if(!collapsed)output.push_back({{ids[v[0]],ids[v[1]],ids[v[2]]},{a,b,c}});
    }
    for(const auto& segment:constraints){
        const auto [lo,hi]=std::minmax(segment.first,segment.second);
        if(incidence[{lo,hi}]!=1)return refuse("boundary_constraint_incidence",true);
    }
    if(output.empty())return refuse("chart_empty");
    return true;
}

} // namespace paacs::geo::detail
