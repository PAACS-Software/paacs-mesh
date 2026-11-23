// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "OcctCurveMetric.hpp"
#include "OcctChartTrim.hpp"
#include "paacs/algo/Delaunay2D.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/geometry/CadConstructionBudget.hpp"
#include <limits>
#include <stdexcept>
#include <tuple>

namespace paacs::geo::detail {

// Insert locally sized interior points in a constrained chart. A candidate
// replaces its visible, connected in-circle cavity without crossing a trim.
// All point positions come from the original face. The caller checks the
// finished chart's physical facets and the closed shell.
template<class Position, class Append, class Charge, class AdmitStorage>
bool spaceOcctChart(const OcctModel& model,int faceTag,std::vector<Vec2>& uv,
    std::vector<int>& ids,std::vector<algo::Tri>& triangles,
    const std::set<std::array<int,2>>& protectedEdges,double targetSize,double chord,
    Position position,Append append,Charge charge,AdmitStorage admitStorage,bool regularFront=false,
    unsigned pointsPerTurn=16,OcctChartTrim* trimClassifier=nullptr,
    const OcctSpatialSizing* sizeField=nullptr,const CadConstructionBudget& construction={},
    std::uint64_t maximumAdditionalNodes=(std::numeric_limits<std::uint64_t>::max)(),
    std::size_t maximumChartNodes=4096,CadSurfaceSpacingReport* report=nullptr) {
    using Sign=predicates::Sign;
    using Edge=std::array<Index,2>;
    const auto* face=model.face(faceTag);
    const BRepAdaptor_Surface surface(OcctModelAccess::face(model,faceTag),false);
    const bool curved=surface.GetType()!=GeomAbs_Plane;
    const auto finite=[](Vec2 p){return std::isfinite(p.u)&&std::isfinite(p.v);};
    const auto dot=[](Vec2 a,Vec2 b){return a.u*b.u+a.v*b.v;};
    const auto cross=[](Vec2 a,Vec2 b){return a.u*b.v-a.v*b.u;};
    const auto edgeKey=[](Index a,Index b){return Edge{std::min(a,b),std::max(a,b)};};
    const auto protectedEdge=[&](Edge edge){return protectedEdges.contains({int(edge[0]),int(edge[1])});};
    constexpr auto noTriangle=std::numeric_limits<std::size_t>::max();
    constexpr unsigned maximumCavity=96;
    // Extend the size from shared curves toward the target in the interior.
    std::vector<double> curveSize(uv.size(),targetSize);
    struct CurveSegment { Vec3 origin, direction; double squaredLength, size; };
    std::vector<CurveSegment> trim;
    for (const auto& [a,b] : protectedEdges) {
        charge(1);
        const auto length=(position(ids[a])-position(ids[b])).norm();
        if (length>0&&std::isfinite(length)) {
            curveSize[a]=std::min(curveSize[a],length);
            curveSize[b]=std::min(curveSize[b],length);
            if(regularFront) {
                const auto start=position(ids[a]);
                trim.push_back({start,position(ids[b])-start,length*length,length});
            }
        }
    }
    const auto localSize=[&](Vec2 p,Vec3 point) {
        charge(1);
        auto h=targetSize;
        if (sizeField) {
            charge(sizeField->queryWork());
            h=sizeField->eval(point,faceTag);
        }
        const auto curvature=curved?face->curvatureAt(p):0.;
        if (std::isfinite(curvature)&&curvature>0) {
            h=std::min(h,2*std::numbers::pi/(pointsPerTurn*curvature));
            h=std::min(h,std::sqrt(6*chord/curvature));
        }
        if(regularFront) {
            for(const auto& edge:trim) {
                charge(1);
                if(edge.size>=h)continue;
                const auto delta=point-edge.origin;
                const auto t=std::clamp(delta.dot(edge.direction)/edge.squaredLength,0.,1.);
                h=std::min(h,edge.size+.45*(delta-edge.direction*t).norm());
            }
        } else {
            for (std::size_t i=0;i<curveSize.size();++i) {
                charge(1);
                h=std::min(h,curveSize[i]+.45*(point-position(ids[i])).norm());
            }
        }
        return h;
    };
    std::vector<double> nodeSizes;
    if(regularFront)for(std::size_t i=0;i<uv.size();++i)
        nodeSizes.push_back(std::min(curveSize[i],localSize(uv[i],position(ids[i]))));
    const bool progressive=construction.resourcePolicy==CadResourcePolicy::Normal;
    const bool legacy=construction.effort==CadConstructionEffort::Legacy&&!progressive;
    std::uint64_t maximumPoints=construction.surfacePointsPerFace?construction.surfacePointsPerFace:256;
    if(progressive&&!construction.surfacePointsPerFace)maximumPoints=maximumAdditionalNodes;
    else if(!legacy&&!construction.surfacePointsPerFace) {
        // One regular surface vertex occupies approximately sqrt(3)/2*h^2.
        // This is an effort estimate, not a geometric or sizing certificate.
        // Sample the actual retained triangles and requested local field; a
        // heavily refined small region therefore contributes its own density.
        long double desiredVertices=0;
        const auto additionalAllowance=static_cast<long double>(maximumAdditionalNodes);
        for(const auto& triangle:triangles) {
            charge(1);
            const auto a=triangle.v[0],b=triangle.v[1],c=triangle.v[2];
            const auto p=position(ids[a]),q=position(ids[b]),r=position(ids[c]);
            const auto center=(uv[a]+uv[b]+uv[c])*(1./3.);
            auto h=localSize(center,face->pointAt(center));
            if(regularFront)h=std::min({h,nodeSizes[a],nodeSizes[b],nodeSizes[c]});
            const auto twiceArea=(q-p).cross(r-p).norm();
            if(!(std::isfinite(h)&&h>0&&std::isfinite(twiceArea))) {
                desiredVertices=additionalAllowance;
                break;
            }
            desiredVertices+=static_cast<long double>(twiceArea)/(std::sqrt(3.L)*h*h);
            if(desiredVertices>=additionalAllowance+static_cast<long double>(uv.size()))break;
        }
        const auto scale=construction.effort==CadConstructionEffort::Thorough?2.L:1.L;
        const auto need=std::max(0.L,desiredVertices-static_cast<long double>(uv.size()))*scale;
        const auto floor=construction.effort==CadConstructionEffort::Thorough?512ULL:256ULL;
        maximumPoints=need>=additionalAllowance?maximumAdditionalNodes:
            std::max<std::uint64_t>(floor,static_cast<std::uint64_t>(std::ceil(need)));
    }
    // Preserve the exact old ceilings and refusal path for untouched Legacy.
    if(!legacy||construction.surfacePointsPerFace||construction.surfaceChordPointsPerFace)
        maximumPoints=std::min(maximumPoints,maximumAdditionalNodes);
    if(report) {
        report->faceTag=faceTag;
        report->pointLimit=maximumPoints;
        report->chartNodeLimit=maximumChartNodes;
        report->termination="constructing";
    }
    std::map<Edge,double> edgeDeviations;
    const auto edgeDeviation=[&](Edge edge) {
        if (!curved||protectedEdge(edge)) return 0.;
        const auto [entry,inserted]=edgeDeviations.try_emplace(edge);
        if (inserted) {
            charge(1);
            const auto sample=(position(ids[edge[0]])+position(ids[edge[1]]))*.5;
            const auto midpoint=(uv[edge[0]]+uv[edge[1]])*.5;
            auto deviation=(face->pointAt(midpoint)-sample).norm();
            if (deviation>.85*chord) {
                charge(1);
                Vec2 closestParameter;
                Vec3 closestPoint;
                if (face->projectPoint(sample,closestParameter,closestPoint))
                    deviation=(closestPoint-sample).norm();
            }
            entry->second=deviation;
        }
        return entry->second;
    };
    struct Score { double chord=0, spacing=0; };
    // Cell IDs never move. A cavity replaces only its own cells and updates
    // the front status of their neighbours. IDs also preserve tie order.
    struct ChartCell {
        algo::Tri triangle;
        Score score;
        std::array<std::size_t,3> neighbours;
        bool live=true;
        bool blocked=false;
    };
    using Priority=std::tuple<bool,double,std::size_t>;
    std::vector<ChartCell> chart;
    const auto reserveChart=[&](std::size_t needed){
        if(needed<=chart.capacity())return;
        const auto capacity=std::max(needed,chart.capacity()+chart.capacity()/2);
        admitStorage(capacity+chart.capacity());
        chart.reserve(capacity);
    };
    reserveChart(triangles.size()+2*maximumPoints);
    struct EdgeUse { std::size_t triangle=std::numeric_limits<std::size_t>::max(); unsigned side=0; };
    std::map<Edge,std::array<EdgeUse,2>> edgeUses;
    std::set<Priority> front;
    unsigned blocked=0;
    std::uint64_t livePending=0;
    const auto originalNodeCount=uv.size();
    const auto blockCandidate=[&] {
        ++blocked;
        if(report)report->blockedCandidates=blocked;
        return blocked>128;
    };
    const auto scoreTriangle=[&](const std::array<Index,3>& t){
            const auto a=position(ids[t[0]]),b=position(ids[t[1]]),c=position(ids[t[2]]);
            const auto ab=(b-a).norm(),bc=(c-b).norm(),ca=(a-c).norm();
            const auto area2=(b-a).cross(c-a).norm();
            const auto longest=std::max({ab,bc,ca});
            Score score{};
            if(regularFront) {
                constexpr double radiusLimit=.7071067811865476;
                const auto radius=area2>0?ab*bc*ca/(2*area2):longest;
                auto h=(nodeSizes[t[0]]+nodeSizes[t[1]]+nodeSizes[t[2]])/3;
                if (sizeField && !sizeField->boxes().controls().empty()) {
                    charge(sizeField->queryWork());
                    const auto uvCenter=(uv[t[0]]+uv[t[1]]+uv[t[2]])*(1./3.);
                    h=std::min(h,sizeField->eval(face->pointAt(uvCenter),faceTag));
                }
                score.spacing=radius/(radiusLimit*h);
            } else {
                const auto parameter=(uv[t[0]]+uv[t[1]]+uv[t[2]])*(1./3);
                charge(1);
                const auto center=face->pointAt(parameter);
                const auto h=localSize(parameter,center);
                const auto shape=2*std::sqrt(3.)*area2/((b-a).norm2()+(c-b).norm2()+(a-c).norm2());
                const auto deviation=(center-(a+b+c)/3.).norm();
                const auto edgeError=std::max({edgeDeviation(edgeKey(t[0],t[1])),
                    edgeDeviation(edgeKey(t[1],t[2])),edgeDeviation(edgeKey(t[2],t[0]))});
                score={std::max(deviation/(.8*chord),edgeError/(.85*chord)),
                       std::max(longest/(1.25*h),longest/h>.75?.40/shape:0.)};
            }
            return score;
    };
    const auto addCell=[&](algo::Tri triangle){
        const auto index=chart.size();
        chart.push_back({triangle,scoreTriangle(triangle.v),{noTriangle,noTriangle,noTriangle}});
        livePending+=std::max(chart.back().score.chord,chart.back().score.spacing)>1;
        if(report)report->pendingEntries=livePending;
        for(unsigned j=0;j<3;++j) {
            charge(1);
            const auto edge=edgeKey(triangle.v[j],triangle.v[(j+1)%3]);
            auto& uses=edgeUses[edge];
            if(uses[0].triangle==noTriangle) {
                uses[0]={index,j};
                continue;
            }
            if(uses[1].triangle!=noTriangle)return false;
            uses[1]={index,j};
            const auto [other,side]=uses[0];
            chart[index].neighbours[j]=other;
            chart[other].neighbours[side]=index;
        }
        return true;
    };
    for(const auto triangle:triangles) {
        charge(1);
        const auto [a,b,c]=triangle.v;
        if(predicates::orient2d(uv[a],uv[b],uv[c])!=Sign::Positive||!addCell(triangle))return false;
    }
    const auto around=[&](std::size_t i,unsigned side){
        const auto& t=chart[i].triangle.v;
        return protectedEdge(edgeKey(t[side],t[(side+1)%3]))?noTriangle:chart[i].neighbours[side];
    };
    const auto active=[&](std::size_t i,unsigned side){
        const auto other=around(i,side);
        if(other==noTriangle)return true;
        const auto score=chart[other].score;
        return std::max(score.chord,score.spacing)<=1;
    };
    const auto priority=[&](std::size_t i)->Priority{
        const auto score=chart[i].score;
        const bool spacing=!(score.chord>1);
        return {spacing,-(spacing?score.spacing:score.chord),i};
    };
    const auto refresh=[&](std::size_t i){
        charge(1);
        const auto [spacing,score,index]=priority(i);
        if(!(score<-1))return;
        if(chart[index].live&&!chart[index].blocked&&
            (!regularFront||active(index,0)||active(index,1)||active(index,2))) {
            front.emplace(spacing,score,index);
        } else {
            front.erase({spacing,score,index});
        }
    };
    for(std::size_t i=0;i<chart.size();++i)refresh(i);
    std::uint64_t added=0;
    for(;added<maximumPoints&&uv.size()<maximumChartNodes&&!front.empty();) {
        const auto worst=std::get<2>(*front.begin());
        front.erase(front.begin());
        const auto t=chart[worst].triangle.v;
        const auto center=(uv[t[0]]+uv[t[1]]+uv[t[2]])*(1./3);
        charge(1);
        const auto metric=occtChartMetricAt(surface,center);
        if (!metric) {
            chart[worst].blocked=true;
            if (blockCandidate()) break;
            continue;
        }
        const auto x=metric->physical(uv[t[0]]-center);
        const auto y=metric->physical(uv[t[1]]-center);
        const auto z=metric->physical(uv[t[2]]-center);
        const auto u=y-x,v=z-x;
        const auto u2=dot(u,u),v2=dot(v,v),denominator=2*cross(u,v);
        const bool hasCircumcenter=std::isfinite(denominator)&&denominator>0;
        Vec2 circum{};
        if (hasCircumcenter)
            circum=x+Vec2{(v.v*u2-u.v*v2)/denominator,(u.u*v2-v.u*u2)/denominator};
        const auto h=regularFront?0.:localSize(center,face->pointAt(center));
        const auto admit=[&](Vec2 candidate,double separation,std::optional<Edge> split={}) {
            charge(1);
            if (!finite(candidate)||!(trimClassifier?trimClassifier->contains(candidate):face->containsParam(candidate))) return false;
            const auto p=face->pointAt(candidate);
            const auto spacing=localSize(candidate,p);
            if (!(spacing>0&&std::isfinite(spacing))) return false;
            const auto pointMetric=occtChartMetricAt(surface,candidate);
            if (!pointMetric) return false;
            const auto mapped=[&](Index id){return pointMetric->physical(uv[id]-candidate);};
            for (std::size_t i=0;i<uv.size();++i) {
                charge(1);
                const auto minimum=regularFront?separation*std::min(spacing,nodeSizes[i]):separation*spacing;
                if (mapped(static_cast<Index>(i)).norm()<minimum) return false;
            }
            for (const auto& [a,b] : protectedEdges) {
                charge(1);
                const auto q=mapped(a),r=mapped(b),edge=r-q;
                const auto squared=dot(edge,edge);
                if (!(squared>0)) continue;
                const auto projection=std::clamp(-dot(q,edge)/squared,0.,1.);
                const auto clearance=regularFront?.32*std::min(spacing,std::lerp(nodeSizes[a],nodeSizes[b],projection)):
                    .45*separation*spacing;
                if ((q+edge*projection).norm()<clearance) return false;
            }
            std::vector<std::size_t> cavity;
            if (split) {
                const auto& uses=edgeUses.at(*split);
                if (uses[1].triangle==noTriangle||protectedEdge(*split)) return false;
                const auto [a,b]=std::minmax(uses[0].triangle,uses[1].triangle);
                cavity={a,b};
            } else {
                const auto sideOf=[&](std::size_t i,unsigned j){
                    charge(1);
                    const auto& a=chart[i].triangle.v;
                    return predicates::orient2d(uv[a[j]],uv[a[(j+1)%3]],candidate);
                };
                // A strictly interior hit is unique. Keep stable insertion
                // order for a point exactly on an edge or a blocked walk.
                auto walk=worst;
                for(unsigned step=0;step<32&&walk!=noTriangle;++step) {
                    const std::array signs{sideOf(walk,0),sideOf(walk,1),sideOf(walk,2)};
                    const auto negative=std::ranges::find(signs,Sign::Negative);
                    if(negative==signs.end()) {
                        if(std::ranges::find(signs,Sign::Zero)==signs.end())cavity.push_back(walk);
                        break;
                    }
                    walk=around(walk,static_cast<unsigned>(negative-signs.begin()));
                }
                if(cavity.empty())
                    for(std::size_t i=0;i<chart.size();++i) {
                        charge(1);
                        if(chart[i].live&&sideOf(i,0)!=Sign::Negative&&sideOf(i,1)!=Sign::Negative&&
                            sideOf(i,2)!=Sign::Negative) {
                            cavity.push_back(i);
                            break;
                        }
                    }
                if (cavity.empty()) return false;
                std::vector<bool> visited(chart.size());
                visited[cavity.front()]=true;
                for (std::size_t i=0;i<cavity.size();++i) for (unsigned j=0;j<3;++j) {
                    charge(1);
                    const auto other=around(cavity[i],j);
                    if (other==noTriangle||visited[other]) continue;
                    visited[other]=true;
                    const auto& a=chart[other].triangle.v;
                    if (predicates::incircle(mapped(a[0]),mapped(a[1]),mapped(a[2]),{})!=Sign::Negative) {
                        if (cavity.size()==maximumCavity) return false;
                        cavity.push_back(other);
                    }
                }
            }
            if (std::ranges::find(cavity,worst)==cavity.end()) return false;
            std::vector<bool> removed(chart.size());
            for (const auto i : cavity) removed[i]=true;
            std::vector<algo::Tri> replacement;
            std::set<Index> outgoing,incoming;
            const auto next=static_cast<Index>(uv.size());
            for (const auto i : cavity) for (unsigned j=0;j<3;++j) {
                charge(1);
                const auto other=around(i,j);
                if (other!=noTriangle&&removed[other]) continue;
                const auto& a=chart[i].triangle.v;
                if (predicates::orient2d(uv[a[j]],uv[a[(j+1)%3]],candidate)!=Sign::Positive)
                    return false;
                if (!outgoing.insert(a[j]).second||!incoming.insert(a[(j+1)%3]).second)
                    return false;
                replacement.push_back({{a[j],a[(j+1)%3],next}});
            }
            // A simple planar cavity with no swallowed vertex has k+2 shell
            // edges. Keep each old node and refuse a disconnected or pinched fan.
            if (replacement.size()!=cavity.size()+2) return false;
            reserveChart(chart.size()+replacement.size());
            const auto id=append(candidate,p);
            uv.push_back(candidate);
            ids.push_back(id);
            if(report)report->addedPoints=uv.size()-originalNodeCount;
            if(regularFront)nodeSizes.push_back(spacing);
            std::vector<std::size_t> neighbours;
            for(const auto i:cavity) {
                livePending-=std::max(chart[i].score.chord,chart[i].score.spacing)>1;
                if(report)report->pendingEntries=livePending;
                chart[i].live=false;
                refresh(i);
                const auto& vertices=chart[i].triangle.v;
                for(unsigned j=0;j<3;++j) {
                    charge(1);
                    const auto at=edgeUses.find(edgeKey(vertices[j],vertices[(j+1)%3]));
                    if(at==edgeUses.end())throw std::logic_error("missing chart edge use");
                    auto& uses=at->second;
                    const unsigned slot=uses[0].triangle==i?0:1;
                    if(uses[slot].triangle!=i)throw std::logic_error("missing chart edge owner");
                    const auto other=uses[1-slot];
                    if(other.triangle==noTriangle) {
                        edgeUses.erase(at);
                    } else {
                        chart[other.triangle].neighbours[other.side]=noTriangle;
                        if(!removed[other.triangle])neighbours.push_back(other.triangle);
                        uses={other,EdgeUse{}};
                    }
                }
            }
            const auto first=chart.size();
            for(const auto triangle:replacement)
                if(!addCell(triangle))throw std::logic_error("nonmanifold chart cavity");
            for(const auto i:neighbours)refresh(i);
            for(auto i=first;i<chart.size();++i)refresh(i);
            return true;
        };
        const auto insert=[&]() {
            const std::array xy{x,y,z};
            if(regularFront) {
                constexpr double altitude=.8660254037844386;
                for(unsigned j=0;j<3;++j) {
                    if(!active(worst,j))continue;
                    const auto mid=(xy[j]+xy[(j+1)%3])*.5;
                    const auto direction=circum-mid;
                    const auto length=direction.norm();
                    if(!hasCircumcenter||!(length>0)||cross(xy[(j+1)%3]-xy[j],direction)<=0)continue;
                    const auto size=(nodeSizes[t[j]]+nodeSizes[t[(j+1)%3]])*.5;
                    const auto at=[&](double height){return center+metric->parameter(mid+direction*(height/length));};
                    const auto firstHeight=std::min(length,altitude*size);
                    const auto apex=at(firstHeight);
                    const auto nextSize=localSize(apex,face->pointAt(apex));
                    const auto height=std::min(firstHeight,altitude*nextSize);
                    if(admit(at(height),.52))return true;
                }
                return false;
            }
            // Chord accuracy can require a smaller point spacing than the
            // size field. Split the offending interior edge before grading.
            for (unsigned j=0;j<3;++j) {
                const auto edge=edgeKey(t[j],t[(j+1)%3]);
                if (edgeDeviation(edge)<=.85*chord) continue;
                const auto midpoint=(uv[edge[0]]+uv[edge[1]])*.5;
                if (admit(midpoint,.10)||admit(midpoint,.08,edge))
                    return true;
            }
            for (unsigned j=0;j<3;++j) if (active(worst,j)) {
                const auto a=xy[j],b=xy[(j+1)%3],edge=b-a;
                const auto length=edge.norm();
                if (!(length>0)) continue;
                const auto midpoint=(a+b)*.5;
                const auto normal=Vec2{-edge.v,edge.u}*(1./length);
                const auto height=std::clamp(.8660254037844386*length,.45*h,.95*h);
                const auto towards=hasCircumcenter?dot(circum-midpoint,normal):0.;
                const auto offset=towards>.35*h?std::min(height,towards):height;
                if (admit(center+metric->parameter(midpoint+normal*offset),.24)) return true;
            }
            if (hasCircumcenter&&admit(center+metric->parameter(circum),.24)) return true;
            if (admit(center,.24)) return true;
            for (unsigned j=0;j<3;++j)
                if (admit((uv[t[j]]+uv[t[(j+1)%3]])*.5,.24)) return true;
            return false;
        };
        if (insert()) {
            ++added;
        } else {
            chart[worst].blocked=true;
            if (blockCandidate()) break;
        }
    }
    triangles.clear();
    std::uint64_t pending=0;
    for(const auto& cell:chart) {
        charge(1);
        if(cell.live) {
            triangles.push_back(cell.triangle);
            pending+=std::max(cell.score.chord,cell.score.spacing)>1;
        }
    }
    if(report) {
        report->addedPoints=added;
        report->pendingEntries=pending;
        report->blockedCandidates=blocked;
        report->queueExhausted=pending==0;
        report->budgetLimited=pending&&(added==maximumPoints||uv.size()>=maximumChartNodes||blocked>128);
        report->termination=!pending?"queue_exhausted":uv.size()>=maximumChartNodes?"chart_node_limit":
            added==maximumPoints?"point_limit":blocked>128?"blocked_candidate_limit":"stalled";
    }
    return true;
}

} // namespace paacs::geo::detail
