// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctCurveMetric.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctStarMesher.hpp"
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Line.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Circ.hxx>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numbers>
#include <sstream>
#include <stdexcept>

using namespace paacs;
using namespace paacs::geo::detail;
namespace {
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
struct WorkStop {};
struct CancelStop {};
struct Work {
    std::uint64_t used=0,maximum=10'000'000,cancelAt=0;
    void operator()(std::uint64_t n){
        if(cancelAt&&used>=cancelAt)throw CancelStop{};
        if(n>maximum-used)throw WorkStop{};
        used+=n;
    }
};
struct Relations {OcctCurveUses uses;OcctVertexEdges incident;};
Relations relations(const geo::OcctModel& model){
    Relations result;
    for(const auto* edge:model.edges()){
        result.uses[edge->tag()];result.incident[edge->startVertexTag()].push_back(edge->tag());
        if(edge->endVertexTag()!=edge->startVertexTag())result.incident[edge->endVertexTag()].push_back(edge->tag());
    }
    for(const auto* solid:model.solids())for(const auto& face:model.faceUses(solid->tag()))
        for(const auto& wire:model.boundaryUses(face.faceTag))for(const auto& use:wire.edges)
            result.uses[use.edgeTag].push_back({face.faceTag,use.useId});
    return result;
}
std::string snapshot(const geo::OcctModel& model){
    std::ostringstream output;BRepTools::Write(OcctModelAccess::shape(model),output);return output.str();
}
int edgeOfType(const geo::OcctModel& model,GeomAbs_CurveType type){
    for(const auto* edge:model.edges())if(BRepAdaptor_Curve(OcctModelAccess::edge(model,edge->tag())).GetType()==type)
        return edge->tag();
    throw std::runtime_error("owned curve fixture is missing its expected native edge");
}
OcctCurveSpacing normal(OcctCurveSamplingReport& report){
    OcctCurveSpacing spacing;spacing.maximumChordError=1e-4;
    spacing.resourcePolicy=geo::CadResourcePolicy::Normal;spacing.maximumSegments=4096;
    spacing.maximumIntegrationIntervals=8192;spacing.maximumWorkingBytes=32ULL*1024*1024;spacing.report=&report;
    return spacing;
}
std::vector<double> sample(const geo::OcctModel& model,int edge,double size,
    const Relations& relation,Work& work,OcctCurveSpacing spacing={}){
    return sampleOcctCurveMetric(model,edge,size,relation.uses,relation.incident,
        [&](std::uint64_t n){work(n);},spacing);
}
void ordered(const std::vector<double>& schedule,double first,double last){
    require(schedule.size()>=2&&schedule.front()==first&&schedule.back()==last,"native curve endpoints changed");
    for(std::size_t i=1;i<schedule.size();++i)
        require(std::isfinite(schedule[i])&&schedule[i]>schedule[i-1],"curve parameters are not strictly ordered");
}
void exerciseCircleAndBoundedIdentity(){
    auto imported=OcctModelAccess::fromShape(BRepPrimAPI_MakeCylinder(1.,2.).Shape());
    require(imported.imported(),"owned cylinder import failed");const auto& model=*imported.model;
    const auto before=snapshot(model);const auto relation=relations(model);
    const auto circle=edgeOfType(model,GeomAbs_Circle),line=edgeOfType(model,GeomAbs_Line);
    Work boundedWork;const auto bounded=sample(model,line,.125,relation,boundedWork);
    double first,last;model.edge(line)->parameterBounds(first,last);
    require(bounded.size()==17&&last-first==2.,"bounded line fixture changed its historical segment count");
    for(std::size_t i=0;i<bounded.size();++i)
        require(std::bit_cast<std::uint64_t>(bounded[i])==
            std::bit_cast<std::uint64_t>(std::lerp(first,last,double(i)/16)),
            "bounded historical binary-exact line schedule changed");
    OcctCurveSpacing oversized;oversized.maximumSegments=10000;oversized.maximumIntegrationIntervals=10000;
    oversized.maximumWorkingBytes=64ULL*1024*1024;Work repeatedWork;
    const auto repeated=sample(model,line,.125,relation,repeatedWork,oversized);
    require(repeated==bounded&&repeatedWork.used==boundedWork.used,"new optional ceilings changed Bounded behavior");
    Work refusedWork;require(sample(model,circle,.01,relation,refusedWork).empty(),
        "fine circle did not expose the historical 256-segment ceiling");
    OcctCurveSamplingReport report;auto spacing=normal(report);Work work;
    const auto fine=sample(model,circle,.01,relation,work,spacing);
    model.edge(circle)->parameterBounds(first,last);ordered(fine,first,last);
    require(fine.size()>257&&report.segments==fine.size()-1&&report.completed&&!report.budgetLimited,
        "Normal curve sampling did not create more than 256 actual segments");
    require(report.integrationIntervals>256&&report.evaluations>257&&
        report.maximumModeledWorkingBytes<=spacing.maximumWorkingBytes,"fine curve bypassed adaptive/resource admission");
    require(std::abs(report.metricLength-2*std::numbers::pi/.01)<1e-8,
        "circle metric does not agree with the independent analytic circumference");
    for(std::size_t i=1;i<fine.size();++i){
        require(fine[i]-fine[i-1]<=.01001,"fine circle arc spacing exceeds its actual request");
        const auto a=model.edge(circle)->pointAt(fine[i-1]),b=model.edge(circle)->pointAt(fine[i]);
        const auto middle=model.edge(circle)->pointAt(std::midpoint(fine[i-1],fine[i]));
        require((middle-(a+b)*.5).norm()<spacing.maximumChordError,"fine circle chord exceeds its retained geometry bound");
    }
    require(snapshot(model)==before,"curve sampling mutated the original CAD");
    std::cout<<"adaptive circle: segments="<<report.segments<<" intervals="<<report.integrationIntervals
        <<" evaluations="<<report.evaluations<<" work="<<work.used<<'\n';
}
const std::array<Vec3,6> poles{{{0,1,3},{1,5,3},{2,.3,3},{3,5,3},{4,.3,3},{5,1,3}}};
double analyticSpeed(double at){
    std::array<Vec3,5> derivative;
    for(unsigned i=0;i<5;++i)derivative[i]=(poles[i+1]-poles[i])*5.;
    for(unsigned degree=4;degree>0;--degree)for(unsigned i=0;i<degree;++i)
        derivative[i]=derivative[i]*(1-at)+derivative[i+1]*at;
    return derivative[0].norm();
}
double referenceLength(double first,double last,unsigned intervals){
    double value=analyticSpeed(first)+analyticSpeed(last);
    for(unsigned i=1;i<intervals;++i)value+=(i%2?4.:2.)*analyticSpeed(std::lerp(first,last,double(i)/intervals));
    return value*(last-first)/(3*intervals);
}
void exerciseVaryingCurve(){
    NCollection_Array1<gp_Pnt> points(1,int(poles.size()));
    for(unsigned i=0;i<poles.size();++i)points(int(i)+1)=gp_Pnt(poles[i].x,poles[i].y,poles[i].z);
    occ::handle<Geom_BezierCurve> curve=new Geom_BezierCurve(points);
    BRep_Builder builder;TopoDS_Compound shape;builder.MakeCompound(shape);
    builder.Add(shape,BRepPrimAPI_MakeBox(1,1,1).Shape());builder.Add(shape,BRepBuilderAPI_MakeEdge(curve).Edge());
    geo::NativeCadImportOptions options;options.retainAuxiliaryGeometry=true;
    auto imported=OcctModelAccess::fromShape(shape,options);require(imported.imported(),"varying-curve import failed");
    const auto& model=*imported.model;const auto before=snapshot(model);
    const auto relation=relations(model);const auto edge=edgeOfType(model,GeomAbs_BezierCurve);
    OcctCurveSamplingReport report;auto spacing=normal(report);spacing.pointsPerTurn=1;spacing.maximumChordError=.01;
    Work work;const auto schedule=sample(model,edge,.01,relation,work,spacing);
    double first,last;model.edge(edge)->parameterBounds(first,last);ordered(schedule,first,last);
    require(schedule.size()>257&&report.completed,"varying curve did not exercise fine adaptive sampling");
    const auto independent=referenceLength(0,1,20000);
    require(std::abs(report.metricLength*.01-independent)<independent*2e-4,
        "adaptive length differs from independent analytic Bezier derivative quadrature");
    for(std::size_t i=1;i<schedule.size();++i)
        require(referenceLength(schedule[i-1],schedule[i],20)<=.0101,
            "varying-speed parameter inversion exceeded requested physical spacing");
    require(report.estimatedMetricError<=report.metricLength*1.1e-4+1e-8&&
        report.maximumModeledWorkingBytes<=spacing.maximumWorkingBytes,"adaptive estimates exceeded their bounds");
    require(before==snapshot(model),"varying-curve sampling modified source geometry");
    std::cout<<"adaptive Bezier: segments="<<report.segments<<" intervals="<<report.integrationIntervals
        <<" reference_length="<<independent<<" work="<<work.used<<'\n';
}
void exerciseUnequalEndpointFeatures(){
    // A wide-parameter straight edge joins circular features with unequal
    // scales. Applying those feature sizes only at the two endpoints caused
    // adaptive relative-error subdivision to run out of parameter precision.
    constexpr double first=-3.5,last=4.5,length=last-first,target=1.25;
    const auto start=BRepBuilderAPI_MakeVertex(gp_Pnt(first,0,0)).Vertex();
    const auto end=BRepBuilderAPI_MakeVertex(gp_Pnt(last,0,0)).Vertex();
    occ::handle<Geom_Line> line=new Geom_Line(gp_Pnt(0,0,0),gp_Dir(1,0,0));
    const auto straight=BRepBuilderAPI_MakeEdge(line,start,end,first,last).Edge();
    const auto arc=[&](double x,double radius,const TopoDS_Vertex& vertex){
        occ::handle<Geom_Circle> curve=new Geom_Circle(gp_Circ(gp_Ax2(gp_Pnt(x,-radius,0),gp_Dir(0,0,1)),radius));
        const auto other=BRepBuilderAPI_MakeVertex(gp_Pnt(x-radius,-radius,0)).Vertex();
        return BRepBuilderAPI_MakeEdge(curve,vertex,other,std::numbers::pi/2,std::numbers::pi).Edge();
    };
    BRep_Builder builder;TopoDS_Compound shape;builder.MakeCompound(shape);
    builder.Add(shape,BRepPrimAPI_MakeBox(gp_Pnt(0,5,0),1,1,1).Shape());
    builder.Add(shape,straight);builder.Add(shape,arc(first,.5,start));builder.Add(shape,arc(last,1.,end));
    geo::NativeCadImportOptions options;options.retainAuxiliaryGeometry=true;
    auto imported=OcctModelAccess::fromShape(shape,options);
    require(imported.imported(),"unequal endpoint fixture import failed");
    const auto& model=*imported.model;const auto before=snapshot(model);const auto relation=relations(model);
    int selected=-1;
    for(const auto* candidate:model.edges()){
        double a,b;candidate->parameterBounds(a,b);
        if(BRepAdaptor_Curve(OcctModelAccess::edge(model,candidate->tag())).GetType()==GeomAbs_Line&&
            a==first&&b==last)selected=candidate->tag();
    }
    require(selected>0,"wide native-parameter line was not retained");
    const auto* source=model.edge(selected);
    require(relation.incident.at(source->startVertexTag()).size()==2&&
        relation.incident.at(source->endVertexTag()).size()==2,"endpoint features lost shared topology");
    OcctCurveSamplingReport report;auto spacing=normal(report);
    spacing.pointsPerTurn=16;spacing.maximumChordError=.02;
    Work work;work.maximum=1'000'000;
    const auto schedule=sample(model,selected,target,relation,work,spacing);
    ordered(schedule,first,last);
    require(report.completed&&!report.budgetLimited&&report.integrationIntervals<512,
        "endpoint size jumps still force unbounded adaptive subdivision");
    const auto slope=std::sqrt(1.8)-1;
    const double near=2*std::numbers::pi*.5/16,far=2*std::numbers::pi/16;
    const auto left=(target-near)/slope,right=length-(target-far)/slope;
    require(left>0&&right>left&&right<length,"independent endpoint envelope fixture changed");
    const auto integral=[&](double x){
        if(x<=left)return std::log1p(slope*x/near)/slope;
        const auto initial=std::log(target/near)/slope;
        if(x<=right)return initial+(x-left)/target;
        return initial+(right-left)/target+std::log(target/(far+slope*(length-x)))/slope;
    };
    const auto reference=integral(length);
    require(std::abs(report.metricLength-reference)<=reference*3e-4,
        "continuous endpoint metric differs from independent analytic line integral");
    const auto metricStep=reference/double(schedule.size()-1);
    double previous=0;
    for(std::size_t i=1;i<schedule.size();++i){
        const auto a=schedule[i-1]-first,b=schedule[i]-first,segment=b-a;
        require(segment<=target*(1+1e-8)&&integral(b)-integral(a)<=metricStep*1.003,
            "endpoint feature size/growth envelope was lost in physical placement");
        if(previous>0)require(std::max(segment/previous,previous/segment)<=std::exp(slope)*1.003,
            "achieved adjacent line segments exceed the independent growth bound");
        previous=segment;
    }
    require(before==snapshot(model)&&report.maximumModeledWorkingBytes<=spacing.maximumWorkingBytes,
        "continuous endpoint sampling changed CAD or escaped modeled memory");
    std::cout<<"continuous endpoint features: segments="<<report.segments
        <<" intervals="<<report.integrationIntervals<<" work="<<work.used<<'\n';
}
void exerciseInteriorCurvatureJump(){
    // C1 quadratic spline: x=t; y=0 before t=.5 and 2*(t-.5)^2
    // afterwards. Speed is continuous while curvature jumps from zero to4.
    NCollection_Array1<gp_Pnt> points(1,4);
    points(1)=gp_Pnt(0,0,0);points(2)=gp_Pnt(.25,0,0);
    points(3)=gp_Pnt(.75,0,0);points(4)=gp_Pnt(1,.5,0);
    NCollection_Array1<double> knots(1,3);knots(1)=0;knots(2)=.5;knots(3)=1;
    NCollection_Array1<int> multiplicity(1,3);multiplicity(1)=3;multiplicity(2)=1;multiplicity(3)=3;
    occ::handle<Geom_BSplineCurve> curve=new Geom_BSplineCurve(points,knots,multiplicity,2,false);
    BRep_Builder builder;TopoDS_Compound shape;builder.MakeCompound(shape);
    builder.Add(shape,BRepPrimAPI_MakeBox(gp_Pnt(0,5,0),1,1,1).Shape());
    builder.Add(shape,BRepBuilderAPI_MakeEdge(curve).Edge());
    geo::NativeCadImportOptions options;options.retainAuxiliaryGeometry=true;
    auto imported=OcctModelAccess::fromShape(shape,options);
    require(imported.imported(),"interior curvature-jump fixture import failed");
    const auto& model=*imported.model;const auto before=snapshot(model);const auto relation=relations(model);
    const auto edge=edgeOfType(model,GeomAbs_BSplineCurve);
    require(model.edge(edge)->curvatureAt(.5-1e-6)<1e-6&&model.edge(edge)->curvatureAt(.5+1e-6)>3.99,
        "owned spline did not reproduce a one-sided curvature jump");
    OcctCurveSamplingReport report;auto spacing=normal(report);
    spacing.pointsPerTurn=16;spacing.maximumChordError=.02;
    Work work;work.maximum=2'000'000;
    const auto schedule=sample(model,edge,.4,relation,work,spacing);ordered(schedule,0,1);
    require(report.completed&&!report.budgetLimited&&report.integrationIntervals<2048,
        "interior curvature jump still exhausts parameter precision or integration resources");
    constexpr std::size_t intervals=20000;
    std::vector<double> arclength(intervals+1),sizes(intervals+1),metric(intervals+1);
    const auto speed=[](double t){const auto q=std::max(0.,t-.5);return std::sqrt(1+16*q*q);};
    for(std::size_t i=0;i<=intervals;++i){
        const auto t=double(i)/intervals,q=std::max(0.,t-.5);
        const auto curvature=t<.5?0.:4/std::pow(1+16*q*q,1.5);
        sizes[i]=curvature>0?std::min({.4,2*std::numbers::pi/(16*curvature),std::sqrt(8*.02/curvature)}):.4;
        if(i){const auto a=double(i-1)/intervals,b=t;
            arclength[i]=arclength[i-1]+(b-a)*(speed(a)+4*speed((a+b)*.5)+speed(b))/6;}
    }
    const auto slope=std::sqrt(1.8)-1;
    // Dense analytic geometry gives an independent discrete Lipschitz minorant
    // of the requested field in physical arclength, including its one-sided
    // knot limit. Simpson quadrature avoids the production logarithmic rule.
    for(std::size_t i=1;i<=intervals;++i)sizes[i]=std::min(sizes[i],sizes[i-1]+slope*(arclength[i]-arclength[i-1]));
    for(std::size_t i=intervals;i>0;--i)sizes[i-1]=std::min(sizes[i-1],sizes[i]+slope*(arclength[i]-arclength[i-1]));
    for(std::size_t i=1;i<=intervals;++i)
        metric[i]=metric[i-1]+(arclength[i]-arclength[i-1])*(1/sizes[i-1]+8/(sizes[i-1]+sizes[i])+1/sizes[i])/6;
    require(std::abs(report.metricLength-metric.back())<=metric.back()*1e-3,
        "prepared spline metric differs from independent dense arclength reference");
    const auto integral=[&](double t){
        if(t>=1)return metric.back();
        const auto scaled=t*intervals;
        const auto index=std::size_t(scaled);
        return std::lerp(metric[index],metric[index+1],scaled-double(index));
    };
    const auto step=metric.back()/double(schedule.size()-1);
    for(std::size_t i=1;i<schedule.size();++i)
        require(integral(schedule[i])-integral(schedule[i-1])<=step*1.01,
            "spline placement misses independent prepared-metric spacing");
    if(report.lastRejectedSizes&&report.lastRejectedRawSizes)
        for(std::size_t i=0;i<5;++i)require((*report.lastRejectedSizes)[i]<=(*report.lastRejectedRawSizes)[i],
            "local growth preparation increased a raw requested size");
    require(before==snapshot(model)&&report.maximumModeledWorkingBytes<=spacing.maximumWorkingBytes,
        "interior curvature sampling changed geometry or escaped modeled memory");
    std::cout<<"interior curvature jump: segments="<<report.segments<<" intervals="<<report.integrationIntervals
        <<" metric="<<report.metricLength<<" reference="<<metric.back()<<" work="<<work.used<<'\n';
}
void exerciseResourceRefusal(){
    auto imported=OcctModelAccess::fromShape(BRepPrimAPI_MakeCylinder(1.,2.).Shape());
    require(imported.imported(),"resource cylinder import failed");const auto& model=*imported.model;
    const auto before=snapshot(model);const auto relation=relations(model);const auto edge=edgeOfType(model,GeomAbs_Circle);
    for(unsigned mode=0;mode<3;++mode){
        OcctCurveSamplingReport report;auto spacing=normal(report);Work work;
        if(mode==0)spacing.maximumSegments=64;
        if(mode==1)spacing.maximumIntegrationIntervals=1;
        if(mode==2)spacing.maximumWorkingBytes=4096;
        require(sample(model,edge,.01,relation,work,spacing).empty()&&report.budgetLimited&&!report.completed,
            "curve resource refusal returned a partial schedule");
        require(report.maximumModeledWorkingBytes<=spacing.maximumWorkingBytes,"curve admission exceeded its explicit memory limit");
    }
    OcctCurveSamplingReport report;auto spacing=normal(report);Work limited;limited.maximum=20;bool stopped=false;
    try{(void)sample(model,edge,.01,relation,limited,spacing);}catch(const WorkStop&){stopped=true;}
    require(stopped&&limited.used<=limited.maximum&&!report.completed,"curve work limit was ignored or masked");
    Work cancelled;cancelled.cancelAt=20;stopped=false;
    try{(void)sample(model,edge,.01,relation,cancelled,spacing);}catch(const CancelStop&){stopped=true;}
    require(stopped&&!report.completed&&before==snapshot(model),"curve cancellation published or mutated geometry");
}
void exerciseManyFaces(){
    constexpr unsigned sides=110;
    BRepBuilderAPI_MakePolygon polygon;
    for(unsigned i=0;i<sides;++i){const auto angle=2*std::numbers::pi*i/sides;polygon.Add(gp_Pnt(3*std::cos(angle),3*std::sin(angle),0));}
    polygon.Close();const auto face=BRepBuilderAPI_MakeFace(polygon.Wire()).Face();
    auto imported=OcctModelAccess::fromShape(BRepPrimAPI_MakePrism(face,gp_Vec(0,0,1)).Shape());
    require(imported.imported(),"many-face prism import failed");const auto& model=*imported.model;const auto before=snapshot(model);
    require(model.faceUses(1).size()>32&&model.edges().size()>256,"many-face fixture did not cross historical count ceilings");
    geo::OcctBoundaryOptions options;options.maximumWork=200'000'000;options.maximumNodes=20000;
    options.maximumTriangles=40000;options.maximumOutputBytes=32ULL*1024*1024;options.maximumChordError=.02;
    for(unsigned curvature:std::array{0u,16u}){
        mesh::Mesh held;held.addNode({99,98,97},{3,77});OcctBoundaryAliasReceipt receipt;
        const auto bounded=buildSpacedOcctBoundary(model,1,options,.9,128ULL*1024*1024,held,receipt,true,false,curvature);
        require(!bounded.published&&bounded.reason==geo::OcctBoundaryReason::BudgetExceeded&&held.nodeCount()==1,
            "Bounded face/edge guard or unchanged destination contract changed");
    }
    options.construction.resourcePolicy=geo::CadResourcePolicy::Normal;
    mesh::Mesh boundary;OcctBoundaryAliasReceipt receipt;
    const auto built=buildSpacedOcctBoundary(model,1,options,.9,128ULL*1024*1024,boundary,receipt,true,false,16);
    require(built.published&&built.closed&&built.oriented&&built.vertexManifold&&built.sourceModelUnmodified,
        "Normal many-face boundary failed unchanged conformity/geometry gates");
    require(receipt.rephasedEdges>256&&receipt.rebuiltCharts>32,"Normal did not process actual edges/faces beyond old ceilings");
    OcctStarOptions audit;audit.maximumWork=100'000'000;audit.maximumWorkingBytes=128ULL*1024*1024;
    const auto certificate=certifyOcctBoundarySurface(boundary,audit);
    require(certificate.certified&&certificate.embedded&&before==snapshot(model),
        "Normal many-face boundary lost exact embedding or source identity");
    std::cout<<"Normal many-face boundary: faces="<<receipt.rebuiltCharts<<" edges="<<receipt.rephasedEdges
        <<" nodes="<<boundary.nodeCount()<<" triangles="<<boundary.cellCount()<<'\n';
}
}
int main(){try{
    exerciseCircleAndBoundedIdentity();exerciseVaryingCurve();exerciseUnequalEndpointFeatures();
    exerciseInteriorCurvatureJump();exerciseResourceRefusal();exerciseManyFaces();
    std::cout<<"Normal CAD curve and face resource guards passed\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
catch(...){std::cerr<<"unexpected CAD curve resource exception\n";return 1;}}
