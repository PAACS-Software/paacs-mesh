// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctCurveMetric.hpp"
#include "normal_span_knot_prism.hpp"
#include "curve_span_fixtures.hpp"
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
#include <TopoDS.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Circ.hxx>
#include <gp_Ax1.hxx>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string_view>

using namespace paacs;
using namespace paacs::geo::detail;
namespace {
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
void exerciseStencilSeedCollision(){
    const auto edge=BRepBuilderAPI_MakeEdge(test::normalSpanKnotCurve()).Edge();
    BRepAdaptor_Curve curve(edge);
    const auto knot=std::nextafter(.5,std::numeric_limits<double>::infinity());
    require(curve.FirstParameter()==0.&&curve.LastParameter()==1.,"owned edge retains ordinary parameter range");
    require(curve.NbIntervals(GeomAbs_C2)==2,"owned spline retains mandatory C2 split");
    NCollection_Array1<double> intervals(1,3);curve.Intervals(intervals,GeomAbs_C2);
    require(intervals(1)==0.&&intervals(2)==knot&&intervals(3)==1.,"native continuity intervals retain adjacent knot");
    for(double parameter:std::array{0.,.125,.25,.5,knot,.75,.875,1.}){
        gp_Pnt point;gp_Vec derivative;curve.D1(parameter,point,derivative);
        require(std::abs(point.X()-parameter)<1e-14&&point.Y()==0.&&point.Z()==0.&&
              std::abs(derivative.X()-1.)<1e-13&&derivative.Y()==0.&&derivative.Z()==0.,
              "Greville poles reproduce physical straight edge and ordinary derivative");
    }
    const auto middle=std::midpoint(.5,knot);
    require(!(std::midpoint(.5,middle)>.5&&middle< std::midpoint(middle,knot)&&std::midpoint(middle,knot)<knot),
          "adjacent parameters cannot form five distinct stencil points");
    geo::detail::OcctCurveSamplingReport report;
    geo::detail::OcctCurveSpacing spacing;
    spacing.resourcePolicy=geo::CadResourcePolicy::Normal;spacing.maximumSegments=10000;
    spacing.maximumIntegrationIntervals=10000;spacing.maximumWorkingBytes=16ULL*1024*1024;
    spacing.maximumChordError=.005;spacing.report=&report;
    std::uint64_t work=0;
    const auto result=geo::detail::curve_metric_detail::progressive(curve,0.,1.,false,spacing,
        [&](double parameter){
            gp_Pnt point;gp_Vec derivative;curve.D1(parameter,point,derivative);
            return geo::detail::curve_metric_detail::Sample{parameter,derivative.Magnitude(),.3,
                {point.X(),point.Y(),point.Z()}};
        },[&](std::uint64_t n){work+=n;});
    require(report.completed&&!report.budgetLimited&&std::string_view(report.diagnostic)=="completed"&&result.size()==5,
          "artificial cut next to native knot no longer creates an unrepresentable interval");
    for(std::size_t i=0;i<result.size();++i)require(std::abs(result[i]-double(i)/4.)<1e-12,
          "corrected line schedule retains uniform physical spacing");
    require(!report.numericalInterval&&work>0&&report.evaluations>0&&work<10000,
          "corrected seed selection is bounded and reports no false numerical refusal");

    const auto mandatoryEdge=BRepBuilderAPI_MakeEdge(test::adjacentMandatoryKnotCurve()).Edge();
    BRepAdaptor_Curve mandatory(mandatoryEdge);
    require(mandatory.NbIntervals(GeomAbs_C2)==3,"two native continuity knots remain mandatory");
    NCollection_Array1<double> mandatoryIntervals(1,4);mandatory.Intervals(mandatoryIntervals,GeomAbs_C2);
    const auto upper=mandatoryIntervals(3);
    require(mandatoryIntervals(2)==.5&&upper>.5&&!geo::detail::curve_metric_detail::stencil(.5,upper).resolved(),
          "owned mandatory interval is not stencil-representable");
    work=0;
    const auto refused=geo::detail::curve_metric_detail::progressive(mandatory,0.,1.,false,spacing,
        [&](double parameter){
            gp_Pnt point;gp_Vec derivative;mandatory.D1(parameter,point,derivative);
            return geo::detail::curve_metric_detail::Sample{parameter,derivative.Magnitude(),.3,
                {point.X(),point.Y(),point.Z()}};
        },[&](std::uint64_t n){work+=n;});
    require(refused.empty()&&!report.completed&&!report.budgetLimited&&
          std::string_view(report.diagnostic)=="curve_parameter_resolution","mandatory unresolved interval still refuses honestly");
    require(report.numericalInterval&&(*report.numericalInterval)[0]==.5&&(*report.numericalInterval)[1]==upper,
          "mandatory knot interval is retained without epsilon merging");
    require(work>0&&work<10000,"mandatory numerical refusal retains bounded measured work");
}

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
std::string shapeSnapshot(const TopoDS_Shape& shape){
    std::ostringstream text;BRepTools::Write(shape,text);return text.str();
}
std::vector<double> sampleOwnedSpan(BRepAdaptor_Curve& curve,bool varyingSecond,
    OcctCurveSamplingReport& report,Work& work,OcctCurveSpacing spacing={}){
    if(!spacing.report)spacing=normal(report);
    return curve_metric_detail::progressive(curve,0.,1.,false,spacing,
        [&](double at,const curve_metric_detail::NativeSpan& span){
            gp_Pnt point,secondPoint;gp_Vec derivative,secondFirst,second;
            span.D1(at,point,derivative);span.D2(at,secondPoint,secondFirst,second);
            const bool left=at<.5||(at==.5&&span.last==.5);
            const auto q=at-.5;
            const double x=varyingSecond?(left?.25*at+.5*at*at:.25+q+q*q):
                (left?.5*at:.25+1.5*q);
            const double first=varyingSecond?(left?.25+at:1+2*q):(left?.5:1.5);
            const double dd=varyingSecond?(left?1.:2.):0.;
            gp_Pnt expected(x,0,0);gp_Vec expectedFirst(first,0,0),expectedSecond(dd,0,0);
            expected.Transform(curve.Trsf());expectedFirst.Transform(curve.Trsf());expectedSecond.Transform(curve.Trsf());
            require(point.Distance(expected)<1e-13&&secondPoint.Distance(expected)<1e-13&&
                (derivative-expectedFirst).Magnitude()<1e-13&&(secondFirst-expectedFirst).Magnitude()<1e-13&&
                (second-expectedSecond).Magnitude()<1e-13,"span derivatives or placement disagree with analytic piecewise polynomial");
            return curve_metric_detail::Sample{at,derivative.Magnitude(),.25,{point.X(),point.Y(),point.Z()}};
        },[&](std::uint64_t n){work(n);});
}
void exerciseOneSidedSpans(){
    const auto spline=test::oneSidedSpeedCurve();
    const auto edge=BRepBuilderAPI_MakeEdge(spline).Edge();const auto before=shapeSnapshot(edge);
    BRepAdaptor_Curve curve(edge);OcctCurveSamplingReport report;Work work;
    const auto schedule=sampleOwnedSpan(curve,false,report,work);
    ordered(schedule,0,1);
    require(report.completed&&report.nativeSpans==2&&report.localD1==2&&report.localD2==2&&
        report.adapterD1>0&&report.adapterD2>0&&schedule.size()==5&&std::abs(report.metricLength-4.)<1e-12,
        "side-owned speed metric differs from exact unit arclength or overrides ordinary endpoints");
    for(std::size_t i=0;i<schedule.size();++i)
        require(std::abs(curve.Value(schedule[i]).X()-double(i)/4.)<1e-12,
            "one-sided speed mapping does not give equal physical spacing");
    require(work.used<10000&&report.maximumModeledWorkingBytes<=normal(report).maximumWorkingBytes,
        "one-sided construction escaped admitted work or memory");
    gp_Trsf placement;placement.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1)),.61);
    placement.SetTranslationPart(gp_Vec(2,-3,4));
    const auto moved=TopoDS::Edge(edge.Moved(TopLoc_Location(placement)));const auto movedBefore=shapeSnapshot(moved);
    BRepAdaptor_Curve located(moved);OcctCurveSamplingReport locatedReport;Work locatedWork;
    const auto placed=sampleOwnedSpan(located,false,locatedReport,locatedWork);
    require(locatedReport.completed&&placed==schedule&&locatedReport.metricLength==report.metricLength,
        "native placement changed one-sided metric or parameter schedule");
    const auto varyingEdge=BRepBuilderAPI_MakeEdge(test::oneSidedSpeedCurve(true)).Edge();
    BRepAdaptor_Curve varying(varyingEdge);OcctCurveSamplingReport varyingReport;Work varyingWork;
    const auto secondSchedule=sampleOwnedSpan(varying,true,varyingReport,varyingWork);
    require(varyingReport.completed&&varyingReport.localD2==2&&std::abs(varyingReport.metricLength-4.)<1e-12,
        "one-sided nonzero second derivatives changed exact physical length");
    ordered(secondSchedule,0,1);
    for(unsigned mode=0;mode<3;++mode){
        OcctCurveSamplingReport limitedReport;auto spacing=normal(limitedReport);Work limited;
        if(mode==0)limited.maximum=30;
        if(mode==1)limited.cancelAt=30;
        if(mode==2)spacing.maximumWorkingBytes=4096;
        bool workStopped=false,cancelled=false;
        try{
            const auto refused=sampleOwnedSpan(curve,false,limitedReport,limited,spacing);
            require(mode==2&&refused.empty()&&limitedReport.budgetLimited,
                "span resource admission unexpectedly returned a completed schedule");
        }catch(const WorkStop&){workStopped=true;}catch(const CancelStop&){cancelled=true;}
        require(!limitedReport.completed&&(mode!=0||(workStopped&&limited.used<=limited.maximum))&&
            (mode!=1||cancelled)&&limitedReport.maximumModeledWorkingBytes<=spacing.maximumWorkingBytes,
            "span work/cancellation/memory refusal was hidden or exceeded admission");
    }
    require(before==shapeSnapshot(edge)&&movedBefore==shapeSnapshot(moved),
        "one-sided sampling changed original or placed source geometry");
    std::cout<<"one-sided speed: segments="<<schedule.size()-1<<" metric="<<report.metricLength
        <<" work="<<work.used<<" internal_sides="<<report.localD1<<'\n';
}
void exerciseOrdinarySpanReference(){
    const auto edge=BRepBuilderAPI_MakeEdge(test::smoothSpanCurve()).Edge();const auto before=shapeSnapshot(edge);
    BRepAdaptor_Curve native(edge),owned(edge);OcctCurveSamplingReport reference,candidate;
    auto referenceSpacing=normal(reference),candidateSpacing=normal(candidate);
    const auto scalar=curve_metric_detail::progressive(native,0.,1.,false,referenceSpacing,
        [&](double at){gp_Pnt p;gp_Vec d;native.D1(at,p,d);return curve_metric_detail::Sample{at,d.Magnitude(),.25,{p.X(),p.Y(),p.Z()}};},[](std::uint64_t){});
    const auto result=curve_metric_detail::progressive(owned,0.,1.,false,candidateSpacing,
        [&](double at,const curve_metric_detail::NativeSpan& span){gp_Pnt p;gp_Vec d;span.D1(at,p,d);return curve_metric_detail::Sample{at,d.Magnitude(),.25,{p.X(),p.Y(),p.Z()}};},[](std::uint64_t){});
    require(reference.completed&&candidate.completed&&result==scalar&&candidate.metricLength==reference.metricLength&&
        candidate.nativeSpans==1&&candidate.localD1==0&&before==shapeSnapshot(edge),
        "ordinary smooth native reference changed or received an unnecessary local derivative override");
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
    require(report.localD2>=2,"model-aware curvature evaluation did not use both internal knot sides");
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
    exerciseStencilSeedCollision();exerciseOneSidedSpans();exerciseOrdinarySpanReference();
    exerciseCircleAndBoundedIdentity();exerciseVaryingCurve();exerciseUnequalEndpointFeatures();
    exerciseInteriorCurvatureJump();exerciseResourceRefusal();exerciseManyFaces();
    std::cout<<"Normal CAD curve and face resource guards passed\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
catch(...){std::cerr<<"unexpected CAD curve resource exception\n";return 1;}}
