// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "OcctModelInternal.hpp"
#include "OcctSpatialSizing.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <Geom_BSplineCurve.hxx>
#include <NCollection_Array1.hxx>
#include <Standard_OutOfMemory.hxx>
#include <Standard_Failure.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <limits>
#include <numbers>
#include <new>
#include <numeric>
#include <optional>
#include <set>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace paacs::geo::detail {

using OcctCurveUses = std::map<int, std::vector<std::pair<int, std::uint64_t>>>;
using OcctVertexEdges = std::map<int, std::vector<int>>;

struct OcctCurveSamplingReport {
    std::uint64_t evaluations=0,integrationIntervals=0,segments=0,maximumModeledWorkingBytes=0;
    double metricLength=0,estimatedMetricError=0;
    bool completed=false,budgetLimited=false;
    const char* diagnostic="not_attempted";
    std::optional<std::array<double,2>> numericalInterval,lastRejectedInterval;
    std::optional<std::array<double,5>> lastRejectedSpeeds,lastRejectedSizes,lastRejectedRawSizes;
    // Metric, length and physical-mapping error estimates and their tolerances.
    std::optional<std::array<double,6>> lastRejectedErrors;
    std::uint64_t nativeSpans=0,localD1=0,localD2=0,adapterD1=0,adapterD2=0;
};
struct OcctCurveSpacing {
    unsigned pointsPerTurn = 32;
    double maximumChordError = 0;
    CadResourcePolicy resourcePolicy=CadResourcePolicy::Bounded;
    std::uint64_t maximumSegments=256,maximumIntegrationIntervals=256;
    std::uint64_t maximumWorkingBytes=1024*1024;
    OcctCurveSamplingReport* report=nullptr;
};

namespace curve_metric_detail {
struct Sample {double at,speed,size;Vec3 point;};
struct SamplingStop {};
// One private native adapter per mandatory C2 span. Only exact internal span
// endpoints need knot-local derivatives; ordinary interior evaluation and the
// original edge endpoints retain their native adapter. Placement is preserved.
struct NativeSpan {
    double first=0.,last=0.;
    occ::handle<Geom_BSplineCurve> spline;
    occ::handle<Adaptor3d_Curve> adapter;
    gp_Trsf placement;
    int firstKnot=0,lastKnot=0;
    bool internalFirst=false,internalLast=false;
    OcctCurveSamplingReport* report=nullptr;
    void D1(double at,gp_Pnt& point,gp_Vec& derivative)const {
        const bool atFirst=at==first&&internalFirst,atLast=at==last&&internalLast;
        if(spline.IsNull()||(!atFirst&&!atLast)){if(report)++report->adapterD1;adapter->D1(at,point,derivative);return;}
        if(report)++report->localD1;
        const auto from=atFirst?firstKnot:std::max(firstKnot,lastKnot-1);
        const auto to=atFirst?std::min(lastKnot,firstKnot+1):lastKnot;
        spline->LocalD1(at,from,to,point,derivative);
        point.Transform(placement);derivative.Transform(placement);
    }
    void D2(double at,gp_Pnt& point,gp_Vec& derivative,gp_Vec& second)const {
        const bool atFirst=at==first&&internalFirst,atLast=at==last&&internalLast;
        if(spline.IsNull()||(!atFirst&&!atLast)){if(report)++report->adapterD2;adapter->D2(at,point,derivative,second);return;}
        if(report)++report->localD2;
        const auto from=atFirst?firstKnot:std::max(firstKnot,lastKnot-1);
        const auto to=atFirst?std::min(lastKnot,firstKnot+1):lastKnot;
        spline->LocalD2(at,from,to,point,derivative,second);
        point.Transform(placement);derivative.Transform(placement);second.Transform(placement);
    }
};
struct Stencil {
    double first,quarter,middle,threeQuarter,last;
    [[nodiscard]] bool resolved()const noexcept {
        return first<quarter&&quarter<middle&&middle<threeQuarter&&threeQuarter<last;
    }
};
inline Stencil stencil(double first,double last) noexcept {
    const auto middle=std::midpoint(first,last);
    return {first,std::midpoint(first,middle),middle,std::midpoint(middle,last),last};
}

// A numerical metric estimate, not a CAD-error certificate. Native face/curve
// identity, trim, chord and final mesh audits remain the consumer's obligation.
// Cache adaptive samples once, then prepare the same bidirectional growth field
// along approximate physical arclength. Integrate its piecewise linear size
// exactly with the logarithmic mean rather than a fixed parameter grid.
template<class Evaluate,class Tick>
std::vector<double> progressive(BRepAdaptor_Curve& curve,double first,double last,
    bool closed,const OcctCurveSpacing& spacing,Evaluate evaluate,Tick tick) {
    OcctCurveSamplingReport localReport;
    auto& report=spacing.report?*spacing.report:localReport;report={};report.diagnostic="sampling";
    const auto fail=[&](const char* why,bool budget){
        report.diagnostic=why;report.budgetLimited=budget;throw SamplingStop{};
    };
    const auto segments=std::min<std::uint64_t>(spacing.maximumSegments,
        std::uint64_t((std::numeric_limits<unsigned>::max)())-1);
    const auto intervals=std::min<std::uint64_t>(spacing.maximumIntegrationIntervals,
        (std::uint64_t((std::numeric_limits<std::size_t>::max)())-1)/5);
    try{
        if(!segments||!intervals||spacing.maximumWorkingBytes<4096)fail("curve_resource_limit",true);
        if(!std::isfinite(first)||!std::isfinite(last)||!(last>first)||!std::isfinite(last-first))
            fail("curve_parameter_range",false);
        const auto span=last-first;
        std::map<std::pair<double,std::size_t>,Sample> cache;
        std::vector<double> breaks,nativeBreaks;
        std::vector<NativeSpan> spans;
        std::uint64_t spanSlots=0;
        static_assert(sizeof(Sample)<=64,"adaptive curve sample payload model");
        struct Interval{double first,last;std::size_t span;};
        std::vector<Interval> pending;
        const auto admit=[&](std::uint64_t samples,std::uint64_t stackSlots,
                             std::uint64_t breakSlots,std::uint64_t finalSamples,std::uint64_t outputSlots){
            std::uint64_t bytes=4096;
            for(const auto [count,unit]:std::array{
                std::pair{spanSlots,std::uint64_t(512)},
                std::pair{samples,std::uint64_t(256)},std::pair{stackSlots,std::uint64_t(64)},
                std::pair{breakSlots,std::uint64_t(32)},std::pair{finalSamples,std::uint64_t(192)},
                std::pair{outputSlots,std::uint64_t(32)}}){
                if(count>(spacing.maximumWorkingBytes-bytes)/unit)fail("curve_working_memory",true);
                bytes+=count*unit;
            }
            report.maximumModeledWorkingBytes=std::max(report.maximumModeledWorkingBytes,bytes);
        };
        const auto nativeIntervals=curve.NbIntervals(GeomAbs_C2);
        if(nativeIntervals<1)fail("curve_continuity_intervals",false);
        if(std::uint64_t(nativeIntervals)>intervals||nativeIntervals==(std::numeric_limits<int>::max)())
            fail("curve_integration_limit",true);
        spanSlots=std::uint64_t(nativeIntervals);
        admit(0,0,2*(std::uint64_t(nativeIntervals)+17),0,0);
        {
            NCollection_Array1<double> native(1,nativeIntervals+1);curve.Intervals(native,GeomAbs_C2);
            breaks.reserve(std::size_t(nativeIntervals)+18);
            breaks.push_back(first);breaks.push_back(last);
            for(int i=1;i<=nativeIntervals+1;++i)if(native(i)>first&&native(i)<last)breaks.push_back(native(i));
        }
        std::sort(breaks.begin(),breaks.end());breaks.erase(std::unique(breaks.begin(),breaks.end()),breaks.end());
        nativeBreaks=breaks;
        spans.reserve(nativeBreaks.size()-1);
        const auto rawSpline=curve.Is3DCurve()&&curve.GetType()==GeomAbs_BSplineCurve?
            curve.Curve().BSpline():occ::handle<Geom_BSplineCurve>{};
        for(std::size_t i=1;i<nativeBreaks.size();++i){
            tick(1);NativeSpan local;local.first=nativeBreaks[i-1];local.last=nativeBreaks[i];
            local.internalFirst=local.first>first;local.internalLast=local.last<last;local.report=&report;
            local.placement=curve.Trsf();
            const auto owned=curve.ShallowCopy();
            local.adapter=owned->Trim(local.first,local.last,0.);
            if(local.adapter.IsNull())fail("curve_local_span",false);
            if(!rawSpline.IsNull()&&!rawSpline->HasEvalRepresentation()&&local.first>=rawSpline->Knot(1)&&
               local.last<=rawSpline->Knot(rawSpline->NbKnots())){
                const auto count=rawSpline->NbKnots();
                if(count==(std::numeric_limits<int>::max)())fail("curve_native_knot_limit",true);
                int lo=1,hi=count+1;
                while(lo<hi){tick(1);const auto middle=lo+(hi-lo)/2;
                    if(rawSpline->Knot(middle)<=local.first)lo=middle+1;else hi=middle;}
                local.firstKnot=std::max(1,lo-1);
                lo=1;hi=count+1;
                while(lo<hi){tick(1);const auto middle=lo+(hi-lo)/2;
                    if(rawSpline->Knot(middle)<local.last)lo=middle+1;else hi=middle;}
                local.lastKnot=std::min(count,lo);
                if(local.firstKnot>=local.lastKnot)fail("curve_local_span",false);
                local.spline=rawSpline;
            }
            spans.push_back(std::move(local));
        }
        report.nativeSpans=spans.size();
        // Native continuity breaks and endpoints are mandatory. Artificial
        // integration seeds must not create an unrepresentable interval next
        // to one of them. Preserve every native break, and admit each sorted
        // seed only when both children support the exact integration stencil.
        const auto initial=std::min<std::uint64_t>(16,intervals);
        for(std::uint64_t i=1;i<initial;++i){
            const auto cut=std::lerp(first,last,double(i)/initial);
            const auto at=std::lower_bound(breaks.begin(),breaks.end(),cut,
                [&](double value,double wanted){tick(1);return value<wanted;});
            if(at==breaks.begin()||at==breaks.end()||*at==cut)continue;
            tick(1);
            if(!stencil(*(at-1),cut).resolved()||!stencil(cut,*at).resolved())continue;
            // Charge the bounded search above and every shifted entry before
            // inserting. Initial reserve and its memory admission include all
            // possible seeds; no new allocation policy is introduced here.
            tick(std::uint64_t(breaks.end()-at)+1);
            breaks.insert(at,cut);
        }
        if(breaks.size()-1>intervals)fail("curve_integration_limit",true);
        admit(0,2*breaks.size(),breaks.capacity(),0,0);
        pending.reserve(breaks.size());
        for(std::size_t i=breaks.size()-1;i>0;--i){
            const auto owner=std::upper_bound(nativeBreaks.begin(),nativeBreaks.end(),breaks[i-1],
                [&](double wanted,double value){tick(1);return wanted<value;});
            const auto id=std::size_t(owner-nativeBreaks.begin()-1);
            if(id>=spans.size()||breaks[i]>spans[id].last)fail("curve_local_span",false);
            pending.push_back({breaks[i-1],breaks[i],id});
        }
        const auto sample=[&](double at,std::size_t spanId)->const Sample&{
            const auto key=std::pair{at,spanId};
            if(auto found=cache.find(key);found!=cache.end())return found->second;
            if(cache.size()>=4*intervals+spans.size())fail("curve_integration_limit",true);
            admit(cache.size()+1,2*pending.capacity(),breaks.capacity(),0,0);
            tick(1);const auto value=[&]{
                if constexpr(std::is_invocable_v<Evaluate,double,const NativeSpan&>)return evaluate(at,spans[spanId]);
                else return evaluate(at);
            }();++report.evaluations;
            if(!(std::isfinite(value.speed)&&value.speed>=0&&std::isfinite(value.size)&&value.size>0&&
                 std::isfinite(value.point.x)&&std::isfinite(value.point.y)&&std::isfinite(value.point.z)))
                fail("curve_nonfinite_sample",false);
            return cache.emplace(key,value).first->second;
        };
        const auto density=[](const Sample& p){return static_cast<long double>(p.speed)/p.size;};
        const double slope=std::sqrt(1.8)-1;
        long double estimatedError=0;
        while(!pending.empty()){
            tick(1);const auto interval=pending.back();pending.pop_back();
            const auto points=stencil(interval.first,interval.last);
            const auto [a,q1,m,q3,b]=points;
            if(!points.resolved()){
                report.numericalInterval=std::array{a,b};fail("curve_parameter_resolution",false);
            }
            std::array<Sample,5> stencil{sample(a,interval.span),sample(q1,interval.span),
                sample(m,interval.span),sample(q3,interval.span),sample(b,interval.span)};
            const std::array rawSizes{stencil[0].size,stencil[1].size,stencil[2].size,stencil[3].size,stencil[4].size};
            // Prepare the growth-limited metric before estimating its error.
            // Pointwise curvature/feature changes must not force endless
            // relative-error subdivision. These local passes only lower sizes
            // using chord distances (a conservative arclength bound); raw
            // cached requests and the later global growth pass stay separate.
            // The estimate describes this prepared metric, not a certificate
            // of integration error in a discontinuous raw curvature field.
            for(std::size_t i=1;i<stencil.size();++i){tick(1);stencil[i].size=std::min(stencil[i].size,
                stencil[i-1].size+slope*(stencil[i].point-stencil[i-1].point).norm());}
            for(std::size_t i=stencil.size()-1;i>0;--i){tick(1);stencil[i-1].size=std::min(stencil[i-1].size,
                stencil[i].size+slope*(stencil[i].point-stencil[i-1].point).norm());}
            const auto& s0=stencil[0];const auto& s1=stencil[1];const auto& s2=stencil[2];
            const auto& s3=stencil[3];const auto& s4=stencil[4];
            const auto width=static_cast<long double>(b-a);
            const auto coarse=width*(density(s0)+density(s4))*.5L;
            const auto fine=width*(density(s0)+2*density(s1)+2*density(s2)+2*density(s3)+density(s4))*.125L;
            const auto error=std::abs(fine-coarse)/15;
            const auto coarseLength=width*(s0.speed+s4.speed)*.5L;
            const auto fineLength=width*(s0.speed+2*s1.speed+2*s2.speed+2*s3.speed+s4.speed)*.125L;
            const auto lengthError=std::abs(fineLength-coarseLength)/15;
            const auto h=std::min({s0.size,s1.size,s2.size,s3.size,s4.size});
            const auto mappingError=std::max({(s1.point-(s0.point*.75+s4.point*.25)).norm(),
                (s2.point-(s0.point+s4.point)*.5).norm(),(s3.point-(s0.point*.25+s4.point*.75)).norm()});
            const auto physicalTolerance=spacing.maximumChordError>0?std::min(h*.01,spacing.maximumChordError*.25):h*.01;
            const bool accurate=std::isfinite(fine)&&std::isfinite(fineLength)&&fine<=1&&
                error<=1e-4L*std::abs(fine)+1e-10L*((b-a)/span)&&
                lengthError<=1e-5L*std::abs(fineLength)+1e-12L*h*((b-a)/span)&&
                std::isfinite(mappingError)&&mappingError<=physicalTolerance;
            if(accurate){++report.integrationIntervals;estimatedError+=error;continue;}
            report.lastRejectedInterval=std::array{a,b};
            report.lastRejectedSpeeds=std::array{s0.speed,s1.speed,s2.speed,s3.speed,s4.speed};
            report.lastRejectedSizes=std::array{s0.size,s1.size,s2.size,s3.size,s4.size};
            report.lastRejectedRawSizes=rawSizes;
            report.lastRejectedErrors=std::array{double(error),double(1e-4L*std::abs(fine)+1e-10L*((b-a)/span)),
                double(lengthError),double(1e-5L*std::abs(fineLength)+1e-12L*h*((b-a)/span)),mappingError,physicalTolerance};
            if(report.integrationIntervals+pending.size()+2>intervals)fail("curve_integration_limit",true);
            admit(cache.size(),std::max<std::uint64_t>(2*pending.capacity(),2*(pending.size()+2)),breaks.capacity(),0,0);
            pending.push_back({m,b,interval.span});pending.push_back({a,m,interval.span});
        }
        admit(cache.size(),2*pending.capacity(),breaks.capacity(),cache.size(),0);
        std::vector<Sample> samples;samples.reserve(cache.size());
        for(const auto& [at,value]:cache){(void)at;tick(1);samples.push_back(value);}
        std::vector<double> length(samples.size()-1),metric(samples.size(),0.);
        for(std::size_t i=0;i<length.size();++i){
            tick(1);length[i]=(samples[i+1].at-samples[i].at)*(samples[i].speed+samples[i+1].speed)*.5;
            if(!std::isfinite(length[i])||length[i]<0)fail("curve_nonfinite_length",false);
        }
        for(std::size_t i=0;i<length.size();++i){tick(1);samples[i+1].size=std::min(samples[i+1].size,samples[i].size+slope*length[i]);}
        for(std::size_t i=length.size();i>0;--i){tick(1);samples[i-1].size=std::min(samples[i-1].size,samples[i].size+slope*length[i-1]);}
        const auto logarithm=[](double a,double b){
            const auto ratio=(b-a)/a;
            return std::isfinite(ratio)&&ratio>-1?std::log1p(ratio):std::log(b)-std::log(a);
        };
        for(std::size_t i=0;i<length.size();++i){
            tick(1);const auto a=samples[i].size,b=samples[i+1].size;
            const auto amount=a==b?length[i]/a:length[i]*logarithm(a,b)/(b-a);
            metric[i+1]=metric[i]+amount;
            if(!std::isfinite(metric[i+1])||amount<0)fail("curve_metric_overflow",true);
        }
        report.metricLength=metric.back();report.estimatedMetricError=double(estimatedError);
        if(!(metric.back()>0))fail("curve_zero_metric",false);
        const auto minimum=spacing.maximumChordError>0&&closed?3u:1u;
        const auto count=std::max<std::uint64_t>(minimum,std::uint64_t(std::ceil(std::min(metric.back(),double(segments)+1))));
        if(count>segments)fail("curve_segment_limit",true);
        admit(cache.size(),2*pending.capacity(),breaks.capacity(),samples.size(),count+1);
        std::vector<double> schedule;schedule.reserve(std::size_t(count+1));schedule.push_back(first);
        for(std::uint64_t i=1;i<count;++i){
            tick(1);const auto target=metric.back()*double(i)/count;
            const auto at=std::size_t(std::lower_bound(metric.begin(),metric.end(),target)-metric.begin());
            if(!at||at>=metric.size()||!(metric[at]>metric[at-1]))fail("curve_metric_resolution",false);
            const auto fraction=std::clamp((target-metric[at-1])/(metric[at]-metric[at-1]),0.,1.);
            const auto a=samples[at-1].size,b=samples[at].size;
            double physicalFraction=fraction;
            if(a!=b){
                const auto ratioLog=logarithm(a,b);
                const auto change=std::abs(ratioLog)<.1?a*std::expm1(fraction*ratioLog):
                    std::exp(std::lerp(std::log(a),std::log(b),fraction))-a;
                physicalFraction=std::clamp(change/(b-a),0.,1.);
            }
            const auto parameter=std::lerp(samples[at-1].at,samples[at].at,physicalFraction);
            if(!(parameter>schedule.back()&&parameter<last))fail("curve_parameter_resolution",false);
            schedule.push_back(parameter);
        }
        schedule.push_back(last);report.segments=count;report.completed=true;report.diagnostic="completed";
        return schedule;
    }catch(const SamplingStop&){return {};}
}
} // namespace curve_metric_detail

template<class Tick>
std::vector<double> sampleOcctCurveMetric(const OcctModel& model, int tag, double sizeCap,
    const OcctCurveUses& uses, const OcctVertexEdges& incident, Tick tick,
    OcctCurveSpacing spacing = {},const OcctSpatialSizing* sizeField = nullptr) {
    constexpr std::size_t intervals = 256;
    const auto localCurvature = [&](int edge, double t,const curve_metric_detail::NativeSpan* local=nullptr) {
        tick(1);
        double curvature = 0;
        const auto include = [&](double value) {
            if (std::isfinite(value)) curvature = std::max(curvature, value);
        };
        if(local){
            gp_Pnt point;gp_Vec derivative,second;
            try{
                local->D2(t,point,derivative,second);
                const auto magnitude=derivative.Magnitude();
                if(magnitude>0.)include(derivative.Crossed(second).Magnitude()/(magnitude*magnitude*magnitude));
            }catch(const Standard_OutOfMemory&){throw std::bad_alloc{};}
             catch(const Standard_Failure&){/* Match Edge::curvatureAt's unavailable-curvature convention. */}
        }else include(model.edge(edge)->curvatureAt(t));
        for (const auto& [face, use] : uses.at(edge)) {
            tick(2);
            Vec2 uv;
            if (model.parameterOnUse(use, t, uv)) include(model.face(face)->curvatureAt(uv));
        }
        return curvature;
    };
    const auto vertexCurvature = [&](int vertex) {
        double curvature = 0;
        for (int edgeTag : incident.at(vertex)) {
            tick(1);
            const auto* edge = model.edge(edgeTag);
            double first, last;
            edge->parameterBounds(first, last);
            if (edge->startVertexTag() == vertex)
                curvature = std::max(curvature, localCurvature(edgeTag, first));
            if (edge->endVertexTag() == vertex)
                curvature = std::max(curvature, localCurvature(edgeTag, last));
        }
        return curvature;
    };
    const auto* edge = model.edge(tag);
    double first, last;
    edge->parameterBounds(first, last);
    BRepAdaptor_Curve curve(OcctModelAccess::edge(model, tag));
    if(spacing.resourcePolicy==CadResourcePolicy::Normal){
        const auto refuse=[&](const char* diagnostic,bool budget){
            if(spacing.report){*spacing.report={};spacing.report->diagnostic=diagnostic;
                spacing.report->budgetLimited=budget;}
            return std::vector<double>{};
        };
        if(!spacing.maximumSegments||!spacing.maximumIntegrationIntervals||spacing.maximumWorkingBytes<4096)
            return refuse("curve_resource_limit",true);
        if(!std::isfinite(first)||!std::isfinite(last)||!(last>first)||!std::isfinite(last-first))
            return refuse("curve_parameter_range",false);
        const auto requestedSize=[&](double curvature,Vec3 position){
            auto h=curvature>0?std::min(sizeCap,2*std::numbers::pi/(spacing.pointsPerTurn*curvature)):sizeCap;
            if(curvature>0&&spacing.maximumChordError>0)
                h=std::min(h,std::sqrt(8*spacing.maximumChordError/curvature));
            if(sizeField){tick(sizeField->queryWork());h=std::min(h,sizeField->eval(position));}
            return h;
        };
        tick(1);const auto start=curve.Value(first);tick(1);const auto end=curve.Value(last);
        const Vec3 startPoint{start.X(),start.Y(),start.Z()},endPoint{end.X(),end.Y(),end.Z()};
        auto startSize=requestedSize(vertexCurvature(edge->startVertexTag()),startPoint);
        auto endSize=requestedSize(vertexCurvature(edge->endVertexTag()),endPoint);
        const bool closed=edge->startVertexTag()==edge->endVertexTag();
        if(closed)startSize=endSize=std::min(startSize,endSize);
        const auto finitePoint=[](Vec3 p){return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);};
        if(!finitePoint(startPoint)||!finitePoint(endPoint)||!(std::isfinite(startSize)&&startSize>0&&
            std::isfinite(endSize)&&endSize>0))return refuse("curve_nonfinite_sample",false);
        // Endpoint feature sizes must influence a continuous neighbourhood
        // before quadrature. An isolated endpoint override makes relative
        // quadrature error persist through every subdivision. Euclidean
        // distance is a lower bound on physical arclength, so this existing
        // growth-slope envelope is conservative; final arclength smoothing
        // and all native geometric audits still apply.
        const double slope=std::sqrt(1.8)-1;
        const auto evaluate=[&](double at,const curve_metric_detail::NativeSpan& local){
            gp_Pnt point;gp_Vec derivative;local.D1(at,point,derivative);
            const Vec3 position{point.X(),point.Y(),point.Z()};
            const auto h=std::min({requestedSize(localCurvature(tag,at,&local),position),
                startSize+slope*(position-startPoint).norm(),endSize+slope*(position-endPoint).norm()});
            return curve_metric_detail::Sample{at,derivative.Magnitude(),h,position};
        };
        return curve_metric_detail::progressive(curve,first,last,closed,
            spacing,evaluate,tick);
    }
    std::array<double, intervals+1> t{}, speed{}, size{}, metric{};
    std::array<double, intervals> length{};
    for (std::size_t i = 0; i <= intervals; ++i) {
        tick(1);
        t[i] = i == intervals ? last : std::lerp(first, last, double(i)/intervals);
        gp_Pnt point;
        gp_Vec derivative;
        curve.D1(t[i], point, derivative);
        speed[i] = derivative.Magnitude();
        auto curvature = localCurvature(tag, t[i]);
        if (i == 0) curvature = std::max(curvature, vertexCurvature(edge->startVertexTag()));
        if (i == intervals) curvature = std::max(curvature, vertexCurvature(edge->endVertexTag()));
        size[i] = curvature > 0 ?
            std::min(sizeCap, 2*std::numbers::pi/(spacing.pointsPerTurn*curvature)) : sizeCap;
        if (curvature > 0 && spacing.maximumChordError > 0)
            size[i] = std::min(size[i], std::sqrt(8*spacing.maximumChordError/curvature));
        if (sizeField) {
            tick(sizeField->queryWork());
            size[i] = std::min(size[i], sizeField->eval({point.X(), point.Y(), point.Z()}));
        }
        if (!(std::isfinite(speed[i]) && std::isfinite(size[i]) && size[i] > 0)) return {};
        if (i) length[i-1] = (t[i]-t[i-1])*(speed[i]+speed[i-1])*.5;
    }
    const double slope = std::sqrt(1.8)-1;
    for (std::size_t i = 0; i < intervals; ++i)
        size[i+1] = std::min(size[i+1], size[i]+slope*length[i]);
    for (std::size_t i = intervals; i > 0; --i)
        size[i-1] = std::min(size[i-1], size[i]+slope*length[i-1]);
    for (std::size_t i = 1; i <= intervals; ++i)
        metric[i] = metric[i-1]+(t[i]-t[i-1])*(speed[i]/size[i]+speed[i-1]/size[i-1])*.5;
    if (!(std::isfinite(metric.back()) && metric.back() > 0 && metric.back() <= 256)) return {};
    const auto minimum = spacing.maximumChordError > 0 && edge->startVertexTag()==edge->endVertexTag()?3u:1u;
    const auto count = std::max(minimum, static_cast<unsigned>(std::ceil(metric.back()-.01)));
    std::vector<double> schedule{first};
    for (unsigned i = 1; i < count; ++i) {
        const auto target = metric.back()*i/count;
        const auto j = std::size_t(std::lower_bound(metric.begin(), metric.end(), target)-metric.begin());
        const auto parameter = std::lerp(t[j-1], t[j], (target-metric[j-1])/(metric[j]-metric[j-1]));
        if (!(parameter > schedule.back() && parameter < last)) return {};
        schedule.push_back(parameter);
    }
    schedule.push_back(last);
    return schedule;
}

struct OcctChartMetric {
    double u, shear, v;

    Vec2 physical(Vec2 uv) const { return {u*uv.u+shear*uv.v,v*uv.v}; }
    Vec2 parameter(Vec2 p) const {
        const double y=p.v/v;
        return {(p.u-shear*y)/u,y};
    }
};

inline std::optional<OcctChartMetric> occtChartMetric(double a,double b,double c) {
    if (!(std::isfinite(a)&&std::isfinite(b)&&std::isfinite(c)&&a>0&&c>0)) return {};
    const auto u=std::sqrt(a),shear=b/u,v=std::sqrt(c-shear*shear);
    if (!(std::isfinite(v)&&v>0)) return {};
    return OcctChartMetric{u,shear,v};
}

inline std::optional<OcctChartMetric> occtChartMetricAt(const BRepAdaptor_Surface& surface,Vec2 uv) {
    gp_Pnt point;
    gp_Vec u,v;
    surface.D1(uv.u,uv.v,point,u,v);
    return occtChartMetric(u.SquareMagnitude(),u.Dot(v),v.SquareMagnitude());
}

template<class Tick>
bool metricOcctChart(const OcctModel& model, int face, std::span<const Vec2> samples,
                     Tick tick, std::vector<Vec2>& coordinates) {
    BRepAdaptor_Surface surface(OcctModelAccess::face(model, face));
    double a = 0, b = 0, c = 0;
    for (const auto uv : samples) {
        tick(1);
        gp_Pnt point;
        gp_Vec u, v;
        surface.D1(uv.u, uv.v, point, u, v);
        a += u.SquareMagnitude();
        b += u.Dot(v);
        c += v.SquareMagnitude();
    }
    if (samples.empty() || !(std::isfinite(a) && std::isfinite(b) && std::isfinite(c) &&
                             a > 0 && c > 0 && a*c > b*b)) return false;
    const auto scaleU=std::sqrt(a/samples.size());
    const auto shearV=b/(samples.size()*scaleU);
    const auto scaleV=std::sqrt(c/samples.size()-shearV*shearV);
    if (!(std::isfinite(scaleV)&&scaleV>0)) return false;
    const OcctChartMetric metric{scaleU,shearV,scaleV};
    coordinates.reserve(samples.size());
    for (const auto uv : samples) {
        const Vec2 transformed=metric.physical(uv);
        if (!(std::isfinite(transformed.u) && std::isfinite(transformed.v))) return false;
        coordinates.push_back(transformed);
    }
    return true;
}

template<class Triangles, class Tick>
std::optional<std::vector<int>> labelOcctChartRegions(const Triangles& triangles,
    const std::set<std::array<Index,2>>& constraints, Tick tick) {
    using Edge = std::array<Index,2>;
    std::map<Edge, std::vector<std::size_t>> owners;
    for (std::size_t i = 0; i < triangles.size(); ++i) {
        tick(3);
        for (unsigned side = 0; side < 3; ++side) {
            const auto [lo,hi] = std::minmax(triangles[i].v[side], triangles[i].v[(side+1)%3]);
            owners[{lo,hi}].push_back(i);
        }
    }
    for (const auto& edge : constraints) {
        tick(1);
        if (!owners.contains(edge)) return {};
    }
    std::vector<std::vector<std::pair<std::size_t, int>>> neighbours(triangles.size());
    std::vector<int> parity(triangles.size(), -1);
    std::vector<std::size_t> pending;
    const auto assign = [&](std::size_t i, int label) {
        if (parity[i] < 0) {
            parity[i] = label;
            pending.push_back(i);
        }
        return parity[i] == label;
    };
    for (const auto& [edge, uses] : owners) {
        tick(1);
        const int crosses = constraints.contains(edge);
        if (uses.size() == 1) {
            if (!assign(uses.front(), crosses)) return {};
        } else if (uses.size() == 2) {
            neighbours[uses[0]].push_back({uses[1], crosses});
            neighbours[uses[1]].push_back({uses[0], crosses});
        } else return {};
    }
    for (std::size_t next = 0; next < pending.size(); ++next)
        for (const auto& [other, crosses] : neighbours[pending[next]]) {
            tick(1);
            if (!assign(other, parity[pending[next]] ^ crosses)) return {};
        }
    if (pending.size() != triangles.size()) return {};
    return parity;
}

} // namespace paacs::geo::detail
