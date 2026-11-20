// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/geometry/OcctQueryContext.hpp"
#include "../src/geometry/OcctModelInternal.hpp"
#include "paacs/geometry/OcctBoundary.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <gp_Trsf.hxx>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <thread>

using namespace paacs;
using namespace paacs::geo::detail;
namespace {
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
bool same(double a,double b){return std::bit_cast<std::uint64_t>(a)==std::bit_cast<std::uint64_t>(b);}
bool same(Vec2 a,Vec2 b){return same(a.u,b.u)&&same(a.v,b.v);}
bool same(Vec3 a,Vec3 b){return same(a.x,b.x)&&same(a.y,b.y)&&same(a.z,b.z);}
struct Query { const geo::GFace* face;Vec3 point; };
struct Result {bool accepted;Vec2 uv;Vec3 point;};
Result project(const Query& query){Result r{false,{91,92},{93,94,95}};r.accepted=query.face->projectPoint(query.point,r.uv,r.point);return r;}
Result trimmed(const geo::OcctModel& model,const Query& query){
    Result result{false,{91,92},{93,94,95}};
    result.accepted=geo::closestPointOnOcctFace(model,query.face->tag(),query.point,result.point);
    return result;
}
void equal(Result a,Result b){
    check(a.accepted==b.accepted&&same(a.uv,b.uv)&&same(a.point,b.point),"reused status/UV/physical point is bitwise identical to fresh query");
    if(!b.accepted)check(same(b.uv,{91,92})&&same(b.point,{93,94,95}),"refused projection keeps outputs unchanged");
}
std::string snapshot(const geo::OcctModel& model){std::ostringstream s;BRepTools::Write(OcctModelAccess::shape(model),s);return s.str();}
void exercise(const TopoDS_Shape& shape){
    auto imported=OcctModelAccess::fromShape(shape);check(imported.imported(),"valid native query fixture");
    const auto& model=*imported.model;const auto before=snapshot(model);std::vector<Query> queries;
    for(const auto* face:model.faces()){
        Vec2 lo,hi;face->paramBounds(lo,hi);
        for(unsigned i=0;i<=6;++i)for(unsigned j=0;j<=6;++j){
            Vec2 uv{lo.u+(hi.u-lo.u)*i/6.,lo.v+(hi.v-lo.v)*j/6.};const auto p=face->pointAt(uv),n=face->normalAt(uv);
            if(!std::isfinite(p.x+p.y+p.z))continue;
            queries.push_back({face,p});
            if(std::isfinite(n.x+n.y+n.z)){queries.push_back({face,p+n*.01});queries.push_back({face,p-n*.01});}
        }
        // Includes UV-chart exterior and nonfinite refusal paths.
        queries.push_back({face,face->pointAt({lo.u-(hi.u-lo.u)*.3,lo.v-(hi.v-lo.v)*.3})});
        queries.push_back({face,{std::numeric_limits<double>::quiet_NaN(),0,0}});
    }
    std::vector<Result> reference,trimmedReference;
    for(const auto& q:queries){reference.push_back(project(q));trimmedReference.push_back(trimmed(model,q));}
    check(OcctQueryContext::active()==nullptr,"default uses fresh native queries");
    {
        OcctQueryContext context(true);check(OcctQueryContext::active()==&context,"explicit current scope");
        for(std::size_t i=0;i<queries.size();++i){
            {OcctQueryContext disabled;check(OcctQueryContext::active()==nullptr,"nested disabled context uses original path");equal(reference[i],project(queries[i]));}
            check(OcctQueryContext::active()==&context,"caller context restored after nested disabled scope");
            equal(reference[i],project(queries[i]));
            equal(trimmedReference[i],trimmed(model,queries[i]));
        }
        const auto count=context.counters();
        check(count.surfaceProjectionQueries>0&&count.surfaceProjectionReuses>0&&
              count.surfaceProjectionInitializations>=(model.faces().size()>1?2u:1u),
              "real native per-face reuse and face replacement counted");
        check(count.surfaceProjectionQueries==count.surfaceProjectionInitializations+count.surfaceProjectionReuses,"each cached native query separately counted");
        check(count.trimmedProjectionQueries>0&&count.trimmedProjectionMemoHits==0,"default reuse context does not memoize results");
        {OcctQueryContext nested(true);equal(reference.front(),project(queries.front()));check(nested.counters().surfaceProjectionQueries==1,"nested contexts own separate initialized engines");}
        check(context.counters().surfaceProjectionQueries==count.surfaceProjectionQueries,"nested context cannot alter caller counters/cache");
    }
    check(OcctQueryContext::active()==nullptr,"query cache released at scope exit");
    // Use one sample across every face to expose aliases in the face key.
    const auto probe=queries.front().point;
    std::vector<Result> faceResults;
    for(const auto* face:model.faces())faceResults.push_back(trimmed(model,{face,probe}));
    {
        OcctQueryContext memo(true,true);
        for(std::size_t i=0;i<queries.size();++i){
            equal(trimmedReference[i],trimmed(model,queries[i]));
            if(trimmedReference[i].accepted){
                const auto beforeHit=memo.counters().trimmedProjectionMemoHits;
                equal(trimmedReference[i],trimmed(model,queries[i]));
                check(memo.counters().trimmedProjectionMemoHits==beforeHit+1,"immediate repetition returns the exact retained trimmed result");
            }
        }
        std::size_t index=0;
        for(const auto* face:model.faces())equal(faceResults[index++],trimmed(model,{face,probe}));
        const auto accepted=std::find_if(trimmedReference.begin(),trimmedReference.end(),[](const auto& r){return r.accepted;});
        check(accepted!=trimmedReference.end(),"fixture contains a successful trimmed projection");
        const auto& query=queries[static_cast<std::size_t>(accepted-trimmedReference.begin())];
        equal(*accepted,trimmed(model,query));
        auto opposite=OcctModelAccess::face(model,query.face->tag());
        opposite.Reverse();
        Vec3 unchanged{93,94,95};
        check(!memo.findTrimmedFacePoint(opposite,query.point,unchanged)&&same(unchanged,{93,94,95}),
              "opposite native orientation cannot reuse an answer");
        auto translated=OcctModelAccess::face(model,query.face->tag());
        gp_Trsf translation;translation.SetTranslation(gp_Vec(1,2,3));
        translated.Move(TopLoc_Location(translation));
        check(!memo.findTrimmedFacePoint(translated,query.point,unchanged)&&same(unchanged,{93,94,95}),
              "different native location cannot reuse an answer");
        const auto beforeNested=memo.counters().trimmedProjectionMemoHits;
        {
            OcctQueryContext nested(true,true);
            equal(*accepted,trimmed(model,query));
            check(nested.counters().trimmedProjectionMemoHits==0,"nested memo starts empty");
        }
        check(memo.counters().trimmedProjectionMemoHits==beforeNested,"nested context leaves caller memo counters unchanged");
        equal(*accepted,trimmed(model,query));
        check(memo.counters().trimmedProjectionMemoHits==beforeNested+1,"nested context restores caller memo");
    }
    check(OcctQueryContext::active()==nullptr,"trimmed projection memo released at scope exit");
    std::array<std::exception_ptr,4> errors;std::array<OcctQueryCounters,4> counts;
    std::array<std::jthread,4> threads;
    for(unsigned worker=0;worker<threads.size();++worker)threads[worker]=std::jthread([&,worker]{try{
        check(OcctQueryContext::active()==nullptr,"new worker starts without another thread's context");
        OcctQueryContext own(true);for(std::size_t i=0;i<queries.size();++i){
            equal(reference[i],project(queries[i]));
            equal(trimmedReference[i],trimmed(model,queries[i]));
        }counts[worker]=own.counters();
    }catch(...){errors[worker]=std::current_exception();}});
    for(auto& thread:threads)thread.join();for(const auto& error:errors)if(error)std::rethrow_exception(error);
    for(const auto& count:counts)check(count.surfaceProjectionQueries>0&&count.surfaceProjectionReuses>0,"each thread has its own bounded query context");
    check(before==snapshot(model),"fresh/reused/seam/refusal/thread queries preserve authoritative CAD bytes");
}
}
int main(){try{
    const auto outer=BRepPrimAPI_MakeCylinder(2,3).Shape(),inner=BRepPrimAPI_MakeCylinder(1,3).Shape();
    BRepAlgoAPI_Cut cut(outer,inner);check(cut.IsDone(),"real trimmed cap/hole fixture");exercise(cut.Shape());
    gp_Trsf rotation;rotation.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(1,2,3)),.73);rotation.SetTranslationPart(gp_Vec(13,-4,7));
    BRepBuilderAPI_Transform moved(cut.Shape(),rotation,true);check(moved.IsDone(),"located tube fixture");exercise(moved.Shape());
    BRepBuilderAPI_NurbsConvert spline(moved.Shape(),true);check(spline.IsDone(),"spline surface fixture");exercise(spline.Shape());
    exercise(BRepPrimAPI_MakeSphere(2).Shape());
    check(OcctQueryContext::maximumRetainedPayloadBytes()>0&&
          OcctQueryContext::maximumRetainedPayloadBytes(true)>OcctQueryContext::maximumRetainedPayloadBytes(),
          "retained logical query payload includes optional memo table");
    std::cout<<"native scoped query reuse: bitwise fresh equivalence, spline/located/seam/pole/trim refusal, nested scopes, thread isolation and source guards passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
