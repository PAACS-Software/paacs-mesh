// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/geometry/OcctLocalDelaunay.hpp"
#include "../src/geometry/OcctNativeValidation.hpp"
#include "../src/geometry/OcctModelInternal.hpp"
#include "../src/geometry/OcctSolidClassifier.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <gp_Trsf.hxx>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace {std::atomic_bool failNextAllocation{false};}
void* operator new(std::size_t n){
    if(failNextAllocation.exchange(false,std::memory_order_relaxed))throw std::bad_alloc{};
    if(void* p=std::malloc(n?n:1))return p;throw std::bad_alloc{};
}
void* operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void* p)noexcept{std::free(p);}
void operator delete[](void* p)noexcept{std::free(p);}
void operator delete(void* p,std::size_t)noexcept{std::free(p);}
void operator delete[](void* p,std::size_t)noexcept{std::free(p);}
using namespace paacs;
using namespace paacs::geo::detail;
namespace {
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}

void exerciseLocalDelaunay(){
    mesh::Mesh source;
    const std::array vertices{Vec3{0,0,0},Vec3{2,0,0},Vec3{0,2,0},Vec3{0,0,2},Vec3{0,0,-2}};
    for(std::size_t i=0;i<vertices.size();++i)source.addNode(vertices[i],{0,int(i)+1});
    constexpr mesh::Classification region{3,11};
    const std::array<bnd::TetFacetCell,2> cells{{{{0,1,2,3},0,region},{{0,2,1,4},0,region}}};
    for(const auto& cell:cells)for(unsigned side=0;side<4;++side){
        const auto face=bnd::tet_facet_detail::outwardFace(cell,side);
        if(*std::max_element(face.begin(),face.end())<3)continue;
        source.addCell(mesh::CellType::Tri3,{face[0],face[1],face[2]},{2,int(source.cellCount())+1});
    }
    const Index boundaryCells=source.cellCount();
    for(const auto& cell:cells){
        const auto& v=cell.vertices;
        source.addCell(mesh::CellType::Tet4,{v[0],v[1],v[2],v[3]},region);
    }
    check(mesh::auditLinearVolumeBoundary(source).valid,"source has a closed oriented shell");

    const auto samePrefix=[&](const mesh::Mesh& candidate,Index cellCount){
        for(Index i=0;i<source.nodeCount();++i){
            const auto a=candidate.node(i),b=source.node(i);
            check(a.x==b.x&&a.y==b.y&&a.z==b.z&&candidate.nodeClass(i)==source.nodeClass(i),
                  "original coordinates and node tags unchanged");
        }
        for(Index i=0;i<cellCount;++i){
            const auto a=candidate.cell(i),b=source.cell(i);
            const auto* v=source.cellNodes(i);
            check(a.type==b.type&&a.onEntity==b.onEntity&&
                  std::equal(v,v+(b.type==mesh::CellType::Tri3?3:4),candidate.cellNodes(i)),
                  "original cell connectivity, winding and CAD tag unchanged");
        }
    };
    const auto material=[](Vec3 p){return p.x>0&&p.y>0&&p.x+p.y+std::abs(p.z)<2;};
    const auto charge=[](std::uint64_t){};
    const LocalDelaunayLimits limits{100'000,64,16,16,1024*1024,1.};
    const std::array sites{Vec3{.25,.25,.25}};
    auto refined=source;
    check(insertOcctDelaunaySites(refined,sites,limits,material,material,charge)==1,
          "admitted support site inserts into the live mesh");
    check(refined.nodeCount()==6&&refined.nodeClass(5)==region&&
          refined.countCells(mesh::CellType::Tet4)==6&&refined.countCells(mesh::CellType::Tri3)==boundaryCells,
          "in-sphere cavity replaces both adjacent cells and retains the shell");
    samePrefix(refined,boundaryCells);
    check(mesh::auditLinearVolumeBoundary(refined).valid,"replacement volume retains a conforming shell");
    for(Index i=boundaryCells;i<refined.cellCount();++i){
        const auto* v=refined.cellNodes(i);
        const auto a=refined.node(v[0]),b=refined.node(v[1]),c=refined.node(v[2]),d=refined.node(v[3]);
        check(predicates::orient3d(a,b,c,d)==predicates::Sign::Positive&&
              material((a+b+c+d)*.25)&&refined.cell(i).onEntity==region,
              "replacement tetrahedron has positive orientation, material and body tag");
    }

    auto refused=source;
    unsigned centroidChecks=0;
    const auto refuseCentroid=[&](Vec3){++centroidChecks;return false;};
    check(insertOcctDelaunaySites(refused,sites,limits,material,refuseCentroid,charge)==0&&centroidChecks>0,
          "CAD centroid rejection refuses the staged cavity");
    check(refused.nodeCount()==source.nodeCount()&&refused.cellCount()==source.cellCount(),
          "CAD refusal leaves all original cells and nodes");
    samePrefix(refused,source.cellCount());

    auto shortBudget=limits;
    shortBudget.queries=1;
    unsigned nodeChecks=0;
    const auto countNode=[&](Vec3 p){++nodeChecks;return material(p);};
    check(insertOcctDelaunaySites(refused,sites,shortBudget,countNode,material,charge)==0&&nodeChecks==1,
          "query exhaustion after staging leaves the cavity uncommitted");
    check(refused.nodeCount()==source.nodeCount()&&refused.cellCount()==source.cellCount(),
          "budget refusal leaves all original cells and nodes");
    samePrefix(refused,source.cellCount());
}

std::string snapshot(const geo::OcctModel& model){std::ostringstream s;BRepTools::Write(OcctModelAccess::shape(model),s);return s.str();}
std::vector<Vec3> points(const gp_Trsf& transform){
    std::vector<Vec3> result;
    for(unsigned i=0;i<128;++i){const auto angle=2*std::numbers::pi*i/128.;
        gp_Pnt p(1.5*std::cos(angle),1.5*std::sin(angle),.1+2.8*(i%17)/16.);p.Transform(transform);result.push_back({p.X(),p.Y(),p.Z()});}
    for(double z:{0.,3.}){gp_Pnt p(1.5,0,z);p.Transform(transform);result.push_back({p.X(),p.Y(),p.Z()});}
    return result;
}
void exercise(const TopoDS_Shape& shape,const gp_Trsf& transform){
    auto imported=OcctModelAccess::fromShape(shape);check(imported.imported(),"valid tube native snapshot");
    const auto& model=*imported.model;const auto before=snapshot(model);auto input=points(transform);
    for(unsigned threads:{1u,2u,4u,8u}){
        OcctNativeValidationOptions options;options.maximumThreads=threads;
        const auto result=validateOcctNativeInteriorPoints(model,1,input,options);
        if(!result.completed)std::cerr<<"worker reason "<<occtNativeValidationReasonName(result.reason)<<" threads="<<threads<<'\n';
        check(result.completed&&result.reason==OcctNativeValidationReason::None&&result.sourceModelUnmodified,"all actual tube material points accepted");
        check(result.nativeQueries==input.size()&&result.points==input.size()&&result.threadsUsed==threads,"every point queried exactly once");
        check(result.workUsed<=options.maximumWork&&result.maximumWorkingPayloadEstimate<=options.maximumWorkingBytes,"native query admission measured");
        if(threads>1)check(result.originalReferenceQueries==input.size()&&result.copiedGeometryValuePairs>0,"entire small fixture independently compared with original classifier");
        check(before==snapshot(model),"native source bytes unchanged after copied query workers");
    }
    // A faceted concave CAD boundary can place a cell centroid slightly in
    // the native void. The explicit distance policy records that allowance;
    // strict inclusion and a genuine hole remain separate cases.
    gp_Pnt nearHole(.99,0,1.5);nearHole.Transform(transform);input[11]={nearHole.X(),nearHole.Y(),nearHole.Z()};
    OcctNativeValidationOptions strict;strict.maximumThreads=4;
    auto strictResult=validateOcctNativeInteriorPoints(model,1,input,strict);
    check(!strictResult.completed&&strictResult.reason==OcctNativeValidationReason::OutsideMaterial,"default validation still refuses near-boundary OUT point");
    std::vector<unsigned char> mask;
    const auto classified=classifyOcctNativeInteriorPoints(model,1,input,strict,mask);
    check(classified.completed&&classified.sourceModelUnmodified&&classified.firstFailedPoint==11&&
          mask.size()==input.size()&&mask[11]==0&&mask[12]==1,"batch classifies through an outside point");
    const auto savedMask=mask;
    const OcctNativeValidationOptions serial;
    auto rayOptions=serial;
    rayOptions.useTrimmedRayWitnesses=true;
    std::vector<Vec3> rayInput;
    for(unsigned repeat=0;repeat<4;++repeat)rayInput.insert(rayInput.end(),input.begin(),input.end());
    std::vector<unsigned char> rayMask;
    OcctNativeStrictAdmissions rayReceipt(4ULL*1024*1024);
    const auto rayResult=classifyOcctNativeInteriorPoints(model,1,rayInput,rayOptions,rayMask,&rayReceipt);
    check(rayResult.completed&&rayResult.sourceModelUnmodified&&rayMask.size()==rayInput.size()&&
          rayResult.nativeQueries+rayResult.trimmedRayAdmissions==rayInput.size(),
          "opt-in ray witnesses account for every original-solid query");
    for(std::size_t i=0;i<rayMask.size();++i)
        check(rayMask[i]==savedMask[i%savedMask.size()],"ray witnesses preserve the tube material mask");
    if(rayResult.trimmedRayAdmissions){
        const std::array point{input[12]};
        const auto ordinary=validateOcctNativeInteriorPoints(model,1,point,serial,&rayReceipt);
        check(ordinary.completed&&ordinary.nativeQueries==1&&ordinary.strictAdmissionReuses==0,
              "the default classifier cannot reuse an opt-in ray receipt");
    }
    auto tinyBatch=strict;tinyBatch.maximumWork=1;
    const auto stoppedBatch=classifyOcctNativeInteriorPoints(model,1,input,tinyBatch,mask);
    check(!stoppedBatch.completed&&stoppedBatch.reason==OcctNativeValidationReason::BudgetExceeded&&mask==savedMask,
          "failed batch preserves the destination mask");
    OcctNativeStrictAdmissions prior(4ULL*1024*1024);
    const auto captured=classifyOcctNativeInteriorPoints(model,1,input,serial,mask,&prior);
    check(captured.completed&&captured.sourceModelUnmodified&&mask==savedMask&&prior.retainedPayloadBytes()>0,
          "only a complete audited serial classification retains strict admissions");
    gp_Pnt freshPoint(1.7,0,1.25);freshPoint.Transform(transform);
    const std::array finalPoints{input[12],input[4],Vec3{freshPoint.X(),freshPoint.Y(),freshPoint.Z()}};
    const auto reused=validateOcctNativeInteriorPoints(model,1,finalPoints,serial,&prior);
    check(reused.completed&&reused.sourceModelUnmodified&&reused.strictAdmissionReuses==2&&reused.nativeQueries==1,
          "reordered exact admissions reuse prior checks and unseen coordinates query OCCT");
    const std::array outsidePoints{input[12],input[11]};
    OcctNativeStrictAdmissions updated(4ULL*1024*1024);
    const auto remasked=classifyOcctNativeInteriorPoints(model,1,outsidePoints,serial,mask,&updated,&prior);
    check(remasked.completed&&mask==std::vector<unsigned char>{1,0}&&
          remasked.nativeQueries==1&&remasked.strictAdmissionReuses==1,
          "masked replay keeps exact positive rows and queries previously refused rows");
    const auto refreshed=classifyOcctNativeInteriorPoints(model,1,finalPoints,serial,mask,&prior,&prior);
    check(refreshed.completed&&mask==std::vector<unsigned char>{1,1,1}&&
          refreshed.nativeQueries==1&&refreshed.strictAdmissionReuses==2,
          "one receipt can be refreshed from its earlier positive rows");
    const auto afterRefresh=validateOcctNativeInteriorPoints(model,1,finalPoints,serial,&prior);
    check(afterRefresh.completed&&afterRefresh.nativeQueries==0&&afterRefresh.strictAdmissionReuses==3,
          "the refreshed receipt keeps reused and newly queried admissions");
    const auto strictOutside=validateOcctNativeInteriorPoints(model,1,outsidePoints,serial,&prior);
    check(!strictOutside.completed&&strictOutside.reason==OcctNativeValidationReason::OutsideMaterial&&
          strictOutside.firstFailedPoint==1&&strictOutside.strictAdmissionReuses==1&&strictOutside.nativeQueries==1,
          "earlier outside result does not become a reusable strict admission");
    auto otherModel=OcctModelAccess::fromShape(shape);
    check(otherModel.imported(),"second native model identity fixture");
    const auto uncached=validateOcctNativeInteriorPoints(*otherModel.model,1,finalPoints,serial,&prior);
    check(uncached.completed&&uncached.nativeQueries==finalPoints.size()&&uncached.strictAdmissionReuses==0,
          "strict admissions cannot be shared across model identities");
    auto wider=serial;wider.maximumOutsideDistance=.02;
    const auto distanceBand=validateOcctNativeInteriorPoints(model,1,outsidePoints,wider,&prior);
    check(distanceBand.completed&&distanceBand.nativeQueries==outsidePoints.size()&&distanceBand.strictAdmissionReuses==0&&
          distanceBand.pointsWithinDistanceBand==1,"non-strict policy does its ordinary native queries");
    for(unsigned threads:{1u,4u}){
        OcctNativeValidationOptions band;band.maximumThreads=threads;band.maximumOutsideDistance=.02;
        const auto result=validateOcctNativeInteriorPoints(model,1,input,band);
        check(result.completed&&result.sourceModelUnmodified&&result.pointsWithinDistanceBand==1&&result.nativeDistanceQueries==1,"explicit native distance band records one admitted OUT point");
        check(std::abs(result.maximumOutsideDistance-.01)<1e-6,"actual distance to the original tube measured");
        band.maximumOutsideDistance=.005;const auto refused=validateOcctNativeInteriorPoints(model,1,input,band);
        check(!refused.completed&&refused.reason==OcctNativeValidationReason::OutsideMaterial&&refused.firstFailedPoint==11,"a point beyond the requested distance still refuses");
    }
    gp_Pnt hole(0,0,1.5);hole.Transform(transform);input[11]={hole.X(),hole.Y(),hole.Z()};
    OcctNativeValidationOptions options;options.maximumThreads=4;
    auto result=validateOcctNativeInteriorPoints(model,1,input,options);
    check(!result.completed&&result.reason==OcctNativeValidationReason::OutsideMaterial&&result.firstFailedPoint==11,"a point in the genuine tube hole remains outside");
    check(result.nativeQueries==input.size(),"parallel refusal still accounts for every performed query");
    options.maximumOutsideDistance=.02;result=validateOcctNativeInteriorPoints(model,1,input,options);
    check(!result.completed&&result.reason==OcctNativeValidationReason::OutsideMaterial,"distance band cannot fill the genuine tube hole");
    options.maximumOutsideDistance=-1;result=validateOcctNativeInteriorPoints(model,1,input,options);
    check(!result.completed&&result.reason==OcctNativeValidationReason::InvalidRequest,"negative native distance allowance refused");
    options.maximumOutsideDistance=std::numeric_limits<double>::quiet_NaN();result=validateOcctNativeInteriorPoints(model,1,input,options);
    check(!result.completed&&result.reason==OcctNativeValidationReason::InvalidRequest,"nonfinite native distance allowance refused");
    input=points(transform);
    options=OcctNativeValidationOptions{};options.maximumThreads=4;options.maximumWork=1;result=validateOcctNativeInteriorPoints(model,1,input,options);
    check(!result.completed&&result.reason==OcctNativeValidationReason::BudgetExceeded&&result.nativeQueries==0,"insufficient query work refuses before launch");
    options=OcctNativeValidationOptions{};options.maximumThreads=4;options.maximumSnapshotBytes=1;
    result=validateOcctNativeInteriorPoints(model,1,input,options);
    check(!result.completed&&result.reason==OcctNativeValidationReason::BudgetExceeded,"bounded native BREP snapshot refuses tiny cap");
    std::atomic_bool cancel{true};options=OcctNativeValidationOptions{};options.cancellation=&cancel;
    result=validateOcctNativeInteriorPoints(model,1,input,options);
    check(!result.completed&&result.reason==OcctNativeValidationReason::Cancelled&&result.nativeQueries==0,"cancelled native query batch refuses before native work");
    options=OcctNativeValidationOptions{};options.maximumThreads=9;
    result=validateOcctNativeInteriorPoints(model,1,input,options);
    check(!result.completed&&result.reason==OcctNativeValidationReason::InvalidRequest,"worker hard cap checked");
    options.maximumThreads=4;input[4].x=std::numeric_limits<double>::quiet_NaN();
    result=validateOcctNativeInteriorPoints(model,1,input,options);
    check(!result.completed&&result.reason==OcctNativeValidationReason::InvalidRequest&&result.firstFailedPoint==4,"nonfinite native query input refused");
    input=points(transform);failNextAllocation.store(true,std::memory_order_relaxed);
    result=validateOcctNativeInteriorPoints(model,1,input,options);failNextAllocation.store(false,std::memory_order_relaxed);
    check(!result.completed&&result.reason==OcctNativeValidationReason::AllocationFailed,"snapshot allocation failure remains distinct");
    check(before==snapshot(model),"all failure modes preserve native source bytes");
}
}
int main(int argc,char** argv){try{
    exerciseLocalDelaunay();
    const auto outer=BRepPrimAPI_MakeCylinder(2,3).Shape(),inner=BRepPrimAPI_MakeCylinder(1,3).Shape();
    BRepAlgoAPI_Cut cut(outer,inner);check(cut.IsDone(),"native hollow cylinder fixture");
    gp_Trsf identity;exercise(cut.Shape(),identity);
    gp_Trsf rotation;rotation.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(1,2,3)),.73);rotation.SetTranslationPart(gp_Vec(13,-4,7));
    BRepBuilderAPI_Transform relocated(cut.Shape(),rotation,true);check(relocated.IsDone(),"located native tube");exercise(relocated.Shape(),rotation);
    BRepBuilderAPI_NurbsConvert converted(relocated.Shape(),true);check(converted.IsDone(),"native spline tube fixture");exercise(converted.Shape(),rotation);
    if(argc>1){
        auto imported=geo::importOcctCad(argv[1],{});check(imported.imported(),"local classifier regression import");const auto& model=*imported.model;const auto before=snapshot(model);
        check(model.solids().size()==1,"local regression expected one original body");const auto tag=model.solids().front()->tag();
        OcctSolidClassifier classifier(OcctModelAccess::solid(model,tag));
        const Vec3 witness{-6.7038004590372005,13.716460113414932,29.738673362698854};
        classifier.Perform(gp_Pnt(witness.x,witness.y,witness.z),model.tolerance(3,tag));check(classifier.State()==TopAbs_IN,"known material point falsely classified outside");
        classifier.Perform(gp_Pnt(.13,20,32),model.tolerance(3,tag));check(classifier.State()==TopAbs_OUT,"real CAD gap was filled by classifier recovery");
        classifier.Perform(gp_Pnt(100,100,100),model.tolerance(3,tag));check(classifier.State()==TopAbs_OUT,"exterior point was admitted");
        check(before==snapshot(model),"local classifier recovery changed source geometry");
    }
    std::cout<<"native validation: analytic/located/spline geometry, original source guards, local Delaunay cavities and resource failures passed\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}catch(...){std::cerr<<"unexpected native validation test failure\n";return 2;}}
