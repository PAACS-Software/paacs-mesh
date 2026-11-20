// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/geometry/OcctContactAssembly.hpp"
#include "../src/geometry/OcctModelInternal.hpp"
#include "paacs/geometry/OcctBoundary.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Trsf.hxx>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <new>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
std::atomic_bool failEveryAllocation{false},failNextAllocation{false};
}
void* operator new(std::size_t n){
    if(failEveryAllocation.load(std::memory_order_relaxed)||failNextAllocation.exchange(false,std::memory_order_relaxed))throw std::bad_alloc{};
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
void check(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
using Position=std::array<std::uint64_t,3>;
Position bits(Vec3 p){return {std::bit_cast<std::uint64_t>(p.x),std::bit_cast<std::uint64_t>(p.y),std::bit_cast<std::uint64_t>(p.z)};}
std::string snapshot(const geo::OcctModel& model){std::ostringstream s;BRepTools::Write(OcctModelAccess::shape(model),s);return s.str();}
TopoDS_Shape compound(const TopoDS_Shape& first,const TopoDS_Shape& second){
    BRepBuilderAPI_Copy a(first,true,false),b(second,true,false);check(a.IsDone()&&b.IsDone(),"independent copied body fixtures");
    TopoDS_Compound c;BRep_Builder builder;builder.MakeCompound(c);builder.Add(c,a.Shape());builder.Add(c,b.Shape());return c;
}
geo::NativeCadMeshOptions settings(){
    geo::NativeCadMeshOptions o;o.elements=algo::PrimitiveElements::Tet4;o.separateMatchedNativeFaceContacts=true;
    o.geometry.targetSize=.4;o.maximumChordError=.05;o.geometry.minimumElementQuality=0;o.geometry.warningElementQuality=.25;
    o.maximumValidationWork=10'000'000;o.geometry.maximumMeshingWork=10'000'000;
    o.geometry.maximumNodes=20000;o.geometry.maximumCells=100000;o.geometry.maximumOutputBytes=32ULL*1024*1024;
    o.maximumPrimitiveWorkingBytes=128ULL*1024*1024;o.geometry.maximumAuditTemporaryBytes=64ULL*1024*1024;
    return o;
}
void verify(const geo::OcctModel& model,const mesh::Mesh& mesh,const geo::NativeCadMeshOutcome& result){
    check(result.published&&result.reason==geo::NativeCadMeshReason::None&&result.audit.valid&&result.quality.valid,
          "shared native assembly publishes only valid positive audited mesh");
    check(result.nativeContactSourceGuard&&result.matchedNativeContactTriangles>0&&result.separateContactPairs==1,
          "actual shared interface and complete native source guard recorded");
    std::array<std::set<Index>,3> bodyNodes;std::array<std::set<int>,3> retainedVertices,retainedFaces;
    for(Index i=0;i<mesh.cellCount();++i){
        const auto& c=mesh.cell(i);const auto* ids=mesh.cellNodes(i);
        if(c.type==mesh::CellType::Tet4){
            check(c.onEntity.dim==3&&c.onEntity.tag>=1&&c.onEntity.tag<=2,"original native body class preserved");
            const int body=c.onEntity.tag;Vec3 centroid{};
            for(unsigned j=0;j<4;++j){bodyNodes[body].insert(ids[j]);centroid=centroid+mesh.node(ids[j]);}
            check(model.solid(body)->isInside(centroid*.25),"every actual Tet4 centroid inside original material");
        } else {
            check(c.type==mesh::CellType::Tri3&&c.onEntity.dim==2,"original native Tri3 boundary class");
            const auto owners=model.faceOwners(c.onEntity.tag);check(owners.size()==1,"fixture source faces owned by exactly one independent body");
            retainedFaces[owners.front().solidTag].insert(c.onEntity.tag);
            const auto a=mesh.node(ids[0]),b=mesh.node(ids[1]),d=mesh.node(ids[2]);
            for(const auto p:std::array{(a+b+d)/3.,(a+b)*.5,(b+d)*.5,(d+a)*.5}){
                Vec3 q;check(geo::closestPointOnOcctFace(model,c.onEntity.tag,p,q),"native trimmed face projection");
                check((p-q).norm()<=settings().maximumChordError+model.tolerance(2,c.onEntity.tag),"all actual original CAD facet samples within chord bound");
            }
        }
    }
    std::set<Position> firstPositions;std::uint64_t duplicates=0;
    for(Index id:bodyNodes[1])firstPositions.insert(bits(mesh.node(id)));
    for(Index id:bodyNodes[2]){check(!bodyNodes[1].contains(id),"contact bodies retain distinct node IDs");duplicates+=firstPositions.contains(bits(mesh.node(id)));}
    check(duplicates>=3,"actual interface has identical positions on distinct body nodes");
    for(int body=1;body<=2;++body){
        std::set<int> expectedVertices,expectedFaces;
        for(const auto& face:model.faceUses(body)){expectedFaces.insert(face.faceTag);
            for(const auto& wire:model.boundaryUses(face.faceTag))for(const auto& edge:wire.edges){expectedVertices.insert(edge.startVertexTag);expectedVertices.insert(edge.endVertexTag);}}
        for(Index id:bodyNodes[body]){const auto cls=mesh.nodeClass(id);const auto p=mesh.node(id);Vec3 q;
            if(cls.dim==3){check(cls.tag==body&&model.solid(body)->isInside(p),"interior node native body identity and inclusion");continue;}
            if(cls.dim==0){retainedVertices[body].insert(cls.tag);check(expectedVertices.contains(cls.tag),"original native vertex provenance");q=model.vertex(cls.tag)->position();}
            else if(cls.dim==1)q=model.edge(cls.tag)->closestPoint(p);
            else{check(cls.dim==2&&expectedFaces.contains(cls.tag)&&geo::closestPointOnOcctFace(model,cls.tag,p,q),"original native face provenance");}
            check((p-q).norm()<=std::max(8*model.tolerance(cls.dim,cls.tag),1e-9),"every classified node matches original native entity");
        }
        check(retainedVertices[body]==expectedVertices&&retainedFaces[body]==expectedFaces,"every original native vertex and face retained");
    }
    check(mesh::auditLinearVolumeBoundary(mesh).valid,"independent complete boundary incidence audit");
}
void refusals(const geo::OcctModel& model,const geo::NativeCadMeshOptions& options,
              const mesh::Mesh& accepted,const geo::NativeCadMeshOutcome& success){
    const auto before=snapshot(model);mesh::Mesh sentinel;sentinel.addNode({-7,-8,-9},{3,999});const auto saved=mesh::meshHash(sentinel);
    const auto expect=[&](const geo::NativeCadMeshOptions& o,geo::NativeCadMeshReason reason,const char* message){
        const auto r=buildOcctContactAssembly(model,o,sentinel);
        if(r.reason!=reason)std::cerr<<message<<" actual="<<static_cast<unsigned>(r.reason)<<" diagnostic="<<r.nativeContactDiagnostic<<'\n';
        check(!r.published&&r.destinationUnchanged&&r.reason==reason,message);
        check(mesh::meshHash(sentinel)==saved&&snapshot(model)==before,"failed assembly source/destination rollback");
    };
    auto o=options;check(success.validationWork>1,"measured cumulative work available");
    o.maximumValidationWork=success.validationWork-1;o.geometry.maximumMeshingWork=o.maximumValidationWork;
    expect(o,geo::NativeCadMeshReason::BudgetExceeded,"one-less cumulative measured work refuses");
    o=options;o.geometry.maximumNodes=accepted.nodeCount()-1;
    expect(o,geo::NativeCadMeshReason::BudgetExceeded,"one-less accepted total node count refuses");
    o=options;const auto bytes=static_cast<std::uint64_t>(accepted.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+
        static_cast<std::uint64_t>(accepted.cellCount())*sizeof(mesh::Cell)+static_cast<std::uint64_t>(accepted.connectivityCount())*sizeof(Index);
    o.geometry.maximumOutputBytes=bytes-1;expect(o,geo::NativeCadMeshReason::BudgetExceeded,"one-less exact accepted output payload refuses");
    o=options;check(success.nativeContactMaximumWorkingPayloadEstimate>1,"measured logical working high-water available");
    o.maximumPrimitiveWorkingBytes=success.nativeContactMaximumWorkingPayloadEstimate-1;
    expect(o,geo::NativeCadMeshReason::BudgetExceeded,"one-less measured logical working high-water refuses");
    std::atomic_bool cancel{true};o=options;o.geometry.cancellation=&cancel;
    expect(o,geo::NativeCadMeshReason::Cancelled,"cancelled assembly refuses before publication");
    o=options;o.geometry.minimumElementQuality=-.1;expect(o,geo::NativeCadMeshReason::InvalidRequest,"invalid quality request refused");
    o=options;o.geometry.minimumElementQuality=1;expect(o,geo::NativeCadMeshReason::QualityRequirementNotMet,"strict quality target refuses poor native star fan");
    failNextAllocation.store(true,std::memory_order_relaxed);
    const auto one=buildOcctContactAssembly(model,options,sentinel);failNextAllocation.store(false,std::memory_order_relaxed);
    check(!one.published&&one.reason==geo::NativeCadMeshReason::AllocationFailed&&mesh::meshHash(sentinel)==saved,"single allocation failure rolls back distinctly");
    failEveryAllocation.store(true,std::memory_order_relaxed);
    const auto persistent=buildOcctContactAssembly(model,options,sentinel);failEveryAllocation.store(false,std::memory_order_relaxed);
    check(!persistent.published&&persistent.reason==geo::NativeCadMeshReason::AllocationFailed&&mesh::meshHash(sentinel)==saved,"persistent allocation failure cannot terminate in diagnostic/source-audit cleanup");
    check(snapshot(model)==before,"allocation refusal source bytes unchanged");
}
void run(const TopoDS_Shape& shape,bool budgets){
    auto imported=OcctModelAccess::fromShape(shape);check(imported.imported()&&imported.model->solids().size()==2,"two valid independent native source bodies");
    const auto before=snapshot(*imported.model);const auto options=settings();mesh::Mesh mesh;
    const auto result=buildOcctContactAssembly(*imported.model,options,mesh);
    if(!result.published)std::cerr<<"assembly failure "<<static_cast<unsigned>(result.reason)<<' '<<result.nativeContactDiagnostic<<'\n';
    verify(*imported.model,mesh,result);check(before==snapshot(*imported.model),"entire native assembly source snapshot unchanged");
    if(budgets)refusals(*imported.model,options,mesh,result);
}
#ifndef PAACS_CONTACT_HELPER_ONLY
void publicContract(const TopoDS_Shape& shape,bool limits){
    auto imported=OcctModelAccess::fromShape(shape);check(imported.imported(),"public native contact source imports");
    const auto before=snapshot(*imported.model);auto options=settings();mesh::Mesh accepted;
    const auto success=geo::meshNativeCad(*imported.model,options,accepted);
    if(!success.published)std::cerr<<"public assembly failure "<<static_cast<unsigned>(success.reason)<<' '<<success.nativeContactDiagnostic<<'\n';
    verify(*imported.model,accepted,success);check(snapshot(*imported.model)==before,"public native source guard spans legacy and generic contact fallback");
    if(!limits)return;
    mesh::Mesh sentinel;sentinel.addNode({-3,-4,-5},{3,888});const auto saved=mesh::meshHash(sentinel);
    const auto expect=[&](const geo::NativeCadMeshOptions& o,geo::NativeCadMeshReason reason,const char* message){
        const auto r=geo::meshNativeCad(*imported.model,o,sentinel);
        if(r.reason!=reason)std::cerr<<message<<" actual="<<static_cast<unsigned>(r.reason)<<" diagnostic="<<r.nativeContactDiagnostic<<'\n';
        check(!r.published&&r.destinationUnchanged&&r.reason==reason,message);
        check(mesh::meshHash(sentinel)==saved&&snapshot(*imported.model)==before,"public legacy plus contact fallback refusal is atomic");
    };
    auto o=options;check(success.validationWork>1,"public exact cumulative work receipt");
    o.maximumValidationWork=success.validationWork-1;o.geometry.maximumMeshingWork=o.maximumValidationWork;
    expect(o,geo::NativeCadMeshReason::BudgetExceeded,"public one-less combined legacy/helper work refuses");
    o=options;const auto bytes=static_cast<std::uint64_t>(accepted.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+
        static_cast<std::uint64_t>(accepted.cellCount())*sizeof(mesh::Cell)+static_cast<std::uint64_t>(accepted.connectivityCount())*sizeof(Index);
    o.geometry.maximumOutputBytes=bytes-1;expect(o,geo::NativeCadMeshReason::BudgetExceeded,"public one-less joined mesh output refuses");
    o=options;o.minimumSize=.1;expect(o,geo::NativeCadMeshReason::UnsupportedControl,"generic contact unsupported minimum size is explicit");
    std::atomic_bool cancel{true};o=options;o.geometry.cancellation=&cancel;
    expect(o,geo::NativeCadMeshReason::Cancelled,"public generic contact cancellation refuses");
}
#endif
}
int main(){try{
    const auto boxA=BRepPrimAPI_MakeBox(1,1,1).Shape(),boxB=BRepPrimAPI_MakeBox(gp_Pnt(1,0,0),1,1,1).Shape();run(compound(boxA,boxB),true);
    const auto outer=BRepPrimAPI_MakeCylinder(2,2).Shape(),inner=BRepPrimAPI_MakeCylinder(1,2).Shape();
    BRepAlgoAPI_Cut tube(outer,inner);check(tube.IsDone(),"actual native tube-hole fixture");
    gp_Trsf translation;translation.SetTranslation(gp_Vec(0,0,2));BRepBuilderAPI_Transform upper(tube.Shape(),translation,true);
    check(upper.IsDone(),"independent coaxial cap-touch tube");run(compound(tube.Shape(),upper.Shape()),false);
#ifndef PAACS_CONTACT_HELPER_ONLY
    publicContract(compound(boxA,boxB),true);publicContract(compound(tube.Shape(),upper.Shape()),false);
#endif
    auto overlap=OcctModelAccess::fromShape(compound(boxA,BRepPrimAPI_MakeBox(gp_Pnt(.5,0,0),1,1,1).Shape()));check(overlap.imported(),"valid overlapping source BREP");
    const auto before=snapshot(*overlap.model);mesh::Mesh sentinel;sentinel.addNode({-1,-2,-3},{3,99});const auto saved=mesh::meshHash(sentinel);
    const auto refused=buildOcctContactAssembly(*overlap.model,settings(),sentinel);
    check(!refused.published&&refused.reason==geo::NativeCadMeshReason::UnresolvedBodyInterface&&mesh::meshHash(sentinel)==saved&&snapshot(*overlap.model)==before,"real native material overlap remains refused with atomic rollback");
    std::cout<<"native contact assembly: distinct body nodes/actual holes/native identity/all-point geometry/audits; overlap/work/output/working/cancel/quality/persistent-OOM rollback passed\n";return 0;
}catch(const std::exception& e){failEveryAllocation=false;failNextAllocation=false;std::cerr<<e.what()<<'\n';return 1;}catch(...){failEveryAllocation=false;failNextAllocation=false;std::cerr<<"unexpected native contact test failure\n";return 2;}}
