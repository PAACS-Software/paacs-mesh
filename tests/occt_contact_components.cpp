// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctContactComponents.hpp"
#include "OcctContactAssembly.hpp"
#include "OcctModelInternal.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Pln.hxx>
#include <atomic>
#include <array>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <new>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {std::atomic_bool failNext{false},failEvery{false};}
void* operator new(std::size_t n){if(failEvery.load()||failNext.exchange(false))throw std::bad_alloc{};
    if(void* p=std::malloc(n?n:1))return p;throw std::bad_alloc{};}
void* operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void* p)noexcept{std::free(p);}void operator delete[](void* p)noexcept{std::free(p);}
void operator delete(void* p,std::size_t)noexcept{std::free(p);}void operator delete[](void* p,std::size_t)noexcept{std::free(p);}
using namespace paacs;
using namespace paacs::geo::detail;
namespace {
void check(bool b,const char* s){if(!b)throw std::runtime_error(s);}
std::string snapshot(const geo::OcctModel& model){std::ostringstream s;BRepTools::Write(OcctModelAccess::shape(model),s);return s.str();}
geo::NativeCadMeshOptions options(){geo::NativeCadMeshOptions o;o.geometry.targetSize=.4;o.geometry.warningElementQuality=.25;
    o.maximumChordError=.05;o.separateMatchedNativeFaceContacts=true;o.separateLowerDimensionalContacts=true;o.useNativeKernelAdmission=true;
    o.maximumValidationWork=20000000;o.geometry.maximumMeshingWork=20000000;o.maximumNativePartitionPieces=512;
    o.maximumNativeValidationThreads=2;o.maximumPrimitiveWorkingBytes=128ULL*1024*1024;
    o.geometry.maximumNodes=20000;o.geometry.maximumCells=100000;o.geometry.maximumOutputBytes=32ULL*1024*1024;
    o.geometry.maximumAuditTemporaryBytes=32ULL*1024*1024;return o;}
TopoDS_Shape source(bool shared=false,bool overlap=false){TopoDS_Compound c;BRep_Builder b;b.MakeCompound(c);
    if(shared){NCollection_List<TopoDS_Shape> args,tools;args.Append(BRepPrimAPI_MakeBox(2,1,1).Shape());
        tools.Append(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(1,0,0),gp_Dir(1,0,0)),-3,3,-3,3).Shape());
        BRepAlgoAPI_Splitter split;split.SetArguments(args);split.SetTools(tools);split.SetNonDestructive(true);split.Build();check(split.IsDone()&&!split.HasErrors(),"shared original native interface fixture");b.Add(c,split.Shape());
    }else {b.Add(c,BRepPrimAPI_MakeBox(1,1,1).Shape());b.Add(c,BRepPrimAPI_MakeBox(gp_Pnt(overlap?.75:1,0,0),1,1,1).Shape());}
    b.Add(c,BRepPrimAPI_MakeBox(gp_Pnt(4,0,0),1,1,1).Shape());return c;
}
void verify(const geo::OcctModel& model,const mesh::Mesh& value,const geo::NativeCadMeshOutcome& r,bool component=true){
    if(!r.published)std::cerr<<"component refusal "<<geo::nativeCadMeshReasonName(r.reason)<<" "<<r.nativeContactDiagnostic<<'\n';
    check(r.published&&r.reason==geo::NativeCadMeshReason::None&&r.audit.valid&&r.quality.complete&&r.quality.valid,"complete original assembly is valid");
    check(r.nativeContactSourceGuard&&r.nativeContactSourceAuditAvailable&&r.matchedNativeContactTriangles>0&&r.bodies.size()==3,"complete source and interface proof");
    std::set<int> bodyTags;for(const auto& body:r.bodies){bodyTags.insert(body.solidTag);
        if(component)check(body.nativeCentroidSourceGuard&&body.nativeCentroidsChecked>0,"every original body gets independent native-worker proof");
        check(body.boundaryEmbeddingCertified,"every original body gets exterior embedding proof");
        check(std::abs(body.meshVolume-body.cadVolume)<body.cadVolume*.01,"native original-body volume");}
    check(bodyTags==std::set<int>{1,2,3},"three original body identities");
    std::array<std::set<Index>,4> bodyNodes;for(Index i=0;i<value.cellCount();++i){const auto& cell=value.cell(i);if(cell.type!=mesh::CellType::Tet4)continue;
        check(cell.onEntity.dim==3&&cell.onEntity.tag>=1&&cell.onEntity.tag<=3,"tet original body owner");Vec3 p{};
        for(int j=0;j<4;++j){const auto id=value.cellNodes(i)[j];p=p+value.node(id);bodyNodes[cell.onEntity.tag].insert(id);}
        check(model.solid(cell.onEntity.tag)->isInside(p*.25),"every original tet centroid inside");}
    for(int body=1;body<=3;++body){check(!bodyNodes[body].empty(),"all original material bodies represented");
        for(Index id:bodyNodes[body]){check(model.solid(body)->isInside(value.node(id)),"every original node including hubs inside");
            for(int other=body+1;other<=3;++other)check(!bodyNodes[other].contains(id),"contact and isolated bodies keep distinct node IDs");}}
    check(mesh::auditLinearVolumeBoundary(value).valid,"independent final complete incidence audit");
}
TopoDS_Shape fallbackSource(bool overlap=false){
    BRepAlgoAPI_Cut tube(BRepPrimAPI_MakeCylinder(1.1,1.).Shape(),BRepPrimAPI_MakeCylinder(1.,1.).Shape());
    check(tube.IsDone()&&!tube.HasErrors(),"thin native tube fixture");
    TopoDS_Compound value;BRep_Builder builder;builder.MakeCompound(value);builder.Add(value,tube.Shape());
    builder.Add(value,BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0,0,overlap?.9:1.),gp_Dir(0,0,1)),1.1,.25).Shape());
    builder.Add(value,BRepPrimAPI_MakeBox(gp_Pnt(4,4,0),.5,.5,.5).Shape());return value;
}
void preservesWholeAssemblyAfterComponentRefusal(){
    auto imported=OcctModelAccess::fromShape(fallbackSource());check(imported.imported(),"tube, touching cap and isolated body import");
    const auto& model=*imported.model;const auto before=snapshot(model);auto o=options();
    o.geometry.targetSize=.5;o.maximumChordError=.005;o.maximumNativePartitionPieces=64;
    mesh::Mesh sentinel;sentinel.addNode({-7,-8,-9},{3,777});const auto saved=mesh::meshHash(sentinel);
    const auto component=buildOcctContactComponents(model,o,sentinel);
    check(!component.published&&component.destinationUnchanged&&component.reason==geo::NativeCadMeshReason::BudgetExceeded&&
          std::string_view(component.nativeContactDiagnostic)=="component_child_not_published",
          "component strategy reaches its native piece limit");
    check(mesh::meshHash(sentinel)==saved&&snapshot(model)==before,"failed strategy leaves source and destination unchanged");
    auto previous=o;previous.useNativeKernelAdmission=false;mesh::Mesh expected;
    const auto whole=buildOcctContactAssembly(model,previous,expected);verify(model,expected,whole,false);
    check(whole.nativeContactPieces<=o.maximumNativePartitionPieces,"established whole strategy fits the same piece limit");
    mesh::Mesh normal;const auto published=geo::meshNativeCad(model,o,normal);verify(model,normal,published,false);
    check(mesh::meshHash(normal)==mesh::meshHash(expected)&&published.validationWork>component.validationWork+whole.validationWork,
          "advanced failure retains the established mesh and charges both attempts");
    check(snapshot(model)==before,"fallback preserves entire original source");
    auto request=o;request.maximumValidationWork=published.validationWork-1;request.geometry.maximumMeshingWork=request.maximumValidationWork;
    auto refused=geo::meshNativeCad(model,request,sentinel);
    check(!refused.published&&refused.destinationUnchanged&&refused.reason==geo::NativeCadMeshReason::BudgetExceeded&&
          mesh::meshHash(sentinel)==saved&&snapshot(model)==before,"one-less cumulative work cannot publish a fallback");
    request=o;request.maximumNativePartitionPieces=32;refused=geo::meshNativeCad(model,request,sentinel);
    check(!refused.published&&refused.reason==geo::NativeCadMeshReason::BudgetExceeded&&mesh::meshHash(sentinel)==saved&&snapshot(model)==before,
          "two strategy failures preserve destination and source");
    request=o;request.geometry.minimumElementQuality=1;refused=geo::meshNativeCad(model,request,sentinel);
    check(!refused.published&&refused.reason==geo::NativeCadMeshReason::QualityRequirementNotMet&&mesh::meshHash(sentinel)==saved,
          "fallback cannot bypass requested strict quality");
    std::atomic_bool cancelled{true};request=o;request.geometry.cancellation=&cancelled;refused=geo::meshNativeCad(model,request,sentinel);
    check(!refused.published&&refused.reason==geo::NativeCadMeshReason::Cancelled&&mesh::meshHash(sentinel)==saved&&snapshot(model)==before,
          "cancelled fallback request preserves source and destination");
    auto overlapping=OcctModelAccess::fromShape(fallbackSource(true));check(overlapping.imported(),"overlap refusal fixture import");
    const auto overlappingBefore=snapshot(*overlapping.model);refused=geo::meshNativeCad(*overlapping.model,o,sentinel);
    check(!refused.published&&refused.reason==geo::NativeCadMeshReason::UnresolvedBodyInterface&&mesh::meshHash(sentinel)==saved&&
          snapshot(*overlapping.model)==overlappingBefore,"strategy fallback cannot admit original material overlap");
}
void run(bool shared,bool limits){auto imported=OcctModelAccess::fromShape(source(shared));check(imported.imported()&&imported.model->solids().size()==3,"three native source bodies");
    std::cerr<<"component fixture shared="<<shared<<" limits="<<limits<<'\n';
    const auto before=snapshot(*imported.model);auto o=options();mesh::Mesh value;const auto success=buildOcctContactComponents(*imported.model,o,value);verify(*imported.model,value,success);
    mesh::Mesh normal;const auto published=geo::meshNativeCad(*imported.model,o,normal);verify(*imported.model,normal,published);
    check(mesh::meshHash(normal)==mesh::meshHash(value)&&published.validationWork>success.validationWork,
          "normal entry point preserves complete component output and charges prior work");
    check(snapshot(*imported.model)==before,"entire original geometry stays byte-identical");if(!limits)return;
    mesh::Mesh sentinel;sentinel.addNode({-4,-5,-6},{3,888});const auto saved=mesh::meshHash(sentinel);
    const auto expect=[&](const geo::NativeCadMeshOptions& request,geo::NativeCadMeshReason reason,const char* label){std::cerr<<label<<'\n';auto r=buildOcctContactComponents(*imported.model,request,sentinel);
        if(r.reason!=reason)std::cerr<<label<<" actual="<<geo::nativeCadMeshReasonName(r.reason)<<" "<<r.nativeContactDiagnostic<<'\n';
        check(!r.published&&r.destinationUnchanged&&r.reason==reason,label);check(mesh::meshHash(sentinel)==saved&&snapshot(*imported.model)==before,"typed refusal preserves destination and original");};
    auto request=o;request.maximumValidationWork=success.validationWork-1;request.geometry.maximumMeshingWork=request.maximumValidationWork;
    expect(request,geo::NativeCadMeshReason::BudgetExceeded,"one-less measured global work");
    request=o;request.geometry.maximumNodes=value.nodeCount()-1;expect(request,geo::NativeCadMeshReason::BudgetExceeded,"one-less whole node count");
    request=o;request.geometry.maximumOutputBytes=static_cast<std::uint64_t>(value.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+
        static_cast<std::uint64_t>(value.cellCount())*sizeof(mesh::Cell)+static_cast<std::uint64_t>(value.connectivityCount())*sizeof(Index)-1;
    expect(request,geo::NativeCadMeshReason::BudgetExceeded,"one-less whole output payload");
    request=o;request.maximumPrimitiveWorkingBytes=1024;expect(request,geo::NativeCadMeshReason::BudgetExceeded,"bounded native snapshots and logical payload");
    request=o;request.geometry.minimumElementQuality=-.1;expect(request,geo::NativeCadMeshReason::InvalidRequest,"invalid quality request");
    request=o;request.geometry.minimumElementQuality=1;expect(request,geo::NativeCadMeshReason::QualityRequirementNotMet,"strict preferred quality");
    std::atomic_bool cancel{true};request=o;request.geometry.cancellation=&cancel;expect(request,geo::NativeCadMeshReason::Cancelled,"cancel before component work");
    request=o;request.minimumSize=.1;expect(request,geo::NativeCadMeshReason::UnsupportedControl,"unhandled controls explicit");
    const auto wrong=before+"changed";auto r=buildOcctContactComponents(*imported.model,o,sentinel,&wrong);
    check(!r.published&&r.reason==geo::NativeCadMeshReason::InternalError&&mesh::meshHash(sentinel)==saved,"caller source snapshot cannot be bypassed");
    std::cerr<<"single allocation failure\n";failNext.store(true);r=buildOcctContactComponents(*imported.model,o,sentinel);failNext.store(false);
    check(!r.published&&r.reason==geo::NativeCadMeshReason::AllocationFailed&&mesh::meshHash(sentinel)==saved,"single allocation failure typed and atomic");
    std::cerr<<"persistent allocation failure\n";failEvery.store(true);r=buildOcctContactComponents(*imported.model,o,sentinel);failEvery.store(false);
    check(!r.published&&r.reason==geo::NativeCadMeshReason::AllocationFailed&&mesh::meshHash(sentinel)==saved,"persistent allocation failure cleanup safe");
    check(snapshot(*imported.model)==before,"fault tests keep original bytes");
    request=o;request.maximumValidationWork=published.validationWork-1;request.geometry.maximumMeshingWork=request.maximumValidationWork;
    r=geo::meshNativeCad(*imported.model,request,sentinel);
    check(!r.published&&r.reason==geo::NativeCadMeshReason::BudgetExceeded&&mesh::meshHash(sentinel)==saved&&snapshot(*imported.model)==before,
          "normal entry point reserves and charges final original snapshot before publication");
}
}
int main(){std::setvbuf(stdout,nullptr,_IONBF,0);try{run(false,true);run(true,false);preservesWholeAssemblyAfterComponentRefusal();auto imported=OcctModelAccess::fromShape(source(false,true));check(imported.imported(),"overlapping individually valid bodies");
    mesh::Mesh value;value.addNode({-1,-2,-3},{3,99});const auto saved=mesh::meshHash(value);const auto before=snapshot(*imported.model);
    auto r=buildOcctContactComponents(*imported.model,options(),value);check(!r.published&&r.reason==geo::NativeCadMeshReason::UnresolvedBodyInterface&&mesh::meshHash(value)==saved&&snapshot(*imported.model)==before,"actual original material overlap refused");
    std::cout<<"Native contact components and original assembly guards passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
