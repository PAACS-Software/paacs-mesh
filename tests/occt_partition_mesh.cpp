// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctPartitionMesher.hpp"
#include "OcctModelInternal.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <sstream>
#include <stdexcept>

static std::atomic_bool failNextAllocation{false};
void* operator new(std::size_t n){if(failNextAllocation.exchange(false))throw std::bad_alloc{};if(void* p=std::malloc(n?n:1))return p;throw std::bad_alloc{};}
void* operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void* p)noexcept{std::free(p);}
void operator delete[](void* p)noexcept{std::free(p);}
void operator delete(void* p,std::size_t)noexcept{std::free(p);}
void operator delete[](void* p,std::size_t)noexcept{std::free(p);}
using namespace paacs;using namespace paacs::geo;using namespace paacs::geo::detail;
namespace {
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
std::string bytes(const OcctModel& m){std::ostringstream out;BRepTools::Write(OcctModelAccess::shape(m),out);return out.str();}
mesh::Mesh sentinel(){mesh::Mesh m;m.addNode({91,92,93},{0,94});return m;}
void unchanged(const mesh::Mesh& m){require(m.nodeCount()==1&&m.cellCount()==0&&m.node(0).x==91&&m.nodeClass(0)==mesh::Classification{0,94},"destination mutated on refusal");}
void refused(const OcctModel& m,OcctPartitionMeshOptions o,OcctPartitionMeshReason reason){auto d=sentinel();const auto result=buildOcctPartitionVolume(m,1,o,d);require(!result.published&&result.destinationUnchanged&&result.reason==reason,result.diagnostic);unchanged(d);}
}
int main(){try{
    const auto outer=BRepPrimAPI_MakeBox(gp_Pnt(0,0,0),gp_Pnt(3,3,1)).Shape();
    const auto hole=BRepPrimAPI_MakeBox(gp_Pnt(1,1,-1),gp_Pnt(2,2,2)).Shape();
    BRepAlgoAPI_Cut cut(outer,hole);require(cut.IsDone(),"construct genus-one native solid");
    auto imported=OcctModelAccess::fromShape(cut.Shape());require(imported.imported()&&imported.solids==1,"snapshot genus-one native solid");
    const auto before=bytes(*imported.model);auto d=sentinel();
    OcctPartitionMeshOptions o;o.boundary.maximumChordError=.05;o.boundary.linearDeflection=.02;
    const auto result=buildOcctPartitionVolume(*imported.model,1,o,d);
    require(result.published&&!result.destinationUnchanged&&result.reason==OcctPartitionMeshReason::None,result.diagnostic);
    require(result.pieces>1&&result.interfaces>0&&result.pairedInterfaceTriangles>0,"partition/interface receipt");
    require(result.exteriorCertificate.certified&&result.exteriorCertificate.eulerCharacteristic==0&&result.exteriorCertificate.connectedShells==1,"original through-hole preserved");
    require(result.audit.valid&&result.quality.valid&&std::abs(result.meshVolume-8)<1e-10,"full original native mesh audit");
    require(result.sourceImmutabilityVerified&&before==bytes(*imported.model),"source native BRep mutated");
    {
        auto liveOptions=o;
        liveOptions.useNativeKernelAdmission=true;
        OcctPartitionRecoveryOptions recovery;
        recovery.enabled=true;
        recovery.useLiveConstraints=true;
        mesh::Mesh reference,withWitnesses;
        const auto slow=buildOcctPartitionVolume(*imported.model,1,liveOptions,recovery,reference);
        recovery.useAdaptiveRefill=true;
        const auto fast=buildOcctPartitionVolume(*imported.model,1,liveOptions,recovery,withWitnesses);
        require(slow.published&&fast.published&&fast.audit.valid&&fast.sourceImmutabilityVerified&&
                fast.sampleWork==4*fast.surfaceTriangles&&slow.sampleWork==fast.sampleWork&&
                fast.maximumBoundarySampleError<=o.boundary.maximumChordError&&
                mesh::meshHash(reference)==mesh::meshHash(withWitnesses)&&before==bytes(*imported.model),
                "original-face witnesses preserve all facet samples and the native mesh");
    }
    for(Index i=0;i<d.nodeCount();++i)if(d.nodeClass(i).dim==3)require(imported.model->solid(1)->isInside(d.node(i)),"volume node fills native hole");
    for(Index i=0;i<d.cellCount();++i)if(d.cell(i).type==mesh::CellType::Tet4){const auto* ids=d.cellNodes(i);const auto c=(d.node(ids[0])+d.node(ids[1])+d.node(ids[2])+d.node(ids[3]))*.25;require(imported.model->solid(1)->isInside(c),"tet fills native hole");}
    {
        auto parallelOptions=o;parallelOptions.maximumNativeValidationThreads=4;
        auto parallelMesh=sentinel();const auto parallel=buildOcctPartitionVolume(*imported.model,1,parallelOptions,parallelMesh);
        require(parallel.published&&parallel.nativeCentroidValidation.completed&&
                parallel.nativeCentroidValidation.sourceModelUnmodified&&
                parallel.nativeCentroidValidation.threadsUsed==4&&
                parallel.nativeCentroidValidation.nativeQueries==result.cells-result.surfaceTriangles&&
                parallel.nativeCentroidValidation.originalReferenceQueries>0&&
                parallel.nativeCentroidValidation.copiedGeometryValuePairs>0,
                "parallel native validation must check every Tet4 and compare owned copies to original CAD");
        require(parallelMesh.nodeCount()==d.nodeCount()&&parallelMesh.cellCount()==d.cellCount()&&
                parallel.meshVolume==result.meshVolume&&parallel.exteriorCertificate.certified&&parallel.audit.valid,
                "native validation workers changed completed mesh or audits");
        for(Index i=0;i<d.nodeCount();++i){const auto a=d.node(i),b=parallelMesh.node(i);
            require(a.x==b.x&&a.y==b.y&&a.z==b.z&&d.nodeClass(i)==parallelMesh.nodeClass(i),
                    "native validation workers changed node payload");}
        for(Index i=0;i<d.cellCount();++i){const auto& a=d.cell(i);const auto& b=parallelMesh.cell(i);
            require(a.type==b.type&&a.onEntity==b.onEntity,"native validation workers changed cell identity");
            for(unsigned j=0;j<mesh::nodesPerCell(a.type);++j)
                require(d.cellNodes(i)[j]==parallelMesh.cellNodes(i)[j],"native validation workers changed connectivity");}
        require(before==bytes(*imported.model),"native validation workers mutated source CAD");
        parallelOptions.maximumWork=parallel.workUsed-1;
        refused(*imported.model,parallelOptions,OcctPartitionMeshReason::BudgetExceeded);
        parallelOptions.maximumNativeValidationThreads=9;
        refused(*imported.model,parallelOptions,OcctPartitionMeshReason::InvalidRequest);
        parallelOptions.maximumNativeValidationThreads=0;
        refused(*imported.model,parallelOptions,OcctPartitionMeshReason::InvalidRequest);
    }
    o.maximumWork=0;refused(*imported.model,o,OcctPartitionMeshReason::BudgetExceeded);
    o={};o.boundary.maximumWork=0;refused(*imported.model,o,OcctPartitionMeshReason::BudgetExceeded);
    o={};o.maximumNodes=1;refused(*imported.model,o,OcctPartitionMeshReason::BudgetExceeded);
    o={};o.maximumOutputBytes=1;refused(*imported.model,o,OcctPartitionMeshReason::BudgetExceeded);
    o={};o.maximumWorkingBytes=1;refused(*imported.model,o,OcctPartitionMeshReason::BudgetExceeded);
    o={};o.partition.maximumPieces=1;refused(*imported.model,o,OcctPartitionMeshReason::BudgetExceeded);
    o={};o.maximumPieceStarWork=1;d=sentinel();
    const auto capped=buildOcctPartitionVolume(*imported.model,1,o,d);
    require(!capped.published&&capped.destinationUnchanged&&capped.reason==OcctPartitionMeshReason::BudgetExceeded
        &&capped.starReason==OcctStarReason::WorkBudget&&capped.starWork<=o.maximumPieceStarWork
        &&capped.workUsed<o.maximumWork,"per-piece star admission must preserve unused global work");unchanged(d);
    std::atomic_bool cancelled{true};o={};o.cancellation=&cancelled;refused(*imported.model,o,OcctPartitionMeshReason::Cancelled);
    d=sentinel();failNextAllocation=true;const auto allocation=buildOcctPartitionVolume(*imported.model,1,{},d);require(allocation.reason==OcctPartitionMeshReason::AllocationFailed&&!allocation.published,"allocation taxonomy");unchanged(d);
    require(before==bytes(*imported.model),"failed attempts mutated source");
    // A curved bore surrounded by convex corner fillets is a different
    // admission problem from the square-hole fixture: cutting every cylinder
    // axis exhausts this small plane allowance, while the concave-only policy
    // retains the actual bore and joins its native sector interfaces.
    {
        TopoDS_Shape rounded=BRepPrimAPI_MakeBox(gp_Pnt(.5,0,0),gp_Pnt(5.5,6,2)).Shape();
        auto unite=[&](const TopoDS_Shape& shape){
            BRepAlgoAPI_Fuse f(rounded,shape);f.SetNonDestructive(true);f.Build();
            require(f.IsDone(),"rounded-bore fixture fuse");f.SimplifyResult();rounded=f.Shape();
        };
        unite(BRepPrimAPI_MakeBox(gp_Pnt(0,.5,0),gp_Pnt(6,5.5,2)).Shape());
        for(double x:{.5,5.5})for(double y:{.5,5.5})
            unite(BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x,y,0),gp_Dir(0,0,1)),.5,2).Shape());
        BRepAlgoAPI_Cut bore(rounded,BRepPrimAPI_MakeCylinder(
            gp_Ax2(gp_Pnt(3,3,-1),gp_Dir(0,0,1)),1.2,4).Shape());
        bore.SetNonDestructive(true);bore.Build();require(bore.IsDone(),"rounded-bore fixture cut");
        bore.SimplifyResult();auto native=OcctModelAccess::fromShape(bore.Shape());
        require(native.imported()&&native.solids==1,"rounded-bore native snapshot");
        const auto original=bytes(*native.model);
        OcctPartitionMeshOptions selective;selective.boundary.linearDeflection=.02;
        selective.boundary.maximumChordError=.05;selective.partition.maximumPlanes=12;
        selective.partition.featureAxisPlanes=8;
        selective.partition.policy=OcctPartitionPolicy::NativeFeatureAxes;
        refused(*native.model,selective,OcctPartitionMeshReason::BudgetExceeded);
        selective.partition.policy=OcctPartitionPolicy::NativeConcaveFeatureAxes;
        auto selectiveMesh=sentinel();const auto built=buildOcctPartitionVolume(*native.model,1,selective,selectiveMesh);
        require(built.published&&built.reason==OcctPartitionMeshReason::None,built.diagnostic);
        require(built.planes<=12&&built.pieces>1&&built.pairedInterfaceTriangles>0,
                "concave-only native sectors must share audited interfaces");
        require(built.exteriorCertificate.certified&&built.exteriorCertificate.eulerCharacteristic==0&&
                built.audit.valid&&built.quality.valid&&built.relativeVolumeError<.01,
                "concave-only mesh preserves original bore and full volume audit");
        require(built.sourceImmutabilityVerified&&original==bytes(*native.model),
                "concave-only publication mutated original rounded solid");
        for(Index i=0;i<selectiveMesh.cellCount();++i)if(selectiveMesh.cell(i).type==mesh::CellType::Tet4){
            const auto* ids=selectiveMesh.cellNodes(i);const auto p=(selectiveMesh.node(ids[0])+selectiveMesh.node(ids[1])+selectiveMesh.node(ids[2])+selectiveMesh.node(ids[3]))*.25;
            require(native.model->solid(1)->isInside(p),"concave-only fan centroid enters curved bore");
        }
        selective.maximumAdaptiveCuts=4;auto adaptiveMesh=sentinel();
        const auto adaptive=buildOcctPartitionVolume(*native.model,1,selective,adaptiveMesh);
        require(adaptive.published&&adaptive.adaptiveAttempts==1&&adaptive.adaptiveCutReports.empty()&&
                adaptive.adaptivePrincipalFaceTag>0&&adaptive.adaptiveWork>0,
                "adaptive admission retains a successful first attempt");
        require(adaptiveMesh.nodeCount()==selectiveMesh.nodeCount()&&adaptiveMesh.cellCount()==selectiveMesh.cellCount(),
                "adaptive admission changed first-attempt mesh counts");
        for(Index i=0;i<selectiveMesh.nodeCount();++i){
            const auto a=selectiveMesh.node(i),b=adaptiveMesh.node(i);
            require(a.x==b.x&&a.y==b.y&&a.z==b.z&&selectiveMesh.nodeClass(i)==adaptiveMesh.nodeClass(i),
                    "adaptive admission changed first-attempt node identity or coordinates");
        }
        for(Index i=0;i<selectiveMesh.cellCount();++i){
            const auto& a=selectiveMesh.cell(i);const auto& b=adaptiveMesh.cell(i);
            require(a.type==b.type&&a.onEntity==b.onEntity,"adaptive admission changed first-attempt cell identity");
            for(unsigned j=0;j<mesh::nodesPerCell(a.type);++j)
                require(selectiveMesh.cellNodes(i)[j]==adaptiveMesh.cellNodes(i)[j],
                        "adaptive admission changed first-attempt connectivity");
        }
        selective.maximumWork=adaptive.workUsed-1;d=sentinel();
        const auto oneLess=buildOcctPartitionVolume(*native.model,1,selective,d);
        require(!oneLess.published&&oneLess.destinationUnchanged&&oneLess.reason==OcctPartitionMeshReason::BudgetExceeded&&
                oneLess.workUsed<=selective.maximumWork,"adaptive work admission is cumulative and transactional");unchanged(d);
        selective.maximumWork=5'000'000;selective.maximumAdaptiveCuts=17;
        refused(*native.model,selective,OcctPartitionMeshReason::InvalidRequest);
        selective.maximumAdaptiveCuts=16;selective.useNativeKernelAdmission=true;
        auto nativeAdmittedMesh=sentinel();const auto nativeAdmitted=buildOcctPartitionVolume(*native.model,1,selective,nativeAdmittedMesh);
        require(nativeAdmitted.published&&nativeAdmitted.adaptiveAttempts==1&&nativeAdmitted.nativeCandidateQueries>0&&nativeAdmitted.nativeCandidateQueries<=nativeAdmitted.starWork&&nativeAdmitted.audit.valid&&nativeAdmitted.sourceImmutabilityVerified,"native candidate admission remains fully audited within adaptive16");
        selective.maximumWork=nativeAdmitted.workUsed-1;refused(*native.model,selective,OcctPartitionMeshReason::BudgetExceeded);
        selective.maximumWork=5'000'000;selective.useNativeKernelAdmission=false;
        selective.maximumAdaptiveCuts=4;selective.cancellation=&cancelled;
        refused(*native.model,selective,OcctPartitionMeshReason::Cancelled);selective.cancellation=nullptr;
        d=sentinel();failNextAllocation=true;
        const auto adaptiveAllocation=buildOcctPartitionVolume(*native.model,1,selective,d);
        require(!adaptiveAllocation.published&&adaptiveAllocation.reason==OcctPartitionMeshReason::AllocationFailed,
                "adaptive allocation failure must preserve publication boundary");unchanged(d);
        selective.maximumWork=0;refused(*native.model,selective,OcctPartitionMeshReason::BudgetExceeded);
        require(original==bytes(*native.model),"concave-only refusal mutated original rounded solid");
    }
    std::puts("Native partition-volume interface, hole, budget, cancellation and rollback checks passed");return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
