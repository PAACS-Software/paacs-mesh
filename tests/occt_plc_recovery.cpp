// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctRecoveryMesher.hpp"
#include "OcctPlcFrontRecovery.hpp"
#include "OcctPartitionMesher.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctModelInternal.hpp"
#include "OcctNativeWireAuthority.hpp"
#include "../src/bnd/TetEdgeGeneralCavityCadPolicy.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/geometry/OcctMesher.hpp"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <gp_Pnt.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepTools.hxx>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <sstream>
#include <map>
#include <limits>
#include <set>
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
void require(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
std::string bytes(const OcctModel& m){std::ostringstream out;BRepTools::Write(OcctModelAccess::shape(m),out);return out.str();}
std::vector<std::uint8_t> archive(const mesh::Mesh& m){std::vector<std::uint8_t> b;require(sdk::serializeMeshArchive(m,b).completed,"serialize fixture");return b;}
mesh::Mesh sentinel(){mesh::Mesh m;m.addNode({91,92,93},{0,94});return m;}
mesh::Mesh boundary(const OcctModel& m){OcctBoundaryOptions o;o.linearDeflection=.08;o.maximumChordError=.2;OcctBoundaryOutcome prepared;auto session=OcctBoundarySession::prepare(m,o,prepared);require(bool(session),occtBoundaryReasonName(prepared.reason));mesh::Mesh b;auto x=session->boundary(1,b);require(x.published,occtBoundaryReasonName(x.reason));return b;}
std::vector<mesh::Classification> parents(const OcctModel& m){std::vector<mesh::Classification> p(1);for(int tag:m.solid(1)->boundaryFaceTags()){if(std::size_t(tag)>=p.size())p.resize(tag+1);p[tag]={2,tag};}return p;}
OcctStarOptions options(){OcctStarOptions o;o.maximumWork=30'000'000;o.maximumWorkingBytes=512ULL*1024*1024;return o;}
void nativeWireAuthority(){
    auto cad=OcctModelAccess::fromShape(BRepPrimAPI_MakeCylinder(1.,2.).Shape());require(cad.imported(),"wire authority cylinder");
    const auto source=bytes(*cad.model);const auto b=boundary(*cad.model);
    std::map<std::array<Index,2>,std::set<int>> faces;
    for(Index i=0;i<b.cellCount();++i){const auto* ids=b.cellNodes(i);for(unsigned j=0;j<3;++j){std::array<Index,2> key{ids[j],ids[(j+1)%3]};std::sort(key.begin(),key.end());faces[key].insert(b.cell(i).onEntity.tag);}}
    std::uint64_t visits=0;unsigned diagonals=0,rims=0;const auto tick=[&]{++visits;};
    for(const auto& [edge,owners]:faces){const auto a=b.nodeClass(edge[0]),c=b.nodeClass(edge[1]);
        if(owners.size()==1&&a.dim==1&&c.dim==1&&a.tag==c.tag&&cad.model->face(*owners.begin())->kind()==SurfaceKind::Plane){
            require(actualOcctWire(*cad.model,a,c,owners,tick)==0,"cap diagonal incorrectly claims native rim authority");++diagonals;}
        if(owners.size()==2&&a.dim==1&&c.dim==1&&a.tag==c.tag){
            require(actualOcctWire(*cad.model,a,c,owners,tick)==a.tag,"actual cross-face rim lost authority");++rims;}
    }
    require(diagonals>0&&rims>0&&visits>0,"wire authority fixture exercises internal cap chords and actual rims");
    unsigned seams=0;
    for(const auto* face:cad.model->faces())for(const auto& wire:cad.model->boundaryUses(face->tag()))for(const auto& use:wire.edges)if(use.seam&&!use.degenerated){
        const std::array owners{face->tag()};require(actualOcctWire(*cad.model,{0,use.startVertexTag},{0,use.endVertexTag},owners,tick)==use.edgeTag,"periodic same-face seam lost repeated-use authority");++seams;
    }
    require(seams==2&&bytes(*cad.model)==source,"periodic seam/source authority checks");
}
void recoveryWithoutBooleanCuts(){
    const auto shape=BRepAlgoAPI_Cut(BRepPrimAPI_MakeCone(2.,1.8,2.).Shape(),
        BRepPrimAPI_MakeCone(1.9,1.7,2.).Shape()).Shape();
    auto cad=OcctModelAccess::fromShape(shape);require(cad.imported(),"uncut hollow frustum snapshot");
    const auto source=bytes(*cad.model);
    NativeCadMeshOptions o;o.geometry.targetSize=.3;o.maximumChordError=.01;
    o.geometry.maximumMeshingWork=100000000;o.maximumValidationWork=100000000;
    o.maximumPrimitiveWorkingBytes=512ULL*1024*1024;o.useNativeKernelAdmission=true;
    o.useNativePartitionFallback=false;o.useLongitudinalRecovery=true;o.reuseNativeQueries=true;
    mesh::Mesh volume;const auto recovered=meshNativeCad(*cad.model,o,volume);
    require(recovered.published&&recovered.audit.valid&&recovered.quality.valid,
            "native recovery was disabled together with Boolean partitioning");
    require(recovered.bodies.size()==1&&recovered.bodies[0].partitionAttempts==0,
            "uncut recovery attempted a Boolean partition");
    require(bytes(*cad.model)==source,"uncut recovery mutated the original CAD");
}
void refuse(const OcctModel& model,const mesh::Mesh& originalBoundary,OcctStarOptions o,OcctPartitionRecoveryOptions r,OcctPlcRecoveryReason reason){auto b=originalBoundary,d=sentinel();auto oldBoundary=archive(b),oldDestination=archive(d);const auto x=recoverOcctPlcVolume(b,1,model,model,1,parents(model),o,r,d);require(!x.published&&x.reason==reason,x.diagnostic);require(oldBoundary==archive(b)&&oldDestination==archive(d),"failed recovery published boundary or volume");}
void smallNonstarFront(){
    const auto shape=BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(4.,3.,1.).Shape(),
        BRepPrimAPI_MakeBox(gp_Pnt(1.,1.,-.5),2.,3.,2.).Shape()).Shape();
    auto cad=OcctModelAccess::fromShape(shape);require(cad.imported(),"U-shaped piece snapshot");
    const auto original=bytes(*cad.model);const auto b=boundary(*cad.model);auto o=options();
    mesh::Mesh star;const auto unsupported=buildOcctStarVolume(b,1,o,star);
    require(!unsupported.published&&unsupported.reason==OcctStarReason::KernelNotFound,"U-piece should require several tetrahedral fronts");
    OcctPartitionRecoveryOptions recovery;recovery.enabled=true;recovery.useAdvancingFront=true;
    auto boundaryCopy=b,d=sentinel();const auto source=archive(b);
    const auto result=recoverOcctPlcVolume(boundaryCopy,1,*cad.model,*cad.model,1,parents(*cad.model),o,recovery,d);
    require(result.published&&result.piecesRecovered==1&&result.coreAttempts==0,result.diagnostic);
    const auto quality=algo::analyzeMixedMeshQuality(d);
    require(mesh::auditLinearVolumeBoundary(d).valid&&quality.valid&&std::abs(quality.totalVolume-8.)<1e-10,"nonstar front boundary/volume audit");
    require(archive(boundaryCopy)==source&&bytes(*cad.model)==original,"nonstar front mutated source");
    auto longer=recovery;longer.maximumFrontSeconds=15.;boundaryCopy=b;
    const auto extended=recoverOcctPlcVolume(boundaryCopy,1,*cad.model,*cad.model,1,parents(*cad.model),o,longer,d);
    require(extended.published&&mesh::auditLinearVolumeBoundary(d).valid&&
            std::abs(algo::analyzeMixedMeshQuality(d).totalVolume-8.)<1e-10,"extended front changed the U-piece domain");
    auto invalid=recovery;invalid.maximumFrontSeconds=0.;refuse(*cad.model,b,o,invalid,OcctPlcRecoveryReason::InvalidRequest);
    invalid.maximumFrontSeconds=std::numeric_limits<double>::quiet_NaN();refuse(*cad.model,b,o,invalid,OcctPlcRecoveryReason::InvalidRequest);
    // The axis-free route must join artificial subdivisions without bridging
    // the open slot of a nonstar solid or changing the original CAD.
    OcctPartitionMeshOptions divided;divided.partition.policy=OcctPartitionPolicy::CartesianGrid;
    divided.partition.cartesianDivisions=3;divided.maximumWork=30'000'000;
    divided.maximumPieceStarWork=5'000'000;divided.maximumWorkingBytes=o.maximumWorkingBytes;
    divided.useNativeKernelAdmission=true;divided.repairSurfaceAliases=true;
    divided.boundary.linearDeflection=.08;divided.boundary.maximumChordError=.2;
    for(unsigned i=1;i<8;++i)divided.partition.extraPlanes.push_back({{.5*i,0,0},{1,0,0}});
    mesh::Mesh joined;const auto partitioned=buildOcctPartitionVolume(*cad.model,1,divided,longer,joined);
    require(partitioned.published&&partitioned.interfaces>0&&partitioned.exteriorCertificate.certified,
            partitioned.diagnostic);
    require(mesh::auditLinearVolumeBoundary(joined).valid&&
            std::abs(algo::analyzeMixedMeshQuality(joined).totalVolume-8.)<1e-10,
            "Cartesian recovery filled the U-piece slot or left an interface gap");
    require(bytes(*cad.model)==original&&archive(boundaryCopy)==source,"Cartesian/front route changed source");
    auto limited=o;limited.maximumOutputBytes=1;refuse(*cad.model,b,limited,recovery,OcctPlcRecoveryReason::BudgetExceeded);
    limited=o;limited.maximumWorkingBytes=1;refuse(*cad.model,b,limited,recovery,OcctPlcRecoveryReason::BudgetExceeded);
    limited=o;limited.maximumWork=1;refuse(*cad.model,b,limited,recovery,OcctPlcRecoveryReason::BudgetExceeded);
    std::atomic_bool cancelled{true};limited=o;limited.cancellation=&cancelled;refuse(*cad.model,b,limited,recovery,OcctPlcRecoveryReason::Cancelled);
    auto broken=b;const auto* f=b.cellNodes(0);broken.addCell(mesh::CellType::Tri3,std::span<const Index>(f,3),b.cell(0).onEntity);
    refuse(*cad.model,broken,o,recovery,OcctPlcRecoveryReason::UnsupportedGeometry);
}
}
int main(){try{
    nativeWireAuthority();recoveryWithoutBooleanCuts();smallNonstarFront();
    auto imported=OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(1,1,1).Shape());require(imported.imported(),"box snapshot");
    const auto source=bytes(*imported.model);const auto initial=boundary(*imported.model);auto b=initial,d=sentinel();
    auto o=options();OcctPartitionRecoveryOptions r;r.enabled=true;
    const auto result=recoverOcctPlcVolume(b,1,*imported.model,*imported.model,1,parents(*imported.model),o,r,d);
    require(result.published&&result.reason==OcctPlcRecoveryReason::None&&result.piecesRecovered==1&&result.coreAttempts>0&&result.generalWorkReserved>0,result.diagnostic);
    require(mesh::auditLinearVolumeBoundary(d).valid&&algo::analyzeMixedMeshQuality(d).valid,"fixed PLC full audit");
    require(bytes(*imported.model)==source,"native source changed");
    const auto first=archive(d);b=initial;d=sentinel();const auto replay=recoverOcctPlcVolume(b,1,*imported.model,*imported.model,1,parents(*imported.model),o,r,d);
    require(replay.published&&replay.workUsed==result.workUsed&&first==archive(d),"recovery replay changed mesh or work");
    auto staged=r;staged.maximumSteinerCandidates=8;b=initial;d=sentinel();
    const auto optional=recoverOcctPlcVolume(b,1,*imported.model,*imported.model,1,parents(*imported.model),o,staged,d);
    require(optional.published&&optional.steinerNodes==0&&optional.workUsed==result.workUsed&&archive(d)==first,
            "enabling the staged fallback changed an existing fixed-PLC success");
    staged.maximumSteinerCandidates=65;refuse(*imported.model,initial,o,staged,OcctPlcRecoveryReason::InvalidRequest);
    auto limited=o;limited.maximumWork=result.workUsed-1;refuse(*imported.model,initial,limited,r,OcctPlcRecoveryReason::BudgetExceeded);
    limited=o;limited.maximumWorkingBytes=1;refuse(*imported.model,initial,limited,r,OcctPlcRecoveryReason::BudgetExceeded);
    limited=o;limited.maximumOutputBytes=1;refuse(*imported.model,initial,limited,r,OcctPlcRecoveryReason::BudgetExceeded);
    std::atomic_bool cancelled{true};limited=o;limited.cancellation=&cancelled;refuse(*imported.model,initial,limited,r,OcctPlcRecoveryReason::Cancelled);
    auto invalid=r;invalid.maximumCoreAttempts=0;refuse(*imported.model,initial,o,invalid,OcctPlcRecoveryReason::InvalidRequest);
    b=initial;d=sentinel();const auto beforeBoundary=archive(b),beforeDestination=archive(d);const auto map=parents(*imported.model);failNextAllocation=true;
    const auto failed=recoverOcctPlcVolume(b,1,*imported.model,*imported.model,1,map,o,r,d);
    require(!failed.published&&failed.reason==OcctPlcRecoveryReason::AllocationFailed,"typed recovery allocation failure");require(beforeBoundary==archive(b)&&beforeDestination==archive(d),"allocation failure published partial payload");
    require(!bnd::detail::exchangeCadGeneralPruning(false),"CAD pruning scope escaped recovery");
    {bnd::detail::ScopedCadGeneralPruning outer(true);{bnd::detail::ScopedCadGeneralPruning inner(false);}require(bnd::detail::exchangeCadGeneralPruning(true),"nested pruning scope failed to restore");}
    require(!bnd::detail::exchangeCadGeneralPruning(false),"nested pruning scope escaped");
    OcctPartitionMeshOptions partition;auto oldMesh=sentinel(),disabledMesh=sentinel();const auto old=buildOcctPartitionVolume(*imported.model,1,partition,oldMesh);OcctPartitionRecoveryReceipt receipt;
    const auto disabled=buildOcctPartitionVolume(*imported.model,1,partition,OcctPartitionRecoveryOptions{},disabledMesh,&receipt);
    require(old.published&&disabled.published&&!receipt.attempted&&old.workUsed==disabled.workUsed&&archive(oldMesh)==archive(disabledMesh),"disabled additive recovery changed legacy bytes or work");
    require(bytes(*imported.model)==source,"failed/replayed recovery mutated source");
    std::puts("Optional native fixed-PLC recovery, replay, budget, cancellation, allocation and legacy checks passed");return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
