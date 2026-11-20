// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctStarMesher.hpp"
#include "OcctModelInternal.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <stdexcept>
#include <sstream>
#include <string_view>

static std::atomic_bool failNextAllocation{false};
static std::atomic_size_t failAllocationSize{0};
static std::atomic_size_t allocationSkips{0},countAllocationSize{0},matchedAllocations{0};
void* operator new(std::size_t n) {
    if(failNextAllocation.exchange(false))throw std::bad_alloc{};
    if(const auto target=countAllocationSize.load();target!=0&&target==n)++matchedAllocations;
    if(const auto target=failAllocationSize.load();target!=0&&target==n){
        if(allocationSkips.load()!=0)--allocationSkips;
        else{failAllocationSize=0;throw std::bad_alloc{};}
    }
    if(void* p=std::malloc(n?n:1))return p;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void* p) noexcept {std::free(p);}
void operator delete[](void* p) noexcept {std::free(p);}
void operator delete(void* p,std::size_t) noexcept {std::free(p);}
void operator delete[](void* p,std::size_t) noexcept {std::free(p);}

using namespace paacs;
using namespace paacs::geo::detail;
namespace {
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
mesh::Mesh sentinel(){mesh::Mesh m;m.addNode({91,92,93},{0,94});return m;}
void unchanged(const mesh::Mesh& m){require(m.nodeCount()==1&&m.cellCount()==0&&m.node(0).x==91&&m.nodeClass(0)==mesh::Classification{0,94},"destination changed on refusal");}
bool samePayload(const mesh::Mesh& a,const mesh::Mesh& b){
    if(a.nodeCount()!=b.nodeCount()||a.cellCount()!=b.cellCount())return false;
    for(Index i=0;i<a.nodeCount();++i){const auto x=a.node(i),y=b.node(i);if(x.x!=y.x||x.y!=y.y||x.z!=y.z||a.nodeClass(i)!=b.nodeClass(i))return false;}
    for(Index i=0;i<a.cellCount();++i){const auto x=a.cell(i),y=b.cell(i);if(x.type!=y.type||x.onEntity!=y.onEntity)return false;for(int j=0;j<mesh::nodesPerCell(x.type);++j)if(a.cellNodes(i)[j]!=b.cellNodes(i)[j])return false;}
    return true;
}
mesh::Mesh cube(){
    mesh::Mesh m;
    for(Vec3 p:std::array<Vec3,8>{{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}}})
        m.addNode(p,{0,static_cast<int>(m.nodeCount()+101)});
    const std::array<std::array<Index,3>,12> triangles{{{0,3,2},{0,2,1},{4,5,6},{4,6,7},
        {0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}}};
    for(auto t:triangles)m.addCell(mesh::CellType::Tri3,t,{2,static_cast<int>(m.cellCount()+201)});
    return m;
}
mesh::Mesh prism(std::span<const Vec3> polygon,std::span<const std::array<Index,3>> cap){
    mesh::Mesh m;const auto n=static_cast<Index>(polygon.size());
    for(double z:{0.,1.})for(auto p:polygon)m.addNode({p.x,p.y,z},{1,static_cast<int>(m.nodeCount()+1)});
    for(auto t:cap){m.addCell(mesh::CellType::Tri3,{t[0],t[2],t[1]},{2,1});m.addCell(mesh::CellType::Tri3,{t[0]+n,t[1]+n,t[2]+n},{2,2});}
    for(Index i=0;i<n;++i){const Index j=(i+1)%n;m.addCell(mesh::CellType::Tri3,{i,j,j+n},{2,static_cast<int>(i+3)});m.addCell(mesh::CellType::Tri3,{i,j+n,i+n},{2,static_cast<int>(i+3)});}
    return m;
}
void refusal(const mesh::Mesh& b,const OcctStarOptions& opts,OcctStarReason reason){
    auto d=sentinel();const auto o=buildOcctStarVolume(b,51,opts,d);
    require(!o.published&&o.destinationUnchanged&&o.reason==reason,o.diagnostic);unchanged(d);
}
void accepted(const mesh::Mesh& b,double volume){
    auto d=sentinel();const auto o=buildOcctStarVolume(b,51,{},d);
    require(o.published&&!o.destinationUnchanged&&o.reason==OcctStarReason::None,o.diagnostic);
    require(o.embedded&&o.kernelFeasible,"missing star certificate");
    require(d.nodeCount()==b.nodeCount()+1&&d.cellCount()==2*b.cellCount(),"fan counts");
    require(d.nodeClass(b.nodeCount())==mesh::Classification{3,51},"body node identity");
    for(Index i=0;i<b.nodeCount();++i)require(d.node(i).x==b.node(i).x&&d.node(i).y==b.node(i).y&&d.node(i).z==b.node(i).z&&d.nodeClass(i)==b.nodeClass(i),"source node changed");
    for(Index i=0;i<b.cellCount();++i){const auto j=i+b.cellCount();require(d.cell(j).type==mesh::CellType::Tri3&&d.cell(j).onEntity==b.cell(i).onEntity,"source face changed");for(unsigned k=0;k<3;++k)require(d.cellNodes(j)[k]==b.cellNodes(i)[k],"source connectivity changed");}
    const auto audit=mesh::auditLinearVolumeBoundary(d);const auto q=algo::analyzeMixedMeshQuality(d);
    require(audit.valid&&q.valid&&std::abs(q.totalVolume-volume)<1e-12,"final fan audit");
}
}
int main(){try{
    const auto box=cube();accepted(box,1);

    // This candidate-selection fixture intentionally uses a smaller native
    // admission domain: boundary ownership remains a separate caller gate.
    // The former cube center produces left-side Tet centroids outside it,
    // while another strict PLC kernel makes every fan centroid native IN.
    const auto admissionShape=BRepPrimAPI_MakeBox(gp_Pnt(.175,0,0),gp_Pnt(1,1,1)).Shape();
    auto admission=OcctModelAccess::fromShape(admissionShape);
    require(admission.imported()&&admission.solids==1,"native admission fixture");
    std::ostringstream nativeBefore;BRepTools::Write(OcctModelAccess::shape(*admission.model),nativeBefore);
    auto ordinary=sentinel();const auto ordinaryResult=buildOcctStarVolume(box,51,{},ordinary);
    require(ordinaryResult.published&&ordinaryResult.nativeCandidateQueries==0,"default fan admission changed");
    bool ordinaryOutside=false;
    for(Index i=0;i<ordinary.cellCount();++i)if(ordinary.cell(i).type==mesh::CellType::Tet4){
        const auto* t=ordinary.cellNodes(i);const auto p=(ordinary.node(t[0])+ordinary.node(t[1])+ordinary.node(t[2])+ordinary.node(t[3]))*.25;
        ordinaryOutside|=!admission.model->solid(1)->isInside(p);
    }
    require(ordinaryOutside,"fixture must reject the former fan point");
    OcctStarOptions nativeOptions;nativeOptions.maximumProposalSweeps=0;
    nativeOptions.nativeAdmissionModel=admission.model.get();nativeOptions.nativeAdmissionBodyTag=1;
    auto admitted=sentinel();const auto admittedResult=buildOcctStarVolume(box,51,nativeOptions,admitted);
    require(admittedResult.published&&admittedResult.kernelFeasible&&admittedResult.nativeCandidateQueries>box.cellCount(),admittedResult.diagnostic);
    require(admittedResult.kernel.x>=.7&&mesh::auditLinearVolumeBoundary(admitted).valid,"alternate native kernel and exact audit");
    require(admission.model->solid(1)->isInside(admittedResult.kernel),"accepted kernel outside native body");
    for(Index i=0;i<admitted.cellCount();++i)if(admitted.cell(i).type==mesh::CellType::Tet4){
        const auto* t=admitted.cellNodes(i);const auto p=(admitted.node(t[0])+admitted.node(t[1])+admitted.node(t[2])+admitted.node(t[3]))*.25;
        require(admission.model->solid(1)->isInside(p),"accepted native fan centroid outside");
    }
    for(Index i=0;i<box.nodeCount();++i)require(admitted.node(i).x==box.node(i).x&&admitted.node(i).y==box.node(i).y&&admitted.node(i).z==box.node(i).z&&admitted.nodeClass(i)==box.nodeClass(i),"native admission moved boundary");
    for(Index i=0;i<box.cellCount();++i){const auto j=i+box.cellCount();require(admitted.cell(j).onEntity==box.cell(i).onEntity,"native admission changed face ownership");for(unsigned k=0;k<3;++k)require(admitted.cellNodes(j)[k]==box.cellNodes(i)[k],"native admission changed directed facet");}
    nativeOptions.maximumWork=admittedResult.workUsed;auto nativeReplay=sentinel();
    const auto nativeRepeated=buildOcctStarVolume(box,51,nativeOptions,nativeReplay);
    require(nativeRepeated.published&&samePayload(admitted,nativeReplay)&&nativeRepeated.workUsed==admittedResult.workUsed&&nativeRepeated.nativeCandidateQueries==admittedResult.nativeCandidateQueries,"native admission deterministic exact-budget replay");
    --nativeOptions.maximumWork;refusal(box,nativeOptions,OcctStarReason::WorkBudget);
    nativeOptions.maximumWork=5'000'000;std::atomic_bool nativeCancelled{true};nativeOptions.cancellation=&nativeCancelled;
    refusal(box,nativeOptions,OcctStarReason::Cancelled);nativeOptions.cancellation=nullptr;
    auto nativeAllocationDestination=sentinel();failNextAllocation=true;
    const auto nativeAllocation=buildOcctStarVolume(box,51,nativeOptions,nativeAllocationDestination);
    require(nativeAllocation.reason==OcctStarReason::AllocationFailed&&!nativeAllocation.published,"native admission allocation rollback");unchanged(nativeAllocationDestination);
    nativeOptions.nativeAdmissionBodyTag=999;refusal(box,nativeOptions,OcctStarReason::InvalidBoundary);
    std::ostringstream nativeAfter;BRepTools::Write(OcctModelAccess::shape(*admission.model),nativeAfter);
    require(nativeBefore.str()==nativeAfter.str(),"native admission changed source CAD");
    auto ordinaryReplay=sentinel();const auto ordinaryRepeated=buildOcctStarVolume(box,51,{},ordinaryReplay);
    require(ordinaryRepeated.published&&samePayload(ordinary,ordinaryReplay)&&ordinaryRepeated.workUsed==ordinaryResult.workUsed,"native opt-in altered default payload/work");

    // This thin concave prism has a strict kernel missed by the former 16-sweep
    // proposal search. Its fixed boundary is retained verbatim by the fan.
    const std::array<Vec3,6> thinPoints{{{0,0,0},{4,0,0},{1004,.001,0},{1001,.001,0},{4001,.004,0},{4000,.004,0}}};
    const std::array<std::array<Index,3>,4> thinCap{{{0,1,3},{1,2,3},{0,3,5},{3,4,5}}};
    const auto thin=prism(thinPoints,thinCap);auto thinOptions=OcctStarOptions{};thinOptions.maximumProposalSweeps=16;
    auto thinDestination=sentinel();const auto thinResult=buildOcctStarVolume(thin,51,thinOptions,thinDestination);
    require(thinResult.published&&thinResult.embedded&&thinResult.kernelFeasible,"thin strict kernel proposal");
    require(thinResult.minimumQuality>0&&std::abs(thinResult.boundaryVolume-.007)<1e-12,"thin fan volume/quality");
    require(mesh::auditLinearVolumeBoundary(thinDestination).valid,"thin complete boundary incidence");
    for(Index i=0;i<thin.nodeCount();++i)require(thinDestination.node(i).x==thin.node(i).x&&thinDestination.node(i).y==thin.node(i).y&&thinDestination.node(i).z==thin.node(i).z&&thinDestination.nodeClass(i)==thin.nodeClass(i),"thin boundary node changed");
    for(Index i=0;i<thin.cellCount();++i){const auto j=i+thin.cellCount();require(thinDestination.cell(j).onEntity==thin.cell(i).onEntity,"thin boundary face identity changed");for(unsigned k=0;k<3;++k)require(thinDestination.cellNodes(j)[k]==thin.cellNodes(i)[k],"thin directed boundary changed");}
    thinOptions.maximumWork=thinResult.workUsed-1;refusal(thin,thinOptions,OcctStarReason::WorkBudget);
    thinOptions.maximumWork=thinResult.workUsed;auto exactAdmission=sentinel();const auto replay=buildOcctStarVolume(thin,51,thinOptions,exactAdmission);
    require(replay.published&&replay.workUsed==thinResult.workUsed,"thin exact work admission");
    require(replay.kernel.x==thinResult.kernel.x&&replay.kernel.y==thinResult.kernel.y&&replay.kernel.z==thinResult.kernel.z,"thin deterministic proposal");
    thinOptions={};thinOptions.maximumProposalSweeps=16;
    // Admit preparation but deny the additional bounded proposal-plane payload.
    thinOptions.maximumWorkingBytes=static_cast<std::uint64_t>(thin.nodeCount()+thin.cellCount())*1024;
    auto planeLimitedDestination=sentinel();const auto planeLimited=buildOcctStarVolume(thin,51,thinOptions,planeLimitedDestination);
    require(!planeLimited.published&&planeLimited.reason==OcctStarReason::WorkBudget&&std::string_view(planeLimited.diagnostic)=="proposal-working-byte-budget","thin plane working admission");unchanged(planeLimitedDestination);
    thinOptions={};thinOptions.maximumProposalSweeps=16;thinOptions.maximumWork=planeLimited.workUsed+thin.cellCount();
    auto sortLimitedDestination=sentinel();const auto sortLimited=buildOcctStarVolume(thin,51,thinOptions,sortLimitedDestination);
    require(!sortLimited.published&&sortLimited.embedded&&!sortLimited.kernelFeasible&&sortLimited.reason==OcctStarReason::WorkBudget&&sortLimited.workUsed==thinOptions.maximumWork,"thin proposal sort work admission");unchanged(sortLimitedDestination);
    thinOptions={};thinOptions.maximumProposalSweeps=16;
    struct PlaneAllocationLayout {Vec3 normal;double offset;std::uint64_t order;};
    auto planeAllocationDestination=sentinel();failAllocationSize=(thin.cellCount()+6)*sizeof(PlaneAllocationLayout);
    const auto planeAllocation=buildOcctStarVolume(thin,51,thinOptions,planeAllocationDestination);const auto failedSize=failAllocationSize.exchange(0);
    require(failedSize==0&&!planeAllocation.published&&planeAllocation.embedded&&!planeAllocation.kernelFeasible&&planeAllocation.reason==OcctStarReason::AllocationFailed&&planeAllocation.workUsed==planeLimited.workUsed,"thin proposal plane allocation rollback");unchanged(planeAllocationDestination);
    thinOptions={};thinOptions.maximumBoundaryTriangles=thin.cellCount()-1;refusal(thin,thinOptions,OcctStarReason::WorkBudget);

    // A determinant-one shear exposes cancellation between nearly parallel
    // projected planes. Numerical retries still require every exact side test.
    auto parallel=thin;
    for(Index i=0;i<parallel.nodeCount();++i){auto p=parallel.node(i);p.x+=p.z/16;p.z+=3*p.y/8;parallel.setNode(i,p);}
    OcctStarOptions parallelOptions;parallelOptions.maximumProposalSweeps=16;
    const auto planeArrayBytes=(parallel.cellCount()+6)*sizeof(PlaneAllocationLayout);
    matchedAllocations=0;countAllocationSize=planeArrayBytes;
    auto parallelDestination=sentinel();const auto parallelResult=buildOcctStarVolume(parallel,51,parallelOptions,parallelDestination);
    countAllocationSize=0;const auto planeAllocations=matchedAllocations.load();
    require(parallelResult.published&&parallelResult.kernelFeasible&&parallelResult.embedded&&parallelResult.minimumQuality>0,"parallel strict kernel retry");
    require(planeAllocations>0&&planeAllocations<=6,"bounded plane proposal allocations");
    const auto parallelAudit=mesh::auditLinearVolumeBoundary(parallelDestination);const auto parallelQuality=algo::analyzeMixedMeshQuality(parallelDestination);
    require(parallelAudit.valid&&parallelQuality.valid&&std::abs(parallelQuality.totalVolume-.007)<1e-12,"parallel complete fan volume audit");
    parallelOptions.maximumWork=parallelResult.workUsed;auto parallelReplay=sentinel();const auto parallelRepeated=buildOcctStarVolume(parallel,51,parallelOptions,parallelReplay);
    require(parallelRepeated.published&&parallelRepeated.workUsed==parallelResult.workUsed&&samePayload(parallelReplay,parallelDestination),"parallel exact work replay");
    parallelOptions.maximumWork=parallelResult.workUsed-1;refusal(parallel,parallelOptions,OcctStarReason::WorkBudget);
    parallelOptions={};parallelOptions.maximumProposalSweeps=16;
    // Inject failure at every observed proposal-plane allocation, including
    // the retry after earlier numerical attempts fail, without assuming which
    // ordinary attempts fail on a different floating-point implementation.
    std::uint64_t lastAllocationWork=0;
    for(std::size_t skip=0;skip<planeAllocations;++skip){
        auto failedDestination=sentinel();allocationSkips=skip;failAllocationSize=planeArrayBytes;
        const auto failedProposal=buildOcctStarVolume(parallel,51,parallelOptions,failedDestination);
        const auto pending=failAllocationSize.exchange(0);allocationSkips=0;
        require(pending==0&&!failedProposal.published&&failedProposal.embedded&&!failedProposal.kernelFeasible&&failedProposal.reason==OcctStarReason::AllocationFailed,"parallel proposal allocation rollback");
        unchanged(failedDestination);lastAllocationWork=failedProposal.workUsed;
    }
    parallelOptions.maximumWork=lastAllocationWork+parallel.cellCount();auto retryLimitedDestination=sentinel();const auto retryLimited=buildOcctStarVolume(parallel,51,parallelOptions,retryLimitedDestination);
    require(!retryLimited.published&&retryLimited.embedded&&!retryLimited.kernelFeasible&&retryLimited.reason==OcctStarReason::WorkBudget&&retryLimited.workUsed==parallelOptions.maximumWork,"parallel retry sorting work limit");unchanged(retryLimitedDestination);
    std::atomic_bool cancelParallel{true};parallelOptions={};parallelOptions.cancellation=&cancelParallel;refusal(parallel,parallelOptions,OcctStarReason::Cancelled);
    const auto boxCertificate=certifyOcctBoundarySurface(box,{});
    require(boxCertificate.certified&&boxCertificate.closedOrientedVertexManifold&&boxCertificate.embedded&&boxCertificate.connectedShells==1&&boxCertificate.eulerCharacteristic==2&&boxCertificate.exactBoundaryPairs==54,"cube surface certificate");
    auto cavity=box;
    for(Index i=0;i<box.nodeCount();++i)cavity.addNode(box.node(i)*.5+Vec3{.25,.25,.25},{0,static_cast<int>(301+i)});
    for(Index i=0;i<box.cellCount();++i){const auto* t=box.cellNodes(i);cavity.addCell(mesh::CellType::Tri3,{t[0]+8,t[2]+8,t[1]+8},{2,static_cast<int>(401+i)});}
    const auto cavityCertificate=certifyOcctBoundarySurface(cavity,{});
    require(cavityCertificate.certified&&cavityCertificate.connectedShells==2&&cavityCertificate.eulerCharacteristic==4&&std::abs(cavityCertificate.signedVolume-.875)<1e-12&&cavityCertificate.exactBoundaryPairs==108,"void shell surface certificate");
    refusal(cavity,{},OcctStarReason::UnsupportedTopology);
    mesh::Mesh ring;
    const std::array<Vec3,8> ringPoints{{{0,0,0},{3,0,0},{3,3,0},{0,3,0},{1,1,0},{2,1,0},{2,2,0},{1,2,0}}};
    for(double z:{0.,1.})for(auto p:ringPoints)ring.addNode({p.x,p.y,z},{0,static_cast<int>(ring.nodeCount()+1)});
    for(auto t:std::array<std::array<Index,3>,8>{{{0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}}}){ring.addCell(mesh::CellType::Tri3,{t[0],t[2],t[1]},{2,1});ring.addCell(mesh::CellType::Tri3,{t[0]+8,t[1]+8,t[2]+8},{2,2});}
    for(Index i=0;i<4;++i){const auto j=(i+1)%4;ring.addCell(mesh::CellType::Tri3,{i,j,j+8},{2,3});ring.addCell(mesh::CellType::Tri3,{i,j+8,i+8},{2,3});const auto a=i+4,b=j+4;ring.addCell(mesh::CellType::Tri3,{a,b+8,b},{2,4});ring.addCell(mesh::CellType::Tri3,{a,a+8,b+8},{2,4});}
    const auto ringCertificate=certifyOcctBoundarySurface(ring,{});
    require(ringCertificate.certified&&ringCertificate.connectedShells==1&&ringCertificate.eulerCharacteristic==0&&std::abs(ringCertificate.signedVolume-8)<1e-12,"genus-one boundary certificate");
    refusal(ring,{},OcctStarReason::UnsupportedTopology);
    const std::array<Vec3,6> l{{{0,0,0},{2,0,0},{2,1,0},{1,1,0},{1,2,0},{0,2,0}}};
    const std::array<std::array<Index,3>,4> lc{{{0,1,3},{1,2,3},{0,3,5},{3,4,5}}};
    accepted(prism(l,lc),3);
    const std::array<Vec3,8> u{{{0,0,0},{3,0,0},{3,3,0},{2,3,0},{2,1,0},{1,1,0},{1,3,0},{0,3,0}}};
    const std::array<std::array<Index,3>,6> uc{{{0,1,4},{1,2,3},{1,3,4},{0,4,5},{0,5,7},{5,6,7}}};
    refusal(prism(u,uc),{},OcctStarReason::KernelNotFound);
    auto opts=OcctStarOptions{};opts.maximumWork=0;refusal(box,opts,OcctStarReason::WorkBudget);
    opts={};opts.maximumNodes=8;refusal(box,opts,OcctStarReason::OutputBudget);
    opts={};opts.maximumCells=23;refusal(box,opts,OcctStarReason::OutputBudget);
    opts={};opts.maximumOutputBytes=1;refusal(box,opts,OcctStarReason::OutputBudget);
    opts={};opts.maximumWorkingBytes=1;refusal(box,opts,OcctStarReason::WorkBudget);
    opts={};opts.auditLimits.maximumFaceUses=0;refusal(box,opts,OcctStarReason::WorkBudget);
    opts={};opts.auditLimits.maximumTemporaryBytes=0;refusal(box,opts,OcctStarReason::WorkBudget);
    std::atomic_bool cancelled{true};opts={};opts.cancellation=&cancelled;refusal(box,opts,OcctStarReason::Cancelled);
    auto open=box;open=mesh::Mesh{};for(Index i=0;i<box.nodeCount();++i)open.addNode(box.node(i),box.nodeClass(i));for(Index i=0;i<box.cellCount()-1;++i)open.addCell(box.cell(i).type,{box.cellNodes(i),3},box.cell(i).onEntity);refusal(open,{},OcctStarReason::UnsupportedTopology);
    auto crossing=box;crossing.setNode(6,{.2,.2,-.2});refusal(crossing,{},OcctStarReason::BoundaryContact);
    auto d=sentinel();failNextAllocation=true;const auto failed=buildOcctStarVolume(box,51,{},d);require(failed.reason==OcctStarReason::AllocationFailed&&!failed.published,"allocation refusal");unchanged(d);
    std::puts("Exact star volume regression checks passed");return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
