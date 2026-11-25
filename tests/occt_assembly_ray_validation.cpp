// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctContactComponents.hpp"
#include "OcctModelInternal.hpp"
#include "OcctNativeValidation.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Trsf.hxx>
#include <atomic>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace paacs;
using namespace paacs::geo;
using namespace paacs::geo::detail;
namespace {
void check(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
std::string snapshot(const OcctModel& model) {
    std::ostringstream stream; BRepTools::Write(OcctModelAccess::shape(model), stream); return stream.str();
}
TopoDS_Compound compound(std::initializer_list<TopoDS_Shape> shapes) {
    TopoDS_Compound result; BRep_Builder builder; builder.MakeCompound(result);
    for (const auto& shape : shapes) builder.Add(result, shape);
    return result;
}
bool same(const mesh::Mesh& a, const mesh::Mesh& b) {
    if (a.nodeCount()!=b.nodeCount() || a.cellCount()!=b.cellCount() || a.connectivityCount()!=b.connectivityCount()) return false;
    for (Index i=0; i<a.nodeCount(); ++i) {
        const auto p=a.node(i),q=b.node(i);
        if (a.nodeClass(i)!=b.nodeClass(i) || std::bit_cast<std::uint64_t>(p.x)!=std::bit_cast<std::uint64_t>(q.x) ||
            std::bit_cast<std::uint64_t>(p.y)!=std::bit_cast<std::uint64_t>(q.y) ||
            std::bit_cast<std::uint64_t>(p.z)!=std::bit_cast<std::uint64_t>(q.z)) return false;
    }
    for (Index i=0; i<a.cellCount(); ++i) {
        const auto& p=a.cell(i); const auto& q=b.cell(i);
        if (p.type!=q.type || p.onEntity!=q.onEntity) return false;
        for (int j=0; j<mesh::nodesPerCell(p.type); ++j) if (a.cellNodes(i)[j]!=b.cellNodes(i)[j]) return false;
    }
    return true;
}
NativeCadMeshOptions options(bool witnesses) {
    NativeCadMeshOptions result;
    result.geometry.targetSize=.35;
    result.geometry.minimumElementQuality=0;
    result.maximumValidationWork=100'000'000;
    result.geometry.maximumMeshingWork=100'000'000;
    result.maximumPrimitiveWorkingBytes=512ULL*1024*1024;
    result.geometry.maximumAuditTemporaryBytes=128ULL*1024*1024;
    result.geometry.maximumOutputBytes=128ULL*1024*1024;
    result.geometry.maximumNodes=100'000;
    result.geometry.maximumCells=500'000;
    result.maximumNativeValidationThreads=1;
    result.useAdvancingInteriorFront=false;
    result.useTrimmedRayWitnesses=witnesses;
    // This test does not enable source substitution, boundary-first selection,
    // or a distance band. Only the existing explicit classifier toggle differs.
    return result;
}
void accounting(const OcctModel& model, const mesh::Mesh& mesh, const NativeCadMeshOutcome& result,
                bool witnesses, bool requireEveryBodyAboveThreshold) {
    if (!result.published) std::cerr << nativeCadMeshReasonName(result.reason) << ' ' << result.nativeContactDiagnostic << '\n';
    check(result.published && result.reason==NativeCadMeshReason::None && result.audit.valid &&
          result.quality.complete && result.quality.valid && result.nativeContactSourceGuard,
          "whole owned assembly retains native publication and original-source guard");
    std::map<int,std::uint64_t> tets;
    std::map<int,std::set<Index>> nodes;
    for (Index i=0;i<mesh.cellCount();++i) if(mesh.cell(i).type==mesh::CellType::Tet4) {
        const auto body=mesh.cell(i).onEntity;
        check(body.dim==3 && model.solid(body.tag),"original body tag retained");
        ++tets[body.tag]; for(int j=0;j<4;++j) nodes[body.tag].insert(mesh.cellNodes(i)[j]);
    }
    std::uint64_t totalRays=0;
    check(result.bodies.size()==model.solids().size(),"every original body has a report");
    for(const auto& body:result.bodies) {
        const auto count=tets.at(body.solidTag);
        check(count>0 && (!requireEveryBodyAboveThreshold || count>32),"fixture exercises a full centroid batch");
        check(body.nativeCentroidSourceGuard && body.nativeCentroidThreads==1,
              "strict serial final batch belongs to original source");
        check(body.nativeCentroidExactReuses==0 && body.nativeCentroidReferenceQueries==0 &&
              body.nativeCentroidGeometryPairs==0,"no copied-geometry or prior-model reuse is counted");
        check(body.nativeCentroidsChecked+body.nativeFinalRayAdmissions==count,
              "ordinary queries plus ray admissions count every final original-body centroid exactly once");
        check(witnesses || body.nativeFinalRayAdmissions==0,"ray-off original batch has no stale child ray count");
        totalRays+=body.nativeFinalRayAdmissions;
    }
    check(!witnesses || totalRays>0,"explicit ray toggle reaches original-body final validation");
    for(auto first=nodes.begin();first!=nodes.end();++first) for(auto second=nodes.begin();second!=first;++second)
        for(auto node:first->second) check(!second->second.contains(node),"independent body nodes were not welded");
}
void independentParts() {
    gp_Trsf placement; placement.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,1,0)),.31);
    placement.SetTranslationPart(gp_Vec(8,1,2));
    const auto shape=compound({BRepPrimAPI_MakeBox(1.,1.,1.).Shape(),
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(3,0,0),gp_Dir(0,0,1)),1.,1.5).Shape(),
        BRepPrimAPI_MakeCylinder(1.,1.5).Shape().Moved(TopLoc_Location(placement))});
    auto imported=OcctModelAccess::fromShape(shape); check(imported.imported(),"owned disconnected assembly imports");
    const auto& model=*imported.model; const auto original=snapshot(model);
    auto off=options(false); off.assemblyPolicy=NativeCadAssemblyPolicy::IndependentParts;
    auto on=off; on.useTrimmedRayWitnesses=true;
    mesh::Mesh ordinary,ray;
    const auto a=meshNativeCad(model,off,ordinary); accounting(model,ordinary,a,false,true);
    const auto b=meshNativeCad(model,on,ray); accounting(model,ray,b,true,true);
    check(same(ordinary,ray),"ray on/off retains exact coordinates, node/cell order, connectivity and CAD tags");
    check(a.interferences.size()==b.interferences.size(),"interference pair count unchanged");
    for(std::size_t i=0;i<a.interferences.size();++i) {
        const auto& p=a.interferences[i];const auto& q=b.interferences[i];
        check(p.firstSolidTag==q.firstSolidTag && p.secondSolidTag==q.secondSolidTag &&
              p.state==q.state && p.nativeOverlapVolume==q.nativeOverlapVolume,"interference evidence unchanged");
    }
    check(snapshot(model)==original,"ray on/off preserves original BREP bytes");
    mesh::Mesh sentinel; sentinel.addNode({9,8,7},{3,99}); const auto saved=mesh::meshHash(sentinel);
    auto bounded=on; bounded.maximumValidationWork=1; bounded.geometry.maximumMeshingWork=1;
    const auto exhausted=meshNativeCad(model,bounded,sentinel);
    check(!exhausted.published && exhausted.reason==NativeCadMeshReason::BudgetExceeded &&
          exhausted.validationWork<=1 && mesh::meshHash(sentinel)==saved,"explicit ray mode retains cumulative-budget rollback");
    std::atomic_bool cancel{true}; bounded=on; bounded.geometry.cancellation=&cancel;
    const auto cancelled=meshNativeCad(model,bounded,sentinel);
    check(!cancelled.published && cancelled.reason==NativeCadMeshReason::Cancelled &&
          mesh::meshHash(sentinel)==saved && snapshot(model)==original,"explicit ray mode retains cancellation/source rollback");
}
void contactComponents() {
    const auto source=compound({BRepPrimAPI_MakeBox(1.,1.,1.).Shape(),
        BRepPrimAPI_MakeBox(gp_Pnt(1,0,0),1.,1.,1.).Shape(),
        BRepPrimAPI_MakeBox(gp_Pnt(4,0,0),1.,1.,1.).Shape()});
    auto imported=OcctModelAccess::fromShape(source); check(imported.imported(),"owned contact-plus-isolated assembly imports");
    const auto& model=*imported.model;const auto original=snapshot(model);
    auto off=options(false);off.separateMatchedNativeFaceContacts=true;
    off.separateLowerDimensionalContacts=true;off.useNativeKernelAdmission=true;
    auto on=off;on.useTrimmedRayWitnesses=true;
    mesh::Mesh ordinary,ray;
    const auto a=buildOcctContactComponents(model,off,ordinary);accounting(model,ordinary,a,false,false);
    const auto b=buildOcctContactComponents(model,on,ray);accounting(model,ray,b,true,false);
    check(a.matchedNativeContactTriangles>0 && a.matchedNativeContactTriangles==b.matchedNativeContactTriangles &&
          a.separateContactPairs==b.separateContactPairs,"exact separate matched interfaces remain certified");
    check(same(ordinary,ray) && snapshot(model)==original,
          "contact original-body validation toggle changes no mesh payload, identity or source bytes");
}
void outsideFallback() {
    auto imported=OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(1.,1.,1.).Shape());
    check(imported.imported(),"owned outside-fallback body imports");
    const auto original=snapshot(*imported.model);
    std::vector<Vec3> points(40,Vec3{.37,.43,.59}); points.push_back({2,2,2});
    OcctNativeValidationOptions query;query.maximumThreads=1;query.useTrimmedRayWitnesses=true;
    const auto result=validateOcctNativeInteriorPoints(*imported.model,1,points,query);
    if(result.completed || result.reason!=OcctNativeValidationReason::OutsideMaterial ||
       result.firstFailedPoint!=40 || !result.trimmedRayAdmissions || !result.nativeQueries) {
        std::cerr<<"outside fixture: completed="<<result.completed<<" reason="<<static_cast<unsigned>(result.reason)
                 <<" first_failed="<<result.firstFailedPoint<<" rays="<<result.trimmedRayAdmissions
                 <<" ordinary="<<result.nativeQueries<<" work="<<result.workUsed
                 <<" working="<<result.maximumWorkingPayloadEstimate<<" snapshot="<<result.sourceSnapshotBytes
                 <<" source_guard="<<result.sourceModelUnmodified<<'\n';
    }
    check(!result.completed && result.reason==OcctNativeValidationReason::OutsideMaterial &&
          result.firstFailedPoint==40 && result.trimmedRayAdmissions>0 && result.nativeQueries>0 &&
          snapshot(*imported.model)==original,
          "outside point is rejected by ordinary fallback; independent source bytes remain unchanged");
    points.back()={0,.43,.59};
    const auto onBoundary=validateOcctNativeInteriorPoints(*imported.model,1,points,query);
    check(onBoundary.completed && onBoundary.nativeQueries>0 && onBoundary.sourceModelUnmodified,
          "ambiguous origin-on-boundary ray uses ordinary IN/ON fallback without changing native tolerance");
}
} // namespace
int main() {
    try {
        independentParts(); contactComponents(); outsideFallback();
        std::cout << "Original-body ray on/off: exact assembly payload, final centroid accounting, matched interfaces and native failure gates passed\n";
        return 0;
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
