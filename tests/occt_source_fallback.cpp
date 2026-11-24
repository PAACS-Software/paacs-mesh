// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "normal_span_knot_prism.hpp"
#include "paacs/geometry/OcctBoundary.hpp"
#include "paacs/geometry/OcctMesher.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctModelInternal.hpp"
#include "OcctStarMesher.hpp"
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools.hxx>
#include <Standard_Failure.hxx>
#include <atomic>
#include <cmath>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>

using namespace paacs;
namespace {
unsigned checks=0;
void check(bool value,const char* why){if(!value)throw std::runtime_error(why);++checks;}
bool named(const char* value,std::string_view expected){return value&&std::string_view(value)==expected;}
std::string snapshot(const geo::OcctModel& model){
    std::ostringstream text;BRepTools::Write(geo::detail::OcctModelAccess::shape(model),text);return text.str();
}
mesh::Mesh sentinel(){mesh::Mesh mesh;mesh.addNode({-3,-4,-5},{3,999});return mesh;}

geo::NativeCadMeshOptions options(){
    auto o=geo::normalNativeCadMeshOptions();
    o.geometry.targetSize=.3;o.geometry.threads=1;o.maximumNativeValidationThreads=1;
    o.geometry.maximumMeshingWork=100'000'000;o.maximumValidationWork=100'000'000;
    o.geometry.maximumNodes=500'000;o.geometry.maximumCells=3'000'000;
    o.geometry.maximumOutputBytes=256ULL*1024*1024;o.geometry.maximumAuditTemporaryBytes=512ULL*1024*1024;
    o.maximumPrimitiveWorkingBytes=1024ULL*1024*1024;
    o.useAdvancingInteriorFront=true;o.recoverBoundaryBeforeInteriorSites=true;
    return o;
}

void checkWitness(const geo::NativeCadSourceShellReport& shell){
    check(shell.count>0&&shell.count<=shell.attempts.size(),"source attempt count fits fixed allowance");
    const auto& first=shell.attempts[0];
    check(named(first.provider,"regular_chart")&&!first.published&&
          named(first.surfaceDiagnostic,"nonconforming_edge"),"original regular-source refusal retained");
    check(named(first.witness.stage,"shared_edge_sampling")&&
          named(first.witness.detail,"curve_parameter_resolution")&&first.witness.geometricRefusal&&
          first.witness.edgeTag>0&&!first.constructionBudgetLimited,"original sampling witness retains geometric cause and edge identity");
}

void checkMesh(const geo::OcctModel& model,const geo::NativeCadMeshOutcome& result,const mesh::Mesh& mesh){
    check(result.published&&!result.destinationUnchanged&&result.reason==geo::NativeCadMeshReason::None,
          "enabled source substitution publishes valid candidate");
    check(result.audit.complete&&result.audit.valid&&result.quality.valid,"native topology and cell-quality gates pass");
    check(mesh::auditLinearVolumeBoundary(mesh).valid,"independent exact boundary incidence audit passes");
    check(mesh.countCells(mesh::CellType::Tet4)>0&&result.bodies.size()==1,"published owned solid contains tetrahedra");
    const auto& body=result.bodies.front();
    checkWitness(body.sourceShell);
    check(body.sourceShell.nativeFallbackRequested&&body.sourceShell.nativeFallbackAttempted&&body.sourceShell.accepted,
          "published mesh actually used requested Native source");
    check(body.sourceShell.count>=2,"source fallback consumes another shared source attempt");
    const auto& source=body.sourceShell.attempts[1];
    check(named(source.provider,"native_tessellation")&&source.published&&source.sourceModelUnmodified&&
          source.volumeCriterionMet&&!source.globalSizingApplied&&!source.witness.geometricRefusal,
          "Native provider has separate truthful receipt and retains volume/source gates");
    check(body.nativeCentroidSourceGuard&&body.nativeCentroidsChecked>0,
          "final Native centroid audit guards original source");
    check(body.relativeVolumeError<=.01&&std::abs(result.quality.totalVolume-1.)<.01,
          "published volume agrees with owned unit prism");
    std::set<int> faces,bodies,expectedFaces,expectedBodies;
    for(Index i=0;i<mesh.cellCount();++i){
        const auto& cell=mesh.cell(i);
        if(cell.onEntity.dim==2)faces.insert(cell.onEntity.tag);
        if(cell.onEntity.dim==3)bodies.insert(cell.onEntity.tag);
    }
    for(const auto* face:model.faces())expectedFaces.insert(face->tag());
    for(const auto* solid:model.solids())expectedBodies.insert(solid->tag());
    check(faces==expectedFaces&&bodies==expectedBodies,"all Native face and body identities survive source substitution");
}

void exercise(const geo::OcctModel& model){
    const auto before=snapshot(model);
    geo::OcctBoundaryOptions sourceOptions;sourceOptions.linearDeflection=.002;sourceOptions.maximumChordError=.005;
    sourceOptions.maximumWork=50'000'000;sourceOptions.construction.resourcePolicy=geo::CadResourcePolicy::Normal;
    mesh::Mesh boundary;
    const auto native=geo::buildOcctBoundary(model,1,sourceOptions,boundary);
    check(native.published&&native.sourceModelUnmodified&&native.relativeVolumeError<.01,
          "owned Native boundary independently qualifies before volume fallback test");
    const auto certificate=geo::detail::certifyOcctBoundarySurface(boundary,geo::detail::OcctStarOptions{});
    check(certificate.certified&&certificate.embedded&&certificate.closedOrientedVertexManifold,
          "owned Native source is an embedded closed oriented manifold");
    auto disabledOptions=options();
    auto disabledMesh=sentinel();const auto disabledHash=mesh::meshHash(disabledMesh);
    const auto disabled=geo::meshNativeCad(model,disabledOptions,disabledMesh);
    check(!disabled.published&&disabled.destinationUnchanged&&disabled.reason==geo::NativeCadMeshReason::GeometryError&&
          mesh::meshHash(disabledMesh)==disabledHash,"disabled fallback refuses original source transactionally");
    check(disabled.bodies.size()==1,"disabled attempt retains body diagnostic");
    checkWitness(disabled.bodies.front().sourceShell);
    check(!disabled.bodies.front().sourceShell.nativeFallbackAttempted,"disabled request never substitutes source");
    auto enabledOptions=options();enabledOptions.allowNativeSourceFallback=true;
    mesh::Mesh enabledMesh;
    const auto enabled=geo::meshNativeCad(model,enabledOptions,enabledMesh);
    checkMesh(model,enabled,enabledMesh);
    auto capped=enabledOptions;capped.construction.interiorCandidates=1;capped.construction.interiorSites=1;
    auto cappedMesh=sentinel();const auto cappedHash=mesh::meshHash(cappedMesh);
    const auto cappedResult=geo::meshNativeCad(model,capped,cappedMesh);
    check(!cappedResult.bodies.empty(),"explicit insertion caps retain body receipt");
    check(cappedResult.bodies.front().interiorFrontCandidates<=1&&
          cappedResult.bodies.front().interiorFrontNodes<=1,"boundary-first constructor respects explicit candidate/site upper bounds");
    if(cappedResult.published)checkMesh(model,cappedResult,cappedMesh);
    else check(cappedResult.destinationUnchanged&&mesh::meshHash(cappedMesh)==cappedHash,
               "capped refusal retains destination");
    auto surfaceCapped=enabledOptions;surfaceCapped.construction.surfacePointsPerFace=1;
    auto surfaceMesh=sentinel();const auto surfaceHash=mesh::meshHash(surfaceMesh);
    const auto surfaceResult=geo::meshNativeCad(model,surfaceCapped,surfaceMesh);
    check(!surfaceResult.published&&surfaceResult.destinationUnchanged&&mesh::meshHash(surfaceMesh)==surfaceHash,
          "explicit surface cap cannot be evaded by source substitution");
    for(const auto& body:surfaceResult.bodies)check(!body.sourceShell.nativeFallbackAttempted,
          "explicit surface cap never enters Native source provider");
    std::atomic_bool cancel=true;auto cancelled=enabledOptions;cancelled.geometry.cancellation=&cancel;
    auto cancelledMesh=sentinel();const auto cancelHash=mesh::meshHash(cancelledMesh);
    const auto cancelResult=geo::meshNativeCad(model,cancelled,cancelledMesh);
    check(!cancelResult.published&&cancelResult.reason==geo::NativeCadMeshReason::Cancelled&&
          cancelResult.destinationUnchanged&&mesh::meshHash(cancelledMesh)==cancelHash,"cancelled fallback request publishes nothing");
    check(before==snapshot(model),"all successful/refused/capped/cancelled attempts preserve original CAD bytes");
}
}

int main(){
    try{
        const auto shape=test::normalSpanKnotPrism(true);
        check(BRepCheck_Analyzer(shape).IsValid(),"owned unit prism is valid Native BRep before import");
        auto imported=geo::detail::OcctModelAccess::fromShape(shape);
        check(imported.imported(),"exact owned fixture imports without healing or file normalization");
        exercise(*imported.model);
        std::cout<<"OCCT source fallback checks="<<checks<<" PASS\n";return 0;
    }catch(const Standard_Failure& e){std::cerr<<"OCCT source fallback kernel: "<<e.what()<<'\n';}
     catch(const std::exception& e){std::cerr<<"OCCT source fallback: "<<e.what()<<'\n';}
    return 1;
}
