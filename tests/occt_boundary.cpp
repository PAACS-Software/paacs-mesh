// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctBoundary.hpp"
#include "paacs/algo/Delaunay2D.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "../src/geometry/OcctBoundaryInternal.hpp"
#include "../src/geometry/OcctCurveMetric.hpp"
#include "../src/geometry/OcctModelInternal.hpp"
#include "../src/geometry/OcctStarMesher.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepLib.hxx>
#include <Geom_BezierSurface.hxx>
#include <NCollection_Array2.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <TopExp.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <NCollection_List.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <algorithm>
#include <bit>
#include <filesystem>
#include <iostream>
#include <map>
#include <sstream>
#include <set>

using namespace paacs;
namespace {
void check(bool value,const char* message) { if(!value) throw std::runtime_error(message); }
std::string snapshot(const geo::OcctModel& model) {
    std::ostringstream stream; BRepTools::Write(geo::detail::OcctModelAccess::shape(model),stream); return stream.str();
}
geo::NativeCadImportOutcome load(const TopoDS_Shape& shape,const std::filesystem::path& path) {
    check(BRepTools::Write(shape,path.string().c_str()),"fixture BREP write");
    geo::NativeCadImportOptions options; options.format=geo::CadFileFormat::Brep;
    auto imported=geo::importOcctCad(path,options); check(imported.imported(),"fixture import"); return imported;
}
bool same(Vec3 a,Vec3 b) {
    return std::bit_cast<std::uint64_t>(a.x)==std::bit_cast<std::uint64_t>(b.x)&&
           std::bit_cast<std::uint64_t>(a.y)==std::bit_cast<std::uint64_t>(b.y)&&
           std::bit_cast<std::uint64_t>(a.z)==std::bit_cast<std::uint64_t>(b.z);
}
int constrainedChartRegions() {
    // Each strip consists of two triangles between the outer and inner wire.
    // The final two triangles fill the inner loop and must be marked exterior.
    const std::vector<algo::Tri> triangles{
        {{0,1,5}},{{0,5,4}},{{1,2,6}},{{1,6,5}},
        {{2,3,7}},{{2,7,6}},{{3,0,4}},{{3,4,7}},
        {{4,5,6}},{{4,6,7}}
    };
    std::set<std::array<Index,2>> wires{
        {0,1},{1,2},{2,3},{0,3},{4,5},{5,6},{6,7},{4,7}
    };
    const auto charge=[](std::uint64_t){};
    const auto labels=geo::detail::labelOcctChartRegions(triangles,wires,charge);
    check(labels&&*labels==std::vector<int>{1,1,1,1,1,1,1,1,0,0},
          "chart parity respects outer and inner constraints");
    wires.insert({0,5});
    check(!geo::detail::labelOcctChartRegions(triangles,wires,charge),
          "chart parity refuses contradictory barriers");
    wires.erase({0,5});
    wires.insert({0,2});
    check(!geo::detail::labelOcctChartRegions(triangles,wires,charge),
          "chart parity refuses missing recovered constraints");
    return 3;
}
int regularPrismSurface(const std::filesystem::path& folder) {
    using namespace geo::detail;
    constexpr double side=6.,height=1.8,target=1.01;
    const auto apexY=side*std::sqrt(3.)/2.;
    const auto trim=BRepBuilderAPI_MakePolygon(
        gp_Pnt(0,0,0),gp_Pnt(side,0,0),gp_Pnt(side/2.,apexY,0),true).Wire();
    const auto cap=BRepBuilderAPI_MakeFace(trim).Face();
    const auto prism=BRepPrimAPI_MakePrism(cap,gp_Vec(0,0,height)).Shape();
    const auto imported=load(prism,folder/"regular-prism.brep");
    const auto source=snapshot(*imported.model);
    geo::OcctBoundaryOptions options;
    options.maximumChordError=.01;
    constexpr std::uint64_t storage=64ULL*1024*1024;
    OcctBoundaryAliasReceipt receipt;
    mesh::Mesh regular;
    const auto result=buildSpacedOcctBoundary(*imported.model,1,options,target,storage,regular,receipt,true);
    check(result.published&&result.sourceModelUnmodified&&receipt.rebuiltCharts==5&&
          result.relativeVolumeError<1e-8&&certifyOcctBoundarySurface(regular,OcctStarOptions{}).certified,
          "regular prism triangulation retains all five CAD faces and encloses its volume");

    std::vector<double> capQuality;
    double capArea=0,maximumCapEdge=0;
    for(Index i=0;i<regular.cellCount();++i) {
        const auto* nodes=regular.cellNodes(i);
        const auto a=regular.node(nodes[0]),b=regular.node(nodes[1]),c=regular.node(nodes[2]);
        if(std::max({std::abs(a.z),std::abs(b.z),std::abs(c.z)})>1e-8)continue;
        const auto ab=(b-a).norm2(),bc=(c-b).norm2(),ca=(a-c).norm2();
        const auto twiceArea=(b-a).cross(c-a).norm();
        capQuality.push_back(2*std::sqrt(3.)*twiceArea/(ab+bc+ca));
        capArea+=twiceArea/2.;
        maximumCapEdge=std::max(maximumCapEdge,std::sqrt(std::max({ab,bc,ca})));
    }
    std::ranges::sort(capQuality);
    check(capQuality.size()>=20&&std::abs(capArea-side*apexY/2.)<1e-8,
          "regular prism cap resolves the requested spacing over the entire trim");
    const auto p10=capQuality[capQuality.size()/10];
    const auto median=capQuality[capQuality.size()/2];
    std::cout<<"regular prism cap: triangles="<<capQuality.size()<<" mean ratio p10="<<p10
             <<" median="<<median<<" maximum edge="<<maximumCapEdge<<'\n';
    check(p10>=.95&&maximumCapEdge<=1.18*target,
          "regular prism cap uses consistently well-shaped triangles near the target size");
    const auto hash=mesh::meshHash(regular);
    const auto replay=buildSpacedOcctBoundary(*imported.model,1,options,target,storage,regular,receipt,true);
    check(replay.published&&mesh::meshHash(regular)==hash&&snapshot(*imported.model)==source,
          "regular prism triangulation is deterministic and keeps the original CAD");
    mesh::Mesh improved;
    const auto polished=buildSpacedOcctBoundary(*imported.model,1,options,target,storage,improved,receipt,true,true);
    check(polished.published&&polished.sourceModelUnmodified&&receipt.rebuiltCharts==5&&
          improved.cellCount()==regular.cellCount()&&improved.nodeCount()==regular.nodeCount()&&
          certifyOcctBoundarySurface(improved,OcctStarOptions{}).certified,
          "regular prism quality pass preserves CAD wires and closes the shell at the same size");
    return 5;
}
int recoverCurvedRoof(const std::filesystem::path& folder) {
    BRepBuilderAPI_Sewing sewing(1e-7);
    const auto plane=[&](std::array<gp_Pnt,4> corners){BRepBuilderAPI_MakePolygon wire;for(auto p:corners)wire.Add(p);wire.Close();sewing.Add(BRepBuilderAPI_MakeFace(wire.Wire()).Face());};
    plane({gp_Pnt(0,0,0),gp_Pnt(0,1,0),gp_Pnt(1,1,0),gp_Pnt(1,0,0)});
    plane({gp_Pnt(0,0,0),gp_Pnt(0,0,1),gp_Pnt(0,1,1),gp_Pnt(0,1,0)});
    plane({gp_Pnt(1,0,0),gp_Pnt(1,1,0),gp_Pnt(1,1,1),gp_Pnt(1,0,1)});
    plane({gp_Pnt(0,0,0),gp_Pnt(1,0,0),gp_Pnt(1,0,1),gp_Pnt(0,0,1)});
    plane({gp_Pnt(0,1,0),gp_Pnt(0,1,1),gp_Pnt(1,1,1),gp_Pnt(1,1,0)});
    NCollection_Array2<gp_Pnt> poles(1,3,1,3);
    for(int i=1;i<=3;++i)for(int j=1;j<=3;++j)poles.SetValue(i,j,gp_Pnt((i-1)*.5,(j-1)*.5,1.+(i==2&&j==2?.6:0.)));
    occ::handle<Geom_BezierSurface> surface=new Geom_BezierSurface(poles);
    sewing.Add(BRepBuilderAPI_MakeFace(surface,1e-7).Face());sewing.Perform();
    check(sewing.SewedShape().ShapeType()==TopAbs_SHELL,"curved roof has one closed shell");
    auto solid=BRepBuilderAPI_MakeSolid(TopoDS::Shell(sewing.SewedShape())).Solid();check(BRepLib::OrientClosedSolid(solid),"curved roof orientation");
    auto cad=load(solid,folder/"curved-roof.brep");const auto before=snapshot(*cad.model);
    geo::OcctBoundaryOptions options;options.linearDeflection=.5;options.angularDeflection=1.5;options.maximumChordError=.01;options.maximumWork=50000000;
    mesh::Mesh coarse;const auto staged=geo::detail::buildOcctBoundaryCandidate(*cad.model,1,options,coarse);
    check(staged.published,"coarse curved roof staging");
    mesh::Mesh repaired;const auto result=geo::buildOcctBoundary(*cad.model,1,options,repaired);
    std::cout<<"curved roof: "<<geo::occtBoundaryReasonName(result.reason)<<" nodes="<<result.nodes<<" triangles="<<result.triangles<<" face="<<result.maximumErrorFace<<" sample_error="<<result.maximumFacetSampleError<<" closed="<<result.closed<<" oriented="<<result.oriented<<" manifold="<<result.vertexManifold<<'\n';
    check(result.published&&result.closed&&result.oriented&&result.vertexManifold&&result.sourceModelUnmodified,"curved roof surface recovered");
    check(repaired.nodeCount()>coarse.nodeCount()&&result.maximumFacetSampleError<=options.maximumChordError+1e-7&&result.relativeVolumeError<.01,"curved roof measured fidelity improved");
    const auto audited=geo::detail::auditOcctBoundaryFacets(*cad.model,1,repaired,options);
    check(audited.completed&&audited.sourceModelUnmodified&&snapshot(*cad.model)==before,"curved roof independent facet audit and source guard");
    std::map<geo::OcctBoundaryNodeKey,Vec3> original;
    for(Index i=0;i<coarse.nodeCount();++i)original.emplace(staged.nodeKeys[i],coarse.node(i));
    for(Index i=0;i<repaired.nodeCount();++i)if(repaired.nodeClass(i).dim<2){const auto at=original.find(result.nodeKeys[i]);check(at!=original.end()&&same(at->second,repaired.node(i)),"curved roof retains canonical wire positions and keys");}
    return 7;
}
int localSharedCurvedPreparation(const std::filesystem::path& folder) {
    NCollection_List<TopoDS_Shape> arguments;
    arguments.Append(BRepPrimAPI_MakeBox(gp_Pnt(-2,-2,0),4.,4.,2.).Shape());
    arguments.Append(BRepPrimAPI_MakeCylinder(1.,2.).Shape());
    BRepAlgoAPI_Splitter splitter;splitter.SetArguments(arguments);splitter.SetNonDestructive(true);splitter.Build();
    check(splitter.IsDone()&&!splitter.HasErrors(),"shared curved fixture split");
    auto cad=load(splitter.Shape(),folder/"local-shared-curved.brep");
    check(cad.model->solids().size()==2,"shared curved fixture has both material bodies");
    const auto before=snapshot(*cad.model);int shared=-1,selected=-1;
    for(const auto* face:cad.model->faces())if(face->kind()==geo::SurfaceKind::Cylinder&&cad.model->faceOwners(face->tag()).size()==2)shared=face->tag();
    check(shared>0,"actual shared native curved face");
    const std::array chartPoints{Vec2{.1,.2},Vec2{.2,.2},Vec2{.1,.3}};
    std::vector<Vec2> chartMetric;
    check(geo::detail::metricOcctChart(*cad.model,shared,chartPoints,[](std::uint64_t){},chartMetric)&&
          chartMetric.size()==3&&std::abs(chartMetric[1].u-chartMetric[0].u-.1)<1e-12&&
          std::abs(chartMetric[2].v-chartMetric[0].v-.1)<1e-12&&
          std::abs(chartMetric[2].u-chartMetric[0].u)<1e-12,"unit cylinder physical chart metric");
    for(const auto& wire:cad.model->boundaryUses(shared))for(const auto& use:wire.edges)
        if(!use.degenerated&&!use.seam){selected=use.edgeTag;break;}
    check(selected>0,"selected actual shared curved wire");
    geo::OcctBoundaryOptions options;options.linearDeflection=.01;options.maximumChordError=.021;
    geo::detail::OcctBoundaryConstructionOptions construction;
    const std::array tags{selected};construction.refinedEdges=tags;construction.localLinearDeflection=.001;
    std::map<int,std::uint64_t> refinedHashes,refinedWork;
    std::map<int,std::vector<geo::OcctBoundaryNodeKey>> refinedKeys;
    std::uint64_t refinedPreparationWork=0;
    double coarseError=0,refinedError=0;Index coarseNodes=0,refinedNodes=0;
    const auto sampleError=[&](const mesh::Mesh& boundary){double maximum=0;
        for(Index i=0;i<boundary.cellCount();++i)if(boundary.cell(i).onEntity.tag==shared){
            const auto* t=boundary.cellNodes(i);const auto a=boundary.node(t[0]),b=boundary.node(t[1]),c=boundary.node(t[2]);
            for(const auto p:std::array{(a+b+c)/3.,(a+b)*.5,(b+c)*.5,(c+a)*.5}){Vec3 q;
                check(geo::closestPointOnOcctFace(*cad.model,shared,p,q),"shared curved facet sample retains trimmed CAD authority");maximum=std::max(maximum,(p-q).norm());}}
        return maximum;};
    for(unsigned pass=0;pass<3;++pass){
        geo::OcctBoundaryOutcome preparation;
        auto session=geo::detail::OcctBoundarySession::prepare(*cad.model,options,preparation,
            geo::detail::OcctBoundaryAuditPolicy::PerBoundary,geo::detail::OcctBoundaryFacetPolicy::Immediate,true,pass?&construction:nullptr);
        check(session&&preparation.sourceModelUnmodified,"local shared refinement prepares an immutable owned compound");
        if(pass==1)refinedPreparationWork=preparation.work;
        if(pass==2)check(preparation.work==refinedPreparationWork,"repeated local preparation retains work accounting");
        std::map<geo::OcctBoundaryNodeKey,Vec3> sharedPoints;
        std::map<std::array<geo::OcctBoundaryNodeKey,3>,int> winding;
        std::map<std::array<geo::OcctBoundaryNodeKey,3>,unsigned> uses;
        for(const auto* body:cad.model->solids()){
            mesh::Mesh boundary;const auto receipt=session->boundary(body->tag(),boundary);
            check(receipt.published&&receipt.closed&&receipt.oriented&&receipt.vertexManifold&&receipt.sourceModelUnmodified,"local refinement publishes audited body boundaries");
            for(Index i=0;i<boundary.nodeCount();++i){auto [it,inserted]=sharedPoints.emplace(receipt.nodeKeys[i],boundary.node(i));check(inserted||same(it->second,boundary.node(i)),"locally refined shared keys retain bitwise positions");}
            for(Index i=0;i<boundary.cellCount();++i)if(boundary.cell(i).onEntity.tag==shared){
                const auto* t=boundary.cellNodes(i);std::array key{receipt.nodeKeys[t[0]],receipt.nodeKeys[t[1]],receipt.nodeKeys[t[2]]};
                int sign=1;for(unsigned a=0;a<3;++a)for(unsigned b=a+1;b<3;++b)if(key[b]<key[a]){std::swap(key[a],key[b]);sign=-sign;}
                winding[key]+=sign;++uses[key];}
            if(pass==0){coarseError=std::max(coarseError,sampleError(boundary));coarseNodes+=boundary.nodeCount();}
            if(pass==1){refinedError=std::max(refinedError,sampleError(boundary));refinedNodes+=boundary.nodeCount();refinedHashes[body->tag()]=mesh::meshHash(boundary);refinedKeys[body->tag()]=receipt.nodeKeys;refinedWork[body->tag()]=receipt.work;}
            if(pass==2)check(mesh::meshHash(boundary)==refinedHashes.at(body->tag())&&receipt.nodeKeys==refinedKeys.at(body->tag())&&receipt.work==refinedWork.at(body->tag()),"repeated local refinement retains exact boundary and keys");
        }
        check(!winding.empty(),"locally refined shared interface retains triangles");
        for(const auto& [key,sign]:winding)check(sign==0&&uses.at(key)==2,"both coowners retain the same curved triangles with opposite winding");
        check(snapshot(*cad.model)==before,"local shared preparation preserves exact Source BREP");
    }
    check(refinedNodes>coarseNodes&&refinedError<coarseError*.5,"local edge preparation improves measured shared curved surface error");
    std::cout<<"local shared curved preparation: nodes="<<coarseNodes<<" -> "<<refinedNodes<<" sample_error="<<coarseError<<" -> "<<refinedError<<'\n';
    geo::OcctBoundaryOutcome preparation;auto bad=construction;
    const auto refuses=[&](const geo::detail::OcctBoundaryConstructionOptions& c,geo::OcctBoundaryReason reason){
        auto session=geo::detail::OcctBoundarySession::prepare(*cad.model,options,preparation,
            geo::detail::OcctBoundaryAuditPolicy::PerBoundary,geo::detail::OcctBoundaryFacetPolicy::Immediate,true,&c);
        check(!session&&preparation.reason==reason&&!preparation.published&&snapshot(*cad.model)==before,"invalid or exhausted local preparation preserves Source");};
    bad.localLinearDeflection=0;refuses(bad,geo::OcctBoundaryReason::InvalidRequest);
    bad=construction;bad.localLinearDeflection=.02;refuses(bad,geo::OcctBoundaryReason::InvalidRequest);
    bad=construction;bad.maximumIncidentFaces=1;refuses(bad,geo::OcctBoundaryReason::BudgetExceeded);
    const std::array duplicate{selected,selected};bad=construction;bad.refinedEdges=duplicate;refuses(bad,geo::OcctBoundaryReason::InvalidRequest);
    bad=construction;bad.refinedEdges=duplicate;bad.maximumRefinedEdges=1;refuses(bad,geo::OcctBoundaryReason::BudgetExceeded);
    const std::array missing{-1};bad=construction;bad.refinedEdges=missing;refuses(bad,geo::OcctBoundaryReason::InvalidRequest);
    auto limited=options;limited.maximumWork=refinedPreparationWork+refinedWork.begin()->second-1;
    auto session=geo::detail::OcctBoundarySession::prepare(*cad.model,limited,preparation,
        geo::detail::OcctBoundaryAuditPolicy::PerBoundary,geo::detail::OcctBoundaryFacetPolicy::Immediate,true,&construction);
    check(session!=nullptr,"local one-less work reaches extraction");mesh::Mesh sentinel;sentinel.addNode({7,8,9},{3,999});const auto hash=mesh::meshHash(sentinel);
    const auto refused=session->boundary(refinedWork.begin()->first,sentinel);
    check(!refused.published&&refused.reason==geo::OcctBoundaryReason::BudgetExceeded&&mesh::meshHash(sentinel)==hash&&snapshot(*cad.model)==before,"local preparation cumulative work rolls back the destination");
    std::atomic<bool> cancel{true};auto cancelled=options;cancelled.cancellation=&cancel;
    session=geo::detail::OcctBoundarySession::prepare(*cad.model,cancelled,preparation,
        geo::detail::OcctBoundaryAuditPolicy::PerBoundary,geo::detail::OcctBoundaryFacetPolicy::Immediate,true,&construction);
    check(!session&&preparation.reason==geo::OcctBoundaryReason::Cancelled&&preparation.work==0&&snapshot(*cad.model)==before,"local preparation cancellation preserves Source");
    return 29;
}
}
int main(int argc,char**argv) {
    try {
        const auto folder=argc==2?std::filesystem::path(argv[1]):std::filesystem::path("out/occt-boundary-tests");
        std::filesystem::create_directories(folder);
        const auto cylinder=BRepPrimAPI_MakeCylinder(1.,2.).Shape();
        gp_Trsf transform;transform.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(1,2,3)),.63);
        transform.SetTranslationPart(gp_Vec(9,-4,12));
        const std::vector<std::pair<const char*,TopoDS_Shape>> fixtures={
            {"cylinder",cylinder},{"rotated",BRepBuilderAPI_Transform(cylinder,transform,true).Shape()},
            {"sphere",BRepPrimAPI_MakeSphere(1.).Shape()},{"torus",BRepPrimAPI_MakeTorus(2.,.6).Shape()},
            {"tube",BRepAlgoAPI_Cut(cylinder,BRepPrimAPI_MakeCylinder(.55,2.).Shape()).Shape()}};
        geo::OcctBoundaryOptions options;options.linearDeflection=.01;options.maximumChordError=.021;
        int checks=0;
        for(const auto& [name,shape]:fixtures) {
            auto imported=load(shape,folder/(std::string(name)+".brep")); const auto before=snapshot(*imported.model);
            mesh::Mesh boundary;boundary.addNode({-7,-8,-9},{3,999});
            const auto result=geo::buildOcctBoundary(*imported.model,1,options,boundary);
            std::cout<<name<<": "<<geo::occtBoundaryReasonName(result.reason)<<" nodes="<<result.nodes<<" triangles="<<result.triangles
                     <<" sample_error="<<result.maximumFacetSampleError<<" seams="<<result.seamUses<<" poles="<<result.degeneratedUses<<'\n';
            check(result.published&&result.closed&&result.oriented&&result.vertexManifold&&result.sourceModelUnmodified,"valid boundary publication");
            check(result.connectedShells==1&&result.maximumFacetSampleError<=.021&&result.maximumNodeProjectionError<1e-6,"boundary geometry/topology receipt");
            if(std::string_view(name)=="sphere"||std::string_view(name)=="tube"){
                geo::OcctBoundaryOutcome preparation;
                auto session=geo::detail::OcctBoundarySession::prepare(*imported.model,options,preparation);
                check(session!=nullptr,"advanced chart session preparation");
                geo::detail::OcctBoundaryAliasOptions repair;
                repair.maximumFlips=0;
                repair.maximumRebuiltCharts=32;
                geo::detail::OcctBoundaryAliasReceipt receipt;
                mesh::Mesh advanced;
                const auto repaired=session->boundary(1,advanced,repair,receipt);
                check(repaired.published&&repaired.sourceModelUnmodified&&receipt.rebuiltCharts==0&&
                      mesh::meshHash(advanced)==mesh::meshHash(boundary),"advanced chart checks preserve valid seam, pole and hole meshes");
            }
            if(std::string_view(name)=="tube"){
                using namespace geo::detail;
                auto spacedOptions=options;
                spacedOptions.maximumChordError=.009;
                constexpr std::uint64_t storage=64ULL*1024*1024;
                OcctBoundaryAliasReceipt receipt;
                mesh::Mesh spaced;
                const auto rebuilt=buildSpacedOcctBoundary(*imported.model,1,spacedOptions,.35,storage,spaced,receipt);
                std::cout<<"spaced tube: "<<geo::occtBoundaryReasonName(rebuilt.reason)<<" faces="<<receipt.rebuiltCharts
                         <<" curves="<<receipt.rephasedEdges<<" chord="<<rebuilt.maximumFacetSampleError<<'\n';
                check(rebuilt.published&&rebuilt.closed&&rebuilt.oriented&&rebuilt.vertexManifold&&
                      rebuilt.sourceModelUnmodified&&receipt.rebuiltCharts==4&&receipt.rephasedEdges>0,
                      "spaced tube uses shared curves across both periodic faces and annular caps");
                const auto facets=auditOcctBoundaryFacets(*imported.model,1,spaced,spacedOptions);
                check(facets.completed&&rebuilt.relativeVolumeError<.01&&
                      certifyOcctBoundarySurface(spaced,OcctStarOptions{}).certified,
                      "spaced tube passes independent chord, volume and embedding checks");
                const auto hash=mesh::meshHash(spaced);
                auto limited=spacedOptions;
                limited.maximumWork=1;
                const auto capped=buildSpacedOcctBoundary(*imported.model,1,limited,.35,storage,spaced,receipt);
                check(!capped.published&&capped.reason==geo::OcctBoundaryReason::BudgetExceeded&&mesh::meshHash(spaced)==hash,
                      "spaced tube work refusal preserves published destination");
                check(buildSpacedOcctBoundary(*imported.model,1,spacedOptions,.35,storage,spaced,receipt).published&&
                      mesh::meshHash(spaced)==hash&&snapshot(*imported.model)==before,
                      "spaced tube replay is deterministic and preserves source");
                checks+=4;
            }
            if(std::string_view(name)=="cylinder"){
                std::vector<int> faces;for(const auto* f:imported.model->faces())faces.push_back(f->tag());
                geo::OcctBoundaryOutcome prepared;auto session=geo::detail::OcctBoundarySession::prepare(*imported.model,options,prepared);
                check(session!=nullptr,"alias session preparation");geo::detail::OcctBoundaryAliasOptions alias;alias.faceTags=faces;geo::detail::OcctBoundaryAliasReceipt receipt;
                mesh::Mesh value;auto extracted=session->boundary(1,value,alias,receipt);
                check(extracted.published&&receipt.flips==0&&receipt.maximumWorkingPayloadEstimate>0&&mesh::meshHash(value)==mesh::meshHash(boundary),"alias admission preserves an already valid native surface");
                mesh::Mesh sentinel;sentinel.addNode({7,8,9},{3,999});value.swap(sentinel);const auto hash=mesh::meshHash(value);
                alias.maximumWorkingBytes=0;extracted=session->boundary(1,value,alias,receipt);
                check(!extracted.published&&extracted.reason==geo::OcctBoundaryReason::BudgetExceeded&&mesh::meshHash(value)==hash,"alias storage refusal is transactional");
                alias.maximumWorkingBytes=64ULL*1024*1024;alias.maximumFlips=257;extracted=session->boundary(1,value,alias,receipt);
                check(!extracted.published&&extracted.reason==geo::OcctBoundaryReason::InvalidRequest&&mesh::meshHash(value)==hash&&snapshot(*imported.model)==before,"alias invalid request preserves native source and destination");checks+=3;
            }
            check(result.nodeKeys.size()==static_cast<std::size_t>(boundary.nodeCount()),"parallel native keys");
            check(snapshot(*imported.model)==before,"source exact BREP immutability");
            check(boundary.countCells(mesh::CellType::Tri3)==boundary.cellCount(),"Tri3-only boundary");
            for(Index i=0;i<boundary.nodeCount();++i) {
                const auto c=boundary.nodeClass(i);const auto key=result.nodeKeys[static_cast<std::size_t>(i)];
                check(c.dim==key.dimension&&c.tag==key.tag,"native key classification");
                check(c.dim==0?imported.model->vertex(c.tag)!=nullptr:c.dim==1?imported.model->edge(c.tag)!=nullptr:c.dim==2&&imported.model->face(c.tag)!=nullptr,"native boundary node identity");
            }
            checks+=6;
            if(std::string_view(name)=="cylinder") {
                check(result.seamUses==2,"cylinder seam branches stitched");++checks;
                const auto hash=mesh::meshHash(boundary);
                auto refusedOptions=options;refusedOptions.maximumNodes=2;
                auto refused=geo::buildOcctBoundary(*imported.model,1,refusedOptions,boundary);
                check(!refused.published&&refused.reason==geo::OcctBoundaryReason::BudgetExceeded&&refused.destinationUnchanged&&mesh::meshHash(boundary)==hash,"node budget rollback");++checks;
                refusedOptions=options;refusedOptions.maximumChordError=1e-8;
                refused=geo::buildOcctBoundary(*imported.model,1,refusedOptions,boundary);
                check(!refused.published&&refused.reason==geo::OcctBoundaryReason::GeometryError&&refused.maximumFacetSampleError>1e-8&&mesh::meshHash(boundary)==hash,"measured excessive chord error rollback");++checks;
                std::atomic<bool> cancel{true};refusedOptions=options;refusedOptions.cancellation=&cancel;
                refused=geo::buildOcctBoundary(*imported.model,1,refusedOptions,boundary);
                check(!refused.published&&refused.reason==geo::OcctBoundaryReason::Cancelled&&mesh::meshHash(boundary)==hash,"cancellation rollback");++checks;
                Vec3 closest{7,8,9};check(!geo::closestPointOnOcctFace(*imported.model,-1,{0,0,0},closest)&&same(closest,{7,8,9}),"failed closest-face query preserves output");++checks;
            }
            if(std::string_view(name)=="sphere") {check(result.degeneratedUses>0,"sphere pole chart uses handled");++checks;}
            if(std::string_view(name)=="tube") {
                int cap=-1;for(const auto* face:imported.model->faces())if(face->kind()==geo::SurfaceKind::Plane){cap=face->tag();break;}
                check(cap>0,"annular cap found");Vec2 uv;Vec3 projected;const auto* face=imported.model->face(cap);
                Vec3 lo,hi;imported.model->solid(1)->bbox(lo,hi);const double z=face->pointAt({0,0}).z;
                check(!face->projectPoint({0,0,z},uv,projected),"stationary projection still refuses annular hole");
                check(geo::closestPointOnOcctFace(*imported.model,cap,{0,0,z},projected)&&std::abs(projected.norm2()-z*z-.55*.55)<1e-8,"closest trimmed face resolves hole rim");
                check(snapshot(*imported.model)==before,"closest trimmed-face distance preserves source");checks+=3;
            }
        }
        // A point just inside a planar trim has an exact perpendicular closest
        // point. OCCT also admits a tolerance-equivalent edge support here;
        // that support must not mask the closer face point.
        auto planar=load(BRepPrimAPI_MakeBox(1.,1.,1.).Shape(),folder/"near-trim-closest.brep");
        const auto planarBefore=snapshot(*planar.model);int lowerCap=-1;
        for(const auto* face:planar.model->faces()) {
            const auto* descriptor=planar.model->surfaceDescriptor(face->tag());
            if(descriptor&&descriptor->kind==geo::SurfaceKind::Plane&&std::abs(descriptor->axis.z)>.99&&
               std::abs(descriptor->origin.z)<1e-12)lowerCap=face->tag();
        }
        check(lowerCap>0,"native planar lower cap");++checks;
        const Vec3 nearTrim{3e-6,.5,1e-4};Vec3 planarClosest{7,8,9};
        check(geo::closestPointOnOcctFace(*planar.model,lowerCap,nearTrim,planarClosest)&&
              (planarClosest-Vec3{nearTrim.x,nearTrim.y,0}).norm()<1e-12&&
              std::abs((nearTrim-planarClosest).norm()-1e-4)<1e-12,
              "near-trim closest point retains actual perpendicular face minimum");++checks;
        check(snapshot(*planar.model)==planarBefore,"multi-support closest query preserves source");++checks;
        NCollection_List<TopoDS_Shape> arguments,tools;arguments.Append(BRepPrimAPI_MakeBox(2.,1.,1.).Shape());
        tools.Append(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(1,0,0),gp_Dir(1,0,0)),-3.,3.,-3.,3.).Shape());
        BRepAlgoAPI_Splitter splitter;splitter.SetArguments(arguments);splitter.SetTools(tools);splitter.SetNonDestructive(true);splitter.Build();
        check(splitter.IsDone()&&!splitter.HasErrors(),"shared-interface fixture split");
        auto split=load(splitter.Shape(),folder/"shared-interface.brep");check(split.model->solids().size()==2,"two split bodies");
        const auto splitBefore=snapshot(*split.model);
        geo::OcctBoundaryOutcome preparation;
        auto session=geo::detail::OcctBoundarySession::prepare(*split.model,options,preparation);check(session!=nullptr,"single native tessellation session");
        std::map<geo::OcctBoundaryNodeKey,Vec3> allKeys;std::size_t common=0;
        std::map<int,std::uint64_t> bodyHashes;
        std::map<int,std::uint64_t> bodyWork;
        std::map<int,std::vector<geo::OcctBoundaryNodeKey>> bodyKeys;
        std::map<std::array<geo::OcctBoundaryNodeKey,3>,int> interfaceWinding;
        int sharedFace=-1;for(const auto* face:split.model->faces())if(split.model->faceOwners(face->tag()).size()==2)sharedFace=face->tag();
        check(sharedFace>0,"actual native shared face");
        for(const auto* solid:split.model->solids()) {
            mesh::Mesh boundary;const auto result=session->boundary(solid->tag(),boundary);
            check(result.published&&result.sourceModelUnmodified,"shared body boundary");
            bodyHashes[solid->tag()]=mesh::meshHash(boundary);bodyKeys[solid->tag()]=result.nodeKeys;
            bodyWork[solid->tag()]=result.work;
            for(Index i=0;i<boundary.nodeCount();++i) {
                const auto key=result.nodeKeys[static_cast<std::size_t>(i)];auto [found,inserted]=allKeys.emplace(key,boundary.node(i));
                if(!inserted){++common;check(same(found->second,boundary.node(i)),"shared native key has bitwise coordinates");}
            }
            for(Index i=0;i<boundary.cellCount();++i)if(boundary.cell(i).onEntity.tag==sharedFace) {
                const auto* nodes=boundary.cellNodes(i);std::array<geo::OcctBoundaryNodeKey,3> keys{result.nodeKeys[nodes[0]],result.nodeKeys[nodes[1]],result.nodeKeys[nodes[2]]};
                int sign=1;for(int a=0;a<3;++a)for(int b=a+1;b<3;++b)if(keys[b]<keys[a]){std::swap(keys[a],keys[b]);sign=-sign;}
                interfaceWinding[keys]+=sign;
            }
        }
        check(common>=4&&!interfaceWinding.empty(),"shared interface uses identical session identities");
        for(const auto& [triangle,winding]:interfaceWinding){(void)triangle;check(winding==0,"shared interface triangulation has opposite body winding");}
        check(snapshot(*split.model)==splitBefore,"session preserves source BREP");checks+=6;
        geo::OcctBoundaryOutcome explicitPreparation;
        geo::detail::OcctBoundaryConstructionOptions emptyConstruction;
        auto emptySession=geo::detail::OcctBoundarySession::prepare(*split.model,options,explicitPreparation,
            geo::detail::OcctBoundaryAuditPolicy::PerBoundary,geo::detail::OcctBoundaryFacetPolicy::Immediate,false,&emptyConstruction);
        check(emptySession&&explicitPreparation.work==preparation.work,"empty construction retains default preparation work");++checks;
        for(const auto* body:split.model->solids()){
            mesh::Mesh boundary;const auto receipt=emptySession->boundary(body->tag(),boundary);
            check(receipt.published&&receipt.work==bodyWork.at(body->tag())&&mesh::meshHash(boundary)==bodyHashes.at(body->tag())&&receipt.nodeKeys==bodyKeys.at(body->tag()),"empty construction retains exact default boundary keys and work");++checks;
        }
        auto explicitDefault=geo::detail::OcctBoundarySession::prepare(*split.model,options,explicitPreparation,
            geo::detail::OcctBoundaryAuditPolicy::PerBoundary,geo::detail::OcctBoundaryFacetPolicy::Immediate,false);
        check(explicitDefault&&explicitPreparation.work==preparation.work,"explicit default preparation retains exact work");++checks;
        for(const auto* body:split.model->solids()){
            mesh::Mesh boundary;const auto report=explicitDefault->boundary(body->tag(),boundary);
            check(report.published&&report.sourceModelUnmodified&&report.work==bodyWork.at(body->tag())&&
                  mesh::meshHash(boundary)==bodyHashes.at(body->tag())&&report.nodeKeys==bodyKeys.at(body->tag()),
                  "explicit default retains exact native mesh keys and work");++checks;
        }
        auto controlledShared=geo::detail::OcctBoundarySession::prepare(*split.model,options,explicitPreparation,
            geo::detail::OcctBoundaryAuditPolicy::Deferred,geo::detail::OcctBoundaryFacetPolicy::Immediate,true);
        check(controlledShared!=nullptr,"controlled shared source preparation");++checks;
        std::map<geo::OcctBoundaryNodeKey,Vec3> controlledKeys;
        for(const auto* body:split.model->solids()){
            mesh::Mesh boundary;const auto report=controlledShared->boundary(body->tag(),boundary);
            check(report.published&&report.closed&&report.vertexManifold&&!report.sourceModelUnmodified,"controlled shared boundary remains staging");
            for(Index i=0;i<boundary.nodeCount();++i){const auto key=report.nodeKeys[i];auto [it,added]=controlledKeys.emplace(key,boundary.node(i));check(added||same(it->second,boundary.node(i)),"controlled shared identities retain bitwise points");}++checks;
        }
        check(controlledShared->auditSource().sourceModelUnmodified&&snapshot(*split.model)==splitBefore,"controlled whole-source immutability gate");++checks;
        auto analytic=load(BRepPrimAPI_MakeCone(1.,.3,1.7).Shape(),folder/"controlled-cone.brep");const auto analyticBefore=snapshot(*analytic.model);
        auto controlled=geo::detail::OcctBoundarySession::prepare(*analytic.model,options,explicitPreparation,
            geo::detail::OcctBoundaryAuditPolicy::PerBoundary,geo::detail::OcctBoundaryFacetPolicy::Immediate,true);
        check(controlled!=nullptr,"controlled analytic cone preparation");++checks;
        mesh::Mesh controlledBoundary;const auto controlledReceipt=controlled->boundary(1,controlledBoundary);
        check(controlledReceipt.published&&controlledReceipt.closed&&controlledReceipt.vertexManifold&&controlledReceipt.sourceModelUnmodified&&
              controlledReceipt.maximumFacetSampleError<=options.maximumChordError&&snapshot(*analytic.model)==analyticBefore,"controlled cone retains full topology geometry source gates");++checks;
        auto controlledCap=options;controlledCap.maximumWork=explicitPreparation.work+controlledReceipt.work-1;
        geo::OcctBoundaryOutcome capPreparation;
        auto oneLess=geo::detail::OcctBoundarySession::prepare(*analytic.model,controlledCap,capPreparation,
            geo::detail::OcctBoundaryAuditPolicy::PerBoundary,geo::detail::OcctBoundaryFacetPolicy::Immediate,true);
        check(oneLess!=nullptr,"controlled one-less work reaches extraction");++checks;
        const auto controlledHash=mesh::meshHash(controlledBoundary);const auto cappedControlled=oneLess->boundary(1,controlledBoundary);
        check(!cappedControlled.published&&cappedControlled.reason==geo::OcctBoundaryReason::BudgetExceeded&&
              mesh::meshHash(controlledBoundary)==controlledHash&&snapshot(*analytic.model)==analyticBefore,"controlled one-less cumulative work rollback");++checks;
        std::atomic<bool> controlledCancel{true};auto controlledCancelledOptions=options;controlledCancelledOptions.cancellation=&controlledCancel;
        auto cancelledControl=geo::detail::OcctBoundarySession::prepare(*analytic.model,controlledCancelledOptions,capPreparation,
            geo::detail::OcctBoundaryAuditPolicy::PerBoundary,geo::detail::OcctBoundaryFacetPolicy::Immediate,true);
        check(!cancelledControl&&capPreparation.reason==geo::OcctBoundaryReason::Cancelled&&capPreparation.work==0&&snapshot(*analytic.model)==analyticBefore,"controlled preparation cancellation preserves source");++checks;
        checks+=constrainedChartRegions();
        checks+=regularPrismSurface(folder);
        checks+=localSharedCurvedPreparation(folder);
        checks+=recoverCurvedRoof(folder);
        mesh::Mesh sentinel;sentinel.addNode({-7,-8,-9},{3,999});const auto sentinelHash=mesh::meshHash(sentinel);
        const auto capped=session->boundary(1,sentinel,0);
        check(!capped.published&&capped.reason==geo::OcctBoundaryReason::BudgetExceeded&&capped.work==0&&mesh::meshHash(sentinel)==sentinelHash,"interleaved per-call work ceiling rollback");++checks;

        // Deferred piece boundaries are staging receipts, never an assertion
        // that the original model has already passed a final source audit.
        std::atomic<bool> auditCancellation{false};auto deferredOptions=options;
        deferredOptions.cancellation=&auditCancellation;
        auto deferred=geo::detail::OcctBoundarySession::prepare(*split.model,deferredOptions,preparation,
            geo::detail::OcctBoundaryAuditPolicy::Deferred);
        check(deferred&&preparation.sourceModelUnmodified,"deferred session starts with full source guard");++checks;
        for(const auto* solid:split.model->solids()) {
            mesh::Mesh boundary;const auto result=deferred->boundary(solid->tag(),boundary);
            check(result.published&&!result.sourceModelUnmodified,"deferred piece requires final source audit");
            check(mesh::meshHash(boundary)==bodyHashes.at(solid->tag())&&result.nodeKeys==bodyKeys.at(solid->tag()),"deferred batching preserves exact boundary payload and keys");checks+=2;
        }
        const auto zeroAudit=deferred->auditSource(0);
        check(zeroAudit.reason==geo::OcctBoundaryReason::BudgetExceeded&&!zeroAudit.sourceModelUnmodified&&zeroAudit.work==0,"source audit respects shared remaining work");++checks;
        auditCancellation.store(true);const auto cancelledAudit=deferred->auditSource();
        check(cancelledAudit.reason==geo::OcctBoundaryReason::Cancelled&&!cancelledAudit.sourceModelUnmodified&&cancelledAudit.work==0,"final source audit cancellation");++checks;
        auditCancellation.store(false);const auto finalAudit=deferred->auditSource();
        check(finalAudit.reason==geo::OcctBoundaryReason::None&&finalAudit.sourceModelUnmodified&&!finalAudit.published&&finalAudit.work==2,"explicit final whole-source immutability audit");++checks;

        // A deliberate private native mutation after staging must be detected
        // before the transactional consumer can publish. Default sessions
        // continue to refuse it during the individual boundary operation.
        auto mutableSource=load(splitter.Shape(),folder/"source-audit-mutation.brep");
        auto guarded=geo::detail::OcctBoundarySession::prepare(*mutableSource.model,options,preparation);
        auto pending=geo::detail::OcctBoundarySession::prepare(*mutableSource.model,options,preparation,
            geo::detail::OcctBoundaryAuditPolicy::Deferred);
        check(guarded&&pending,"source mutation guard sessions prepared");++checks;
        mesh::Mesh staged;const auto stagedReceipt=pending->boundary(1,staged);
        check(stagedReceipt.published&&!stagedReceipt.sourceModelUnmodified,"piece staged before deliberate source mutation");++checks;
        NCollection_IndexedMap<TopoDS_Shape,TopTools_ShapeMapHasher> nativeVertices;
        TopExp::MapShapes(geo::detail::OcctModelAccess::shape(*mutableSource.model),TopAbs_VERTEX,nativeVertices);
        const auto changedVertex=TopoDS::Vertex(nativeVertices(1));BRep_Builder mutation;
        mutation.UpdateVertex(changedVertex,BRep_Tool::Tolerance(changedVertex)+.001);
        const auto rejectedAudit=pending->auditSource();
        check(rejectedAudit.reason==geo::OcctBoundaryReason::SourceChanged&&!rejectedAudit.sourceModelUnmodified&&!rejectedAudit.published,"deferred final audit refuses changed original source");++checks;
        const auto rejectedBoundary=guarded->boundary(1,sentinel);
        check(!rejectedBoundary.published&&rejectedBoundary.reason==geo::OcctBoundaryReason::SourceChanged&&rejectedBoundary.destinationUnchanged&&mesh::meshHash(sentinel)==sentinelHash,"default boundary still audits source before publication");++checks;

        // OCCT accepts this bounded BRep with positive box volume. Its extra
        // zero-volume shell has two distinct coincident faces sharing native
        // edges. A linear surface must refuse the resulting two-arc links,
        // rather than publish a fictitious vertex-manifold boundary.
        const auto outer=BRepPrimAPI_MakeBox(2.,2.,2.).Solid();
        const auto triangularFace=BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakePolygon(
            gp_Pnt(.5,.5,1),gp_Pnt(1.5,.5,1),gp_Pnt(1,1.5,1),true).Wire()).Face();
        const auto oppositeFace=BRepBuilderAPI_MakeFace(BRepTools::OuterWire(triangularFace)).Face();
        BRep_Builder builder;TopoDS_Shell doubled;builder.MakeShell(doubled);
        builder.Add(doubled,triangularFace);builder.Add(doubled,oppositeFace.Reversed());doubled.Closed(true);
        TopoDS_Solid doubledShellBody;builder.MakeSolid(doubledShellBody);
        for(TopoDS_Iterator iterator(outer);iterator.More();iterator.Next())builder.Add(doubledShellBody,iterator.Value());
        builder.Add(doubledShellBody,doubled);
        auto degenerate=load(doubledShellBody,folder/"doubled-interior-shell.brep");const auto degenerateBefore=snapshot(*degenerate.model);
        const auto refused=geo::buildOcctBoundary(*degenerate.model,1,options,sentinel);
        check(!refused.published&&refused.reason==geo::OcctBoundaryReason::Nonmanifold&&refused.destinationUnchanged,"native doubled-shell simplicial link refusal");
        check(mesh::meshHash(sentinel)==sentinelHash&&snapshot(*degenerate.model)==degenerateBefore,"doubled-shell source and destination unchanged");checks+=2;
        std::cout<<"OCCT boundary checks="<<checks<<" PASS\n";return 0;
    } catch(const Standard_Failure& error){std::cerr<<"OCCT boundary kernel: "<<error.GetMessageString()<<'\n';}
      catch(const std::exception& error){std::cerr<<"OCCT boundary: "<<error.what()<<'\n';}
    return 1;
}
