// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctFrontMesher.hpp"
#include "OcctInwardSites.hpp"
#include "OcctShellProximity.hpp"
#include "OcctLiveCavityRelocation.hpp"
#include "OcctModelInternal.hpp"
#include "OcctPlcMaterialFilter.hpp"
#include "../src/bnd/TetConstraintRecovery.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <sstream>
#include <iostream>
using namespace paacs;using namespace paacs::geo;using namespace paacs::geo::detail;
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
std::string bytes(const OcctModel& model){std::ostringstream stream;BRepTools::Write(OcctModelAccess::shape(model),stream);return stream.str();}
void checkInterruptedSpatialIndex(){
    bnd::detail::PlcGrowingAabbIndex index;
    using Box=bnd::detail::PlcGrowingAabbIndex::Box;
    const auto bounds=[](std::size_t i){const Vec3 p{double(i),0,0};return Box{p,p};};
    const Box all{{0,0,0},{299,0,0}};
    unsigned calls=0;
    const auto interrupted=[&](std::uint64_t){if(++calls==4)throw std::runtime_error("interrupted index build");};
    const auto ignore=[](std::size_t){return true;};
    const auto noCharge=[](std::uint64_t){};
    const auto verify=[&](std::size_t count){
        std::set<std::size_t> found;
        check(index.allAppending(count,all,bounds,noCharge,[&](auto id){return found.insert(id).second;})&&
              found.size()==count,"resumed spatial index visits every cell once");
    };
    bool stopped=false;
    try{index.allAppending(150,all,bounds,interrupted,ignore);}catch(const std::runtime_error&){stopped=true;}
    check(stopped,"initial spatial build interruption");
    verify(150);
    stopped=false;
    calls=0;
    try{index.allAppending(300,all,bounds,interrupted,ignore);}catch(const std::runtime_error&){stopped=true;}
    check(stopped,"spatial append interruption");
    verify(300);
}
void checkLiveSteinerCavity(){
    using namespace bnd;
    using bnd::detail::BistellarComplex;
    BistellarComplex live;
    live.walkConstraints=true;
    live.useLookahead=false;
    // The old non-star witness requires a Steiner point. Unrelated vertex
    // slots make passing the full mesh to a 14-vertex solver fail.
    for(unsigned i=0;i<20;++i)live.vertices.push_back({{10.+i,10,10},i+1});
    const std::array<Vec3,7> points{{{1,-2.6,-2},{2,-4,-4},{2,-4,-1.3333333333333335},
        {2,-1.9000000000000004,-1},{3,-4,-4},{3,-4,-1.3333333333333335},{3,-1.9000000000000004,-1}}};
    for(auto p:points)live.vertices.push_back({p,live.vertices.size()+1});
    for(auto ids:std::array<std::array<TetFacetVertexId,4>,4>{
        {{0,4,5,1},{0,2,6,3},{0,5,6,2},{4,5,6,0}}}){
        for(auto& id:ids)id+=20;
        live.insert({ids,0,{3,0}});
    }
    const BistellarComplex::Edge target{21,23};
    live.requiredEdges.insert(target);
    bnd::detail::TetConstraintRecoveryReport report;
    bnd::detail::constraint_detail::recoverEdgeCavity(live,target,12000,2'000'000,report);
    check(report.steinerNodes==1&&report.steinerAttempts>0&&live.vertices.size()==28&&
          live.edges.contains(target)&&live.protectedEdges.contains(target),"live one-Steiner edge recovery");
    std::vector<TetFacetCell> cells;
    for(std::size_t i=0;i<live.cells.size();++i)if(live.alive[i]){
        cells.push_back(live.cells[i]);
        const auto faces=live.faceKeys(live.cells[i]);
        for(unsigned side=0;side<4;++side){
            const auto next=live.neighbours[i][side];
            const auto& owners=live.faces.at(faces[side]);
            if(next==live.noCell){check(owners.size()==1,"unlinked internal tet face");continue;}
            check(next<live.alive.size()&&live.alive[next]&&owners.contains(next)&&owners.size()==2&&
                  std::ranges::count(live.neighbours[next],i)==1,"reciprocal tet neighbours after Steiner replacement");
        }
    }
    const TetFacetComplexView view{live.vertices,cells,{},{},{},{3,0}};
    std::vector<std::size_t> all(cells.size());
    std::iota(all.begin(),all.end(),0);
    check(auditTetFacetBall(view,all).valid,"live Steiner replacement preserves a positive ball");
}
void checkLiveFacetWalk(){
    using bnd::detail::BistellarComplex;
    BistellarComplex live;
    live.walkConstraints=true;
    live.useLookahead=false;
    for(auto p:std::array<Vec3,5>{{{-1,-1,0},{1,-1,0},{0,1,0},{0,0,-1},{0,0,1}}})
        live.vertices.push_back({p,live.vertices.size()+1});
    for(unsigned side=0;side<3;++side){
        bnd::TetFacetCell cell{{3,4,side,(side+1)%3},0,{3,1}};
        check(live.positive(cell),"facet walk fixture orientation");
        live.insert(cell);
    }
    const BistellarComplex::Face target{0,1,2};
    live.requiredFaces.insert(target);
    auto piercing=live;
    check(bnd::detail::constraint_detail::recoverPiercedFace(piercing,target,{3,1})==1&&
          piercing.protectedFaces.contains(target)&&!piercing.edges.contains({3,4}),
          "piercing-edge removal restores the target facet");
    bnd::detail::constraint_detail::recoverFace(live,target,{3,1},true);
    check(live.facetWalks==1&&live.faces.contains(target)&&live.faces.at(target).size()==2&&
          live.protectedFaces.contains(target)&&live.metrics.convexBipyramids==1,
          "face walk recovers the triangle through an exact convex flip");
    check(!live.flip23(target,{3,4}),"convex flip retains a protected base face");
    live.protectedFaces.clear();
    const auto inverse=live.flip23(target,{3,4});
    check(inverse&&inverse->convexBipyramid&&live.preservesSurvivors(*inverse),"inverse convex flip certificate");
    live.walkConstraints=false;
    const auto checked=live.flip23(target,{3,4});
    const auto sameTet=[](const auto& a,const auto& b) {
        return a.vertices==b.vertices&&a.classification==b.classification;
    };
    check(checked&&inverse->delta==-1&&inverse->delta==checked->delta&&
          std::ranges::equal(inverse->added,checked->added,sameTet),
          "certified flip matches the general cavity check");
    const bnd::TetFacetComplexView view{live.vertices,inverse->added,{},{},{},{3,1}};
    const std::array<std::size_t,3> indices{0,1,2};
    check(bnd::auditTetFacetBall(view,indices).valid,"inverse convex flip produces a positive ball");
}
void checkLiveMaterialRelocation() {
    bnd::detail::BistellarComplex live;
    live.walkConstraints=true;
    for(const auto p:std::array<Vec3,5>{{{0,0,0},{4,0,0},{0,4,0},{0,0,4},{.05,.05,.05}}})
        live.vertices.push_back({p,live.vertices.size()+1});
    const std::array<std::size_t,4> star{0,1,2,3};
    for(unsigned i=0;i<4;++i){
        bnd::TetFacetCell cell{{4,i,(i+1)%4,(i+2)%4},0,{3,1}};
        check(live.positive(cell),"positive interior star");
        live.insert(cell);
    }
    const auto charge=[](std::uint64_t){};
    const auto original=live.vertices[4].position;
    const auto accepted=relocateOcctLiveCavityApex(live,4,star,0,[](Vec3){return true;},charge);
    check(accepted&&accepted->x>.5&&accepted->y>.5&&accepted->z>.5,
          "a closed interior star admits relocation into its kernel");
    check(!relocateOcctLiveCavityApex(live,4,star,0,[](Vec3){return false;},charge)&&
          live.vertices[4].position.x==original.x,"CAD refusal leaves the original apex in place");
    check(!relocateOcctLiveCavityApex(live,0,star,0,[](Vec3){return true;},charge),
          "an open boundary vertex star cannot be relocated");
}
void checkInwardSiteUsesOriginalFacet() {
    mesh::Mesh boundary;
    for(const auto point:std::array<Vec3,3>{{{0,0,0},{0,1,0},{1,0,0}}})
        boundary.addNode(point,{2,1});
    boundary.addCell(mesh::CellType::Tri3,{0,1,2},{2,1});
    auto filled=boundary;
    filled.addNode({0,0,1},{3,1});
    filled.addCell(mesh::CellType::Tet4,{0,1,2,3},{3,1});
    const std::array<unsigned char,1> rejected{0}, accepted{1};
    const auto inside=[](Vec3 p){return p.z>0.;};
    const auto charge=[](std::uint64_t){};
    const auto sites=occtInwardSites(boundary,filled,rejected,.1,1,inside,inside,charge);
    check(sites.size()==1&&sites.front().z>0.&&boundary.nodeCount()==3,
          "an inward site comes from a rejected tetrahedron's original facet");
    check(occtInwardSites(boundary,filled,accepted,.1,1,inside,inside,charge).empty(),
          "accepted tetrahedra do not request new interior sites");
    check(occtInwardSites(boundary,filled,rejected,.1,1,inside,
                          [](Vec3){return false;},charge).empty(),
          "a proposed site also needs a CAD-admitted facet centroid");
    const std::array<Vec3,4> candidates{{{.3,.3,.05},{.8,.8,0},{.3,.3,.5},{.51,.51,0}}};
    const auto near=occtCentroidsNearShell(boundary,candidates,.1,1000,charge);
    check(near==std::vector<std::size_t>{0,3},
          "shell proximity includes the face and nearby edge and excludes empty box corners");
    check(occtCentroidsNearShell(boundary,candidates,.1,1,charge).empty(),
          "a limited shell query leaves the full-classification fallback available");
}
void checkConstraintEscalation() {
    using namespace bnd;
    using bnd::detail::BistellarComplex;
    BistellarComplex live;
    live.walkConstraints=true;
    live.useLookahead=false;
    const mesh::Classification region{3,9};
    const std::array<Vec3,8> points{{{5,-8,0},{0,-4,13},{4,0,0},{-10,3,15},
        {-17,-14,-20},{-18,5,-10},{-19,-14,3},{-2,-16,-8}}};
    for (const auto point : points) live.vertices.push_back({point,live.vertices.size()+1});
    const std::array<std::array<TetFacetVertexId,4>,6> cells{{{4,5,6,7},{6,3,2,5},
        {7,6,2,5},{6,1,2,3},{7,0,2,6},{6,0,2,1}}};
    for (const auto tet : cells) live.insert({tet,0,region});

    const BistellarComplex::Edge target{3,4};
    live.requiredEdges.insert(target);
    auto cheap=live;
    check(!cheap.recoverDirectEdge(target)&&!cheap.recoverEdge(target,4,false)&&cheap.metrics.accepted==0,
          "short-link recovery stalls before its first move");
    const TetFacetComplexView before{live.vertices,live.cells,{},{},{},region};
    const auto shell=tet_facet_detail::boundarySenses(tet_facet_detail::candidateFaceMap(before,live.cells));

    bnd::detail::TetConstraintRecoveryReport report;
    const std::array<BistellarComplex::Edge,1> missingEdges{target};
    bnd::detail::constraint_detail::recoverConstraintFlips(live,missingEdges,{},region,12000,report);
    check(live.protectedEdges.contains(target)&&live.metrics.flipsNM>0,
          "stalled short-link pass reaches a long-star replacement");
    std::vector<TetFacetCell> replaced;
    for (std::size_t i=0;i<live.cells.size();++i) if (live.alive[i]) replaced.push_back(live.cells[i]);
    const TetFacetComplexView after{live.vertices,replaced,{},{},{},region};
    std::vector<std::size_t> ids(replaced.size());
    std::iota(ids.begin(),ids.end(),0);
    check(auditTetFacetBall(after,ids).valid&&
          shell==tet_facet_detail::boundarySenses(tet_facet_detail::candidateFaceMap(after,replaced)),
          "stronger recovery keeps a valid oriented cavity shell");
}
void checkAdjacentConstraintRecovery(bool pair) {
    using namespace bnd;
    using bnd::detail::BistellarComplex;
    using bnd::detail::TetConstraintRecoveryReport;
    using namespace bnd::detail::constraint_detail;
    using Complex=BistellarComplex;
    const std::vector<Vec3> points=pair?std::vector<Vec3>{
        {5,-10.300000000000001,3},{5,-1.7763568394002505e-15,3},
        {7.6384650723747658,1.0387160589922075,0},
        {12.574999999999989,-10.300000000000001,.40000000000000002},
        {12.575000000000014,-10.300000000000001,3},
        {14.38440328417316,3.7309177010787886,.40000000000000002},
        {15.765596715826838,-3.6534855830943727,0}}:std::vector<Vec3>{
        {-25.149999999999999,-7.2999999999999998,10.800000000000001},
        {-23.150000000000002,-10.300000000000001,8.4000000000000004},
        {-23.150000000000002,-10.300000000000001,10.200000000000001},
        {-23.150000000000002,-7.2999999999999998,8.4000000000000004},
        {-23.150000000000002,-7.2999999999999998,10.200000000000001},
        {-20.150000000000002,-10.300000000000001,8.4000000000000004},
        {-20.150000000000002,-10.300000000000001,10.200000000000001},
        {-20.150000000000002,-7.2999999999999998,8.4000000000000004}};
    using Tet=std::array<TetFacetVertexId,4>;
    const std::vector<Tet> connectivity=pair?
        std::vector<Tet>{{4,1,2,0},{3,4,2,0},{1,2,6,4},{2,6,5,1},{3,4,6,2}}:
        std::vector<Tet>{{2,0,4,1},{3,1,4,0},{1,6,4,2},{1,3,4,7},{1,5,4,6},{1,7,4,5},{7,6,4,5}};
    const Complex::Edge target=pair?Complex::Edge{3,5}:Complex::Edge{2,3};
    const std::array<Complex::Edge,1> requested{target};
    const mesh::Classification region{3,9};
    Complex original;
    original.walkConstraints=true;
    original.useLookahead=false;
    for (const auto p:points) original.vertices.push_back({p,original.vertices.size()+1});
    for (const auto& ids:connectivity) {
        TetFacetCell cell{ids,0,region};
        check(original.positive(cell),"positive ABC5454 source tet");
        original.insert(cell);
    }
    original.requiredEdges.insert(target);
    const TetFacetComplexView source{original.vertices,original.cells,{},{},{},region};
    const auto expected=tet_facet_detail::boundarySenses(tet_facet_detail::candidateFaceMap(source,original.cells));
    auto baseline=original;
    TetConstraintRecoveryReport baselineReport;
    const auto selected=baseline.edgeCavity(target);
    recoverEdgeCavity(baseline,target,residualMaximumAttemptWork,0,baselineReport,&selected);
    check(!baseline.edges.contains(target),"ABC5454 crossed strip alone cannot recover its edge");

    auto interrupted=original;
    TetConstraintRecoveryReport interruptedReport;
    unsigned charges=0;
    interrupted.chargeWork=[&](std::uint64_t n){
        if (++charges==20) throw Cancelled{};
        interruptedReport.work+=n;
    };
    bool cancelled=false;
    try {
        recoverResidualCavities(interrupted,requested,{},region,1'000'000,interruptedReport);
    } catch (const Cancelled&) {
        cancelled=true;
    }
    check(cancelled&&!interrupted.edges.contains(target)&&
          interrupted.cells.size()==original.cells.size()&&interrupted.alive==original.alive,
          "cancellation before publication preserves the live patch");

    auto live=original;
    TetConstraintRecoveryReport report;
    constexpr std::uint64_t allowance=10'000;
    live.chargeWork=[&](std::uint64_t n){
        if (n>allowance-report.work) throw WorkLimit{};
        report.work+=n;
    };
    recoverResidualCavities(live,requested,{},region,allowance,report);
    check(!live.edges.contains(target)&&report.residualAttempts==0&&report.work<=allowance,
          "insufficient expansion allowance skips the full cavity solver");

    constexpr std::uint64_t limit=1'000'000;
    live.chargeWork=[&](std::uint64_t n){
        if (n>limit-report.work) throw WorkLimit{};
        report.work+=n;
    };
    recoverResidualCavities(live,requested,{},region,limit,report);
    std::vector<TetFacetCell> output;
    for (std::size_t i=0;i<live.cells.size();++i) if (live.alive[i]) output.push_back(live.cells[i]);
    const TetFacetComplexView after{live.vertices,output,{},{},{},region};
    std::vector<std::size_t> ids(output.size());
    std::iota(ids.begin(),ids.end(),0);
    check(live.edges.contains(target)&&live.protectedEdges.contains(target)&&report.residualExpanded==1&&
          report.residualAttempts==1&&report.residualLinkRejected>0&&report.work<=limit,
          "bounded adjacent cavity search recovers the edge after link filtering");
    check(auditTetFacetBall(after,ids).valid&&live.vertices.size()==points.size()&&
          expected==tet_facet_detail::boundarySenses(tet_facet_detail::candidateFaceMap(after,output)),
          "adjacent cavity recovery keeps a valid oriented shell and the existing nodes");
}

void checkEnclosedVoid(){
    mesh::Mesh surface,seed;
    const std::array<Vec3,8> points{{{0,0,0},{8,0,0},{0,8,0},{0,0,8},
                                      {1,1,1},{2,1,1},{1,2,1},{1,1,2}}};
    for(auto point:points)surface.addNode(point,{2,1});
    for(unsigned shell=0;shell<2;++shell){
        const auto offset=4*shell;
        bnd::TetFacetCell cell{{offset,offset+1,offset+2,offset+3},1,{3,1}};
        if(predicates::orient3d(points[offset],points[offset+1],points[offset+2],points[offset+3])!=predicates::Sign::Positive)
            std::swap(cell.vertices[0],cell.vertices[1]);
        for(unsigned side=0;side<4;++side){
            const auto face=bnd::tet_facet_detail::outwardFace(cell,side);
            std::array<Index,3> ids{Index(face[0]),Index(face[1]),Index(face[2])};
            if(shell==1)std::swap(ids[1],ids[2]);
            surface.addCell(mesh::CellType::Tri3,ids,{2,int(shell+1)});
        }
    }
    const auto before=mesh::meshHash(surface);
    const auto result=bnd::detail::recoverTetConstraints(surface,{3,1},2'000'000,128ULL*1024*1024,1000,nullptr,seed);
    check(std::string_view(result.diagnostic)=="complete_seed"&&seed.cellCount()>0,"closed void constraint recovery");
    const auto quality=algo::analyzeMixedMeshQuality(seed);
    check(quality.valid&&std::abs(quality.totalVolume-511./6)<1e-10&&mesh::meshHash(surface)==before,
          "oriented carving preserves enclosed void volume and source identity");
}
int main(){try{
 checkInterruptedSpatialIndex();
 checkLiveSteinerCavity();
 checkLiveFacetWalk();
 checkLiveMaterialRelocation();
 checkInwardSiteUsesOriginalFacet();
 checkConstraintEscalation();
 checkAdjacentConstraintRecovery(false);
 checkAdjacentConstraintRecovery(true);
 checkEnclosedVoid();
 mesh::Mesh provisional;for(auto p:std::array<Vec3,7>{{{0,0,0},{1,0,0},{0,1,0},{0,0,1},{.1,.1,-.1},{.1,.1,.1},{2,2,-1}}})provisional.addNode(p,{3,1});
 provisional.addCell(mesh::CellType::Tet4,std::array<Index,4>{0,1,2,5},{3,1});provisional.addCell(mesh::CellType::Tet4,std::array<Index,4>{0,1,4,5},{3,1});provisional.addCell(mesh::CellType::Tet4,std::array<Index,4>{0,2,4,5},{3,1});provisional.addCell(mesh::CellType::Tet4,std::array<Index,4>{0,1,2,3},{3,1});
 const auto originalHash=mesh::meshHash(provisional);mesh::Mesh retained;std::uint64_t filterWork=0;
 const auto admits=[](Vec3 p){return p.z>=0&&!(p.x==.25&&p.y==.25&&p.z==.25);};
 const auto filtered=retainOcctPlcMaterial(provisional,4,admits,[&](std::uint64_t n){filterWork+=n;},retained);
 check(filtered.removedInteriorNodes==2&&filtered.removedIncidentCells==2&&filtered.removedCentroidCells==1,"outside-node incident star/centroid removal");check(retained.nodeCount()==5&&retained.cellCount()==1&&retained.cellNodes(0)[3]==4,"native partial compaction");for(Index i=0;i<4;++i)check(retained.node(i).x==provisional.node(i).x&&retained.node(i).y==provisional.node(i).y&&retained.node(i).z==provisional.node(i).z&&retained.nodeClass(i)==provisional.nodeClass(i),"native partial boundary prefix");check(filterWork>0&&mesh::meshHash(provisional)==originalHash,"native partial source transaction");
 const auto retainedHash=mesh::meshHash(retained);bool refused=false;try{retainOcctPlcMaterial(provisional,4,admits,[](std::uint64_t){throw std::runtime_error("cancelled");},retained);}catch(const std::runtime_error&){refused=true;}check(refused&&mesh::meshHash(retained)==retainedHash&&mesh::meshHash(provisional)==originalHash,"native partial cancellation transaction");
 mesh::Mesh separate;const auto distinct=retainOcctPlcMaterial(provisional,4,admits,[](Vec3){return true;},[](std::uint64_t){},separate);
 check(distinct.removedInteriorNodes==2&&distinct.removedIncidentCells==2&&distinct.removedCentroidCells==0&&separate.cellCount()==2,"centroid allowance cannot retain a refused interior node or its incident cells");
 auto imported=OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(1.,1.,1.).Shape());check(imported.imported(),"box import");const auto& model=*imported.model;const auto before=bytes(model);NativeCadMeshOptions options;options.geometry.targetSize=.3;
 mesh::Mesh result;const auto built=buildOcctFrontVolume(model,1,options,5000000,result);check(built.published&&built.reason==NativeCadMeshReason::None,"direct native recovery");check(mesh::auditLinearVolumeBoundary(result).valid&&built.nativeValidation.completed&&built.nativeValidation.sourceModelUnmodified,"complete boundary and native certificate");check(bytes(model)==before,"source BREP changed");const auto saved=mesh::meshHash(result);
 check(std::string_view(built.constraintRecovery.diagnostic)=="not_attempted","default skips constraint recovery");
 options.recoveryMode=NativeCadRecoveryMode::Advanced;options.maximumNativeValidationThreads=2;
 mesh::Mesh advanced;const auto recovered=buildOcctFrontVolume(model,1,options,5000000,advanced);
 if(!recovered.published)std::cerr<<recovered.diagnostic<<" / "<<recovered.constraintRecovery.diagnostic<<'\n';
 check(recovered.published&&recovered.constraintRecovery.candidateCells>0&&recovered.constraintRecovery.seedCells>0,
       "advanced route constructs and uses a constraint seed");
 check(mesh::auditLinearVolumeBoundary(advanced).valid&&recovered.nativeValidation.completed&&
       recovered.nativeValidation.sourceModelUnmodified&&bytes(model)==before,"advanced output and native source checked");
 options.recoveryMode=NativeCadRecoveryMode::Standard;options.maximumNativeValidationThreads=1;
 options.maximumRecoveryPatchStates=0;mesh::Mesh disabled;const auto withoutRegion=buildOcctFrontVolume(model,1,options,5000000,disabled);check(withoutRegion.published&&mesh::meshHash(disabled)==saved&&built.regionRecovery.attempts==0,"accepted front should bypass region fallback");options.maximumRecoveryPatchStates=65537;auto invalid=meshNativeCad(model,options,result);check(!invalid.published&&invalid.reason==NativeCadMeshReason::InvalidRequest&&mesh::meshHash(result)==saved,"invalid region recovery limit publication");options.maximumRecoveryPatchStates=2048;
 options.allowCentroidsWithinChordTolerance=true;mesh::Mesh sampled;const auto band=buildOcctFrontVolume(model,1,options,5000000,sampled);
 check(band.published&&mesh::meshHash(sampled)==saved&&band.nativeValidation.pointsWithinDistanceBand==0&&band.chordSampleValidation.completed&&band.chordSampleValidation.nativeQueries==4*band.nativeValidation.nativeQueries,"explicit chord policy preserves exact box mesh and checks four additional points per tetrahedron");options.allowCentroidsWithinChordTolerance=false;
 auto budget=buildOcctFrontVolume(model,1,options,1,result);check(!budget.published&&budget.reason==NativeCadMeshReason::BudgetExceeded&&mesh::meshHash(result)==saved,"budget publication");
 std::atomic_bool cancel{true};options.geometry.cancellation=&cancel;auto stopped=buildOcctFrontVolume(model,1,options,5000000,result);check(!stopped.published&&stopped.reason==NativeCadMeshReason::Cancelled&&mesh::meshHash(result)==saved,"cancelled publication");
 options.geometry.cancellation=nullptr;options.maximumPrimitiveWorkingBytes=1;auto storage=buildOcctFrontVolume(model,1,options,5000000,result);check(!storage.published&&storage.reason==NativeCadMeshReason::BudgetExceeded&&mesh::meshHash(result)==saved&&bytes(model)==before,"storage budget publication");
 std::cout<<"Native front completion, source identity, budget and cancellation checks passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
