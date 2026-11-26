// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctFacetRefinement.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctModelInternal.hpp"
#include "OcctNativeWireAuthority.hpp"
#include "OcctRecoveryMesher.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/predicates/Predicates.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <NCollection_List.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <algorithm>
#include <bit>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
using namespace paacs;using namespace paacs::geo;using namespace paacs::geo::detail;
namespace {
int checks=0;
void check(bool b,const char* message){++checks;if(!b)throw std::runtime_error(message);}
std::string bytes(const OcctModel& m){std::ostringstream s;BRepTools::Write(OcctModelAccess::shape(m),s);return s.str();}
bool same(Vec3 a,Vec3 b){return std::bit_cast<std::uint64_t>(a.x)==std::bit_cast<std::uint64_t>(b.x)&&std::bit_cast<std::uint64_t>(a.y)==std::bit_cast<std::uint64_t>(b.y)&&std::bit_cast<std::uint64_t>(a.z)==std::bit_cast<std::uint64_t>(b.z);}
void prefix(const mesh::Mesh& a,const mesh::Mesh& b){
 check(b.nodeCount()>=a.nodeCount(),"prefix size");
 for(Index i=0;i<a.nodeCount();++i)check(same(a.node(i),b.node(i))&&a.nodeClass(i)==b.nodeClass(i),"bitwise original node/class prefix");
}
void committed(const mesh::Mesh& before,const mesh::Mesh& after,int body,const OcctFacetRefinementOutcome& outcome){
 check(outcome.published&&!outcome.committedPatches.empty()&&outcome.committedPatches.size()==outcome.addedNodes,
       "published facet receipt covers every committed point");
 std::set<Index> inserted;
 for(const auto& patch:outcome.committedPatches){
  if(patch.bodyTag!=body)continue;
  check(patch.originalTriangle>=0&&patch.originalTriangle<before.cellCount()&&
        patch.insertedVertex>=before.nodeCount()&&patch.insertedVertex<after.nodeCount()&&
        inserted.insert(patch.insertedVertex).second,"committed patch IDs are body-local and unique");
  const auto* original=before.cellNodes(patch.originalTriangle);
  check(std::equal(original,original+3,patch.originalVertices.begin())&&
        before.cell(patch.originalTriangle).onEntity==mesh::Classification{2,patch.faceTag}&&
        after.nodeClass(patch.insertedVertex)==mesh::Classification{2,patch.faceTag},
        "committed patch preserves original directed triangle and actual face tag");
  unsigned children=0;
  for(Index ci=0;ci<after.cellCount();++ci){const auto* p=after.cellNodes(ci);
   if(std::find(p,p+3,patch.insertedVertex)==p+3)continue;
   check(after.cell(ci).onEntity==mesh::Classification{2,patch.faceTag},"committed child keeps actual face identity");
   bool expected=false;
   for(unsigned j=0;j<3;++j){const std::array child{original[j],original[(j+1)%3],patch.insertedVertex};
    for(unsigned shift=0;shift<3;++shift)
     expected|=p[0]==child[shift]&&p[1]==child[(shift+1)%3]&&p[2]==child[(shift+2)%3];}
   check(expected,"committed child has original edge and preserved winding");++children;
  }
  check(children==3,"one committed patch describes exactly three child triangles");
 }
 check(inserted.size()==std::uint64_t(after.nodeCount()-before.nodeCount()),"body-local receipt covers full source-prefix growth");
}
void noCommitted(const OcctFacetRefinementOutcome& outcome){
 check(!outcome.published&&outcome.committedPatches.empty(),"refused batch exposes no committed patch records");
}
int faceAt(const OcctModel& model,Vec3 point){
 int selected=-1;
 for(const auto* face:model.faces()){
  Vec3 closest;if(!closestPointOnOcctFace(model,face->tag(),point,closest))continue;
  if((closest-point).norm()>model.tolerance(2,face->tag()))continue;
  check(selected<0,"owned fixture face probe is unambiguous");selected=face->tag();
 }
 check(selected>0,"owned fixture face probe belongs to actual CAD");return selected;
}
mesh::Classification curveClassAt(const OcctModel& model,Vec3 point){
 for(const auto* vertex:model.vertices())
  if((vertex->position()-point).norm()<=model.tolerance(0,vertex->tag()))return {0,vertex->tag()};
 for(const auto* edge:model.edges())
  if((edge->closestPoint(point)-point).norm()<=model.tolerance(1,edge->tag()))return {1,edge->tag()};
 throw std::runtime_error("owned fixture boundary node lacks actual CAD curve identity");
}
using SourceEdge=std::array<Index,2>;
std::map<SourceEdge,int> sourceCurveEdges(const OcctModel& model,const mesh::Mesh& boundary){
 std::map<SourceEdge,std::set<int>> owners;
 for(Index i=0;i<boundary.cellCount();++i){const auto* p=boundary.cellNodes(i);
  for(unsigned j=0;j<3;++j){SourceEdge edge{p[j],p[(j+1)%3]};std::sort(edge.begin(),edge.end());owners[edge].insert(boundary.cell(i).onEntity.tag);}}
 std::map<SourceEdge,int> result;
 for(const auto& [edge,faces]:owners){
  const auto tag=actualOcctWire(model,boundary.nodeClass(edge[0]),boundary.nodeClass(edge[1]),faces,[]{});
  if(tag>0)result.emplace(edge,tag);
 }
 return result;
}
void addPositiveTet(mesh::Mesh& volume,std::array<Index,4> ids,int body){
 auto sign=predicates::orient3d(volume.node(ids[0]),volume.node(ids[1]),volume.node(ids[2]),volume.node(ids[3]));
 check(sign!=predicates::Sign::Zero,"owned geometric cavity tet is nondegenerate");
 if(sign==predicates::Sign::Negative)std::swap(ids[1],ids[2]);
 check(predicates::orient3d(volume.node(ids[0]),volume.node(ids[1]),volume.node(ids[2]),volume.node(ids[3]))==predicates::Sign::Positive,
       "owned geometric cavity tet has exact positive orientation");
 volume.addCell(mesh::CellType::Tet4,ids,{3,body});
}
void exerciseConcaveGeometryRefinement(){
 // A closed coarse bore can meet its facet/volume bounds while a tetrahedron
 // just outside the faceted bore still has its centroid in the true CAD void.
 // This is geometry approximation, not a request for a better quality score.
 const auto solid=BRepPrimAPI_MakeBox(gp_Pnt(-4,-4,0),8.,8.,1.).Shape();
 // Place the periodic seam on an existing hexagon vertex/vertical edge.
 const auto tool=BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0,0,0),gp_Dir(0,0,1),gp_Dir(0,1,0)),1.,1.).Shape();
 NCollection_List<TopoDS_Shape> arguments,tools;arguments.Append(solid);tools.Append(tool);
 BRepAlgoAPI_Cut cut;cut.SetArguments(arguments);cut.SetTools(tools);cut.SetNonDestructive(true);cut.Build();
 check(cut.IsDone()&&BRepCheck_Analyzer(cut.Shape()).IsValid(),"owned bored box is valid native CAD");
 auto imported=OcctModelAccess::fromShape(cut.Shape());check(imported.imported()&&imported.model->solids().size()==1,"owned bored box import");
 const auto& model=*imported.model;const int body=model.solids().front()->tag();const auto before=bytes(model);
 const int bottom=faceAt(model,{3,0,0}),top=faceAt(model,{3,0,1}),bore=faceAt(model,{1,0,.5});
 const double c=std::sqrt(3.)/2.,d=4./std::sqrt(3.);
 const std::array<Vec3,6> inner{{{c,.5,0},{0,1,0},{-c,.5,0},{-c,-.5,0},{0,-1,0},{c,-.5,0}}};
 const std::array<Vec3,10> outer{{{4,d,0},{4,4,0},{0,4,0},{-4,4,0},{-4,d,0},
                                {-4,-d,0},{-4,-4,0},{0,-4,0},{4,-4,0},{4,-d,0}}};
 constexpr std::array<unsigned,6> anchors{0,2,4,5,7,9};
 mesh::Mesh boundary;std::array<std::array<Index,6>,2> innerIds{};std::array<std::array<Index,10>,2> outerIds{};
 for(unsigned layer=0;layer<2;++layer){
  for(unsigned i=0;i<inner.size();++i){auto p=inner[i];p.z=layer;innerIds[layer][i]=boundary.addNode(p,curveClassAt(model,p));}
  for(unsigned i=0;i<outer.size();++i){auto p=outer[i];p.z=layer;outerIds[layer][i]=boundary.addNode(p,curveClassAt(model,p));}
 }
 for(unsigned i=0;i<inner.size();++i){const auto j=(i+1)%inner.size();
  std::vector<unsigned> chain{anchors[i]};while(chain.back()!=anchors[j])chain.push_back((chain.back()+1)%outer.size());
  for(unsigned layer=0;layer<2;++layer){
   for(std::size_t k=0;k+1<chain.size();++k){std::array tri{innerIds[layer][i],outerIds[layer][chain[k]],outerIds[layer][chain[k+1]]};
    if(!layer)std::swap(tri[1],tri[2]);boundary.addCell(mesh::CellType::Tri3,tri,{2,layer?top:bottom});}
   std::array tri{innerIds[layer][i],outerIds[layer][anchors[j]],innerIds[layer][j]};
   if(!layer)std::swap(tri[1],tri[2]);boundary.addCell(mesh::CellType::Tri3,tri,{2,layer?top:bottom});
  }
 }
 Index selected=-1;
 for(unsigned i=0;i<inner.size();++i){const auto j=(i+1)%inner.size();
  boundary.addCell(mesh::CellType::Tri3,{innerIds[0][i],innerIds[1][i],innerIds[1][j]},{2,bore});
  if(i==5)selected=boundary.cellCount();
  boundary.addCell(mesh::CellType::Tri3,{innerIds[0][i],innerIds[1][j],innerIds[0][j]},{2,bore});
 }
 for(unsigned i=0;i<outer.size();++i){const auto j=(i+1)%outer.size();
  const auto midpoint=(boundary.node(outerIds[0][i])+boundary.node(outerIds[1][j]))*.5;
  const int face=faceAt(model,midpoint);
  boundary.addCell(mesh::CellType::Tri3,{outerIds[0][i],outerIds[0][j],outerIds[1][j]},{2,face});
  boundary.addCell(mesh::CellType::Tri3,{outerIds[0][i],outerIds[1][j],outerIds[1][i]},{2,face});
 }
 check(boundary.nodeCount()==32&&boundary.cellCount()==64&&selected>=0,"owned closed coarse bore populations");
 const auto certificate=certifyOcctBoundarySurface(boundary,OcctStarOptions{});
 check(certificate.certified&&certificate.embedded&&certificate.closedOrientedVertexManifold,
       "owned coarse bore PLC is embedded and exactly closed");
 const double relativeVolumeError=std::abs(certificate.signedVolume-model.volume(body))/model.volume(body);
 check(relativeVolumeError>0&&relativeVolumeError<.01,"coarse bore PLC meets one-percent volume criterion");
 OcctBoundaryOptions audit;audit.maximumChordError=.14;
 const auto initialAudit=auditOcctBoundaryFacets(model,body,boundary,audit);
 check(initialAudit.completed&&initialAudit.sourceModelUnmodified&&initialAudit.maximumFacetSampleError<=.14,
       "coarse bore PLC meets original CAD facet bound");
 const auto original=boundary;const auto curves=sourceCurveEdges(model,boundary);check(!curves.empty(),"owned source has actual CAD curve IDs");
 const Vec3 apex{1.33,0,.5};const auto* face=boundary.cellNodes(selected);
 const Vec3 center=(boundary.node(face[0])+boundary.node(face[1])+boundary.node(face[2])+apex)*.25;
 check(std::hypot(center.x,center.y)<.991,"analytic coarse centroid is in true cylindrical void");
 BRepClass3d_SolidClassifier classifier(OcctModelAccess::solid(model,body));const double tolerance=model.tolerance(3,body);
 classifier.Perform(gp_Pnt(center.x,center.y,center.z),tolerance);check(classifier.State()==TopAbs_OUT,"strict native CAD rejects coarse tet centroid");
 classifier.Perform(gp_Pnt(apex.x,apex.y,apex.z),tolerance);check(classifier.State()==TopAbs_IN,"same owned cavity apex is strictly material");
 mesh::Mesh coarse;for(Index i=0;i<boundary.nodeCount();++i)coarse.addNode(boundary.node(i),boundary.nodeClass(i));
 const auto coarseApex=coarse.addNode(apex,{3,body});addPositiveTet(coarse,{face[0],face[1],face[2],coarseApex},body);
 check(algo::analyzeMixedMeshQuality(coarse).valid,"coarse tet is oriented and nondegenerate despite wrong CAD centroid");
 const std::array selections{OcctFacetSelection{body,selected}};OcctFacetRefinementOptions options;
 options.selections=selections;options.maximumProjectionDistance=.14;options.maximumAddedNodes=1;
 const std::array input{OcctFacetBoundary{body,&boundary}};const auto refined=refineOcctFacets(model,input,options);
 check(refined.published&&refined.addedNodes==1&&refined.addedTriangles==2&&refined.maximumDisplacement>.11&&refined.maximumDisplacement<.14,
       "bounded native curved-facet projection addresses geometric mismatch");
 committed(original,boundary,body,refined);
 prefix(original,boundary);check(sourceCurveEdges(model,boundary)==curves,"all actual source-curve node IDs and CAD tags remain unchanged");
 check(auditOcctBoundaryFacets(model,body,boundary,audit).completed,"refined bore still meets unchanged facet bound");
 const auto newNode=original.nodeCount();check(boundary.nodeClass(newNode)==mesh::Classification{2,bore},"new point retains actual bore face identity");
 mesh::Mesh repaired;for(Index i=0;i<boundary.nodeCount();++i)repaired.addNode(boundary.node(i),boundary.nodeClass(i));
 const auto repairedApex=repaired.addNode(apex,{3,body});unsigned cells=0;
 for(Index i=0;i<boundary.cellCount();++i){const auto* p=boundary.cellNodes(i);
  if(std::find(p,p+3,newNode)==p+3)continue;
  const Vec3 centroid=(boundary.node(p[0])+boundary.node(p[1])+boundary.node(p[2])+apex)*.25;
  check(centroid.x>1.01,"analytic refined centroid is beyond true bore radius");
  classifier.Perform(gp_Pnt(centroid.x,centroid.y,centroid.z),tolerance);
  check(classifier.State()==TopAbs_IN,"strict native CAD admits every refined tet centroid");
  addPositiveTet(repaired,{p[0],p[1],p[2],repairedApex},body);++cells;
 }
 check(cells==3&&algo::analyzeMixedMeshQuality(repaired).valid,"actual facet children yield three positive Tet4 cells");
 check(bytes(model)==before,"geometric refinement preserves original CAD and native tolerances");
}
void exerciseMixedFacetCandidates(){
 // This deliberately coarse annular PLC is an embedded projection-candidate
 // fixture, not a claim that its volume meets a Native meshing tolerance.
 const gp_Ax2 axes(gp_Pnt(0,0,0),gp_Dir(0,0,1),gp_Dir(0,1,0));
 NCollection_List<TopoDS_Shape> arguments,tools;
 arguments.Append(BRepPrimAPI_MakeCylinder(axes,1.3,1.).Shape());tools.Append(BRepPrimAPI_MakeCylinder(axes,1.,1.).Shape());
 BRepAlgoAPI_Cut cut;cut.SetArguments(arguments);cut.SetTools(tools);cut.SetNonDestructive(true);cut.Build();
 check(cut.IsDone()&&BRepCheck_Analyzer(cut.Shape()).IsValid(),"mixed-candidate annulus is valid native CAD");
 auto imported=OcctModelAccess::fromShape(cut.Shape());check(imported.imported()&&imported.model->solids().size()==1,"mixed-candidate annulus import");
 const auto& model=*imported.model;const auto before=bytes(model);const int body=model.solids().front()->tag();
 const int bottom=faceAt(model,{1.15,0,0}),top=faceAt(model,{1.15,0,1}),innerFace=faceAt(model,{1,0,.5}),outerFace=faceAt(model,{1.3,0,.5});
 const double c=std::sqrt(3.)/2.;const std::array<Vec3,6> ring{{{c,.5,0},{0,1,0},{-c,.5,0},{-c,-.5,0},{0,-1,0},{c,-.5,0}}};
 mesh::Mesh source;std::array<std::array<Index,6>,2> inner{},outer{};
 for(unsigned layer=0;layer<2;++layer){
  for(unsigned i=0;i<ring.size();++i){auto p=ring[i];p.z=layer;inner[layer][i]=source.addNode(p,curveClassAt(model,p));}
  for(unsigned i=0;i<ring.size();++i){auto p=ring[i]*1.3;p.z=layer;outer[layer][i]=source.addNode(p,curveClassAt(model,p));}
 }
 Index unavailable=-1,planar=-1,curved=-1;
 for(unsigned i=0;i<ring.size();++i){const auto j=(i+1)%ring.size();
  for(unsigned layer=0;layer<2;++layer){
   std::array first{inner[layer][i],outer[layer][i],outer[layer][j]},second{inner[layer][i],outer[layer][j],inner[layer][j]};
   if(!layer){std::swap(first[1],first[2]);std::swap(second[1],second[2]);}
   if(i==5&&layer==1)planar=source.cellCount();
   source.addCell(mesh::CellType::Tri3,first,{2,layer?top:bottom});
   if(i==5&&layer==1)unavailable=source.cellCount();
   source.addCell(mesh::CellType::Tri3,second,{2,layer?top:bottom});
  }
  source.addCell(mesh::CellType::Tri3,{inner[0][i],inner[1][i],inner[1][j]},{2,innerFace});
  if(i==5)curved=source.cellCount();
  source.addCell(mesh::CellType::Tri3,{inner[0][i],inner[1][j],inner[0][j]},{2,innerFace});
  source.addCell(mesh::CellType::Tri3,{outer[0][i],outer[0][j],outer[1][j]},{2,outerFace});
  source.addCell(mesh::CellType::Tri3,{outer[0][i],outer[1][j],outer[1][i]},{2,outerFace});
 }
 check(source.nodeCount()==24&&source.cellCount()==48&&unavailable>=0&&planar>=0&&curved>=0,"mixed-candidate source populations");
 check(certifyOcctBoundarySurface(source,OcctStarOptions{}).certified,"mixed-candidate source is an exactly closed embedded PLC");
 const auto* bad=source.cellNodes(unavailable);const auto badCenter=(source.node(bad[0])+source.node(bad[1])+source.node(bad[2]))/3.;Vec2 uv;Vec3 projected;
 check(std::hypot(badCenter.x,badCenter.y)<.98&&!model.face(top)->projectPoint(badCenter,uv,projected),
       "actual annular cap trim refuses its coarse chord triangle center");
 const auto hash=mesh::meshHash(source);const auto curves=sourceCurveEdges(model,source);
 check(curves.size()>=24,"mixed-candidate source retains every physical circle segment identity");
 const std::array mixed{OcctFacetSelection{body,curved},OcctFacetSelection{body,unavailable}};
 OcctFacetRefinementOptions options;options.selections=mixed;options.maximumProjectionDistance=.14;options.maximumAddedNodes=1;
 auto current=source;const std::array input{OcctFacetBoundary{body,&current}};
 const auto refused=refineOcctFacets(model,input,options);
 check(!refused.published&&refused.reason==OcctFacetRefinementReason::UnsupportedGeometry&&
       std::string_view(refused.diagnostic)=="native-trim-projection-failed"&&refused.failedTriangle==unavailable&&
       refused.unavailableCandidates==0&&mesh::meshHash(current)==hash,"default mixed projection batch rolls back all candidates");
 noCommitted(refused);
 options.skipUnavailableCandidates=true;const auto admitted=refineOcctFacets(model,input,options);
 check(admitted.published&&admitted.addedNodes==1&&admitted.addedTriangles==2&&admitted.unavailableCandidates==1&&
       admitted.unchangedCandidates==0,"opt-in omission retains a usable curved facet");
 committed(source,current,body,admitted);
 prefix(source,current);check(sourceCurveEdges(model,current)==curves&&certifyOcctBoundarySurface(current,OcctStarOptions{}).certified,
       "mixed-candidate publication preserves curves and complete embedding certificate");
 const std::array withPlanar{OcctFacetSelection{body,curved},OcctFacetSelection{body,unavailable},OcctFacetSelection{body,planar}};
 current=source;options.selections=withPlanar;options.requireGeometricChange=true;
 const auto geometric=refineOcctFacets(model,input,options);
 check(geometric.published&&geometric.addedNodes==1&&geometric.unchangedCandidates==2&&geometric.unavailableCandidates==0,
       "geometry mode skips planar candidates without projecting them");
 committed(source,current,body,geometric);
 const std::array onlyPlanar{OcctFacetSelection{body,planar}};current=source;options.selections=onlyPlanar;
 const auto unchanged=refineOcctFacets(model,input,options);
 check(!unchanged.published&&std::string_view(unchanged.diagnostic)=="no-usable-facet-candidates"&&
       unchanged.unchangedCandidates==1&&unchanged.nativeQueries==0&&mesh::meshHash(current)==hash,
       "all-no-op geometry selection refuses atomically");
 noCommitted(unchanged);
 const std::array onlyUnavailable{OcctFacetSelection{body,unavailable}};options.selections=onlyUnavailable;options.requireGeometricChange=false;
 const auto empty=refineOcctFacets(model,input,options);
 check(!empty.published&&std::string_view(empty.diagnostic)=="no-usable-facet-candidates"&&
       empty.unavailableCandidates==1&&mesh::meshHash(current)==hash,"all unavailable projections refuse atomically");
 noCommitted(empty);
 options.selections=mixed;
 const auto limited=[&](OcctFacetRefinementOptions limits,OcctFacetRefinementReason reason){
  current=source;const auto result=refineOcctFacets(model,input,limits);
  check(!result.published&&result.reason==reason&&mesh::meshHash(current)==hash,"candidate omission preserves hard stops and rollback");
  noCommitted(result);
 };
 auto budget=options;budget.limits.maximumWork=1;limited(budget,OcctFacetRefinementReason::BudgetExceeded);
 budget=options;budget.limits.maximumWorkingBytes=1;limited(budget,OcctFacetRefinementReason::BudgetExceeded);
 budget=options;budget.limits.maximumBoundaryNodes=source.nodeCount();limited(budget,OcctFacetRefinementReason::BudgetExceeded);
 std::atomic<bool> cancelled{true};budget=options;budget.limits.cancellation=&cancelled;limited(budget,OcctFacetRefinementReason::Cancelled);
 // A route-selection cap counts usable patches, not attempted raw IDs.
 // More than 32 planar no-ops precede the useful curved proposal.
 std::vector<OcctFacetSelection> delayed(40,OcctFacetSelection{body,planar});
 delayed.push_back({body,unavailable});delayed.push_back({body,curved});delayed.push_back({body,curved-1});
 current=source;auto capped=options;capped.selections=delayed;capped.requireGeometricChange=true;
 capped.maximumCandidatePatches=1;capped.maximumAddedNodes=1;
 const auto selected=refineOcctFacets(model,input,capped);
 check(selected.published&&selected.addedNodes==1&&selected.unchangedCandidates==41&&
       selected.committedPatches.size()==1&&selected.committedPatches.front().originalTriangle==curved,
       "usable patch cap scans beyond raw first 32 no-ops without consuming their insertion allowance");
 committed(source,current,body,selected);
 auto hard=capped;hard.maximumAddedNodes=0;limited(hard,OcctFacetRefinementReason::BudgetExceeded);
 hard=capped;hard.limits.maximumWork=1;limited(hard,OcctFacetRefinementReason::BudgetExceeded);
 hard=capped;hard.limits.maximumWorkingBytes=1;limited(hard,OcctFacetRefinementReason::BudgetExceeded);
 hard=capped;hard.limits.cancellation=&cancelled;limited(hard,OcctFacetRefinementReason::Cancelled);
 current=source;capped.selections=onlyPlanar;
 const auto noUseful=refineOcctFacets(model,input,capped);
 check(!noUseful.published&&std::string_view(noUseful.diagnostic)=="no-usable-facet-candidates"&&
       mesh::meshHash(current)==hash,"usable cap with no geometric change leaves fallback source untouched");
 noCommitted(noUseful);
 check(bytes(model)==before,"candidate selection never changes Native source or tolerances");
}
}
int main(){try{
 auto box=OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(1.,1.,1.).Shape());check(box.imported(),"box import");
 const auto before=bytes(*box.model);OcctBoundaryOptions bo;bo.linearDeflection=.03;bo.maximumChordError=.05;
 mesh::Mesh original;check(buildOcctBoundary(*box.model,1,bo,original).published,"box boundary");
 const std::array selections{OcctFacetSelection{1,0},OcctFacetSelection{1,1},OcctFacetSelection{1,0}};
 OcctFacetRefinementOptions o;o.selections=selections;o.maximumProjectionDistance=.1;
 mesh::Mesh refined=original;const std::array one{OcctFacetBoundary{1,&refined}};
 const auto success=refineOcctFacets(*box.model,one,o);
 check(success.published&&success.addedNodes==2&&success.addedTriangles==4&&success.nativeQueries==4&&success.maximumDisplacement==0,"native planar facet split and duplicate selection deduplication");
 committed(original,refined,1,success);
 prefix(original,refined);check(bytes(*box.model)==before,"box source BREP unchanged");
 // Only the first two triangles were replaced; every other directed triangle
 // and face class must remain identical, in its original relative order.
 for(Index i=2;i<original.cellCount();++i){const auto* a=original.cellNodes(i);const auto* b=refined.cellNodes(i+4);
  check(original.cell(i).onEntity==refined.cell(i+4).onEntity&&std::equal(a,a+3,b),"unaffected directed triangles unchanged");}
 check(auditOcctBoundaryFacets(*box.model,1,refined,bo).completed,"refined box CAD facet checks");
 mesh::Mesh volume;OcctStarOptions star;star.nativeAdmissionModel=box.model.get();star.nativeAdmissionBodyTag=1;
 const auto fan=buildOcctStarVolume(refined,1,star,volume);const auto quality=algo::analyzeMixedMeshQuality(volume);
 check(fan.published&&quality.valid&&volume.countCells(mesh::CellType::Tet4)==16&&std::abs(quality.totalVolume-1.)<1e-12&&mesh::auditLinearVolumeBoundary(volume).valid,"facet split participates in a valid native-admitted volume");
 mesh::Mesh repeated=original;const std::array repeat{OcctFacetBoundary{1,&repeated}};const auto deterministic=refineOcctFacets(*box.model,repeat,o);
 check(deterministic.published&&mesh::meshHash(repeated)==mesh::meshHash(refined)&&deterministic.workUsed==success.workUsed,"deterministic mesh and work");
 committed(original,repeated,1,deterministic);
 std::vector<mesh::Classification> parent(1,{-1,-1});for(const auto* face:box.model->faces()){if(std::size_t(face->tag())>=parent.size())parent.resize(face->tag()+1,{-1,-1});parent[face->tag()]={2,face->tag()};}
 OcctPartitionRecoveryOptions recovery;recovery.enabled=true;recovery.maximumCoreAttempts=4;
 OcctStarOptions coreLimits;coreLimits.maximumWork=20000000;
 mesh::Mesh unseededSource=refined,unseededVolume;unseededVolume.addNode({9,8,7},{3,99});const auto sentinel=mesh::meshHash(unseededVolume),inputHash=mesh::meshHash(unseededSource);
 const auto unseeded=recoverOcctPlcVolume(unseededSource,1,*box.model,*box.model,1,parent,coreLimits,recovery,unseededVolume);
 check(!unseeded.published&&mesh::meshHash(unseededVolume)==sentinel&&mesh::meshHash(unseededSource)==inputHash,"minimal unseeded recovery failure remains atomic");
 recovery.nativeKernelSeed=true;mesh::Mesh seededSource=refined,seededVolume;
 const auto seeded=recoverOcctPlcVolume(seededSource,1,*box.model,*box.model,1,parent,coreLimits,recovery,seededVolume);
 const auto seededQuality=algo::analyzeMixedMeshQuality(seededVolume);
 check(seeded.published&&seeded.coreAttempts==2&&seeded.quadFlips>0&&seeded.nativeQueries>0&&seededVolume.nodeCount()==refined.nodeCount()+1&&seededVolume.countCells(mesh::CellType::Tet4)==16&&seededQuality.valid&&std::abs(seededQuality.totalVolume-1.)<1e-12&&mesh::auditLinearVolumeBoundary(seededVolume).valid,"native kernel point plus diagonal retry resolves the minimal recovery failure");
 prefix(refined,seededVolume);check(auditOcctBoundaryFacets(*box.model,1,seededSource,bo).completed&&bytes(*box.model)==before,"recovered facet changes retain original CAD and source gates");
 coreLimits.maximumWork=seeded.workUsed-1;mesh::Mesh refusedSource=refined,refusedVolume;refusedVolume.addNode({9,8,7},{3,99});
 const auto exhausted=recoverOcctPlcVolume(refusedSource,1,*box.model,*box.model,1,parent,coreLimits,recovery,refusedVolume);
 check(!exhausted.published&&exhausted.reason==OcctPlcRecoveryReason::BudgetExceeded&&mesh::meshHash(refusedSource)==inputHash&&mesh::meshHash(refusedVolume)==sentinel,"seeded recovery exhausted budget preserves both destinations");
 const auto refuses=[&](OcctFacetRefinementOptions limited,OcctFacetRefinementReason reason){
  mesh::Mesh sentinel=original;const auto hash=mesh::meshHash(sentinel);const std::array input{OcctFacetBoundary{1,&sentinel}};
  const auto r=refineOcctFacets(*box.model,input,limited);
  check(!r.published&&r.reason==reason&&mesh::meshHash(sentinel)==hash&&bytes(*box.model)==before,"refusal preserves destination and source");
  noCommitted(r);
 };
 auto limited=o;limited.limits.maximumWork=success.workUsed-1;refuses(limited,OcctFacetRefinementReason::BudgetExceeded);
 limited=o;limited.maximumAddedNodes=1;refuses(limited,OcctFacetRefinementReason::BudgetExceeded);
 limited=o;limited.limits.maximumBoundaryNodes=9;refuses(limited,OcctFacetRefinementReason::BudgetExceeded);
 limited=o;limited.limits.maximumWorkingBytes=1;refuses(limited,OcctFacetRefinementReason::BudgetExceeded);
 limited=o;limited.limits.maximumOutputBytes=1;refuses(limited,OcctFacetRefinementReason::BudgetExceeded);
 limited=o;limited.maximumProjectionDistance=0;refuses(limited,OcctFacetRefinementReason::InvalidRequest);
 const std::array badSelection{OcctFacetSelection{1,original.cellCount()}};limited=o;limited.selections=badSelection;refuses(limited,OcctFacetRefinementReason::InvalidRequest);
 std::atomic<bool> cancelled{true};limited=o;limited.limits.cancellation=&cancelled;refuses(limited,OcctFacetRefinementReason::Cancelled);
 NCollection_List<TopoDS_Shape> args;args.Append(BRepPrimAPI_MakeBox(gp_Pnt(-1,0,0),1.,1.,1.).Shape());args.Append(BRepPrimAPI_MakeBox(1.,1.,1.).Shape());
 BRepAlgoAPI_Splitter splitter;splitter.SetArguments(args);splitter.SetNonDestructive(true);splitter.Build();check(splitter.IsDone()&&!splitter.HasErrors(),"shared fixture");
 auto shared=OcctModelAccess::fromShape(splitter.Shape());check(shared.imported()&&shared.model->solids().size()==2,"shared bodies");const auto sharedBefore=bytes(*shared.model);
 int sharedFace=-1;for(const auto* face:shared.model->faces())if(shared.model->faceOwners(face->tag()).size()==2)sharedFace=face->tag();check(sharedFace>0,"actual shared face");
 OcctBoundaryOutcome prep;auto session=OcctBoundarySession::prepare(*shared.model,bo,prep);check(bool(session),"shared preparation");
 std::array<mesh::Mesh,2> boundaries;const std::array both{OcctFacetBoundary{shared.model->solids()[0]->tag(),&boundaries[0]},OcctFacetBoundary{shared.model->solids()[1]->tag(),&boundaries[1]}};
 for(unsigned i=0;i<2;++i)check(session->boundary(both[i].bodyTag,boundaries[i]).published,"shared body extraction");
 Index selected=-1;for(Index i=0;i<boundaries[0].cellCount();++i)if(boundaries[0].cell(i).onEntity.tag==sharedFace){selected=i;break;}
 check(selected>=0,"shared triangle selection");const std::array sharedSelection{OcctFacetSelection{both[0].bodyTag,selected}};o.selections=sharedSelection;
 const auto saved=boundaries;const auto hash0=mesh::meshHash(boundaries[0]),hash1=mesh::meshHash(boundaries[1]);
 const auto missing=refineOcctFacets(*shared.model,std::span(both).first(1),o);
 check(!missing.published&&std::string_view(missing.diagnostic)=="missing-shared-face-coowner"&&mesh::meshHash(boundaries[0])==hash0,"missing coowner refuses atomically");
 noCommitted(missing);
 auto candidateSelection=o;candidateSelection.skipUnavailableCandidates=true;candidateSelection.requireGeometricChange=true;
 const auto stillMissing=refineOcctFacets(*shared.model,std::span(both).first(1),candidateSelection);
 check(!stillMissing.published&&std::string_view(stillMissing.diagnostic)=="missing-shared-face-coowner"&&
       stillMissing.unavailableCandidates==0&&stillMissing.unchangedCandidates==0&&mesh::meshHash(boundaries[0])==hash0,
       "candidate omission and plane skipping cannot bypass actual coowner requirements");
 noCommitted(stillMissing);
 const auto coupled=refineOcctFacets(*shared.model,both,o);check(coupled.published&&coupled.addedNodes==2&&coupled.addedTriangles==4,"both owners split together");
 committed(saved[0],boundaries[0],both[0].bodyTag,coupled);committed(saved[1],boundaries[1],both[1].bodyTag,coupled);
 using Point=std::array<std::uint64_t,3>;using Triangle=std::array<Point,3>;std::map<Triangle,std::pair<int,unsigned>> winding;
 for(unsigned owner=0;owner<2;++owner){prefix(saved[owner],boundaries[owner]);check(auditOcctBoundaryFacets(*shared.model,both[owner].bodyTag,boundaries[owner],bo).completed,"shared facet CAD checks");
  for(Index i=0;i<boundaries[owner].cellCount();++i)if(boundaries[owner].cell(i).onEntity.tag==sharedFace){const auto* ids=boundaries[owner].cellNodes(i);Triangle key;
   for(unsigned j=0;j<3;++j){const auto p=boundaries[owner].node(ids[j]);key[j]={std::bit_cast<std::uint64_t>(p.x),std::bit_cast<std::uint64_t>(p.y),std::bit_cast<std::uint64_t>(p.z)};}
   int sign=1;for(unsigned a=0;a<3;++a)for(unsigned b=a+1;b<3;++b)if(key[b]<key[a]){std::swap(key[a],key[b]);sign=-sign;}auto& entry=winding[key];entry.first+=sign;++entry.second;
  }
 }
 check(!winding.empty(),"refined shared facets exist");for(const auto& [key,value]:winding)check(value.first==0&&value.second==2,"bitwise matching coowner facets and opposite winding");
 check(bytes(*shared.model)==sharedBefore,"shared source BREP unchanged");
 boundaries=saved;limited=o;limited.maximumAddedNodes=1;const auto late=refineOcctFacets(*shared.model,both,limited);
 check(!late.published&&late.reason==OcctFacetRefinementReason::BudgetExceeded&&mesh::meshHash(boundaries[0])==hash0&&mesh::meshHash(boundaries[1])==hash1,"late second-owner budget refusal rolls back both");
 noCommitted(late);
 boundaries=saved;limited=o;limited.maximumCandidatePatches=1;limited.maximumAddedNodes=1;
 const auto cappedCoowner=refineOcctFacets(*shared.model,both,limited);
 check(!cappedCoowner.published&&cappedCoowner.reason==OcctFacetRefinementReason::BudgetExceeded&&
       mesh::meshHash(boundaries[0])==hash0&&mesh::meshHash(boundaries[1])==hash1,
       "one canonical candidate still counts both actual coowner insertions against hard node budget");
 noCommitted(cappedCoowner);
 // A changed coowner triangle is a refusal, even when its surface remains
 // topologically closed and embedded. No coordinate welding is attempted.
 boundaries=saved;const auto* ids=boundaries[0].cellNodes(selected);const auto chosen=boundaries[0].node(ids[0]);
 for(Index i=0;i<boundaries[1].nodeCount();++i)if(same(boundaries[1].node(i),chosen)){auto p=boundaries[1].node(i);p.y+=.01;boundaries[1].setNode(i,p);break;}
 const auto mismatch0=mesh::meshHash(boundaries[0]),mismatch1=mesh::meshHash(boundaries[1]);const auto mismatch=refineOcctFacets(*shared.model,both,o);
 check(!mismatch.published&&std::string_view(mismatch.diagnostic)=="missing-bitwise-shared-facet"&&mesh::meshHash(boundaries[0])==mismatch0&&mesh::meshHash(boundaries[1])==mismatch1,"mismatched shared facet is not silently welded");
 noCommitted(mismatch);
 auto cylinder=OcctModelAccess::fromShape(BRepPrimAPI_MakeCylinder(1.,2.).Shape());check(cylinder.imported(),"cylinder import");mesh::Mesh round;
 check(buildOcctBoundary(*cylinder.model,1,bo,round).published,"curved cylinder boundary");const auto roundBefore=bytes(*cylinder.model);Index curved=-1;
 for(Index i=0;i<round.cellCount();++i)if(cylinder.model->face(round.cell(i).onEntity.tag)->kind()==SurfaceKind::Cylinder){curved=i;break;}
 check(curved>=0,"native cylinder facet");const auto roundOriginal=round;const std::array curvedSelection{OcctFacetSelection{1,curved}};o.selections=curvedSelection;o.maximumProjectionDistance=.1;
 const std::array roundInput{OcctFacetBoundary{1,&round}};const auto curvedResult=refineOcctFacets(*cylinder.model,roundInput,o);
 check(curvedResult.published&&curvedResult.addedNodes==1&&curvedResult.maximumDisplacement>1e-6,"curved facet insertion uses native surface rather than flat centroid");
 committed(roundOriginal,round,1,curvedResult);
 prefix(roundOriginal,round);Vec3 target;const auto cls=round.nodeClass(round.nodeCount()-1);
 check(cls.dim==2&&closestPointOnOcctFace(*cylinder.model,cls.tag,round.node(round.nodeCount()-1),target)&&(target-round.node(round.nodeCount()-1)).norm()<=cylinder.model->tolerance(2,cls.tag)&&auditOcctBoundaryFacets(*cylinder.model,1,round,bo).completed,"inserted point and curved facets pass original CAD gates");
 check(bytes(*cylinder.model)==roundBefore,"curved source unchanged");round=roundOriginal;o.maximumProjectionDistance=1e-10;const auto roundHash=mesh::meshHash(round);const auto jump=refineOcctFacets(*cylinder.model,roundInput,o);
 check(!jump.published&&std::string_view(jump.diagnostic)=="projection-displacement-limit"&&mesh::meshHash(round)==roundHash,"excess projection displacement rolls back");
 noCommitted(jump);
 exerciseConcaveGeometryRefinement();
 exerciseMixedFacetCandidates();
 std::cout<<checks<<" facet-refinement checks passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
