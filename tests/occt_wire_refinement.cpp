// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctWireRefinement.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctModelInternal.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <NCollection_List.hxx>
#include <bit>
#include <iostream>
#include <map>
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
}
int main(){try{
 auto box=OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(1.,1.,1.).Shape());check(box.imported(),"box import");
 const auto before=bytes(*box.model);OcctBoundaryOptions bo;bo.linearDeflection=.01;bo.maximumChordError=.03;
 mesh::Mesh original;check(buildOcctBoundary(*box.model,1,bo,original).published,"box boundary");
 OcctWireRefinementOptions o;o.edgeTag=1;o.segmentsPerExistingSegment=4;
 mesh::Mesh refined=original;const std::array one{OcctWireBoundary{1,&refined}};
 const auto success=refineOcctWireRefinement(*box.model,one,o);
 check(success.published&&success.addedNodes==3&&success.addedTriangles==6&&success.nativeUseChecks==6,"four-way actual wire split");
 prefix(original,refined);check(bytes(*box.model)==before,"box source BREP unchanged");
 check(auditOcctBoundaryFacets(*box.model,1,refined,bo).completed,"refined box CAD facet checks");
 mesh::Mesh volume;OcctStarOptions star;star.nativeAdmissionModel=box.model.get();star.nativeAdmissionBodyTag=1;
 const auto fan=buildOcctStarVolume(refined,1,star,volume);
 const auto volumeQuality=algo::analyzeMixedMeshQuality(volume);
 check(fan.published&&volumeQuality.valid&&volume.countCells(mesh::CellType::Tet4)==18&&
       std::abs(volumeQuality.totalVolume-1.)<1e-12&&mesh::auditLinearVolumeBoundary(volume).valid,
       "refined actual wire participates in a valid native-admitted cube volume");
 mesh::Mesh repeated=original;const std::array repeat{OcctWireBoundary{1,&repeated}};
 const auto deterministic=refineOcctWireRefinement(*box.model,repeat,o);
 check(deterministic.published&&mesh::meshHash(repeated)==mesh::meshHash(refined)&&deterministic.workUsed==success.workUsed,"deterministic mesh and work");
 const auto refuses=[&](OcctWireRefinementOptions limited,OcctWireRefinementReason reason){
  mesh::Mesh sentinel=original;const auto hash=mesh::meshHash(sentinel);const std::array item{OcctWireBoundary{1,&sentinel}};
  const auto r=refineOcctWireRefinement(*box.model,item,limited);
  check(!r.published&&r.reason==reason&&mesh::meshHash(sentinel)==hash&&bytes(*box.model)==before,"refusal preserves exact destination and source");
 };
 auto limited=o;limited.limits.maximumWork=success.workUsed-1;refuses(limited,OcctWireRefinementReason::BudgetExceeded);
 limited=o;limited.maximumAddedNodes=2;refuses(limited,OcctWireRefinementReason::BudgetExceeded);
 limited=o;limited.limits.maximumBoundaryNodes=10;refuses(limited,OcctWireRefinementReason::BudgetExceeded);
 limited=o;limited.limits.maximumWorkingBytes=1;refuses(limited,OcctWireRefinementReason::BudgetExceeded);
 limited=o;limited.edgeTag=-1;refuses(limited,OcctWireRefinementReason::InvalidRequest);
 std::atomic<bool> cancelled{true};limited=o;limited.limits.cancellation=&cancelled;refuses(limited,OcctWireRefinementReason::Cancelled);
 NCollection_List<TopoDS_Shape> args;
 args.Append(BRepPrimAPI_MakeBox(gp_Pnt(-1,0,0),1.,1.,1.).Shape());
 args.Append(BRepPrimAPI_MakeBox(gp_Pnt(0,0,0),1.,1.,1.).Shape());
 BRepAlgoAPI_Splitter splitter;splitter.SetArguments(args);splitter.SetNonDestructive(true);splitter.Build();
 check(splitter.IsDone()&&!splitter.HasErrors(),"shared planar fixture split");
 auto shared=OcctModelAccess::fromShape(splitter.Shape());check(shared.imported()&&shared.model->solids().size()==2,"both shared bodies");
 const auto sharedBefore=bytes(*shared.model);int sharedFace=-1,wire=-1;
 for(const auto* face:shared.model->faces())if(shared.model->faceOwners(face->tag()).size()==2)sharedFace=face->tag();
 check(sharedFace>0,"actual shared face authority");
 for(const auto& ring:shared.model->boundaryUses(sharedFace))for(const auto& use:ring.edges)
  if(!use.degenerated&&shared.model->edge(use.edgeTag)->kind()==CurveKind::Line){wire=use.edgeTag;break;}
 check(wire>0,"shared linear wire authority");
 OcctBoundaryOutcome prep;auto session=OcctBoundarySession::prepare(*shared.model,bo,prep);
 check(bool(session),"shared owned preparation");
 std::array<mesh::Mesh,2> boundaries;
 const std::array both{OcctWireBoundary{shared.model->solids()[0]->tag(),&boundaries[0]},OcctWireBoundary{shared.model->solids()[1]->tag(),&boundaries[1]}};
 for(unsigned i=0;i<2;++i)check(session->boundary(both[i].bodyTag,boundaries[i]).published,"shared body extraction");
 const auto saved=boundaries;const auto hash0=mesh::meshHash(boundaries[0]),hash1=mesh::meshHash(boundaries[1]);
 o.edgeTag=wire;const auto missing=refineOcctWireRefinement(*shared.model,std::span(both).first(1),o);
 check(!missing.published&&std::string_view(missing.diagnostic)=="missing-shared-face-coowner"&&mesh::meshHash(boundaries[0])==hash0&&mesh::meshHash(boundaries[1])==hash1,"missing coowner refuses atomically");
 const auto coupled=refineOcctWireRefinement(*shared.model,both,o);
 check(coupled.published&&coupled.addedNodes==6&&coupled.addedTriangles==12,"both owners split together");
 // Compare every directed shared triangle by exact point bits. Local node
 // indices may differ; coordinate evaluation and winding must agree.
 using Point=std::array<std::uint64_t,3>;using Triangle=std::array<Point,3>;
 std::map<Triangle,std::pair<int,unsigned>> winding;
 for(unsigned owner=0;owner<2;++owner){prefix(saved[owner],boundaries[owner]);
  check(auditOcctBoundaryFacets(*shared.model,both[owner].bodyTag,boundaries[owner],bo).completed,"shared CAD facet checks");
  for(Index i=0;i<boundaries[owner].cellCount();++i)if(boundaries[owner].cell(i).onEntity.tag==sharedFace){
   const auto* ids=boundaries[owner].cellNodes(i);Triangle key;
   for(unsigned j=0;j<3;++j){const auto p=boundaries[owner].node(ids[j]);key[j]={std::bit_cast<std::uint64_t>(p.x),std::bit_cast<std::uint64_t>(p.y),std::bit_cast<std::uint64_t>(p.z)};}
   int sign=1;for(unsigned a=0;a<3;++a)for(unsigned b=a+1;b<3;++b)if(key[b]<key[a]){std::swap(key[a],key[b]);sign=-sign;}
   auto& entry=winding[key];entry.first+=sign;++entry.second;
  }
 }
 check(!winding.empty(),"shared refined triangles exist");
 for(const auto& [key,value]:winding)check(value.first==0&&value.second==2,"bitwise matching shared facets have opposite winding");
 check(bytes(*shared.model)==sharedBefore,"shared exact source unchanged");
 // Refusal after private construction of the first body must preserve BOTH.
 boundaries=saved;limited=o;limited.maximumAddedNodes=4;
 const auto partial=refineOcctWireRefinement(*shared.model,both,limited);
 check(!partial.published&&partial.reason==OcctWireRefinementReason::BudgetExceeded&&mesh::meshHash(boundaries[0])==hash0&&mesh::meshHash(boundaries[1])==hash1,"late second-body failure rolls back both bodies");
 auto cylinder=OcctModelAccess::fromShape(BRepPrimAPI_MakeCylinder(1.,2.).Shape());check(cylinder.imported(),"cylinder import");
 mesh::Mesh round;check(buildOcctBoundary(*cylinder.model,1,bo,round).published,"cylinder surface");
 int circular=-1;for(const auto* edge:cylinder.model->edges())if(edge->periodic()){circular=edge->tag();break;}
 check(circular>0,"periodic circular curve found");o.edgeTag=circular;o.allowNonlinearWireCandidates=true;
 const auto roundHash=mesh::meshHash(round);const std::array roundInput{OcctWireBoundary{1,&round}};
 const auto periodic=refineOcctWireRefinement(*cylinder.model,roundInput,o);
 check(!periodic.published&&periodic.reason==OcctWireRefinementReason::UnsupportedGeometry&&mesh::meshHash(round)==roundHash,"periodic candidate preserves the original boundary");
 const auto roundOriginal=round;const auto roundBytes=bytes(*cylinder.model);o.allowClosedWireCandidates=true;o.segmentsPerExistingSegment=2;
 const auto closed=refineOcctWireRefinement(*cylinder.model,roundInput,o);
 check(closed.published&&closed.addedNodes>0&&closed.nativeUseChecks>0,"closed native wire refines across its actual seam");prefix(roundOriginal,round);
 check(bytes(*cylinder.model)==roundBytes&&auditOcctBoundaryFacets(*cylinder.model,1,round,bo).completed,"closed wire keeps original CAD and facet fidelity");
 mesh::Mesh closedRepeat=roundOriginal;const std::array closedAgain{OcctWireBoundary{1,&closedRepeat}};const auto repeatedClosed=refineOcctWireRefinement(*cylinder.model,closedAgain,o);
 check(repeatedClosed.published&&mesh::meshHash(closedRepeat)==mesh::meshHash(round)&&repeatedClosed.workUsed==closed.workUsed,"closed seam and cap subdivision are deterministic");
 auto closedLimited=o;closedLimited.limits.maximumWork=closed.workUsed-1;closedRepeat=roundOriginal;const auto closedBound=refineOcctWireRefinement(*cylinder.model,closedAgain,closedLimited);
 check(!closedBound.published&&closedBound.reason==OcctWireRefinementReason::BudgetExceeded&&mesh::meshHash(closedRepeat)==roundHash&&bytes(*cylinder.model)==roundBytes,"closed loop one-less work rolls back all changes");
 o.allowClosedWireCandidates=false;
 // A full circular trim shared by two body owners exercises both cap ears
 // and opposite face windings. The interface must receive identical facets.
 NCollection_List<TopoDS_Shape> roundArgs;roundArgs.Append(BRepPrimAPI_MakeCylinder(1.,1.).Shape());roundArgs.Append(BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0,0,1),gp_Dir(0,0,1)),1.,1.).Shape());
 BRepAlgoAPI_Splitter roundSplit;roundSplit.SetArguments(roundArgs);roundSplit.SetNonDestructive(true);roundSplit.Build();check(roundSplit.IsDone()&&!roundSplit.HasErrors(),"shared circular fixture split");
 auto joined=OcctModelAccess::fromShape(roundSplit.Shape());check(joined.imported()&&joined.model->solids().size()==2,"shared circular owners");int cap=-1,loop=-1;
 for(const auto* face:joined.model->faces())if(joined.model->faceOwners(face->tag()).size()==2)cap=face->tag();check(cap>0,"shared disk found");
 for(const auto& ring:joined.model->boundaryUses(cap))for(const auto& use:ring.edges)if(joined.model->edge(use.edgeTag)->startVertexTag()==joined.model->edge(use.edgeTag)->endVertexTag())loop=use.edgeTag;check(loop>0,"shared closed wire found");
 OcctBoundaryOutcome circularPreparation;auto circularSession=OcctBoundarySession::prepare(*joined.model,bo,circularPreparation);check(bool(circularSession),"shared circular source session");
 std::array<mesh::Mesh,2> closedBodies;for(int b=1;b<=2;++b)check(circularSession->boundary(b,closedBodies[b-1]).published,"shared circular source body");const auto circularSaved=closedBodies;
 const std::array circularOwners{OcctWireBoundary{1,&closedBodies[0]},OcctWireBoundary{2,&closedBodies[1]}};o.edgeTag=loop;o.allowClosedWireCandidates=true;
 const auto sharedClosed=refineOcctWireRefinement(*joined.model,circularOwners,o);check(sharedClosed.published,"closed curve updates both owners atomically");
 using PointKey=std::array<std::uint64_t,3>;using TriKey=std::array<PointKey,3>;std::map<TriKey,int> closedWinding;std::map<TriKey,unsigned> closedCounts;
 for(unsigned b=0;b<2;++b){prefix(circularSaved[b],closedBodies[b]);for(Index i=0;i<closedBodies[b].cellCount();++i)if(closedBodies[b].cell(i).onEntity.tag==cap){const auto* ids=closedBodies[b].cellNodes(i);TriKey key;for(unsigned j=0;j<3;++j){const auto p=closedBodies[b].node(ids[j]);key[j]={std::bit_cast<std::uint64_t>(p.x),std::bit_cast<std::uint64_t>(p.y),std::bit_cast<std::uint64_t>(p.z)};}int sign=1;for(unsigned j=0;j<3;++j)for(unsigned k=j+1;k<3;++k)if(key[k]<key[j]){std::swap(key[k],key[j]);sign=-sign;}closedWinding[key]+=sign;++closedCounts[key];}}
 check(!closedWinding.empty(),"shared closed interface retains triangles");for(const auto& [key,sign]:closedWinding)check(sign==0&&closedCounts.at(key)==2,"shared disk has exactly matching triangles and opposite winding");
 o.allowClosedWireCandidates=false;
 auto sector=OcctModelAccess::fromShape(BRepPrimAPI_MakeCylinder(1.,2.,1.5707963267948966).Shape());check(sector.imported(),"quarter cylinder import");
 mesh::Mesh arc;check(buildOcctBoundary(*sector.model,1,bo,arc).published,"quarter cylinder boundary");const auto unrefined=arc;const auto sectorBytes=bytes(*sector.model);int trimmedArc=-1;
 for(const auto* edge:sector.model->edges())if(edge->periodic()&&edge->startVertexTag()!=edge->endVertexTag()){trimmedArc=edge->tag();break;}
 check(trimmedArc>0,"trimmed periodic arc has distinct endpoints");o.edgeTag=trimmedArc;const std::array arcInput{OcctWireBoundary{1,&arc}};const auto arcResult=refineOcctWireRefinement(*sector.model,arcInput,o);
 check(arcResult.published&&arcResult.addedNodes>0&&arcResult.nativeUseChecks>0,"bounded periodic arc refinement");prefix(unrefined,arc);check(bytes(*sector.model)==sectorBytes,"trimmed arc source unchanged");check(auditOcctBoundaryFacets(*sector.model,1,arc,bo).completed,"trimmed arc CAD facet checks");
 std::cout<<checks<<" wire-refinement checks passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
