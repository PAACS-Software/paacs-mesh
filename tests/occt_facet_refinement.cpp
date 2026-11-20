// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctFacetRefinement.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctModelInternal.hpp"
#include "OcctRecoveryMesher.hpp"
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
 const auto before=bytes(*box.model);OcctBoundaryOptions bo;bo.linearDeflection=.03;bo.maximumChordError=.05;
 mesh::Mesh original;check(buildOcctBoundary(*box.model,1,bo,original).published,"box boundary");
 const std::array selections{OcctFacetSelection{1,0},OcctFacetSelection{1,1},OcctFacetSelection{1,0}};
 OcctFacetRefinementOptions o;o.selections=selections;o.maximumProjectionDistance=.1;
 mesh::Mesh refined=original;const std::array one{OcctFacetBoundary{1,&refined}};
 const auto success=refineOcctFacets(*box.model,one,o);
 check(success.published&&success.addedNodes==2&&success.addedTriangles==4&&success.nativeQueries==4&&success.maximumDisplacement==0,"native planar facet split and duplicate selection deduplication");
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
 const auto coupled=refineOcctFacets(*shared.model,both,o);check(coupled.published&&coupled.addedNodes==2&&coupled.addedTriangles==4,"both owners split together");
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
 // A changed coowner triangle is a refusal, even when its surface remains
 // topologically closed and embedded. No coordinate welding is attempted.
 boundaries=saved;const auto* ids=boundaries[0].cellNodes(selected);const auto chosen=boundaries[0].node(ids[0]);
 for(Index i=0;i<boundaries[1].nodeCount();++i)if(same(boundaries[1].node(i),chosen)){auto p=boundaries[1].node(i);p.y+=.01;boundaries[1].setNode(i,p);break;}
 const auto mismatch0=mesh::meshHash(boundaries[0]),mismatch1=mesh::meshHash(boundaries[1]);const auto mismatch=refineOcctFacets(*shared.model,both,o);
 check(!mismatch.published&&std::string_view(mismatch.diagnostic)=="missing-bitwise-shared-facet"&&mesh::meshHash(boundaries[0])==mismatch0&&mesh::meshHash(boundaries[1])==mismatch1,"mismatched shared facet is not silently welded");
 auto cylinder=OcctModelAccess::fromShape(BRepPrimAPI_MakeCylinder(1.,2.).Shape());check(cylinder.imported(),"cylinder import");mesh::Mesh round;
 check(buildOcctBoundary(*cylinder.model,1,bo,round).published,"curved cylinder boundary");const auto roundBefore=bytes(*cylinder.model);Index curved=-1;
 for(Index i=0;i<round.cellCount();++i)if(cylinder.model->face(round.cell(i).onEntity.tag)->kind()==SurfaceKind::Cylinder){curved=i;break;}
 check(curved>=0,"native cylinder facet");const auto roundOriginal=round;const std::array curvedSelection{OcctFacetSelection{1,curved}};o.selections=curvedSelection;o.maximumProjectionDistance=.1;
 const std::array roundInput{OcctFacetBoundary{1,&round}};const auto curvedResult=refineOcctFacets(*cylinder.model,roundInput,o);
 check(curvedResult.published&&curvedResult.addedNodes==1&&curvedResult.maximumDisplacement>1e-6,"curved facet insertion uses native surface rather than flat centroid");
 prefix(roundOriginal,round);Vec3 target;const auto cls=round.nodeClass(round.nodeCount()-1);
 check(cls.dim==2&&closestPointOnOcctFace(*cylinder.model,cls.tag,round.node(round.nodeCount()-1),target)&&(target-round.node(round.nodeCount()-1)).norm()<=cylinder.model->tolerance(2,cls.tag)&&auditOcctBoundaryFacets(*cylinder.model,1,round,bo).completed,"inserted point and curved facets pass original CAD gates");
 check(bytes(*cylinder.model)==roundBefore,"curved source unchanged");round=roundOriginal;o.maximumProjectionDistance=1e-10;const auto roundHash=mesh::meshHash(round);const auto jump=refineOcctFacets(*cylinder.model,roundInput,o);
 check(!jump.published&&std::string_view(jump.diagnostic)=="projection-displacement-limit"&&mesh::meshHash(round)==roundHash,"excess projection displacement rolls back");
 std::cout<<checks<<" facet-refinement checks passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
