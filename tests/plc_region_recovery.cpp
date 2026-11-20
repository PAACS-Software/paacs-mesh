// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/bnd/TetPlcCavityRecovery.hpp"
#include "paacs/bnd/TetFacetGeometry.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <iostream>
#include <map>
#include <cstdlib>
#include <new>
#include <stdexcept>
#include "allocation_nothrow.hpp"
bool countAllocations=false;std::size_t allocationCount=0,faultOrdinal=SIZE_MAX;
void* operator new(std::size_t size){if(countAllocations&&allocationCount++==faultOrdinal)throw std::bad_alloc();if(auto* p=std::malloc(size?size:1))return p;throw std::bad_alloc();}
void operator delete(void* p)noexcept{std::free(p);}void operator delete(void* p,std::size_t)noexcept{std::free(p);}
using namespace paacs;
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
std::vector<std::uint8_t> bytes(const mesh::Mesh& m){std::vector<std::uint8_t> out;check(sdk::serializeMeshArchive(m,out).completed,"archive");return out;}
struct Fixture {mesh::Mesh source,partial;};
Fixture concave(bool keep=false){
 std::map<std::array<int,3>,unsigned> nodes;std::vector<bnd::TetFacetVertex> vertices;std::vector<bnd::TetFacetCell> cells;Fixture f;
 const std::array<std::array<int,3>,8> offsets{{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}}};
 for(int y=0;y<2;++y)for(int x=0;x<3;++x){if(x==1&&y==1)continue;std::array<unsigned,8> ids;
  for(unsigned i=0;i<8;++i){const auto o=offsets[i];const std::array<int,3> key{x+o[0],y+o[1],o[2]};auto [at,added]=nodes.emplace(key,unsigned(vertices.size()));
   if(added)vertices.push_back({{double(key[0]),double(key[1]),double(key[2])},vertices.size()+1});ids[i]=at->second;}
  for(auto local:std::array<std::array<unsigned,4>,6>{{{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}}}){
   std::array<unsigned,4> c;for(unsigned j=0;j<4;++j)c[j]=ids[local[j]];cells.push_back({c,cells.size()+1,{3,7}});}}
 for(auto v:vertices){f.source.addNode(v.position,{0,int(v.sequence)});f.partial.addNode(v.position,{0,int(v.sequence)});}
 bnd::TetFacetComplexView view{vertices,cells,{},{},{},{3,7}};const auto faces=bnd::tet_facet_detail::candidateFaceMap(view,cells);int tag=1;
 for(const auto& [key,uses]:faces){(void)key;if(uses.size()==1){std::array<Index,3> ids;for(unsigned j=0;j<3;++j)ids[j]=uses[0].directed[j];f.source.addCell(mesh::CellType::Tri3,ids,{2,tag++});}}
 if(keep)for(unsigned i=0;i<6;++i){std::array<Index,4> ids;for(unsigned j=0;j<4;++j)ids[j]=cells[i].vertices[j];f.partial.addCell(mesh::CellType::Tet4,ids,{3,7});}
 return f;
}
Fixture blockedGroup(bool backtracking=false){
 const auto base=concave();Fixture f;std::vector<Index> inverse(base.source.nodeCount());
 // Captured seed-2 advancing-front witness. Freeze coordinates/connectivity;
 // the regression does not depend on a library's random shuffle algorithm.
 const std::array<Vec3,24> first{{{2,1,0},{0,2,0},{0,1,1},{3,2,1},{1,2,1},{3,1,1},{0,2,1},{1,0,1},{0,0,1},{1,1,1},{3,0,1},{2,0,1},{2,0,0},{2,2,1},{1,1,0},{3,2,0},{2,1,1},{0,1,0},{3,0,0},{1,2,0},{2,2,0},{1,0,0},{0,0,0},{3,1,0}}};
 const std::array<Vec3,24> second{{{0,2,1},{1,2,1},{2,1,0},{1,0,1},{2,0,0},{3,1,0},{3,0,0},{0,1,1},{0,1,0},{1,0,0},{2,2,0},{0,2,0},{1,2,0},{3,0,1},{2,2,1},{0,0,1},{1,1,0},{0,0,0},{2,0,1},{3,2,1},{3,1,1},{2,1,1},{1,1,1},{3,2,0}}};
 const auto& points=backtracking?second:first;
 for(auto p:points){Index old=kInvalidIndex;for(Index i=0;i<base.source.nodeCount();++i){auto q=base.source.node(i);if(q.x==p.x&&q.y==p.y&&q.z==p.z){old=i;break;}}check(old!=kInvalidIndex,"frozen witness node");inverse[old]=f.source.nodeCount();f.source.addNode(p,base.source.nodeClass(old));f.partial.addNode(p,base.source.nodeClass(old));}
 for(Index i=0;i<base.source.cellCount();++i){std::array<Index,3> ids;for(unsigned j=0;j<3;++j)ids[j]=inverse[base.source.cellNodes(i)[j]];f.source.addCell(mesh::CellType::Tri3,ids,base.source.cell(i).onEntity);}
 if(backtracking){for(auto ids:std::array<std::array<Index,4>,8>{{{0,1,7,8},{0,1,8,11},{1,7,8,16},{1,7,16,22},{1,8,11,12},{1,8,12,16},{4,5,2,18},{9,4,2,18}}})f.partial.addCell(mesh::CellType::Tet4,ids,{3,7});}
 else for(auto ids:std::array<std::array<Index,4>,6>{{{21,12,0,11},{12,0,11,5},{0,11,5,16},{12,0,5,23},{0,16,5,3},{0,3,5,23}}})f.partial.addCell(mesh::CellType::Tet4,ids,{3,7});
 return f;
}
void exactAudit(const Fixture& f,const mesh::Mesh& mesh){
 check(mesh::auditLinearVolumeBoundary(mesh).valid,"boundary incidence");
 std::vector<bnd::TetFacetVertex> vertices;std::vector<bnd::TetFacetCell> cells;
 for(Index i=0;i<mesh.nodeCount();++i)vertices.push_back({mesh.node(i),std::uint64_t(i)+1});
 for(Index i=f.source.cellCount();i<mesh.cellCount();++i){std::array<unsigned,4> ids;for(unsigned j=0;j<4;++j)ids[j]=unsigned(mesh.cellNodes(i)[j]);cells.push_back({ids,cells.size()+1,{3,7}});}
 bnd::TetFacetComplexView view{vertices,cells,{},{},{},{3,7}};
 for(std::size_t i=0;i<cells.size();++i){const auto c=cells[i];check(predicates::orient3d(vertices[c.vertices[0]].position,vertices[c.vertices[1]].position,vertices[c.vertices[2]].position,vertices[c.vertices[3]].position)==predicates::Sign::Positive,"orientation");
  for(std::size_t j=0;j<i;++j)check(!bnd::tet_facet_detail::tetImproperlyIntersects(view,c,cells[j]),"tet overlap");
  std::array<Vec3,4> corners;for(unsigned j=0;j<4;++j)corners[j]=vertices[c.vertices[j]].position;check(!bnd::tet_facet_detail::pointStrictlyInsideTet({1.5,1.5,.5},corners),"physical gap filled");}
 for(Index i=0;i<f.source.cellCount();++i){check(mesh.cell(i).onEntity==f.source.cell(i).onEntity,"face identity");for(unsigned j=0;j<3;++j)check(mesh.cellNodes(i)[j]==f.source.cellNodes(i)[j],"source triangles changed");}
}
int main(int argc,char** argv){try{
 unsigned passed=0;auto f=concave();const auto source=bytes(f.source),partial=bytes(f.partial);
 bnd::TetPlcCavityOptions options;options.region={3,7};options.maximumPatches=1;mesh::Mesh result=f.source;
 auto old=bnd::recoverTetPlcCavities(f.source,f.partial,options,result);check(!old.published&&bytes(result)==source,"single-apex unexpectedly fills U");++passed;
 options.patchSearch.maximumStates=4096;options.maximumNodes=f.source.nodeCount();
 const auto out=bnd::recoverTetPlcCavities(f.source,f.partial,options,result);
 std::cout<<"{\"stage\":\"concave-region\",\"published\":"<<out.published<<",\"diagnostic\":\""<<out.diagnostic<<"\",\"states\":"<<out.regionSearchStates<<",\"backtracks\":"<<out.regionSearchBacktracks<<",\"candidates\":"<<out.regionSearchCandidates<<",\"cache_hits\":"<<out.regionSearchPairCacheHits<<",\"bound_hits\":"<<out.regionSearchBoundHits<<",\"work\":"<<out.workUsed<<",\"search_seconds\":"<<out.regionSearchSeconds<<",\"preparation_seconds\":"<<out.preparationSeconds<<",\"publication_seconds\":"<<out.publicationSeconds<<",\"validation_seconds\":"<<out.validationSeconds<<"}\n"<<std::flush;
 check(out.published&&out.regionSearchPatches==1&&result.nodeCount()==f.source.nodeCount()&&std::abs(out.meshVolume-5)<1e-11,"concave region completion");exactAudit(f,result);check(bytes(f.source)==source&&bytes(f.partial)==partial,"source mutation");++passed;
 mesh::Mesh again;const auto repeat=bnd::recoverTetPlcCavities(f.source,f.partial,options,again);check(repeat.published&&bytes(result)==bytes(again),"repeatable local region search");++passed;
 auto limited=options;limited.patchSearch.maximumStates=1;mesh::Mesh refused=f.source;auto no=bnd::recoverTetPlcCavities(f.source,f.partial,limited,refused);check(!no.published&&no.reason==bnd::TetPlcCavityReason::BudgetExceeded&&bytes(refused)==source&&no.regionSearchBoundHits,"search allowance");++passed;
 limited=options;limited.patchSearch.maximumCandidateTets=1;no=bnd::recoverTetPlcCavities(f.source,f.partial,limited,refused);check(!no.published&&no.reason==bnd::TetPlcCavityReason::BudgetExceeded&&bytes(refused)==source,"candidate allowance");++passed;
 limited=options;limited.patchSearch.maximumVertices=4;no=bnd::recoverTetPlcCavities(f.source,f.partial,limited,refused);check(!no.published&&no.reason==bnd::TetPlcCavityReason::BudgetExceeded&&bytes(refused)==source,"vertex allowance");++passed;
 limited=options;limited.interiorAdmission=[](Vec3){return false;};no=bnd::recoverTetPlcCavities(f.source,f.partial,limited,refused);check(!no.published&&bytes(refused)==source,"material predicate bypass");++passed;
 limited=options;std::stop_source stop;stop.request_stop();limited.cancellation=stop.get_token();no=bnd::recoverTetPlcCavities(f.source,f.partial,limited,refused);check(!no.published&&no.reason==bnd::TetPlcCavityReason::Cancelled&&bytes(refused)==source,"cancelled region publication");++passed;
 {auto larger=options;larger.patchSearch.maximumVertices=512;mesh::Mesh rebuilt=f.source;const auto recovered=bnd::recoverTetPlcCavities(f.source,f.partial,larger,rebuilt);check(recovered.published,"larger declared region allowance");exactAudit(f,rebuilt);check(bytes(f.source)==source&&bytes(f.partial)==partial,"larger allowance source mutation");++passed;}
 limited=options;limited.patchSearch.maximumVertices=513;no=bnd::recoverTetPlcCavities(f.source,f.partial,limited,refused);check(!no.published&&no.reason==bnd::TetPlcCavityReason::InvalidRequest&&bytes(refused)==source,"invalid region limits");++passed;
 {auto group=blockedGroup();const auto beforeSource=bytes(group.source),beforePartial=bytes(group.partial);auto repair=options;repair.maximumPatches=32;repair.useVisibilityExpansion=false;repair.patchSearch.maximumShells=0;mesh::Mesh rebuilt=group.source;
  // Separating angular sectors can now close this former no-shell refusal.
  // Qualify either outcome against the original geometry, not a search trace.
  auto blocked=bnd::recoverTetPlcCavities(group.source,group.partial,repair,rebuilt);
  if(blocked.published)exactAudit(group,rebuilt);
  else check(blocked.reason==bnd::TetPlcCavityReason::Unavailable&&bytes(rebuilt)==beforeSource,"fixed chamber refusal rollback");
  repair.patchSearch.maximumShells=3;auto recovered=bnd::recoverTetPlcCavities(group.source,group.partial,repair,rebuilt);if(!recovered.published)std::cerr<<"group refill "<<unsigned(recovered.reason)<<' '<<recovered.diagnostic<<" removed "<<recovered.removedCells<<" added "<<recovered.addedCells<<'\n';check(recovered.published&&recovered.regionSearchPatches>=1&&rebuilt.nodeCount()<=repair.maximumNodes,"connected group replacement");exactAudit(group,rebuilt);check(bytes(group.source)==beforeSource&&bytes(group.partial)==beforePartial,"group recovery source mutation");++passed;
  std::cout<<"{\"stage\":\"connected-cell-repair\",\"removed\":"<<recovered.removedCells<<",\"added\":"<<recovered.addedCells<<",\"search_attempts\":"<<recovered.regionSearchAttempts<<",\"states\":"<<recovered.regionSearchStates<<",\"search_seconds\":"<<recovered.regionSearchSeconds<<"}\n";
  if(argc==3){const std::string prefix=argv[2];check(sdk::writeMeshArchive(group.source,prefix+"-source.paacsmesh").completed&&sdk::writeMeshArchive(group.partial,prefix+"-partial.paacsmesh").completed&&sdk::writeMeshArchive(rebuilt,prefix+"-result.paacsmesh").completed,"export group witness");}
 }
 {auto group=blockedGroup(true);auto repair=options;repair.maximumPatches=32;repair.useVisibilityExpansion=false;repair.patchSearch.maximumShells=0;mesh::Mesh rebuilt=group.source;const auto before=bytes(rebuilt);const auto refused=bnd::recoverTetPlcCavities(group.source,group.partial,repair,rebuilt);if(refused.published)exactAudit(group,rebuilt);else check(bytes(rebuilt)==before,"backtracking refusal rollback");repair.patchSearch.maximumShells=3;const auto out=bnd::recoverTetPlcCavities(group.source,group.partial,repair,rebuilt);check(out.published&&rebuilt.nodeCount()<=repair.maximumNodes,"region recovery with exact contacts");exactAudit(group,rebuilt);++passed;
  std::cout<<"{\"stage\":\"backtracking-repair\",\"backtracks\":"<<out.regionSearchBacktracks+refused.regionSearchBacktracks<<",\"pair_cache_hits\":"<<out.regionSearchPairCacheHits+refused.regionSearchPairCacheHits<<"}\n";
 }
 {limited=options;std::stop_source stop;unsigned queries=0;limited.cancellation=stop.get_token();limited.interiorAdmission=[&](Vec3){if(++queries==8)stop.request_stop();return true;};no=bnd::recoverTetPlcCavities(f.source,f.partial,limited,refused);check(!no.published&&no.reason==bnd::TetPlcCavityReason::Cancelled&&bytes(refused)==source&&queries>=8,"cancel during patch search");++passed;}
 {mesh::Mesh healthy;allocationCount=0;faultOrdinal=SIZE_MAX;countAllocations=true;const auto healthyOut=bnd::recoverTetPlcCavities(f.source,f.partial,options,healthy);countAllocations=false;const auto total=allocationCount;check(healthyOut.published&&total>32,"allocation campaign setup");
  const std::set<std::size_t> samples{0,1,2,total/8,total/4,total/2,3*total/4,7*total/8,total-3,total-2,total-1};
  for(auto ordinal:samples){auto target=f.source;allocationCount=0;faultOrdinal=ordinal;countAllocations=true;const auto failure=bnd::recoverTetPlcCavities(f.source,f.partial,options,target);countAllocations=false;check(!failure.published&&failure.reason==bnd::TetPlcCavityReason::AllocationFailed&&bytes(target)==source,"sampled allocation rollback");}
  faultOrdinal=SIZE_MAX;std::cout<<"{\"allocation_sites\":"<<total<<",\"sampled_faults\":"<<samples.size()<<",\"all_unchanged\":true}\n";++passed;
 }
 if(argc==2)check(sdk::writeMeshArchive(result,argv[1]).completed,"write demonstration mesh");
 std::cout<<"{\"passed\":"<<passed<<",\"failed\":0}\n";return 0;
 }catch(const std::exception& e){countAllocations=false;std::cerr<<e.what()<<'\n';return 1;}}
