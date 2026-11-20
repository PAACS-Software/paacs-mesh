// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/bnd/TetPlcCavityRecovery.hpp"
#include "paacs/bnd/TetFacetGeometry.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <map>
#include <iostream>
#include <stdexcept>
#include <atomic>
#include <cstdlib>
#include <new>
#include "allocation_nothrow.hpp"
std::atomic_bool failAllocation{false};
void* operator new(std::size_t n){if(failAllocation.exchange(false))throw std::bad_alloc();if(auto* p=std::malloc(n?n:1))return p;throw std::bad_alloc();}
void operator delete(void* p)noexcept{std::free(p);}void operator delete(void* p,std::size_t)noexcept{std::free(p);}
using namespace paacs;using bnd::TetPlcCavityReason;
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
std::vector<std::uint8_t> bytes(const mesh::Mesh& m){std::vector<std::uint8_t> result;check(sdk::serializeMeshArchive(m,result).completed,"fixture encoding");return result;}
struct Fixture {mesh::Mesh source,partial;std::vector<bnd::TetFacetCell> all;};
Fixture fixture(std::vector<Vec3> points,std::vector<std::array<unsigned,4>> cells,std::vector<unsigned> keep={}){
 Fixture f;std::vector<bnd::TetFacetVertex> vertices;for(auto p:points){const auto id=f.source.nodeCount();f.source.addNode(p,{0,static_cast<int>(id)+1});f.partial.addNode(p,{0,static_cast<int>(id)+1});vertices.push_back({p,vertices.size()+1});}
 for(auto ids:cells)f.all.push_back({ids,f.all.size()+1,{3,9}});bnd::TetFacetComplexView view{vertices,f.all,{},{},{},{3,9}};auto faces=bnd::tet_facet_detail::candidateFaceMap(view,f.all);int tag=100;
 for(const auto& [key,rows]:faces){(void)key;if(rows.size()!=1)continue;std::array<Index,3> ids;for(unsigned i=0;i<3;++i)ids[i]=rows[0].directed[i];f.source.addCell(mesh::CellType::Tri3,ids,{2,tag++});}
 for(auto i:keep){std::array<Index,4> ids;for(unsigned j=0;j<4;++j)ids[j]=cells.at(i)[j];f.partial.addCell(mesh::CellType::Tet4,ids,{3,9});}return f;
}
Fixture tetra(){return fixture({{0,0,0},{1,0,0},{0,1,0},{0,0,1}},{{0,1,2,3}});}
Fixture brick(bool uShape=false){std::map<std::array<int,3>,unsigned> nodes;std::vector<Vec3> points;std::vector<std::array<unsigned,4>> cells;std::vector<unsigned> keep;
 for(int y=0;y<2;++y)for(int x=0;x<3;++x){if(uShape&&x==1&&y==1)continue;std::array<unsigned,8> ids;const std::array<std::array<int,3>,8> offsets{{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}}};for(unsigned i=0;i<8;++i){const auto o=offsets[i];std::array<int,3> key{x+o[0],y+o[1],o[2]};auto [at,inserted]=nodes.emplace(key,unsigned(points.size()));if(inserted)points.push_back({double(key[0]),double(key[1]),double(key[2])});ids[i]=at->second;}
  for(auto local:std::array<std::array<unsigned,4>,6>{{{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}}}){std::array<unsigned,4> c;for(unsigned j=0;j<4;++j)c[j]=ids[local[j]];if((!uShape&&x==1&&y==1)||(uShape&&!(x==0&&y==0)))keep.push_back(unsigned(cells.size()));cells.push_back(c);}}
 return fixture(points,cells,keep);
}
void verify(const Fixture& f,const mesh::Mesh& result,const bnd::TetPlcCavityOutcome& out,double volume){check(out.published&&out.reason==TetPlcCavityReason::None&&!out.destinationUnchanged,"completion");check(mesh::auditLinearVolumeBoundary(result).valid,"completed incidence");check(std::abs(out.meshVolume-volume)<1e-11,"completed volume");
 for(Index i=0;i<f.source.nodeCount();++i){auto a=f.source.node(i),b=result.node(i);check(a.x==b.x&&a.y==b.y&&a.z==b.z&&f.source.nodeClass(i)==result.nodeClass(i),"source node identity/class changed");}
 for(Index i=0;i<f.source.cellCount();++i){check(f.source.cell(i).onEntity==result.cell(i).onEntity,"face classification changed");for(unsigned j=0;j<3;++j)check(f.source.cellNodes(i)[j]==result.cellNodes(i)[j],"oriented source face changed");}
 for(Index i=f.source.cellCount();i<result.cellCount();++i)check(result.cell(i).onEntity==mesh::Classification{3,9},"body identity changed");
}
Fixture ringBridge(){std::map<std::array<int,3>,unsigned> nodes;std::vector<Vec3> points;std::vector<std::array<unsigned,4>> cells;std::vector<unsigned> keep;
 for(int y=0;y<5;++y)for(int x=0;x<5;++x){std::array<unsigned,8> ids;const std::array<std::array<int,3>,8> offsets{{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}}};for(unsigned i=0;i<8;++i){const auto o=offsets[i];std::array<int,3> key{x+o[0],y+o[1],o[2]};auto [at,inserted]=nodes.emplace(key,unsigned(points.size()));if(inserted)points.push_back({double(key[0]),double(key[1]),double(key[2])});ids[i]=at->second;}
  for(auto local:std::array<std::array<unsigned,4>,6>{{{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}}}){std::array<unsigned,4> c;for(unsigned j=0;j<4;++j)c[j]=ids[local[j]];if(x>0&&x<4&&y>0&&y<4&&!(x==2&&y==2))keep.push_back(unsigned(cells.size()));cells.push_back(c);}}
 return fixture(points,cells,keep);
}
int main(){try{unsigned passed=0;bnd::TetPlcCavityOptions options;options.region={3,9};
 {auto f=fixture({{0,0,0},{1,0,0},{0,1,0},{0,0,1},{0,0,-1}},{{0,1,2,3},{0,2,1,4}});
  auto local=options;local.maximumLocalFrontTets=2;local.useVisibilityExpansion=false;
  mesh::Mesh result;auto out=bnd::recoverTetPlcCavities(f.source,f.partial,local,result);
  verify(f,result,out,1./3);check(out.localFrontTets==2&&out.localFrontClosed==1&&out.removedCells==0,"three-face advance followed by four-face closure");
  const auto before=bytes(result);local.centroidAdmission=[](Vec3 p){return p.z>=0;};
  out=bnd::recoverTetPlcCavities(f.source,f.partial,local,result);
  check(!out.published&&out.localFrontTets==1&&bytes(result)==before,"local front respects centroid admission and transactional publication");++passed;}
 {auto f=fixture({{0,0,-1},{0,0,1},{1,0,0},{0,1,0},{-1,0,0},{0,-1,0}},
                 {{0,1,2,3},{0,1,3,4},{0,1,4,5},{0,1,5,2}});
  auto local=options;local.maximumLocalFrontTets=1;local.useVisibilityExpansion=false;
  mesh::Mesh result=f.source;const auto before=bytes(result);
  auto out=bnd::recoverTetPlcCavities(f.source,f.partial,local,result);
  check(!out.published&&out.localFrontTets==0&&bytes(result)==before,"two-face advance requires a pair budget");
  local.maximumLocalFrontTets=4;local.centroidAdmission=[](Vec3 p){return p.x>0&&p.y>0;};
  out=bnd::recoverTetPlcCavities(f.source,f.partial,local,result);
  check(!out.published&&out.localFrontTets==0&&bytes(result)==before,"refused second tet prevents a partial pair commit");
  local.centroidAdmission={};out=bnd::recoverTetPlcCavities(f.source,f.partial,local,result);
  verify(f,result,out,4./3);check(out.localFrontTets==4&&out.localFrontClosed==1,"two-face pair unlocks the octahedron");++passed;}
 {auto f=tetra();mesh::Mesh empty,result;auto out=bnd::recoverTetPlcCavities(f.source,empty,options,result);verify(f,result,out,1./6);++passed;}
 {auto f=brick();mesh::Mesh result;const auto before=bytes(f.partial);auto out=bnd::recoverTetPlcCavities(f.source,f.partial,options,result);verify(f,result,out,6);check(out.removedCells>0&&bytes(f.partial)==before,"targeted expansion was not used transactionally");++passed;}
 {auto f=brick(true);mesh::Mesh result;auto out=bnd::recoverTetPlcCavities(f.source,f.partial,options,result);verify(f,result,out,5);for(Index i=f.source.cellCount();i<result.cellCount();++i){std::array<Vec3,4> p;for(unsigned j=0;j<4;++j)p[j]=result.node(result.cellNodes(i)[j]);check(!bnd::tet_facet_detail::pointStrictlyInsideTet({1.5,1.5,.5},p),"filled physical U-shaped gap");}++passed;}
 {auto f=fixture({{0,0,0},{1,0,0},{0,1,0},{0,0,1},{4,0,0},{5,0,0},{4,1,0},{4,0,1}},{{0,1,2,3},{4,5,6,7}});mesh::Mesh empty,result;auto out=bnd::recoverTetPlcCavities(f.source,empty,options,result);verify(f,result,out,1./3);++passed;}
 {auto f=tetra();mesh::Mesh open;for(Index i=0;i<f.source.nodeCount();++i)open.addNode(f.source.node(i),f.source.nodeClass(i));for(Index i=1;i<f.source.cellCount();++i)open.addCell(mesh::CellType::Tri3,{f.source.cellNodes(i),3},f.source.cell(i).onEntity);mesh::Mesh result=f.source;const auto before=bytes(result);auto out=bnd::recoverTetPlcCavities(open,f.partial,options,result);check(!out.published&&out.reason==TetPlcCavityReason::InvalidSource&&bytes(result)==before,"open source refusal");++passed;}
 {auto f=tetra();std::array<Index,4> ids{0,1,2,3};f.partial.addCell(mesh::CellType::Tet4,ids,{3,9});f.partial.addCell(mesh::CellType::Tet4,ids,{3,9});mesh::Mesh result=f.source;const auto before=bytes(result);auto out=bnd::recoverTetPlcCavities(f.source,f.partial,options,result);check(!out.published&&out.reason==TetPlcCavityReason::InvalidPartialMesh&&bytes(result)==before,"overlapping partial input");++passed;}
 {auto f=tetra();f.partial.addCell(mesh::CellType::Tet4,{0,1,2,3},{3,10});mesh::Mesh result;auto out=bnd::recoverTetPlcCavities(f.source,f.partial,options,result);check(!out.published&&out.reason==TetPlcCavityReason::InvalidPartialMesh,"joined differing region tags");++passed;}
 {auto f=tetra();auto result=f.source;const auto before=bytes(result);auto limited=options;limited.maximumWork=1;auto out=bnd::recoverTetPlcCavities(f.source,f.partial,limited,result);check(!out.published&&out.reason==TetPlcCavityReason::BudgetExceeded&&bytes(result)==before,"work budget publication");limited=options;limited.maximumWorkingBytes=1;out=bnd::recoverTetPlcCavities(f.source,f.partial,limited,result);check(!out.published&&out.reason==TetPlcCavityReason::BudgetExceeded&&bytes(result)==before,"storage budget publication");++passed;}
 {auto f=tetra();auto result=f.source;const auto before=bytes(result);auto limited=options;std::stop_source stop;stop.request_stop();limited.cancellation=stop.get_token();auto out=bnd::recoverTetPlcCavities(f.source,f.partial,limited,result);check(!out.published&&out.reason==TetPlcCavityReason::Cancelled&&bytes(result)==before,"cancelled publication");++passed;}
 {auto f=tetra();const auto before=bytes(f.source);auto out=bnd::recoverTetPlcCavities(f.source,f.partial,options,f.source);check(!out.published&&out.reason==TetPlcCavityReason::InvalidRequest&&bytes(f.source)==before,"aliased publication");++passed;}
 {auto f=brick();auto result=f.source;const auto before=bytes(result);auto limited=options;limited.maximumExpansionCells=0;auto out=bnd::recoverTetPlcCavities(f.source,f.partial,limited,result);check(!out.published&&out.reason==TetPlcCavityReason::BudgetExceeded&&bytes(result)==before,"expansion budget publication");++passed;}
 {auto f=tetra();mesh::Mesh empty,result;auto warning=options;warning.warningQuality=1;auto out=bnd::recoverTetPlcCavities(f.source,empty,warning,result);verify(f,result,out,1./6);check(out.qualityWarning,"quality warning was discarded");++passed;}
 {auto f=tetra();mesh::Mesh result=f.source;const auto before=bytes(result);failAllocation.store(true);auto out=bnd::recoverTetPlcCavities(f.source,f.partial,options,result);failAllocation.store(false);check(!out.published&&out.reason==TetPlcCavityReason::AllocationFailed&&bytes(result)==before,"allocation failure publication");++passed;}
 {auto f=tetra();mesh::Mesh empty,result=f.source;const auto before=bytes(result);auto admitted=options;admitted.interiorAdmission=[](Vec3){return false;};auto out=bnd::recoverTetPlcCavities(f.source,empty,admitted,result);check(!out.published&&out.reason==TetPlcCavityReason::Unavailable&&bytes(result)==before&&!out.unresolvedSourceTriangles.empty(),"external material admission bypassed");++passed;}
 {auto f=tetra();f.partial.addCell(mesh::CellType::Tet4,{0,1,2,3},{3,9});mesh::Mesh result=f.source;const auto before=bytes(result);auto admitted=options;admitted.interiorAdmission=[](Vec3){return false;};auto out=bnd::recoverTetPlcCavities(f.source,f.partial,admitted,result);check(!out.published&&out.reason==TetPlcCavityReason::InvalidPartialMesh&&bytes(result)==before,"invalid partial material admission");++passed;}
 {auto f=tetra();mesh::Mesh empty,result=f.source;const auto before=bytes(result);auto admitted=options;admitted.interiorAdmission=[](Vec3)->bool{throw std::runtime_error("provider refused");};auto out=bnd::recoverTetPlcCavities(f.source,empty,admitted,result);check(!out.published&&out.reason==TetPlcCavityReason::InternalError&&bytes(result)==before,"throwing material provider publication");++passed;}
 {auto f=tetra();mesh::Mesh empty,refused=f.source;std::vector<Vec3> oldPoints;auto baseline=options;baseline.maximumCandidatePoints=128;baseline.interiorAdmission=[&](Vec3 p){oldPoints.push_back(p);return false;};auto no=bnd::recoverTetPlcCavities(f.source,empty,baseline,refused);check(!no.published&&!oldPoints.empty(),"deterministic proposal witness");
  auto randomized=options;randomized.maximumCandidatePoints=128;randomized.maximumRandomCandidatePoints=64;randomized.randomSeed=17;randomized.interiorAdmission=[&](Vec3 p){return std::none_of(oldPoints.begin(),oldPoints.end(),[&](Vec3 q){return p.x==q.x&&p.y==q.y&&p.z==q.z;});};mesh::Mesh a,b;const auto sourceBefore=bytes(f.source);auto first=bnd::recoverTetPlcCavities(f.source,empty,randomized,a);auto second=bnd::recoverTetPlcCavities(f.source,empty,randomized,b);verify(f,a,first,1./6);verify(f,b,second,1./6);check(first.randomCandidatePointsProposed==64&&first.candidatePointsTested>0&&bytes(a)==bytes(b)&&bytes(f.source)==sourceBefore,"seeded retry reproducibility/source identity");
  randomized.randomSeed=29;mesh::Mesh c;auto other=bnd::recoverTetPlcCavities(f.source,empty,randomized,c);verify(f,c,other,1./6);check(bytes(a)!=bytes(c),"seed did not diversify interior proposals");++passed;}
 {auto f=tetra();mesh::Mesh empty,result=f.source;const auto before=bytes(result);auto bad=options;bad.maximumRandomCandidatePoints=bad.maximumCandidatePoints+1;auto out=bnd::recoverTetPlcCavities(f.source,empty,bad,result);check(!out.published&&out.reason==TetPlcCavityReason::InvalidRequest&&bytes(result)==before,"unbounded random candidate request");bad=options;bad.maximumRandomCandidatePoints=16;bad.interiorAdmission=[](Vec3){return false;};out=bnd::recoverTetPlcCavities(f.source,empty,bad,result);check(!out.published&&out.reason==TetPlcCavityReason::Unavailable&&bytes(result)==before&&out.randomCandidatePointsProposed==16,"random proposal bypassed material admission");++passed;}
 {auto f=ringBridge();const auto source=bytes(f.source),partial=bytes(f.partial);auto limited=options;limited.maximumPatches=1;limited.maximumCandidatePoints=128;mesh::Mesh result=f.source;auto blocked=bnd::recoverTetPlcCavities(f.source,f.partial,limited,result);check(!blocked.published&&blocked.reason==TetPlcCavityReason::BudgetExceeded&&bytes(result)==source,"separate-cavity patch bound");limited.maximumMergedComponents=1;auto merged=bnd::recoverTetPlcCavities(f.source,f.partial,limited,result);if(!merged.published)std::cerr<<"merge fixture: "<<merged.diagnostic<<" remaining "<<merged.remainingFrontFaces<<" other "<<merged.otherComponentBlockedCandidates<<'\n';verify(f,result,merged,25);check(merged.mergedComponents==1&&merged.cavitiesClosed==2&&merged.patchesCommitted==1&&bytes(f.source)==source&&bytes(f.partial)==partial,"bounded component merge/source transaction");++passed;}
 std::cout<<"{\"passed\":"<<passed<<",\"failed\":0}\n";return 0;
 }catch(const std::exception& e){failAllocation.store(false);std::cerr<<e.what()<<'\n';return 1;}}
