// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Dense-source differential regression: scheduling must preserve construction,
// including after a bounded run resumes, while avoiding repeated full scans.
#include "../src/bnd/TetPlcGiftWrapping.hpp"
#include "../src/geometry/OcctSeedSampling.hpp"
#include <iostream>
using namespace paacs;
using namespace paacs::bnd;
using namespace paacs::bnd::detail;
int main(){try{
 PlcGiftWrapping prototype;std::map<std::array<unsigned,3>,unsigned> ids;
 constexpr unsigned n=27;
 const auto vertex=[&](std::array<unsigned,3> p){auto [at,added]=ids.try_emplace(p,unsigned(prototype.vertices.size()));if(added)prototype.vertices.push_back({{double(p[0])/n,double(p[1])/n,double(p[2])/n},prototype.vertices.size()+1});return at->second;};
 for(unsigned axis=0;axis<3;++axis)for(unsigned side=0;side<2;++side)for(unsigned i=0;i<n;++i)for(unsigned j=0;j<n;++j){
  const auto point=[&](unsigned u,unsigned v){std::array<unsigned,3> p{};p[axis]=side*n;p[(axis+1)%3]=u;p[(axis+2)%3]=v;return vertex(p);};
  const auto a=point(i,j),b=point(i+1,j),c=point(i+1,j+1),d=point(i,j+1);
  for(auto f:std::array<PlcGiftWrapping::Face,2>{{{a,b,c},{a,c,d}}}){if(predicates::orient3d(prototype.vertices[f[0]].position,prototype.vertices[f[1]].position,prototype.vertices[f[2]].position,{.5,.5,.5})==predicates::Sign::Positive)std::swap(f[1],f[2]);prototype.boundary.push_back(f);}
 }
 prototype.vertices.push_back({{.5,.5,.5},prototype.vertices.size()+1});prototype.useLocalCandidates=true;
 auto legacy=prototype,incremental=prototype;legacy.useIncrementalFrontQueue=false;
 bool okay=true;double seconds[2]{};
 for(unsigned pass=0;pass<2;++pass){auto& g=pass?incremental:legacy;const auto start=std::chrono::steady_clock::now();
  g.run(32,120,{},2000000);g.run(64,120,{},2000000);seconds[pass]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 }
 okay&=legacy.cells.size()==64&&incremental.cells.size()==64&&legacy.front.size()==incremental.front.size();
 for(std::size_t i=0;i<legacy.cells.size()&&i<incremental.cells.size();++i){const auto a=legacy.cells[i],b=incremental.cells[i];okay&=a.vertices==b.vertices&&a.sequence==b.sequence;
  const auto& v=incremental.vertices;okay&=predicates::orient3d(v[b.vertices[0]].position,v[b.vertices[1]].position,v[b.vertices[2]].position,v[b.vertices[3]].position)==predicates::Sign::Positive;
 }
 for(const auto& [key,a]:legacy.front){auto at=incremental.front.find(key);if(at==incremental.front.end()){okay=false;break;}const auto& b=at->second;okay&=a.directed==b.directed&&a.next==b.next&&a.prepared==b.prepared&&a.full==b.full&&a.apex==b.apex;}
 okay&=legacy.metrics.apexTests==incremental.metrics.apexTests&&legacy.metrics.boundaryTests==incremental.metrics.boundaryTests&&legacy.metrics.pairTests==incremental.metrics.pairTests&&legacy.metrics.blocked==incremental.metrics.blocked;
 TetFacetComplexView view{incremental.vertices,incremental.cells,{},{},{},{3,0}};for(std::size_t i=0;i<incremental.cells.size();++i)for(std::size_t j=0;j<i;++j)okay&=!tet_facet_detail::tetImproperlyIntersects(view,incremental.cells[i],incremental.cells[j]);
 std::map<PlcGiftWrapping::Face,int> ledger;
 const auto record=[&](PlcGiftWrapping::Face f,int delta){int sense=1;for(unsigned i=0;i<3;++i)for(unsigned j=i+1;j<3;++j)if(f[j]<f[i]){std::swap(f[i],f[j]);sense=-sense;}ledger[f]+=sense*delta;};
 for(auto f:prototype.boundary)record(f,1);
 for(auto cell:incremental.cells)for(unsigned side=0;side<4;++side)record(tet_facet_detail::outwardFace(cell,side),-1);
 for(const auto& [key,row]:incremental.front){(void)key;record(row.directed,-1);}
 for(const auto& [key,count]:ledger){(void)key;okay&=count==0;}
 okay&=incremental.metrics.queueRowsScanned<legacy.metrics.queueRowsScanned/8;
 std::uint64_t samplingWork=0;const auto charge=[&](std::uint64_t n){samplingWork+=n;};
 std::vector<Vec3> dense;for(unsigned i=0;i<25;++i)for(unsigned j=0;j<25;++j)for(unsigned k=0;k<25;++k)dense.push_back({double(i)/24,double(j)/24,double(k)/24});
 const auto sampled=geo::detail::occtSeedBoundarySample(dense,{0,0,0},{1,1,1},charge);
 okay&=sampled.size()==2048&&std::set<Index>(sampled.begin(),sampled.end()).size()==2048&&sampled==geo::detail::occtSeedBoundarySample(dense,{0,0,0},{1,1,1},charge);
 std::array<unsigned,8> octants{};for(auto i:sampled){okay&=i>=0&&std::size_t(i)<dense.size();const auto p=dense[i];++octants[(p.x>=.5?1:0)|(p.y>=.5?2:0)|(p.z>=.5?4:0)];}for(auto count:octants)okay&=count>100;
 const std::array<Vec3,4> small{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}};const auto unchanged=geo::detail::occtSeedBoundarySample(small,{0,0,0},{1,1,1},charge);okay&=unchanged==std::vector<Index>{0,1,2,3}&&samplingWork>0;
 std::cout<<"{\"passed\":"<<(okay?"true":"false")<<",\"source_faces\":"<<prototype.boundary.size()<<",\"identical_tetrahedra\":"<<incremental.cells.size()<<",\"apex_tests\":"<<incremental.metrics.apexTests<<",\"legacy_queue_rows\":"<<legacy.metrics.queueRowsScanned<<",\"incremental_queue_rows\":"<<incremental.metrics.queueRowsScanned<<",\"legacy_seconds\":"<<seconds[0]<<",\"incremental_seconds\":"<<seconds[1]<<"}\n";return okay?0:1;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
