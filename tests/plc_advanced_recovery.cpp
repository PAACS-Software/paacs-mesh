// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Reproductions for pinched fronts, optional apices, and retained repair progress.

#include "../src/bnd/TetPlcCavityExpansion.hpp"
#include <iostream>
using namespace paacs;
int sectors(){
 using namespace bnd::detail;
 PlcGiftWrapping g;
 const std::array<Vec3,6> p{{{0,0,0},{0,0,1},{1,0,.1},{0,1,.2},{-1,0,.1},{0,-1,.2}}};
 for(auto v:p)g.vertices.push_back({v,g.vertices.size()+1});
 for(auto ids:std::array<std::array<unsigned,4>,2>{{{0,1,2,3},{0,1,4,5}}}){
  bnd::TetFacetCell c;c.vertices=ids;
  for(unsigned j=0;j<4;++j){const auto f=bnd::tet_facet_detail::outwardFace(c,j);g.front.emplace(PlcGiftWrapping::key(f),PlcGiftWrapping::Front{f});}
 }
 
 auto revised=plcFrontComponents(g);
 bool okay=revised.size()==2;
 for(const auto& component:revised){okay&=component.size()==4;std::map<std::array<unsigned,2>,std::pair<unsigned,int>> edges;
  for(auto k:component){const auto f=g.front.at(k).directed;for(unsigned j=0;j<3;++j){const auto a=f[j],b=f[(j+1)%3];auto& e=edges[{std::min(a,b),std::max(a,b)}];++e.first;e.second+=a<b?1:-1;}}
  for(const auto& [k,e]:edges){(void)k;okay&=e.first==2&&e.second==0;}
 }
 std::cout<<"{\"passed\":"<<(okay?"true":"false")<<",\"previous_components\":"<<1<<",\"recovered_chambers\":"<<revised.size()<<"}\n";
 return okay?0:1;
}

#include "../src/bnd/TetPlcCavityExpansion.hpp"
#include <iostream>
using namespace paacs;
int topology(){
 using namespace bnd::detail;
 std::map<PlcExpansionFace,PlcExpansionFace> region;
 for(auto ids:std::array<std::array<unsigned,4>,2>{{{0,1,2,3},{0,4,5,6}}}){bnd::TetFacetCell cell;cell.vertices=ids;for(unsigned side=0;side<4;++side){const auto f=bnd::tet_facet_detail::outwardFace(cell,side);region.emplace(PlcGiftWrapping::key(f),f);}}
 std::uint64_t work=0;PlcFrontExpansionLimits budget;budget.work=&work;budget.maximumWork=1000000;const auto pinched=plcPinchedRegionVertices(region,budget);bool okay=pinched==std::set<unsigned>{0};
 std::vector<bnd::TetFacetVertex> vertices;for(auto p:std::array<Vec3,10>{{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1},{.5,.5,.5},{.2,.3,.4}}})vertices.push_back({p,vertices.size()+1});
 std::map<PlcExpansionFace,std::pair<unsigned,PlcExpansionFace>> uses;
 for(auto ids:std::array<std::array<unsigned,4>,6>{{{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}}}){bnd::TetFacetCell cell;cell.vertices=ids;for(unsigned side=0;side<4;++side){const auto f=bnd::tet_facet_detail::outwardFace(cell,side);auto& row=uses[PlcGiftWrapping::key(f)];++row.first;row.second=f;}}
 std::vector<PlcExpansionFace> faces;for(const auto& [key,row]:uses){(void)key;if(row.first==1)faces.push_back(row.second);}
 bnd::TetPlcPatchSearchOptions options;options.maximumStates=4096;options.maximumVertices=64;options.maximumCandidateTets=32768;PlcPatchSearchMetrics metrics;
 PlcPatchSearch search(vertices,faces,options,512,metrics,[&](std::uint64_t n){budget.tick(n);},[](Vec3 p){return !(p.x==.2&&p.y==.3&&p.z==.4);});
 std::vector<bnd::TetFacetCell> cells;const bool complete=search.run(cells);okay&=complete;for(auto cell:cells)okay&=std::find(cell.vertices.begin(),cell.vertices.end(),9)==cell.vertices.end();
 PlcGiftWrapping g;g.vertices=vertices;g.boundary=faces;g.interiorAdmission=[](Vec3 p){return p.x>=0&&p.x<=1&&p.y>=0&&p.y<=1&&p.z>=0&&p.z<=1;};okay&=g.initialize();std::vector<PlcExpansionFace> component;for(const auto& [k,row]:g.front){(void)row;component.push_back(k);}
 budget.patchSearch=options;budget.patchSearch.maximumStates=1;PlcFrontExpansionMetrics kernelMetrics;const bool kernelComplete=closePlcFrontRegion(g,component,kernelMetrics,budget);okay&=kernelComplete&&g.front.empty()&&g.cells.size()==12&&kernelMetrics.regionStates==0;
 for(auto cell:g.cells){const auto a=g.vertices[cell.vertices[0]].position,b=g.vertices[cell.vertices[1]].position,c=g.vertices[cell.vertices[2]].position,d=g.vertices[cell.vertices[3]].position;okay&=predicates::orient3d(a,b,c,d)==predicates::Sign::Positive&&g.interiorAdmission((a+b+c+d)*.25);}
 std::cout<<"{\"passed\":"<<(okay?"true":"false")<<",\"detected_pinched_vertices\":"<<pinched.size()<<",\"optional_candidate_rejected\":true,\"region_completed\":"<<(complete?"true":"false")<<",\"tetrahedra\":"<<cells.size()<<",\"states\":"<<metrics.states<<",\"kernel_region_completed\":"<<(kernelComplete?"true":"false")<<",\"kernel_tetrahedra\":"<<g.cells.size()<<",\"work\":"<<work<<"}\n";return okay?0:1;
}

#include "../src/bnd/TetPlcRetainedRecovery.hpp"
#include "paacs/bnd/TetFacetGeometry.hpp"
#include <iostream>
using namespace paacs;
int continuation(){
 mesh::Mesh source;
 const std::array<Vec3,8> corners{{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}}};
 const std::array<std::array<unsigned,4>,6> tets{{{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}}};
 int faceTag=1;
 for(unsigned cube=0;cube<2;++cube){const unsigned offset=unsigned(source.nodeCount());for(auto p:corners)source.addNode(p+Vec3{3.*cube,0,0},{0,int(source.nodeCount())});
  using Face=std::array<unsigned,3>;std::map<Face,std::pair<unsigned,Face>> faces;
  for(auto ids:tets){bnd::TetFacetCell c;for(unsigned j=0;j<4;++j)c.vertices[j]=ids[j]+offset;for(unsigned side=0;side<4;++side){const auto f=bnd::tet_facet_detail::outwardFace(c,side);auto k=f;std::sort(k.begin(),k.end());auto& row=faces[k];++row.first;row.second=f;}}
  for(const auto& [key,row]:faces)if(row.first==1){(void)key;source.addCell(mesh::CellType::Tri3,std::array<Index,3>{row.second[0],row.second[1],row.second[2]},{2,faceTag++});}
 }
 mesh::Mesh partial,destination,progress;destination.addNode({7,8,9},{0,0});destination.addNode({8,8,9},{0,1});destination.addCell(mesh::CellType::Line2,std::array<Index,2>{0,1},{1,1});
 bnd::TetPlcCavityOptions o;o.region={3,1};o.maximumNodes=1000;o.maximumCells=1000;o.maximumSourceTriangles=1000;o.maximumWork=5000000;o.maximumWorkingBytes=128ULL*1024*1024;o.interiorAdmission=[](Vec3 p){return p.x<1.1;};
 const auto first=bnd::recoverTetPlcCavitiesRetaining(source,partial,o,destination,&progress);
 bool okay=!first.published&&destination.nodeCount()==2&&destination.cellCount()==1&&destination.node(0).x==7&&progress.cellCount()==12&&progress.nodeCount()>source.nodeCount();
 for(Index i=0;i<source.nodeCount();++i){const auto a=source.node(i),b=progress.node(i);okay&=a.x==b.x&&a.y==b.y&&a.z==b.z&&source.nodeClass(i)==progress.nodeClass(i);}
 for(Index i=0;i<progress.cellCount();++i){const auto* p=progress.cellNodes(i);okay&=predicates::orient3d(progress.node(p[0]),progress.node(p[1]),progress.node(p[2]),progress.node(p[3]))==predicates::Sign::Positive;const auto center=(progress.node(p[0])+progress.node(p[1])+progress.node(p[2])+progress.node(p[3]))*.25;okay&=center.x<1.1;}
 o.interiorAdmission={};mesh::Mesh completed;const auto second=bnd::recoverTetPlcCavitiesRetaining(source,progress,o,completed,nullptr);
 okay&=second.published&&second.outputCells==24;
 if(second.published)for(Index i=0;i<progress.cellCount();++i){const auto* a=progress.cellNodes(i);const auto* b=completed.cellNodes(source.cellCount()+i);for(unsigned j=0;j<4;++j)okay&=a[j]==b[j];}
 std::cout<<"{\"passed\":"<<(okay?"true":"false")<<",\"partial_published\":"<<(first.published?"true":"false")<<",\"retained_tetrahedra\":"<<progress.cellCount()<<",\"remaining_front\":"<<first.remainingFrontFaces<<",\"resume_published\":"<<(second.published?"true":"false")<<",\"final_tetrahedra\":"<<second.outputCells<<",\"diagnostic\":\""<<first.diagnostic<<"\"}\n";
 return okay?0:1;
}

int main(){try{return sectors()|topology()|continuation();}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
