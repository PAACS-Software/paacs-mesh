// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// A pyramid with a frozen, triangulated base needs one tetrahedron per base
// triangle: all other points are coplanar. 325 vertices therefore require
// exactly 578 tetrahedra, despite fitting the existing region vertex budget.
#include "../src/bnd/TetPlcPatchSearch.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include <iostream>
using namespace paacs;
int main(){try{
 using namespace bnd;using namespace bnd::detail;constexpr unsigned n=17,apex=(n+1)*(n+1);
 std::vector<TetFacetVertex> vertices;for(unsigned j=0;j<=n;++j)for(unsigned i=0;i<=n;++i)vertices.push_back({{double(i)/n,double(j)/n,0},vertices.size()+1});vertices.push_back({{.5,.5,1},vertices.size()+1});
 const auto id=[](unsigned i,unsigned j){return j*(n+1)+i;};
 std::vector<PlcGiftWrapping::Face> faces;
 for(unsigned j=0;j<n;++j)for(unsigned i=0;i<n;++i){faces.push_back({id(i,j),id(i+1,j+1),id(i+1,j)});faces.push_back({id(i,j),id(i,j+1),id(i+1,j+1)});}
 for(unsigned i=0;i<n;++i){faces.push_back({id(i,0),id(i+1,0),apex});faces.push_back({id(n,i),id(n,i+1),apex});faces.push_back({id(n-i,n),id(n-i-1,n),apex});faces.push_back({id(0,n-i),id(0,n-i-1),apex});}
 TetPlcPatchSearchOptions options;options.maximumVertices=512;options.maximumStates=4096;options.maximumCandidateTets=32768;
 std::uint64_t work=0;const auto charge=[&](std::uint64_t count){work+=count;if(work>100000000)throw std::runtime_error("deep refill work budget");};
 PlcPatchSearchMetrics limitedMetrics,expandedMetrics;std::vector<TetFacetCell> limited,completed;
 PlcPatchSearch small(vertices,faces,options,512,limitedMetrics,charge,{});const bool smallDone=small.run(limited);
 PlcPatchSearch large(vertices,faces,options,2048,expandedMetrics,charge,{});const bool done=large.run(completed);
 bool okay=!smallDone&&limited.empty()&&limitedMetrics.bounded&&limitedMetrics.states>512&&done&&completed.size()==2*n*n;
 mesh::Mesh result;for(const auto& v:vertices)result.addNode(v.position,{0,0});for(auto f:faces)result.addCell(mesh::CellType::Tri3,std::array<Index,3>{f[0],f[1],f[2]},{2,1});
 long double volume=0;for(auto c:completed){const auto a=vertices[c.vertices[0]].position,b=vertices[c.vertices[1]].position,c0=vertices[c.vertices[2]].position,d=vertices[c.vertices[3]].position;okay&=predicates::orient3d(a,b,c0,d)==predicates::Sign::Positive;volume+=(b-a).cross(c0-a).dot(d-a)/6.;result.addCell(mesh::CellType::Tet4,std::array<Index,4>{c.vertices[0],c.vertices[1],c.vertices[2],c.vertices[3]},{3,1});}
 okay&=std::abs(double(volume)-1./3.)<1e-12&&mesh::auditLinearVolumeBoundary(result).valid;
 TetFacetComplexView view{vertices,completed,{},{},{},{3,1}};for(std::size_t i=0;i<completed.size();++i)for(std::size_t j=0;j<i;++j)okay&=!tet_facet_detail::tetImproperlyIntersects(view,completed[i],completed[j]);
 std::cout<<"{\"passed\":"<<(okay?"true":"false")<<",\"vertices\":"<<vertices.size()<<",\"boundary_faces\":"<<faces.size()<<",\"tetrahedra\":"<<completed.size()<<",\"old_cell_bound_reproduced\":"<<(!smallDone&&limitedMetrics.bounded?"true":"false")<<",\"states\":"<<expandedMetrics.states<<",\"work\":"<<work<<"}\n";return okay?0:1;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
