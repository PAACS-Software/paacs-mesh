// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Standalone checks for experimental recovery. These do not enable a CAD route.
#include "TetPlcGiftWrapping.hpp"
#include "TetBistellarRecovery.hpp"
#include <stdexcept>
#include <iomanip>
using namespace paacs;
using Front=bnd::detail::PlcGiftWrapping;
void require(bool condition,const char* reason){if(!condition)throw std::runtime_error(reason);}
std::vector<Front::Face> exterior(const std::vector<bnd::TetFacetVertex>& vertices,const std::vector<bnd::TetFacetCell>& cells){
 bnd::TetFacetComplexView view{vertices,cells,{},{},{},{3,0}};
 auto uses=bnd::tet_facet_detail::candidateFaceMap(view,cells);std::vector<Front::Face> faces;
 for(const auto& [key,rows]:uses){(void)key;if(rows.size()==1)faces.push_back(rows[0].directed);else require(rows.size()==2&&rows[0].sense==-rows[1].sense,"fixture interface");}return faces;
}
double volume(const Front& g){double total=0;for(auto cell:g.cells){const auto& v=cell.vertices;auto a=g.vertices[v[0]].position,b=g.vertices[v[1]].position,c=g.vertices[v[2]].position,d=g.vertices[v[3]].position;total+=(b-a).dot((c-a).cross(d-a))/6.;}return total;}
void audit(const Front& g){require(g.front.empty(),"unfinished front");bnd::TetFacetComplexView view{g.vertices,g.cells,{},{},{},{3,0}};
 std::map<bnd::TetFacetFaceKey,int> expected;for(auto f:g.boundary)expected.emplace(bnd::canonicalTetFacetFace(f),bnd::tet_facet_detail::faceSense(view,f));
 require(bnd::tet_facet_detail::boundarySenses(bnd::tet_facet_detail::candidateFaceMap(view,g.cells))==expected,"source facet mismatch");
 for(std::size_t i=0;i<g.cells.size();++i){auto c=g.cells[i];auto ids=c.vertices;require(predicates::orient3d(g.vertices[ids[0]].position,g.vertices[ids[1]].position,g.vertices[ids[2]].position,g.vertices[ids[3]].position)==predicates::Sign::Positive,"nonpositive cell");for(std::size_t j=0;j<i;++j)require(!bnd::tet_facet_detail::tetImproperlyIntersects(view,c,g.cells[j]),"overlapping cells");}
}
Front tetra(){Front g;g.vertices={{{0,0,0},1},{{1,0,0},2},{{0,1,0},3},{{0,0,1},4}};g.boundary=exterior(g.vertices,{{{0,1,2,3},1,{3,0}}});return g;}
int main(){try{unsigned passed=0;
 {auto g=tetra();const auto source=g.boundary;require(g.run(100,5),"tet recovery");audit(g);require(std::abs(volume(g)-1./6)<1e-14,"tet volume");require(g.boundary==source,"tet source changed");++passed;}
 {Front g;std::vector<Vec3> p{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};for(auto q:p)g.vertices.push_back({q,g.vertices.size()+1});std::vector<bnd::TetFacetCell> cells;for(auto c:std::array<std::array<unsigned,4>,6>{{{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}}})cells.push_back({c,cells.size()+1,{3,0}});g.boundary=exterior(g.vertices,cells);require(g.run(100,5),"cube recovery");audit(g);require(std::abs(volume(g)-1)<1e-14,"cube volume");++passed;}
 {auto g=tetra();auto other=g.vertices;for(auto v:other){v.position.x+=4;v.sequence=g.vertices.size()+1;g.vertices.push_back(v);}g.boundary=exterior(g.vertices,{{{0,1,2,3},1,{3,0}},{{4,5,6,7},2,{3,0}}});require(g.run(100,5),"disconnected recovery");audit(g);require(std::abs(volume(g)-1./3)<1e-14,"disconnected volume");for(auto c:g.cells){std::array<Vec3,4> p;for(unsigned j=0;j<4;++j)p[j]=g.vertices[c.vertices[j]].position;require(!bnd::tet_facet_detail::pointStrictlyInsideTet({2,.1,.1},p),"filled required gap");}++passed;}
 {auto g=tetra();g.boundary.pop_back();require(!g.run(100,5)&&g.cells.empty()&&g.front.empty(),"open boundary refusal is not transactional");++passed;}
 {auto g=tetra();g.boundary[0][0]=100;require(!g.run(100,5)&&g.cells.empty(),"invalid node reference");++passed;}
 {auto g=tetra();std::stop_source stop;stop.request_stop();require(!g.run(100,5,stop.get_token())&&g.cells.empty(),"cancelled run");++passed;}
 {auto g=tetra();require(!g.run(0,5)&&g.cells.empty(),"cell budget");require(!g.run(100,0)&&g.cells.empty(),"time budget");++passed;}
 {auto g=tetra();auto other=g.vertices;for(auto v:other){v.position.x+=4;v.sequence=g.vertices.size()+1;g.vertices.push_back(v);}g.boundary=exterior(g.vertices,{{{0,1,2,3},1,{3,0}},{{4,5,6,7},2,{3,0}}});require(!g.run(1,5)&&g.cells.size()<=1,"kernel exceeded cell budget");++passed;}
 {auto g=tetra();const auto source=g.boundary;require(g.initialize(),"seed fixture");require(g.seedUnfilledFront(1)==1&&g.vertices.size()==5,"front seed cap");require(g.run(100,5),"seeded tetrahedron recovery");audit(g);require(g.boundary==source&&std::abs(volume(g)-1./6)<1e-14,"seeding changed source or volume");++passed;}
 {auto g=tetra();require(g.initialize(),"cancelled seed fixture");const auto size=g.vertices.size();std::stop_source stop;stop.request_stop();g.cancellation=stop.get_token();require(g.seedUnfilledFront(10)==0&&g.vertices.size()==size,"cancelled seeding changed vertices");++passed;}
 {auto g=tetra();auto other=g.vertices;for(auto v:other){v.position.x+=4;v.sequence=g.vertices.size()+1;g.vertices.push_back(v);}g.boundary=exterior(g.vertices,{{{0,1,2,3},1,{3,0}},{{4,5,6,7},2,{3,0}}});require(g.initialize(),"disconnected seed fixture");require(g.seedUnfilledFront(8)>0,"disconnected seeds");require(g.run(100,5),"seeded disconnected recovery");audit(g);require(std::abs(volume(g)-1./3)<1e-14,"seeding bridged physical gap");++passed;}
 for(bool protectShared:std::array{false,true}){bnd::detail::BistellarComplex g;g.vertices={{{-1,-1,0},1},{{1,-1,0},2},{{0,1,0},3},{{0,0,-1},4},{{0,0,1},5}};std::vector<bnd::TetFacetCell> initial{{{0,1,2,4},1,{3,0}},{{0,2,1,3},2,{3,0}}};for(auto c:initial)g.insert(c);const auto surface=exterior(g.vertices,initial);for(auto f:surface)g.protectedFaces.insert(g.key(f));if(protectShared)g.protectedFaces.insert({0,1,2});const bool result=g.recoverEdge({3,4},16);require(result!=protectShared,"protected flip recovery");Front checked;checked.vertices=g.vertices;checked.boundary=surface;require(checked.initialize(),"flip fixture initialization");for(std::size_t i=0;i<g.cells.size();++i)if(g.alive[i])checked.insert(g.cells[i]);audit(checked);if(protectShared)require(checked.cells.size()==2,"protected refusal changed mesh");else require(g.edges.contains({3,4})&&checked.cells.size()==3,"missing recovered edge");++passed;}
 {bnd::detail::BistellarComplex g;g.vertices={{{-1,-1,0},1},{{1,-1,0},2},{{0,1,0},3},{{0,0,-1},4},{{0,0,1},5}};g.insert({{0,1,2,4},1,{3,1}});g.insert({{0,2,1,3},2,{3,2}});require(!g.recoverEdge({3,4},16)&&g.cells.size()==2&&std::count(g.alive.begin(),g.alive.end(),1)==2,"flip joined separately classified bodies");++passed;}
 {auto g=tetra();g.requireManifoldFront=true;require(g.run(100,5),"manifold tetrahedron recovery");audit(g);++passed;}
 {auto g=tetra();auto other=g.vertices;for(auto v:other){v.position.x+=4;v.sequence=g.vertices.size()+1;g.vertices.push_back(v);}g.boundary=exterior(g.vertices,{{{0,1,2,3},1,{3,0}},{{4,5,6,7},2,{3,0}}});g.requireManifoldFront=true;require(g.initialize(),"manifold gap fixture");const bnd::TetFacetCell bridge{{0,1,2,7},3,{3,0}};require(!g.manifoldAfter(bridge),"admitted a pinched vertex link");require(g.manifoldAfter({{0,1,2,3},1,{3,0}}),"refused legitimate component closure");g.insert({{0,1,2,3},1,{3,0}});g.insert({{4,5,6,7},2,{3,0}});require(g.manifoldAfter(bridge),"topology refusal was incorrectly treated as permanent");++passed;}
 std::cout<<"{\"scope\":\"Standalone experimental recovery checks; no CAD corpus credit\",\"passed\":"<<passed<<",\"failed\":0}\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
