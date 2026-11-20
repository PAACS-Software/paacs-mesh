// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "OcctNativeWireAuthority.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/bnd/TetFacetGeometry.hpp"
#include <map>
#include <set>
namespace paacs::geo::detail {
struct OcctBoundaryDiagonalClosureResult {unsigned flips=0;double maximumHeight=0;};
struct AllOcctBoundaryDiagonalPairs {
 bool operator()(Index,Index,mesh::Classification)const noexcept{return true;}
};
template<class Charge,class AllowPair=AllOcctBoundaryDiagonalPairs>
OcctBoundaryDiagonalClosureResult closeOcctBoundaryDiagonals(const OcctModel& model,
 const mesh::Mesh& boundary,const mesh::Mesh& partial,unsigned maximumFlips,mesh::Mesh& candidate,
 Charge charge,AllowPair allowPair={}){
 using Face=std::array<Index,3>;using Edge=std::array<Index,2>;
 const auto key=[](Face f){std::sort(f.begin(),f.end());return f;};
 const auto edge=[](Index a,Index b){return Edge{std::min(a,b),std::max(a,b)};};
 const auto cyclic=[](Face a,Face b){for(unsigned j=0;j<3;++j)if(a[0]==b[j]&&a[1]==b[(j+1)%3]&&a[2]==b[(j+2)%3])return true;return false;};
 struct Use{unsigned count=0;Face directed;};std::map<Face,Use> cellFaces;
 for(Index i=0;i<partial.cellCount();++i){charge(1);bnd::TetFacetCell c;const auto* p=partial.cellNodes(i);for(unsigned j=0;j<4;++j)c.vertices[j]=unsigned(p[j]);for(unsigned side=0;side<4;++side){const auto f=bnd::tet_facet_detail::outwardFace(c,side);Face row{f[0],f[1],f[2]};auto& use=cellFaces[key(row)];++use.count;use.directed=row;}}
 std::map<Edge,std::vector<Index>> sourceEdges;std::set<Face> sourceFaces;
 for(Index i=0;i<boundary.cellCount();++i){const auto* p=boundary.cellNodes(i);Face f{p[0],p[1],p[2]};sourceFaces.insert(key(f));for(unsigned j=0;j<3;++j){charge(1);sourceEdges[edge(p[j],p[(j+1)%3])].push_back(i);}}
 std::map<Index,Face> replacement;std::set<Edge> newDiagonals;OcctBoundaryDiagonalClosureResult result;
 for(const auto& [diagonal,uses]:sourceEdges){
  charge(1);if(result.flips==maximumFlips)break;if(uses.size()!=2||replacement.contains(uses[0])||replacement.contains(uses[1]))continue;
  const auto face=boundary.cell(uses[0]).onEntity;if(face!=boundary.cell(uses[1]).onEntity||face.dim!=2||!allowPair(uses[0],uses[1],face))continue;
  const auto* p=boundary.cellNodes(uses[0]);const auto* q=boundary.cellNodes(uses[1]);Face first{p[0],p[1],p[2]},second{q[0],q[1],q[2]};
  if(cellFaces.contains(key(first))||cellFaces.contains(key(second)))continue;
  Index a=-1,b=-1,c=-1,d=-1;for(unsigned j=0;j<3;++j)if(edge(p[j],p[(j+1)%3])==diagonal){a=p[j];b=p[(j+1)%3];c=p[(j+2)%3];break;}for(auto id:second)if(id!=a&&id!=b)d=id;
  if(a<0||d<0||c==d||sourceEdges.contains(edge(c,d))||newDiagonals.contains(edge(c,d)))continue;
  const std::set<int> owners{face.tag};if(actualOcctWire(model,boundary.nodeClass(a),boundary.nodeClass(b),owners,[&]{charge(1);})!=0)continue;
  const Face x{c,a,d},y{c,d,b};if(sourceFaces.contains(key(x))||sourceFaces.contains(key(y)))continue;
  const auto foundX=cellFaces.find(key(x)),foundY=cellFaces.find(key(y));
  if(foundX==cellFaces.end()||foundY==cellFaces.end()||foundX->second.count!=1||foundY->second.count!=1||!cyclic(x,foundX->second.directed)||!cyclic(y,foundY->second.directed))continue;
  const auto A=boundary.node(a),B=boundary.node(b),C=boundary.node(c),D=boundary.node(d);
  const double area=std::max({(B-A).cross(C-A).norm(),(A-B).cross(D-B).norm(),(A-C).cross(D-C).norm(),(D-C).cross(B-C).norm()});
  if(!(area>0&&std::isfinite(area)))continue;
  const double height=std::abs(predicates::orient3dDeterminantEstimate(A,B,C,D))/area;
  if(!(std::isfinite(height)&&height<=8*std::max(model.tolerance(2,face.tag),1e-10)))continue;
  replacement.emplace(uses[0],x);replacement.emplace(uses[1],y);newDiagonals.insert(edge(c,d));++result.flips;result.maximumHeight=std::max(result.maximumHeight,height);
 }
 if(result.flips){mesh::Mesh staged;for(Index i=0;i<boundary.nodeCount();++i){charge(1);staged.addNode(boundary.node(i),boundary.nodeClass(i));}for(Index i=0;i<boundary.cellCount();++i){charge(1);const auto at=replacement.find(i);if(at!=replacement.end())staged.addCell(mesh::CellType::Tri3,at->second,boundary.cell(i).onEntity);else staged.addCell(mesh::CellType::Tri3,std::span<const Index>(boundary.cellNodes(i),3),boundary.cell(i).onEntity);}candidate.swap(staged);}
 return result;
}
}
