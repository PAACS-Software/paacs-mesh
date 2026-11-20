// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private value-object visibility expansion. Source PLC triangles are immutable;
// all proposed cells are checked before replacing any part of the input complex.
#include "TetPlcGiftWrapping.hpp"
#include "TetPlcPatchSearch.hpp"
#include <cstdint>
#include <stop_token>
namespace paacs::bnd::detail {
using PlcExpansionFace=PlcGiftWrapping::Face;
using PlcCellOwners=std::map<PlcExpansionFace,std::vector<std::size_t>>;
struct PlcExpansionStop {bool budget;};
struct PlcFrontExpansionLimits {
 unsigned maximumRemovedCells=128,maximumCandidatePoints=32,maximumNewCells=512,maximumComponentFaces=512;
 unsigned maximumRandomCandidatePoints=0;std::uint64_t randomSeed=1;
 unsigned maximumMergedComponents=0;
 std::uint64_t maximumNodes=UINT64_MAX;
 TetPlcPatchSearchOptions patchSearch;
 std::uint64_t* work=nullptr;std::uint64_t maximumWork=5'000'000;std::stop_token cancellation;
 void tick(std::uint64_t n=1)const {if(cancellation.stop_requested())throw PlcExpansionStop{false};if(!work||*work>maximumWork||n>maximumWork-*work)throw PlcExpansionStop{true};*work+=n;}
};
inline PlcCellOwners plcCellOwners(const PlcGiftWrapping& g,const PlcFrontExpansionLimits& limits){
    PlcCellOwners owners;
    for(std::size_t i=0;i<g.cells.size();++i)for(unsigned side=0;side<4;++side){
        limits.tick();
        owners[PlcGiftWrapping::key(tet_facet_detail::outwardFace(g.cells[i],side))].push_back(i);
    }
    return owners;
}
inline bool plcPatchClearOfSurvivors(PlcGiftWrapping& g,const TetFacetComplexView& view,
    const TetFacetCell& cell,PlcGiftWrapping::Box box,const std::set<std::size_t>& removed,
    const PlcFrontExpansionLimits& limits){
    const auto bounds=[&](std::size_t i){const auto b=g.cellBoxes[i];return PlcGrowingAabbIndex::Box{b.lo,b.hi};};
    return g.occupiedIndex.all(g.cells.size(),{box.lo,box.hi},bounds,
        [&](std::uint64_t work){limits.tick(work);},[&](std::size_t i){
            return removed.contains(i)||!g.cellsConflict(view,cell,g.cells[i]);
        });
}
// MSVC's empty map allocates a sentinel. Defer it until guarded recovery work,
// so constructing metrics before the noexcept entry's try cannot terminate.
struct PlcFrontExpansionMetrics{unsigned points=0,randomProposed=0,sourceBlocked=0,otherHoleBlocked=0,expansionBudget=0,kernelFailures=0,pairFailures=0,closed=0,removed=0,added=0,merged=0;std::uint64_t regionAttempts=0,regionStates=0,regionBacktracks=0,regionCandidates=0,regionPairCacheHits=0,regionPatches=0,regionBoundHits=0,regionVertexBoundHits=0,regionStateBoundHits=0,largestRegionVertices=0,largestRegionFaces=0;std::optional<std::map<PlcExpansionFace,std::uint64_t>> blockedSources;std::optional<std::set<PlcExpansionFace>> lastRegionDependencies,lastRegionTouched;};
inline std::vector<std::vector<PlcExpansionFace>> plcFrontComponents(const PlcGiftWrapping& g){std::map<std::array<unsigned,2>,std::vector<PlcExpansionFace>> byEdge;for(const auto& [k,row]:g.front)for(unsigned i=0;i<3;++i){auto a=k[i],b=k[(i+1)%3];byEdge[{std::min(a,b),std::max(a,b)}].push_back(k);}std::set<PlcExpansionFace> seen;std::vector<std::vector<PlcExpansionFace>> out;
  std::map<PlcExpansionFace,std::vector<PlcExpansionFace>> adjacent;
  for(const auto& [edge,faces]:byEdge){
   const auto join=[&](PlcExpansionFace a,PlcExpansionFace b){adjacent[a].push_back(b);adjacent[b].push_back(a);};
   if(faces.size()==2){join(faces[0],faces[1]);continue;}
   bool valid=faces.size()>2&&faces.size()%2==0;
   struct Ray{double angle;PlcExpansionFace face;int sense;unsigned third;};std::vector<Ray> rays;
   const auto origin=g.vertices[edge[0]].position;auto axis=g.vertices[edge[1]].position-origin;const double length=axis.norm();valid&=length>0&&std::isfinite(length);
   if(valid){axis=axis/length;const Vec3 reference=std::abs(axis.x)<=std::abs(axis.y)&&std::abs(axis.x)<=std::abs(axis.z)?Vec3{1,0,0}:std::abs(axis.y)<=std::abs(axis.z)?Vec3{0,1,0}:Vec3{0,0,1};auto u=axis.cross(reference);u=u/u.norm();const auto v=axis.cross(u);
    for(auto k:faces){const auto f=g.front.at(k).directed;int sense=0;unsigned third=0;for(unsigned j=0;j<3;++j){if(f[j]==edge[0]&&f[(j+1)%3]==edge[1])sense=1;if(f[j]==edge[1]&&f[(j+1)%3]==edge[0])sense=-1;if(f[j]!=edge[0]&&f[j]!=edge[1])third=f[j];}const auto p=g.vertices[third].position-origin;const double x=p.dot(u),y=p.dot(v);valid&=sense!=0&&std::isfinite(x)&&std::isfinite(y)&&(x!=0||y!=0);rays.push_back({std::atan2(y,x),k,sense,third});}
    std::sort(rays.begin(),rays.end(),[](const auto& a,const auto& b){return a.angle!=b.angle?a.angle<b.angle:a.face<b.face;});
    for(std::size_t j=0;j<rays.size();++j){const auto& a=rays[j];const auto& b=rays[(j+1)%rays.size()];valid&=a.sense==-b.sense;
     if(a.angle==b.angle)valid=false;
     // Exact signs certify each floating-point angular ordering decision.
     const auto sign=predicates::orient3d(origin,g.vertices[edge[1]].position,g.vertices[a.third].position,g.vertices[b.third].position);
     double delta=b.angle-a.angle;if(delta<0)delta+=6.2831853071795864769;
     if(delta<3.14159265358979323846)valid&=sign==predicates::Sign::Positive;else if(delta>3.14159265358979323846)valid&=sign==predicates::Sign::Negative;
    }
   }
   if(valid){for(std::size_t j=0;j<rays.size();++j)if(rays[j].sense<0)join(rays[j].face,rays[(j+1)%rays.size()].face);}
   else {for(auto a:faces)for(auto b:faces)if(a!=b)adjacent[a].push_back(b);}
  }
  for(const auto& [first,row]:g.front){(void)row;if(seen.contains(first))continue;std::vector<PlcExpansionFace> queue{first};seen.insert(first);for(std::size_t at=0;at<queue.size();++at)for(auto next:adjacent[queue[at]])if(seen.insert(next).second)queue.push_back(next);out.push_back(std::move(queue));}
  return out;
}
inline bool plcInsideFront(const PlcGiftWrapping& g,const std::vector<PlcExpansionFace>& component,Vec3 p,const PlcFrontExpansionLimits& limits){double sum=0;for(auto k:component){limits.tick();auto f=g.front.at(k).directed;const auto u=g.vertices[f[0]].position-p,v=g.vertices[f[1]].position-p,w=g.vertices[f[2]].position-p;sum+=2*std::atan2(u.dot(v.cross(w)),u.norm()*v.norm()*w.norm()+u.dot(v)*w.norm()+v.dot(w)*u.norm()+w.dot(u)*v.norm());}return std::abs(sum)>6.283185307179586;}
inline std::set<unsigned> plcPinchedRegionVertices(const std::map<PlcExpansionFace,PlcExpansionFace>& region,const PlcFrontExpansionLimits& limits){
 std::map<unsigned,std::map<unsigned,std::set<unsigned>>> links;
 for(const auto& [k,f]:region){(void)k;for(unsigned j=0;j<3;++j){limits.tick();const auto a=f[j],b=f[(j+1)%3],c=f[(j+2)%3];links[a][b].insert(c);links[a][c].insert(b);}}
 std::set<unsigned> bad;
 for(const auto& [vertex,link]:links){bool valid=true;for(const auto& [other,next]:link){(void)other;limits.tick();valid&=next.size()==2;}
  std::set<unsigned> seen;std::vector<unsigned> queue{link.begin()->first};seen.insert(queue.front());for(std::size_t i=0;i<queue.size();++i)for(auto next:link.at(queue[i])){limits.tick();if(seen.insert(next).second)queue.push_back(next);}valid&=seen.size()==link.size();if(!valid)bad.insert(vertex);
 }
 return bad;
}
inline bool closePlcFrontComponent(PlcGiftWrapping& g,const std::vector<PlcExpansionFace>& component,PlcFrontExpansionMetrics& metrics,const PlcFrontExpansionLimits& limits,
 const PlcCellOwners* sharedOwners=nullptr){
 const auto maxRemoved=limits.maximumRemovedCells;
 if(component.size()>limits.maximumComponentFaces)throw PlcExpansionStop{true};
 const auto localOwners=sharedOwners?PlcCellOwners{}:plcCellOwners(g,limits);
 const auto& owners=sharedOwners?*sharedOwners:localOwners;
 const std::set<PlcExpansionFace> baseOriginal(component.begin(),component.end());std::set<PlcExpansionFace> source;for(auto f:g.boundary)source.insert(PlcGiftWrapping::key(f));std::vector<Vec3> points;
 std::vector<std::vector<PlcExpansionFace>> allComponents;std::map<PlcExpansionFace,unsigned> componentOf;
 if(limits.maximumMergedComponents){allComponents=plcFrontComponents(g);for(unsigned i=0;i<allComponents.size();++i)for(auto f:allComponents[i]){limits.tick();componentOf.emplace(f,i);}}
 auto addPoint=[&](Vec3 p,bool adjacentCell=false){limits.tick();if(points.size()>=4096)throw PlcExpansionStop{true};if(!bnd::tet_facet_detail::finite(p)||(!adjacentCell&&!plcInsideFront(g,component,p,limits)))return;for(auto q:points)if(q.x==p.x&&q.y==p.y&&q.z==p.z)return;points.push_back(p);};
 const auto origin=g.vertices[component.front()[0]].position;long double volume=0;Vec3 centroid{};
 for(auto k:component){limits.tick();auto f=g.front.at(k).directed;auto a=g.vertices[f[0]].position,b=g.vertices[f[1]].position,c=g.vertices[f[2]].position;const auto v=(a-origin).dot((b-origin).cross(c-origin))/6.;volume+=v;centroid=centroid+(origin+a+b+c)*(v*.25);}if(volume>0)addPoint(centroid/double(volume));
 std::set<unsigned> used;for(auto c:g.cells)used.insert(c.vertices.begin(),c.vertices.end());for(unsigned i=g.fixedVertexCount;i<g.vertices.size();++i)if(!used.contains(i))addPoint(g.vertices[i].position);
 auto bounds=g.box(g.front.at(component.front()).directed);for(auto k:component){limits.tick();auto q=g.box(g.front.at(k).directed);bounds.lo={std::min(bounds.lo.x,q.lo.x),std::min(bounds.lo.y,q.lo.y),std::min(bounds.lo.z,q.lo.z)};bounds.hi={std::max(bounds.hi.x,q.hi.x),std::max(bounds.hi.y,q.hi.y),std::max(bounds.hi.z,q.hi.z)};}const auto span=(bounds.hi-bounds.lo).norm();addPoint((bounds.lo+bounds.hi)*.5);
 for(auto k:component){limits.tick();auto f=g.front.at(k).directed;auto a=g.vertices[f[0]].position,b=g.vertices[f[1]].position,c=g.vertices[f[2]].position;auto n=(b-a).cross(c-a);const auto length=n.norm();if(!(length>0))continue;n=n*(-1./length);const auto center=(a+b+c)/3.;
  for(double fraction:std::array{.01,.001,.0001,.05})addPoint(center+n*(span*fraction));}
 // A material-valid recovery point can lie in an adjacent accepted cell,
 // especially when the entire initial sliver cavity is outside curved CAD.
 // Visibility expansion must remove the intervening live cells; unchanged
 // source barriers and the final exact contact/front tests still apply.
 std::set<std::size_t> adjacent;for(auto k:component)if(auto at=owners.find(k);at!=owners.end())for(auto id:at->second)adjacent.insert(id);
 for(auto id:adjacent){limits.tick();Vec3 p{};for(auto v:g.cells[id].vertices)p=p+g.vertices[v].position*.25;addPoint(p,true);}
 const auto deterministicCount=points.size();
 if(limits.maximumRandomCandidatePoints){
  auto state=PlcGiftWrapping::mix(limits.randomSeed);for(auto k:component)for(auto v:k)state=PlcGiftWrapping::mix(state^v);
  const auto uniform=[&](){state=PlcGiftWrapping::mix(state);return double(state>>11)*0x1.0p-53;};
  std::vector<std::size_t> adjacentIds(adjacent.begin(),adjacent.end());
  for(unsigned attempt=0;attempt<limits.maximumRandomCandidatePoints;++attempt){limits.tick();++metrics.randomProposed;
   if(!adjacentIds.empty()&&attempt%2==0){const auto id=adjacentIds[std::min(adjacentIds.size()-1,std::size_t(uniform()*adjacentIds.size()))];std::array<double,3> u{uniform(),uniform(),uniform()};std::sort(u.begin(),u.end());const std::array<double,4> weights{u[0],u[1]-u[0],u[2]-u[1],1-u[2]};Vec3 p{};for(unsigned j=0;j<4;++j)p=p+g.vertices[g.cells[id].vertices[j]].position*weights[j];addPoint(p,true);}
   else addPoint({bounds.lo.x+(bounds.hi.x-bounds.lo.x)*uniform(),bounds.lo.y+(bounds.hi.y-bounds.lo.y)*uniform(),bounds.lo.z+(bounds.hi.z-bounds.lo.z)*uniform()});
  }
 }
 auto score=[&](Vec3 p){unsigned count=0;for(auto k:component){limits.tick();auto f=g.front.at(k).directed;count+=predicates::orient3d(g.vertices[f[0]].position,g.vertices[f[1]].position,g.vertices[f[2]].position,p)!=predicates::Sign::Negative;}return count;};
 std::vector<std::pair<unsigned,std::size_t>> ranked,randomRanked;ranked.reserve(deterministicCount);for(std::size_t i=0;i<points.size();++i)(i<deterministicCount?ranked:randomRanked).emplace_back(score(points[i]),i);std::stable_sort(ranked.begin(),ranked.end());std::stable_sort(randomRanked.begin(),randomRanked.end());
 const auto randomCount=std::min<std::size_t>(randomRanked.size(),limits.maximumRandomCandidatePoints);if(ranked.size()>limits.maximumCandidatePoints-randomCount)ranked.resize(limits.maximumCandidatePoints-randomCount);randomRanked.resize(randomCount);std::vector<Vec3> ordered;for(auto [rank,index]:ranked){(void)rank;ordered.push_back(points[index]);}for(auto [rank,index]:randomRanked){(void)rank;ordered.push_back(points[index]);}points=std::move(ordered);
 for(auto point:points){limits.tick();++metrics.points;auto original=baseOriginal;std::map<PlcExpansionFace,PlcExpansionFace> front;for(auto k:component)front.emplace(k,g.front.at(k).directed);std::set<std::size_t> removed;std::set<unsigned> merged;bool failed=false;
  if(g.interiorAdmission&&!g.interiorAdmission(point))continue;
  while(true){limits.tick();std::optional<PlcExpansionFace> bad;for(const auto& [k,f]:front)if(predicates::orient3d(g.vertices[f[0]].position,g.vertices[f[1]].position,g.vertices[f[2]].position,point)!=predicates::Sign::Negative){bad=k;break;}if(!bad)break;
   if(source.contains(*bad)){++metrics.sourceBlocked;if(!metrics.blockedSources)metrics.blockedSources.emplace();++(*metrics.blockedSources)[*bad];failed=true;break;}auto at=owners.find(*bad);std::optional<std::size_t> next;unsigned remainingOwners=0;if(at!=owners.end())for(auto owner:at->second)if(!removed.contains(owner)){next=owner;++remainingOwners;}if(remainingOwners!=1){++metrics.kernelFailures;failed=true;break;}
   if(removed.size()==maxRemoved){++metrics.expansionBudget;failed=true;break;}const auto id=*next;std::set<unsigned> additions;
   for(unsigned side=0;side<4;++side){auto f=bnd::tet_facet_detail::outwardFace(g.cells[id],side);auto k=PlcGiftWrapping::key(f);if(g.front.contains(k)&&!original.contains(k)){if(!limits.maximumMergedComponents){++metrics.otherHoleBlocked;failed=true;break;}additions.insert(componentOf.at(k));}}
   if(failed)break;
   if(merged.size()+additions.size()>limits.maximumMergedComponents){++metrics.otherHoleBlocked;failed=true;break;}
   for(auto index:additions){limits.tick();if(original.size()+allComponents[index].size()>limits.maximumComponentFaces){++metrics.expansionBudget;failed=true;break;}merged.insert(index);
    for(auto k:allComponents[index]){limits.tick();original.insert(k);const auto f=g.front.at(k).directed;if(auto frontAt=front.find(k);frontAt!=front.end())front.erase(frontAt);else front.emplace(k,f);}}
   if(failed)break;
   removed.insert(id);
   for(unsigned side=0;side<4;++side){auto f=bnd::tet_facet_detail::outwardFace(g.cells[id],side);auto k=PlcGiftWrapping::key(f);if(auto found=front.find(k);found!=front.end())front.erase(found);else front.emplace(k,f);}
  }
  if(failed)continue;
  if(front.size()>limits.maximumNewCells)throw PlcExpansionStop{true};
  auto vertices=g.vertices;const auto apex=unsigned(vertices.size());bool duplicate=false;for(unsigned i=0;i<vertices.size();++i){const auto q=vertices[i].position;duplicate|=q.x==point.x&&q.y==point.y&&q.z==point.z;}if(duplicate){++metrics.kernelFailures;continue;}vertices.push_back({point,vertices.size()+1});
  std::vector<bnd::TetFacetCell> proposed,kept;std::vector<PlcGiftWrapping::Box> boxes;for(std::size_t i=0;i<g.cells.size();++i)if(!removed.contains(i))kept.push_back(g.cells[i]);for(const auto& [k,f]:front){(void)k;proposed.push_back({{f[0],f[2],f[1],apex},proposed.size()+1,{3,0}});}
  bnd::TetFacetComplexView view{vertices,kept,{},{},{},{3,0}};bool valid=true;
  for(auto c:proposed){limits.tick();const auto ids=c.vertices;if(predicates::orient3d(vertices[ids[0]].position,vertices[ids[1]].position,vertices[ids[2]].position,vertices[ids[3]].position)!=predicates::Sign::Positive){valid=false;break;}if(!g.admitCentroid((vertices[ids[0]].position+vertices[ids[1]].position+vertices[ids[2]].position+vertices[ids[3]].position)*.25)){valid=false;break;}}
  if(!valid){++metrics.kernelFailures;continue;}
  // Use a separate immutable source index; no original candidate is changed.
  PlcGiftWrapping validator;validator.vertices=vertices;validator.boundary=g.boundary;validator.cancellation=limits.cancellation;validator.interiorAdmission=g.interiorAdmission;validator.centroidAdmission=g.centroidAdmission;validator.useExactContactShortcut=g.useExactContactShortcut;if(!validator.initialize())return false;
  for(std::size_t i=0;i<proposed.size()&&valid;++i){const auto c=proposed[i];const auto box=validator.box(c.vertices);if(!validator.boundaryAccepts(c,box,0)){valid=false;break;}
   if(g.useExactContactShortcut)valid=plcPatchClearOfSurvivors(g,view,c,box,removed,limits);
   else for(std::size_t j=0;j<kept.size();++j)if((limits.tick(),true)&&PlcGiftWrapping::overlap(box,validator.box(kept[j].vertices))&&g.cellsConflict(view,c,kept[j])){valid=false;break;}
   for(std::size_t j=0;j<i&&valid;++j)if((limits.tick(),true)&&PlcGiftWrapping::overlap(box,validator.box(proposed[j].vertices))&&g.cellsConflict(view,c,proposed[j]))valid=false;
  }
  if(!valid){++metrics.pairFailures;continue;}
  for(auto c:kept)validator.insert(c);
  for(auto c:proposed)validator.insert(c);
  if(validator.front.size()!=g.front.size()-original.size()){++metrics.pairFailures;continue;}
  for(const auto& [k,row]:validator.front){auto old=g.front.find(k);if(original.contains(k)||old==g.front.end()||bnd::tet_facet_detail::faceSense(view,row.directed)!=bnd::tet_facet_detail::faceSense(view,old->second.directed)){valid=false;break;}}
  if(!valid){++metrics.pairFailures;continue;}
  ++metrics.closed;metrics.removed+=unsigned(removed.size());metrics.added+=unsigned(proposed.size());metrics.merged+=unsigned(merged.size());g=std::move(validator);return true;
 }
 return false;
}

// A completed local search is checked against the unchanged surrounding mesh
// before one value-object replacement. Preparation/copying happens once per
// completed proposal, never at each backtracking state.
inline bool publishPlcRegion(PlcGiftWrapping& g,const std::set<PlcExpansionFace>& original,
 const std::set<std::size_t>& removed,const std::vector<TetFacetCell>& proposed,
 const PlcFrontExpansionLimits& limits){
 limits.tick(g.vertices.size()+g.boundary.size()+g.cells.size());
 PlcGiftWrapping validator;validator.vertices=g.vertices;validator.boundary=g.boundary;
 validator.cancellation=limits.cancellation;validator.interiorAdmission=g.interiorAdmission;validator.centroidAdmission=g.centroidAdmission;
 validator.useExactContactShortcut=g.useExactContactShortcut;
 validator.chargeWork=[&](std::uint64_t n){limits.tick(n);};if(!validator.initialize())return false;
 TetFacetComplexView view{g.vertices,g.cells,{},{},{},{3,0}};
 std::vector<std::pair<double,std::size_t>> kept;
 if(!g.useExactContactShortcut){
  kept.reserve(g.cells.size()-removed.size());
  for(std::size_t i=0;i<g.cells.size();++i)if(!removed.contains(i))kept.emplace_back(g.cellBoxes[i].lo.x,i);
 }
 std::sort(kept.begin(),kept.end());std::vector<PlcGiftWrapping::Box> boxes;
 for(std::size_t i=0;i<proposed.size();++i){limits.tick();const auto c=proposed[i];const auto box=validator.box(c.vertices);
  if(!validator.boundaryAccepts(c,box,0))return false;
  if(g.useExactContactShortcut&&!plcPatchClearOfSurvivors(g,view,c,box,removed,limits))return false;
  for(auto [lower,j]:kept){limits.tick();if(lower>box.hi.x)break;if(PlcGiftWrapping::overlap(box,g.cellBoxes[j])&&g.cellsConflict(view,c,g.cells[j]))return false;}
  for(std::size_t j=0;j<i;++j){limits.tick();if(PlcGiftWrapping::overlap(box,boxes[j])&&g.cellsConflict(view,c,proposed[j]))return false;}boxes.push_back(box);
 }
 for(std::size_t i=0;i<g.cells.size();++i)if(!removed.contains(i)){limits.tick();validator.insert(g.cells[i]);}
 for(auto c:proposed){limits.tick();validator.insert(c);}
 if(validator.front.size()!=g.front.size()-original.size())return false;
 for(const auto& [k,row]:validator.front){limits.tick();const auto at=g.front.find(k);if(original.contains(k)||at==g.front.end()||tet_facet_detail::faceSense(view,row.directed)!=tet_facet_detail::faceSense(view,at->second.directed))return false;}
 // Do not retain a callback capturing this stack frame in the value object.
 validator.chargeWork={};g=std::move(validator);return true;
}

inline bool closePlcFrontRegion(PlcGiftWrapping& g,const std::vector<PlcExpansionFace>& component,
 PlcFrontExpansionMetrics& metrics,const PlcFrontExpansionLimits& limits,const PlcCellOwners* sharedOwners=nullptr){
 if(!metrics.lastRegionDependencies)metrics.lastRegionDependencies.emplace();else metrics.lastRegionDependencies->clear();
 if(!metrics.lastRegionTouched)metrics.lastRegionTouched.emplace();else metrics.lastRegionTouched->clear();
 metrics.lastRegionDependencies->insert(component.begin(),component.end());
 if(!limits.patchSearch.maximumStates)return false;
 if(component.size()>limits.maximumComponentFaces){++metrics.regionBoundHits;return false;}
 const auto localOwners=sharedOwners?PlcCellOwners{}:plcCellOwners(g,limits);
 const auto& owners=sharedOwners?*sharedOwners:localOwners;
 std::set<PlcExpansionFace> original(component.begin(),component.end()),source;
 std::map<PlcExpansionFace,PlcExpansionFace> region;for(auto k:component)region.emplace(k,g.front.at(k).directed);
 for(auto f:g.boundary){limits.tick();source.insert(PlcGiftWrapping::key(f));}
 const auto components=plcFrontComponents(g);std::map<PlcExpansionFace,unsigned> componentOf;
 for(unsigned i=0;i<components.size();++i)for(auto f:components[i]){limits.tick();componentOf.emplace(f,i);}
 std::set<std::size_t> removed;std::set<unsigned> merged;
 for(unsigned shell=0;shell<=limits.patchSearch.maximumShells;++shell){
  const auto previousVertices=g.vertices.size();
  std::set<unsigned> ids;for(const auto& [k,f]:region){(void)f;limits.tick();ids.insert(k.begin(),k.end());}
  for(const auto& [k,f]:region){(void)f;metrics.lastRegionDependencies->insert(k);}for(auto id:removed)for(unsigned side=0;side<4;++side)metrics.lastRegionDependencies->insert(PlcGiftWrapping::key(tet_facet_detail::outwardFace(g.cells[id],side)));
  if(shell&&!removed.empty()){
   std::vector<std::pair<double,Vec3>> proposed;
   for(auto id:removed){limits.tick();const auto c=g.cells[id].vertices;const auto a=g.vertices[c[0]].position,b=g.vertices[c[1]].position,d=g.vertices[c[2]].position,e=g.vertices[c[3]].position;proposed.emplace_back(std::abs((b-a).cross(d-a).dot(e-a)),(a+b+d+e)*.25);}
   std::stable_sort(proposed.begin(),proposed.end(),[](const auto& a,const auto& b){return a.first>b.first;});
   unsigned added=0;
   for(auto [volume,p]:proposed){(void)volume;limits.tick();if(added==8||ids.size()>=limits.patchSearch.maximumVertices||g.vertices.size()>=limits.maximumNodes)break;if(g.interiorAdmission&&!g.interiorAdmission(p))continue;bool duplicate=false;for(auto id:ids){limits.tick();const auto q=g.vertices[id].position;if(q.x==p.x&&q.y==p.y&&q.z==p.z){duplicate=true;break;}}if(duplicate)continue;
    double winding=0;for(const auto& [k,f]:region){(void)k;limits.tick();const auto a=g.vertices[f[0]].position-p,b=g.vertices[f[1]].position-p,c=g.vertices[f[2]].position-p;winding+=2*std::atan2(a.dot(b.cross(c)),a.norm()*b.norm()*c.norm()+a.dot(b)*c.norm()+a.dot(c)*b.norm()+b.dot(c)*a.norm());}
    if(!(std::abs(winding)>6.283185307179586))continue;
    const auto index=unsigned(g.vertices.size());
    g.vertices.push_back({p,g.vertices.size()+1});
    ids.insert(index);
    ++added;
   }
   
  }
  metrics.largestRegionVertices=std::max(metrics.largestRegionVertices,std::uint64_t(ids.size()));metrics.largestRegionFaces=std::max(metrics.largestRegionFaces,std::uint64_t(region.size()));
  if(ids.size()>limits.patchSearch.maximumVertices||region.size()>limits.maximumComponentFaces){g.vertices.resize(previousVertices);++metrics.regionBoundHits;metrics.regionVertexBoundHits+=ids.size()>limits.patchSearch.maximumVertices;break;}
  std::vector<unsigned> global;
  {std::set<unsigned> fixed;for(const auto& [k,f]:region){(void)f;fixed.insert(k.begin(),k.end());}global.assign(fixed.begin(),fixed.end());for(auto id:ids)if(!fixed.contains(id))global.push_back(id);}
  std::vector<TetFacetVertex> vertices;std::map<unsigned,unsigned> local;
  for(auto id:global){limits.tick();local.emplace(id,unsigned(vertices.size()));vertices.push_back({g.vertices[id].position,vertices.size()+1});}
  std::vector<PlcExpansionFace> faces;for(const auto& [k,f]:region){(void)k;faces.push_back({local.at(f[0]),local.at(f[1]),local.at(f[2])});}
  PlcPatchSearchMetrics counters;++metrics.regionAttempts;std::vector<TetFacetCell> candidate;
  const auto pinchedVertices=plcPinchedRegionVertices(region,limits);
  bool complete=false;
  if(pinchedVertices.empty()&&g.vertices.size()<limits.maximumNodes){
   PlcGiftWrapping kernel;kernel.vertices=vertices;kernel.boundary=faces;kernel.interiorAdmission=g.interiorAdmission;kernel.centroidAdmission=g.centroidAdmission;kernel.cancellation=limits.cancellation;kernel.chargeWork=[&](std::uint64_t n){limits.tick(n);};
   if(kernel.initialize()&&kernel.finishKernels(limits.maximumNewCells)&&kernel.front.empty()){
    candidate=kernel.cells;for(std::size_t i=vertices.size();i<kernel.vertices.size();++i){global.push_back(unsigned(g.vertices.size()));g.vertices.push_back({kernel.vertices[i].position,g.vertices.size()+1});}complete=true;
    
   }
  }
  if(!complete&&pinchedVertices.empty()){PlcPatchSearch search(std::move(vertices),std::move(faces),limits.patchSearch,limits.maximumNewCells,counters,[&](std::uint64_t n){limits.tick(n);},g.interiorAdmission,g.centroidAdmission);complete=search.run(candidate);}
  metrics.regionStates+=counters.states;metrics.regionBacktracks+=counters.backtracks;
  
  metrics.regionCandidates+=counters.candidates;metrics.regionPairCacheHits+=counters.pairCacheHits;metrics.regionBoundHits+=counters.bounded;metrics.regionStateBoundHits+=counters.stateBound;
  if(complete){for(auto& c:candidate)for(auto& id:c.vertices)id=global[id];
   std::set<PlcExpansionFace> touched;for(auto id:removed)for(unsigned side=0;side<4;++side)touched.insert(PlcGiftWrapping::key(tet_facet_detail::outwardFace(g.cells[id],side)));for(auto c:candidate)for(unsigned side=0;side<4;++side)touched.insert(PlcGiftWrapping::key(tet_facet_detail::outwardFace(c,side)));
   if(publishPlcRegion(g,original,removed,candidate,limits)){*metrics.lastRegionTouched=std::move(touched);++metrics.closed;++metrics.regionPatches;metrics.removed+=unsigned(removed.size());metrics.added+=unsigned(candidate.size());metrics.merged+=unsigned(merged.size());return true;}
   ++metrics.pairFailures;
  }
  g.vertices.resize(previousVertices);
  if(shell==limits.patchSearch.maximumShells)break;
  // Grow a connected cell layer across artificial cavity faces. Native source
  // triangles are barriers, and other unfilled components require admission.
  std::set<std::size_t> layer;
  using Edge=std::array<unsigned,2>;std::map<Edge,std::pair<unsigned,int>> edgeUses;
  for(const auto& [k,f]:region){(void)k;for(unsigned j=0;j<3;++j){limits.tick();const auto a=f[j],b=f[(j+1)%3];auto& use=edgeUses[{std::min(a,b),std::max(a,b)}];++use.first;use.second+=a<b?1:-1;}}
  std::set<Edge> pinched;for(const auto& [edge,use]:edgeUses)if(use.first!=2||use.second)pinched.insert(edge);
  if(!pinchedVertices.empty()){
   for(std::size_t id=0;id<g.cells.size();++id)if(!removed.contains(id))for(auto vertex:g.cells[id].vertices){limits.tick();if(pinchedVertices.contains(vertex)){layer.insert(id);break;}}
   
  }else if(!pinched.empty()){
   // A greedy front can leave chambers meeting only along an edge. A single
   // adjacent-cell addition cannot cure that topology: reopen its entire live
   // edge star before trying another fill. This is a targeted group, not a
   // global enlargement of every face of the cavity.
   for(std::size_t id=0;id<g.cells.size();++id)if(!removed.contains(id)){
    const auto v=g.cells[id].vertices;for(auto e:tet_facet_detail::kTetEdges){limits.tick();if(pinched.contains({std::min(v[e[0]],v[e[1]]),std::max(v[e[0]],v[e[1]])})){layer.insert(id);break;}}
   }
  }else for(const auto& [k,f]:region){(void)f;if(source.contains(k))continue;auto at=owners.find(k);if(at==owners.end())continue;
   for(auto id:at->second)if(!removed.contains(id))layer.insert(id);}
  bool progress=false;
  for(auto id:layer){limits.tick();if(removed.size()>=limits.maximumRemovedCells){++metrics.regionBoundHits;break;}
   std::set<unsigned> absorb;for(unsigned side=0;side<4;++side){const auto k=PlcGiftWrapping::key(tet_facet_detail::outwardFace(g.cells[id],side));
    if(g.front.contains(k)&&!original.contains(k))absorb.insert(componentOf.at(k));}
   if(merged.size()+absorb.size()>limits.maximumMergedComponents)continue;
   std::size_t addedFaces=0;for(auto part:absorb)addedFaces+=components[part].size();
   if(original.size()+addedFaces>limits.maximumComponentFaces){++metrics.regionBoundHits;continue;}
   for(auto part:absorb){merged.insert(part);for(auto k:components[part]){limits.tick();original.insert(k);region.emplace(k,g.front.at(k).directed);}}
   removed.insert(id);progress=true;
   for(unsigned side=0;side<4;++side){auto f=tet_facet_detail::outwardFace(g.cells[id],side);const auto k=PlcGiftWrapping::key(f);if(auto at=region.find(k);at!=region.end())region.erase(at);else region.emplace(k,f);}
  }
  if(!progress)break;
 }
 return false;
}

} // namespace paacs::bnd::detail
