// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Plain-value constrained advancing front. Source facets are immutable;
// no cell may cross them or overlap a previously accepted cell. Publication
// requires the oriented front to close completely and independent auditing.
#include "paacs/bnd/TetFacetGeometry.hpp"
#include <chrono>
#include <functional>
#include <iostream>
#include <numeric>
#include <stop_token>
#include <tuple>
#include "PlcSpatialIndex.hpp"

namespace paacs::bnd::detail {
class PlcGiftWrapping {
public:
    using Face=std::array<TetFacetVertexId,3>;
    struct Box {Vec3 lo,hi;};
    struct Metrics {std::uint64_t apexTests=0,boundaryTests=0,pairTests=0,inserted=0,blocked=0,kernelAttempts=0,kernelCached=0,frontSeedAttempts=0,frontSeeds=0,nonmanifoldProposals=0,queueRebuilds=0,queueRowsScanned=0;};
    struct Front {
        Face directed;
        std::vector<std::pair<double,TetFacetVertexId>> apex{};
        std::size_t next=0;
        bool prepared=false;
        bool full=true;
        std::set<TetFacetVertexId> geometryRejected{};
        bool topologyBlocked=false;
    };
    std::vector<TetFacetVertex> vertices;
    std::vector<unsigned char> sourceVertex;
    std::vector<Face> boundary;
    std::vector<TetFacetCell> cells;
    std::vector<Box> cellBoxes,boundaryBoxes;
    std::map<Face,Front> front;
    bool requireManifoldFront=false;
    bool useLocalCandidates=false;
    bool useExactContactShortcut=false;
    bool useIncrementalFrontQueue=true;
    bool manifoldIndexValid=false;
    bool lastTopologyRefusal=false;
    bool retryTopologyRequired=false;
    std::map<TetFacetVertexId,std::set<Face>> frontVertexFaces;
    std::set<std::vector<Face>> exhaustedKernels;
    Metrics metrics;
    static Face key(Face f){std::sort(f.begin(),f.end());return f;}
    static bool overlap(Box a,Box b){return a.lo.x<=b.hi.x&&b.lo.x<=a.hi.x&&a.lo.y<=b.hi.y&&b.lo.y<=a.hi.y&&a.lo.z<=b.hi.z&&b.lo.z<=a.hi.z;}
    template<std::size_t N> Box box(std::array<TetFacetVertexId,N> ids)const{
        auto p=vertices[ids[0]].position;Box b{p,p};for(auto id:ids){p=vertices[id].position;
            b.lo={std::min(b.lo.x,p.x),std::min(b.lo.y,p.y),std::min(b.lo.z,p.z)};b.hi={std::max(b.hi.x,p.x),std::max(b.hi.y,p.y),std::max(b.hi.z,p.z)};}return b;
    }
    struct Branch {Box bounds{};unsigned first=0,last=0,left=0,right=0;};
    std::vector<Branch> tree;std::vector<unsigned> order;
    bool initialized=false;
    std::stop_token cancellation;
    // Optional external material admission, in addition to exact PLC tests.
    // It can reject a proposal; it cannot bypass any geometric check.
    std::function<bool(Vec3)> interiorAdmission,centroidAdmission;
    bool admitCentroid(Vec3 point)const{return centroidAdmission?centroidAdmission(point):!interiorAdmission||interiorAdmission(point);}
    std::function<void(std::uint64_t)> chargeWork;
    std::function<void(std::uint64_t)> admitApexStorage;
    std::uint64_t apexCapacityBytes=0;
    std::ostream* progressStream=nullptr;
    void charge(std::uint64_t count=1)const{if(chargeWork)chargeWork(count);}
    PlcGrowingAabbIndex fixedIndex,occupiedIndex;
    TetFacetVertexId fixedVertexCount=0;
    unsigned build(unsigned first,unsigned last){
        const auto id=unsigned(tree.size());tree.push_back({});Box b=boundaryBoxes[order[first]];
        for(unsigned i=first+1;i<last;++i){const auto q=boundaryBoxes[order[i]];b.lo={std::min(b.lo.x,q.lo.x),std::min(b.lo.y,q.lo.y),std::min(b.lo.z,q.lo.z)};b.hi={std::max(b.hi.x,q.hi.x),std::max(b.hi.y,q.hi.y),std::max(b.hi.z,q.hi.z)};}
        tree[id].bounds=b;tree[id].first=first;tree[id].last=last;if(last-first<=8)return id;
        const auto span=b.hi-b.lo;const unsigned axis=span.x>=span.y&&span.x>=span.z?0:(span.y>=span.z?1:2);
        const auto coordinate=[&](unsigned i){auto q=boundaryBoxes[i];return axis==0?q.lo.x+q.hi.x:axis==1?q.lo.y+q.hi.y:q.lo.z+q.hi.z;};const auto mid=(first+last)/2;
        std::nth_element(order.begin()+first,order.begin()+mid,order.begin()+last,[&](auto a,auto c){return std::pair(coordinate(a),a)<std::pair(coordinate(c),c);});
        const auto left=build(first,mid),right=build(mid,last);tree[id].left=left;tree[id].right=right;return id;
    }
    bool boundaryAccepts(const TetFacetCell& cell,Box bounds,unsigned branch,bool canonicalVertices=false){
        charge();
        if(cancellation.stop_requested())return false;
        const auto& b=tree[branch];if(!overlap(bounds,b.bounds))return true;
        if(b.left||b.right)return boundaryAccepts(cell,bounds,b.left,canonicalVertices)&&
                                 boundaryAccepts(cell,bounds,b.right,canonicalVertices);
        TetFacetComplexView view{vertices,{}, {},{}, {},{3,0}};
        for(unsigned j=b.first;j<b.last;++j){const auto i=order[j];if(!overlap(bounds,boundaryBoxes[i]))continue;
            charge();
            ++metrics.boundaryTests;view.target.vertices=boundary[i];
            if(useExactContactShortcut&&tet_facet_detail::separatedByTetFace(view,cell,boundary[i],[]{}))continue;
            const auto contact=canonicalVertices?
                tet_facet_detail::classifyTetFacetContactWithCanonicalVertices(view,cell):
                classifyTetFacetContact(view,cell);
            if(contact!=TetFacetContact::Disjoint&&contact!=TetFacetContact::CoplanarExactFace&&contact!=TetFacetContact::DeclaredBoundaryTouch)return false;}
        return true;
    }
    bool cellsConflict(const TetFacetComplexView& view,const TetFacetCell& a,const TetFacetCell& b)const{
        if(useExactContactShortcut&&(tet_facet_detail::separatedByTetFace(view,a,b.vertices,[]{})||
            tet_facet_detail::separatedByTetFace(view,b,a.vertices,[]{})))return false;
        return tet_facet_detail::tetImproperlyIntersects(view,a,b);
    }
    void prepare(Front& row,bool forceFull=false){
        row.full=forceFull||requireManifoldFront||!useLocalCandidates;row.apex.clear();row.next=0;
        row.prepared=true;const auto f=row.directed;const auto a=vertices[f[0]].position,b=vertices[f[1]].position,c=vertices[f[2]].position;
        const auto u=b-a,v=c-a,n=u.cross(v);const auto norm=n.norm(),den=2*n.dot(n);if(!(norm>0&&std::isfinite(den)&&den>0))return;
        const auto center=a+(v.cross(n)*u.dot(u)+n.cross(u)*v.dot(v))/den;const auto radius=(center-a).dot(center-a);const auto inward=n*(-1./norm);
        std::vector<TetFacetVertexId> selected;
        if(row.full){selected.resize(vertices.size());std::iota(selected.begin(),selected.end(),0);}
        else{
            const auto middle=(a+b+c)/3.;const double reach=2.5*std::max({(a-b).norm(),(b-c).norm(),(c-a).norm()});
            const PlcGrowingAabbIndex::Box query{middle-Vec3{reach,reach,reach},middle+Vec3{reach,reach,reach}};
            const auto bounds=[&](std::size_t i){const auto p=vertices[i].position;return PlcGrowingAabbIndex::Box{p,p};};
            fixedIndex.all(fixedVertexCount,query,bounds,[&](std::uint64_t n){charge(n);},[&](std::size_t id){selected.push_back(static_cast<TetFacetVertexId>(id));return true;});
            if(selected.size()>256){std::nth_element(selected.begin(),selected.begin()+256,selected.end(),[&](auto x,auto y){charge(1);const auto a=(vertices[x].position-middle).norm2(),b=(vertices[y].position-middle).norm2();return a!=b?a<b:x<y;});selected.resize(256);}
            for(TetFacetVertexId id=fixedVertexCount;id<vertices.size();++id)selected.push_back(id);
        }
        charge(selected.size());
        const auto oldCapacity=row.apex.capacity();
        if(admitApexStorage)admitApexStorage(apexCapacityBytes+std::uint64_t(selected.size())*sizeof(decltype(row.apex)::value_type));
        row.apex.reserve(selected.size());apexCapacityBytes+=std::uint64_t(row.apex.capacity()-oldCapacity)*sizeof(decltype(row.apex)::value_type);
        for(TetFacetVertexId id:selected){if(std::find(f.begin(),f.end(),id)!=f.end())continue;const auto p=vertices[id].position;
            if(predicates::orient3d(a,b,c,p)!=predicates::Sign::Negative)continue;
            const auto height=(p-center).dot(inward);
            if(!(height>0))continue;
            const auto priority=((p-center).dot(p-center)-radius)/(2*height);if(std::isfinite(priority))row.apex.emplace_back(priority,id);}
        std::sort(row.apex.begin(),row.apex.end());
    }
    void indexManifoldFront(){
        if(manifoldIndexValid)return;
        frontVertexFaces.clear();for(const auto& [k,row]:front){(void)row;for(auto v:k)frontVertexFaces[v].insert(k);}
        manifoldIndexValid=true;
    }
    // Prevent an insertion from leaving chambers that meet only at an edge or
    // vertex. The remaining front must have a single cycle around each vertex.
    // This is a topology admission check; exact geometry checks still follow.
    bool manifoldAfter(const TetFacetCell& cell){
        indexManifoldFront();std::array<Face,4> changed;for(unsigned side=0;side<4;++side)changed[side]=key(tet_facet_detail::outwardFace(cell,side));
        for(auto vertex:cell.vertices){
            std::map<TetFacetVertexId,std::set<TetFacetVertexId>> link;
            auto add=[&](Face f){std::array<TetFacetVertexId,2> other;unsigned at=0;for(auto v:f)if(v!=vertex)other[at++]=v;if(at!=2)return false;return link[other[0]].insert(other[1]).second&&link[other[1]].insert(other[0]).second;};
            if(auto found=frontVertexFaces.find(vertex);found!=frontVertexFaces.end())for(auto f:found->second){if(std::find(changed.begin(),changed.end(),f)!=changed.end())continue;if(!add(f))return false;}
            for(auto f:changed)if(!front.contains(f)&&std::find(f.begin(),f.end(),vertex)!=f.end()&&!add(f))return false;
            if(link.empty())continue;
            for(const auto& [v,neighbors]:link){(void)v;if(neighbors.size()!=2)return false;}
            std::set<TetFacetVertexId> seen;std::vector<TetFacetVertexId> queue{link.begin()->first};seen.insert(queue.front());
            for(std::size_t i=0;i<queue.size();++i)for(auto v:link.at(queue[i]))if(seen.insert(v).second)queue.push_back(v);
            if(seen.size()!=link.size())return false;
        }
        return true;
    }
    bool valid(const Front& row,TetFacetVertexId apex,TetFacetCell& cell){
        charge();
        lastTopologyRefusal=false;
        if(cancellation.stop_requested())return false;
        ++metrics.apexTests;const auto f=row.directed;cell={{f[0],f[2],f[1],apex},cells.size()+1,{3,0}};
        if(predicates::orient3d(vertices[cell.vertices[0]].position,vertices[cell.vertices[1]].position,vertices[cell.vertices[2]].position,vertices[apex].position)!=predicates::Sign::Positive)return false;
        if(requireManifoldFront&&!manifoldAfter(cell)){lastTopologyRefusal=true;++metrics.nonmanifoldProposals;return false;}
        TetFacetComplexView view{vertices,cells,{}, {},{}, {3,0}};
        for(unsigned side=0;side<4;++side){const auto directed=tet_facet_detail::outwardFace(cell,side);const auto k=key(directed);
            if(auto i=front.find(k);i!=front.end()&&tet_facet_detail::faceSense(view,i->second.directed)!=tet_facet_detail::faceSense(view,directed))return false;}
        const auto bounds=box(cell.vertices);
        const std::array<Vec3,4> corners{vertices[cell.vertices[0]].position,vertices[cell.vertices[1]].position,vertices[cell.vertices[2]].position,vertices[apex].position};
        const PlcGrowingAabbIndex::Box query{bounds.lo,bounds.hi};
        const auto work=[&](std::uint64_t n){charge(n);};
        const auto fixedBounds=[&](std::size_t i){const auto p=vertices[i].position;return PlcGrowingAabbIndex::Box{p,p};};
        if(!fixedIndex.all(fixedVertexCount,query,fixedBounds,work,[&](std::size_t i){
            if(std::find(cell.vertices.begin(),cell.vertices.end(),i)!=cell.vertices.end())return true;
            if(i<sourceVertex.size()&&sourceVertex[i]){
                // A frozen boundary vertex cannot disappear inside a cell or
                // become an undeclared point on its edge/face. Check this
                // cheap witness before the complete source-facet contacts.
                for(unsigned side=0;side<4;++side){const auto f=tet_facet_detail::outwardFace(cell,side);
                    if(predicates::orient3d(vertices[f[0]].position,vertices[f[1]].position,vertices[f[2]].position,vertices[i].position)==predicates::Sign::Positive)return true;}
                return false;
            }
            return !tet_facet_detail::pointStrictlyInsideTet(vertices[i].position,corners);
        }))return false;
        if(!boundaryAccepts(cell,bounds,0))return false;
        const auto occupiedBounds=[&](std::size_t i){return PlcGrowingAabbIndex::Box{cellBoxes[i].lo,cellBoxes[i].hi};};
        if(!occupiedIndex.all(cells.size(),query,occupiedBounds,work,[&](std::size_t i){++metrics.pairTests;return !cellsConflict(view,cell,cells[i]);}))return false;
        if(interiorAdmission){
            if(apex>=fixedVertexCount&&!interiorAdmission(vertices[apex].position))return false;
        }
        if(!admitCentroid((corners[0]+corners[1]+corners[2]+corners[3])*.25))return false;
        return true;
    }
    void insert(const TetFacetCell& cell){
        if(requireManifoldFront)indexManifoldFront();else {manifoldIndexValid=false;frontVertexFaces.clear();}
        for(unsigned side=0;side<4;++side){auto directed=tet_facet_detail::outwardFace(cell,side);const auto k=key(directed);
            if(auto i=front.find(k);i!=front.end()){apexCapacityBytes-=std::uint64_t(i->second.apex.capacity())*sizeof(decltype(i->second.apex)::value_type);front.erase(i);if(requireManifoldFront)for(auto v:k){auto at=frontVertexFaces.find(v);at->second.erase(k);if(at->second.empty())frontVertexFaces.erase(at);}}
            else {std::swap(directed[1],directed[2]);front.emplace(k,Front{directed});if(requireManifoldFront)for(auto v:k)frontVertexFaces[v].insert(k);}}
        cells.push_back(cell);cellBoxes.push_back(box(cell.vertices));++metrics.inserted;
        // Geometry refusals are monotone as material is occupied. Topology
        // refusals are not: an adjacent insertion can reconnect a vertex link.
        if(requireManifoldFront&&retryTopologyRequired){for(auto& [k,row]:front){(void)k;if(row.topologyBlocked){row.next=0;row.topologyBlocked=false;}}retryTopologyRequired=false;}
    }
    bool initialize(){
        if(initialized)return true;
        if(boundary.empty()||vertices.size()<4||vertices.size()>=kInvalidTetFacetVertex)return false;
        for(std::size_t i=0;i<vertices.size();++i)if(!tet_facet_detail::finite(vertices[i].position)||(i&&vertices[i-1].sequence>=vertices[i].sequence))return false;
        std::set<Face> unique;
        std::map<std::array<TetFacetVertexId,2>,std::pair<unsigned,int>> edges;
        for(auto f:boundary){
            if(f[0]>=vertices.size()||f[1]>=vertices.size()||f[2]>=vertices.size()||f[0]==f[1]||f[1]==f[2]||f[2]==f[0]||!unique.insert(key(f)).second)return false;
            if(tet_facet_detail::projectionForTriangle(vertices[f[0]].position,vertices[f[1]].position,vertices[f[2]].position)==tet_facet_detail::Projection::Invalid)return false;
            for(unsigned j=0;j<3;++j){auto a=f[j],b=f[(j+1)%3];auto& edge=edges[{std::min(a,b),std::max(a,b)}];++edge.first;edge.second+=a<b?1:-1;}
        }
        for(const auto& [key,row]:edges){(void)key;if(row.first!=2||row.second!=0)return false;}
        sourceVertex.assign(vertices.size(),0);
        for(auto f:boundary){if(!front.emplace(key(f),Front{f}).second)return false;boundaryBoxes.push_back(box(f));for(auto v:f){sourceVertex[v]=1;fixedVertexCount=std::max(fixedVertexCount,TetFacetVertexId(v+1));}}
        order.resize(boundary.size());std::iota(order.begin(),order.end(),0);build(0,unsigned(order.size()));
        initialized=true;return true;
    }
    // Propose points from the remaining front's local clearance rather than a
    // global grid. Each point must form an exact positive, noncrossing cone
    // inside the unfilled complex; it is not a completed recovery certificate.
    // Only value-object vertices are appended. The original PLC is immutable.
    unsigned seedUnfilledFront(unsigned maximumSeeds=128){
        if(!initialized||front.empty()||cancellation.stop_requested()||maximumSeeds==0)return 0;
        std::vector<Face> snapshot;for(const auto& [k,row]:front){charge();(void)k;snapshot.push_back(row.directed);}
        Box bounds=box(snapshot.front());for(auto f:snapshot){const auto b=box(f);bounds.lo={std::min(bounds.lo.x,b.lo.x),std::min(bounds.lo.y,b.lo.y),std::min(bounds.lo.z,b.lo.z)};bounds.hi={std::max(bounds.hi.x,b.hi.x),std::max(bounds.hi.y,b.hi.y),std::max(bounds.hi.z,b.hi.z)};}
        const auto span=(bounds.hi-bounds.lo).norm();if(!(span>0&&std::isfinite(span)))return 0;
        unsigned inserted=0;
        for(auto f:snapshot){
            if(inserted>=maximumSeeds||cancellation.stop_requested())break;
            const auto a=vertices[f[0]].position,b=vertices[f[1]].position,c=vertices[f[2]].position;
            auto normal=(b-a).cross(c-a);const auto length=normal.norm();if(!(length>0&&std::isfinite(length)))continue;
            normal=normal*(-1./length);const auto center=(a+b+c)/3.;const auto end=center+normal*(2*span);double clearance=2*span;
            for(auto other:snapshot){
                charge();
                if(key(other)==key(f))continue;
                const auto p=vertices[other[0]].position,q=vertices[other[1]].position,r=vertices[other[2]].position;
                if(!tet_facet_detail::segmentCrossesOpenTriangle(center,end,p,q,r)&&!tet_facet_detail::segmentTouchesTriangleBoundary(center,end,p,q,r))continue;
                const auto n=(q-p).cross(r-p);const auto den=n.dot(normal);if(den==0)continue;
                const auto distance=n.dot(p-center)/den;
                if(distance>span*1e-12&&distance<clearance)clearance=distance;
            }
            for(double fraction:std::array{.25,.5,.1,.01,.001}){
                if(cancellation.stop_requested())break;
                ++metrics.frontSeedAttempts;const auto p=center+normal*(clearance*fraction);if(!tet_facet_detail::finite(p))continue;
                bool duplicate=false;for(const auto& v:vertices){charge();duplicate|=v.position.x==p.x&&v.position.y==p.y&&v.position.z==p.z;}if(duplicate)continue;
                const auto id=static_cast<TetFacetVertexId>(vertices.size());if(id==kInvalidTetFacetVertex)return inserted;
                vertices.push_back({p,vertices.size()+1});TetFacetCell cone;bool admissible=valid(Front{f},id,cone);
                if(admissible){TetFacetComplexView view{vertices,cells,{},{},{},{3,0}};
                    for(auto other:snapshot){charge();view.target.vertices=other;const auto contact=classifyTetFacetContact(view,cone);
                        if(contact!=TetFacetContact::Disjoint&&contact!=TetFacetContact::CoplanarExactFace&&contact!=TetFacetContact::DeclaredBoundaryTouch){admissible=false;break;}}
                }
                if(!admissible){vertices.pop_back();continue;}
                ++inserted;++metrics.frontSeeds;break;
            }
        }
        if(inserted){for(auto& [k,row]:front){(void)k;row.apex.clear();row.next=0;row.prepared=false;}exhaustedKernels.clear();}
        return inserted;
    }
    struct Plane {Vec3 normal;double offset;std::uint64_t order;};
    static std::uint64_t mix(std::uint64_t v){v+=0x9e3779b97f4a7c15ULL;v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;v=(v^(v>>27))*0x94d049bb133111ebULL;return v^(v>>31);}
    // Numerical feasibility is only a point proposal. All orientation,
    // boundary contact and overlap decisions remain exact on stored doubles.
    bool kernel(const std::vector<Face>& component,Vec3& result){
        Box bounds=box(component.front());for(auto f:component){const auto b=box(f);bounds.lo={std::min(bounds.lo.x,b.lo.x),std::min(bounds.lo.y,b.lo.y),std::min(bounds.lo.z,b.lo.z)};bounds.hi={std::max(bounds.hi.x,b.hi.x),std::max(bounds.hi.y,b.hi.y),std::max(bounds.hi.z,b.hi.z)};}
        const auto origin=(bounds.lo+bounds.hi)*.5;const auto scale=(bounds.hi-bounds.lo).norm();if(!(scale>0))return false;
        for(double margin:std::array{1e-8,1e-10,1e-12}){
            std::vector<Plane> planes{{{1,0,0},1,0},{{-1,0,0},1,1},{{0,1,0},1,2},{{0,-1,0},1,3},{{0,0,1},1,4},{{0,0,-1},1,5}};
            for(auto f:component){const auto a=vertices[f[0]].position,b=vertices[f[1]].position,c=vertices[f[2]].position;const auto n=(b-a).cross(c-a);const auto length=n.norm();if(!(length>0))return false;
                const auto unit=n/length;auto k=key(f);planes.push_back({unit,unit.dot((a-origin)/scale)-margin,mix(k[0])^mix(k[1]+1)^mix(k[2]+2)});}
            std::sort(planes.begin()+6,planes.end(),[](auto a,auto b){return a.order<b.order;});Vec3 point{};bool feasible=true;
            for(std::size_t i=0;i<planes.size()&&feasible;++i){const auto plane=planes[i];const auto excess=plane.normal.dot(point)-plane.offset;if(!(excess>0))continue;
                const Vec3 axis=std::abs(plane.normal.x)<=std::abs(plane.normal.y)&&std::abs(plane.normal.x)<=std::abs(plane.normal.z)?Vec3{1,0,0}:std::abs(plane.normal.y)<=std::abs(plane.normal.z)?Vec3{0,1,0}:Vec3{0,0,1};
                auto u=plane.normal.cross(axis);u=u/u.norm();const auto v=plane.normal.cross(u),base=point-plane.normal*excess;double s=0,t=0;
                for(std::size_t j=0;j<i&&feasible;++j){const auto previous=planes[j];const auto a=previous.normal.dot(u),b=previous.normal.dot(v),c=previous.offset-previous.normal.dot(base);if(!(a*s+b*t>c))continue;
                    const auto length=std::hypot(a,b);if(!(length>0)){feasible=false;break;}const auto nx=a/length,ny=b/length,height=c/length,sx=nx*height,sy=ny*height,dx=-ny,dy=nx;
                    double lower=-std::numeric_limits<double>::infinity(),upper=std::numeric_limits<double>::infinity();
                    for(std::size_t k=0;k<j;++k){const auto prior=planes[k];const auto ra=prior.normal.dot(u),rb=prior.normal.dot(v),rc=prior.offset-prior.normal.dot(base),coefficient=ra*dx+rb*dy,rhs=rc-ra*sx-rb*sy;
                        if(std::abs(coefficient)<=32*std::numeric_limits<double>::epsilon()){const auto uncertainty=64*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(rc),std::abs(ra*sx),std::abs(rb*sy)});if(rhs < -uncertainty){feasible=false;break;}continue;}
                        const auto bound=rhs/coefficient;if(!std::isfinite(bound)){feasible=false;break;}if(coefficient>0)upper=std::min(upper,bound);else lower=std::max(lower,bound);if(lower>upper){feasible=false;break;}}
                    if(!feasible)break;
                    const auto parameter=std::clamp(0.,lower,upper);
                    s=sx+dx*parameter;
                    t=sy+dy*parameter;
                }point=base+u*s+v*t;
            }
            if(!feasible)continue;
            const auto candidate=origin+point*scale;
            if(!tet_facet_detail::finite(candidate))continue;
            bool exact=true;for(auto f:component)if(predicates::orient3d(vertices[f[0]].position,vertices[f[1]].position,vertices[f[2]].position,candidate)!=predicates::Sign::Negative){exact=false;break;}
            if(exact){result=candidate;return true;}
        }return false;
    }
    bool finishKernels(unsigned maxCells=30000){
        if(cells.size()>maxCells||cancellation.stop_requested())return false;
        if(front.empty())return true;
        std::map<std::array<TetFacetVertexId,2>,std::vector<Face>> edgeFaces;for(const auto& [k,row]:front){for(unsigned i=0;i<3;++i){std::array<TetFacetVertexId,2> e{row.directed[i],row.directed[(i+1)%3]};std::sort(e.begin(),e.end());edgeFaces[e].push_back(k);}}
        std::set<Face> seen;std::vector<std::vector<Face>> components;
        for(const auto& [first,row]:front){(void)row;if(seen.contains(first))continue;std::vector<Face> keys{first};seen.insert(first);
            for(std::size_t at=0;at<keys.size();++at){const auto f=keys[at];for(unsigned i=0;i<3;++i){std::array<TetFacetVertexId,2> e{f[i],f[(i+1)%3]};std::sort(e.begin(),e.end());for(auto next:edgeFaces[e])if(seen.insert(next).second)keys.push_back(next);}}
            std::vector<Face> directed;for(auto k:keys)directed.push_back(front.at(k).directed);components.push_back(std::move(directed));}
        bool progress=false;
        for(const auto& component:components){if(cancellation.stop_requested())return false;if(component.size()<4||component.size()>maxCells-cells.size())continue;if(!exhaustedKernels.insert(component).second){++metrics.kernelCached;continue;}++metrics.kernelAttempts;Vec3 point;if(!kernel(component,point))continue;
            const auto id=TetFacetVertexId(vertices.size());vertices.push_back({point,vertices.size()+1});std::vector<TetFacetCell> proposed;bool accepted=true;
            for(auto f:component){TetFacetCell c;if(!valid(Front{f},id,c)){accepted=false;break;}proposed.push_back(c);}
            TetFacetComplexView view{vertices,{}, {},{}, {},{3,0}};
            if(accepted)for(std::size_t i=0;i<proposed.size()&&accepted;++i)for(std::size_t j=0;j<i;++j)if(tet_facet_detail::tetImproperlyIntersects(view,proposed[i],proposed[j])){accepted=false;break;}
            if(!accepted){vertices.pop_back();continue;}
            for(auto c:proposed)insert(c);
            progress=true;
            if(progressStream)*progressStream<<"{\"kernel_front_faces\":"<<component.size()<<",\"front_remaining\":"<<front.size()<<"}\n"<<std::flush;
        }return progress;
    }
    bool run(unsigned maxCells=30000,double maxSeconds=90,std::stop_token token={},std::uint64_t maxApexTests=1'000'000){
        cancellation=token;
        const auto start=std::chrono::steady_clock::now();
        if(!maxCells||!std::isfinite(maxSeconds)||maxSeconds<=0||!maxApexTests||cancellation.stop_requested())return false;
        if(!initialized&&!initialize())return false;
        // Rebuilding a large heap after every tetrahedron revisits the whole
        // surface. Retain unaffected rows in the same (next, face) order and
        // invalidate only the four changed faces. A failed batch rejoins at
        // the next round, exactly as in the rebuilding implementation.
        const bool incremental=useIncrementalFrontQueue&&!requireManifoldFront&&front.size()>=8192;
        using QueueEntry=std::tuple<std::size_t,Face,std::uint64_t>;
        std::vector<QueueEntry> pending;
        std::map<Face,std::uint64_t> tickets;
        std::vector<Face> deferred;
        std::uint64_t serial=0;
        const auto eligible=[](const Front& row){return !row.prepared||row.next<row.apex.size()||!row.full;};
        const auto enqueue=[&](Face f){
            tickets.erase(f);auto found=front.find(f);if(found==front.end()||!eligible(found->second))return;
            const auto ticket=++serial;tickets.emplace(f,ticket);pending.emplace_back(found->second.next,f,ticket);
            std::push_heap(pending.begin(),pending.end(),std::greater<>{});
        };
        const auto rebuild=[&](){
            pending.clear();tickets.clear();deferred.clear();pending.reserve(front.size());
            ++metrics.queueRebuilds;metrics.queueRowsScanned+=front.size();
            for(const auto& [f,row]:front)if(eligible(row)){const auto ticket=incremental?++serial:0;if(incremental)tickets.emplace(f,ticket);pending.emplace_back(row.next,f,ticket);}
            std::make_heap(pending.begin(),pending.end(),std::greater<>{});
        };
        if(incremental)rebuild();
        while(!front.empty()&&cells.size()<maxCells){bool progress=false,searchProgress=false;
            if(cancellation.stop_requested()||metrics.apexTests>=maxApexTests)return false;
            if(front.size()<400&&finishKernels(maxCells)){if(front.empty())return true;if(incremental)rebuild();}
            // Prefer unexamined local fronts. One difficult face must not
            // force a global apex scan before other faces can make progress.
            if(!incremental||pending.size()>3*front.size()+1024)rebuild();
            while(!pending.empty()){std::pop_heap(pending.begin(),pending.end(),std::greater<>{});const auto [priority,f,ticket]=pending.back();pending.pop_back();(void)priority;
                auto current=front.find(f);if(current==front.end())continue;
                if(incremental){const auto known=tickets.find(f);if(known==tickets.end()||known->second!=ticket)continue;}
                auto& row=current->second;if(!row.prepared)prepare(row);else if(!row.full&&row.next==row.apex.size())prepare(row,true);
                TetFacetCell cell;bool found=false;
                // Failed apices cannot become admissible after more material
                // is occupied; cache the exhausted prefix of each front row.
                unsigned batch=0;while(row.next<row.apex.size()&&batch++<16){if(cancellation.stop_requested()||metrics.apexTests>=maxApexTests)return false;searchProgress=true;const auto apex=row.apex[row.next++].second;if(requireManifoldFront&&row.geometryRejected.contains(apex))continue;if(valid(row,apex,cell)){found=true;break;}if(requireManifoldFront){if(lastTopologyRefusal){row.topologyBlocked=true;retryTopologyRequired=true;}else row.geometryRejected.insert(apex);}
                    if(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>maxSeconds)return false;}
                if(found){insert(cell);progress=true;
                    if(incremental)for(unsigned side=0;side<4;++side)enqueue(key(tet_facet_detail::outwardFace(cell,side)));
                    break;
                }++metrics.blocked;if(incremental)deferred.push_back(f);
            }
            if(incremental){for(auto f:deferred)enqueue(f);deferred.clear();}
            if(!progress&&!searchProgress)return false;
            if(progressStream&&cells.size()%100==0)*progressStream<<"{\"gift_cells\":"<<cells.size()<<",\"front\":"<<front.size()<<",\"apex_tests\":"<<metrics.apexTests<<"}\n"<<std::flush;
            if(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>maxSeconds)return false;
        }
        return front.empty();
    }
};
}
