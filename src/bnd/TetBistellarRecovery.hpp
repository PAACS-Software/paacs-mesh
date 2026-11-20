// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private, bounded boundary-recovery prototype. Bistellar changes retain the
// oriented patch boundary; they cannot carve material or change a PLC facet.
#include "paacs/bnd/TetFacetGeometry.hpp"
#include "PlcSpatialIndex.hpp"
#include <chrono>
#include <functional>
#include <stdexcept>
#include <unordered_set>

namespace paacs::bnd::detail {
struct BistellarMetrics {
    std::uint64_t proposed=0, positive=0, accepted=0, flips23=0, flips32=0,
        flips44=0, flipsNM=0, predicateCalls=0, neutral=0, rejectedProtected=0;
    std::uint64_t convexBipyramids=0;
    double longStarSeconds=0,validationSeconds=0,survivorSeconds=0;
};
struct BistellarTimer {
    double* seconds;
    std::chrono::steady_clock::time_point began;
    explicit BistellarTimer(double* target)
        : seconds(target),began(target?std::chrono::steady_clock::now():std::chrono::steady_clock::time_point{}) {}
    ~BistellarTimer(){if(seconds)*seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();}
};
class BistellarComplex {
public:
    using Face=std::array<TetFacetVertexId,3>;
    using Edge=std::array<TetFacetVertexId,2>;
    std::vector<TetFacetVertex> vertices;
    std::vector<TetFacetCell> cells;
    std::vector<unsigned char> alive;
    std::map<Face,std::set<std::size_t>> faces;
    std::map<Edge,std::set<std::size_t>> edges;
    std::set<Face> protectedFaces;
    std::set<Edge> protectedEdges;
    std::set<Face> requiredFaces;
    std::set<Edge> requiredEdges;
    std::function<void(std::uint64_t)> chargeWork;
    std::size_t maximumStoredCells=std::numeric_limits<std::size_t>::max();
    bool useLookahead=true;
    bool walkConstraints=false;
    bool profile=false;
    std::uint64_t edgeWalks=0,facetWalks=0,walkFallbacks=0;
    BistellarMetrics metrics;
    using Box=PlcGrowingAabbIndex::Box;
    std::vector<Box> cellBoxes;
    PlcGrowingAabbIndex index;
    static constexpr std::size_t noCell=std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> vertexCell;
    std::vector<std::array<std::size_t,4>> neighbours;
    void charge(std::uint64_t n=1) const {if(chargeWork)chargeWork(n);}
    template<std::size_t N> Box bounds(const std::array<TetFacetVertexId,N>& ids) const {
        auto p=vertices[ids[0]].position;Box box{p,p};
        for(auto id:ids){p=vertices[id].position;box.lo={std::min(box.lo.x,p.x),std::min(box.lo.y,p.y),std::min(box.lo.z,p.z)};
            box.hi={std::max(box.hi.x,p.x),std::max(box.hi.y,p.y),std::max(box.hi.z,p.z)};}
        return box;
    }
    static bool overlaps(Box a,Box b){
        return a.lo.x<=b.hi.x&&b.lo.x<=a.hi.x&&a.lo.y<=b.hi.y&&b.lo.y<=a.hi.y&&
               a.lo.z<=b.hi.z&&b.lo.z<=a.hi.z;
    }
    template<class Visitor> bool visit(Box box,Visitor visitor){
        const auto boxAt=[&](std::size_t i){return cellBoxes[i];};
        const auto work=[&](std::uint64_t n){charge(n);};
        const auto live=[&](std::size_t i){return !alive[i]||visitor(i);};
        return walkConstraints?index.allAppending(cells.size(),box,boxAt,work,live):
                               index.all(cells.size(),box,boxAt,work,live,2048);
    }
    static Face key(Face a){std::sort(a.begin(),a.end());return a;}
    static Edge key(Edge a){if(a[1]<a[0])std::swap(a[0],a[1]);return a;}
    static std::array<Face,4> faceKeys(const TetFacetCell& c){
        std::array<Face,4> r;for(unsigned i=0;i<4;++i)r[i]=key(tet_facet_detail::outwardFace(c,i));return r;
    }
    static std::array<Edge,6> edgeKeys(const TetFacetCell& c){
        std::array<Edge,6> r;unsigned k=0;for(unsigned i=0;i<4;++i)for(unsigned j=i+1;j<4;++j)r[k++]=key(Edge{c.vertices[i],c.vertices[j]});return r;
    }
    void insert(TetFacetCell c,bool account=true){
        if(cells.size()>=maximumStoredCells)throw std::length_error("constraint cell storage limit");
        if(account)charge(11);c.sequence=cells.size()+1;const auto id=cells.size();cells.push_back(c);alive.push_back(1);cellBoxes.push_back(bounds(c.vertices));
        const auto keys=faceKeys(c);
        if(walkConstraints){
            vertexCell.resize(vertices.size(),noCell);
            for(auto v:c.vertices)vertexCell[v]=id;
            neighbours.push_back({noCell,noCell,noCell,noCell});
        }
        for(unsigned side=0;side<4;++side){
            auto& owners=faces[keys[side]];
            if(walkConstraints&&!owners.empty()){
                const auto other=*owners.begin();
                const auto otherKeys=faceKeys(cells[other]);
                const auto otherSide=std::ranges::find(otherKeys,keys[side])-otherKeys.begin();
                neighbours[id][side]=other;neighbours[other][otherSide]=id;
            }
            owners.insert(id);
        }
        for(auto e:edgeKeys(c))edges[e].insert(id);
    }
    void erase(std::size_t id,bool account=true){
        if(account)charge(11);
        for(auto f:faceKeys(cells[id])){auto i=faces.find(f);i->second.erase(id);if(i->second.empty())faces.erase(i);}
        if(walkConstraints)for(auto other:neighbours[id])if(other!=noCell)
            for(auto& adjacent:neighbours[other])if(adjacent==id){adjacent=noCell;break;}
        for(auto e:edgeKeys(cells[id])){auto i=edges.find(e);i->second.erase(id);if(i->second.empty())edges.erase(i);}
        alive[id]=0;
    }
    void restore(std::size_t id){charge(11);alive[id]=1;for(auto f:faceKeys(cells[id]))faces[f].insert(id);for(auto e:edgeKeys(cells[id]))edges[e].insert(id);}
    bool positive(TetFacetCell& c){
        charge();++metrics.predicateCalls;
        const auto& v=c.vertices;const auto s=predicates::orient3d(vertices[v[0]].position,vertices[v[1]].position,vertices[v[2]].position,vertices[v[3]].position);
        if(s==predicates::Sign::Zero)return false;if(s==predicates::Sign::Negative)std::swap(c.vertices[0],c.vertices[1]);return true;
    }
    bool crosses(Face f,Edge target){
        charge();
        const bool first=std::find(f.begin(),f.end(),target[0])!=f.end(),last=std::find(f.begin(),f.end(),target[1])!=f.end();if(first&&last)return false;
        const auto p=vertices[target[0]].position,q=vertices[target[1]].position;
        const auto a=vertices[f[0]].position,b=vertices[f[1]].position,c=vertices[f[2]].position;
        if(std::max(p.x,q.x)<std::min({a.x,b.x,c.x})||std::min(p.x,q.x)>std::max({a.x,b.x,c.x})||
           std::max(p.y,q.y)<std::min({a.y,b.y,c.y})||std::min(p.y,q.y)>std::max({a.y,b.y,c.y})||
           std::max(p.z,q.z)<std::min({a.z,b.z,c.z})||std::min(p.z,q.z)>std::max({a.z,b.z,c.z}))return false;
        metrics.predicateCalls+=5;
        using predicates::Sign;
        const auto sp=predicates::orient3d(a,b,c,p),sq=predicates::orient3d(a,b,c,q);
        if(!first&&!last&&tet_facet_detail::opposite(sp,sq)){
            const auto ab=predicates::orient3d(p,q,a,b),bc=predicates::orient3d(p,q,b,c),ca=predicates::orient3d(p,q,c,a);
            return (ab!=Sign::Positive&&bc!=Sign::Positive&&ca!=Sign::Positive)||
                   (ab!=Sign::Negative&&bc!=Sign::Negative&&ca!=Sign::Negative);
        }
        if(sp!=Sign::Zero||sq!=Sign::Zero)return false;
        const auto projection=tet_facet_detail::projectionForTriangle(a,b,c);if(projection==tet_facet_detail::Projection::Invalid)return false;
        const auto pp=tet_facet_detail::project(p,projection),qq=tet_facet_detail::project(q,projection);
        for(unsigned i=0;i<3;++i)if(tet_facet_detail::segmentContact2(pp,qq,tet_facet_detail::project(vertices[f[i]].position,projection),tet_facet_detail::project(vertices[f[(i+1)%3]].position,projection))==tet_facet_detail::Segment2Contact::Proper)return true;
        return false;
    }
    bool patchConflict(const TetFacetComplexView& view,const TetFacetCell& a,const TetFacetCell& b) const {
        const auto tick=[&]{charge();};
        if(tet_facet_detail::separatedByPositiveTetFace(view,a,b.vertices,tick)||
           tet_facet_detail::separatedByPositiveTetFace(view,b,a.vertices,tick))return false;
        charge(32);
        return tet_facet_detail::tetImproperlyIntersects(view,a,b);
    }
    bool facetConflict(const TetFacetComplexView& view,const TetFacetCell& cell) const {
        if(tet_facet_detail::separatedByPositiveTetFace(view,cell,view.target.vertices,[&]{charge();}))return false;
        charge(16);
        // Candidate construction and Steiner insertion preserve the finite,
        // monotonically sequenced vertex table throughout recovery.
        const auto contact=tet_facet_detail::classifyTetFacetContactWithCanonicalVertices(view,cell);
        return contact!=TetFacetContact::Disjoint&&contact!=TetFacetContact::DeclaredBoundaryTouch&&
               contact!=TetFacetContact::CoplanarExactFace;
    }
    struct Proposal {
        std::vector<std::size_t> removed;
        std::vector<TetFacetCell> added;
        int delta=0;
        unsigned kind=0;
        bool convexBipyramid=false;
    };
    bool isConvexBipyramid(Edge apices,Face base) const {
        charge(5);
        return tet_facet_detail::segmentCrossesOpenTriangle(vertices[apices[0]].position,
            vertices[apices[1]].position,vertices[base[0]].position,vertices[base[1]].position,
            vertices[base[2]].position);
    }
    bool valid(Proposal& p,std::optional<Edge> target,TetFacetVertexId staged=kInvalidTetFacetVertex){
        const BistellarTimer timer{profile?&metrics.validationSeconds:nullptr};
        charge(1+16*(p.removed.size()+p.added.size()));++metrics.proposed;
        std::vector<TetFacetCell> old;std::set<TetFacetVertexId> oldVertices,newVertices;
        for(auto id:p.removed){if(id>=alive.size()||!alive[id])return false;old.push_back(cells[id]);for(auto v:cells[id].vertices)oldVertices.insert(v);}
        if(old.empty())return false;
        for(const auto& c:old)if(c.classification.dim!=old.front().classification.dim||c.classification.tag!=old.front().classification.tag)return false;
        for(auto& c:p.added){c.classification=old.front().classification;if(!positive(c))return false;for(auto v:c.vertices)newVertices.insert(v);}
        if(staged!=kInvalidTetFacetVertex&&(!newVertices.erase(staged)||oldVertices.contains(staged)))return false;
        if(oldVertices!=newVertices)return false;
        const TetFacetComplexView view{vertices,old,{}, {},{},old.front().classification};
        const auto oldFaces=tet_facet_detail::candidateFaceMap(view,old),newFaces=tet_facet_detail::candidateFaceMap(view,p.added);
        if(tet_facet_detail::boundarySenses(oldFaces)!=tet_facet_detail::boundarySenses(newFaces))return false;
        for(const auto& [f,uses]:newFaces){(void)f;if(uses.size()>2||(uses.size()==2&&uses[0].sense==uses[1].sense))return false;}
        std::set<Edge> newEdges;for(const auto& c:p.added)for(auto e:edgeKeys(c))newEdges.insert(e);
        for(const auto& [f,uses]:oldFaces)if(protectedFaces.contains(f.vertices)&&!newFaces.contains(f)){++metrics.rejectedProtected;return false;}
        for(const auto& c:old)for(auto e:edgeKeys(c))if(protectedEdges.contains(e)&&!newEdges.contains(e)){++metrics.rejectedProtected;return false;}
        if(!walkConstraints)
            for(std::size_t i=0;i<p.added.size();++i)for(std::size_t j=0;j<i;++j)
                if(patchConflict(view,p.added[i],p.added[j]))return false;
        p.delta=0;
        if(target){
            for(const auto& [f,u]:newFaces){(void)u;if(!oldFaces.contains(f))p.delta+=crosses(f.vertices,*target);}
            for(const auto& [f,u]:oldFaces){(void)u;if(!newFaces.contains(f))p.delta-=crosses(f.vertices,*target);}
        }
        ++metrics.positive;return true;
    }
    std::optional<Proposal> flip23(Face f,Edge target){
        return flip23(f,std::optional<Edge>{target});
    }
    std::optional<Proposal> flip23(Face f,std::optional<Edge> target){
        const auto it=faces.find(f);if(it==faces.end()||it->second.size()!=2||protectedFaces.contains(f))return {};
        Proposal p;p.removed.assign(it->second.begin(),it->second.end());p.kind=23;
        std::array<TetFacetVertexId,2> apex;
        for(unsigned i=0;i<2;++i)for(auto v:cells[p.removed[i]].vertices)if(std::find(f.begin(),f.end(),v)==f.end())apex[i]=v;
        if(apex[0]==apex[1]||edges.contains(key(apex)))return {};
        for(unsigned i=0;i<3;++i)p.added.push_back({{apex[0],apex[1],f[i],f[(i+1)%3]},0,{}});
        if(walkConstraints&&isConvexBipyramid(apex,f)){
            // Both triangulations cover this same convex bipyramid. Only the
            // base face disappears, and no existing edge is removed.
            const BistellarTimer timer{profile?&metrics.validationSeconds:nullptr};
            charge(1+16*(p.removed.size()+p.added.size()));
            ++metrics.proposed;
            const auto classification=cells[p.removed[0]].classification;
            if(cells[p.removed[1]].classification!=classification)return {};
            if(target)p.delta=-int(crosses(f,*target));
            for(unsigned i=0;i<3;++i){
                auto& cell=p.added[i];
                cell.classification=classification;
                if(!positive(cell))return {};
                if(target)p.delta+=crosses({apex[0],apex[1],f[i]},*target);
            }
            p.convexBipyramid=true;
            ++metrics.positive;
            return p;
        }
        if(!valid(p,target))return {};
        return p;
    }
    // Cheap direct attempts before considering intermediate flips. The only
    // 2-to-3 proposals that can create target use its endpoints as apices.
    bool recoverDirectEdge(Edge target){
        std::set<Face> localFaces;
        std::set<Edge> localEdges;
        const auto affected=edgeCavity(target);
        if(affected.size()>32)return false;
        for(auto id:affected){
            for(auto face:faceKeys(cells[id]))localFaces.insert(face);
            for(auto edge:edgeKeys(cells[id]))localEdges.insert(edge);
        }
        for(auto face:localFaces){
            charge();
            const auto& owners=faces.at(face);
            if(owners.size()!=2||protectedFaces.contains(face))continue;
            Edge apices{};
            unsigned i=0;
            for(auto id:owners){
                for(auto vertex:cells[id].vertices)
                    if(std::find(face.begin(),face.end(),vertex)==face.end()){apices[i++]=vertex;break;}
            }
            if(i!=2||key(apices)!=target)continue;
            if(auto patch=flip23(face,target);patch&&preservesSurvivors(*patch)){apply(*patch);return true;}
        }
        for(auto edge:localEdges){
            charge();
            const auto& owners=edges.at(edge);
            if(owners.size()!=4||protectedEdges.contains(edge))continue;
            bool first=false,last=false;
            for(auto id:owners)for(auto vertex:cells[id].vertices){
                first|=vertex==target[0];last|=vertex==target[1];
            }
            if(!first||!last)continue;
            for(auto& patch:edgeFlips(edge,target)){
                const bool creates=std::ranges::any_of(patch.added,[&](const TetFacetCell& cell){
                    const auto keys=edgeKeys(cell);return std::ranges::find(keys,target)!=keys.end();
                });
                if(creates&&preservesSurvivors(patch)){apply(patch);return true;}
            }
        }
        return false;
    }
    std::vector<Proposal> edgeFlips(Edge e,Edge target,const TetFacetComplexView* facet=nullptr){
        return edgeFlips(e,std::optional<Edge>{target},facet);
    }
    std::vector<Proposal> edgeFlips(Edge e,std::optional<Edge> target,const TetFacetComplexView* facet=nullptr){
        charge();
        std::vector<Proposal> out;const auto it=edges.find(e);if(it==edges.end()||protectedEdges.contains(e)||it->second.size()<3||it->second.size()>20)return out;
        std::map<TetFacetVertexId,std::set<TetFacetVertexId>> link;
        for(auto id:it->second){std::array<TetFacetVertexId,2> v;unsigned n=0;for(auto x:cells[id].vertices)if(x!=e[0]&&x!=e[1])v[n++]=x;if(n!=2)return out;link[v[0]].insert(v[1]);link[v[1]].insert(v[0]);}
        for(auto& [v,n]:link)if(n.size()!=2)return out;
        std::vector<TetFacetVertexId> ring{link.begin()->first};auto prev=ring.front(),cur=*link.begin()->second.begin();
        while(cur!=ring.front()){if(ring.size()>=link.size())return out;ring.push_back(cur);auto next=*link[cur].begin();if(next==prev)next=*link[cur].rbegin();prev=cur;cur=next;}
        if(ring.size()!=it->second.size())return out;
        Proposal base;base.removed.assign(it->second.begin(),it->second.end());
        if(ring.size()==3){
            base.kind=32;
            const Face triangle{ring[0],ring[1],ring[2]};
            for(auto v:e)base.added.push_back({{triangle[0],triangle[1],triangle[2],v},0,{}});
            if(valid(base,target)){
                if(walkConstraints)base.convexBipyramid=isConvexBipyramid(e,triangle);
                out.push_back(std::move(base));
            }
        }
        else if(ring.size()==4)for(unsigned i=0;i<2;++i){Proposal p=base;p.kind=44;auto a=ring[i],b=ring[(i+2)%4];if(edges.contains(key(Edge{a,b})))continue;
            for(auto v:e){p.added.push_back({{a,b,ring[(i+1)%4],v},0,{}});p.added.push_back({{a,b,ring[(i+3)%4],v},0,{}});}
            if(valid(p,target))out.push_back(std::move(p));}
        else {
            const BistellarTimer timer{profile?&metrics.longStarSeconds:nullptr};
            // Removing an n-valent edge is a triangulation of its cyclic
            // opposite-vertex link, with one cone from either endpoint. A
            // dynamic program examines O(n^3) link triangles, rather than
            // enumerating all tetrahedra on the cavity vertices.
            std::vector<TetFacetCell> old;for(auto id:base.removed)old.push_back(cells[id]);
            const TetFacetComplexView view{vertices,old,{}, {},{},old.front().classification};
            const auto boundary=tet_facet_detail::boundarySenses(tet_facet_detail::candidateFaceMap(view,old));
            std::vector<std::pair<Face,Box>> shell;
            for(const auto& [face,sense]:boundary){(void)sense;shell.emplace_back(face.vertices,bounds(face.vertices));}
            const unsigned n=unsigned(ring.size()),inf=1000000000;
            const auto clearOfShell=[&](const TetFacetCell& cell){
                auto contact=view;
                const auto box=bounds(cell.vertices);
                for(const auto& [face,faceBox]:shell){
                    if(!overlaps(box,faceBox))continue;
                    contact.target.vertices=face;
                    if(facetConflict(contact,cell))return false;
                }
                return true;
            };
            const auto triangle=[&](unsigned i,unsigned j,unsigned k,std::array<TetFacetCell,2>& t){
                const Face f{ring[i],ring[j],ring[k]};
                metrics.predicateCalls+=2;
                const auto a=predicates::orient3d(vertices[f[0]].position,vertices[f[1]].position,vertices[f[2]].position,vertices[e[0]].position),
                           b=predicates::orient3d(vertices[f[0]].position,vertices[f[1]].position,vertices[f[2]].position,vertices[e[1]].position);
                if(!tet_facet_detail::opposite(a,b))return false;
                for(unsigned h=0;h<2;++h){t[h]={{e[h],f[0],f[1],f[2]},0,old.front().classification};if(!positive(t[h]))return false;
                    for(unsigned side=0;side<4;++side){const auto directed=tet_facet_detail::outwardFace(t[h],side);const auto fk=tet_facet_detail::faceKey(directed);
                        if(auto found=boundary.find(fk);found!=boundary.end()&&found->second!=tet_facet_detail::faceSense(view,directed))return false;}
                }return true;
            };
            bool retry=false;
            const auto search=[&](bool checkAllTriangles){
                std::vector<std::vector<unsigned>> cost(n,std::vector<unsigned>(n,inf)),split(n,std::vector<unsigned>(n,n));
                for(unsigned i=0;i+1<n;++i)cost[i][i+1]=0;
                for(unsigned width=2;width<n;++width)for(unsigned i=0;i+width<n;++i){
                    const auto k=i+width;
                    for(unsigned j=i+1;j<k;++j){
                        charge();
                        if(cost[i][j]==inf||cost[j][k]==inf)continue;
                        std::array<TetFacetCell,2> t;
                        if(!triangle(i,j,k,t))continue;
                        if(checkAllTriangles&&(!clearOfShell(t[0])||!clearOfShell(t[1])))continue;
                        unsigned penalty=1;
                        for(const auto& c:t){
                            if(facet)penalty+=100*facetConflict(*facet,c);
                            else if(target)for(auto f:faceKeys(c))penalty+=100*crosses(f,*target);
                        }
                        const auto value=cost[i][j]+cost[j][k]+penalty;
                        if(value<cost[i][k]){cost[i][k]=value;split[i][k]=j;}
                    }
                }
                if(cost[0][n-1]==inf)return;
                Proposal patch=base;
                patch.kind=nmKind;
                const auto append=[&](auto&& self,unsigned i,unsigned k)->bool{
                    if(k-i<2)return true;
                    const auto j=split[i][k];
                    std::array<TetFacetCell,2> t;
                    if(!triangle(i,j,k,t))return false;
                    if(!checkAllTriangles&&(!clearOfShell(t[0])||!clearOfShell(t[1]))){retry=true;return false;}
                    patch.added.insert(patch.added.end(),t.begin(),t.end());
                    return self(self,i,j)&&self(self,j,k);
                };
                if(append(append,0,n-1)&&valid(patch,target))out.push_back(std::move(patch));
            };
            // If the lowest-cost oriented link passes the exact shell check,
            // filtering every other triangulation cannot improve its rank.
            search(!walkConstraints);
            if(retry)search(true);
        }
        return out;
    }
    static constexpr unsigned nmKind=99;
    bool preservesSurvivors(const Proposal& p){
        const BistellarTimer timer{profile?&metrics.survivorSeconds:nullptr};
        // An open crossing through the base proves the two triangulations
        // cover exactly the same convex bipyramid. The positive cells and
        // unchanged oriented boundary were already checked by valid().
        if(p.convexBipyramid){++metrics.convexBipyramids;return true;}
        const TetFacetComplexView view{vertices,{}, {},{},{},{3,0}};
        std::vector<Box> boxes;
        boxes.reserve(p.added.size());
        for(const auto& cell:p.added)boxes.push_back(bounds(cell.vertices));
        // The live scheduler ranks many provisional moves before selecting a
        // patch. Certify its internal contacts here, before any mutation.
        if(walkConstraints)
            for(std::size_t i=0;i<p.added.size();++i)for(std::size_t j=0;j<i;++j)
                if(overlaps(boxes[i],boxes[j])&&patchConflict(view,p.added[i],p.added[j]))return false;
        for(std::size_t i=0;i<p.added.size();++i)if(!visit(boxes[i],[&](std::size_t id){
            if(std::find(p.removed.begin(),p.removed.end(),id)!=p.removed.end())return true;
            return !patchConflict(view,p.added[i],cells[id]);
        }))return false;
        return true;
    }
    void apply(const Proposal& p){
        if(p.added.size()>maximumStoredCells-cells.size())throw std::length_error("constraint cell storage limit");
        charge(11*(p.removed.size()+p.added.size()));
        for(auto id:p.removed)erase(id,false);for(auto c:p.added){insert(c,false);
        for(auto f:faceKeys(c))if(requiredFaces.contains(f))protectedFaces.insert(f);
        for(auto e:edgeKeys(c))if(requiredEdges.contains(e))protectedEdges.insert(e);
    }++metrics.accepted;if(p.kind==23)++metrics.flips23;else if(p.kind==32)++metrics.flips32;else if(p.kind==44)++metrics.flips44;else ++metrics.flipsNM;if(p.delta==0)++metrics.neutral;}
    std::uint64_t signature(const Proposal& p){
        std::uint64_t h=1469598103934665603ULL;std::vector<std::array<TetFacetVertexId,4>> keys;
        for(const auto& c:p.added){auto k=c.vertices;std::sort(k.begin(),k.end());keys.push_back(k);}std::sort(keys.begin(),keys.end());
        for(const auto& k:keys)for(auto v:k){h^=v;h*=1099511628211ULL;}return h;
    }
    std::vector<Proposal> proposals(Edge target,const std::set<std::size_t>& nearby,bool longLinks=true){
        std::set<Face> localFaces;std::set<Edge> localEdges;
        for(auto id:nearby)if(id<alive.size()&&alive[id]){charge(11);for(auto f:faceKeys(cells[id]))localFaces.insert(f);for(auto e:edgeKeys(cells[id]))localEdges.insert(e);}
        std::vector<Proposal> out;for(auto f:localFaces)if(auto p=flip23(f,target))out.push_back(std::move(*p));
        for(auto e:localEdges)if(longLinks||edges.at(e).size()<=4)
            for(auto p:edgeFlips(e,target))out.push_back(std::move(p));
        return out;
    }
    // Every committed proposal preserves the complete temporary complex. The
    // caller owns publication and can discard this value object on refusal.
    std::set<std::size_t> edgeCavityIndex(Edge target){
        std::set<std::size_t> nearby;std::set<Face> visited;
        visit(bounds(target),[&](std::size_t id){for(auto f:faceKeys(cells[id]))if(visited.insert(f).second&&crosses(f,target)){
            const auto& uses=faces.at(f);nearby.insert(uses.begin(),uses.end());}return true;});
        return nearby;
    }
    // A transverse segment has one exit per visited tetrahedron. Shared
    // vertices, coplanar travel and edge hits fall back to exact indexed discovery.
    std::optional<std::set<std::size_t>> walkEdge(Edge target) {
        const auto p=vertices[target[0]].position,q=vertices[target[1]].position;
        const auto cut=[&](const Face& f){
            charge(4);
            return tet_facet_detail::segmentCrossesOpenTriangle(p,q,vertices[f[0]].position,
                vertices[f[1]].position,vertices[f[2]].position);
        };
        const auto anchor=vertexCell[target[0]];
        if(anchor==noCell||!alive[anchor])return {};
        std::vector<std::size_t> star{anchor};
        std::set<std::size_t> seen{anchor},path;
        std::size_t previous=noCell,current=noCell;
        for(std::size_t i=0;i<star.size()&&current==noCell;++i){
            if(star.size()>128)return {};
            const auto id=star[i];
            for(unsigned side=0;side<4;++side){
                charge();
                const auto adjacent=neighbours[id][side];
                if(adjacent==noCell)continue;
                if(cells[id].vertices[side]==target[0]){
                    if(cut(tet_facet_detail::outwardFace(cells[id],side))){
                        previous=id;current=adjacent;break;
                    }
                }else if(seen.insert(adjacent).second)star.push_back(adjacent);
            }
        }
        if(current==noCell)return {};
        path.insert(previous);
        for(unsigned step=0;step<64;++step){
            if(!path.insert(current).second)return {};
            if(std::ranges::find(cells[current].vertices,target[1])!=cells[current].vertices.end())return path;
            std::size_t next=noCell;
            for(unsigned side=0;side<4;++side){
                charge();
                const auto adjacent=neighbours[current][side];
                if(adjacent==noCell||adjacent==previous||!cut(tet_facet_detail::outwardFace(cells[current],side)))continue;
                if(next!=noCell)return {};
                next=adjacent;
            }
            if(next==noCell)return {};
            previous=current;current=next;
        }
        return {};
    }
    std::set<std::size_t> edgeCavity(Edge target){
        if(walkConstraints){
            if(auto walked=walkEdge(target)){++edgeWalks;return std::move(*walked);}
            ++walkFallbacks;
        }
        return edgeCavityIndex(target);
    }
    std::set<std::size_t> facetCavity(const TetFacetComplexView& view) {
        std::set<std::size_t> affected,seen;
        std::vector<std::size_t> queue;
        const auto enqueue=[&](std::size_t id){if(id!=noCell&&seen.insert(id).second)queue.push_back(id);};
        if(walkConstraints){
            const auto& f=view.target.vertices;
            for(unsigned i=0;i<3;++i)
                for(auto id:edges.at(key(Edge{f[i],f[(i+1)%3]})))enqueue(id);
            for(std::size_t i=0;i<queue.size()&&affected.size()<=32;++i){
                charge();
                const auto id=queue[i];
                if(!facetConflict(view,cells[id]))continue;
                affected.insert(id);
                for(auto adjacent:neighbours[id])enqueue(adjacent);
            }
            if(!affected.empty()){++facetWalks;return affected;}
            ++walkFallbacks;
        }
        visit(bounds(view.target.vertices),[&](std::size_t id){
            if(facetConflict(view,cells[id]))affected.insert(id);
            return affected.size()<=32;
        });
        return affected;
    }
    bool recoverEdge(Edge target,unsigned maximumSteps=1000,bool longLinks=true){
        target=key(target);std::unordered_set<std::uint64_t> tried;
        for(unsigned step=0;step<maximumSteps;++step){
            if(edges.contains(target)){protectedEdges.insert(target);return true;}
            const auto nearby=edgeCavity(target);
            if(nearby.empty()||(!useLookahead&&nearby.size()>32))return false;
            auto available=proposals(target,nearby,longLinks);
            std::stable_sort(available.begin(),available.end(),[](const Proposal& a,const Proposal& b){
                return std::tuple(a.delta,a.added.size(),a.kind)<std::tuple(b.delta,b.added.size(),b.kind);
            });
            const auto best=std::find_if(available.begin(),available.end(),[&](const Proposal& p){
                return p.delta<=0&&!tried.contains(signature(p))&&preservesSurvivors(p);
            });
            if(best==available.end()){
                if(!useLookahead||walkConstraints)return false;
                bool committed=false;unsigned attempts=0;for(const auto& p:available){if(p.delta<0||p.delta>3||tried.contains(signature(p))||++attempts>64)continue;
                    const auto begin=cells.size();for(auto id:p.removed)erase(id);for(auto c:p.added)insert(c);
                    std::set<std::size_t> expanded=nearby;for(std::size_t id=begin;id<cells.size();++id)expanded.insert(id);
                    auto second=proposals(target,expanded);std::optional<Proposal> winner;
                    for(auto& q:second)if(q.delta+p.delta<0&&!tried.contains(signature(q))&&(!winner||q.delta<winner->delta))winner=q;
                    for(std::size_t id=begin;id<cells.size();++id)erase(id);cells.resize(begin);alive.resize(begin);cellBoxes.resize(begin);for(auto id:p.removed)restore(id);
                    if(winner){tried.insert(signature(p));apply(p);
                        // Temporary and committed inserts have the same IDs.
                        tried.insert(signature(*winner));apply(*winner);committed=true;break;}
                }
                if(!committed)return false;
            }else {tried.insert(signature(*best));apply(*best);}
        }
        return edges.contains(target);
    }
};
}
