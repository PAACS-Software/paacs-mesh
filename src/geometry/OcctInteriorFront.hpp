// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "OcctDelaunayCavity.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/mesh/Quality.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace paacs::geo::detail {

struct OcctInteriorFrontLimits {
    std::uint64_t work=42'000'000,queries=12'000,nodes=50'000,tets=250'000;
    std::uint64_t workingBytes=256ULL*1024*1024;
    std::size_t maximumSites=1024,maximumCavityCells=64;
};

struct OcctInteriorFrontOutcome {
    bool initialized=false,changed=false,targetsSatisfied=false,budgetLimited=false;
    std::size_t sites=0,remainingOversized=0;
    std::uint64_t workUsed=0,nativeQueries=0,candidates=0,cavities=0;
    std::uint64_t frontActivations=0,frontAdvances=0,voidActivations=0;
    std::uint64_t spacingRejects=0,shapeRejects=0,cavityRejects=0;
    std::uint64_t validityRejects=0,nativeRejects=0;
    std::uint64_t spacingBoxPrunes=0,verifiedWitnesses=0;
    double minimumMeanRatio=0.,minimumScaledJacobian=0.;
};

namespace occt_interior_front_detail {
inline constexpr double separationFactor=.65,maximumSizeRatio=1.25;
inline constexpr double coneMeanRatioFloor=.01,coneJacobianFloor=1e-6;
struct Shape {double mean=0.,jacobian=0.;};
inline Shape shape(const std::array<Vec3,4>& p) {
    auto determinant=(p[1]-p[0]).cross(p[2]-p[0]).dot(p[3]-p[0]);
    if (std::isfinite(determinant)&&!(determinant>0.))
        determinant=predicates::orient3dDeterminantEstimate(p[0],p[1],p[2],p[3]);
    double sum=0.;
    for (unsigned i=0;i<4;++i) for (unsigned j=i+1;j<4;++j) sum+=(p[j]-p[i]).norm2();
    if (!(determinant>0.)||!(sum>0.)) return {};
    const auto mean=12.*std::pow(determinant/2.,2./3.)/sum;
    const auto jacobian=mesh::tetScaledJacobian(p[0],p[1],p[2],p[3]);
    return {std::isfinite(mean)?mean:0.,std::isfinite(jacobian)?jacobian:0.};
}
inline bool finite(Vec3 p) {
    return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);
}
// The edge metric integral for a linearly interpolated positive endpoint size
// is length/L(a,b). Actual CAD fields need not be linear: this is a symmetric
// ranking model, while actual point samples and conservative cell bounds still
// govern spacing and refinement discovery.
inline double logarithmicMeanSize(double a,double b) {
    if (!(a>0.)||!(b>0.)||!std::isfinite(a)||!std::isfinite(b)) return 0.;
    if (b<a) std::swap(a,b);
    if (a==b) return a;
    const auto difference=b-a;
    const auto relative=difference/a;
    const auto logarithm=std::isfinite(relative)?std::log1p(relative):std::log(b)-std::log(a);
    const auto mean=difference/logarithm;
    if (!(mean>0.)||!std::isfinite(mean)) return 0.;
    return std::clamp(mean,a,b);
}
} // namespace occt_interior_front_detail

// Construction on a complete, constraint-certified, positive material scaffold
// of one body. Tri3 cells protect its exact prescribed boundary. All original
// vertices remain fixed. SizeField supplies eval(p), minimumOnBox(lo,hi), and
// queryWork(); no control shapes or primitive types are inspected here.
//
// Front faces separate the boundary, an already-sized cell, or accepted cells
// from an oversized cell ahead. A separate conservative size-violation queue
// activates enclosed regions which the advancing front has not reached.
// Geometry is ranked on the actual connected replacement cavity. Every new
// cone must clear the shape floors above; these are local admission floors,
// not a promised global quality or sizing certificate for fixed boundaries.
//
// External callback exceptions propagate with the input mesh unchanged.
// Internal limits stop between fully admitted edits and return a valid partial
// construction, explicitly reporting targetsSatisfied=false when unresolved.
template<class SizeField,class AdmitNode,class AdmitCentroid,class Charge>
OcctInteriorFrontOutcome advanceOcctInteriorFront(mesh::Mesh& mesh,const SizeField& field,
    const OcctInteriorFrontLimits& limits,AdmitNode admitNode,
    AdmitCentroid admitCentroid,Charge charge) {
    using namespace occt_interior_front_detail;
    using Complex=bnd::detail::BistellarComplex;
    using Cell=bnd::TetFacetCell;
    using Vertex=bnd::TetFacetVertexId;
    using Sign=predicates::Sign;
    using bnd::tet_facet_detail::outwardFace;
    OcctInteriorFrontOutcome out;
    const auto initialTets=std::size_t(mesh.countCells(mesh::CellType::Tet4));
    const auto originalNodes=std::size_t(mesh.nodeCount());
    const auto maximumSites=std::min<std::size_t>(limits.maximumSites,1024);
    const auto cavityLimit=std::min<std::size_t>(limits.maximumCavityCells,64);
    if (!initialTets||!originalNodes||initialTets>limits.tets||originalNodes>limits.nodes||
        originalNodes>=bnd::kInvalidTetFacetVertex||!cavityLimit) return out;
    const auto maximumNodes=std::min({limits.nodes,std::uint64_t(originalNodes+maximumSites),
        std::uint64_t(bnd::kInvalidTetFacetVertex)});
    // Account for the live complex and its indices, two lazy heaps (including
    // rebuild headroom), node/edge caches, candidate patches and spatial search,
    // and the staged output while the caller retains the original input mesh.
    const std::uint64_t fixedBytes=32ULL*1024*1024+maximumNodes*256+
        std::uint64_t(mesh.cellCount())*64+std::uint64_t(cavityLimit)*16*1024;
    if (!limits.work||!limits.queries||fixedBytes>=limits.workingBytes) {
        out.budgetLimited=true;return out;
    }
    const auto maximumStored=std::min<std::uint64_t>({std::uint64_t(initialTets)+
        std::uint64_t(maximumSites)*(2*cavityLimit+2),
        (limits.workingBytes-fixedBytes)/8192,limits.work/128});
    const auto reserve=8*(maximumStored+maximumNodes+std::uint64_t(mesh.cellCount()));
    if (maximumStored<initialTets||reserve>=limits.work/2) {
        out.budgetLimited=true;return out;
    }
    struct Limit {};
    auto workLimit=limits.work-reserve;
    const auto tick=[&](std::uint64_t n=1) {
        if (out.workUsed>workLimit||n>workLimit-out.workUsed) throw Limit{};
        charge(n);out.workUsed+=n;
    };
    const auto sizeAt=[&](Vec3 p) {tick(field.queryWork());return field.eval(p);};
    const auto minimumOnBox=[&](Vec3 lo,Vec3 hi) {
        tick(field.queryWork());return field.minimumOnBox(lo,hi);
    };
    const auto admit=[&](auto& callback,Vec3 p) {
        if (out.nativeQueries>=limits.queries) throw Limit{};
        tick();++out.nativeQueries;
        const auto accepted=callback(p);
        tick(0);
        return accepted;
    };
    Complex live;
    live.walkConstraints=true;
    live.maximumStoredCells=std::size_t(maximumStored);
    live.chargeWork=tick;
    struct Info {Vec3 center{},lo{},hi{};double lower=0.,ratio=0.,centerSize=0.,frontTarget=0.,frontRatio=0.;};
    std::vector<Info> infos;
    std::vector<std::uint64_t> versions;
    std::vector<unsigned char> covered;
    std::vector<double> sizes;
    double maximumKnownSize=0.;
    auto activeTetCount=initialTets;
    std::map<Complex::Face,Complex::Face> boundaryWinding;
    mesh::Classification region;
    Vec3 lower=mesh.node(0),upper=lower;
    const auto positions=[&](const Cell& cell,Vec3 staged={}) {
        std::array<Vec3,4> p{};
        for (unsigned i=0;i<4;++i) p[i]=cell.vertices[i]==live.vertices.size()?
            staged:live.vertices[cell.vertices[i]].position;
        return p;
    };
    const auto vertexMinimum=[&](const Cell& cell,double stagedSize=std::numeric_limits<double>::infinity()) {
        auto h=stagedSize;
        for (const auto v:cell.vertices) if (v!=live.vertices.size()) h=std::min(h,sizes[v]);
        return h;
    };
    const auto infoFor=[&](const std::array<Vec3,4>& p,double vertexSize) {
        Info info;info.lo=info.hi=p[0];
        std::array<double,6> lengths{};unsigned next=0;
        for (unsigned i=0;i<4;++i) {
            info.center=info.center+p[i]*.25;
            info.lo={std::min(info.lo.x,p[i].x),std::min(info.lo.y,p[i].y),std::min(info.lo.z,p[i].z)};
            info.hi={std::max(info.hi.x,p[i].x),std::max(info.hi.y,p[i].y),std::max(info.hi.z,p[i].z)};
            for (unsigned j=0;j<i;++j) lengths[next++]=(p[i]-p[j]).norm();
        }
        info.lower=minimumOnBox(info.lo,info.hi);
        info.centerSize=sizeAt(info.center);
        if (!(info.lower>0.)||!std::isfinite(info.lower)||!(info.centerSize>0.)||!std::isfinite(info.centerSize))
            throw std::runtime_error("invalid interior-front size bound");
        std::ranges::sort(lengths);
        info.ratio=(lengths[2]+lengths[3])/(2*info.lower);
        info.frontTarget=std::min(vertexSize,info.centerSize);
        info.frontRatio=(lengths[2]+lengths[3])/(2*info.frontTarget);
        return info;
    };
    const auto oversized=[](const Info& info) {return info.ratio>maximumSizeRatio;};
    // Sampled readiness schedules the advancing front; it is not a proof that
    // enclosed field features are resolved. The independent void bound remains
    // active even when every existing point sample has a coarse target size.
    const auto frontOversized=[](const Info& info) {return info.frontRatio>maximumSizeRatio;};
    struct Front {
        double target,ratio;std::size_t ahead,behind;std::uint64_t version;Complex::Face face;
        bool operator<(const Front& other) const {
            if (target!=other.target) return target>other.target;
            return ratio==other.ratio?std::tie(ahead,face,behind)>
                std::tie(other.ahead,other.face,other.behind):ratio<other.ratio;
        }
    };
    struct Void {
        double target,ratio;std::size_t cell;std::uint64_t version;
        bool operator<(const Void& other) const {
            if (target!=other.target) return target>other.target;
            return ratio==other.ratio?cell>other.cell:ratio<other.ratio;
        }
    };
    std::priority_queue<Front> front;
    std::priority_queue<Void> voids;
    const auto enqueue=[&](std::size_t id) {
        tick();
        if (!live.alive[id]) return;
        const auto version=++versions[id];
        if (!oversized(infos[id])) return;
        voids.push({infos[id].lower,infos[id].ratio,id,version});
        if (covered[id]||!frontOversized(infos[id])) return;
        for (unsigned side=0;side<4;++side) {
            tick();
            const auto face=Complex::key(outwardFace(live.cells[id],side));
            const auto other=live.neighbours[id][side];
            if (other==Complex::noCell) {
                if (live.protectedFaces.contains(face))
                    front.push({infos[id].frontTarget,infos[id].frontRatio,id,other,version,face});
            } else if (covered[other]||!frontOversized(infos[other]))
                front.push({infos[id].frontTarget,infos[id].frontRatio,id,other,version,face});
        }
    };

    // Independent all-vertex proximity cache. It includes survivors outside
    // every proposed cavity, and uses symmetric .65*min(h(p),h(q)) spacing.
    using Key=std::array<std::int64_t,3>;
    std::map<Key,std::vector<Vertex>> buckets;
    double bucketWidth=0.;
    const auto keyAt=[&](Vec3 p,Key& key) {
        const auto v=(p-lower)/bucketWidth;
        const std::array values{v.x,v.y,v.z};
        const auto bound=double(std::numeric_limits<std::int64_t>::max()/4);
        for (unsigned i=0;i<3;++i) {
            if (!std::isfinite(values[i])||std::abs(values[i])>bound) return false;
            key[i]=std::int64_t(std::floor(values[i]));
        }
        return true;
    };
    const auto separated=[&](Vec3 p,double h) {
        const auto check=[&](Vertex v) {
            tick();
            const auto d=p-live.vertices[v].position;
            return std::hypot(d.x,d.y,d.z)>=separationFactor*std::min(h,sizes[v]);
        };
        Key key{};
        const auto reach=std::ceil(separationFactor*std::min(h,maximumKnownSize)/bucketWidth)+1;
        const auto width=2*reach+1;
        if (!keyAt(p,key)||!std::isfinite(reach)||reach>7||
            width*width*width>double(live.vertices.size())) {
            for (Vertex v=0;v<live.vertices.size();++v) if (!check(v)) return false;
            return true;
        }
        const auto r=int(reach);
        for (int x=-r;x<=r;++x) for (int y=-r;y<=r;++y) for (int z=-r;z<=r;++z) {
            tick();
            const auto it=buckets.find({key[0]+x,key[1]+y,key[2]+z});
            if (it!=buckets.end()) for (const auto v:it->second) if (!check(v)) return false;
        }
        return true;
    };
    const auto boxForbidden=[&](Vec3 lo,Vec3 hi,double lowerSize) {
        tick();
        const auto maximumRadius=separationFactor*std::min(lowerSize,maximumKnownSize);
        const auto diameterUpper=2*std::nextafter(maximumRadius,std::numeric_limits<double>::infinity());
        const auto widthLower=std::max({std::nextafter(hi.x-lo.x,0.),
            std::nextafter(hi.y-lo.y,0.),std::nextafter(hi.z-lo.z,0.)});
        // No ball of this maximum radius can strictly contain a box whose
        // projection on even one axis is already at least its diameter.
        if (widthLower>=diameterUpper) return false;
        const Vec3 center{std::midpoint(lo.x,hi.x),std::midpoint(lo.y,hi.y),std::midpoint(lo.z,hi.z)};
        const auto check=[&](Vertex v) {
            tick();const auto p=live.vertices[v].position;
            // Round the farthest-corner distance outwards and the forbidden
            // radius inwards. A positive strict result covers the entire box,
            // including survivors which are outside the proposed cavity.
            const auto outward=[](double x) {return std::nextafter(x,std::numeric_limits<double>::infinity());};
            const auto x=outward(std::max(std::abs(lo.x-p.x),std::abs(hi.x-p.x)));
            const auto y=outward(std::max(std::abs(lo.y-p.y),std::abs(hi.y-p.y)));
            const auto z=outward(std::max(std::abs(lo.z-p.z),std::abs(hi.z-p.z)));
            const auto radius=std::nextafter(separationFactor*std::min(lowerSize,sizes[v]),0.);
            return outward(std::hypot(x,y,z))<radius;
        };
        Key key{};
        const auto reach=std::ceil(maximumRadius/bucketWidth)+1;
        const auto width=2*reach+1;
        if (!keyAt(center,key)||!std::isfinite(reach)||reach>7||
            width*width*width>double(live.vertices.size())) {
            for (Vertex v=0;v<live.vertices.size();++v) if (check(v)) return true;
            return false;
        }
        const auto r=int(reach);
        for (int x=-r;x<=r;++x) for (int y=-r;y<=r;++y) for (int z=-r;z<=r;++z) {
            tick();const auto it=buckets.find({key[0]+x,key[1]+y,key[2]+z});
            if (it!=buckets.end()) for (const auto v:it->second) if (check(v)) return true;
        }
        return false;
    };
    const auto contains=[&](std::size_t id,Vec3 p) {
        for (unsigned side=0;side<4;++side) {
            tick();const auto f=outwardFace(live.cells[id],side);
            if (predicates::orient3dFiltered(live.vertices[f[0]].position,
                live.vertices[f[1]].position,live.vertices[f[2]].position,p)==Sign::Positive) return false;
        }
        return true;
    };
    const auto locate=[&](Vec3 p) {
        auto owner=Complex::noCell;
        live.visit({p,p},[&](std::size_t id) {
            if (!contains(id,p)) return true;
            owner=id;return false;
        });
        return owner;
    };
    std::map<Complex::Edge,double> edgeSizes;
    const auto edgeSize=[&](Vertex a,Vertex b) {
        const auto key=Complex::key(Complex::Edge{a,b});
        if (const auto it=edgeSizes.find(key);it!=edgeSizes.end()) return it->second;
        const auto h=sizeAt((live.vertices[a].position+live.vertices[b].position)*.5);
        edgeSizes.emplace(key,h);return h;
    };

    try {
        live.vertices.reserve(std::size_t(maximumNodes));
        infos.reserve(std::size_t(maximumStored));versions.reserve(std::size_t(maximumStored));
        covered.reserve(std::size_t(maximumStored));sizes.reserve(std::size_t(maximumNodes));
        for (Index i=0;i<mesh.nodeCount();++i) {
            tick();const auto p=mesh.node(i);
            if (!finite(p)) return out;
            const auto h=sizeAt(p);
            if (!(h>0.)||!std::isfinite(h)) return out;
            live.vertices.push_back({p,std::uint64_t(i)+1});sizes.push_back(h);
            maximumKnownSize=std::max(maximumKnownSize,h);
            lower={std::min(lower.x,p.x),std::min(lower.y,p.y),std::min(lower.z,p.z)};
            upper={std::max(upper.x,p.x),std::max(upper.y,p.y),std::max(upper.z,p.z)};
        }
        bucketWidth=(upper-lower).norm()/32.;
        if (!(bucketWidth>0.)||!std::isfinite(bucketWidth)) return out;
        for (Vertex v=0;v<live.vertices.size();++v) {
            tick();Key key{};if (!keyAt(live.vertices[v].position,key)) return out;
            buckets[key].push_back(v);
        }
        bool foundRegion=false;
        for (Index i=0;i<mesh.cellCount();++i) {
            tick();const auto& cell=mesh.cell(i);const auto* v=mesh.cellNodes(i);
            for (int j=0;j<mesh::nodesPerCell(cell.type);++j)
                if (v[j]<0||std::size_t(v[j])>=originalNodes) return out;
            if (cell.type==mesh::CellType::Tri3) {
                const Complex::Face oriented{Vertex(v[0]),Vertex(v[1]),Vertex(v[2])};
                const auto face=Complex::key(oriented);
                if (!live.protectedFaces.insert(face).second) return out;
                boundaryWinding.emplace(face,oriented);
                for (unsigned j=0;j<3;++j) live.protectedEdges.insert(Complex::key(
                    Complex::Edge{face[j],face[(j+1)%3]}));
            } else if (cell.type==mesh::CellType::Tet4) {
                if (!foundRegion) {region=cell.onEntity;foundRegion=true;}
                if (region.dim!=3||cell.onEntity!=region) return out;
                Cell tet{{Vertex(v[0]),Vertex(v[1]),Vertex(v[2]),Vertex(v[3])},0,region};
                const auto p=positions(tet);
                tick();
                if (predicates::orient3dFiltered(p[0],p[1],p[2],p[3])!=Sign::Positive) return out;
                live.insert(tet);infos.push_back(infoFor(p,vertexMinimum(tet)));versions.push_back(0);covered.push_back(0);
            } else return out;
        }
        for (const auto& [face,owners]:live.faces) {
            tick();if (owners.size()>2||(owners.size()==1&&!live.protectedFaces.contains(face))) return out;
        }
        for (const auto& face:live.protectedFaces) {
            tick();const auto it=live.faces.find(face);
            if (it==live.faces.end()||it->second.size()!=1) return out;
            bool sameWinding=false;
            const auto& owner=live.cells[*it->second.begin()];
            const auto& prescribed=boundaryWinding.at(face);
            for (unsigned side=0;side<4;++side) {
                const auto outward=outwardFace(owner,side);
                if (Complex::key(outward)!=face) continue;
                for (unsigned start=0;start<3;++start)
                    sameWinding=sameWinding || (outward[0]==prescribed[start] &&
                        outward[1]==prescribed[(start+1)%3] && outward[2]==prescribed[(start+2)%3]);
                break;
            }
            if (!sameWinding) return out;
        }
        out.initialized=true;
        for (std::size_t id=0;id<live.cells.size();++id) enqueue(id);

        while (out.sites<maximumSites&&live.vertices.size()<maximumNodes) {
            tick();
            std::size_t ahead=Complex::noCell;
            std::optional<Complex::Face> activeFace;
            while (!front.empty()) {
                tick();const auto item=front.top();
                if (!live.alive[item.ahead]||versions[item.ahead]!=item.version||covered[item.ahead]||
                    !frontOversized(infos[item.ahead])) {front.pop();continue;}
                const auto owners=live.faces.find(item.face);
                if (owners==live.faces.end()||!owners->second.contains(item.ahead)) {front.pop();continue;}
                if (item.behind!=Complex::noCell&&(!live.alive[item.behind]||
                    !owners->second.contains(item.behind)||
                    !(covered[item.behind]||!frontOversized(infos[item.behind])))) {front.pop();continue;}
                break;
            }
            while (!voids.empty()) {
                tick();const auto item=voids.top();
                if (!live.alive[item.cell]||versions[item.cell]!=item.version||!oversized(infos[item.cell])) {
                    voids.pop();continue;
                }
                break;
            }
            // Resolve the finest target first, then the largest normalized
            // defect, with ties favouring front propagation. Otherwise a large
            // coarse region can consume the budget needed by requested fine
            // sizing elsewhere, despite offering easier shape improvements.
            if (!front.empty()&&(voids.empty()||front.top().target<voids.top().target||
                (front.top().target==voids.top().target&&front.top().ratio>=voids.top().ratio))) {
                const auto item=front.top();front.pop();
                ahead=item.ahead;activeFace=item.face;++out.frontActivations;
            } else if (!voids.empty()) {
                ahead=voids.top().cell;voids.pop();++out.voidActivations;
            }
            if (ahead==Complex::noCell) break;
            const auto cellInfo=infos[ahead];
            auto localTarget=cellInfo.frontTarget;
            bool foundFineWitness=false;
            struct Candidate {Vec3 point;bool fineWitness;};
            std::vector<Candidate> candidates;
            const auto addCandidate=[&](Vec3 p,bool inAhead=false,bool fineWitness=false) {
                tick();
                if (candidates.size()>=16||!finite(p)||
                    (p-cellInfo.center).norm()>2*(cellInfo.hi-cellInfo.lo).norm()||
                    (inAhead&&!contains(ahead,p))) return;
                if (activeFace) {
                    for (unsigned side=0;side<4;++side) {
                        const auto f=outwardFace(live.cells[ahead],side);
                        if (Complex::key(f)!=*activeFace) continue;
                        tick();
                        if (predicates::orient3dFiltered(live.vertices[f[0]].position,
                            live.vertices[f[1]].position,live.vertices[f[2]].position,p)!=Sign::Negative) return;
                        break;
                    }
                }
                for (auto& candidate:candidates) if ((p-candidate.point).norm()<=bucketWidth*1e-10) {
                    candidate.fineWitness=candidate.fineWitness||fineWitness;return;
                }
                candidates.push_back({p,fineWitness});
            };
            if (!activeFace && cellInfo.lower<cellInfo.centerSize*.75) {
                // A conservative box bound is potential refinement, not a
                // promise that its value is attained inside the tetrahedron.
                // Only the void queue localizes fine sizing away from its
                // current center. Front propagation follows actual samples.
                // Binary boxes match the size oracle's domain. Exact source
                // halfspaces prune boxes wholly outside the real tetrahedron
                // before any potentially expensive field-bound query.
                struct SearchBox {
                    Vec3 lo,hi;double h;unsigned depth;std::uint64_t order;
                    bool operator<(const SearchBox& b) const {
                        return h==b.h?(depth==b.depth?order>b.order:depth<b.depth):h>b.h;
                    }
                };
                const auto intersectsAhead=[&](Vec3 lo,Vec3 hi) {
                    for (unsigned side=0;side<4;++side) {
                        const auto f=outwardFace(live.cells[ahead],side);
                        bool outside=true;
                        for (unsigned corner=0;corner<8;++corner) {
                            const Vec3 p{corner&1?hi.x:lo.x,corner&2?hi.y:lo.y,corner&4?hi.z:lo.z};
                            tick();
                            if (predicates::orient3dFiltered(live.vertices[f[0]].position,
                                live.vertices[f[1]].position,live.vertices[f[2]].position,p)!=Sign::Positive) {
                                outside=false;break;
                            }
                        }
                        if (outside) return false;
                    }
                    return true;
                };
                std::priority_queue<SearchBox> search;
                search.push({cellInfo.lo,cellInfo.hi,cellInfo.lower,0,0});
                std::uint64_t ordinal=1;unsigned visited=0;
                auto terminalLower=std::numeric_limits<double>::infinity();
                const auto unseenThreshold=cellInfo.centerSize*.75;
                while (!search.empty()&&visited++<96) {
                    tick();const auto part=search.top();search.pop();
                    if (boxForbidden(part.lo,part.hi,part.h)) {++out.spacingBoxPrunes;continue;}
                    if (part.h>=unseenThreshold) {
                        terminalLower=std::min(terminalLower,part.h);break;
                    }
                    const Vec3 center{std::midpoint(part.lo.x,part.hi.x),
                        std::midpoint(part.lo.y,part.hi.y),std::midpoint(part.lo.z,part.hi.z)};
                    const auto currentLower=std::min(terminalLower,part.h);
                    if (contains(ahead,center)) {
                        const auto actualSize=sizeAt(center);
                        if (actualSize<unseenThreshold&&actualSize<=currentLower*1.5) {
                            if (separated(center,actualSize)) {
                                addCandidate(center,true,true);++out.verifiedWitnesses;
                                foundFineWitness=true;localTarget=std::min(localTarget,actualSize);
                                terminalLower=std::min(terminalLower,part.h);break;
                            }
                            ++out.spacingRejects;
                        }
                    }
                    if (part.depth>=48) {
                        terminalLower=std::min(terminalLower,part.h);continue;
                    }
                    const auto extent=part.hi-part.lo;
                    const std::array length{extent.x,extent.y,extent.z};
                    const auto axis=std::size_t(std::max_element(length.begin(),length.end())-length.begin());
                    const std::array mid{center.x,center.y,center.z};
                    const std::array oldLo{part.lo.x,part.lo.y,part.lo.z},oldHi{part.hi.x,part.hi.y,part.hi.z};
                    if (!(mid[axis]>oldLo[axis]&&mid[axis]<oldHi[axis])) {
                        terminalLower=std::min(terminalLower,part.h);continue;
                    }
                    for (unsigned child=0;child<2;++child) {
                        auto lo=oldLo,hi=oldHi;
                        if (child) lo[axis]=mid[axis];else hi[axis]=mid[axis];
                        const Vec3 childLo{lo[0],lo[1],lo[2]},childHi{hi[0],hi[1],hi[2]};
                        if (!intersectsAhead(childLo,childHi)) continue;
                        // Parent lower size is valid throughout both children,
                        // so a covered child can be discarded before its field
                        // oracle is queried.
                        if (boxForbidden(childLo,childHi,part.h)) {++out.spacingBoxPrunes;continue;}
                        const auto h=std::max(part.h,minimumOnBox(childLo,childHi));
                        search.push({childLo,childHi,h,part.depth+1,ordinal++});
                    }
                }
                // Search bounds guide localization only; retain the original
                // input-cell lower bound for all reported target certificates.
            }
            // A low box bound alone does not justify creating a coarse point.
            // With no actual size defect or spaced fine witness, leave this
            // unresolved potential in the diagnostic instead of filling it.
            if (!activeFace&&!foundFineWitness&&!frontOversized(cellInfo)) continue;
            const auto offerFace=[&](Complex::Face key) {
                Complex::Face f{};bool found=false;
                for (unsigned side=0;side<4;++side) {
                    const auto face=outwardFace(live.cells[ahead],side);
                    if (Complex::key(face)==key) {f=face;found=true;break;}
                }
                if (!found) return;
                const auto a=live.vertices[f[0]].position,b=live.vertices[f[1]].position,c=live.vertices[f[2]].position;
                const auto u=b-a,v=c-a,n=u.cross(v);const auto n2=n.norm2();
                if (!(n2>0.)) return;
                auto base=(a+b+c)/3.;
                if (u.dot(v)>0.&&(a-b).dot(c-b)>0.&&(a-c).dot(b-c)>0.)
                    base=a+(v.cross(n)*u.norm2()+n.cross(u)*v.norm2())/(2*n2);
                const auto h=(edgeSize(f[0],f[1])+edgeSize(f[1],f[2])+edgeSize(f[2],f[0]))/3.;
                const auto r2=std::max({(base-a).norm2(),(base-b).norm2(),(base-c).norm2()});
                if (h>0.&&r2<h*h) {
                    const auto height=std::min(std::sqrt(h*h-r2),std::sqrt(2./3.)*h);
                    for (const auto factor:std::array{1.,.7,.45}) addCandidate(base-n*(height*factor/std::sqrt(n2)));
                }
                addCandidate((a+b+c)/3.*.35+cellInfo.center*.65);
            };
            if (activeFace) offerFace(*activeFace);
            else for (const auto face:Complex::faceKeys(live.cells[ahead])) offerFace(face);
            const auto p=positions(live.cells[ahead]);
            const auto a=p[1]-p[0],b=p[2]-p[0],c=p[3]-p[0];
            const auto determinant=2*predicates::orient3dDeterminantEstimate(p[0],p[1],p[2],p[3]);
            if (determinant>0.) addCandidate(p[0]+(b.cross(c)*a.norm2()+c.cross(a)*b.norm2()+a.cross(b)*c.norm2())/determinant);
            addCandidate(cellInfo.center);

            struct Plan {Vec3 point;double size,mean,jacobian,error,score;bool sizeAppropriate,fineWitness;Complex::Proposal patch;};
            std::vector<Plan> plans;
            for (const auto candidate:candidates) {
                const auto site=candidate.point;
                tick();++out.candidates;
                const auto h=sizeAt(site);
                if (!(h>0.)||!std::isfinite(h)) continue;
                if (!separated(site,h)) {++out.spacingRejects;continue;}
                const auto containing=locate(site);
                if (containing==Complex::noCell) {++out.cavityRejects;continue;}
                std::size_t owner=Complex::noCell;
                auto patch=buildOcctDelaunayCavity(live,containing,site,cavityLimit,owner,true,coneJacobianFloor);
                // The new site must actually replace the cell ahead of this
                // front face. Its position may cross several old cells; a
                // skinny initial cell must not cap the healthy front height.
                if (!patch||std::ranges::find(patch->removed,ahead)==patch->removed.end()) {
                    ++out.cavityRejects;continue;
                }
                ++out.cavities;
                if (patch->added.size()>live.maximumStoredCells-live.cells.size()) {out.budgetLimited=true;continue;}
                const auto retainedTets=activeTetCount-patch->removed.size();
                if (retainedTets>limits.tets||patch->added.size()>limits.tets-retainedTets) {
                    out.budgetLimited=true;continue;
                }
                double minMean=1.,minJacobian=1.,error=0.;
                std::map<Vertex,double> radial;
                bool acceptable=true;
                for (const auto& cell:patch->added) {
                    tick(4);const auto q=positions(cell,site);
                    const auto measured=shape(q);
                    if (measured.mean<coneMeanRatioFloor||measured.jacobian<coneJacobianFloor) {acceptable=false;break;}
                    minMean=std::min(minMean,measured.mean);minJacobian=std::min(minJacobian,measured.jacobian);
                    for (const auto v:cell.vertices) if (v!=live.vertices.size()&&!radial.contains(v)) {
                        const auto target=logarithmicMeanSize(h,sizes[v]);
                        if (!(target>0.)||!std::isfinite(target)) {acceptable=false;break;}
                        const auto ratio=(site-live.vertices[v].position).norm()/target;
                        radial.emplace(v,std::abs(std::log(ratio)));
                    }
                    if (!acceptable) break;
                }
                if (!acceptable) {++out.shapeRejects;continue;}
                for (const auto& [v,e]:radial) {(void)v;error+=e;}
                error/=double(std::max<std::size_t>(1,radial.size()));
                // Both terms are dimensionless. The geometric mean of the
                // worst mean ratio and Jacobian measures angular regularity;
                // mean absolute log endpoint-metric edge error penalizes
                // undersizing and oversizing symmetrically. Size influences every rank,
                // rather than only breaking exactly equal shape scores.
                const auto score=std::sqrt(minMean*minJacobian)/(1+error);
                plans.push_back({site,h,minMean,minJacobian,error,score,h<=maximumSizeRatio*localTarget,
                    candidate.fineWitness,std::move(*patch)});
            }
            std::stable_sort(plans.begin(),plans.end(),[](const Plan& a,const Plan& b) {
                // Requested local sizing takes precedence over a coarser point
                // with a higher shape score. Coarse fallback is considered only
                // after every size-appropriate plan fails full admission.
                if (a.sizeAppropriate!=b.sizeAppropriate) return a.sizeAppropriate;
                // A verified unresolved fine-field witness is construction
                // progress. A coarser high-shape-score plan cannot replace it
                // merely because its tetrahedra are easier to make regular.
                if (a.fineWitness!=b.fineWitness) return a.fineWitness;
                if (a.score!=b.score) return a.score>b.score;
                return a.mean==b.mean?(a.jacobian==b.jacobian?a.error<b.error:a.jacobian>b.jacobian):a.mean>b.mean;
            });
            bool committed=false;
            for (auto& plan:plans) {
                tick();
                std::vector<Info> addedInfos;addedInfos.reserve(plan.patch.added.size());
                for (const auto& cell:plan.patch.added)
                    addedInfos.push_back(infoFor(positions(cell,plan.point),vertexMinimum(cell,plan.size)));
                const auto staged=Vertex(live.vertices.size());
                live.vertices.push_back({plan.point,std::uint64_t(staged)+1});
                bool applied=false;
                try {
                    if (!live.valid(plan.patch,std::nullopt,staged)||!live.preservesSurvivors(plan.patch)) {
                        ++out.validityRejects;
                    } else if (!admit(admitNode,plan.point)) {
                        ++out.nativeRejects;
                    } else {
                        bool material=true;
                        for (const auto& cell:plan.patch.added) {
                            Vec3 center{};for (const auto v:cell.vertices) center=center+live.vertices[v].position*.25;
                            if (!admit(admitCentroid,center)) {material=false;break;}
                        }
                        if (material) {
                            const auto first=live.cells.size();
                            live.apply(plan.patch);applied=true;
                            activeTetCount=activeTetCount-plan.patch.removed.size()+plan.patch.added.size();
                            sizes.push_back(plan.size);maximumKnownSize=std::max(maximumKnownSize,plan.size);Key key{};
                            if (!keyAt(plan.point,key)) throw std::runtime_error("invalid interior-front point index");
                            buckets[key].push_back(staged);
                            infos.insert(infos.end(),addedInfos.begin(),addedInfos.end());
                            versions.resize(live.cells.size(),0);covered.resize(live.cells.size(),0);
                            ++out.sites;if (activeFace) ++out.frontAdvances;
                            std::set<std::size_t> affected;
                            for (auto id=first;id<live.cells.size();++id) {
                                covered[id]=1;affected.insert(id);
                                for (const auto other:live.neighbours[id]) if (other!=Complex::noCell) affected.insert(other);
                            }
                            for (const auto id:affected) enqueue(id);
                            committed=true;break;
                        }
                        ++out.nativeRejects;
                    }
                } catch (...) {
                    if (!applied) live.vertices.pop_back();
                    throw;
                }
                live.vertices.pop_back();
            }
            (void)committed;
            if (front.size()+voids.size()>maximumStored*8) {
                front={};voids={};
                for (std::size_t id=0;id<live.cells.size();++id) if (live.alive[id]) enqueue(id);
            }
        }
    } catch (const Limit&) {
        out.budgetLimited=true;
    }
    if (!out.initialized) return out;
    workLimit=limits.work;
    out.minimumMeanRatio=out.minimumScaledJacobian=1.;
    for (std::size_t id=0;id<live.cells.size();++id) {
        tick(4);if (!live.alive[id]) continue;
        out.remainingOversized+=oversized(infos[id]);
        const auto q=shape(positions(live.cells[id]));
        out.minimumMeanRatio=std::min(out.minimumMeanRatio,q.mean);
        out.minimumScaledJacobian=std::min(out.minimumScaledJacobian,q.jacobian);
    }
    out.targetsSatisfied=out.remainingOversized==0;
    if (!out.targetsSatisfied&&(out.sites>=maximumSites||live.vertices.size()>=maximumNodes)) out.budgetLimited=true;
    if (!out.sites) return out;
    mesh::Mesh result;
    for (Index i=0;i<mesh.nodeCount();++i) {tick();result.addNode(mesh.node(i),mesh.nodeClass(i));}
    for (std::size_t i=originalNodes;i<live.vertices.size();++i) {tick();result.addNode(live.vertices[i].position,region);}
    for (Index i=0;i<mesh.cellCount();++i) if (mesh.cell(i).type==mesh::CellType::Tri3) {
        tick();result.addCell(mesh::CellType::Tri3,std::span(mesh.cellNodes(i),3),mesh.cell(i).onEntity);
    }
    for (std::size_t id=0;id<live.cells.size();++id) if (live.alive[id]) {
        tick();const auto& cell=live.cells[id];const auto& v=cell.vertices;
        result.addCell(mesh::CellType::Tet4,{v[0],v[1],v[2],v[3]},cell.classification);
    }
    tick(0);mesh.swap(result);out.changed=true;
    return out;
}

} // namespace paacs::geo::detail
