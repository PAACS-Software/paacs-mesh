// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/geometry/OcctInteriorFront.hpp"
#include "../src/geometry/OcctModelInternal.hpp"
#include "../src/geometry/OcctSpatialSizing.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/mesh/Quality.hpp"
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

using namespace paacs;
using namespace paacs::geo::detail;

namespace {
constexpr mesh::Classification region{3,7};
using Coordinate=std::array<std::uint64_t,3>;

void check(bool value,const char* message) {
    if (!value) throw std::runtime_error(message);
}
Coordinate coordinate(Vec3 p) {
    return {std::bit_cast<std::uint64_t>(p.x),std::bit_cast<std::uint64_t>(p.y),
        std::bit_cast<std::uint64_t>(p.z)};
}
void sameMesh(const mesh::Mesh& a,const mesh::Mesh& b) {
    check(a.nodeCount()==b.nodeCount()&&a.cellCount()==b.cellCount(),"refusal changed mesh counts");
    for (Index i=0;i<a.nodeCount();++i)
        check(coordinate(a.node(i))==coordinate(b.node(i))&&a.nodeClass(i)==b.nodeClass(i),
            "refusal changed an original point or CAD classification");
    for (Index i=0;i<a.cellCount();++i) {
        const auto x=a.cell(i),y=b.cell(i);
        check(x.type==y.type&&x.onEntity==y.onEntity&&std::equal(a.cellNodes(i),
            a.cellNodes(i)+mesh::nodesPerCell(x.type),b.cellNodes(i)),
            "refusal changed original connectivity or entity identity");
    }
}
auto boundarySignature(const mesh::Mesh& mesh) {
    std::vector<std::tuple<int,int,std::array<Coordinate,3>>> triangles;
    std::vector<std::tuple<Coordinate,int,int>> nodes;
    std::set<Index> used;
    for (Index i=0;i<mesh.cellCount();++i) if (mesh.cell(i).type==mesh::CellType::Tri3) {
        const auto* v=mesh.cellNodes(i);
        std::array<Coordinate,3> f{coordinate(mesh.node(v[0])),coordinate(mesh.node(v[1])),coordinate(mesh.node(v[2]))};
        std::rotate(f.begin(),std::min_element(f.begin(),f.end()),f.end());
        triangles.emplace_back(mesh.cell(i).onEntity.dim,mesh.cell(i).onEntity.tag,f);
        used.insert(v,v+3);
    }
    for (const auto id:used) nodes.emplace_back(coordinate(mesh.node(id)),mesh.nodeClass(id).dim,mesh.nodeClass(id).tag);
    std::ranges::sort(triangles);std::ranges::sort(nodes);
    return std::pair{nodes,triangles};
}
void fixedBoundary(const mesh::Mesh& before,const mesh::Mesh& after) {
    check(after.nodeCount()>=before.nodeCount(),"front removed original vertices");
    for (Index i=0;i<before.nodeCount();++i)
        check(coordinate(before.node(i))==coordinate(after.node(i))&&before.nodeClass(i)==after.nodeClass(i),
            "front moved or reclassified an original fixed vertex");
    check(boundarySignature(before)==boundarySignature(after),"front changed the prescribed oriented classified boundary");
    check(mesh::auditLinearVolumeBoundary(after).valid,"front broke the embedded material boundary");
    const auto oldVolume=mesh::totalTetVolume(before),newVolume=mesh::totalTetVolume(after);
    check(std::abs(oldVolume-newVolume)<=1e-11*oldVolume,"front changed material volume");
}
double meanRatio(const std::array<Vec3,4>& p) {
    const auto determinant=predicates::orient3dDeterminantEstimate(p[0],p[1],p[2],p[3]);
    double sum=0.;
    for (unsigned i=0;i<4;++i) for (unsigned j=i+1;j<4;++j) sum+=(p[i]-p[j]).norm2();
    return determinant>0.&&sum>0.?12.*std::pow(determinant/2.,2./3.)/sum:0.;
}
void allAddedCones(const mesh::Mesh& mesh,Index originalNodes) {
    for (Index i=0;i<mesh.cellCount();++i) if (mesh.cell(i).type==mesh::CellType::Tet4) {
        const auto* v=mesh.cellNodes(i);
        const std::array p{mesh.node(v[0]),mesh.node(v[1]),mesh.node(v[2]),mesh.node(v[3])};
        check(mesh.cell(i).onEntity==region&&predicates::orient3d(p[0],p[1],p[2],p[3])==predicates::Sign::Positive,
            "front lost same-body exact positive orientation");
        if (std::ranges::any_of(std::span(v,4),[&](Index id){return id>=originalNodes;}))
            check(meanRatio(p)>=.01-1e-12&&mesh::tetScaledJacobian(p[0],p[1],p[2],p[3])>=1e-6-1e-14,
                "one cone of an admitted whole cavity is nearly coplanar or below its shape floor");
    }
}
template<class Field> void allVertexSeparation(const mesh::Mesh& mesh,Index originalNodes,const Field& field) {
    for (Index i=originalNodes;i<mesh.nodeCount();++i) for (Index j=0;j<i;++j) {
        const auto p=mesh.node(i),q=mesh.node(j),d=p-q;
        const auto distance=std::hypot(d.x,d.y,d.z);
        check(distance>=.65*std::min(field.eval(p),field.eval(q))*(1.-1e-12),
            "front site is too close to an existing live vertex outside or inside its cavity");
    }
}
void addClosedTet(mesh::Mesh& mesh,std::array<Index,4> vertices,int faceTag=101) {
    auto p=[&](unsigned i){return mesh.node(vertices[i]);};
    const auto sign=predicates::orient3d(p(0),p(1),p(2),p(3));
    check(sign!=predicates::Sign::Zero,"fixture tetrahedron is degenerate");
    if (sign==predicates::Sign::Negative) std::swap(vertices[0],vertices[1]);
    const bnd::TetFacetCell tet{{bnd::TetFacetVertexId(vertices[0]),bnd::TetFacetVertexId(vertices[1]),
        bnd::TetFacetVertexId(vertices[2]),bnd::TetFacetVertexId(vertices[3])},0,region};
    for (unsigned side=0;side<4;++side) {
        const auto f=bnd::tet_facet_detail::outwardFace(tet,side);
        mesh.addCell(mesh::CellType::Tri3,{Index(f[0]),Index(f[1]),Index(f[2])},{2,faceTag+int(side)});
    }
    mesh.addCell(mesh::CellType::Tet4,vertices,region);
}
mesh::Mesh regularTet(double scale=1.) {
    mesh::Mesh mesh;
    const std::array points{Vec3{0,0,0},Vec3{1,0,0},Vec3{.5,std::sqrt(3.)/2,0},
        Vec3{.5,std::sqrt(3.)/6,std::sqrt(2./3.)}};
    for (std::size_t i=0;i<points.size();++i) mesh.addNode(points[i]*scale,{0,int(i)+1});
    addClosedTet(mesh,{0,1,2,3});return mesh;
}
mesh::Mesh cubeScaffold() {
    mesh::Mesh mesh;
    std::array<Index,27> ids{};ids.fill(kInvalidIndex);
    const auto at=[](unsigned x,unsigned y,unsigned z){return x*9+y*3+z;};
    for (unsigned x=0;x<3;++x) for (unsigned y=0;y<3;++y) for (unsigned z=0;z<3;++z) {
        if (x==1&&y==1&&z==1) continue;
        ids[at(x,y,z)]=mesh.addNode({(double(x)-1)*.5,(double(y)-1)*.5,(double(z)-1)*.5},
            {0,int(mesh.nodeCount())+1});
    }
    const auto center=mesh.addNode({0,0,0},region);
    for (unsigned axis=0;axis<3;++axis) for (unsigned side=0;side<2;++side)
        for (unsigned u=0;u<2;++u) for (unsigned v=0;v<2;++v) {
            const auto vertex=[&](unsigned a,unsigned b) {
                std::array<unsigned,3> index{};index[axis]=side*2;
                index[(axis+1)%3]=a;index[(axis+2)%3]=b;
                return ids[at(index[0],index[1],index[2])];
            };
            const auto a=vertex(u,v),b=vertex(u+1,v),c=vertex(u+1,v+1),d=vertex(u,v+1);
            for (auto face:std::array{std::array<Index,3>{a,b,c},std::array<Index,3>{a,c,d}}) {
                if (predicates::orient3d(mesh.node(face[0]),mesh.node(face[1]),mesh.node(face[2]),mesh.node(center))==predicates::Sign::Positive)
                    std::swap(face[1],face[2]);
                mesh.addCell(mesh::CellType::Tri3,face,{2,int(axis*2+side)+1});
                mesh.addCell(mesh::CellType::Tet4,{center,face[0],face[1],face[2]},region);
            }
        }
    check(mesh::auditLinearVolumeBoundary(mesh).valid,"cube scaffold is not a conforming closed positive fill");
    return mesh;
}
struct ConstantField {
    double h;
    double eval(Vec3) const {return h;}
    double minimumOnBox(Vec3,Vec3) const {return h;}
    std::uint64_t queryWork() const {return 1;}
};
// An analytic bounded field with no control-container API. Its jump isolates
// enclosed-region discovery from a particular growth formula or CAD primitive.
class EnclosedField {
    Vec3 lo_,hi_;double fine_,coarse_;
public:
    mutable std::uint64_t boundQueries=0;
    mutable std::uint64_t pointQueries=0,finePointQueries=0;
    EnclosedField(Vec3 lo,Vec3 hi,double fine,double coarse):lo_(lo),hi_(hi),fine_(fine),coarse_(coarse) {}
    bool refined(Vec3 p) const {
        return p.x>=lo_.x&&p.x<=hi_.x&&p.y>=lo_.y&&p.y<=hi_.y&&p.z>=lo_.z&&p.z<=hi_.z;
    }
    double eval(Vec3 p) const {
        ++pointQueries;
        if (refined(p)) {++finePointQueries;return fine_;}
        return coarse_;
    }
    double minimumOnBox(Vec3 lo,Vec3 hi) const {
        ++boundQueries;
        return lo.x<=hi_.x&&hi.x>=lo_.x&&lo.y<=hi_.y&&hi.y>=lo_.y&&lo.z<=hi_.z&&hi.z>=lo_.z?fine_:coarse_;
    }
    std::uint64_t queryWork() const {return 1;}
};
class SmoothRegionField {
    sizing::RegionSizeField field_;
public:
    explicit SmoothRegionField(sizing::RegionSizeField field):field_(std::move(field)) {}
    double eval(Vec3 p) const {return field_.eval(p);}
    double minimumOnBox(Vec3 lo,Vec3 hi) const {return field_.minimumOnBox(lo,hi);}
    std::uint64_t queryWork() const {return field_.boxQueryWorkBound();}
};
OcctInteriorFrontLimits testLimits() {
    OcctInteriorFrontLimits limits;
    limits.work=4'000'000;limits.queries=4000;limits.nodes=256;limits.tets=4096;
    limits.workingBytes=64ULL*1024*1024;limits.maximumSites=4;limits.maximumCavityCells=64;
    return limits;
}
void outcomeBounds(const OcctInteriorFrontOutcome& out,const OcctInteriorFrontLimits& limits,
    const mesh::Mesh& before,const mesh::Mesh& after,std::uint64_t charged,std::uint64_t queries) {
    check(out.workUsed==charged&&out.workUsed<=limits.work&&out.nativeQueries==queries&&queries<=limits.queries,
        "front work or native admissions escaped their enclosing accounting");
    check(std::size_t(after.nodeCount()-before.nodeCount())==out.sites&&out.changed==(out.sites!=0)&&
        std::uint64_t(after.nodeCount())<=limits.nodes&&
        std::uint64_t(after.countCells(mesh::CellType::Tet4))<=limits.tets,
        "front site count, publication flag or output limits are inconsistent");
    check(!out.targetsSatisfied||out.remainingOversized==0,"front claimed completed sizing with unresolved cells");
}

void exerciseLogarithmicMetric() {
    using occt_interior_front_detail::logarithmicMeanSize;
    const std::array sizes{1e-250,1e-12,.06,.45,1.,1e12,1e250};
    for (const auto a:sizes) {
        check(logarithmicMeanSize(a,a)==a,"uniform endpoint metric changed its physical size");
        for (const auto b:sizes) {
            const auto mean=logarithmicMeanSize(a,b),reverse=logarithmicMeanSize(b,a);
            check(mean>0.&&std::isfinite(mean)&&mean>=std::min(a,b)*(1.-1e-12)&&
                mean<=std::max(a,b)*(1.+1e-12),"endpoint metric lost positivity, finiteness or endpoint bounds");
            check(std::abs(mean-reverse)<=1e-12*std::max(mean,reverse),
                "endpoint metric depends on edge orientation");
            for (const auto scale:{1e-50,.1,10.,1e50}) {
                const auto scaled=logarithmicMeanSize(a*scale,b*scale),expected=mean*scale;
                check(scaled>0.&&std::isfinite(scaled)&&std::abs(scaled-expected)<=2e-12*expected,
                    "endpoint metric changes under a physical unit conversion");
            }
        }
    }
    const auto adjacent=std::nextafter(1.,2.);
    check(std::abs(logarithmicMeanSize(1.,adjacent)-1.)<=4*std::numeric_limits<double>::epsilon(),
        "near-equal endpoint metric suffers cancellation");
    double previous=0.;
    for (const auto endpoint:{.2,.21,.5,2.,20.}) {
        const auto mean=logarithmicMeanSize(.2,endpoint);
        check(mean>=previous,"increasing an endpoint size increased its normalized metric distance");
        previous=mean;
    }
    // Independent quadrature of ds/h(s), with h varying linearly along this
    // edge. This validates the physical metric interpretation, not a claim
    // that an arbitrary CAD sizing field is linear between its endpoints.
    constexpr unsigned intervals=2048;
    constexpr double a=.15,b=.6,length=2.4;
    double quadrature=0.;
    for (unsigned i=0;i<=intervals;++i) {
        const auto h=std::lerp(a,b,double(i)/intervals);
        quadrature+=(i==0||i==intervals?1.:i%2?4.:2.)/h;
    }
    quadrature*=length/(3*intervals);
    const auto metric=length/logarithmicMeanSize(a,b);
    check(std::abs(metric-quadrature)<=1e-9*quadrature,
        "endpoint metric disagrees with the integrated linear-field distance");
}

void exerciseUniformAndEnclosed() {
    for (const auto scale:{1e-4,1.,1e4}) {
        auto mesh=regularTet(scale);const auto before=mesh;const ConstantField field{2*scale};
        std::uint64_t work=0,queries=0;
        const auto out=advanceOcctInteriorFront(mesh,field,testLimits(),
            [&](Vec3){++queries;return true;},[&](Vec3){++queries;return true;},
            [&](std::uint64_t n){work+=n;});
        check(out.initialized&&out.targetsSatisfied&&!out.changed&&queries==0,
            "a uniformly satisfied regular tetrahedron should need no proposals");
        sameMesh(before,mesh);outcomeBounds(out,testLimits(),before,mesh,work,queries);
    }
    {
        auto mesh=regularTet();const auto before=mesh;const ConstantField field{.2};
        auto limits=testLimits();limits.maximumSites=2;
        const auto out=advanceOcctInteriorFront(mesh,field,limits,
            [](Vec3){return true;},[](Vec3){return true;},[](std::uint64_t){});
        check(out.initialized&&out.sites>0&&out.frontAdvances>0&&out.voidActivations>0,
            "accepted front cells must feed the conservative remaining-volume activation queue");
        fixedBoundary(before,mesh);allAddedCones(mesh,before.nodeCount());
        allVertexSeparation(mesh,before.nodeCount(),field);
    }
    const ConstantField field{.45};
    auto mesh=cubeScaffold();const auto before=mesh;const auto limits=testLimits();
    std::uint64_t work=0,queries=0;
    const auto out=advanceOcctInteriorFront(mesh,field,limits,
        [&](Vec3){++queries;return true;},[&](Vec3){++queries;return true;},
        [&](std::uint64_t n){work+=n;});
    check(out.initialized&&out.sites>=2&&out.frontActivations>0&&out.frontAdvances>0,
        "uniform scaffold must actually advance an interior front through admitted cavities");
    outcomeBounds(out,limits,before,mesh,work,queries);fixedBoundary(before,mesh);
    allAddedCones(mesh,before.nodeCount());allVertexSeparation(mesh,before.nodeCount(),field);

    auto enclosed=cubeScaffold();const auto original=enclosed;
    const EnclosedField hidden({.16,.07,-.13},{.20,.11,-.09},.06,2.);
    for (Index i=0;i<original.cellCount();++i) if (original.cell(i).type==mesh::CellType::Tet4) {
        const auto* v=original.cellNodes(i);Vec3 center{};
        for (unsigned j=0;j<4;++j) center=center+original.node(v[j])*.25;
        check(hidden.eval(center)==2.,"hidden control unexpectedly touches an old cell centroid");
    }
    work=queries=0;
    const auto seeded=advanceOcctInteriorFront(enclosed,hidden,limits,
        [&](Vec3){++queries;return true;},[&](Vec3){++queries;return true;},
        [&](std::uint64_t n){work+=n;});
    bool found=false;
    for (Index i=original.nodeCount();i<enclosed.nodeCount();++i) found|=hidden.refined(enclosed.node(i));
    if (!(seeded.initialized&&seeded.sites>0&&found&&hidden.boundQueries>0&&
        seeded.frontActivations+seeded.voidActivations>0)) {
        std::cerr<<std::setprecision(17)<<"hidden front: sites="<<seeded.sites
            <<" maxSites="<<limits.maximumSites<<" front="<<seeded.frontActivations
            <<" advances="<<seeded.frontAdvances<<" void="<<seeded.voidActivations
            <<" candidates="<<seeded.candidates<<" cavities="<<seeded.cavities
            <<" work="<<seeded.workUsed<<" native="<<seeded.nativeQueries
            <<" budget="<<seeded.budgetLimited<<" remaining="<<seeded.remainingOversized
            <<" shapeRejects="<<seeded.shapeRejects<<" spacingRejects="<<seeded.spacingRejects
            <<" cavityRejects="<<seeded.cavityRejects<<" validityRejects="<<seeded.validityRejects
            <<" bounds="<<hidden.boundQueries<<" evals="<<hidden.pointQueries
            <<" fineEvals="<<hidden.finePointQueries<<'\n';
        for (Index i=original.nodeCount();i<enclosed.nodeCount();++i) {
            const auto p=enclosed.node(i);
            std::cerr<<"hidden inserted "<<i<<": "<<p.x<<','<<p.y<<','<<p.z
                <<" h="<<hidden.eval(p)<<" inControl="<<hidden.refined(p)<<'\n';
        }
        unsigned tet=0;
        for (Index i=0;i<original.cellCount();++i) if (original.cell(i).type==mesh::CellType::Tet4) {
            const auto* v=original.cellNodes(i);Vec3 center{},lo=original.node(v[0]),hi=lo;
            double minimumVertex=2.;bool contains=true;
            const bnd::TetFacetCell cell{{bnd::TetFacetVertexId(v[0]),bnd::TetFacetVertexId(v[1]),
                bnd::TetFacetVertexId(v[2]),bnd::TetFacetVertexId(v[3])},0,region};
            for (unsigned j=0;j<4;++j) {
                const auto p=original.node(v[j]);center=center+p*.25;minimumVertex=std::min(minimumVertex,hidden.eval(p));
                lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};
                hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};
                const auto f=bnd::tet_facet_detail::outwardFace(cell,j);
                if (predicates::orient3d(original.node(f[0]),original.node(f[1]),original.node(f[2]),Vec3{.18,.09,-.11})==predicates::Sign::Positive)
                    contains=false;
            }
            if (contains) std::cerr<<"hidden original owner tet="<<tet<<" meshCell="<<i
                <<" vertexMin="<<minimumVertex<<" centroidH="<<hidden.eval(center)
                <<" bound="<<hidden.minimumOnBox(lo,hi)<<'\n';
            ++tet;
        }
    }
    check(seeded.initialized&&seeded.sites>0&&found&&hidden.boundQueries>0&&
        seeded.frontActivations+seeded.voidActivations>0,
        "conservative-field activation failed to find the enclosed control missed by every old centroid");
    outcomeBounds(seeded,limits,original,enclosed,work,queries);fixedBoundary(original,enclosed);
    allAddedCones(enclosed,original.nodeCount());allVertexSeparation(enclosed,original.nodeCount(),hidden);
}

void exerciseWholeCavityShape() {
    auto mesh=cubeScaffold();const auto before=mesh;
    for (Index i=0;i<before.cellCount();++i) if (before.cell(i).type==mesh::CellType::Tet4) {
        const auto* v=before.cellNodes(i);
        check(meanRatio({before.node(v[0]),before.node(v[1]),before.node(v[2]),before.node(v[3])})>.1,
            "prospective-cavity fixture must begin with regular owner cells");
    }
    // All old cells are regular. A localized target immediately beside a
    // fixed coarse face can only make a nearly flat cone on that face. Scoring
    // the owner alone, or only the best cone, would miss this refusal.
    const EnclosedField field({.20,.11,-.49996},{.24,.15,-.49982},1e-4,2.);
    auto limits=testLimits();limits.work=2'000'000;
    std::uint64_t queries=0,work=0;
    const auto out=advanceOcctInteriorFront(mesh,field,limits,
        [&](Vec3){++queries;return true;},[&](Vec3){++queries;return true;},
        [&](std::uint64_t n){work+=n;});
    check(out.initialized&&out.shapeRejects>0&&out.sites==0&&queries==0&&!out.targetsSatisfied,
        "near-coplanar prospective shell cones must be refused before native admission");
    sameMesh(before,mesh);outcomeBounds(out,limits,before,mesh,work,queries);
}

void exerciseFrontAndVoidFairness() {
    auto mesh=cubeScaffold();const auto before=mesh;
    // The .45 background already drives the ordinary boundary front. The
    // independent .06 interior request must also receive an admitted site
    // within the same small allowance, rather than waiting for that front to
    // drain. The field exposes only the generic point/bound query interface.
    const EnclosedField field({.16,.07,-.13},{.20,.11,-.09},.06,.45);
    auto limits=testLimits();limits.maximumSites=4;
    std::uint64_t work=0,queries=0;
    const auto out=advanceOcctInteriorFront(mesh,field,limits,
        [&](Vec3){++queries;return true;},[&](Vec3){++queries;return true;},
        [&](std::uint64_t n){work+=n;});
    outcomeBounds(out,limits,before,mesh,work,queries);
    fixedBoundary(before,mesh);allAddedCones(mesh,before.nodeCount());
    allVertexSeparation(mesh,before.nodeCount(),field);
    bool found=false;std::size_t firstFine=0;
    for (Index i=before.nodeCount();i<mesh.nodeCount();++i) if (field.refined(mesh.node(i))) {
        found=true;firstFine=std::size_t(i-before.nodeCount())+1;break;
    }
    if (!found) {
        std::cerr<<std::setprecision(17)<<"front/void fairness: sites="<<out.sites
            <<" maxSites="<<limits.maximumSites<<" front="<<out.frontActivations
            <<" advances="<<out.frontAdvances<<" void="<<out.voidActivations
            <<" work="<<out.workUsed<<" native="<<out.nativeQueries
            <<" budget="<<out.budgetLimited<<" remaining="<<out.remainingOversized
            <<" fineEvals="<<field.finePointQueries<<'\n';
        for (Index i=before.nodeCount();i<mesh.nodeCount();++i) {
            const auto p=mesh.node(i);
            std::cerr<<"fairness inserted "<<i<<": "<<p.x<<','<<p.y<<','<<p.z
                <<" h="<<field.eval(p)<<'\n';
        }
    }
    check(out.initialized&&out.sites>0&&out.sites<=limits.maximumSites&&found&&
        firstFine<=limits.maximumSites,
        "active coarse front starved enclosed refinement within the four-site allowance");
}

void exerciseSmoothRegionalCoverage() {
    constexpr double target=.15;
    const Vec3 lower{-.15,-.2,-.2},upper{.25,.2,.2};
    sizing::RegionSizeSpecification specification;
    specification.baseSize=.45;specification.minimumSize=target;
    specification.controls.push_back({sizing::SizeRegionKind::Box,lower,upper,1.,target,1.55});
    auto compiled=sizing::compileRegionSizeField(specification);
    check(compiled.compiled,"smooth regional field failed compilation");
    const SmoothRegionField field(std::move(compiled.field));
    const auto inside=[&](Vec3 p) {
        return p.x>=lower.x&&p.x<=upper.x&&p.y>=lower.y&&p.y<=upper.y&&p.z>=lower.z&&p.z<=upper.z;
    };
    auto mesh=cubeScaffold();const auto before=mesh;auto limits=testLimits();
    limits.maximumSites=64;
    std::uint64_t work=0,queries=0;
    const auto out=advanceOcctInteriorFront(mesh,field,limits,
        [&](Vec3){++queries;return true;},[&](Vec3){++queries;return true;},
        [&](std::uint64_t n){work+=n;});
    outcomeBounds(out,limits,before,mesh,work,queries);
    fixedBoundary(before,mesh);allAddedCones(mesh,before.nodeCount());
    allVertexSeparation(mesh,before.nodeCount(),field);
    std::vector<double> cellSizes;
    for (Index i=0;i<mesh.cellCount();++i) if (mesh.cell(i).type==mesh::CellType::Tet4) {
        const auto* v=mesh.cellNodes(i);Vec3 center{};
        for (unsigned j=0;j<4;++j) center=center+mesh.node(v[j])*.25;
        if (!inside(center)) continue;
        std::array<double,6> edges{};unsigned next=0;
        for (unsigned a=0;a<4;++a) for (unsigned b=a+1;b<4;++b)
            edges[next++]=(mesh.node(v[a])-mesh.node(v[b])).norm();
        std::ranges::sort(edges);cellSizes.push_back((edges[2]+edges[3])*.5);
    }
    std::ranges::sort(cellSizes);
    const auto count=cellSizes.size();
    const auto median=count?(count%2?cellSizes[count/2]:
        (cellSizes[count/2-1]+cellSizes[count/2])*.5):std::numeric_limits<double>::infinity();
    std::size_t regionalNodes=0;
    for (const auto p:mesh.nodes()) regionalNodes+=inside(p);
    std::cout<<std::setprecision(17)<<"smooth regional coverage: sites="<<out.sites
        <<" maxSites="<<limits.maximumSites<<" roiNodes="<<regionalNodes<<" roiTets="<<count
        <<" medianSixEdges="<<median<<" target="<<target<<" ratio="<<median/target
        <<" front="<<out.frontActivations<<" advances="<<out.frontAdvances<<" void="<<out.voidActivations
        <<" work="<<out.workUsed<<" native="<<out.nativeQueries
        <<" budget="<<out.budgetLimited<<" remaining="<<out.remainingOversized<<'\n';
    // This region is at least .25 from the fixed shell: an achieved .2025
    // characteristic size does not require changing its coarse boundary.
    // Measure actual Tet4 edges, independently of node count or field internals.
    check(out.initialized&&count>0&&median<=1.35*target,
        "smooth interior region remains too coarse despite its finer requested size");
}

void exerciseAllLiveVertexSpacing() {
    const ConstantField field{.45};auto limits=testLimits();limits.maximumSites=1;
    auto reference=cubeScaffold();const auto originalNodes=reference.nodeCount();
    const auto first=advanceOcctInteriorFront(reference,field,limits,
        [](Vec3){return true;},[](Vec3){return true;},[](std::uint64_t){});
    check(first.sites==1,"spacing fixture needs one genuine front proposal");
    const auto p=reference.node(originalNodes);
    const std::array coordinates{p.x,p.y,p.z};
    unsigned axis=0;for (unsigned i=1;i<3;++i) if (std::abs(coordinates[i])>std::abs(coordinates[axis])) axis=i;
    const auto sign=coordinates[axis]>=0.?1.:-1.;
    auto q=p;std::array<double*,3> qAxes{&q.x,&q.y,&q.z};*qAxes[axis]=sign*.5001;
    check((p-q).norm()<.65*field.h,"survivor fixture is not near the old accepted proposal");
    auto mesh=cubeScaffold();
    Vec3 n{},u{},v{};std::array<double*,3> na{&n.x,&n.y,&n.z},ua{&u.x,&u.y,&u.z},va{&v.x,&v.y,&v.z};
    *na[axis]=sign;*ua[(axis+1)%3]=1.;*va[(axis+2)%3]=1.;
    const auto a=mesh.addNode(q,{0,201}),b=mesh.addNode(q+n*.02,{0,202}),
        c=mesh.addNode(q+n*.01+u*.02,{0,203}),d=mesh.addNode(q+n*.01+v*.02,{0,204});
    addClosedTet(mesh,{a,b,c,d},201);
    check(mesh::auditLinearVolumeBoundary(mesh).valid,"separate live-vertex component must be embedded");
    const auto before=mesh;bool reachedNative=false;
    const auto out=advanceOcctInteriorFront(mesh,field,limits,
        [&](Vec3 candidate){reachedNative|=(candidate-p).norm()<1e-10;return true;},
        [](Vec3){return true;},[](std::uint64_t){});
    check(out.initialized&&out.spacingRejects>0&&!reachedNative,
        "a nearby live survivor outside the candidate cavity was omitted from spacing checks");
    fixedBoundary(before,mesh);allAddedCones(mesh,before.nodeCount());
    allVertexSeparation(mesh,before.nodeCount(),field);
}

void exerciseNativeAndLimits() {
    const ConstantField field{.45};const auto limits=testLimits();
    for (const bool refuseNode:{false,true}) {
        auto mesh=cubeScaffold();const auto before=mesh;std::uint64_t nodeCalls=0,centroidCalls=0;
        const auto out=advanceOcctInteriorFront(mesh,field,limits,
            [&](Vec3){++nodeCalls;return !refuseNode;},
            [&](Vec3){++centroidCalls;return false;},[](std::uint64_t){});
        check(out.initialized&&out.sites==0&&out.nativeRejects>0&&nodeCalls>0&&
            (refuseNode?centroidCalls==0:centroidCalls>0),"native refusal did not guard the complete prospective cavity");
        sameMesh(before,mesh);
    }
    for (unsigned mode=0;mode<7;++mode) {
        auto mesh=cubeScaffold();const auto before=mesh;auto bounded=limits;
        if (mode==0) bounded.nodes=std::uint64_t(mesh.nodeCount());
        if (mode==1) bounded.tets=std::uint64_t(mesh.countCells(mesh::CellType::Tet4));
        if (mode==2) bounded.queries=0;
        if (mode==3) bounded.work=0;
        if (mode==4) bounded.workingBytes=1;
        if (mode==5) bounded.maximumSites=0;
        if (mode==6) bounded.maximumCavityCells=0;
        std::uint64_t queries=0;
        const auto out=advanceOcctInteriorFront(mesh,field,bounded,
            [&](Vec3){++queries;return true;},[&](Vec3){++queries;return true;},[](std::uint64_t){});
        check(out.sites==0&&queries==0&&out.workUsed<=bounded.work&&out.nativeQueries<=bounded.queries,
            "front exceeded a node, tet, query, work, storage or site limit before refusal");
        sameMesh(before,mesh);
    }
    for (const bool nodeException:{false,true}) {
        auto mesh=cubeScaffold();const auto before=mesh;unsigned nodes=0,centroids=0;bool caught=false;
        struct NativeStop {};
        try {
            (void)advanceOcctInteriorFront(mesh,field,limits,
                [&](Vec3){if (++nodes==2&&nodeException) throw NativeStop{};return true;},
                [&](Vec3){if (++centroids==2&&!nodeException) throw NativeStop{};return true;},[](std::uint64_t){});
        } catch (const NativeStop&) {caught=true;}
        check(caught&&(nodeException?nodes==2:centroids==2),"external native exception was swallowed or never exercised");
        sameMesh(before,mesh);
    }
    {
        auto mesh=cubeScaffold();const auto before=mesh;unsigned nodes=0;bool caught=false;
        struct WorkStop {};
        try {
            (void)advanceOcctInteriorFront(mesh,field,limits,
                [&](Vec3){++nodes;return true;},[](Vec3){return true;},
                [&](std::uint64_t){if (nodes>=1) throw WorkStop{};});
        } catch (const WorkStop&) {caught=true;}
        check(caught&&nodes==1,"external charged-work interruption was not exercised after staging");
        sameMesh(before,mesh);
    }
    {
        auto mesh=cubeScaffold();const auto before=mesh;
        for (Index i=0;i<mesh.cellCount();++i) if (mesh.cell(i).type==mesh::CellType::Tet4) {
            // Rebuild one cell with another original-body identity, preserving
            // the exact geometry. Mixed material must not enter the front.
            mesh::Mesh mixed;
            for (Index n=0;n<mesh.nodeCount();++n) mixed.addNode(mesh.node(n),mesh.nodeClass(n));
            for (Index c=0;c<mesh.cellCount();++c) {
                const auto cell=mesh.cell(c);
                mixed.addCell(cell.type,std::span(mesh.cellNodes(c),mesh::nodesPerCell(cell.type)),
                    c==i?mesh::Classification{3,8}:cell.onEntity);
            }
            const auto original=mixed;
            const auto out=advanceOcctInteriorFront(mixed,field,limits,
                [](Vec3){return true;},[](Vec3){return true;},[](std::uint64_t){});
            check(!out.initialized&&out.sites==0,"mixed original-body scaffold was accepted");
            sameMesh(original,mixed);break;
        }
        sameMesh(before,mesh);
    }
}

std::string snapshot(const geo::OcctModel& model) {
    std::ostringstream stream;BRepTools::Write(OcctModelAccess::shape(model),stream);return stream.str();
}
void exerciseNativeFieldBounds() {
    auto imported=OcctModelAccess::fromShape(BRepPrimAPI_MakeCylinder(1.,2.).Shape());
    check(imported.imported(),"field-bound cylinder import failed");
    int wall=-1;
    for (const auto* face:imported.model->faces()) {
        const auto* surface=imported.model->surfaceDescriptor(face->tag());
        if (surface&&surface->kind==geo::SurfaceKind::Cylinder) wall=face->tag();
    }
    check(wall>0,"field-bound fixture needs its curved wall");
    geo::NativeCadMeshOptions options;options.geometry.targetSize=.7;options.minimumSize=.03;
    options.faceSizes.push_back({wall,.2,1.55});
    sizing::RegionSizeSpecification specification;specification.baseSize=.7;specification.minimumSize=.03;
    specification.controls.push_back({sizing::SizeRegionKind::Box,{.145,.065,.875},{.195,.115,.925},1.,.05,1.55});
    auto compiled=sizing::compileRegionSizeField(specification);check(compiled.compiled,"field-bound control compilation failed");
    const OcctSpatialSizing field(*imported.model,options,std::move(compiled.field));
    const std::array boxes{std::pair{Vec3{-.2,-.2,.6},Vec3{.3,.3,1.2}},
        std::pair{Vec3{.75,-.08,.5},Vec3{1.05,.08,.9}},
        std::pair{Vec3{2.,-.2,.4},Vec3{2.3,.2,.8}}};
    for (const auto& [lo,hi]:boxes) {
        const auto lower=field.minimumOnBox(lo,hi);
        check(lower>0.&&std::isfinite(lower),"native field returned an invalid conservative lower bound");
        for (unsigned i=0;i<3;++i) for (unsigned j=0;j<3;++j) for (unsigned k=0;k<3;++k) {
            const Vec3 p{std::lerp(lo.x,hi.x,i*.5),std::lerp(lo.y,hi.y,j*.5),std::lerp(lo.z,hi.z,k*.5)};
            check(lower<=field.eval(p)+1e-12,"native face/box lower bound exceeds a sampled field value");
        }
    }
    check(field.minimumOnBox(boxes[0].first,boxes[0].second)<=.05+1e-12,
        "native conservative field query missed a tiny enclosed box");
}

void exercisePublicNativeFront() {
    auto imported=OcctModelAccess::fromShape(BRepPrimAPI_MakeCylinder(1.,2.).Shape());
    check(imported.imported(),"production front cylinder import failed");
    const auto& model=*imported.model;const auto originalCad=snapshot(model);
    int wall=-1;
    for (const auto* face:model.faces()) {
        const auto* surface=model.surfaceDescriptor(face->tag());
        if (surface&&surface->kind==geo::SurfaceKind::Cylinder) wall=face->tag();
    }
    check(wall>0,"production front fixture needs its original wall tag");
    geo::NativeCadMeshOptions options;
    options.geometry.targetSize=.7;options.geometry.threads=1;
    options.geometry.maximumNodes=50'000;options.geometry.maximumCells=250'000;
    options.geometry.maximumOutputBytes=128ULL*1024*1024;options.geometry.maximumMeshingWork=100'000'000;
    options.geometry.maximumAuditTemporaryBytes=512ULL*1024*1024;
    options.geometry.minimumElementQuality=0.;options.geometry.warningElementQuality=.25;
    options.maximumValidationWork=100'000'000;options.maximumPrimitiveWorkingBytes=512ULL*1024*1024;
    options.maximumNativeValidationThreads=1;options.maximumChordError=.0086602540378443865;
    options.maximumRelativeVolumeError=.01;options.minimumSize=.12;
    options.faceSizes.push_back({wall,.16,1.55});
    options.sizeBoxes.push_back({{-.3,-.3,.7},{.3,.3,1.3},.15,1.55});
    options.useSpatialSizing=true;options.improveNativeTetQuality=true;
    options.allowCentroidsWithinChordTolerance=false;options.recoveryMode=geo::NativeCadRecoveryMode::Standard;
    mesh::Mesh baseline;
    const auto ordinary=geo::meshNativeCad(model,options,baseline);
    check(ordinary.published&&ordinary.audit.valid&&ordinary.quality.valid&&ordinary.bodies.size()==1,
        "front-disabled production baseline failed");
    check(std::string(ordinary.bodies.front().interiorFrontDiagnostic)=="not_attempted",
        "disabled front changed the established route");
    options.useAdvancingInteriorFront=true;
    mesh::Mesh advanced;const auto built=geo::meshNativeCad(model,options,advanced);
    if (!built.published) std::cerr<<"native front publication: "<<geo::nativeCadMeshReasonName(built.reason)<<'\n';
    check(built.published&&built.audit.complete&&built.audit.valid&&built.quality.complete&&built.quality.valid&&built.bodies.size()==1,
        "public advancing-front route must publish an audited native mesh");
    const auto& body=built.bodies.front();
    const auto tets=std::uint64_t(advanced.countCells(mesh::CellType::Tet4));
    check(body.solidTag==1&&body.boundaryEmbeddingCertified&&body.nativeCentroidsChecked==tets&&
        body.nativeCentroidThreads==1&&body.nativeCentroidSourceGuard&&body.centroidOutsideDistanceLimit==0.&&
        body.maximumOutsideCentroidDistance==0.&&body.centroidsWithinChordTolerance==0.,
        "advancing front lost strict original-body native centroid coverage");
    check(body.interiorFrontNodes>0&&body.interiorFrontCandidates>=body.interiorFrontNodes&&
        body.interiorFrontActivations>0&&body.interiorFrontAdvances>0&&body.interiorFrontWork>0&&
        body.interiorFrontWork<=options.geometry.maximumMeshingWork&&body.interiorFrontQueries>0,
        "production route did not exercise admitted advancing-front cavities");
    const std::string diagnostic=body.interiorFrontDiagnostic;
    check(diagnostic!="not_attempted"&&diagnostic!="scaffold_refused"&&
        (diagnostic!="targets_satisfied"||body.interiorFrontRemainingOversized==0),
        "production front result hides refusal or unresolved size targets");
    check(boundarySignature(baseline)==boundarySignature(advanced),
        "advancing-front opt-in changed exact prescribed boundary coordinates, winding or tags");
    for (Index i=0;i<advanced.cellCount();++i) {
        const auto cell=advanced.cell(i);
        if (cell.type==mesh::CellType::Tet4) check(cell.onEntity==mesh::Classification{3,1},"front lost original body tags");
        else check(cell.type==mesh::CellType::Tri3&&cell.onEntity.dim==2&&model.face(cell.onEntity.tag),
            "front lost original CAD boundary face tags");
    }
    check(snapshot(model)==originalCad,"advancing-front production route modified original CAD");
    std::cout<<"native front cylinder: sites="<<body.interiorFrontNodes<<" advances="<<body.interiorFrontAdvances
        <<" voids="<<body.interiorVoidActivations<<" diagnostic="<<diagnostic<<'\n';
}
}

int main(int argc,char** argv) {
    try {
        const bool pure=argc==2&&std::string_view(argv[1])=="--pure";
        check(argc==1||pure,"usage: interior_front [--pure]");
        exerciseLogarithmicMetric();
        exerciseUniformAndEnclosed();
        exerciseFrontAndVoidFairness();
        exerciseSmoothRegionalCoverage();
        exerciseWholeCavityShape();
        exerciseAllLiveVertexSpacing();
        exerciseNativeAndLimits();
        exerciseNativeFieldBounds();
        if (!pure) exercisePublicNativeFront();
        std::cout<<"native interior front geometry, enclosed fields, live spacing and transactional admissions passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n';return 1;
    } catch (...) {
        std::cerr<<"unexpected interior-front regression failure\n";return 2;
    }
}
