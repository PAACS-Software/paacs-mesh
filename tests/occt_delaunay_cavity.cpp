// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/geometry/OcctLocalDelaunay.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include <BRepPrimAPI_MakeBox.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace paacs;
using namespace paacs::geo::detail;

namespace {
using Complex = bnd::detail::BistellarComplex;
using Vertex = bnd::TetFacetVertexId;
constexpr mesh::Classification region{3,11};
constexpr double shapeFloor = 1e-6;

void check(bool value,const char* message) {
    if (!value) throw std::runtime_error(message);
}

Complex bipyramid(double scale=1.,bool bottom=true) {
    Complex live;
    live.walkConstraints=true;
    const std::array points{Vec3{0,0,0},Vec3{2,0,0},Vec3{0,2,0},
                            Vec3{0,0,2},Vec3{0,0,-2}};
    for (std::size_t i=0;i<(bottom?points.size():points.size()-1);++i)
        live.vertices.push_back({points[i]*scale,i+1});
    live.insert({{0,1,2,3},0,region});
    if (bottom) live.insert({{0,2,1,4},0,region});
    for (const auto& [face,owners]:live.faces) if (owners.size()==1) {
        live.protectedFaces.insert(face);
        for (unsigned i=0;i<3;++i)
            live.protectedEdges.insert(Complex::key(Complex::Edge{face[i],face[(i+1)%3]}));
    }
    return live;
}

mesh::Mesh asMesh(const Complex& live) {
    mesh::Mesh result;
    for (const auto& vertex:live.vertices) result.addNode(vertex.position,region);
    for (std::size_t id=0;id<live.cells.size();++id) if (live.alive[id]) {
        const auto& cell=live.cells[id];
        for (unsigned side=0;side<4;++side) {
            const auto face=bnd::tet_facet_detail::outwardFace(cell,side);
            if (live.protectedFaces.contains(Complex::key(face)))
                result.addCell(mesh::CellType::Tri3,{face[0],face[1],face[2]},{2,1});
        }
    }
    for (std::size_t id=0;id<live.cells.size();++id) if (live.alive[id]) {
        const auto& cell=live.cells[id];
        const auto& v=cell.vertices;
        result.addCell(mesh::CellType::Tet4,{v[0],v[1],v[2],v[3]},cell.classification);
    }
    return result;
}

void sameMesh(const mesh::Mesh& before,const mesh::Mesh& after) {
    check(before.nodeCount()==after.nodeCount() && before.cellCount()==after.cellCount(),
          "refused cavity changed mesh counts");
    for (Index i=0;i<before.nodeCount();++i) {
        const auto a=before.node(i), b=after.node(i);
        check(a.x==b.x && a.y==b.y && a.z==b.z && before.nodeClass(i)==after.nodeClass(i),
              "refused cavity changed an original vertex");
    }
    for (Index i=0;i<before.cellCount();++i) {
        const auto a=before.cell(i), b=after.cell(i);
        const auto count=a.type==mesh::CellType::Tri3?3:4;
        check(a.type==b.type && a.onEntity==b.onEntity &&
              std::equal(before.cellNodes(i),before.cellNodes(i)+count,after.cellNodes(i)),
              "refused cavity changed original connectivity or classification");
    }
}

double shape(const Complex& live,const bnd::TetFacetCell& cell) {
    const auto& v=cell.vertices;
    return mesh::tetScaledJacobian(live.vertices[v[0]].position,live.vertices[v[1]].position,
                                  live.vertices[v[2]].position,live.vertices[v[3]].position);
}

void exerciseExpansion(double scale) {
    auto live=bipyramid(scale);
    const auto source=asMesh(live);
    check(mesh::auditLinearVolumeBoundary(source).valid,"fixture needs an oriented closed shell");
    const Vec3 point{.25*scale,.25*scale,std::ldexp(scale,-52)};
    auto containing=Complex::noCell;

    auto original=buildOcctDelaunayCavity(live,0,point,2,containing,false);
    check(original && original->removed.size()==1 && original->added.size()==4,
          "default shape policy must retain the original local cavity");
    const auto staged=Vertex(live.vertices.size());
    live.vertices.push_back({point,std::uint64_t(staged)+1});
    check(live.valid(*original,std::nullopt,staged),"nearly flat fixture must be exactly admissible");
    check(std::ranges::any_of(original->added,[&](const auto& cell) {
        const auto& v=cell.vertices;
        return shape(live,cell)<shapeFloor && predicates::orient3d(
            live.vertices[v[0]].position,live.vertices[v[1]].position,
            live.vertices[v[2]].position,live.vertices[v[3]].position)==predicates::Sign::Positive;
    }),"default cavity must expose the positive near-flat cone regression");
    live.vertices.pop_back();

    std::uint64_t shapeWork=0;
    live.chargeWork=[&](std::uint64_t amount) { if (amount==4) shapeWork+=amount; };
    auto expanded=buildOcctDelaunayCavity(live,0,point,2,containing,false,shapeFloor);
    check(expanded && expanded->removed.size()==2 && expanded->added.size()==6 && containing==0,
          "poor cone must expand through the same-region interior face");
    check(shapeWork>0,"cone shape evaluation must consume work");
    live.chargeWork={};
    live.vertices.push_back({point,std::uint64_t(staged)+1});
    check(live.valid(*expanded,std::nullopt,staged) && live.preservesSurvivors(*expanded),
          "expanded cavity must preserve vertices, oriented shell and surviving cells");
    for (const auto& cell:expanded->added) {
        const auto& v=cell.vertices;
        check(shape(live,cell)>=shapeFloor && predicates::orient3d(
            live.vertices[v[0]].position,live.vertices[v[1]].position,
            live.vertices[v[2]].position,live.vertices[v[3]].position)==predicates::Sign::Positive,
              "expanded cavity retained a poor or nonpositive cone");
    }
    live.apply(*expanded);
    const auto refined=asMesh(live);
    check(mesh::auditLinearVolumeBoundary(refined).valid,"expanded replacement changed the shell");
    const auto before=mesh::totalTetVolume(source), after=mesh::totalTetVolume(refined);
    check(std::abs(after-before)<=1e-12*before,"expanded cavity changed total volume");
}

void exerciseRefusals() {
    const Vec3 point{.25,.25,std::ldexp(1.,-52)};
    auto live=bipyramid();
    const auto before=asMesh(live);
    auto containing=Complex::noCell;
    check(!buildOcctDelaunayCavity(live,0,point,1,containing,false,shapeFloor),
          "shape expansion must respect the cavity cell limit");
    sameMesh(before,asMesh(live));

    const Complex::Face shared{0,1,2};
    live.protectedFaces.insert(shared);
    check(!buildOcctDelaunayCavity(live,0,point,2,containing,false,shapeFloor),
          "shape expansion crossed a protected face");
    live.protectedFaces.erase(shared);
    sameMesh(before,asMesh(live));

    live.cells[1].classification={3,12};
    check(!buildOcctDelaunayCavity(live,0,point,2,containing,false,shapeFloor),
          "shape expansion crossed a region boundary");
    live.cells[1].classification=region;
    sameMesh(before,asMesh(live));

    struct WorkLimit {};
    bool stopped=false;
    live.chargeWork=[](std::uint64_t amount) { if (amount==4) throw WorkLimit{}; };
    try {
        (void)buildOcctDelaunayCavity(live,0,point,2,containing,false,shapeFloor);
    } catch (const WorkLimit&) { stopped=true; }
    live.chargeWork={};
    check(stopped,"shape evaluation must respect a work refusal");
    sameMesh(before,asMesh(live));
}

void exerciseSpatialAdmission() {
    auto imported=OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(2,2,2).Shape());
    check(imported.imported(),"native spatial field fixture import failed");
    sizing::RegionSizeSpecification specification;
    specification.baseSize=1.;
    specification.minimumSize=.05;
    sizing::RegionSizeControl control;
    control.lower={.2,.2,1e-16};
    control.upper={.3,.3,3e-16};
    control.targetSize=.1;
    specification.controls.push_back(control);
    auto compiled=sizing::compileRegionSizeField(specification);
    check(compiled.compiled,"spatial box field fixture failed");
    const OcctSpatialSizing field(*imported.model,{},std::move(compiled.field));
    LocalDelaunayLimits limits{100'000,64,32,64,1024*1024,1.,&field,true};
    const auto source=asMesh(bipyramid(1.,false));
    auto refused=source;
    const auto charge=[](std::uint64_t) {};
    const auto accept=[](Vec3) { return true; };
    const auto onlyNearBase=[](Vec3 p) { return p.z>0. && p.z<1e-12; };
    check(insertOcctDelaunaySites(refused,{},limits,onlyNearBase,accept,charge)==0,
          "spatial seed near a protected shell face must be refused");
    sameMesh(source,refused);

    auto legacy=source;
    auto legacyLimits=limits;
    legacyLimits.improveSpatialQuality=false;
    check(insertOcctDelaunaySites(legacy,{},legacyLimits,onlyNearBase,accept,charge)>0,
          "quality opt-out must preserve the existing spatial insertion policy");

    const std::array support{Vec3{.25,.25,2e-16}};
    auto supported=source;
    check(insertOcctDelaunaySites(supported,support,limits,accept,accept,charge)==1,
          "supplied support sites must keep their existing admission policy");
    check(mesh::auditLinearVolumeBoundary(supported).valid,"support insertion changed the shell");

    const std::array regularSupport{Vec3{.25,.25,.25}};
    unsigned centroidChecks=0;
    const auto refuseCentroid=[&](Vec3) { ++centroidChecks; return false; };
    check(insertOcctDelaunaySites(refused,regularSupport,limits,accept,refuseCentroid,charge)==0 &&
          centroidChecks>0,"native centroid refusal must leave the cavity uncommitted");
    sameMesh(source,refused);
    limits.queries=1;
    check(insertOcctDelaunaySites(refused,regularSupport,limits,accept,accept,charge)==0,
          "query exhaustion after staging must leave the cavity uncommitted");
    sameMesh(source,refused);
}
}

int main() {
    try {
        for (double scale:{1e-6,1.,1e6}) exerciseExpansion(scale);
        exerciseRefusals();
        exerciseSpatialAdmission();
        std::cout<<"native Delaunay cavity shape, protected shell and rollback regressions passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    } catch (...) {
        std::cerr<<"unexpected native Delaunay cavity regression failure\n";
        return 2;
    }
}
