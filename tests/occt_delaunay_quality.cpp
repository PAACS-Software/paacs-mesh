// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/geometry/OcctDelaunayMeshQuality.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/mesh/MeshAudit.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace paacs;
namespace {
using Complex=bnd::detail::BistellarComplex;
using Vertex=bnd::TetFacetVertexId;
constexpr std::array<unsigned char,5> fixed{1,1,1,1,0};
void require(bool ok,const char* message) {
    if (!ok) throw std::runtime_error(message);
}
Complex closedStar() {
    // A thin, asymmetric material star: its volume centroid remains inside,
    // but moving all the way there makes the worst radius ratio smaller.
    Complex live;
    live.walkConstraints=true;
    for (const auto p:std::array{Vec3{0,0,0},Vec3{1,0,0},Vec3{0,1,0},
                                Vec3{0,0,.2},Vec3{.25,.22,.048}})
        live.vertices.push_back({p,live.vertices.size()+1});
    for (Vertex omitted=0;omitted<4;++omitted) {
        bnd::TetFacetCell cell;
        cell.vertices[0]=4;
        unsigned next=1;
        for (Vertex vertex=0;vertex<4;++vertex)
            if (vertex!=omitted) cell.vertices[next++]=vertex;
        cell.classification={3,17};
        require(live.positive(cell),"fixture must contain positive tetrahedra");
        live.insert(cell);
        const auto face=bnd::tet_facet_detail::outwardFace(cell,
            unsigned(std::ranges::find(cell.vertices,Vertex{4})-cell.vertices.begin()));
        live.protectedFaces.insert(Complex::key(face));
        for (unsigned i=0;i<3;++i)
            live.protectedEdges.insert(Complex::key(Complex::Edge{face[i],face[(i+1)%3]}));
    }
    // Freeze topology through the supported storage bound. The test measures
    // geometric improvement of this existing fill, with its CAD shell fixed.
    live.maximumStoredCells=live.cells.size();
    return live;
}
mesh::Mesh asMesh(const Complex& live) {
    mesh::Mesh out;
    for (std::size_t i=0;i<live.vertices.size();++i)
        out.addNode(live.vertices[i].position,i<4?mesh::Classification{0,int(i)+1}:
                                                   mesh::Classification{3,17});
    int faceTag=1;
    for (const auto& face:live.protectedFaces)
        out.addCell(mesh::CellType::Tri3,{Index(face[0]),Index(face[1]),Index(face[2])},
                    {2,faceTag++});
    for (std::size_t id=0;id<live.cells.size();++id) if (live.alive[id]) {
        const auto& cell=live.cells[id];
        out.addCell(mesh::CellType::Tet4,
            {Index(cell.vertices[0]),Index(cell.vertices[1]),
             Index(cell.vertices[2]),Index(cell.vertices[3])},cell.classification);
    }
    return out;
}
double minimumRadiusRatio(const Complex& live,Vec3 interior) {
    double result=1.;
    for (std::size_t id=0;id<live.cells.size();++id) if (live.alive[id]) {
        std::array<Vec3,4> p{};
        for (unsigned j=0;j<4;++j) {
            const auto v=live.cells[id].vertices[j];
            p[j]=v==4?interior:live.vertices[v].position;
        }
        result=std::min(result,mesh::tetQuality(p[0],p[1],p[2],p[3]));
    }
    return result;
}
std::array<unsigned,2> poorCounts(const Complex& live) {
    std::array<unsigned,2> result{};
    for (std::size_t id=0;id<live.cells.size();++id) if (live.alive[id]) {
        mesh::Mesh single;
        for (const auto v:live.cells[id].vertices)
            single.addNode(live.vertices[v].position,{3,17});
        single.addCell(mesh::CellType::Tet4,{0,1,2,3},{3,17});
        const auto quality=algo::analyzeMixedMeshQuality(single);
        require(quality.complete&&quality.valid,"each final tet is exact positive");
        result[0]+=quality.minimumQuality<.1;
        result[1]+=quality.minimumQuality<.25;
    }
    return result;
}
std::uint64_t signature(const Complex& live) {
    return mesh::meshHash(asMesh(live));
}
void requireUnchanged(const Complex& live,std::uint64_t before,unsigned edits) {
    require(edits==0&&signature(live)==before,"refused move must preserve the mesh exactly");
}
struct NativeStop {};
struct WorkStop {};
}
int main() {
    try {
        auto live=closedStar();
        const auto before=algo::analyzeMixedMeshQuality(asMesh(live));
        const auto beforePoor=poorCounts(live);
        const auto initial=live.vertices[4].position;
        const Vec3 fullCentroid{.25,.25,.05};
        const auto beforeRadius=minimumRadiusRatio(live,initial);
        require(minimumRadiusRatio(live,fullCentroid)<beforeRadius-1e-4,
                "fixture reproduces a full centroid step that degrades radius ratio");
        const auto boundary=live.protectedFaces;
        const auto edges=live.protectedEdges;
        std::array<Vec3,4> boundaryPositions{};
        for (unsigned i=0;i<4;++i) boundaryPositions[i]=live.vertices[i].position;
        std::uint64_t work=0, geometricWork=0;
        unsigned edits=0, nodeQueries=0, centroidQueries=0;
        live.chargeWork=[&](std::uint64_t n){work+=n;};
        geo::detail::optimizeOcctDelaunayQuality(live,fixed,
            [&](Vec3){++centroidQueries;return true;},
            [&](Vec3){++nodeQueries;geometricWork=work;return true;},edits,4);
        const auto after=algo::analyzeMixedMeshQuality(asMesh(live));
        require(edits==1&&after.complete&&after.valid,
                "bounded relocation must improve the positive closed star");
        require(after.minimumQuality>before.minimumQuality+1e-3&&
                minimumRadiusRatio(live,live.vertices[4].position)>=beforeRadius,
                "partial relocation improves reported quality without worsening radius ratio");
        const auto afterPoor=poorCounts(live);
        require(afterPoor[0]<=beforePoor[0]&&afterPoor[1]<=beforePoor[1],
                "relocation cannot increase reported poor-cell counts");
        require(std::abs(after.totalVolume-before.totalVolume)<=before.totalVolume*1e-12,
                "relocation preserves material volume");
        require(nodeQueries==1&&centroidQueries==4,
                "only the selected geometric trial requires native admissions");
        require(live.cells.size()==4&&live.protectedFaces==boundary&&live.protectedEdges==edges,
                "relocation retains the boundary and connectivity");
        for (unsigned i=0;i<4;++i) {
            const auto p=live.vertices[i].position, original=boundaryPositions[i];
            require(p.x==original.x&&p.y==original.y&&p.z==original.z,
                    "CAD boundary vertices stay bitwise fixed");
        }
        for (const auto& cell:live.cells) {
            const auto& v=cell.vertices;
            require(cell.classification==mesh::Classification{3,17}&&
                    predicates::orient3d(live.vertices[v[0]].position,live.vertices[v[1]].position,
                        live.vertices[v[2]].position,live.vertices[v[3]].position)==predicates::Sign::Positive,
                    "body classification and exact positive orientation are preserved");
        }
        const auto audit=mesh::auditMesh(asMesh(live));
        require(audit.passed()&&audit.boundaryMatched==4&&audit.boundaryOrphans==0&&
                audit.boundaryUnexpected==0&&audit.nonManifoldFaces==0,
                "improved star passes the independent conformity audit");

        for (bool rejectNode:{false,true}) {
            auto refused=closedStar();
            const auto original=signature(refused);
            unsigned refusedEdits=0;
            geo::detail::optimizeOcctDelaunayQuality(refused,fixed,
                [&](Vec3){return rejectNode;},[&](Vec3){return !rejectNode;},refusedEdits,4);
            requireUnchanged(refused,original,refusedEdits);
        }
        for (bool throwAtNode:{false,true}) {
            auto refused=closedStar();
            const auto original=signature(refused);
            unsigned refusedEdits=0, queries=0;
            bool caught=false;
            try {
                geo::detail::optimizeOcctDelaunayQuality(refused,fixed,
                    [&](Vec3){if (++queries==2) throw NativeStop{};return true;},
                    [&](Vec3){if (throwAtNode) throw NativeStop{};return true;},refusedEdits,4);
            } catch (const NativeStop&) {caught=true;}
            require(caught,"native interruption must propagate to the caller");
            requireUnchanged(refused,original,refusedEdits);
        }
        require(geometricWork>0,"successful trial records bounded geometric work");
        for (bool afterNative:{false,true}) {
            auto refused=closedStar();
            const auto original=signature(refused);
            unsigned refusedEdits=0, queries=0;
            std::uint64_t spent=0;
            bool allAdmitted=false, caught=false;
            refused.chargeWork=[&](std::uint64_t n) {
                if ((afterNative&&allAdmitted)||(!afterNative&&n>geometricWork-1-spent))
                    throw WorkStop{};
                spent+=n;
            };
            try {
                geo::detail::optimizeOcctDelaunayQuality(refused,fixed,
                    [&](Vec3){allAdmitted=++queries==4;return true;},
                    [](Vec3){return true;},refusedEdits,4);
            } catch (const WorkStop&) {caught=true;}
            require(caught,"work interruption must propagate before coordinate publication");
            requireUnchanged(refused,original,refusedEdits);
        }
        auto frozen=closedStar();
        const auto original=signature(frozen);
        const std::array<unsigned char,5> allFixed{1,1,1,1,1};
        unsigned frozenEdits=0;
        geo::detail::optimizeOcctDelaunayQuality(frozen,allFixed,
            [](Vec3){return true;},[](Vec3){return true;},frozenEdits,4);
        requireUnchanged(frozen,original,frozenEdits);

        auto transaction=asMesh(closedStar());
        const auto originalMesh=mesh::meshHash(transaction);
        for (const auto queries:std::array<std::uint64_t,2>{0,1}) {
            const auto refused=geo::detail::improveOcctDelaunayMeshQuality(transaction,
                200000,queries,4,64*1024*1024,[](Vec3){return true;},
                [](Vec3){return true;},[](std::uint64_t){});
            require(!refused&&mesh::meshHash(transaction)==originalMesh,
                    "postfill query exhaustion cannot publish a partial transaction");
        }
        const auto accepted=geo::detail::improveOcctDelaunayMeshQuality(transaction,
            200000,100,4,64*1024*1024,[](Vec3){return true;},
            [](Vec3){return true;},[](std::uint64_t){});
        require(accepted&&algo::analyzeMixedMeshQuality(transaction).minimumQuality>
                    before.minimumQuality+1e-3&&mesh::auditMesh(transaction).passed(),
                "independent postfill transaction publishes an audited improvement");
        for (Index i=0;i<4;++i) {
            const auto p=transaction.node(i), expected=boundaryPositions[std::size_t(i)];
            require(p.x==expected.x&&p.y==expected.y&&p.z==expected.z&&
                    transaction.nodeClass(i)==mesh::Classification{0,int(i)+1},
                    "postfill transaction preserves boundary positions and tags");
        }
        std::cout<<"Partial closed-star quality improvement preserves exact, native, volume, and boundary gates.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
