// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/geometry/OcctDelaunayMeshQuality.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/mesh/MeshAudit.hpp"

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
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

Complex capBipyramid(double epsilon=std::ldexp(1.,-45),double scale=1.) {
    Complex live;
    live.walkConstraints=true;
    for (const auto p:std::array{Vec3{0,0,0},Vec3{scale,0,0},Vec3{0,scale,0},
                                Vec3{.3*scale,.3*scale,epsilon},Vec3{.3*scale,.3*scale,-scale}})
        live.vertices.push_back({p,live.vertices.size()+1});
    for (const auto vertices:std::array<std::array<Vertex,4>,2>{{{0,1,2,3},{0,2,1,4}}}) {
        bnd::TetFacetCell cell{vertices,0,{3,17}};
        require(live.positive(cell),"cap bipyramid requires positive original cells");
        live.insert(cell);
    }
    for (const auto& [face,owners]:live.faces) if (owners.size()==1) {
        live.protectedFaces.insert(face);
        for (unsigned i=0;i<3;++i)
            live.protectedEdges.insert(Complex::key(Complex::Edge{face[i],face[(i+1)%3]}));
    }
    return live;
}

mesh::Mesh asCapMesh(const Complex& live) {
    mesh::Mesh result;
    for (std::size_t i=0;i<live.vertices.size();++i)
        result.addNode(live.vertices[i].position,{0,int(i)+1});
    int faceTag=1;
    for (const auto& key:live.protectedFaces) {
        const auto owner=*live.faces.at(key).begin();
        const auto& cell=live.cells[owner];
        for (unsigned side=0;side<4;++side) {
            auto face=bnd::tet_facet_detail::outwardFace(cell,side);
            if (Complex::key(face)!=key) continue;
            // Fix a cyclic starting vertex without changing orientation so
            // equivalent oriented patch shells have identical test payloads.
            std::rotate(face.begin(),std::min_element(face.begin(),face.end()),face.end());
            result.addCell(mesh::CellType::Tri3,{face[0],face[1],face[2]},{2,faceTag++});
            break;
        }
    }
    for (std::size_t id=0;id<live.cells.size();++id) if (live.alive[id]) {
        const auto& cell=live.cells[id];
        const auto& v=cell.vertices;
        result.addCell(mesh::CellType::Tet4,{v[0],v[1],v[2],v[3]},cell.classification);
    }
    return result;
}

void sameCapBoundary(const mesh::Mesh& before,const mesh::Mesh& after) {
    require(before.nodeCount()==after.nodeCount() &&
            before.countCells(mesh::CellType::Tri3)==after.countCells(mesh::CellType::Tri3),
            "cap face repair must preserve all boundary nodes and six facets");
    for (Index i=0;i<before.nodeCount();++i) {
        const auto a=before.node(i),b=after.node(i);
        require(std::bit_cast<std::uint64_t>(a.x)==std::bit_cast<std::uint64_t>(b.x) &&
                std::bit_cast<std::uint64_t>(a.y)==std::bit_cast<std::uint64_t>(b.y) &&
                std::bit_cast<std::uint64_t>(a.z)==std::bit_cast<std::uint64_t>(b.z) &&
                before.nodeClass(i)==after.nodeClass(i),
                "cap repair moved or reclassified an original boundary vertex");
    }
    for (Index i=0;i<before.countCells(mesh::CellType::Tri3);++i) {
        const auto a=before.cell(i),b=after.cell(i);
        require(a.type==b.type && a.onEntity==b.onEntity &&
                std::equal(before.cellNodes(i),before.cellNodes(i)+3,after.cellNodes(i)),
                "cap repair changed the exact oriented CAD facet payload");
    }
}

void exerciseCriticalFaceRepair() {
    constexpr std::array<unsigned char,5> allFixed{1,1,1,1,1};
    auto live=capBipyramid();
    const auto source=asCapMesh(live);
    const auto before=algo::analyzeMixedMeshQuality(source);
    const auto beforePoor=poorCounts(live);
    const auto beforeRadius=minimumRadiusRatio(live,live.vertices[4].position);
    require(before.valid && before.minimumQuality<.01 &&
            mesh::auditLinearVolumeBoundary(source).valid,
            "cap fixture must reproduce an exactly positive sliver inside a closed oriented shell");
    const auto boundary=live.protectedFaces;
    const auto edges=live.protectedEdges;
    unsigned edits=0,nodeQueries=0,centroidQueries=0;
    std::uint64_t work=0,firstNativeWork=0;
    live.chargeWork=[&](std::uint64_t n) {work+=n;};
    geo::detail::optimizeOcctDelaunayQuality(live,allFixed,
        [&](Vec3) {if (!centroidQueries) firstNativeWork=work;++centroidQueries;return true;},
        [&](Vec3) {++nodeQueries;return false;},edits,3);
    const auto repaired=asCapMesh(live);
    const auto after=algo::analyzeMixedMeshQuality(repaired);
    require(edits==1 && live.metrics.flips23==1 && after.valid && after.tet4==3 &&
            after.minimumQuality>.25 && minimumRadiusRatio(live,live.vertices[4].position)>beforeRadius,
            "fixed-vertex cap sliver requires an improving two-to-three face repair");
    require(nodeQueries==0 && centroidQueries==3 && live.edges.contains(Complex::Edge{3,4}) &&
            !live.faces.contains(Complex::Face{0,1,2}),
            "cap repair must replace only its unprotected shared face with the interior apex edge");
    require(live.protectedFaces==boundary && live.protectedEdges==edges &&
            mesh::auditLinearVolumeBoundary(repaired).valid &&
            std::abs(after.totalVolume-before.totalVolume)<=before.totalVolume*1e-12,
            "critical face repair must retain the protected shell and material volume");
    sameCapBoundary(source,repaired);
    const auto afterPoor=poorCounts(live);
    require(afterPoor[0]==0 && afterPoor[1]==0 && beforePoor[0]>0,
            "critical face repair must remove the sliver without increasing poor-cell counts");
    for (std::size_t id=0;id<live.cells.size();++id) if (live.alive[id]) {
        const auto& cell=live.cells[id];
        const auto& v=cell.vertices;
        require(cell.classification==mesh::Classification{3,17} &&
                predicates::orient3d(live.vertices[v[0]].position,live.vertices[v[1]].position,
                    live.vertices[v[2]].position,live.vertices[v[3]].position)==predicates::Sign::Positive,
                "cap replacement must preserve original body tags and exact positive orientation");
    }

    // Each refusal leaves the complete fixed-vertex complex untouched. The
    // boundary edges have stars of size two, so no edge move can hide a face
    // repair that crossed its protected/interface or resource boundary.
    for (unsigned refusal=0;refusal<5;++refusal) {
        auto refused=capBipyramid();
        if (refusal==0) refused.protectedFaces.insert({0,1,2});
        if (refusal==1) refused.cells[1].classification={3,18};
        if (refusal==2) refused.maximumStoredCells=refused.cells.size()+2;
        const auto original=mesh::meshHash(asCapMesh(refused));
        unsigned refusedEdits=0,queries=0;
        geo::detail::optimizeOcctDelaunayQuality(refused,allFixed,
            [&](Vec3) {++queries;return refusal!=4;},
            [](Vec3) {return true;},refusedEdits,refusal==3?2:3);
        require(refusedEdits==0 && mesh::meshHash(asCapMesh(refused))==original,
                "cap face repair violated a protected face, solid interface, storage/tet cap or native refusal");
        if (refusal==4) require(queries>0,"cap repair must obtain original-solid centroid admission");
    }
    {
        auto refused=capBipyramid();
        const auto original=mesh::meshHash(asCapMesh(refused));
        unsigned refusedEdits=0,queries=0;
        bool caught=false;
        try {
            geo::detail::optimizeOcctDelaunayQuality(refused,allFixed,
                [&](Vec3) {if (++queries==2) throw NativeStop{};return true;},
                [](Vec3) {return true;},refusedEdits,3);
        } catch (const NativeStop&) {caught=true;}
        require(caught && refusedEdits==0 && mesh::meshHash(asCapMesh(refused))==original,
                "interrupted cap centroid admission must leave the face transaction uncommitted");
    }
    require(firstNativeWork>0,"critical face repair must charge its geometric proposal work");
    for (bool afterNative:{false,true}) {
        auto refused=capBipyramid();
        const auto original=mesh::meshHash(asCapMesh(refused));
        unsigned refusedEdits=0,queries=0;
        std::uint64_t spent=0;
        bool caught=false;
        refused.chargeWork=[&](std::uint64_t n) {
            if ((afterNative && queries==3) || (!afterNative && n>firstNativeWork-1-spent))
                throw WorkStop{};
            spent+=n;
        };
        try {
            geo::detail::optimizeOcctDelaunayQuality(refused,allFixed,
                [&](Vec3) {++queries;return true;},[](Vec3) {return true;},refusedEdits,3);
        } catch (const WorkStop&) {caught=true;}
        require(caught && refusedEdits==0 && mesh::meshHash(asCapMesh(refused))==original,
                "cap face repair must stop before publication when charged work is exhausted");
    }
    for (const auto queries:std::array<std::uint64_t,3>{0,1,2}) {
        auto transaction=source;
        const auto original=mesh::meshHash(transaction);
        const auto refused=geo::detail::improveOcctDelaunayMeshQuality(transaction,
            200000,queries,3,64*1024*1024,[](Vec3) {return true;},
            [](Vec3) {return true;},[](std::uint64_t) {});
        require(!refused && mesh::meshHash(transaction)==original,
                "postfill cap query exhaustion must retain the full original mesh");
    }
    auto transaction=source;
    const auto accepted=geo::detail::improveOcctDelaunayMeshQuality(transaction,
        200000,3,3,64*1024*1024,[](Vec3) {return false;},
        [](Vec3) {return true;},[](std::uint64_t) {});
    require(accepted==1 && algo::analyzeMixedMeshQuality(transaction).minimumQuality>.25 &&
            mesh::auditLinearVolumeBoundary(transaction).valid,
            "postfill transaction must publish the fully admitted fixed-boundary cap repair");
    sameCapBoundary(source,transaction);
}

void exerciseSevereTopologyOnly() {
    constexpr geo::detail::OcctDelaunayQualityOptions options{1e-6,64,true};
    constexpr std::array<unsigned char,5> allFixed{1,1,1,1,1};
    auto live=capBipyramid(1e-9,2.);
    const auto source=asCapMesh(live);
    const auto faces=live.protectedFaces;
    const auto edges=live.protectedEdges;
    const auto before=algo::analyzeMixedMeshQuality(source);
    const auto poor=poorCounts(live);
    require(before.complete&&before.valid&&before.minimumQuality<options.criticalBelowMeanRatio&&
            mesh::auditLinearVolumeBoundary(source).valid,
            "owned severe cap must start below the numeric-cleanup threshold inside a valid shell");
    unsigned edits=0,nodeQueries=0,centroidQueries=0;
    geo::detail::optimizeOcctDelaunayQuality(live,allFixed,
        [&](Vec3){++centroidQueries;return true;},[&](Vec3){++nodeQueries;return false;},edits,3,options);
    const auto result=asCapMesh(live);
    const auto after=algo::analyzeMixedMeshQuality(result);
    require(edits==1&&edits<=options.maxEdits&&nodeQueries==0&&centroidQueries==3&&
            after.complete&&after.valid&&after.minimumQuality>.25,
            "severe-only mode must repair the fixed cap through admitted topology without node moves");
    require(live.protectedFaces==faces&&live.protectedEdges==edges&&
            mesh::auditLinearVolumeBoundary(result).valid&&
            std::abs(after.totalVolume-before.totalVolume)<=before.totalVolume*1e-12,
            "severe-only cleanup changed the protected material domain");
    sameCapBoundary(source,result);
    const auto nextPoor=poorCounts(live);
    require(nextPoor[0]<=poor[0]&&nextPoor[1]<=poor[1],
            "severe-only cleanup increased the low-quality population");

    // This movable interior star is improved by the ordinary default policy.
    // The numeric-only policy must return before that relocation/search stage.
    auto movable=closedStar();const auto original=signature(movable);
    unsigned skipped=0,moveQueries=0;
    geo::detail::optimizeOcctDelaunayQuality(movable,fixed,[](Vec3){return true;},
        [&](Vec3){++moveQueries;return true;},skipped,4,options);
    requireUnchanged(movable,original,skipped);
    require(moveQueries==0,"severe-only policy queried an ordinary smoothing proposal");

    const auto small=asCapMesh(capBipyramid(1e-9,2.));
    for(unsigned refusal=0;refusal<4;++refusal){
        auto transaction=small;const auto signatureBefore=mesh::meshHash(transaction);
        auto bounded=options;if(refusal==3)bounded.maxEdits=0;
        const auto changed=geo::detail::improveOcctDelaunayMeshQuality(transaction,
            refusal==0?1:200000,refusal==1?2:3,3,64*1024*1024,
            [](Vec3){throw std::runtime_error("numeric cleanup attempted a node move");return false;},
            [&](Vec3){return refusal!=2;},[](std::uint64_t){},bounded);
        require(changed==0&&mesh::meshHash(transaction)==signatureBefore,
                "numeric-only resource, Native or edit-cap refusal published a partial transaction");
    }
    auto transaction=small;
    const auto changed=geo::detail::improveOcctDelaunayMeshQuality(transaction,
        200000,3,3,64*1024*1024,[](Vec3){return false;},[](Vec3){return true;},
        [](std::uint64_t){},options);
    require(changed==1&&mesh::auditLinearVolumeBoundary(transaction).valid&&
            algo::analyzeMixedMeshQuality(transaction).minimumQuality>.25,
            "numeric-only mesh wrapper did not publish the complete admitted topology repair");
    sameCapBoundary(small,transaction);

    auto defaultPolicy=closedStar(),explicitDefault=closedStar();
    unsigned first=0,second=0;
    geo::detail::optimizeOcctDelaunayQuality(defaultPolicy,fixed,[](Vec3){return true;},
        [](Vec3){return true;},first,4);
    geo::detail::optimizeOcctDelaunayQuality(explicitDefault,fixed,[](Vec3){return true;},
        [](Vec3){return true;},second,4,{});
    require(first==second&&signature(defaultPolicy)==signature(explicitDefault),
            "explicit default options changed the established optimization policy");
}
}
int main() {
    try {
        exerciseCriticalFaceRepair();
        exerciseSevereTopologyOnly();
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
