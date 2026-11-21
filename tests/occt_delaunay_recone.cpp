// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/geometry/OcctDelaunayRecone.hpp"
#include "../src/geometry/OcctDelaunayQuality.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/mesh/Quality.hpp"
#include <array>
#include <bit>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

using namespace paacs;

namespace {
using Complex=bnd::detail::BistellarComplex;
using Vertex=bnd::TetFacetVertexId;
using geo::detail::buildOcctDelaunayEdgeRecone;
using geo::detail::occtDelaunayEdgeReconeApices;
constexpr mesh::Classification region{3,1};
constexpr Complex::Edge target{1,6};
constexpr Vertex apex=8;
constexpr std::array<unsigned char,9> fixed{1,1,1,1,1,1,1,0,0};
constexpr std::array<mesh::Classification,9> classifications{{
    {1,1},{1,1},{2,1},{2,1},{2,2},{2,2},{2,2},{3,1},{3,1}
}};

void require(bool value,const char* message) {
    if (!value) throw std::runtime_error(message);
}

bool samePoint(Vec3 a,Vec3 b) {
    return std::bit_cast<std::uint64_t>(a.x)==std::bit_cast<std::uint64_t>(b.x)&&
        std::bit_cast<std::uint64_t>(a.y)==std::bit_cast<std::uint64_t>(b.y)&&
        std::bit_cast<std::uint64_t>(a.z)==std::bit_cast<std::uint64_t>(b.z);
}

Complex capCavity() {
    // Nine cells extracted from this project's generated, rigidly transformed
    // cylinder. Canonical source IDs are 16,17,94,403,432,473,480,762,763.
    // The almost-flat cap tet is cell 6; neither diagonal of its four-cell
    // edge star is a valid refill. The existing interior apex needs five more
    // face-adjacent cells before it sees the full unchanged oriented shell.
    const std::array<Vec3,9> points{{
        {4.23193161000574,-3.231931610005458,4.211880524929385},
        {4.128446996391514,-3.3496787384737075,4.2048188472023895},
        {4.234079813076767,-3.1231218958066425,4.098793598361706},
        {4.111951377780863,-3.2317757582966564,4.084291400371732},
        {4.269853127273194,-3.3699108730557903,4.124788266989986},
        {4.167843699656054,-3.465852896378875,4.127792369296103},
        {4.350391838191693,-3.252442777926622,4.1430727990392935},
        {4.092239052086066,-3.2713831409052188,3.930529638354639},
        {4.185652860532872,-3.1734056628882756,3.9548270292027814}
    }};
    constexpr std::array<std::uint64_t,9> ids{16,17,94,403,432,473,480,762,763};
    constexpr std::array<std::array<Vertex,4>,9> cells{{
        {5,3,4,7},{3,8,4,7},{5,3,1,4},{8,3,6,2},{8,4,6,3},
        {6,1,3,4},{6,1,4,0},{2,1,3,6},{2,1,6,0}
    }};
    Complex live;
    live.walkConstraints=true;
    for (std::size_t i=0;i<points.size();++i) live.vertices.push_back({points[i],ids[i]+1});
    for (const auto& vertices:cells) {
        require(predicates::orient3d(points[vertices[0]],points[vertices[1]],
            points[vertices[2]],points[vertices[3]])==predicates::Sign::Positive,
            "owned cap-cavity source must be exactly positive without reorientation");
        live.insert({vertices,0,region});
    }
    for (const auto& [face,owners]:live.faces) if (owners.size()==1) {
        live.protectedFaces.insert(face);
        for (unsigned i=0;i<3;++i)
            live.protectedEdges.insert(Complex::key(Complex::Edge{face[i],face[(i+1)%3]}));
    }
    return live;
}

std::map<Complex::Face,int> prescribedTags() {
    return {{{0,1,2},1},{{0,1,4},2},{{0,4,6},2},{{1,2,3},1},{{1,4,5},2}};
}

mesh::Mesh asMesh(const Complex& live) {
    mesh::Mesh result;
    for (std::size_t i=0;i<live.vertices.size();++i)
        result.addNode(live.vertices[i].position,i<classifications.size()?classifications[i]:region);
    const auto tags=prescribedTags();
    for (std::size_t id=0;id<live.cells.size();++id) if (live.alive[id]) {
        const auto& cell=live.cells[id];
        for (unsigned side=0;side<4;++side) {
            const auto face=bnd::tet_facet_detail::outwardFace(cell,side);
            const auto key=Complex::key(face);
            if (!live.protectedFaces.contains(key)) continue;
            const auto tag=tags.find(key);
            result.addCell(mesh::CellType::Tri3,{face[0],face[1],face[2]},
                {2,tag==tags.end()?90:tag->second});
        }
        result.addCell(mesh::CellType::Tet4,{cell.vertices[0],cell.vertices[1],
            cell.vertices[2],cell.vertices[3]},cell.classification);
    }
    return result;
}

void unchanged(const Complex& before,const Complex& after) {
    require(before.vertices.size()==after.vertices.size()&&before.cells.size()==after.cells.size()&&
        before.alive==after.alive&&before.faces==after.faces&&before.edges==after.edges&&
        before.neighbours==after.neighbours&&before.vertexCell==after.vertexCell&&
        before.protectedFaces==after.protectedFaces&&before.protectedEdges==after.protectedEdges,
        "proposal planning changed the live complex or an index");
    for (std::size_t i=0;i<before.vertices.size();++i)
        require(samePoint(before.vertices[i].position,after.vertices[i].position)&&
            before.vertices[i].sequence==after.vertices[i].sequence,
            "proposal planning changed a canonical point or identity");
    for (std::size_t i=0;i<before.cells.size();++i)
        require(before.cells[i].vertices==after.cells[i].vertices&&
            before.cells[i].sequence==after.cells[i].sequence&&
            before.cells[i].classification==after.cells[i].classification,
            "proposal planning changed source connectivity or body identity");
}

double minimumMeanRatio(const Complex& live,std::span<const bnd::TetFacetCell> cells) {
    double result=1.;
    for (const auto& cell:cells) {
        const auto& v=cell.vertices;
        const auto determinant=predicates::orient3dDeterminantEstimate(
            live.vertices[v[0]].position,live.vertices[v[1]].position,
            live.vertices[v[2]].position,live.vertices[v[3]].position);
        double sum=0.;
        for (const auto edge:Complex::edgeKeys(cell))
            sum+=(live.vertices[edge[0]].position-live.vertices[edge[1]].position).norm2();
        require(determinant>0.&&sum>0.,"quality witness must retain exact positive cells");
        result=std::min(result,12.*std::pow(determinant/2.,2./3.)/sum);
    }
    return result;
}

void exerciseRecone() {
    auto live=capCavity();
    require(live.edges.at(target)==std::set<std::size_t>{5,6,7,8},
        "fixture must contain the complete nonconvex four-cell edge star");
    require(live.edgeFlips(target,std::nullopt).empty()&&
        !live.flip23(Complex::Face{1,4,6},std::nullopt)&&
        !live.flip23(Complex::Face{0,1,6},std::nullopt),
        "fixture must defeat both star diagonals and the critical cell's local face swaps");
    const auto before=live;
    const auto source=asMesh(live);
    require(mesh::auditLinearVolumeBoundary(source).valid,"owned source cavity needs an embedded oriented shell");
    std::uint64_t work=0;
    live.chargeWork=[&](std::uint64_t n) {work+=n;};
    const auto apices=occtDelaunayEdgeReconeApices(live,target,fixed);
    require(apices==std::vector<Vertex>{8,5},
        "one-hop candidates must prefer the existing nonfixed interior apex");
    require(occtDelaunayEdgeReconeApices(live,target,fixed,1)==std::vector<Vertex>{8}&&
        occtDelaunayEdgeReconeApices(live,target,fixed,0).empty(),
        "apex ranking must respect its explicit candidate limit");
    auto patch=buildOcctDelaunayEdgeRecone(live,target,apices.front(),9);
    require(patch&&patch->removed.size()==9&&patch->added.size()==9,
        "existing-apex recone must close the full expanded nine-cell cavity");
    require(std::set<std::size_t>(patch->removed.begin(),patch->removed.end())==
        std::set<std::size_t>{0,1,2,3,4,5,6,7,8},"recone missed a required closure neighbour");
    require(work>0,"discovery and exact visibility checks must consume work");
    unchanged(before,live);
    live.chargeWork={};
    require(live.valid(*patch,std::nullopt)&&live.preservesSurvivors(*patch),
        "expanded recone must pass canonical shell and exact contact certification");
    require(minimumMeanRatio(live,live.cells)<1e-8&&minimumMeanRatio(live,patch->added)>.18,
        "expanded existing-apex refill must remove the near-zero-quality cap tetrahedron");
    for (const auto& cell:patch->added) {
        const auto edges=Complex::edgeKeys(cell);
        require(cell.classification==region&&std::ranges::find(edges,target)==edges.end(),
            "expanded recone must retain body identity and remove the whole target edge");
    }
    const bnd::TetFacetComplexView view{live.vertices,{}, {},{}, {},region};
    const auto oldShell=bnd::tet_facet_detail::boundarySenses(
        bnd::tet_facet_detail::candidateFaceMap(view,live.cells));
    const auto newShell=bnd::tet_facet_detail::boundarySenses(
        bnd::tet_facet_detail::candidateFaceMap(view,patch->added));
    require(oldShell==newShell,"existing-apex proposal changed the canonical oriented shell");
    for (const auto& [face,tag]:prescribedTags()) {
        require(tag>0&&newShell.contains(bnd::tet_facet_detail::faceKey(face)),
            "existing-apex proposal lost an original classified CAD face");
    }
    live.apply(*patch);
    require(!live.edges.contains(target)&&live.vertices.size()==before.vertices.size()&&
        live.protectedFaces==before.protectedFaces&&live.protectedEdges==before.protectedEdges,
        "committed recone changed the point set or fixed constraints");
    for (std::size_t i=0;i<before.vertices.size();++i)
        require(samePoint(before.vertices[i].position,live.vertices[i].position)&&
            before.vertices[i].sequence==live.vertices[i].sequence,
            "committed recone moved or replaced an existing canonical vertex");
    const auto refined=asMesh(live);
    require(mesh::auditLinearVolumeBoundary(refined).valid,"expanded recone changed the material boundary");
    const auto oldVolume=mesh::totalTetVolume(source),newVolume=mesh::totalTetVolume(refined);
    require(std::abs(oldVolume-newVolume)<=oldVolume*1e-12,"expanded recone changed material volume");
}

void exerciseRefusals() {
    for (const auto limit:{std::size_t(0),std::size_t(3),std::size_t(8)}) {
        auto live=capCavity();const auto before=live;
        require(!buildOcctDelaunayEdgeRecone(live,target,apex,limit),
            "recone exceeded the complete-star or expansion cell limit");
        unchanged(before,live);
    }
    {
        auto live=capCavity();live.maximumStoredCells=live.cells.size()+8;const auto before=live;
        require(!buildOcctDelaunayEdgeRecone(live,target,apex),"recone exceeded append storage");
        unchanged(before,live);
    }
    for (const auto face:{Complex::Face{1,3,4},Complex::Face{1,4,6}}) {
        auto live=capCavity();live.protectedFaces.insert(face);const auto before=live;
        require(!buildOcctDelaunayEdgeRecone(live,target,apex),
            "recone crossed or swallowed a protected face");
        unchanged(before,live);
    }
    {
        auto live=capCavity();live.protectedEdges.insert(target);const auto before=live;
        require(!buildOcctDelaunayEdgeRecone(live,target,apex)&&
            occtDelaunayEdgeReconeApices(live,target,fixed).empty(),"recone removed a protected edge");
        unchanged(before,live);
    }
    for (const auto cell:{std::size_t(2),std::size_t(5)}) {
        auto live=capCavity();live.cells[cell].classification={3,2};const auto before=live;
        require(!buildOcctDelaunayEdgeRecone(live,target,apex),"recone crossed a body interface");
        unchanged(before,live);
    }
    {
        auto live=capCavity();const auto before=live;
        require(!buildOcctDelaunayEdgeRecone(live,target,bnd::kInvalidTetFacetVertex)&&
            !buildOcctDelaunayEdgeRecone(live,target,target[0]),"recone accepted an invalid apex");
        unchanged(before,live);
    }
    {
        auto live=capCavity();live.vertices.push_back({live.vertices[apex].position,999});const auto before=live;
        require(!buildOcctDelaunayEdgeRecone(live,target,9),
            "recone used an apex absent from the removed existing vertex set");
        unchanged(before,live);
    }
    for (const bool discover:{false,true}) {
        auto live=capCavity();const auto before=live;
        struct WorkStop {};
        std::uint64_t work=0;bool caught=false;
        live.chargeWork=[&](std::uint64_t n) {if ((work+=n)>30) throw WorkStop{};};
        try {
            if (discover) (void)occtDelaunayEdgeReconeApices(live,target,fixed);
            else (void)buildOcctDelaunayEdgeRecone(live,target,apex);
        } catch (const WorkStop&) {caught=true;}
        require(caught,"recone work interruption was ignored");
        unchanged(before,live);
    }
}

std::vector<bnd::TetFacetCell> activeCells(const Complex& live) {
    std::vector<bnd::TetFacetCell> result;
    for (std::size_t id=0;id<live.cells.size();++id)
        if (live.alive[id]) result.push_back(live.cells[id]);
    return result;
}

void exerciseOptimizerIntegration() {
    using geo::detail::optimizeOcctDelaunayQuality;
    // Freezing every point isolates topology: repairing this cap cannot rely
    // on relocation, jitter, a new point or a changed prescribed boundary.
    constexpr std::array<unsigned char,9> allFixed{1,1,1,1,1,1,1,1,1};
    auto live=capCavity();
    const auto before=live;
    const auto source=asMesh(before);
    unsigned edits=0,nodeQueries=0,centroidQueries=0;
    std::uint64_t work=0,firstNativeWork=0;
    live.chargeWork=[&](std::uint64_t n) {work+=n;};
    optimizeOcctDelaunayQuality(live,allFixed,
        [&](Vec3) {
            if (!centroidQueries) firstNativeWork=work;
            ++centroidQueries;
            return true;
        },
        [&](Vec3) {++nodeQueries;return false;},edits,9);
    const auto cells=activeCells(live);
    require(edits>=1&&centroidQueries>0&&nodeQueries==0&&firstNativeWork>0,
        "critical optimizer must admit an existing-apex topology repair without any node query");
    require(!live.edges.contains(target)&&!cells.empty()&&cells.size()<=9&&
        minimumMeanRatio(live,cells)>.01,
        "critical optimizer must remove the edge and all critical cells within the live-tet limit");
    require(live.vertices.size()==before.vertices.size()&&
        live.protectedFaces==before.protectedFaces&&live.protectedEdges==before.protectedEdges,
        "optimizer changed the fixed point set or prescribed constraints");
    for (std::size_t i=0;i<before.vertices.size();++i)
        require(samePoint(before.vertices[i].position,live.vertices[i].position)&&
            before.vertices[i].sequence==live.vertices[i].sequence,
            "optimizer moved or replaced an all-fixed canonical vertex");
    for (const auto& cell:cells) {
        const auto& v=cell.vertices;
        require(cell.classification==region&&predicates::orient3d(
            live.vertices[v[0]].position,live.vertices[v[1]].position,
            live.vertices[v[2]].position,live.vertices[v[3]].position)==predicates::Sign::Positive,
            "critical optimizer changed body identity or exact positive orientation");
    }
    const bnd::TetFacetComplexView view{live.vertices,{}, {},{}, {},region};
    require(bnd::tet_facet_detail::boundarySenses(
        bnd::tet_facet_detail::candidateFaceMap(view,before.cells))==
        bnd::tet_facet_detail::boundarySenses(
            bnd::tet_facet_detail::candidateFaceMap(view,cells)),
        "critical optimizer changed the original oriented cavity shell");
    const auto refined=asMesh(live);
    require(mesh::auditLinearVolumeBoundary(refined).valid,
        "critical optimizer broke the original fixed material boundary");
    const auto oldVolume=mesh::totalTetVolume(source),newVolume=mesh::totalTetVolume(refined);
    require(std::abs(oldVolume-newVolume)<=oldVolume*1e-12,
        "critical optimizer changed the original material volume");
    for (Index i=0;i<source.nodeCount();++i)
        require(source.nodeClass(i)==refined.nodeClass(i),"critical optimizer changed a CAD node classification");

    {
        auto refused=capCavity();const auto original=refused;
        unsigned refusedEdits=0,refusedNodes=0,refusedCentroids=0;
        optimizeOcctDelaunayQuality(refused,allFixed,
            [&](Vec3) {++refusedCentroids;return false;},
            [&](Vec3) {++refusedNodes;return false;},refusedEdits,9);
        require(refusedEdits==0&&refusedCentroids>0&&refusedNodes==0,
            "native centroid refusal must block the critical optimizer's complete edit");
        unchanged(original,refused);
    }
    {
        auto refused=capCavity();const auto original=refused;
        struct WorkStop {};
        unsigned refusedEdits=0,queries=0;
        std::uint64_t used=0;
        const auto limit=firstNativeWork-1;
        refused.chargeWork=[&](std::uint64_t n) {
            if (n>limit-used) throw WorkStop{};
            used+=n;
        };
        bool caught=false;
        try {
            optimizeOcctDelaunayQuality(refused,allFixed,
                [&](Vec3) {++queries;return true;},
                [&](Vec3) {++queries;return false;},refusedEdits,9);
        } catch (const WorkStop&) {caught=true;}
        require(caught&&queries==0&&refusedEdits==0,
            "pre-admission work exhaustion must leave the critical edit uncommitted");
        unchanged(original,refused);
    }
    for (const auto stopAt:{1u,2u}) {
        auto refused=capCavity();const auto original=refused;
        struct NativeStop {};
        unsigned refusedEdits=0,queries=0,nodeCalls=0;
        bool caught=false;
        try {
            optimizeOcctDelaunayQuality(refused,allFixed,
                [&](Vec3) {if (++queries==stopAt) throw NativeStop{};return true;},
                [&](Vec3) {++nodeCalls;return false;},refusedEdits,9);
        } catch (const NativeStop&) {caught=true;}
        require(caught&&queries==stopAt&&nodeCalls==0&&refusedEdits==0,
            "interrupted native admission must roll back the complete critical topology edit");
        unchanged(original,refused);
    }
}
}

int main() {
    try {
        exerciseRecone();
        exerciseRefusals();
        exerciseOptimizerIntegration();
        std::cout<<"native existing-apex edge recone, expanded shell and refusal regressions passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    } catch (...) {
        std::cerr<<"unexpected native existing-apex recone regression failure\n";
        return 2;
    }
}
