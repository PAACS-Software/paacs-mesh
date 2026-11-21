// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "../bnd/TetBistellarRecovery.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <span>
#include <tuple>
#include <vector>

namespace paacs::geo::detail {
namespace occt_recone_detail {
using Complex=bnd::detail::BistellarComplex;
using Vertex=bnd::TetFacetVertexId;
inline constexpr std::size_t maximumStarCells=64;

inline bool finite(Vec3 p) {
    return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);
}

inline std::optional<std::vector<std::size_t>> edgeStar(
    const Complex& live,Complex::Edge edge,std::size_t maximumCells) {
    edge=Complex::key(edge);
    if (!live.walkConstraints || !maximumCells || edge[0]==edge[1] ||
        edge[1]>=live.vertices.size() || live.protectedEdges.contains(edge)) return {};
    maximumCells=std::min(maximumCells,maximumStarCells);
    const auto owners=live.edges.find(edge);
    if (owners==live.edges.end()) return {};
    std::vector<std::size_t> result;
    std::optional<mesh::Classification> region;
    for (const auto id:owners->second) {
        live.charge();
        if (id>=live.alive.size()||id>=live.cells.size()) return {};
        if (!live.alive[id]) continue;
        if (result.size()>=maximumCells) return {};
        const auto& cell=live.cells[id];
        if (!region) region=cell.classification;
        if (cell.classification!=*region || region->dim!=3 ||
            std::ranges::find(cell.vertices,edge[0])==cell.vertices.end() ||
            std::ranges::find(cell.vertices,edge[1])==cell.vertices.end()) return {};
        result.push_back(id);
    }
    if (result.empty()) return {};
    const std::set<std::size_t> selected(result.begin(),result.end());
    std::set<std::size_t> visited{result.front()};
    std::vector<std::size_t> connected{result.front()};
    for (std::size_t i=0;i<connected.size();++i) {
        const auto id=connected[i];
        if (id>=live.neighbours.size()) return {};
        for (unsigned side=0;side<4;++side) {
            live.charge();
            const auto next=live.neighbours[id][side];
            if (!selected.contains(next)) continue;
            if (live.protectedFaces.contains(Complex::key(
                bnd::tet_facet_detail::outwardFace(live.cells[id],side)))) return {};
            if (visited.insert(next).second) connected.push_back(next);
        }
    }
    if (visited.size()!=selected.size()) return {};
    return result;
}
} // namespace occt_recone_detail

// Candidate apices come from one face-adjacent layer outside the complete
// edge star. Prefer movable interior vertices, then distance and canonical ID.
// This only ranks existing vertices; it changes neither geometry nor topology.
inline std::vector<bnd::TetFacetVertexId> occtDelaunayEdgeReconeApices(
    const bnd::detail::BistellarComplex& live,
    bnd::detail::BistellarComplex::Edge edge,
    std::span<const unsigned char> fixed,std::size_t maximumApices=8) {
    using namespace occt_recone_detail;
    if (!maximumApices) return {};
    const auto star=edgeStar(live,edge,maximumStarCells);
    if (!star) return {};
    maximumApices=std::min<std::size_t>(maximumApices,8);
    const std::set<std::size_t> initial(star->begin(),star->end());
    std::set<Vertex> initialVertices,candidates;
    const auto region=live.cells[star->front()].classification;
    for (const auto id:*star) for (const auto vertex:live.cells[id].vertices) {
        live.charge();
        if (vertex>=live.vertices.size()) return {};
        initialVertices.insert(vertex);
    }
    for (const auto id:*star) {
        if (id>=live.neighbours.size()) return {};
        for (unsigned side=0;side<4;++side) {
            live.charge();
            const auto next=live.neighbours[id][side];
            const auto face=bnd::tet_facet_detail::outwardFace(live.cells[id],side);
            if (next==Complex::noCell || initial.contains(next) ||
                live.protectedFaces.contains(Complex::key(face))) continue;
            if (next>=live.alive.size() || next>=live.cells.size()) return {};
            if (!live.alive[next] || live.cells[next].classification!=region) continue;
            for (const auto vertex:live.cells[next].vertices) {
                live.charge();
                if (vertex>=live.vertices.size()) return {};
                if (!initialVertices.contains(vertex)&&finite(live.vertices[vertex].position))
                    candidates.insert(vertex);
            }
        }
    }
    edge=Complex::key(edge);
    const auto a=live.vertices[edge[0]].position,b=live.vertices[edge[1]].position;
    if (!finite(a)||!finite(b)) return {};
    const Vec3 midpoint{std::midpoint(a.x,b.x),std::midpoint(a.y,b.y),std::midpoint(a.z,b.z)};
    std::vector<std::tuple<bool,double,Vertex>> ranked;
    for (const auto vertex:candidates) {
        live.charge();
        const auto p=live.vertices[vertex].position;
        const auto distance=std::hypot(p.x-midpoint.x,p.y-midpoint.y,p.z-midpoint.z);
        if (std::isfinite(distance))
            ranked.emplace_back(vertex<fixed.size()&&fixed[vertex]!=0,distance,vertex);
    }
    std::ranges::sort(ranked,[&](const auto& a,const auto& b) {live.charge();return a<b;});
    std::vector<Vertex> result;
    for (const auto& item:ranked) {
        if (result.size()>=maximumApices) break;
        live.charge();
        result.push_back(std::get<2>(item));
    }
    return result;
}

// Remove the entire unprotected edge star, then grow only through obstructing
// shell faces until an existing apex sees the remaining shell. Faces already
// containing the apex need no cone. Exact orientation and the old oriented
// shell are retained; protected faces and original-body interfaces block growth.
// This is a proposal only. The caller still certifies quality, material/native
// admissions, live.valid() and preservesSurvivors() before applying it.
inline std::optional<bnd::detail::BistellarComplex::Proposal> buildOcctDelaunayEdgeRecone(
    const bnd::detail::BistellarComplex& live,
    bnd::detail::BistellarComplex::Edge edge,bnd::TetFacetVertexId apex,
    std::size_t maximumCells=64) {
    using namespace occt_recone_detail;
    using predicates::Sign;
    using bnd::tet_facet_detail::outwardFace;
    edge=Complex::key(edge);
    if (apex>=live.vertices.size() || apex==edge[0] || apex==edge[1] ||
        !finite(live.vertices[apex].position)) return {};
    const auto star=edgeStar(live,edge,maximumCells);
    if (!star) return {};
    maximumCells=std::min(maximumCells,maximumStarCells);
    Complex::Proposal patch;
    patch.kind=Complex::nmKind;
    patch.removed=*star;
    std::set<std::size_t> removed(star->begin(),star->end());
    const auto region=live.cells[star->front()].classification;
    const auto position=[&](Vertex vertex) -> const Vec3& {return live.vertices[vertex].position;};
    for (std::size_t i=0;i<patch.removed.size();++i) {
        const auto id=patch.removed[i];
        if (id>=live.neighbours.size()) return {};
        const auto& cell=live.cells[id];
        std::set<Vertex> vertices;
        for (const auto vertex:cell.vertices) {
            live.charge();
            if (vertex>=live.vertices.size()||!finite(position(vertex))) return {};
            vertices.insert(vertex);
        }
        const auto& v=cell.vertices;
        live.charge();
        if (vertices.size()!=4 || predicates::orient3dFiltered(
            position(v[0]),position(v[1]),position(v[2]),position(v[3]))!=Sign::Positive) return {};
        for (unsigned side=0;side<4;++side) {
            live.charge();
            const auto next=live.neighbours[id][side];
            const auto face=outwardFace(cell,side);
            if (removed.contains(next) || std::ranges::find(face,apex)!=face.end()) continue;
            live.charge();
            if (predicates::orient3dFiltered(position(face[0]),position(face[1]),
                position(face[2]),position(apex))==Sign::Negative) continue;
            if (patch.removed.size()>=maximumCells || next==Complex::noCell ||
                next>=live.alive.size() || next>=live.cells.size() || !live.alive[next] ||
                live.cells[next].classification!=region ||
                live.protectedFaces.contains(Complex::key(face))) return {};
            removed.insert(next);
            patch.removed.push_back(next);
        }
    }

    std::vector<bnd::TetFacetCell> old;
    std::set<Vertex> oldVertices,newVertices;
    for (const auto id:patch.removed) {
        live.charge();
        const auto& cell=live.cells[id];
        old.push_back(cell);
        for (const auto vertex:cell.vertices) oldVertices.insert(vertex);
        for (unsigned side=0;side<4;++side) {
            live.charge();
            const auto face=outwardFace(cell,side);
            if (removed.contains(live.neighbours[id][side])) continue;
            if (std::ranges::find(face,apex)!=face.end()) continue;
            const bnd::TetFacetCell cone{{apex,face[0],face[1],face[2]},0,region};
            const auto& v=cone.vertices;
            live.charge();
            if (predicates::orient3dFiltered(position(v[0]),position(v[1]),position(v[2]),
                position(v[3]))!=Sign::Positive) return {};
            if (std::ranges::find(v,edge[0])!=v.end()&&std::ranges::find(v,edge[1])!=v.end()) return {};
            for (const auto vertex:v) newVertices.insert(vertex);
            patch.added.push_back(cone);
        }
    }
    if (patch.added.empty() || !oldVertices.contains(apex) || oldVertices!=newVertices ||
        live.cells.size()>live.maximumStoredCells ||
        patch.added.size()>live.maximumStoredCells-live.cells.size()) return {};
    live.charge(16*(old.size()+patch.added.size()));
    const bnd::TetFacetComplexView view{live.vertices,old,{}, {},{},region};
    const auto oldFaces=bnd::tet_facet_detail::candidateFaceMap(view,old);
    const auto newFaces=bnd::tet_facet_detail::candidateFaceMap(view,patch.added);
    const auto oriented=[](const auto& faces) {
        return std::ranges::all_of(faces,[](const auto& item) {
            const auto& uses=item.second;
            return uses.size()==1||(uses.size()==2&&uses[0].sense!=uses[1].sense);
        });
    };
    if (!oriented(oldFaces)||!oriented(newFaces)||
        bnd::tet_facet_detail::boundarySenses(oldFaces)!=bnd::tet_facet_detail::boundarySenses(newFaces)) return {};
    for (const auto& [face,uses]:oldFaces) {
        (void)uses;
        live.charge();
        if (live.protectedFaces.contains(face.vertices)&&!newFaces.contains(face)) return {};
    }
    std::set<Complex::Edge> newEdges;
    for (const auto& cell:patch.added) for (const auto e:Complex::edgeKeys(cell)) {
        live.charge();
        newEdges.insert(e);
    }
    for (const auto& cell:old) for (const auto e:Complex::edgeKeys(cell)) {
        live.charge();
        if (live.protectedEdges.contains(e)&&!newEdges.contains(e)) return {};
    }
    return patch;
}

} // namespace paacs::geo::detail
