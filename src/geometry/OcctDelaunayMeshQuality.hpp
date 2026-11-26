// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "OcctDelaunayQuality.hpp"
#include "paacs/mesh/Mesh.hpp"
#include <new>

namespace paacs::geo::detail {

struct OcctDelaunayQualityLimit {};

// A separate, optional postfill transaction. Its budget does not reduce the
// insertion budget or change the prescribed surface. Only admitted edits are
// copied back; exhaustion keeps the last fully admitted complex.
template<class AdmitNode, class AdmitCentroid, class Charge>
unsigned improveOcctDelaunayMeshQuality(mesh::Mesh& mesh, std::uint64_t maximumWork,
    std::uint64_t maximumQueries, std::uint64_t maximumTets,
    std::uint64_t maximumWorkingBytes, AdmitNode&& admitNode,
    AdmitCentroid&& admitCentroid, Charge&& charge,
    const OcctDelaunayQualityOptions& options={}) try {
    using Complex=bnd::detail::BistellarComplex;
    using Vertex=bnd::TetFacetVertexId;
    const auto tets=std::uint64_t(mesh.countCells(mesh::CellType::Tet4));
    const auto nodeBytes=std::uint64_t(mesh.nodeCount())*64;
    if (nodeBytes>=maximumWorkingBytes) return 0;
    const auto stored=std::min(tets+8192,(maximumWorkingBytes-nodeBytes)/4096);
    if (!tets || tets>maximumTets || stored<=tets || !maximumQueries || !maximumWork ||
        std::uint64_t(mesh.nodeCount())>=bnd::kInvalidTetFacetVertex) return 0;
    Complex live;
    live.walkConstraints=true;
    live.maximumStoredCells=std::size_t(stored);
    std::uint64_t work=0,queries=0;
    live.chargeWork=[&](std::uint64_t n) {
        if (n>maximumWork-work) throw OcctDelaunayQualityLimit{};
        charge(n);
        work+=n;
    };
    const auto admit=[&](auto& callback,Vec3 p) {
        if (queries>=maximumQueries) throw OcctDelaunayQualityLimit{};
        charge(0);
        ++queries;
        const auto accepted=callback(p);
        charge(0);
        return accepted;
    };
    std::vector<unsigned char> fixed(std::size_t(mesh.nodeCount()));
    unsigned edits=0;
    try {
        live.vertices.reserve(std::size_t(mesh.nodeCount()));
        for (Index i=0;i<mesh.nodeCount();++i) {
            live.charge();
            live.vertices.push_back({mesh.node(i),std::uint64_t(i)+1});
            fixed[std::size_t(i)]=mesh.nodeClass(i).dim!=3;
        }
        for (Index i=0;i<mesh.cellCount();++i) {
            live.charge();
            const auto& cell=mesh.cell(i);
            const auto* v=mesh.cellNodes(i);
            if (cell.type==mesh::CellType::Tet4) {
                live.insert({{Vertex(v[0]),Vertex(v[1]),Vertex(v[2]),Vertex(v[3])},0,cell.onEntity});
            } else if (cell.type==mesh::CellType::Tri3) {
                const Complex::Face face{Vertex(v[0]),Vertex(v[1]),Vertex(v[2])};
                live.protectedFaces.insert(Complex::key(face));
                for (unsigned j=0;j<3;++j) {
                    fixed[std::size_t(face[j])]=1;
                    live.protectedEdges.insert(Complex::key(Complex::Edge{face[j],face[(j+1)%3]}));
                }
            } else return 0;
        }
        optimizeOcctDelaunayQuality(live,fixed,
            [&](Vec3 p){return admit(admitCentroid,p);},
            [&](Vec3 p){return admit(admitNode,p);},edits,maximumTets,options);
    } catch (const OcctDelaunayQualityLimit&) {
        // Initialization has no edits; optimization stops between transactions.
    }
    if (!edits) return 0;
    mesh::Mesh result;
    for (Index i=0;i<mesh.nodeCount();++i)
        result.addNode(live.vertices[std::size_t(i)].position,mesh.nodeClass(i));
    for (Index i=0;i<mesh.cellCount();++i) if (mesh.cell(i).type==mesh::CellType::Tri3)
        result.addCell(mesh::CellType::Tri3,std::span(mesh.cellNodes(i),3),mesh.cell(i).onEntity);
    for (std::size_t i=0;i<live.cells.size();++i) if (live.alive[i]) {
        const auto& cell=live.cells[i];
        const auto& v=cell.vertices;
        result.addCell(mesh::CellType::Tet4,{Index(v[0]),Index(v[1]),Index(v[2]),Index(v[3])},cell.classification);
    }
    try {charge(0);} catch (const OcctDelaunayQualityLimit&) {}
    mesh.swap(result);
    return edits;
} catch (const std::bad_alloc&) {
    // The input is untouched until swap, including allocation during export.
    return 0;
}

} // namespace paacs::geo::detail
