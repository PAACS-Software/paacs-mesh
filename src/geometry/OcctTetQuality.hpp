// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "../bnd/TetBistellarRecovery.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/mesh/Quality.hpp"
#include <cmath>
#include <span>

namespace paacs::geo::detail {

struct OcctTetQualityLimits {
    std::uint64_t maximumWork;
    std::uint64_t maximumNodes;
    std::uint64_t maximumTets;
    double targetSize;
    std::uint64_t& remainingNativeQueries;
};

// Improve a small, already conforming piece. Its original exterior and cut
// triangles are frozen. Replacement cells are checked before the mesh changes.
template<class Admit, class Charge>
bool improveOcctTetPiece(mesh::Mesh& mesh, const OcctTetQualityLimits& limits,
                        Admit&& admit, Charge&& charge) {
    using Complex = bnd::detail::BistellarComplex;
    using Cell = bnd::TetFacetCell;
    using Vertex = bnd::TetFacetVertexId;
    using Sign = predicates::Sign;
    constexpr double searchBelow = .25;
    constexpr double poorQuality = .10;
    constexpr std::uint64_t maximumMoves = 96;
    struct TrialLimit {};
    const auto tetCount = std::uint64_t(mesh.countCells(mesh::CellType::Tet4));
    if (tetCount < 2 || tetCount > 4096 || tetCount > limits.maximumTets ||
        std::uint64_t(mesh.nodeCount()) > limits.maximumNodes) return false;

    std::uint64_t work = 0, moves = 0, live = tetCount, queries = 0;
    const auto nativeAdmits = [&](Vec3 p) {
        if (++queries>192 || !limits.remainingNativeQueries) throw TrialLimit{};
        --limits.remainingNativeQueries;
        return admit(p);
    };
    Complex complex;
    complex.walkConstraints = true;
    complex.maximumStoredCells = static_cast<std::size_t>(tetCount + maximumMoves*64);
    complex.chargeWork = [&](std::uint64_t amount) {
        if (amount > limits.maximumWork-work) throw TrialLimit{};
        charge(amount);
        work += amount;
    };
    const auto point = [&](Vertex vertex) -> const Vec3& { return complex.vertices[vertex].position; };
    const auto quality = [&](const Cell& cell) {
        const auto& v = cell.vertices;
        return mesh::tetScaledJacobian(point(v[0]),point(v[1]),point(v[2]),point(v[3]));
    };
    const auto oriented = [&](const auto& face, Vec3 p) {
        complex.charge();
        return predicates::orient3d(point(face[0]),point(face[1]),point(face[2]),p);
    };
    const auto accept = [&](const Complex::Proposal& patch, double oldQuality) {
        if (live-patch.removed.size()+patch.added.size() > limits.maximumTets) return false;
        // Adding cells below the reporting threshold can improve the local
        // minimum while making the aggregate bad-element count larger.
        const auto poorBefore=std::ranges::count_if(patch.removed,[&](std::size_t id){
            return quality(complex.cells[id])<poorQuality;
        });
        const auto poorAfter=std::ranges::count_if(patch.added,[&](const Cell& cell){
            return quality(cell)<poorQuality;
        });
        if (poorAfter>poorBefore) return false;
        for (const auto& cell : patch.added) if (!(quality(cell) >= oldQuality+1e-5)) return false;
        if (!complex.preservesSurvivors(patch)) return false;
        for (const auto& cell : patch.added) {
            Vec3 centroid{};
            for (const auto id : cell.vertices) centroid=centroid+point(id)*.25;
            complex.charge();
            if (!nativeAdmits(centroid)) return false;
        }
        complex.apply(patch);
        live=live-patch.removed.size()+patch.added.size();
        ++moves;
        return true;
    };
    try {
        complex.vertices.reserve(mesh.nodeCount()+16);
        for (const auto p : mesh.nodes()) {
            complex.charge();
            complex.vertices.push_back({p,complex.vertices.size()+1});
        }
        std::vector<std::pair<double,std::size_t>> queue;
        queue.reserve(tetCount);
        for (Index i=0; i<mesh.cellCount(); ++i) {
            complex.charge();
            const auto& item=mesh.cell(i);
            const auto* v=mesh.cellNodes(i);
            if (item.type==mesh::CellType::Tri3) {
                const Complex::Face face{Vertex(v[0]),Vertex(v[1]),Vertex(v[2])};
                complex.protectedFaces.insert(Complex::key(face));
                for (unsigned j=0;j<3;++j)
                    complex.protectedEdges.insert(Complex::key(Complex::Edge{face[j],face[(j+1)%3]}));
                continue;
            }
            if (item.type!=mesh::CellType::Tet4) return false;
            Cell tet{{Vertex(v[0]),Vertex(v[1]),Vertex(v[2]),Vertex(v[3])},0,item.onEntity};
            const auto score=quality(tet);
            if (score<searchBelow) queue.emplace_back(score,complex.cells.size());
            complex.insert(tet);
        }
        if (queue.empty()) return false;
        std::ranges::sort(queue);
        const auto improve = [&](std::size_t id) {
            const auto cell=complex.cells[id];
            const auto better = [&](const Complex::Proposal& patch) {
                double old=1.;
                for (const auto owner : patch.removed) old=std::min(old,quality(complex.cells[owner]));
                return accept(patch,old);
            };
            for (const auto face : Complex::faceKeys(cell)) {
                if (auto patch=complex.flip23(face,std::nullopt); patch && better(*patch)) return;
            }
            for (const auto edge : Complex::edgeKeys(cell)) {
                const auto found=complex.edges.find(edge);
                if (found==complex.edges.end() || found->second.size()>7) continue;
                for (const auto& patch : complex.edgeFlips(edge,std::nullopt))
                    if (better(patch)) return;
            }
        };
        for (unsigned pass=0;pass<2 && moves<maximumMoves;++pass) {
            const auto previousMoves=moves;
            for (const auto& [score,id] : queue) {
                complex.charge();
                if (complex.alive[id]) improve(id);
                if (moves==maximumMoves) break;
            }
            if (moves==previousMoves) break;
            queue.clear();
            for (std::size_t i=0;i<complex.cells.size();++i) if (complex.alive[i]) {
                complex.charge();
                const auto score=quality(complex.cells[i]);
                if (score<searchBelow) queue.emplace_back(score,i);
            }
            std::ranges::sort(queue);
        }
        // Try the circumcenter of each remaining oversized poor cell in its
        // local in-sphere cavity. A seed outside the material or a cavity
        // that cannot see its oriented shell leaves the original cells in use.
        unsigned trials=0;
        const auto initialNodes=complex.vertices.size();
        for (const auto& [score,id] : queue) {
            if (moves==maximumMoves || trials>=64 ||
                complex.vertices.size()-initialNodes>=16 ||
                complex.vertices.size()>=limits.maximumNodes) break;
            complex.charge();
            if (!complex.alive[id] || !(limits.targetSize>0)) continue;
            const auto& v=complex.cells[id].vertices;
            const auto a=point(v[1])-point(v[0]), b=point(v[2])-point(v[0]), c=point(v[3])-point(v[0]);
            const auto denominator=2*a.dot(b.cross(c));
            if (!(denominator>0)) continue;
            const auto offset=(b.cross(c)*a.norm2()+c.cross(a)*b.norm2()+a.cross(b)*c.norm2())/denominator;
            const auto radius=offset.norm();
            if (!(std::isfinite(radius) && radius>limits.targetSize*.5)) continue;
            ++trials;
            const auto seed=point(v[0])+offset;
            // Locate the seed by walking across any outward face that it
            // violates. Stop at the fixed exterior and at a local walk limit.
            auto containing=id;
            for (unsigned step=0;step<32;++step) {
                auto next=containing;
                for (unsigned side=0;side<4;++side) {
                    if (oriented(bnd::tet_facet_detail::outwardFace(complex.cells[containing],side),seed)!=Sign::Positive)
                        continue;
                    next=complex.neighbours[containing][side];
                    break;
                }
                if (next==containing) break;
                containing=next;
                if (containing==Complex::noCell || step==31) { containing=Complex::noCell; break; }
            }
            if (containing==Complex::noCell) continue;
            Complex::Proposal patch;
            patch.removed.push_back(containing);
            for (std::size_t i=0;i<patch.removed.size()&&patch.removed.size()<=16;++i) {
                for (const auto other : complex.neighbours[patch.removed[i]]) {
                    complex.charge();
                    if (other==Complex::noCell || std::ranges::find(patch.removed,other)!=patch.removed.end()) continue;
                    const auto& w=complex.cells[other].vertices;
                    if (predicates::insphere(point(w[0]),point(w[1]),point(w[2]),point(w[3]),seed)==Sign::Positive)
                        patch.removed.push_back(other);
                }
            }
            if (patch.removed.size()>16) continue;
            const auto staged=Vertex(complex.vertices.size());
            bool visible=true;
            double old=1.;
            for (const auto owner : patch.removed) {
                old=std::min(old,quality(complex.cells[owner]));
                for (unsigned side=0;side<4;++side) {
                    if (std::ranges::find(patch.removed,complex.neighbours[owner][side])!=patch.removed.end()) continue;
                    const auto face=bnd::tet_facet_detail::outwardFace(complex.cells[owner],side);
                    if (oriented(face,seed)!=Sign::Negative) { visible=false; break; }
                    patch.added.push_back({{staged,face[0],face[1],face[2]},0,complex.cells[owner].classification});
                }
                if (!visible) break;
            }
            if (!visible) continue;
            complex.vertices.push_back({seed,std::uint64_t(staged)+1});
            try {
                const bool better=std::ranges::all_of(patch.added,[&](const Cell& cell){
                    return quality(cell)>=old+1e-5;
                });
                if (!better || !complex.valid(patch,std::nullopt,staged)) complex.vertices.pop_back();
                else {
                    complex.charge();
                    if (!nativeAdmits(seed) || !accept(patch,old)) complex.vertices.pop_back();
                }
            } catch (const TrialLimit&) {
                complex.vertices.pop_back();
                throw;
            }
        }
    } catch (const TrialLimit&) {
        // Work is charged before a transaction. Earlier committed cavities
        // are still conforming and can be used.
    }
    if (!moves) return false;
    mesh::Mesh result;
    const auto region=complex.cells.front().classification;
    for (Index i=0;i<mesh.nodeCount();++i) result.addNode(mesh.node(i),mesh.nodeClass(i));
    for (std::size_t i=mesh.nodeCount();i<complex.vertices.size();++i)
        result.addNode(complex.vertices[i].position,region);
    for (Index i=0;i<mesh.cellCount();++i) if (mesh.cell(i).type==mesh::CellType::Tri3)
        result.addCell(mesh::CellType::Tri3,std::span(mesh.cellNodes(i),3),mesh.cell(i).onEntity);
    for (std::size_t i=0;i<complex.cells.size();++i) if (complex.alive[i]) {
        const auto& tet=complex.cells[i];
        const auto& v=tet.vertices;
        result.addCell(mesh::CellType::Tet4,{Index(v[0]),Index(v[1]),Index(v[2]),Index(v[3])},tet.classification);
    }
    mesh.swap(result);
    return true;
}

} // namespace paacs::geo::detail
