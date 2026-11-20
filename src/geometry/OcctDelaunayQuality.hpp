// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "../bnd/TetBistellarRecovery.hpp"
#include "paacs/mesh/Quality.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <span>
#include <vector>

namespace paacs::geo::detail {

// Improve radius ratio in an already embedded, positive Tet4 complex.
// The caller supplies a fixed-vertex mask and original-solid witnesses.
template<class AdmitCentroid,class AdmitNode>
void optimizeOcctDelaunayQuality(bnd::detail::BistellarComplex& live,
    std::span<const unsigned char> fixedVertices, AdmitCentroid&& admitCentroid,
    AdmitNode&& admitNode, unsigned& edits, std::uint64_t maximumTets=UINT64_MAX) {
    using Complex=bnd::detail::BistellarComplex;
    using Cell=bnd::TetFacetCell;
    using Vertex=bnd::TetFacetVertexId;
    using Sign=predicates::Sign;
    using bnd::tet_facet_detail::outwardFace;
    constexpr double searchBelow=.3;
    constexpr unsigned maximumEdits=1024;
    std::size_t tetCount=0;
    for (const auto alive:live.alive) {
        live.charge();
        tetCount+=alive!=0;
    }

    const auto points=[&](const Cell& cell) {
        const auto& v=cell.vertices;
        return std::array{live.vertices[v[0]].position,live.vertices[v[1]].position,
                          live.vertices[v[2]].position,live.vertices[v[3]].position};
    };
    const auto quality=[](const std::array<Vec3,4>& p) {
        const auto q=mesh::tetQuality(p[0],p[1],p[2],p[3]);
        return std::isfinite(q)?q:0.;
    };
    const auto centroid=[](const std::array<Vec3,4>& p) {
        Vec3 center{};
        for (const auto point:p) center=center+point*.25;
        return center;
    };
    const auto signedVolume=[](const std::array<Vec3,4>& p) {
        return mesh::tetVolume(p[0],p[1],p[2],p[3]);
    };
    const auto swapEdge=[&](Complex::Edge edge) {
        live.charge();
        const auto owners=live.edges.find(edge);
        if (owners==live.edges.end() || live.protectedEdges.contains(edge)) return false;
        const auto count=owners->second.size();
        if (count<3 || count>7 || tetCount-count+2*(count-2)>maximumTets ||
            2*(count-2)>live.maximumStoredCells-live.cells.size()) return false;

        // Each tet contributes one link edge. A valid removable edge has one
        // closed link, with a different ring vertex between each pair of tets.
        std::map<Vertex,std::vector<Vertex>> link;
        Complex::Proposal patch;
        patch.removed.assign(owners->second.begin(),owners->second.end());
        patch.kind=count==3?32:count==4?44:Complex::nmKind;
        double beforeQuality=1.;
        double beforeVolume=0.;
        for (const auto owner:patch.removed) {
            live.charge(4);
            const auto& cell=live.cells[owner];
            const auto p=points(cell);
            beforeQuality=std::min(beforeQuality,quality(p));
            beforeVolume+=signedVolume(p);
            Complex::Edge opposite{};
            unsigned next=0;
            for (const auto vertex:cell.vertices)
                if (vertex!=edge[0] && vertex!=edge[1]) opposite[next++]=vertex;
            link[opposite[0]].push_back(opposite[1]);
            link[opposite[1]].push_back(opposite[0]);
        }
        if (!(beforeVolume>0.) || link.size()!=count) return false;
        for (const auto& [vertex,adjacent]:link) {
            (void)vertex;
            if (adjacent.size()!=2 || adjacent.front()==adjacent.back()) return false;
        }
        std::vector<Vertex> ring{link.begin()->first};
        auto previous=ring.front();
        auto current=link.begin()->second.front();
        while (current!=ring.front() && ring.size()<count) {
            ring.push_back(current);
            const auto& adjacent=link.at(current);
            const auto next=adjacent.front()==previous?adjacent.back():adjacent.front();
            previous=current;
            current=next;
        }
        if (ring.size()!=count || current!=ring.front()) return false;

        struct Cone {
            std::array<Cell,2> cells{};
            double quality=0.,volume=0.;
        };
        std::array<std::array<std::array<Cone,7>,7>,7> cones{};
        for (std::size_t i=0;i+2<count;++i) for (std::size_t j=i+1;j+1<count;++j)
            for (std::size_t k=j+1;k<count;++k) {
                live.charge(2);
                auto& cone=cones[i][j][k];
                const auto side=[&](Vertex v) {
                    return predicates::orient3dFiltered(live.vertices[ring[i]].position,
                        live.vertices[ring[j]].position,live.vertices[ring[k]].position,
                        live.vertices[v].position);
                };
                const auto first=side(edge[0]), second=side(edge[1]);
                if (first==Sign::Zero || second==Sign::Zero || first==second) continue;
                cone.cells[0].vertices={ring[i],ring[j],ring[k],edge[0]};
                cone.cells[1].vertices={ring[i],ring[j],ring[k],edge[1]};
                auto& reversed=cone.cells[first==Sign::Negative?0:1].vertices;
                std::swap(reversed[1],reversed[2]);
                const auto a=points(cone.cells[0]), b=points(cone.cells[1]);
                cone.quality=std::min(quality(a),quality(b));
                cone.volume=signedVolume(a)+signedVolume(b);
            }

        // A seven-vertex ring has at most 42 triangulations. Rank all geometric
        // improvements and certify the strongest admissible one before editing.
        struct Candidate {double quality; std::vector<Cell> cells;};
        std::vector<Candidate> candidates;
        using Interval=std::array<std::size_t,2>;
        const auto enumerate=[&](auto&& self,std::array<Interval,7> remaining,
                                 unsigned pending,double score,double volume)->void {
            if (!pending) {
                if (std::abs(volume-beforeVolume)<=1e-10*(volume+beforeVolume))
                    candidates.push_back({score,patch.added});
                return;
            }
            const auto [i,k]=remaining[--pending];
            for (auto j=i+1;j<k;++j) {
                live.charge();
                const auto& cone=cones[i][j][k];
                const auto nextScore=std::min(score,cone.quality);
                if (!(nextScore>beforeQuality+1e-12)) continue;
                auto next=pending;
                if (j-i>1) remaining[next++]={i,j};
                if (k-j>1) remaining[next++]={j,k};
                patch.added.insert(patch.added.end(),cone.cells.begin(),cone.cells.end());
                self(self,remaining,next,nextScore,volume+cone.volume);
                patch.added.resize(patch.added.size()-2);
            }
        };
        std::array<Interval,7> remaining{{{0,count-1}}};
        enumerate(enumerate,remaining,1,1.,0.);
        std::stable_sort(candidates.begin(),candidates.end(),[](const auto& a,const auto& b) {
            return a.quality>b.quality;
        });
        for (auto& candidate:candidates) {
            patch.added=std::move(candidate.cells);
            if (!live.valid(patch,std::nullopt) || !live.preservesSurvivors(patch)) continue;
            if (!std::ranges::all_of(patch.added,[&](const Cell& cell) {
                live.charge();
                return admitCentroid(centroid(points(cell)));
            })) continue;
            live.apply(patch);
            tetCount=tetCount-patch.removed.size()+patch.added.size();
            ++edits;
            return true;
        }
        return false;
    };

    for (unsigned pass=0;pass<2 && edits<maximumEdits;++pass) {
        const auto previous=edits;
        const auto limit=live.cells.size();
        for (std::size_t id=0;id<limit && edits<maximumEdits;++id) {
            live.charge();
            if (!live.alive[id] || quality(points(live.cells[id]))>=searchBelow) continue;
            for (const auto edge:Complex::edgeKeys(live.cells[id]))
                if (swapEdge(edge)) break;
        }
        if (edits==previous) break;
    }

    // After topology changes, gather interior stars once. A closed positive
    // star with the same oriented link bounds the same volume after any move
    // into the link's strict kernel. No spatial query is needed for that move.
    if (!live.walkConstraints || edits>=maximumEdits) return;
    std::vector<std::vector<std::size_t>> stars(live.vertices.size());
    std::vector<unsigned char> fixed(live.vertices.size());
    std::copy_n(fixedVertices.begin(),std::min(fixed.size(),fixedVertices.size()),fixed.begin());
    for (const auto& face:live.protectedFaces) for (const auto v:face) fixed[v]=1;
    for (const auto& edge:live.protectedEdges) for (const auto v:edge) fixed[v]=1;
    for (std::size_t id=0;id<live.cells.size();++id) {
        live.charge();
        if (live.alive[id]) for (const auto v:live.cells[id].vertices)
            if (!fixed[v]) stars[v].push_back(id);
    }
    live.index={};
    for (Vertex vertex=0;vertex<stars.size() && edits<maximumEdits;++vertex) {
        live.charge();
        const auto& star=stars[vertex];
        if (star.size()<4 || star.size()>96) continue;
        const auto position=live.vertices[vertex].position;
        double beforeVolume=0., beforeQuality=1.;
        Vec3 weightedOffset{};
        bool closed=true;
        for (const auto id:star) {
            live.charge(4);
            const auto& cell=live.cells[id];
            const auto side=std::ranges::find(cell.vertices,vertex)-cell.vertices.begin();
            for (unsigned f=0;f<4;++f) {
                if (f==side) continue;
                const auto other=live.neighbours[id][f];
                if (other==Complex::noCell || !live.alive[other] ||
                    live.cells[other].classification!=cell.classification ||
                    !std::ranges::binary_search(star,other)) closed=false;
            }
            if (!closed) break;
            const auto p=points(cell);
            const auto volume=signedVolume(p);
            if (!(volume>0.)) {closed=false;break;}
            beforeQuality=std::min(beforeQuality,quality(p));
            beforeVolume+=volume;
            weightedOffset=weightedOffset+(centroid(p)-position)*volume;
        }
        if (!closed || beforeQuality>=searchBelow || !(beforeVolume>0.)) continue;
        const auto moved=position+weightedOffset/beforeVolume;
        if (!bnd::tet_facet_detail::finite(moved) ||
            (moved.x==position.x && moved.y==position.y && moved.z==position.z)) continue;
        double afterVolume=0.;
        std::vector<Vec3> centroids;
        for (const auto id:star) {
            live.charge();
            const auto& cell=live.cells[id];
            auto p=points(cell);
            p[std::ranges::find(cell.vertices,vertex)-cell.vertices.begin()]=moved;
            if (predicates::orient3dFiltered(p[0],p[1],p[2],p[3])!=Sign::Positive ||
                quality(p)<beforeQuality) {closed=false;break;}
            afterVolume+=signedVolume(p);
            centroids.push_back(centroid(p));
        }
        if (!closed || !(std::abs(afterVolume-beforeVolume)<=1e-10*beforeVolume) ||
            !admitNode(moved) || !std::ranges::all_of(centroids,[&](Vec3 p){return admitCentroid(p);})) continue;
        live.charge(star.size());
        live.vertices[vertex].position=moved;
        for (const auto id:star) live.cellBoxes[id]=live.bounds(live.cells[id].vertices);
        ++edits;
    }
}

} // namespace paacs::geo::detail
