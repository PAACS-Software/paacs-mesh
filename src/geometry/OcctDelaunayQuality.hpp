// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "../bnd/TetBistellarRecovery.hpp"
#include "OcctDelaunayRecone.hpp"
#include "paacs/mesh/Quality.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <span>
#include <vector>

namespace paacs::geo::detail {

// Improve mean ratio without degrading radius ratio in an embedded Tet4 complex.
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
    const auto meanRatio=[&](const std::array<Vec3,4>& p) {
        // Match the reported subtet mean ratio, including its determinant
        // fallback when the ordinary volume of a positive cell cancels.
        auto volume=(p[1]-p[0]).cross(p[2]-p[0]).dot(p[3]-p[0])/6.;
        if (std::isfinite(volume) && !(volume>0.))
            volume=predicates::orient3dDeterminantEstimate(p[0],p[1],p[2],p[3])/6.;
        double sum=0.;
        for (unsigned i=0;i<4;++i) for (unsigned j=i+1;j<4;++j)
            sum+=(p[j]-p[i]).norm2();
        if (!(volume>0.) || !(sum>0.)) return 0.;
        const auto q=12.*std::pow(3.*volume,2./3.)/sum;
        return std::isfinite(q)?q:0.;
    };
    const auto swapFace=[&](Complex::Face face) {
        live.charge();
        if (tetCount>=maximumTets || live.cells.size()>live.maximumStoredCells ||
            live.maximumStoredCells-live.cells.size()<3) return false;
        // The existing bistellar builder certifies a convex two-cell union
        // or validates its complete oriented replacement shell. Its vertices
        // remain fixed, including all original CAD boundary vertices.
        auto patch=live.flip23(Complex::key(face),std::nullopt);
        if (!patch) return false;
        double beforeQuality=1.,beforeMeanRatio=1.,beforeVolume=0.;
        double afterQuality=1.,afterMeanRatio=1.,afterVolume=0.;
        std::array<unsigned,2> beforePoor{},afterPoor{};
        for (const auto id:patch->removed) {
            live.charge(4);
            const auto p=points(live.cells[id]);
            const auto mean=meanRatio(p);
            beforeQuality=std::min(beforeQuality,quality(p));
            beforeMeanRatio=std::min(beforeMeanRatio,mean);
            beforePoor[0]+=mean<.1;
            beforePoor[1]+=mean<.25;
            beforeVolume+=signedVolume(p);
        }
        for (const auto& cell:patch->added) {
            live.charge(4);
            const auto p=points(cell);
            const auto mean=meanRatio(p);
            afterQuality=std::min(afterQuality,quality(p));
            afterMeanRatio=std::min(afterMeanRatio,mean);
            afterPoor[0]+=mean<.1;
            afterPoor[1]+=mean<.25;
            afterVolume+=signedVolume(p);
        }
        if (!(beforeVolume>0.) || !(afterVolume>0.) ||
            !(afterQuality>beforeQuality+1e-12) || afterMeanRatio<beforeMeanRatio ||
            afterPoor[0]>beforePoor[0] || afterPoor[1]>beforePoor[1] ||
            !(std::abs(afterVolume-beforeVolume)<=1e-10*(afterVolume+beforeVolume)) ||
            !live.valid(*patch,std::nullopt) || !live.preservesSurvivors(*patch)) return false;
        for (const auto& cell:patch->added) {
            live.charge();
            if (!admitCentroid(centroid(points(cell)))) return false;
        }
        live.apply(*patch);
        tetCount=tetCount-patch->removed.size()+patch->added.size();
        ++edits;
        return true;
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
        double beforeMeanRatio=1.;
        std::array<unsigned,2> beforePoor{};
        double beforeVolume=0.;
        for (const auto owner:patch.removed) {
            live.charge(4);
            const auto& cell=live.cells[owner];
            const auto p=points(cell);
            const auto mean=meanRatio(p);
            beforeQuality=std::min(beforeQuality,quality(p));
            beforeMeanRatio=std::min(beforeMeanRatio,mean);
            beforePoor[0]+=mean<.1;
            beforePoor[1]+=mean<.25;
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
            std::array<unsigned,2> poor{};
            double quality=0.,meanRatio=0.,volume=0.;
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
                const auto meanA=meanRatio(a), meanB=meanRatio(b);
                cone.quality=std::min(quality(a),quality(b));
                cone.meanRatio=std::min(meanA,meanB);
                cone.poor={unsigned(meanA<.1)+unsigned(meanB<.1),
                           unsigned(meanA<.25)+unsigned(meanB<.25)};
                cone.volume=signedVolume(a)+signedVolume(b);
            }

        // A seven-vertex ring has at most 42 triangulations. Rank all geometric
        // improvements and certify the strongest admissible one before editing.
        struct Candidate {double meanRatio,quality; std::vector<Cell> cells;};
        std::vector<Candidate> candidates;
        using Interval=std::array<std::size_t,2>;
        const auto enumerate=[&](auto&& self,std::array<Interval,7> remaining,
                                 unsigned pending,double score,double meanScore,double volume,
                                 std::array<unsigned,2> poor)->void {
            if (!pending) {
                if (poor[0]<=beforePoor[0] && poor[1]<=beforePoor[1] &&
                    std::abs(volume-beforeVolume)<=1e-10*(volume+beforeVolume))
                    candidates.push_back({meanScore,score,patch.added});
                return;
            }
            const auto [i,k]=remaining[--pending];
            for (auto j=i+1;j<k;++j) {
                live.charge();
                const auto& cone=cones[i][j][k];
                const auto nextScore=std::min(score,cone.quality);
                const auto nextMeanScore=std::min(meanScore,cone.meanRatio);
                if (!(nextScore>beforeQuality+1e-12) || nextMeanScore<beforeMeanRatio) continue;
                auto next=pending;
                if (j-i>1) remaining[next++]={i,j};
                if (k-j>1) remaining[next++]={j,k};
                patch.added.insert(patch.added.end(),cone.cells.begin(),cone.cells.end());
                self(self,remaining,next,nextScore,nextMeanScore,volume+cone.volume,
                    {poor[0]+cone.poor[0],poor[1]+cone.poor[1]});
                patch.added.resize(patch.added.size()-2);
            }
        };
        std::array<Interval,7> remaining{{{0,count-1}}};
        enumerate(enumerate,remaining,1,1.,1.,0.,{});
        std::stable_sort(candidates.begin(),candidates.end(),[](const auto& a,const auto& b) {
            return a.meanRatio==b.meanRatio?a.quality>b.quality:a.meanRatio>b.meanRatio;
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

    const auto reconeEdge=[&](Complex::Edge edge) {
        live.charge();
        if (live.protectedEdges.contains(Complex::key(edge))) return false;
        for (const auto apex:occtDelaunayEdgeReconeApices(live,edge,fixedVertices,8)) {
            live.charge();
            auto patch=buildOcctDelaunayEdgeRecone(live,edge,apex,64);
            if (!patch || patch->removed.size()>tetCount ||
                live.cells.size()>live.maximumStoredCells ||
                patch->added.size()>live.maximumStoredCells-live.cells.size()) continue;
            const auto retainedTets=tetCount-patch->removed.size();
            if (retainedTets>maximumTets || patch->added.size()>maximumTets-retainedTets) continue;
            double beforeQuality=1.,beforeMeanRatio=1.,beforeVolume=0.;
            double afterQuality=1.,afterMeanRatio=1.,afterVolume=0.;
            std::array<unsigned,2> beforePoor{},afterPoor{};
            for (const auto id:patch->removed) {
                live.charge(4);
                const auto p=points(live.cells[id]);
                const auto mean=meanRatio(p);
                beforeQuality=std::min(beforeQuality,quality(p));
                beforeMeanRatio=std::min(beforeMeanRatio,mean);
                beforePoor[0]+=mean<.1;
                beforePoor[1]+=mean<.25;
                beforeVolume+=signedVolume(p);
            }
            bool clearsCriticalFloor=true;
            for (const auto& cell:patch->added) {
                live.charge(4);
                const auto p=points(cell);
                const auto mean=meanRatio(p);
                clearsCriticalFloor=clearsCriticalFloor && mean>=.01;
                afterQuality=std::min(afterQuality,quality(p));
                afterMeanRatio=std::min(afterMeanRatio,mean);
                afterPoor[0]+=mean<.1;
                afterPoor[1]+=mean<.25;
                afterVolume+=signedVolume(p);
            }
            if (!clearsCriticalFloor || !(beforeVolume>0.) || !(afterVolume>0.) ||
                !(afterQuality>beforeQuality+1e-12) || afterMeanRatio<beforeMeanRatio ||
                afterPoor[0]>beforePoor[0] || afterPoor[1]>beforePoor[1] ||
                !(std::abs(afterVolume-beforeVolume)<=1e-10*(afterVolume+beforeVolume)) ||
                !live.valid(*patch,std::nullopt) || !live.preservesSurvivors(*patch)) continue;
            if (!std::ranges::all_of(patch->added,[&](const Cell& cell) {
                live.charge();
                return admitCentroid(centroid(points(cell)));
            })) continue;
            live.apply(*patch);
            tetCount=retainedTets+patch->added.size();
            ++edits;
            return true;
        }
        return false;
    };

    // A boundary-only sliver cannot be improved by moving an interior vertex.
    // Repair the worst critical cells before relocation or ordinary edge
    // searches consume their shared work and original-solid query allowance.
    for (unsigned pass=0;pass<2 && edits<maximumEdits;++pass) {
        std::vector<std::pair<double,std::size_t>> critical;
        for (std::size_t id=0;id<live.cells.size();++id) {
            live.charge();
            if (!live.alive[id]) continue;
            const auto mean=meanRatio(points(live.cells[id]));
            if (mean<.01) critical.emplace_back(mean,id);
        }
        std::sort(critical.begin(),critical.end());
        const auto previous=edits;
        for (const auto& [mean,id]:critical) {
            (void)mean;
            if (edits>=maximumEdits) break;
            live.charge();
            if (!live.alive[id] || meanRatio(points(live.cells[id]))>=.01) continue;
            bool changed=false;
            for (const auto face:Complex::faceKeys(live.cells[id])) {
                if (swapFace(face)) {changed=true;break;}
            }
            if (!changed) for (const auto edge:Complex::edgeKeys(live.cells[id]))
                if (swapEdge(edge)) {changed=true;break;}
            // A nonconvex link can block every direct flip. Reconstruct its
            // complete edge star and bounded visibility closure from an
            // existing adjacent apex, with every point and PLC facet fixed.
            if (!changed) for (const auto edge:Complex::edgeKeys(live.cells[id]))
                if (reconeEdge(edge)) break;
        }
        if (edits==previous) break;
    }

    // Give interior spacing a bounded sweep before the ordinary edge search spends
    // its work allowance. Repeat only after topology changes make new stars.
    // A closed positive star retains its oriented link and material volume.
    const auto relocate=[&] {
        if (!live.walkConstraints || edits>=maximumEdits) return;
        // Bound successful moves separately so relocation leaves room for
        // topology edits within the shared native-query and work budgets.
        constexpr unsigned maximumMovesPerSweep=64;
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
        // Rank the eligible stars once per sweep, worst reported quality
        // first. The live star is measured again before each proposed move.
        std::vector<std::pair<double,Vertex>> ordered;
        for (Vertex vertex=0;vertex<stars.size();++vertex) {
            live.charge();
            const auto& star=stars[vertex];
            if (star.size()<4 || star.size()>96) continue;
            double minimumMeanRatio=1., minimumQuality=1.;
            for (const auto id:star) {
                live.charge();
                const auto p=points(live.cells[id]);
                minimumMeanRatio=std::min(minimumMeanRatio,meanRatio(p));
                minimumQuality=std::min(minimumQuality,quality(p));
            }
            if (minimumQuality<searchBelow) ordered.emplace_back(minimumMeanRatio,vertex);
        }
        std::sort(ordered.begin(),ordered.end());
        live.index={};
        unsigned moves=0;
        for (const auto& entry:ordered) {
            if (moves>=maximumMovesPerSweep || edits>=maximumEdits) break;
            live.charge();
            const auto vertex=entry.second;
            const auto& star=stars[vertex];
            const auto position=live.vertices[vertex].position;
            double beforeVolume=0., beforeQuality=1., beforeMeanRatio=1.;
            std::array<unsigned,2> beforePoor{};
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
                const auto mean=meanRatio(p);
                beforeQuality=std::min(beforeQuality,quality(p));
                beforeMeanRatio=std::min(beforeMeanRatio,mean);
                beforePoor[0]+=mean<.1;
                beforePoor[1]+=mean<.25;
                beforeVolume+=volume;
                weightedOffset=weightedOffset+(centroid(p)-position)*volume;
            }
            if (!closed || beforeQuality>=searchBelow || !(beforeVolume>0.)) continue;
            const auto offset=weightedOffset/beforeVolume;
            Vec3 moved{};
            double bestMeanRatio=beforeMeanRatio, bestQuality=beforeQuality;
            std::vector<Vec3> bestCentroids;
            // A full centroid step can overshoot the positive star's useful
            // region. Rank a bounded set of geometric trials before native queries.
            for (const auto alpha:std::array{1.,.5,.25,.125,.0625}) {
                const auto candidate=position+offset*alpha;
                if (!bnd::tet_facet_detail::finite(candidate) ||
                    (candidate.x==position.x && candidate.y==position.y && candidate.z==position.z)) continue;
                double afterVolume=0., afterQuality=1., afterMeanRatio=1.;
                std::array<unsigned,2> afterPoor{};
                std::vector<Vec3> centroids;
                centroids.reserve(star.size());
                bool admissible=true;
                for (const auto id:star) {
                    live.charge();
                    const auto& cell=live.cells[id];
                    auto p=points(cell);
                    p[std::ranges::find(cell.vertices,vertex)-cell.vertices.begin()]=candidate;
                    if (predicates::orient3dFiltered(p[0],p[1],p[2],p[3])!=Sign::Positive) {
                        admissible=false;break;
                    }
                    const auto mean=meanRatio(p);
                    afterQuality=std::min(afterQuality,quality(p));
                    afterMeanRatio=std::min(afterMeanRatio,mean);
                    afterPoor[0]+=mean<.1;
                    afterPoor[1]+=mean<.25;
                    if (afterQuality<beforeQuality || afterMeanRatio<beforeMeanRatio) {
                        admissible=false;break;
                    }
                    afterVolume+=signedVolume(p);
                    centroids.push_back(centroid(p));
                }
                if (!admissible || afterPoor[0]>beforePoor[0] || afterPoor[1]>beforePoor[1] ||
                    !(std::abs(afterVolume-beforeVolume)<=1e-10*beforeVolume) ||
                    afterMeanRatio<bestMeanRatio ||
                    (afterMeanRatio==bestMeanRatio && afterQuality<=bestQuality)) continue;
                moved=candidate;
                bestMeanRatio=afterMeanRatio;
                bestQuality=afterQuality;
                bestCentroids=std::move(centroids);
            }
            if (bestCentroids.empty() || !admitNode(moved) ||
                !std::ranges::all_of(bestCentroids,[&](Vec3 p){return admitCentroid(p);})) continue;
            live.charge(star.size());
            live.vertices[vertex].position=moved;
            for (const auto id:star) live.cellBoxes[id]=live.bounds(live.cells[id].vertices);
            ++edits;
            ++moves;
        }
    };
    relocate();
    const auto beforeSwaps=edits;

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
    if (edits!=beforeSwaps) relocate();
}

} // namespace paacs::geo::detail
