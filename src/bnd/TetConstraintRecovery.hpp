// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "TetBistellarRecovery.hpp"
#include "TetEdgeGeneralCavityCadPolicy.hpp"
#include "paacs/algo/TetCandidate.hpp"
#include "paacs/bnd/TetEdgeGeneralCavity.hpp"
#include <atomic>
#include <chrono>

namespace paacs::bnd::detail {

struct TetConstraintRecoveryReport {
    const char* diagnostic = "not_attempted";
    std::uint64_t candidateCells=0, missingEdgesBefore=0, missingFacesBefore=0;
    std::uint64_t missingEdgesAfter=0, missingFacesAfter=0, patches=0, seedCells=0, work=0;
    double seconds=0;
    std::uint64_t edgeWalks=0,facetWalks=0,walkFallbacks=0;
    std::uint64_t steinerAttempts=0,steinerWork=0,steinerNodes=0,retainedSteinerNodes=0;
    double initialSeconds=0,recoverySeconds=0,cavitySeconds=0,steinerSeconds=0,carveSeconds=0;
    double longStarSeconds=0,validationSeconds=0,survivorSeconds=0;
    std::uint64_t piercingPatches=0,convexBipyramids=0;
    std::uint64_t residualAttempts=0,residualLinkRejected=0,residualRecovered=0,residualExpanded=0;
};

namespace constraint_detail {
using Complex = BistellarComplex;
using Face = Complex::Face;
using Edge = Complex::Edge;
struct WorkLimit {};
struct Cancelled {};
using Clock=std::chrono::steady_clock;
struct StageTime {
    double& seconds;
    Clock::time_point began=Clock::now();
    ~StageTime(){seconds+=std::chrono::duration<double>(Clock::now()-began).count();}
};

inline constexpr std::size_t residualMaximumCells=40;
inline constexpr std::size_t residualMaximumVertices=18;
inline constexpr std::uint64_t residualMinimumAttemptWork=10'000;
inline constexpr std::uint64_t residualMaximumAttemptWork=180'000;

inline void recoverEdgeCavity(Complex& complex, Edge target, std::uint64_t workLimit,
                              std::uint64_t steinerWork, TetConstraintRecoveryReport& report,
                              const std::set<std::size_t>* prepared=nullptr) {
    if (!workLimit&&!steinerWork) return;
    StageTime timer{steinerWork?report.steinerSeconds:report.cavitySeconds};
    const auto discovered=prepared?std::set<std::size_t>{}:complex.edgeCavity(target);
    const auto& affected=prepared?*prepared:discovered;
    if (affected.empty()||affected.size()>(prepared?residualMaximumCells:32)) return;
    Complex::Proposal patch;
    patch.removed.assign(affected.begin(),affected.end());
    std::vector<TetFacetCell> old;
    std::vector<TetFacetVertexId> ids;
    std::set<Edge> edges;
    std::set<Face> faces;
    for (auto id:affected) {
        const auto& cell=complex.cells[id];
        old.push_back(cell);
        ids.insert(ids.end(),cell.vertices.begin(),cell.vertices.end());
        for (auto edge:complex.edgeKeys(cell)) if (complex.protectedEdges.contains(edge)) edges.insert(edge);
        for (auto face:complex.faceKeys(cell)) if (complex.protectedFaces.contains(face)) faces.insert(face);
    }
    std::sort(ids.begin(),ids.end());
    ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
    if(ids.size()>(prepared?residualMaximumVertices:14)||!std::ranges::binary_search(ids,target[0])||!std::ranges::binary_search(ids,target[1]))return;
    const auto local=[&](TetFacetVertexId id) {
        return static_cast<TetFacetVertexId>(std::ranges::lower_bound(ids,id)-ids.begin());
    };
    std::vector<TetFacetVertex> points;
    for(auto id:ids)points.push_back(complex.vertices[id]);
    for(auto& cell:old)for(auto& id:cell.vertices)id=local(id);
    std::vector<TetFacetEdgeKey> protectedEdges;
    std::vector<TetFacetProtectedFace> protectedFaces;
    for (auto edge:edges) protectedEdges.push_back({local(edge[0]),local(edge[1])});
    for (auto face:faces) {
        for(auto& id:face)id=local(id);
        protectedFaces.push_back({{face},protectedFaces.size()+1,{2,0}});
    }
    const TetEdgeGeneralCavityInput input{points,old,{local(target[0]),local(target[1])},protectedEdges,
        protectedFaces,old.front().classification};
    auto commit=[&](std::vector<TetFacetCell> cells,TetFacetVertexId staged=kInvalidTetFacetVertex) {
        for(auto& cell:cells)for(auto& id:cell.vertices)id=id==ids.size()?staged:ids[id];
        patch.added=std::move(cells);
        if(!complex.valid(patch,target,staged)||!complex.preservesSurvivors(patch))return false;
        complex.apply(patch);
        return true;
    };
    if(workLimit){
        const auto solve=[&] {
            if (!prepared) return solveTetEdgeGeneralCavity(input,{32,14,1001,1500,workLimit});
            // The existing solver has its own fixed work cap. Check the live
            // cancellation hook on both sides of that bounded call.
            complex.charge(0);
            ++report.residualAttempts;
            const ScopedCadGeneralPruning pruning(true);
            return solveTetEdgeGeneralCavity(input,{residualMaximumCells,residualMaximumVertices,
                                                    3060,8000,workLimit});
        };
        auto solved=solve();
        complex.charge(solved.metrics.workUnits);
        if(solved.status==TetEdgeGeneralCavityStatus::AllocationFailed)throw std::bad_alloc{};
        if(solved.solved()&&commit(std::move(solved.cells)))return;
    }
    if(!steinerWork||ids.size()>=14)return;
    const ScopedCadGeneralPruning pruning(true);
    auto solved=solveTetEdgeGeneralCavityOneSteiner(input,{32,14,1001,1500,steinerWork},20,steinerWork);
    report.steinerAttempts+=solved.candidateAttempts;
    report.steinerWork+=solved.cumulativeWorkUnits;
    complex.charge(solved.cumulativeWorkUnits);
    if(solved.status==TetEdgeGeneralCavityStatus::AllocationFailed)throw std::bad_alloc{};
    if(!solved.solved())return;
    const auto staged=static_cast<TetFacetVertexId>(complex.vertices.size());
    complex.vertices.push_back({solved.steinerPosition,std::uint64_t(staged)+1});
    try{
        if(commit(std::move(solved.cavity.cells),staged)){++report.steinerNodes;return;}
    }catch(...){complex.vertices.pop_back();throw;}
    complex.vertices.pop_back();
}

inline void recoverFace(Complex& complex, Face face, mesh::Classification region,
                        bool direct, bool longLinks=false) {
    const Edge guide{face[0],face[1]};
    if (!complex.edges.contains(guide)||!complex.edges.contains(Edge{face[0],face[2]})||
        !complex.edges.contains(Edge{face[1],face[2]})) return;
    TetFacetComplexView view{complex.vertices,{}, {face}, {},{},region};
    const auto conflicts=[&](const TetFacetCell& cell) {return complex.facetConflict(view,cell);};
    for (unsigned step=0;step<8&&!complex.faces.contains(face);++step) {
        const auto affected=complex.facetCavity(view);
        if (affected.empty()||affected.size()>32) return;
        std::vector<Complex::Proposal> proposals;
        if (direct) {
            // With all three target edges present, the direct 3-to-2 move
            // uses the target triangle as the complete link of an edge.
            std::set<Edge> nearbyEdges;
            for (auto id:affected) for (auto edge:complex.edgeKeys(complex.cells[id])) nearbyEdges.insert(edge);
            for (auto edge:nearbyEdges) {
                complex.charge();
                const auto& owners=complex.edges.at(edge);
                if (owners.size()!=3||complex.protectedEdges.contains(edge)) continue;
                std::set<TetFacetVertexId> link;
                for (auto id:owners) for (auto vertex:complex.cells[id].vertices)
                    if (vertex!=edge[0]&&vertex!=edge[1]) link.insert(vertex);
                if (!std::ranges::equal(link,face)) continue;
                proposals=complex.edgeFlips(edge,guide);
                if (!proposals.empty()&&complex.preservesSurvivors(proposals.front())) {
                    complex.apply(proposals.front());return;
                }
            }
            return;
        }
        proposals=complex.proposals(guide,affected,longLinks);
        int bestDelta=0;
        const Complex::Proposal* best=nullptr;
        for (const auto& patch:proposals) {
            int delta=0;
            bool createsFace=false;
            for (auto id:patch.removed) delta-=int(affected.contains(id));
            for (const auto& cell:patch.added) {
                delta+=int(conflicts(cell));
                for (auto key:complex.faceKeys(cell)) createsFace|=key==face;
            }
            if (!createsFace&&delta>=bestDelta) continue;
            if (!complex.preservesSurvivors(patch)) continue;
            best=&patch;
            if (createsFace) break;
            bestDelta=delta;
        }
        if (!best) return;
        complex.apply(*best);
    }
}

inline void recoverConstraintFlips(Complex& complex, std::span<const Edge> missingEdges,
    std::span<const Face> missingFaces, mesh::Classification region,
    std::uint64_t workLimit, TetConstraintRecoveryReport& report) {
    // Run the stronger pass even if the short-link pass gets stuck.
    for (const bool longLinks : {false,true}) {
        for (const auto edge : missingEdges) {
            complex.charge();
            if (complex.edges.contains(edge)) continue;
            if (!complex.recoverEdge(edge,longLinks?6:4,longLinks) && longLinks)
                recoverEdgeCavity(complex,edge,std::min<std::uint64_t>(12000,workLimit-report.work),0,report);
        }
        for (const auto face : missingFaces) {
            complex.charge();
            if (!complex.faces.contains(face)) recoverFace(complex,face,region,false,longLinks);
        }
    }
}

// Work on the complete star of an edge piercing the target. Generic facet
// moves rank large edge links against a guide segment already in the mesh.
// Here the link triangulation is scored against the triangle itself.
inline std::uint64_t recoverPiercedFace(Complex& complex, Face face, mesh::Classification region) {
    const Edge guide{face[0],face[1]};
    if(!complex.edges.contains(guide)||!complex.edges.contains({face[0],face[2]})||
       !complex.edges.contains({face[1],face[2]}))return 0;
    const TetFacetComplexView view{complex.vertices,{}, {face}, {},{},region};
    const auto points=tet_facet_detail::trianglePoints(view,view.target);
    std::uint64_t accepted=0;
    for(unsigned step=0;step<4&&!complex.faces.contains(face);++step){
        const auto affected=complex.facetCavity(view);
        if(affected.empty()||affected.size()>32)break;
        std::set<Edge> tested;
        std::optional<Complex::Proposal> best;
        int bestDelta=0;
        for(auto id:affected)for(auto edge:complex.edgeKeys(complex.cells[id])){
            complex.charge();
            if(complex.protectedEdges.contains(edge)||!tested.insert(edge).second)continue;
            if(!tet_facet_detail::segmentCrossesOpenTriangle(complex.vertices[edge[0]].position,
                    complex.vertices[edge[1]].position,points[0],points[1],points[2]))continue;
            for(auto& patch:complex.edgeFlips(edge,guide,&view)){
                int delta=0;
                for(auto old:patch.removed)delta-=complex.facetConflict(view,complex.cells[old]);
                bool creates=false;
                for(const auto& cell:patch.added){
                    delta+=complex.facetConflict(view,cell);
                    const auto faces=complex.faceKeys(cell);
                    creates|=std::ranges::find(faces,face)!=faces.end();
                }
                if(!creates&&delta>=bestDelta)continue;
                if(!complex.preservesSurvivors(patch))continue;
                if(creates){
                    complex.apply(patch);
                    return accepted+1;
                }
                best=std::move(patch);
                bestDelta=delta;
            }
        }
        if(!best)break;
        complex.apply(*best);
        ++accepted;
    }
    return accepted;
}

// Outward source facets label both sides, including interior void shells.
// Missing facets block intersecting cells entirely. Unlabelled components
// remain private and excluded from the partial seed.
inline bool carve(Complex& complex, const std::map<Face,Face>& source, const mesh::Mesh& boundary,
                  mesh::Classification region, bool complete, std::uint64_t maximumCells,
                  mesh::Mesh& destination, TetConstraintRecoveryReport& report) {
    enum Label : std::uint8_t { Unknown, Exterior, Material, Blocked };
    std::vector<Label> labels(complex.cells.size(),Unknown);
    std::vector<std::size_t> queue;
    TetFacetComplexView view{complex.vertices,{}, {},{}, {},region};
    for (const auto& [key,directed]:source) if (!complex.faces.contains(key)) {
        view.target={directed};
        complex.visit(complex.bounds(key),[&](std::size_t id) {
            if (labels[id]!=Blocked&&complex.facetConflict(view,complex.cells[id])) labels[id]=Blocked;
            return true;
        });
    }
    const auto label=[&](std::size_t id,Label value) {
        if (labels[id]==Blocked) return true;
        if (labels[id]==Unknown) {labels[id]=value;queue.push_back(id);}
        return labels[id]==value;
    };
    for (const auto& [key,directed]:source) {
        complex.charge();
        const auto owners=complex.faces.find(key);
        if (owners==complex.faces.end()) continue;
        if (owners->second.size()!=2) return false;
        for (auto id:owners->second) {
            for (unsigned side=0;side<4;++side) {
                const auto outward=tet_facet_detail::outwardFace(complex.cells[id],side);
                if (complex.key(outward)!=key) continue;
                if (!label(id,tet_facet_detail::faceSense(view,outward)==
                    tet_facet_detail::faceSense(view,directed)?Material:Exterior)) return false;
            }
        }
    }
    const auto sourceNodes=static_cast<TetFacetVertexId>(boundary.nodeCount());
    const auto isGuard=[&](TetFacetVertexId id){return id>=sourceNodes&&id<sourceNodes+8;};
    for (std::size_t id=0;id<complex.cells.size();++id) {
        complex.charge();
        if (complex.alive[id]&&std::ranges::any_of(complex.cells[id].vertices,isGuard)&&!label(id,Exterior)) return false;
    }
    for (std::size_t i=0;i<queue.size();++i) {
        const auto id=queue[i];
        complex.charge(5);
        for(unsigned side=0;side<4;++side){
            const auto neighbour=complex.neighbours[id][side];
            if(neighbour==Complex::noCell)continue;
            const auto key=complex.key(tet_facet_detail::outwardFace(complex.cells[id],side));
            if(!source.contains(key)&&!label(neighbour,labels[id]))return false;
        }
    }
    mesh::Mesh seed;
    std::vector<Index> nodeMap(complex.vertices.size(),kInvalidIndex);
    for (Index i=0;i<boundary.nodeCount();++i) {
        complex.charge();
        nodeMap[i]=seed.addNode(boundary.node(i),boundary.nodeClass(i));
    }
    for (std::size_t id=0;id<labels.size();++id) {
        complex.charge();
        if (!complex.alive[id]) continue;
        if (complete&&labels[id]==Unknown) return false;
        if (labels[id]!=Material) continue;
        if (std::uint64_t(seed.cellCount())>=maximumCells) throw WorkLimit{};
        std::array<Index,4> ids;
        for(unsigned j=0;j<4;++j){
            const auto v=complex.cells[id].vertices[j];
            if(isGuard(v))return false;
            if(nodeMap[v]==kInvalidIndex)nodeMap[v]=seed.addNode(complex.vertices[v].position,region);
            ids[j]=nodeMap[v];
        }
        seed.addCell(mesh::CellType::Tet4,ids,region);
    }
    if (!seed.cellCount()) return false;
    report.retainedSteinerNodes=seed.nodeCount()-boundary.nodeCount();
    destination.swap(seed);
    return true;
}
inline std::set<TetFacetVertexId> residualCavityNodes(
    Complex& complex, const std::set<std::size_t>& cells) {
    std::set<TetFacetVertexId> ids;
    for (const auto cell : cells) {
        complex.charge(4);
        ids.insert(complex.cells[cell].vertices.begin(),complex.cells[cell].vertices.end());
    }
    return ids;
}

// Only expand through unprotected shell triangles. Prefer neighbors that
// close a gap along two or more faces without enlarging the vertex set.
inline std::vector<std::size_t> residualNeighbours(
    Complex& complex, const std::set<std::size_t>& cells) {
    struct Contact { unsigned shared=0; bool protectedFace=false; };
    std::map<std::size_t,Contact> rows;
    const auto ids=residualCavityNodes(complex,cells);
    for (const auto cell : cells) {
        const auto keys=complex.faceKeys(complex.cells[cell]);
        for (unsigned side=0;side<4;++side) {
            complex.charge();
            const auto adjacent=complex.neighbours[cell][side];
            if (adjacent==Complex::noCell||cells.contains(adjacent)) continue;
            auto& row=rows[adjacent];
            ++row.shared;
            row.protectedFace|=complex.protectedFaces.contains(keys[side]);
        }
    }
    std::vector<std::tuple<unsigned,unsigned,std::size_t>> ranked;
    for (const auto& [cell,row] : rows) {
        complex.charge(5);
        if (row.protectedFace) continue;
        unsigned added=0;
        for (const auto vertex : complex.cells[cell].vertices) added+=!ids.contains(vertex);
        if (ids.size()+added<=residualMaximumVertices) ranked.emplace_back(4-row.shared,added,cell);
    }
    std::sort(ranked.begin(),ranked.end());
    std::vector<std::size_t> result;
    for (const auto& rank : ranked) result.push_back(std::get<2>(rank));
    return result;
}

// In a valid internal target-edge star, each link vertex has one entering and
// one leaving edge under the orientation about the target. A shell admitting
// no directed cycle cannot contain that star.
inline bool residualLinkPossible(
    Complex& complex, Edge target, const std::set<std::size_t>& cells,
    std::map<std::pair<Face,Edge>,bool>& contacts) {
    const auto ids=residualCavityNodes(complex,cells);
    if (!ids.contains(target[0])||!ids.contains(target[1])||ids.size()>residualMaximumVertices) return false;
    std::vector<TetFacetVertexId> rim;
    for (const auto id : ids) if (id!=target[0]&&id!=target[1]) rim.push_back(id);
    if (rim.size()<3) return false;
    std::map<Face,unsigned> uses;
    for (const auto id : cells) {
        complex.charge(4);
        for (const auto face : complex.faceKeys(complex.cells[id])) ++uses[face];
    }
    std::vector<std::pair<Face,Complex::Box>> shell;
    for (const auto& [face,count] : uses) {
        complex.charge();
        if (count==1||complex.protectedFaces.contains(face)) shell.emplace_back(face,complex.bounds(face));
    }
    const auto region=complex.cells[*cells.begin()].classification;
    TetFacetComplexView view{complex.vertices,{}, {},{},{},region};
    std::vector<std::vector<unsigned>> arcs(rim.size());
    std::vector<unsigned> indegrees(rim.size(),0);
    const auto point=[&](TetFacetVertexId id){return complex.vertices[id].position;};
    for (unsigned i=0;i<rim.size();++i) for(unsigned j=i+1;j<rim.size();++j) {
        complex.charge();
        const auto side=predicates::orient3d(point(target[0]),point(target[1]),point(rim[i]),point(rim[j]));
        if (side==predicates::Sign::Zero) continue;
        unsigned from=i,to=j;
        if (side==predicates::Sign::Negative) std::swap(from,to);
        const TetFacetCell tet{{target[0],target[1],rim[from],rim[to]},0,region};
        const auto box=complex.bounds(tet.vertices);
        bool valid=true;
        for (const auto& [face,bounds] : shell) {
            complex.charge();
            if (!complex.overlaps(box,bounds)) continue;
            const auto key=std::pair{face,Edge{rim[i],rim[j]}};
            const auto previous=contacts.find(key);
            if (previous!=contacts.end()) {
                if (previous->second) {valid=false; break;}
                continue;
            }
            view.target.vertices=face;
            const auto conflict=complex.facetConflict(view,tet);
            contacts.emplace(key,conflict);
            if (conflict) {valid=false; break;}
        }
        if (valid) {arcs[from].push_back(to); ++indegrees[to];}
    }
    std::vector<unsigned> queue;
    for (unsigned i=0;i<indegrees.size();++i) if (!indegrees[i]) queue.push_back(i);
    for (std::size_t i=0;i<queue.size();++i)
        for (const auto to : arcs[queue[i]]) if (!--indegrees[to]) queue.push_back(to);
    return queue.size()!=rim.size();
}

inline void recoverResidualCavities(
    Complex& complex, std::span<const Edge> missingEdges, std::span<const Face> missingFaces,
    mesh::Classification region, std::uint64_t workLimit, TetConstraintRecoveryReport& report) {
    const auto facetsReady=[&](Face face) {
        return complex.edges.contains(Edge{face[0],face[1]})&&
               complex.edges.contains(Edge{face[0],face[2]})&&
               complex.edges.contains(Edge{face[1],face[2]});
    };
    const auto remaining=[&]{return workLimit-report.work;};
    std::vector<std::tuple<std::size_t,std::size_t,Edge>> ranked;
    for (const auto edge : missingEdges) {
        complex.charge();
        if (remaining()<residualMinimumAttemptWork) break;
        if (complex.edges.contains(edge)) continue;
        const auto selected=complex.edgeCavity(edge);
        if (selected.empty()||selected.size()>residualMaximumCells) continue;
        const auto ids=residualCavityNodes(complex,selected);
        if (ids.size()<=residualMaximumVertices&&ids.contains(edge[0])&&ids.contains(edge[1]))
            ranked.emplace_back(ids.size(),selected.size(),edge);
    }
    std::sort(ranked.begin(),ranked.end());
    for (const auto& row : ranked) {
        const auto edge=std::get<2>(row);
        if (remaining()<residualMinimumAttemptWork) break;
        if (complex.edges.contains(edge)) continue;
        const auto selected=complex.edgeCavity(edge);
        if (selected.empty()||selected.size()>residualMaximumCells) continue;
        std::map<std::pair<Face,Edge>,bool> contacts;
        const auto attempt=[&](const std::set<std::size_t>& cells) {
            if (remaining()<residualMinimumAttemptWork||cells.size()>residualMaximumCells)
                return false;
            complex.charge(16*cells.size());
            if (!residualLinkPossible(complex,edge,cells,contacts)) {
                ++report.residualLinkRejected;
                return false;
            }
            recoverEdgeCavity(complex,edge,std::min<std::uint64_t>(residualMaximumAttemptWork,remaining()),0,report,&cells);
            return complex.edges.contains(edge);
        };
        if (attempt(selected)) {++report.residualRecovered; continue;}
        auto neighbours=residualNeighbours(complex,selected);
        if (neighbours.size()>8) neighbours.resize(8);
        std::set<std::set<std::size_t>> seen;
        bool recovered=false;
        for (const auto neighbour : neighbours) {
            auto expanded=selected;
            expanded.insert(neighbour);
            if (seen.insert(expanded).second&&attempt(expanded)) {recovered=true; break;}
        }
        if (!recovered) {
            unsigned pairs=0;
            for (const auto neighbour : neighbours) {
                if (pairs==8||remaining()<residualMinimumAttemptWork) break;
                auto expanded=selected;
                expanded.insert(neighbour);
                auto second=residualNeighbours(complex,expanded);
                if (second.size()>2) second.resize(2);
                for (const auto extra : second) {
                    auto two=expanded;
                    two.insert(extra);
                    if (!seen.insert(two).second) continue;
                    ++pairs;
                    if (attempt(two)) {recovered=true; break;}
                    if (pairs==8) break;
                }
                if (recovered) break;
            }
        }
        if (recovered) {++report.residualRecovered; ++report.residualExpanded;}
    }
    if (report.residualRecovered)
        for (const auto face : missingFaces) {
            complex.charge();
            if (remaining()<residualMinimumAttemptWork) break;
            if (!complex.faces.contains(face)&&facetsReady(face)) recoverFace(complex,face,region,true);
            if (!complex.faces.contains(face)&&facetsReady(face)) recoverFace(complex,face,region,false,true);
        }
}
}

// The caller supplies a closed, consistently oriented, embedded Tri3
// boundary. Output is a private Tet4 seed, never an audited publication.
inline TetConstraintRecoveryReport recoverTetConstraints(const mesh::Mesh& boundary,
    mesh::Classification region, std::uint64_t maximumWork, std::uint64_t maximumBytes,
    std::uint64_t maximumCells, const std::atomic_bool* cancellation, mesh::Mesh& destination,
    std::uint64_t maximumSteinerNodes=64, bool profile=false, bool expandResidualEdges=true) {
    using namespace constraint_detail;
    TetConstraintRecoveryReport report;
    const auto start=Clock::now();
    std::uint64_t workLimit=maximumWork;
    const auto charge=[&](std::uint64_t n) {
        if (cancellation&&cancellation->load()) throw Cancelled{};
        if (n>workLimit-report.work) throw WorkLimit{};
        report.work+=n;
    };
    const auto run=[&]() -> const char* {
        charge(0);
        const auto ns=boundary.nodeCount(),nf=boundary.cellCount();
        if (ns<4||nf<4||&boundary==&destination) return "invalid_source";
        if (ns>12000||nf>24000||maximumWork<100000||maximumBytes<64ULL*1024*1024)
            return "admission_budget";
        Complex complex;
        complex.chargeWork=charge;
        complex.maximumStoredCells=static_cast<std::size_t>((maximumBytes-32ULL*1024*1024)/1536);
        complex.useLookahead=false;
        complex.walkConstraints=true;
        complex.profile=profile;
        std::map<Face,Face> source;
        for (Index i=0;i<nf;++i) {
            charge(4);
            if (boundary.cell(i).type!=mesh::CellType::Tri3) return "invalid_source";
            const auto* v=boundary.cellNodes(i);
            if (std::any_of(v,v+3,[&](Index id){return id<0||id>=ns;})) return "invalid_source";
            const Face face{TetFacetVertexId(v[0]),TetFacetVertexId(v[1]),TetFacetVertexId(v[2])};
            const auto key=complex.key(face);
            if (!source.emplace(key,face).second) return "invalid_source";
            complex.requiredFaces.insert(key);
            for (unsigned j=0;j<3;++j) complex.requiredEdges.insert(complex.key(Edge{face[j],face[(j+1)%3]}));
        }
        std::vector<Vec3> points(boundary.nodes().begin(),boundary.nodes().end());
        Vec3 lo=points.front(),hi=lo;
        for (const auto p:points) {
            charge(1);
            lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};
            hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};
        }
        const auto span=std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z});
        if (!(std::isfinite(span)&&span>0)) return "invalid_source";
        for (unsigned i=0;i<8;++i) points.push_back({(i&1)?hi.x+span:lo.x-span,
            (i&2)?hi.y+span:lo.y-span,(i&4)?hi.z+span:lo.z-span});
        auto built=algo::buildTetCandidate(points);
        const auto& m=built.metrics;
        const auto candidateWork=m.construction.conflictInsphereTests+m.construction.conflictOrientationTests+
            m.validationFaceUses+m.validationHullTests+m.candidateCells+m.discardedCandidateWork;
        // The shared Delaunay builder returns its counters on completion.
        // Charge the reserved allowance even when its construction exhausts it.
        const auto charged=std::min(candidateWork,maximumWork-report.work);
        charge(charged);
        if (built.reason==algo::TetCandidateReason::AllocationFailed) throw std::bad_alloc{};
        if (charged!=candidateWork) throw WorkLimit{};
        if (!built.built()||built.candidate->points.size()!=points.size()) return "delaunay_refused";
        report.candidateCells=m.candidateCells;
        if (report.candidateCells>complex.maximumStoredCells||
            m.constructionCandidateConservativeBytes>maximumBytes) return "admission_budget";
        for (auto p:points) complex.vertices.push_back({p,complex.vertices.size()+1});
        for (const auto& src:built.candidate->cells) {
            TetFacetCell tet{{},0,region};
            for (unsigned j=0;j<4;++j) tet.vertices[j]=TetFacetVertexId(built.candidate->pointToInput[src.vertices[j]]);
            complex.insert(tet);
        }
        built.candidate.reset();
        std::vector<Edge> missingEdges;
        std::vector<Face> missingFaces;
        for (auto edge:complex.requiredEdges) {
            charge(1);
            if (complex.edges.contains(edge)) complex.protectedEdges.insert(edge); else missingEdges.push_back(edge);
        }
        for (auto face:complex.requiredFaces) {
            charge(1);
            if (complex.faces.contains(face)) complex.protectedFaces.insert(face); else missingFaces.push_back(face);
        }
        report.missingEdgesBefore=missingEdges.size();
        report.missingFacesBefore=missingFaces.size();
        const auto lengthSquared=[&](Edge edge) {
            const auto d=points[edge[0]]-points[edge[1]];
            return d.dot(d);
        };
        std::sort(missingEdges.begin(),missingEdges.end(),[&](Edge a,Edge b) {
            charge(1);
            return std::tuple(lengthSquared(a),a)<std::tuple(lengthSquared(b),b);
        });
        bool limited=false;
        workLimit=std::max(report.work,maximumWork-maximumWork/3);
        report.initialSeconds=std::chrono::duration<double>(Clock::now()-start).count();
        const auto recoveryStart=Clock::now();
        try {
            for (auto edge:missingEdges) {
                charge(1);
                if (!complex.edges.contains(edge)) complex.recoverDirectEdge(edge);
            }
            for (auto face:missingFaces) {
                charge(1);
                if (!complex.faces.contains(face)) recoverFace(complex,face,region,true);
            }
            recoverConstraintFlips(complex,missingEdges,missingFaces,region,workLimit,report);
            // Reserve the broad Steiner search for the constraints that
            // survived the flip and existing-node passes.
            const auto steinerLimit=std::min<std::uint64_t>(maximumWork/8,20'000'000);
            for (auto edge:missingEdges) {
                if(report.steinerNodes>=maximumSteinerNodes||report.steinerWork>=steinerLimit)break;
                charge(1);
                if(complex.edges.contains(edge))continue;
                recoverEdgeCavity(complex,edge,0,std::min({std::uint64_t{160000},
                    workLimit-report.work,steinerLimit-report.steinerWork}),report);
            }
            if(report.steinerNodes)
                for(auto face:missingFaces)if(!complex.faces.contains(face))recoverFace(complex,face,region,true);
            for(auto face:missingFaces){
                charge(1);
                if(!complex.faces.contains(face))
                    report.piercingPatches+=recoverPiercedFace(complex,face,region);
            }
            if (expandResidualEdges) {
                const auto previousLimit=workLimit;
                const auto allowance=std::min({previousLimit-report.work,maximumWork/8,
                                                std::uint64_t{24'000'000}});
                workLimit=report.work+allowance;
                try {
                    recoverResidualCavities(complex,missingEdges,missingFaces,region,workLimit,report);
                } catch (const WorkLimit&) {
                    // Keep any already certified patches and preserve the carve reserve.
                }
                workLimit=previousLimit;
            }
        } catch (const WorkLimit&) { limited=true; }
          catch (const std::length_error&) { limited=true; }
        report.recoverySeconds=std::chrono::duration<double>(Clock::now()-recoveryStart).count();
        workLimit=maximumWork;
        report.patches=complex.metrics.accepted;
        report.convexBipyramids=complex.metrics.convexBipyramids;
        report.longStarSeconds=complex.metrics.longStarSeconds;
        report.validationSeconds=complex.metrics.validationSeconds;
        report.survivorSeconds=complex.metrics.survivorSeconds;
        report.edgeWalks=complex.edgeWalks;report.facetWalks=complex.facetWalks;report.walkFallbacks=complex.walkFallbacks;
        for (auto edge:missingEdges) report.missingEdgesAfter+=!complex.edges.contains(edge);
        for (auto face:missingFaces) report.missingFacesAfter+=!complex.faces.contains(face);
        const bool complete=report.missingFacesAfter==0;
        {
            StageTime timer{report.carveSeconds};
            if (!carve(complex,source,boundary,region,complete,maximumCells,destination,report)) return "material_carve_refused";
        }
        report.seedCells=destination.cellCount();
        return complete?"complete_seed":limited?"partial_seed_budget":"partial_seed";
    };
    try {report.diagnostic=run();}
    catch (const Cancelled&) {report.diagnostic="cancelled";}
    catch (const WorkLimit&) {report.diagnostic="work_budget";}
    catch (const std::length_error&) {report.diagnostic="storage_budget";}
    catch (const std::bad_alloc&) {report.diagnostic="allocation_failed";}
    report.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return report;
}
}
