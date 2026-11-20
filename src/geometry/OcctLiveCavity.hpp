// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "../bnd/TetBistellarRecovery.hpp"
#include "OcctLiveCavityRelocation.hpp"
#include "paacs/mesh/Mesh.hpp"
#include <span>

namespace paacs::geo::detail {

struct LiveCavityReceipt {
    std::uint64_t initialRejected=0, remainingRejected=0, patched=0, flips=0,
        seeds=0, relocations=0, proposals=0, exactPatches=0, rejectedByCad=0, queries=0, work=0;
};

// The caller supplies a certified discrete fill and the original-CAD mask
// in Tet4 order. Local changes keep the oriented cavity boundary fixed.
template<class Admit, class Charge>
LiveCavityReceipt repairLiveCavities(mesh::Mesh& mesh, std::span<const unsigned char> admitted,
    Admit admit, Charge charge, std::uint64_t maxWork, std::uint64_t maxQueries,
    std::uint64_t maxNodes, std::uint64_t maxTets, bool tryFeatureEdges=false) {
    using Complex=bnd::detail::BistellarComplex;
    using Cell=bnd::TetFacetCell;
    using Vertex=bnd::TetFacetVertexId;
    using Sign=predicates::Sign;
    using bnd::tet_facet_detail::outwardFace;
    struct Limit {};
    constexpr std::size_t maxPatch=18;
    LiveCavityReceipt out;
    const auto count=std::size_t(mesh.countCells(mesh::CellType::Tet4));
    if(count!=admitted.size()||count>maxTets||std::uint64_t(mesh.nodeCount())>maxNodes)return out;
    Complex live;
    live.walkConstraints=true;
    live.maximumStoredCells=count+4096;
    live.chargeWork=[&](std::uint64_t n){
        if(n>maxWork-out.work)throw Limit{};
        charge(n);out.work+=n;
    };
    const auto point=[&](Vertex id){return live.vertices[id].position;};
    const auto center=[&](const Cell& cell){
        Vec3 result{};
        for(const auto id:cell.vertices)result=result+point(id)*.25;
        return result;
    };
    const auto state=[&](Vec3 p){
        if(out.queries==maxQueries)throw Limit{};
        ++out.queries;
        return admit(p);
    };
    std::vector<unsigned char> accepted(admitted.begin(),admitted.end());
    std::vector<std::size_t> rejected;
    for(std::size_t i=0;i<count;++i)if(!accepted[i])rejected.push_back(i);
    out.initialRejected=out.remainingRejected=rejected.size();
    if(rejected.empty())return out;
    auto tetCount=count;
    try {
        live.vertices.reserve(mesh.nodeCount());
        for(const auto p:mesh.nodes()){
            live.charge();
            live.vertices.push_back({p,live.vertices.size()+1});
        }
        for(Index i=0;i<mesh.cellCount();++i){
            live.charge();
            const auto& cell=mesh.cell(i);
            const auto* v=mesh.cellNodes(i);
            if(cell.type==mesh::CellType::Tri3){
                const Complex::Face face{Vertex(v[0]),Vertex(v[1]),Vertex(v[2])};
                live.protectedFaces.insert(Complex::key(face));
                for(unsigned k=0;k<3;++k)
                    live.protectedEdges.insert(Complex::key(Complex::Edge{face[k],face[(k+1)%3]}));
            }else if(cell.type==mesh::CellType::Tet4){
                live.insert(Cell{{Vertex(v[0]),Vertex(v[1]),Vertex(v[2]),Vertex(v[3])},0,cell.onEntity});
            }else return out;
        }
        const auto orient=[&](const Complex::Face& f,Vec3 p){
            live.charge();
            return predicates::orient3d(point(f[0]),point(f[1]),point(f[2]),p);
        };
        const auto commit=[&](const Complex::Proposal& patch, bool seed){
            ++out.proposals;
            if(tetCount-patch.removed.size()+patch.added.size()>maxTets||
               patch.added.size()>live.maximumStoredCells-live.cells.size())return false;
            if(!live.preservesSurvivors(patch))return false;
            ++out.exactPatches;
            // Visit cells against the CAD shell first. They have the highest
            // chance of rejecting a proposal on an inward curved boundary.
            for(bool atBoundary:{true,false})for(const auto& cell:patch.added){
                live.charge();
                const auto faces=Complex::faceKeys(cell);
                const bool touches=std::ranges::any_of(faces,[&](const auto& key){
                    return live.protectedFaces.contains(key);
                });
                if(touches==atBoundary&&!state(center(cell))){
                    ++out.rejectedByCad;
                    return false;
                }
            }
            if(seed&&!state(point(Vertex(live.vertices.size()-1))))return false;
            live.apply(patch);
            for(const auto id:patch.removed)out.remainingRejected-=!accepted[id];
            accepted.resize(live.cells.size(),1);
            tetCount=tetCount-patch.removed.size()+patch.added.size();
            ++out.patched;
            return true;
        };
        const auto insert=[&](std::size_t bad,Vec3 seed,std::size_t patchLimit){
            Complex::Proposal patch;
            patch.removed.push_back(bad);
            for(std::size_t j=0;j<patch.removed.size();++j){
                const auto id=patch.removed[j];
                for(unsigned side=0;side<4;++side){
                    const auto neighbour=live.neighbours[id][side];
                    if(std::ranges::find(patch.removed,neighbour)!=patch.removed.end())continue;
                    const auto face=outwardFace(live.cells[id],side);
                    if(orient(face,seed)==Sign::Negative)continue;
                    if(neighbour==Complex::noCell||patch.removed.size()==patchLimit||
                       live.protectedFaces.contains(Complex::key(face)))return false;
                    patch.removed.push_back(neighbour);
                }
            }
            const auto staged=Vertex(live.vertices.size());
            for(const auto id:patch.removed)for(unsigned side=0;side<4;++side){
                if(std::ranges::find(patch.removed,live.neighbours[id][side])!=patch.removed.end())continue;
                const auto face=outwardFace(live.cells[id],side);
                patch.added.push_back({{staged,face[0],face[1],face[2]},0,live.cells[id].classification});
            }
            live.vertices.push_back({seed,std::uint64_t(staged)+1});
            try {
                if(live.valid(patch,std::nullopt,staged)&&commit(patch,true)){
                    ++out.seeds;
                    return true;
                }
            }catch(const Limit&){
                live.vertices.pop_back();
                throw;
            }
            live.vertices.pop_back();
            return false;
        };
        for(unsigned pass=0;pass<2;++pass){
            const auto old=out.remainingRejected;
            for(const auto id:rejected){
                live.charge();
                if(!live.alive[id])continue;
                const auto cell=live.cells[id];
                bool repaired=false;
                for(const auto edge:Complex::edgeKeys(cell)){
                    for(const auto& patch:live.edgeFlips(edge,std::nullopt))
                        if(commit(patch,false)){++out.flips;repaired=true;break;}
                    if(repaired)break;
                }
                if(repaired||live.vertices.size()>=maxNodes)continue;
                // Expand only across faces the apex cannot see. The seed
                // starts toward a nearby admitted cell and all resulting
                // faces must still close the same local region.
                std::vector<std::size_t> near{id};
                const auto origin=center(cell);
                for(std::size_t j=0;j<near.size()&&!repaired;++j){
                    for(const auto next:live.neighbours[near[j]]){
                        if(next==Complex::noCell||std::ranges::find(near,next)!=near.end())continue;
                        if(near.size()<16)near.push_back(next);
                        if(!accepted[next])continue;
                        const auto target=center(live.cells[next]);
                        for(double fraction:{.65,1.,.4}){
                            if(insert(id,origin+(target-origin)*fraction,maxPatch)){repaired=true;break;}
                        }
                        if(repaired)break;
                    }
                }
            }
            if(out.remainingRejected==0||out.remainingRejected==old)break;
        }
        // At a feature edge, each of the two protected triangles needs
        // inward clearance. Move along their common wedge, keeping both
        // triangles fixed and checking the entire replacement cavity.
        if(tryFeatureEdges&&out.remainingRejected>0&&out.remainingRejected<=4){
            unsigned trials=0;
            for(const auto id:rejected){
                live.charge();
                if(!live.alive[id]||accepted[id]||live.vertices.size()>=maxNodes)continue;
                const auto cell=live.cells[id];
                std::array<unsigned,4> exterior{};
                unsigned faces=0;
                for(unsigned side=0;side<4;++side)
                    if(live.protectedFaces.contains(Complex::key(outwardFace(cell,side))))
                        exterior[faces++]=side;
                if(faces!=2)continue;
                const auto inward=[&](unsigned side){
                    const auto f=outwardFace(cell,side);
                    return (point(f[2])-point(f[0])).cross(point(f[1])-point(f[0])).normalized();
                };
                const auto first=inward(exterior[0]),second=inward(exterior[1]);
                std::array<Vec3,2> shared{};
                unsigned endpoints=0;
                for(unsigned i=0;i<4;++i)
                    if(i!=exterior[0]&&i!=exterior[1])shared[endpoints++]=point(cell.vertices[i]);
                const auto edge=shared[1]-shared[0];
                const auto length=edge.norm();
                if(!(length>0.)||!std::isfinite(length))continue;
                const auto propose=[&](Vec3 origin){
                    for(double weight:{.5,.25,.75,.125,.875}){
                        const auto direction=(first*weight+second*(1.-weight)).normalized();
                        if(!(direction.dot(first)>0.&&direction.dot(second)>0.))continue;
                        for(double fraction:{.02,.01,.005,.04,.08,.12,.16,.24,.4}){
                            if(trials>=160)return false;
                            ++trials;
                            if(insert(id,origin+direction*(length*fraction),32))return true;
                        }
                    }
                    return false;
                };
                for(double along:{.5,.25,.75,.125,.875}){
                    if(trials>=160||propose(shared[0]+edge*along))break;
                }
            }
        }
        // Relocation is last because it changes coordinates recorded by
        // the spatial index used for connectivity patches above.
        std::map<Vertex,std::vector<std::size_t>> stars;
        for(const auto id:rejected){
            live.charge();
            if(!live.alive[id]||accepted[id])continue;
            for(const auto vertex:live.cells[id].vertices)
                if(vertex>=std::uint64_t(mesh.nodeCount())||mesh.nodeClass(Index(vertex)).dim==3)
                    stars.try_emplace(vertex);
        }
        for(std::size_t id=0;!stars.empty()&&id<live.cells.size();++id){
            live.charge();
            if(!live.alive[id])continue;
            for(const auto vertex:live.cells[id].vertices)
                if(auto found=stars.find(vertex);found!=stars.end()&&found->second.size()<=16)
                    found->second.push_back(id);
        }
        for(const auto id:rejected){
            live.charge();
            if(!live.alive[id]||accepted[id])continue;
            for(const auto vertex:live.cells[id].vertices){
                const auto found=stars.find(vertex);
                if(found==stars.end())continue;
                const auto& incident=found->second;
                if(const auto moved=relocateOcctLiveCavityApex(live,vertex,incident,id,state,live.chargeWork)){
                    live.vertices[vertex].position=*moved;
                    for(const auto owner:incident)if(!accepted[owner]){
                        --out.remainingRejected;accepted[owner]=1;
                    }
                    ++out.relocations;
                    ++out.patched;
                    break;
                }
            }
        }
    }catch(const Limit&){
        // Earlier exact transactions remain available within the budget.
    }
    if(!out.patched)return out;
    mesh::Mesh repaired;
    for(Index i=0;i<mesh.nodeCount();++i)repaired.addNode(point(Vertex(i)),mesh.nodeClass(i));
    const auto region=live.cells.front().classification;
    for(std::size_t i=mesh.nodeCount();i<live.vertices.size();++i)
        repaired.addNode(live.vertices[i].position,region);
    for(Index i=0;i<mesh.cellCount();++i)if(mesh.cell(i).type==mesh::CellType::Tri3)
        repaired.addCell(mesh::CellType::Tri3,std::span(mesh.cellNodes(i),3),mesh.cell(i).onEntity);
    for(std::size_t i=0;i<live.cells.size();++i)if(live.alive[i]){
        const auto& cell=live.cells[i];
        const auto& v=cell.vertices;
        repaired.addCell(mesh::CellType::Tet4,{Index(v[0]),Index(v[1]),Index(v[2]),Index(v[3])},cell.classification);
    }
    charge(0);
    mesh.swap(repaired);
    return out;
}
}

