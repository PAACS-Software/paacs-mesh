// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/bnd/TetFacetGeometry.hpp"
#include "TetFacetSourcePatch.hpp"
#include <atomic>
#include <bit>
#include <new>

namespace paacs::bnd::detail {
struct FacetPatchRemapLimits {
    std::uint64_t maximumWork=10'000'000, maximumWorkingBytes=256ULL*1024*1024;
    std::uint64_t maximumNodes=1'000'000, maximumCells=5'000'000;
    std::uint64_t maximumPatches=128,maximumRemovedCells=500'000;
    const std::atomic<bool>* cancellation=nullptr;
};
enum class FacetPatchRemapReason {None,InvalidInput,BudgetExceeded,Cancelled,AllocationFailed,InternalError};
struct FacetPatchRemapOutcome {
    bool prepared=false;
    FacetPatchRemapReason reason=FacetPatchRemapReason::InternalError;
    const char* diagnostic="not_attempted";
    std::uint64_t workUsed=0,modeledBytes=0,broadPairs=0,narrowPairs=0;
    std::uint64_t inputCells=0,retainedCells=0,removedCells=0,degenerateSweeps=0;
    // Input cell indices, available only after successful preparation. These
    // are not stable final-mesh element IDs or proof of volume completion.
    std::vector<Index> removedInputCells;
};
namespace facet_patch_remap_detail {
struct Stop {FacetPatchRemapReason reason;const char* diagnostic;};
struct Box {Vec3 lower{},upper{};};
struct Sweep {Box bounds{};TetFacetCell cell{};bool degenerate=false;};
inline bool same(Vec3 a,Vec3 b){
    return std::bit_cast<std::uint64_t>(a.x)==std::bit_cast<std::uint64_t>(b.x)&&
        std::bit_cast<std::uint64_t>(a.y)==std::bit_cast<std::uint64_t>(b.y)&&
        std::bit_cast<std::uint64_t>(a.z)==std::bit_cast<std::uint64_t>(b.z);
}
inline Box box(std::span<const TetFacetVertex> points,const std::array<TetFacetVertexId,4>& ids){
    Box result{points[ids[0]].position,points[ids[0]].position};
    for(auto id:ids){const auto p=points[id].position;
        result.lower={std::min(result.lower.x,p.x),std::min(result.lower.y,p.y),std::min(result.lower.z,p.z)};
        result.upper={std::max(result.upper.x,p.x),std::max(result.upper.y,p.y),std::max(result.upper.z,p.z)};
    }
    return result;
}
inline bool overlaps(Box a,Box b){
    return a.lower.x<=b.upper.x&&b.lower.x<=a.upper.x&&
        a.lower.y<=b.upper.y&&b.lower.y<=a.upper.y&&
        a.lower.z<=b.upper.z&&b.lower.z<=a.upper.z;
}
}

// Trusted provenance contract: patches come only from a successfully committed
// facet-only ABC -> ABP/BCP/CAP edit, preserving every old source node. The caller
// has independently certified the new embedded surface. This prepares a partial
// candidate ONLY. It never replaces recoverTetPlcCavitiesRetaining's complete
// source/contact/centroid/pair checks or any final Native validity gate.
inline FacetPatchRemapOutcome remapFacetPatchPartial(const mesh::Mesh& boundary,
    const mesh::Mesh& partial,Index previousBoundaryNodes,
    std::span<const FacetSourcePatch> patches,mesh::Classification region,
    const FacetPatchRemapLimits& limits,mesh::Mesh& destination) noexcept {
    using Sweep=facet_patch_remap_detail::Sweep;
    using RemapStop=facet_patch_remap_detail::Stop;
    FacetPatchRemapOutcome out;
    const auto need=[](bool ok,FacetPatchRemapReason why,const char* diagnostic){if(!ok)throw RemapStop{why,diagnostic};};
    const auto tick=[&](std::uint64_t n=1){
        if(limits.cancellation&&limits.cancellation->load())throw RemapStop{FacetPatchRemapReason::Cancelled,"cancelled"};
        need(out.workUsed<=limits.maximumWork&&n<=limits.maximumWork-out.workUsed,
             FacetPatchRemapReason::BudgetExceeded,"facet-remap-work-budget");out.workUsed+=n;
    };
    try {
        tick(0);
        need(&boundary!=&destination && &partial!=&destination&&previousBoundaryNodes>=0&&
             previousBoundaryNodes<=boundary.nodeCount()&&previousBoundaryNodes<=partial.nodeCount()&&region.dim==3&&region.tag>=0,
             FacetPatchRemapReason::InvalidInput,"facet-remap-prefix-or-alias");
        const auto sourceNodes=std::uint64_t(boundary.nodeCount());
        const auto added=sourceNodes-std::uint64_t(previousBoundaryNodes);
        const auto interiors=std::uint64_t(partial.nodeCount()-previousBoundaryNodes);
        need(sourceNodes<=UINT64_MAX-interiors,FacetPatchRemapReason::BudgetExceeded,"facet-remap-node-overflow");
        const auto nodes=sourceNodes+interiors,cells=std::uint64_t(partial.cellCount());out.inputCells=cells;
        need(nodes<kInvalidTetFacetVertex&&nodes<=limits.maximumNodes&&patches.size()<=limits.maximumPatches,
             FacetPatchRemapReason::BudgetExceeded,"facet-remap-count-budget");
        need(added==patches.size(),FacetPatchRemapReason::InvalidInput,"facet-remap-patch-coverage");
        const auto addBytes=[&](std::uint64_t count,std::uint64_t unit){
            need(!unit||count<=(UINT64_MAX-out.modeledBytes)/unit,
                 FacetPatchRemapReason::BudgetExceeded,"facet-remap-byte-overflow");
            out.modeledBytes+=count*unit;
        };
        addBytes(nodes,sizeof(TetFacetVertex)+sizeof(Vec3)+sizeof(mesh::Classification));
        addBytes(std::min(cells,limits.maximumCells),sizeof(mesh::Cell)+4*sizeof(Index));
        addBytes(cells,sizeof(Index));addBytes(patches.size(),sizeof(Sweep)+1);addBytes(1,16*1024);
        need(out.modeledBytes<=limits.maximumWorkingBytes,FacetPatchRemapReason::BudgetExceeded,"facet-remap-working-budget");
        std::vector<TetFacetVertex> vertices;vertices.reserve(static_cast<std::size_t>(nodes));
        for(Index i=0;i<boundary.nodeCount();++i){tick();const auto p=boundary.node(i);
            need(tet_facet_detail::finite(p)&&boundary.nodeClass(i).dim>=0&&boundary.nodeClass(i).dim<=2&&boundary.nodeClass(i).tag>=0,
                 FacetPatchRemapReason::InvalidInput,"facet-remap-source-node");
            if(i<previousBoundaryNodes)need(facet_patch_remap_detail::same(p,partial.node(i))&&boundary.nodeClass(i)==partial.nodeClass(i),
                FacetPatchRemapReason::InvalidInput,"facet-remap-prefix-changed");
            vertices.push_back({p,vertices.size()+1});
        }
        for(Index i=previousBoundaryNodes;i<partial.nodeCount();++i){tick();
            need(tet_facet_detail::finite(partial.node(i))&&partial.nodeClass(i)==region,
                 FacetPatchRemapReason::InvalidInput,"facet-remap-interior-node");
            vertices.push_back({partial.node(i),vertices.size()+1});
        }
        std::vector<unsigned char> seen(static_cast<std::size_t>(added),0);
        std::vector<Sweep> sweeps;sweeps.reserve(patches.size());
        for(const auto& patch:patches){tick();
            need(patch.bodyTag==region.tag&&patch.faceTag>=0&&patch.originalTriangle>=0&&patch.insertedVertex>=previousBoundaryNodes&&
                 patch.insertedVertex<boundary.nodeCount()&&boundary.nodeClass(patch.insertedVertex)==mesh::Classification{2,patch.faceTag},
                 FacetPatchRemapReason::InvalidInput,"facet-remap-patch-identity");
            const auto inserted=static_cast<std::size_t>(patch.insertedVertex-previousBoundaryNodes);
            need(!seen[inserted]++,FacetPatchRemapReason::InvalidInput,"facet-remap-duplicate-patch");
            Sweep sweep;
            for(unsigned j=0;j<3;++j){const auto id=patch.originalVertices[j];
                need(id>=0&&id<previousBoundaryNodes,FacetPatchRemapReason::InvalidInput,"facet-remap-patch-index");
                sweep.cell.vertices[j]=static_cast<TetFacetVertexId>(id);
            }
            need(patch.originalVertices[0]!=patch.originalVertices[1]&&patch.originalVertices[0]!=patch.originalVertices[2]&&
                 patch.originalVertices[1]!=patch.originalVertices[2],FacetPatchRemapReason::InvalidInput,"facet-remap-repeated-source-vertex");
            sweep.cell.vertices[3]=static_cast<TetFacetVertexId>(patch.insertedVertex);
            const auto& p=sweep.cell.vertices;
            const auto sign=predicates::orient3d(vertices[p[0]].position,vertices[p[1]].position,vertices[p[2]].position,vertices[p[3]].position);
            sweep.degenerate=sign==predicates::Sign::Zero;out.degenerateSweeps+=sweep.degenerate;
            if(sign==predicates::Sign::Negative)std::swap(sweep.cell.vertices[1],sweep.cell.vertices[2]);
            sweep.bounds=facet_patch_remap_detail::box(vertices,sweep.cell.vertices);sweeps.push_back(sweep);
        }
        need(std::find(seen.begin(),seen.end(),0)==seen.end(),FacetPatchRemapReason::InvalidInput,"facet-remap-unrecorded-new-node");
        const auto mapped=[&](Index id){
            need(id>=0&&id<partial.nodeCount(),FacetPatchRemapReason::InvalidInput,"facet-remap-tet-index");
            return static_cast<TetFacetVertexId>(std::uint64_t(id)+(id>=previousBoundaryNodes?added:0));
        };
        const TetFacetComplexView view{vertices,{},{},{},{},region};
        std::vector<Index> removed;removed.reserve(static_cast<std::size_t>(std::min(cells,limits.maximumRemovedCells)));
        for(Index ci=0;ci<partial.cellCount();++ci){tick();
            need(partial.cell(ci).type==mesh::CellType::Tet4&&partial.cell(ci).onEntity==region,
                 FacetPatchRemapReason::InvalidInput,"facet-remap-tet-region");
            TetFacetCell cell;const auto* source=partial.cellNodes(ci);
            for(unsigned j=0;j<4;++j)cell.vertices[j]=mapped(source[j]);
            const auto& p=cell.vertices;
            need(predicates::orient3d(vertices[p[0]].position,vertices[p[1]].position,vertices[p[2]].position,vertices[p[3]].position)==predicates::Sign::Positive,
                 FacetPatchRemapReason::InvalidInput,"facet-remap-nonpositive-input");
            const auto bounds=facet_patch_remap_detail::box(vertices,cell.vertices);bool drop=false;
            for(const auto& sweep:sweeps){tick();++out.broadPairs;if(!facet_patch_remap_detail::overlaps(bounds,sweep.bounds))continue;
                if(sweep.degenerate){drop=true;break;}
                tick();++out.narrowPairs;
                if(tet_facet_detail::tetImproperlyIntersects(view,cell,sweep.cell)){drop=true;break;}
            }
            if(drop){need(removed.size()<limits.maximumRemovedCells,FacetPatchRemapReason::BudgetExceeded,"facet-remap-removal-budget");removed.push_back(ci);}
        }
        out.removedCells=removed.size();out.retainedCells=cells-out.removedCells;
        need(out.retainedCells<=limits.maximumCells,FacetPatchRemapReason::BudgetExceeded,"facet-remap-cell-budget");
        mesh::Mesh staged;staged.reserveAppend(static_cast<std::size_t>(nodes),static_cast<std::size_t>(out.retainedCells),static_cast<std::size_t>(out.retainedCells*4));
        for(Index i=0;i<boundary.nodeCount();++i){tick();staged.addNode(boundary.node(i),boundary.nodeClass(i));}
        for(Index i=previousBoundaryNodes;i<partial.nodeCount();++i){tick();staged.addNode(partial.node(i),partial.nodeClass(i));}
        std::size_t nextRemoved=0;
        for(Index ci=0;ci<partial.cellCount();++ci){tick();
            if(nextRemoved<removed.size()&&removed[nextRemoved]==ci){++nextRemoved;continue;}
            std::array<Index,4> ids;const auto* source=partial.cellNodes(ci);
            for(unsigned j=0;j<4;++j)ids[j]=mapped(source[j]);
            staged.addCell(mesh::CellType::Tet4,ids,partial.cell(ci).onEntity);
        }
        tick(0);destination.swap(staged);out.removedInputCells.swap(removed);
        out.prepared=true;out.reason=FacetPatchRemapReason::None;out.diagnostic="partial_remap_prepared";
    }catch(const RemapStop& e){out.reason=e.reason;out.diagnostic=e.diagnostic;}
    catch(const std::bad_alloc&){out.reason=FacetPatchRemapReason::AllocationFailed;out.diagnostic="allocation-failed";}
    catch(...){out.reason=FacetPatchRemapReason::InternalError;out.diagnostic="internal-error";}
    return out;
}
}
