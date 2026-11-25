// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/bnd/TetFacetPatchRemap.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <iostream>
#include <stdexcept>

using namespace paacs;
using namespace paacs::bnd::detail;
namespace {
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
std::array<Index,4> positive(const mesh::Mesh& mesh,std::array<Index,4> ids){
    if(predicates::orient3d(mesh.node(ids[0]),mesh.node(ids[1]),mesh.node(ids[2]),mesh.node(ids[3]))==predicates::Sign::Negative)
        std::swap(ids[1],ids[2]);
    require(predicates::orient3d(mesh.node(ids[0]),mesh.node(ids[1]),mesh.node(ids[2]),mesh.node(ids[3]))==predicates::Sign::Positive,"owned tet positive");
    return ids;
}
struct Split {mesh::Mesh boundary;FacetSourcePatch patch;};
Split split(const mesh::Mesh& old,Index selected,Vec3 point){
    Split result;for(Index i=0;i<old.nodeCount();++i)result.boundary.addNode(old.node(i),old.nodeClass(i));
    const auto entity=old.cell(selected).onEntity;const auto inserted=result.boundary.addNode(point,entity);
    const auto* source=old.cellNodes(selected);
    result.patch={1,entity.tag,selected,{source[0],source[1],source[2]},inserted};
    for(Index i=0;i<old.cellCount();++i){const auto* p=old.cellNodes(i);
        if(i!=selected)result.boundary.addCell(mesh::CellType::Tri3,{p[0],p[1],p[2]},old.cell(i).onEntity);
        else for(unsigned j=0;j<3;++j)result.boundary.addCell(mesh::CellType::Tri3,{p[j],p[(j+1)%3],inserted},entity);
    }
    return result;
}
Index triangle(const mesh::Mesh& source,std::array<Index,3> ids){
    std::sort(ids.begin(),ids.end());
    for(Index i=0;i<source.cellCount();++i){const auto* p=source.cellNodes(i);std::array actual{p[0],p[1],p[2]};std::sort(actual.begin(),actual.end());if(actual==ids)return i;}
    throw std::runtime_error("owned source triangle absent");
}
void exercise(){
    mesh::Mesh source;for(const auto p:std::array{Vec3{0,0,0},Vec3{1,0,0},Vec3{0,1,0},Vec3{0,0,1}})
        source.addNode(p,{0,int(source.nodeCount()+1)});
    const auto parent=positive(source,{0,1,2,3});paacs::bnd::TetFacetCell tet;
    for(unsigned j=0;j<4;++j)tet.vertices[j]=static_cast<paacs::bnd::TetFacetVertexId>(parent[j]);
    for(unsigned side=0;side<4;++side){const auto face=paacs::bnd::tet_facet_detail::outwardFace(tet,side);
        source.addCell(mesh::CellType::Tri3,{face[0],face[1],face[2]},{2,int(side)+1});}
    mesh::Mesh partial;for(Index i=0;i<source.nodeCount();++i)partial.addNode(source.node(i),source.nodeClass(i));
    const auto apex=partial.addNode({.2,.2,.2},{3,1});
    partial.addCell(mesh::CellType::Tet4,positive(partial,{0,1,2,apex}),{3,1});
    // Same AABB as the swept simplex, but beyond its x+y supporting plane.
    const auto near=partial.nodeCount();
    for(const auto p:std::array{Vec3{.85,.10,.02},Vec3{.86,.10,.02},Vec3{.85,.11,.02},Vec3{.85,.10,.03}})partial.addNode(p,{3,1});
    partial.addCell(mesh::CellType::Tet4,positive(partial,{near,near+1,near+2,near+3}),{3,1});
    const auto far=partial.nodeCount();
    for(const auto p:std::array{Vec3{.1,.1,.5},Vec3{.2,.1,.5},Vec3{.1,.2,.5},Vec3{.1,.1,.6}})partial.addNode(p,{3,1});
    partial.addCell(mesh::CellType::Tet4,positive(partial,{far,far+1,far+2,far+3}),{3,1});
    const auto originalHash=mesh::meshHash(partial);
    auto edited=split(source,triangle(source,{0,1,2}),{.25,.25,.1});
    const std::array changes{edited.patch};FacetPatchRemapLimits limits;mesh::Mesh remapped;
    const auto result=remapFacetPatchPartial(edited.boundary,partial,source.nodeCount(),changes,{3,1},limits,remapped);
    require(result.prepared&&result.removedInputCells==std::vector<Index>{0}&&result.retainedCells==2,"only changed tet removed");
    require(result.narrowPairs>=2,"near disjoint tet reaches exact narrow phase");
    require(remapped.nodeCount()==partial.nodeCount()+1&&remapped.cellCount()==2,"partial populations");
    for(Index i=0;i<partial.nodeCount();++i){const auto mapped=i<source.nodeCount()?i:i+1;
        require(facet_patch_remap_detail::same(partial.node(i),remapped.node(mapped))&&partial.nodeClass(i)==remapped.nodeClass(mapped),"coordinates and classifications preserved");}
    for(Index i=0;i<2;++i){const auto* a=partial.cellNodes(i+1);const auto* b=remapped.cellNodes(i);
        for(unsigned j=0;j<4;++j)require(b[j]==a[j]+1,"retained connectivity remapped in order");}
    // A second committed patch preserves the previous inserted source node,
    // then shifts each old interior index once more; no whole-volume front.
    auto second=split(edited.boundary,triangle(edited.boundary,{0,1,4}),{.1,.05,.03});
    const std::array secondChange{second.patch};mesh::Mesh twice;
    const auto again=remapFacetPatchPartial(second.boundary,remapped,edited.boundary.nodeCount(),secondChange,{3,1},limits,twice);
    require(again.prepared&&again.removedCells==0&&again.retainedCells==2,"second local patch retains both disjoint tets");
    for(Index i=0;i<remapped.nodeCount();++i){const auto mapped=i<edited.boundary.nodeCount()?i:i+1;
        require(facet_patch_remap_detail::same(remapped.node(i),twice.node(mapped))&&remapped.nodeClass(i)==twice.nodeClass(mapped),"two-round prefix and interior coordinates");}
    auto planar=split(source,triangle(source,{0,1,2}),{.2,.2,0});const std::array planarChange{planar.patch};mesh::Mesh degenerate;
    const auto flat=remapFacetPatchPartial(planar.boundary,partial,source.nodeCount(),planarChange,{3,1},limits,degenerate);
    require(flat.prepared&&flat.degenerateSweeps==1&&flat.removedInputCells==std::vector<Index>{0},"zero-volume sweep uses conservative closed AABB");
    // An outward sweep shares the original boundary face legally. It expands
    // the domain, so the unchanged old tetrahedron must survive the remap.
    mesh::Mesh unit;for(Index i=0;i<source.nodeCount();++i)unit.addNode(source.node(i),source.nodeClass(i));
    unit.addCell(mesh::CellType::Tet4,parent,{3,1});
    auto outward=split(source,triangle(source,{0,1,2}),{.25,.25,-.1});
    const std::array outwardChange{outward.patch};mesh::Mesh expanded;
    const auto legal=remapFacetPatchPartial(outward.boundary,unit,source.nodeCount(),outwardChange,{3,1},limits,expanded);
    require(legal.prepared&&legal.narrowPairs==1&&legal.removedCells==0&&expanded.cellCount()==1,
            "legal old ABC shared-face contact retains the old tet");
    // P has a new identity. A geometric contact on an existing tet face or
    // edge is hanging even when the swept simplex has no volume overlap.
    for(const auto points:std::array{
            std::array{Vec3{.2,.2,.1},Vec3{.4,.2,.1},Vec3{.2,.4,.1},Vec3{.2,.2,.2}},
            std::array{Vec3{.2,.25,.1},Vec3{.3,.25,.1},Vec3{.2,.35,.2},Vec3{.2,.25,.2}}}){
        mesh::Mesh contact;for(Index i=0;i<source.nodeCount();++i)contact.addNode(source.node(i),source.nodeClass(i));
        for(const auto p:points)contact.addNode(p,{3,1});
        contact.addCell(mesh::CellType::Tet4,positive(contact,{4,5,6,7}),{3,1});
        mesh::Mesh clipped;
        const auto result=remapFacetPatchPartial(edited.boundary,contact,source.nodeCount(),changes,{3,1},limits,clipped);
        require(result.prepared&&result.degenerateSweeps==0&&result.narrowPairs==1&&
                result.removedInputCells==std::vector<Index>{0}&&clipped.cellCount()==0,
                "new source P touching an old tet face or edge removes the hanging contact");
    }
    // Classification errors are rejected before a partial candidate can be
    // published; the helper does not infer or rewrite CAD ownership.
    for(const auto wrong:std::array{mesh::Classification{3,1},mesh::Classification{2,-1}}){
        mesh::Mesh invalid;
        for(Index i=0;i<edited.boundary.nodeCount();++i)
            invalid.addNode(edited.boundary.node(i),i==0?wrong:edited.boundary.nodeClass(i));
        for(Index i=0;i<edited.boundary.cellCount();++i)
            invalid.addCell(mesh::CellType::Tri3,std::span<const Index>(edited.boundary.cellNodes(i),3),edited.boundary.cell(i).onEntity);
        mesh::Mesh sentinel;sentinel.addNode({9,8,7},{3,99});const auto hash=mesh::meshHash(sentinel);
        const auto result=remapFacetPatchPartial(invalid,partial,source.nodeCount(),changes,{3,1},limits,sentinel);
        require(!result.prepared&&result.reason==FacetPatchRemapReason::InvalidInput&&mesh::meshHash(sentinel)==hash,
                "invalid source classification refuses atomically");
    }
    const auto refusal=[&](FacetPatchRemapLimits bound,FacetPatchRemapReason reason){
        mesh::Mesh sentinel;sentinel.addNode({9,8,7},{3,99});const auto hash=mesh::meshHash(sentinel);
        const auto stopped=remapFacetPatchPartial(edited.boundary,partial,source.nodeCount(),changes,{3,1},bound,sentinel);
        require(!stopped.prepared&&stopped.reason==reason&&mesh::meshHash(sentinel)==hash&&mesh::meshHash(partial)==originalHash,"atomic bounded refusal");
    };
    auto bound=limits;bound.maximumWorkingBytes=1;refusal(bound,FacetPatchRemapReason::BudgetExceeded);
    bound=limits;bound.maximumWork=1;refusal(bound,FacetPatchRemapReason::BudgetExceeded);
    bound=limits;bound.maximumRemovedCells=0;refusal(bound,FacetPatchRemapReason::BudgetExceeded);
    bound=limits;bound.maximumNodes=partial.nodeCount();refusal(bound,FacetPatchRemapReason::BudgetExceeded);
    std::atomic<bool> cancel{true};bound=limits;bound.cancellation=&cancel;refusal(bound,FacetPatchRemapReason::Cancelled);
    require(mesh::meshHash(partial)==originalHash,"original partial unchanged after both successful rounds");
}
}
int main(){try{exercise();std::cout<<"local patch remap checks passed\n";return 0;}
catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
