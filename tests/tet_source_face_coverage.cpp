// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/bnd/TetSourceFaceCoverage.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <iostream>
#include <stdexcept>
#include <string_view>

using namespace paacs;
using namespace paacs::bnd;
using namespace paacs::bnd::detail;
namespace {
unsigned checks=0;
void require(bool value,const char* why){++checks;if(!value)throw std::runtime_error(why);}
std::array<Index,4> positive(const mesh::Mesh& mesh,std::array<Index,4> ids){
    if(predicates::orient3d(mesh.node(ids[0]),mesh.node(ids[1]),mesh.node(ids[2]),mesh.node(ids[3]))==predicates::Sign::Negative)
        std::swap(ids[1],ids[2]);
    require(predicates::orient3d(mesh.node(ids[0]),mesh.node(ids[1]),mesh.node(ids[2]),mesh.node(ids[3]))==predicates::Sign::Positive,"owned fixture tet is positive");
    return ids;
}
void appendSurface(mesh::Mesh& mesh,std::array<Index,4> ids){
    TetFacetCell cell;for(unsigned j=0;j<4;++j)cell.vertices[j]=static_cast<TetFacetVertexId>(ids[j]);
    for(unsigned side=0;side<4;++side){const auto f=tet_facet_detail::outwardFace(cell,side);
        mesh.addCell(mesh::CellType::Tri3,{f[0],f[1],f[2]},{2,int(mesh.cellCount())+1});}
}
mesh::Mesh changedSurface(const mesh::Mesh& old,Index reversed,bool duplicate=false){
    mesh::Mesh result;for(Index i=0;i<old.nodeCount();++i)result.addNode(old.node(i),old.nodeClass(i));
    for(Index ci=0;ci<old.cellCount();++ci){const auto* p=old.cellNodes(ci);std::array ids{p[0],p[1],p[2]};
        if(ci==reversed)std::swap(ids[1],ids[2]);result.addCell(mesh::CellType::Tri3,ids,old.cell(ci).onEntity);}
    if(duplicate)result.addCell(mesh::CellType::Tri3,std::span<const Index>(old.cellNodes(0),3),old.cell(0).onEntity);
    return result;
}
void exercise(){
    mesh::Mesh source;
    for(auto p:std::array{Vec3{0,0,0},Vec3{1,0,0},Vec3{0,1,0},Vec3{0,0,1}})
        source.addNode(p,{0,int(source.nodeCount())+1});
    const auto parent=positive(source,{0,1,2,3});appendSurface(source,parent);
    mesh::Mesh full;for(Index i=0;i<source.nodeCount();++i)full.addNode(source.node(i),source.nodeClass(i));
    full.addCell(mesh::CellType::Tet4,parent,{3,1});
    SourceFaceCoverageLimits limits;
    const auto complete=collectMissingSourceTriangles(source,full,{3,1},limits);
    require(complete.complete&&complete.sourceFacesCovered==4&&complete.tetFacesVisited==4&&complete.missingSourceTriangles.empty(),"all four directed source faces covered");
    const auto repeated=collectMissingSourceTriangles(source,full,{3,1},limits);
    require(repeated.complete&&repeated.workUsed==complete.workUsed,"coverage work is deterministic");
    mesh::Mesh partial;for(Index i=0;i<source.nodeCount();++i)partial.addNode(source.node(i),source.nodeClass(i));
    const auto apex=partial.addNode({.2,.2,.2},{3,1});
    partial.addCell(mesh::CellType::Tet4,positive(partial,{0,1,2,apex}),{3,1});
    const auto some=collectMissingSourceTriangles(source,partial,{3,1},limits);
    require(some.complete&&some.sourceFacesCovered==1&&some.missingSourceTriangles.size()==3&&
            std::is_sorted(some.missingSourceTriangles.begin(),some.missingSourceTriangles.end()),"missing IDs belong to original source and are sorted");
    const auto sourceHash=mesh::meshHash(source),partialHash=mesh::meshHash(partial);
    const auto refused=[&](const mesh::Mesh& boundary,const mesh::Mesh& mesh,SourceFaceCoverageLimits cap,
                           SourceFaceCoverageReason reason,std::string_view diagnostic={}){
        const auto result=collectMissingSourceTriangles(boundary,mesh,{3,1},cap);
        require(!result.complete&&result.reason==reason&&result.missingSourceTriangles.empty(),"refused diagnosis publishes no partial missing-face receipt");
        if(!diagnostic.empty())require(result.diagnostic==diagnostic,"precise structural refusal diagnostic");
    };
    const auto reversed=changedSurface(source,0);
    refused(reversed,full,limits,SourceFaceCoverageReason::InvalidInput,"source-coverage-source-winding");
    const auto duplicated=changedSurface(source,-1,true);
    refused(duplicated,full,limits,SourceFaceCoverageReason::InvalidInput,"source-coverage-duplicate-source-face");
    auto doubleOwner=full;doubleOwner.addCell(mesh::CellType::Tet4,parent,{3,1});
    refused(source,doubleOwner,limits,SourceFaceCoverageReason::InvalidInput,"source-coverage-multiple-source-owners");
    auto changed=full;changed.setNode(0,{.001,0,0});
    refused(source,changed,limits,SourceFaceCoverageReason::InvalidInput,"source-coverage-prefix");
    auto cap=limits;cap.maximumWorkingBytes=1;refused(source,partial,cap,SourceFaceCoverageReason::BudgetExceeded);
    cap=limits;cap.maximumWork=1;refused(source,partial,cap,SourceFaceCoverageReason::BudgetExceeded);
    cap=limits;cap.maximumSourceTriangles=3;refused(source,partial,cap,SourceFaceCoverageReason::BudgetExceeded);
    cap=limits;cap.maximumPartialCells=0;refused(source,partial,cap,SourceFaceCoverageReason::BudgetExceeded);
    std::atomic<bool> cancelled{true};cap=limits;cap.cancellation=&cancelled;
    refused(source,partial,cap,SourceFaceCoverageReason::Cancelled);
    require(mesh::meshHash(source)==sourceHash&&mesh::meshHash(partial)==partialHash,"all diagnostics preserve source and partial bytes");
    // Several independent components put genuinely missing source IDs beyond
    // the first 32 entries; no file name, CAD shape recognition or Native call.
    mesh::Mesh many;
    std::vector<std::array<Index,4>> parents;
    for(unsigned component=0;component<12;++component){const auto begin=many.nodeCount();
        for(auto p:std::array{Vec3{0,0,0},Vec3{1,0,0},Vec3{0,1,0},Vec3{0,0,1}}){p.x+=component*2.;many.addNode(p,{0,int(many.nodeCount())+1});}
        parents.push_back(positive(many,{begin,begin+1,begin+2,begin+3}));appendSurface(many,parents.back());
    }
    mesh::Mesh eleven;for(Index i=0;i<many.nodeCount();++i)eleven.addNode(many.node(i),many.nodeClass(i));
    for(unsigned i=0;i<11;++i)eleven.addCell(mesh::CellType::Tet4,parents[i],{3,1});
    const auto later=collectMissingSourceTriangles(many,eleven,{3,1},limits);
    require(later.complete&&later.sourceFacesCovered==44&&later.missingSourceTriangles==std::vector<Index>({44,45,46,47}),"all components inspected without a first-32 raw-ID cutoff");
}
}
int main(){try{exercise();std::cout<<checks<<" source face coverage checks passed\n";return 0;}
catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
