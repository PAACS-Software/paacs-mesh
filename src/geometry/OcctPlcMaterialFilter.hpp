// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private adapter step: remove every incident cell before deleting an interior
// point rejected by native material. The original boundary prefix is immutable.
#include "paacs/mesh/Mesh.hpp"
#include <array>
#include <stdexcept>
#include <vector>
namespace paacs::geo::detail {
struct OcctPlcMaterialFilterStats {
    std::uint64_t removedInteriorNodes=0,removedIncidentCells=0,removedCentroidCells=0;
};
template<class AdmitNode,class AdmitCentroid,class Charge>
OcctPlcMaterialFilterStats retainOcctPlcMaterial(const mesh::Mesh& input,Index boundaryNodes,
    AdmitNode&& admitNode,AdmitCentroid&& admitCentroid,Charge&& charge,mesh::Mesh& output){
    OcctPlcMaterialFilterStats stats;mesh::Mesh retained;
    if(&input==&output||boundaryNodes<0||boundaryNodes>input.nodeCount())throw std::invalid_argument("native PLC material filter prefix/alias");
    std::vector<Index> map(input.nodeCount(),kInvalidIndex);
    for(Index i=0;i<input.nodeCount();++i){
        charge(1);if(i>=boundaryNodes&&!admitNode(input.node(i))){++stats.removedInteriorNodes;continue;}
        map[i]=retained.addNode(input.node(i),input.nodeClass(i));
    }
    for(Index i=0;i<input.cellCount();++i){
        charge(1);const auto& cell=input.cell(i);if(cell.type==mesh::CellType::Tri3)continue;
        if(cell.type!=mesh::CellType::Tet4)throw std::invalid_argument("native PLC material filter cell type");
        const auto* ids=input.cellNodes(i);std::array<Index,4> mapped;Vec3 center{};bool incident=false;
        for(unsigned j=0;j<4;++j){if(ids[j]<0||ids[j]>=input.nodeCount())throw std::invalid_argument("native PLC material filter cell index");mapped[j]=map[ids[j]];incident|=mapped[j]==kInvalidIndex;center=center+input.node(ids[j])*.25;}
        if(incident){++stats.removedIncidentCells;continue;}
        if(!admitCentroid(center)){++stats.removedCentroidCells;continue;}
        retained.addCell(mesh::CellType::Tet4,mapped,cell.onEntity);
    }
    charge(0);output.swap(retained);return stats;
}
template<class Admit,class Charge>
OcctPlcMaterialFilterStats retainOcctPlcMaterial(const mesh::Mesh& input,Index boundaryNodes,
    Admit&& admit,Charge&& charge,mesh::Mesh& output){
    return retainOcctPlcMaterial(input,boundaryNodes,admit,admit,charge,output);
}
}
