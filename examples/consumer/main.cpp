// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <paacs/sdk/MeshArchive.hpp>
#include <paacs/sdk/Version.hpp>
#include <array>
#include <iostream>
int main() {
    paacs::mesh::Mesh m;
    (void)m.addNode({0,0,0},{0,1});(void)m.addNode({1,0,0},{0,2});
    const std::array<paacs::Index,2> conn{0,1};
    (void)m.addCell(paacs::mesh::CellType::Line2,std::span<const paacs::Index>(conn),{1,1});
    std::vector<std::uint8_t> bytes;
    if(!paacs::sdk::serializeMeshArchive(m,bytes).completed)return 1;
    paacs::mesh::Mesh read;
    if(!paacs::sdk::deserializeMeshArchive(bytes,read).completed || read.cellCount()!=1)return 1;
    std::cout<<paacs::sdk::kApiVersion<<" C++ consumer passed\n";
}
