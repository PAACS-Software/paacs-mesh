// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <paacs/geometry/OcctImport.hpp>
#include <paacs/algo/GeometryMesher.hpp>
#include <iostream>
int main(int argc,char **argv) {
    if(argc!=2)return 2;
    auto imported=paacs::geo::importPlanarCad(argv[1]);
    if(!imported.imported())return 1;
    paacs::algo::GeometryMeshingOptions options;options.targetSize=.5;
    paacs::mesh::Mesh mesh;
    const auto result=paacs::algo::meshGeometryModel(*imported.model,options,mesh);
    if(!result.published || !result.audit.valid)return 1;
    std::cout<<"Installed CAD consumer passed: "<<mesh.nodeCount()<<" nodes\n";
}
