// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <paacs/geometry/OcctMesher.hpp>
#include <paacs/sdk/MeshArchive.hpp>
#include <paacs/mesh/MeshAudit.hpp>
#include <iostream>

int main(int argc,char** argv) {
    if(argc!=2)return 2;
    const auto imported=paacs::geo::importOcctCad(argv[1]);
    if(!imported.imported()||imported.curvedFaces!=1||imported.solids!=1)return 1;
    paacs::geo::NativeCadMeshOptions options;options.geometry.targetSize=.3;
    paacs::mesh::Mesh mesh;
    const auto result=paacs::geo::meshNativeCad(*imported.model,options,mesh);
    if(!result.published||!result.audit.valid||!result.quality.valid)return 1;
    const auto written=paacs::sdk::writeMeshArchive(mesh,"native-consumer.paacsmesh");
    paacs::mesh::Mesh readback;
    if(!written.completed||!paacs::sdk::readMeshArchive("native-consumer.paacsmesh",readback).completed||
       paacs::mesh::meshHash(mesh)!=paacs::mesh::meshHash(readback))return 1;
    const auto audit=paacs::mesh::auditLinearVolumeBoundary(readback);
    if(!audit.valid)return 1;
    std::cout<<"Installed native CAD and archive consumer passed: "<<readback.nodeCount()<<" nodes\n";
    return 0;
}
