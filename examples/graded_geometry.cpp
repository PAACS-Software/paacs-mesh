// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/geometry/PolyhedralModel.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <iostream>

int main() {
    using namespace paacs;
    const std::vector<std::vector<Vec2>> profile{
        {{0,0},{3,0},{3,2},{0,2}}, {{.1,.1},{.1,1.9},{2.9,1.9},{2.9,.1}}};
    auto geometry=geo::makePolyhedralModel(geo::extrudedPolygon(profile,{.3,0,2}));
    if(!geometry.built()) return 1;
    algo::GeometryMeshingOptions options;
    options.targetSize=.25;
    options.maximumOutputBytes=64ULL*1024*1024;
    options.maximumAuditTemporaryBytes=32ULL*1024*1024;
    algo::GeometryBodyControl control;
    control.solidTag=geometry.model->solids().front()->tag();
    control.layers=12;
    control.endToStartLayerWidthRatio=.25;
    options.bodyControls={control};
    mesh::Mesh mesh;
    const auto result=algo::meshGeometryModel(*geometry.model,options,mesh);
    if(!result.published) {
        std::cerr<<algo::geometryReasonName(result.reason)<<'\n';return 2;
    }
    const auto &body=result.analysis.bodies.front();
    std::cout<<algo::geometryStrategyName(body.strategy)<<": "<<body.layers
             <<" layers; last/first width="<<body.endToStartLayerWidthRatio
             <<"; audit temporary bytes="<<result.audit.temporaryBytes<<'\n';
    return sdk::writeMeshArchive(mesh,"graded-example.paacsmesh").completed?0:3;
}
