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
    // Rectangular through-hole: the two loops have opposite winding.
    const std::vector<std::vector<Vec2>> loops{
        {{0,0},{3,0},{3,3},{0,3}}, {{1,1},{1,2},{2,2},{2,1}}};
    auto geometry=geo::makePolyhedralModel(geo::extrudedPolygon(loops,{.2,0,1}));
    if(!geometry.built())return 1;
    algo::GeometryMeshingOptions options;
    options.targetSize=.25;
    options.maximumNodes=100000;
    options.maximumOutputBytes=64ULL*1024*1024;
    options.minimumElementQuality=.05;
    mesh::Mesh result;
    const auto outcome=algo::meshGeometryModel(*geometry.model,options,result);
    if(!outcome.published) {
        std::cerr<<algo::geometryReasonName(outcome.reason)<<'\n';return 2;
    }
    std::cout<<algo::geometryStrategyName(outcome.analysis.bodies[0].strategy)<<" "
             <<result.nodeCount()<<" nodes, "<<result.cellCount()<<" cells\n";
    return sdk::writeMeshArchive(result,"example.paacsmesh").completed?0:3;
}
