// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/AdaptiveBoxMesher.hpp"
#include "paacs/algo/LayeredPatchMesher.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <iostream>

int main() {
    using namespace paacs;
    sizing::RegionSizeSpecification sizes;
    sizes.baseSize = .25;
    // A finite rectangle on the selected cube face; its influence decays
    // with physical distance through the volume.
    sizes.controls.push_back(
        sizing::boxFaceSizeControl({0, 0, 0}, {1, 1, 1}, sizing::BoxFace::XMin, .0625, 1.3));
    auto compiled = sizing::compileRegionSizeField(sizes);
    if (!compiled.compiled)
        return 1;
    algo::AdaptiveBoxRequest request;
    request.field = &compiled.field;
    request.maximumAdjacentSizeRatio = 1.3; // Optional achieved-size acceptance gate.
    mesh::Mesh volume;
    const auto result = algo::meshAdaptiveBox(request, volume);
    if (!result.published)
        return 2;
    std::cout << result.tetrahedra << " tets; measured adjacent size ratio "
              << result.audit.maximumAdjacentSizeRatio << '\n';
    if (!sdk::writeMeshArchive(volume, "local-controls-example.paacsmesh").completed)
        return 3;

    // A prepared Tri3 patch can be structured or unstructured. This minimal
    // example uses two triangles; all shared nodes remain shared through the stack.
    mesh::Mesh patch;
    for (Vec3 p : std::array<Vec3, 4>{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}})
        patch.addNode(p, {2, 1});
    patch.addCell(mesh::CellType::Tri3, {0, 1, 2}, {2, 1});
    patch.addCell(mesh::CellType::Tri3, {0, 2, 3}, {2, 1});
    algo::LayeredPatchRequest layers;
    layers.schedule.firstLayerThickness = .01;
    layers.schedule.wallLayers = 5;
    layers.schedule.growthRatio = 1.2;
    layers.schedule.totalDepth = 1;
    layers.schedule.maximumCoreThickness = .15;
    mesh::Mesh wedges;
    const auto stack = algo::meshLayeredPatch(patch, layers, wedges);
    if (!stack.published)
        return 4;
    std::cout << stack.sweep.metrics.wedges << " wedges; " << stack.schedule.wallLayers
              << " wall + " << stack.schedule.coreLayers << " core layers\n";
    return sdk::writeMeshArchive(wedges, "layer-controls-example.paacsmesh").completed ? 0 : 5;
}
