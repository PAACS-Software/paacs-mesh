// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <iostream>
#include <paacs/algo/LayeredPatchMesher.hpp>
#include <paacs/algo/PrimitiveSizeControls.hpp>
int main() {
    using namespace paacs;
    algo::PrimitiveShape cylinder;
    cylinder.radius = .6;
    cylinder.height = 4.8;
    sizing::CappedSizeSpecification sizes;
    sizes.length = cylinder.height;
    sizes.globalSize = .36;
    sizes.controls = {algo::primitiveFaceSizeControl(cylinder, 2, .12, 1.25)};
    const auto compiled = sizing::compileCappedSizeProfile(sizes);
    if (!compiled.compiled)
        return 1;
    for (auto family : {algo::PrimitiveElements::Tet4, algo::PrimitiveElements::Tet10,
                        algo::PrimitiveElements::HexDominant}) {
        algo::PrimitiveMeshingOptions request;
        request.elements = family;
        request.sizeProfile = &compiled.profile;
        request.profileStrategy = algo::PrimitiveProfileStrategy::Volume;
        request.maximumChordError = .015;
        request.maximumRelativeVolumeError = .03;
        mesh::Mesh output;
        const auto result = algo::meshPrimitive(cylinder, request, output);
        if (!result.published) {
            std::cerr << algo::primitiveReasonName(result.reason) << '\n';
            return 2;
        }
        std::cout << algo::primitiveElementsName(family) << ": " << output.nodeCount() << " nodes; "
                  << result.quality.lowQualityCells << " low-quality warnings\n";
    }
    // The same physical profile drives prism thickness on any validated planar patch.
    mesh::Mesh patch;
    patch.addNode({0, 0, 0}, {2, 1});
    patch.addNode({1, 0, 0}, {2, 1});
    patch.addNode({0, 1, 0}, {2, 1});
    patch.addCell(mesh::CellType::Tri3, {0, 1, 2}, {2, 1});
    algo::LayeredPatchRequest layers;
    layers.sizeProfile = &compiled.profile;
    mesh::Mesh prisms;
    return algo::meshLayeredPatch(patch, layers, prisms).published ? 0 : 3;
}
