// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <iostream>
#include <paacs/algo/PrimitiveMesher.hpp>
#include <paacs/sdk/MeshArchive.hpp>
int main() {
    for (auto kind : {paacs::algo::PrimitiveKind::Cylinder, paacs::algo::PrimitiveKind::Cone,
                      paacs::algo::PrimitiveKind::Torus}) {
        paacs::algo::PrimitiveShape shape;
        shape.kind = kind;
        paacs::algo::PrimitiveMeshingOptions options;
        options.targetSize = .15;
        options.maximumChordError = .001;
        options.threads = 4;
        paacs::mesh::Mesh mesh;
        const auto result = paacs::algo::meshPrimitive(shape, options, mesh);
        if (!result.published) {
            std::cerr << paacs::algo::primitiveReasonName(result.reason) << '\n';
            return 1;
        }
        if (!paacs::sdk::writeMeshArchive(mesh, std::string(paacs::algo::primitiveKindName(kind)) +
                                                    ".paacsmesh")
                 .completed)
            return 1;
        std::cout << paacs::algo::primitiveKindName(kind) << ": " << result.hexes + result.pyramids
                  << " volume cells, " << 100 * result.relativeVolumeError << "% volume error\n";
    }
}
