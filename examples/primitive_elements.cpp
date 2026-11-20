// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <iostream>
#include <paacs/algo/PrimitiveMesher.hpp>
#include <paacs/sdk/MeshArchive.hpp>
int main() {
    using namespace paacs;
    for (auto family : {algo::PrimitiveElements::Tet4, algo::PrimitiveElements::Tet10}) {
        algo::PrimitiveShape shape; // Unit-radius cylinder, height two.
        algo::PrimitiveMeshingOptions options;
        options.elements = family;
        options.cylinderStrategy = algo::PrimitiveCylinderStrategy::Uniform;
        options.targetSize = .25;
        options.threads = 4;
        options.warningQuality = .1;
        options.qualityPolicy = algo::PrimitiveQualityPolicy::Warn;
        mesh::Mesh mesh;
        const auto result = algo::meshPrimitive(shape, options, mesh);
        if (!result.published) {
            std::cerr << algo::primitiveReasonName(result.reason) << '\n';
            return 1;
        }
        if (result.qualityWarning)
            std::cerr << result.quality.lowQualityCells << " poor cells; worst cell "
                      << result.quality.worstCell << '\n';
        const std::string path =
            std::string("cylinder-") + algo::primitiveElementsName(family) + ".paacsmesh";
        if (!sdk::writeMeshArchive(mesh, path).completed)
            return 1;
        std::cout << path << ": " << result.tet4s + result.tet10s << " cells, " << result.nodes
                  << " nodes\n";
    }
}
