// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <iostream>
#include <paacs/algo/PrimitiveMesher.hpp>
#include <paacs/algo/SurfaceFeatures.hpp>

int main() {
    using namespace paacs;
    algo::PrimitiveShape shape;
    shape.kind = algo::PrimitiveKind::Elbow;
    algo::PrimitiveMeshingOptions options;
    options.elements = algo::PrimitiveElements::Tet10;
    options.circularStrategy = algo::PrimitiveCircularStrategy::Uniform;
    options.targetSize = .2;
    options.threads = 4;
    mesh::Mesh mesh;
    const auto result = algo::meshPrimitive(shape, options, mesh);
    if (!result.published)
        return 1;

    // A separate read-only request: ordinary meshing has no feature-analysis cost.
    algo::SurfaceFeatureOptions features;
    features.enabled = true;
    features.angleDegrees = 30;
    features.angularToleranceDegrees = 1e-8;
    features.minimumLength = .01; // Output filter, not a mesh-sizing constraint.
    const auto detected = algo::detectSurfaceFeatures(mesh, features);
    if (!detected.completed || detected.edges.size() != 8 * result.segmentsPerQuarter)
        return 2;
    std::cout << result.tet10s << " certified Tet10 cells; " << detected.edges.size()
              << " rim edges; minimum mean ratio " << result.minimumQuality << '\n';
    return 0;
}
