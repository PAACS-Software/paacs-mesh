// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/geometry/PolyhedralModel.hpp"
#include <iostream>
int main() {
    using namespace paacs;
    // Rectangular duct with two holes; each loop is a polygonal CAD boundary.
    const std::vector<std::vector<Vec2>> rings{{{0, 0}, {4, 0}, {4, 3}, {0, 3}},
                                               {{.4, .4}, {1.7, .4}, {1.7, 2.6}, {.4, 2.6}},
                                               {{2.3, .4}, {3.6, .4}, {3.6, 2.6}, {2.3, 2.6}}};
    auto model = geo::makePolyhedralModel(geo::extrudedPolygon(rings, {0, 0, 4}));
    if (!model.built())
        return 1;
    sizing::CappedSizeSpecification specification;
    specification.length = 4;
    specification.globalSize = .4;
    specification.controls = {{0, 0, .1, 1.25}, {2, 2.4, .16, 1.2}};
    auto profile = sizing::compileCappedSizeProfile(specification);
    if (!profile.compiled)
        return 2;
    algo::GeometryBodyControl body;
    body.solidTag = 1;
    body.sourceFaceTag = 1; // extrudedPolygon assigns bottom=1, top=2.
    body.sizeProfile = &profile.profile;
    body.elements = algo::GeometryElements::Tet4; // Native retains Prism6 here.
    algo::GeometryMeshingOptions options;
    options.targetSize = .3; // Independent tangential source-face resolution.
    options.warningElementQuality = .15;
    options.bodyControls = {body};
    mesh::Mesh mesh;
    const auto outcome = algo::meshGeometryModel(*model.model, options, mesh);
    if (!outcome.published) {
        std::cerr << algo::geometryReasonName(outcome.reason) << '\n';
        return 3;
    }
    const auto &report = outcome.analysis.bodies.front();
    std::cout << mesh.nodeCount() << " nodes, " << report.layers << " layers, quality "
              << outcome.quality.minimumQuality << ", warning=" << outcome.qualityWarning << '\n';
}
