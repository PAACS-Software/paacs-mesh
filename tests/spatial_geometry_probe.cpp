// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/geometry/PolyhedralModel.hpp"
#include "paacs/sizing/ControlField.hpp"
#include <iomanip>
#include <iostream>
using namespace paacs;
int main() {
    const std::vector<std::vector<Vec2>> loops{{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
    auto model = geo::makePolyhedralModel(geo::extrudedPolygon(loops, {0, 0, 1}));
    if (!model.built())
        return 2;
    for (double fine : {.5, .4, .3, .25})
        for (double ratio : {0., 1.5, 2., 3.}) {
            sizing::ControlFieldSpecification s;
            s.maximumSize = .6;
            s.nodes = {sizing::field::distancePoint({.5, .5, .5}),
                       sizing::field::threshold(0, 0, .7, fine, .6)};
            s.output = 1;
            auto f = sizing::compileControlField(s);
            if (!f.compiled)
                return 3;
            algo::GeometryMeshingOptions o;
            o.targetSize = .6;
            o.spatialField = &f.field;
            o.threads = 1;
            o.warningElementSizeRatio = 1.35;
            if (ratio > 0) {
                algo::GeometryBodyControl c;
                c.solidTag = 1;
                c.tetRefinement.enabled = true;
                c.tetRefinement.maximumEdgeToTarget = ratio;
                o.bodyControls = {c};
            }
            mesh::Mesh m;
            const auto r = algo::meshGeometryModel(*model.model, o, m);
            std::cout << std::setprecision(17) << "{\"fine\":" << fine
                      << ",\"refinement_ratio\":" << ratio
                      << ",\"published\":" << (r.published ? "true" : "false") << ",\"reason\":\""
                      << algo::geometryReasonName(r.reason) << "\",\"nodes\":" << m.nodeCount()
                      << ",\"cells\":" << m.cellCount();
            if (!r.analysis.bodies.empty()) {
                const auto &b = r.analysis.bodies[0];
                std::cout << ",\"boundary_reason\":" << b.boundaryReason
                          << ",\"refinement_reason\":" << b.refinementReason
                          << ",\"insertions\":" << b.refinementInsertions
                          << ",\"ratio\":" << b.sizeAudit.maximumConservativeRatio
                          << ",\"warn_cells\":" << b.sizeAudit.cellsAboveWarning;
            }
            std::cout << "}\n" << std::flush;
        }
}
