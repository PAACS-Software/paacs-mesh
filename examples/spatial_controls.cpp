// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/AdaptiveBoxMesher.hpp"
#include "paacs/algo/MeshSizeAudit.hpp"
#include "paacs/sizing/ControlField.hpp"
#include <iostream>

int main() {
    using namespace paacs;
    sizing::ControlFieldSpecification specification;
    specification.minimumSize = .03;
    specification.maximumSize = .4;
    specification.sizeFactor = 1;
    // Dependencies refer to earlier node indices. Fine near a segment or ball;
    // Minimum means either influence can request refinement.
    specification.nodes = {sizing::field::distanceSegment({.15, .5, .1}, {.85, .5, .9}),
                           sizing::field::threshold(0, .015, .3, .07, .4),
                           sizing::field::distanceBall({.7, .5, .25}, .08),
                           sizing::field::threshold(2, 0, .2, .07, .4, true),
                           sizing::field::minimum(1, 3)};
    specification.output = 4;
    auto compiled = sizing::compileControlField(specification);
    if (!compiled.compiled) {
        std::cerr << sizing::controlFieldReasonName(compiled.reason) << '\n';
        return 1;
    }
    algo::AdaptiveBoxRequest request;
    request.field = &compiled.field; // Keep this immutable field alive for the call.
    request.threads = 1;
    request.maximumSizingWork = 1'000'000;
    mesh::Mesh mesh;
    auto result = algo::meshAdaptiveBox(request, mesh);
    if (!result.published) {
        std::cerr << algo::adaptiveBoxReasonName(result.reason) << '\n';
        return 2;
    }
    algo::MeshSizeAuditOptions policy;
    policy.maximumEdgeToTargetRatio = 1.05;
    auto size = algo::auditMeshSize(mesh, compiled.field, policy);
    if (!size.complete || !size.requirementMet)
        return 3;
    std::cout << result.tetrahedra << " Tet4 cells; conservative size ratio "
              << size.maximumConservativeRatio << "; field payload "
              << compiled.field.logicalBytes() << " bytes\n";
}
