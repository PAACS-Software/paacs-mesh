// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/StraightPrismLayerMesher.hpp"
#include "paacs/algo/StructuredGrowth.hpp"
#include "paacs/algo/StructuredBlockControlPlanner.hpp"
#include "paacs/algo/StructuredPolyhedralMesher.hpp"
#include "paacs/mesh/PolyhedralMesh.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sizing/MeshControls.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <vector>

extern int gFailures;
static const char* gControlsCurrent = "";
#define CONTROLS_TEST(name)                 \
    do {                                    \
        gControlsCurrent = name;            \
        std::printf("[test] %s\n", name); \
    } while (0)
#define CONTROLS_CHECK(cond)                                      \
    do {                                                          \
        if (!(cond)) {                                            \
            std::printf("  FAIL %s:%d  %s  (in %s)\n", __FILE__,\
                        __LINE__, #cond, gControlsCurrent);        \
            ++gFailures;                                          \
        }                                                         \
    } while (0)

namespace {
namespace algo = paacs::algo;
namespace mesh = paacs::mesh;
namespace sizing = paacs::sizing;
using paacs::Index;
using paacs::Vec3;

algo::StructuredBlockGeometry block(std::uint32_t u, std::uint32_t v,
                                    std::uint32_t w) {
    algo::StructuredBlockGeometry out;
    out.origin = {-1.0, 2.0, 3.0};
    out.axisU = {2.0, 0.0, 0.0};
    out.axisV = {0.2, 3.0, 0.0};
    out.axisW = {0.1, 0.15, 4.0};
    out.cellsU = u;
    out.cellsV = v;
    out.cellsW = w;
    out.volumeClassification = {3, 41};
    for (std::size_t i = 0; i < out.boundaryClassifications.size(); ++i)
        out.boundaryClassifications[i] =
            {2, static_cast<int>(101 + i)};
    return out;
}

mesh::Mesh squareSurface(std::uint32_t n) {
    mesh::Mesh out;
    const std::uint64_t side = static_cast<std::uint64_t>(n) + 1U;
    out.reserveAppend(static_cast<std::size_t>(side * side),
                      static_cast<std::size_t>(2U * n * n),
                      static_cast<std::size_t>(6U * n * n));
    for (std::uint32_t j = 0; j <= n; ++j)
        for (std::uint32_t i = 0; i <= n; ++i) {
            const bool boundary = i == 0 || j == 0 || i == n || j == n;
            (void)out.addNode(
                {static_cast<double>(i) / n,
                 static_cast<double>(j) / n, 0.0},
                boundary ? mesh::Classification{1, 10}
                         : mesh::Classification{2, 20});
        }
    const auto node = [n](std::uint32_t i, std::uint32_t j) {
        return static_cast<Index>(
            static_cast<std::uint64_t>(j) * (n + 1U) + i);
    };
    for (std::uint32_t j = 0; j < n; ++j)
        for (std::uint32_t i = 0; i < n; ++i) {
            const Index a = node(i, j);
            const Index b = node(i + 1U, j);
            const Index c = node(i + 1U, j + 1U);
            const Index d = node(i, j + 1U);
            (void)out.addCell(mesh::CellType::Tri3, {a, b, c}, {2, 20});
            (void)out.addCell(mesh::CellType::Tri3, {a, c, d}, {2, 20});
        }
    return out;
}

void compiledControls() {
    CONTROLS_TEST("controls.v6 compile once and apply deterministic precedence");
    sizing::MeshControlSpecification spec;
    spec.minimumSize = 0.05;
    spec.maximumSize = 2.0;
    spec.maximumGrowthRate = 1.4;
    spec.entitySizes = {{3, 7, 0.8}, {2, 12, 0.25}, {3, 9, 0.5},
                        {2, 12, 0.25}};
    spec.bodyMethods = {{9, sizing::VolumeEngineKind::MappedHex8},
                        {7, sizing::VolumeEngineKind::Tet4}};
    spec.prismLayers = {{12, 5, 0.01, 1.2, 0.2}};
    const auto out = sizing::compileMeshControls(spec);
    CONTROLS_CHECK(out.compiled);
    CONTROLS_CHECK(out.reason == sizing::MeshControlCompileReason::None);
    CONTROLS_CHECK(out.controls.entitySizeCount() == 3);
    CONTROLS_CHECK(out.controls.bodyMethodCount() == 2);
    CONTROLS_CHECK(out.controls.prismLayerCount() == 1);
    CONTROLS_CHECK(std::abs(out.controls.effectiveSize(2, 12, 1.0) -
                            0.25) < 1.0e-15);
    CONTROLS_CHECK(std::abs(out.controls.effectiveInterfaceSize(
                            12, 7, 9, 1.0) - 0.25) < 1.0e-15);
    CONTROLS_CHECK(std::abs(out.controls.effectiveSize(2, 99, 0.01) -
                            0.05) < 1.0e-15);
    CONTROLS_CHECK(out.controls.engineForBody(9) ==
                   sizing::VolumeEngineKind::MappedHex8);
    CONTROLS_CHECK(out.controls.engineForBody(100) ==
                   sizing::VolumeEngineKind::Auto);
    const auto* layer = out.controls.prismLayerForFace(12);
    CONTROLS_CHECK(layer != nullptr);
    if (layer) CONTROLS_CHECK(layer->layers == 5);

    spec.entitySizes.push_back({2, 12, 0.3});
    const auto conflict = sizing::compileMeshControls(spec);
    CONTROLS_CHECK(!conflict.compiled);
    CONTROLS_CHECK(conflict.reason ==
                   sizing::MeshControlCompileReason::ConflictingDuplicateSize);
}


void structuredBodyPlanner() {
    CONTROLS_TEST(
        "controls.v6 body/face sizing compiles direct structured counts");
    sizing::MeshControlSpecification spec;
    spec.minimumSize = 0.01;
    spec.maximumSize = 10.0;
    spec.maximumGrowthRate = 1.3;
    spec.entitySizes = {{3, 41, 0.5}, {2, 101, 0.25}};
    spec.bodyMethods = {{41, sizing::VolumeEngineKind::MappedHex8}};
    spec.prismLayers = {{101, 6, 0.002, 1.2, 0.03}};
    const auto compiled = sizing::compileMeshControls(spec);
    CONTROLS_CHECK(compiled.compiled);

    algo::StructuredBlockSizingPolicy policy;
    policy.baseSize = 1.0;
    policy.maximumCells = {128, 128, 128};
    policy.growth.u.mode = algo::StructuredAxisGrowthMode::
        EndToStartCellWidthRatio;
    policy.growth.u.endToStartCellWidthRatio = 4.0;
    const auto plan = algo::planStructuredBlockFromControls(
        block(1, 1, 1), 41, compiled.controls, policy);
    CONTROLS_CHECK(plan.planned);
    CONTROLS_CHECK(plan.reason == algo::StructuredBlockPlanReason::None);
    CONTROLS_CHECK(plan.requestedEngine ==
                   sizing::VolumeEngineKind::MappedHex8);
    CONTROLS_CHECK(plan.geometry.cellsU == 8);
    CONTROLS_CHECK(plan.geometry.cellsV == 13);
    CONTROLS_CHECK(plan.geometry.cellsW == 9);
    CONTROLS_CHECK(plan.controlLookups == 16);

    algo::MappedHexBlockRequest request;
    request.geometry = plan.geometry;
    request.growth = plan.growth;
    request.validation = algo::MappedHexValidationMode::FullTopologyAudit;
    mesh::Mesh generated;
    const auto outcome = algo::generateMappedHexBlock(request, generated);
    CONTROLS_CHECK(outcome.published());
    CONTROLS_CHECK(outcome.metrics.hexes == 936);
    CONTROLS_CHECK(outcome.metrics.audit.valid);

    algo::StructuredPolyBlockRequest polyRequest;
    polyRequest.geometry = plan.geometry;
    polyRequest.growth = plan.growth;
    polyRequest.validation = algo::StructuredPolyValidationMode::FullAudit;
    mesh::PolyhedralMesh poly;
    const auto polyOutcome = algo::generateStructuredPolyhedralBlock(
        polyRequest, poly);
    CONTROLS_CHECK(polyOutcome.published());
    CONTROLS_CHECK(polyOutcome.metrics.cells == 936);
    CONTROLS_CHECK(polyOutcome.metrics.audit.valid);

    algo::StraightPrismLayerRequest layerRequest;
    CONTROLS_CHECK(algo::applyCompiledPrismLayerControl(
        101, compiled.controls, layerRequest));
    CONTROLS_CHECK(layerRequest.layers == 6);
    CONTROLS_CHECK(std::abs(layerRequest.firstLayerThickness - 0.002) <
                   1.0e-15);
    CONTROLS_CHECK(std::abs(layerRequest.cellToCellGrowth - 1.2) <
                   1.0e-15);
    CONTROLS_CHECK(std::abs(layerRequest.maximumAllowedGrowth - 1.3) <
                   1.0e-15);
    CONTROLS_CHECK(!algo::applyCompiledPrismLayerControl(
        999, compiled.controls, layerRequest));
}

void geometricAxisGrowth() {
    CONTROLS_TEST("growth.v6 geometric schedule is bounded and normalized");
    algo::StructuredAxisGrowth growth;
    growth.mode = algo::StructuredAxisGrowthMode::
        EndToStartCellWidthRatio;
    growth.endToStartCellWidthRatio = 4.0;
    growth.maximumAdjacentCellWidthRatio = 1.3;
    const auto out = algo::buildStructuredAxisDistribution(12, growth);
    CONTROLS_CHECK(out.valid);
    CONTROLS_CHECK(out.boundaries.size() == 13);
    CONTROLS_CHECK(out.boundaries.front() == 0.0);
    CONTROLS_CHECK(out.boundaries.back() == 1.0);
    CONTROLS_CHECK(out.maximumAdjacentRatio <= 1.3);
    CONTROLS_CHECK(std::abs(out.lastCellFraction /
                            out.firstCellFraction - 4.0) < 1.0e-12);
    for (std::size_t i = 1; i < out.boundaries.size(); ++i)
        CONTROLS_CHECK(out.boundaries[i] > out.boundaries[i - 1]);

    growth.endToStartCellWidthRatio = 1000.0;
    growth.maximumAdjacentCellWidthRatio = 1.05;
    const auto rejected = algo::buildStructuredAxisDistribution(4, growth);
    CONTROLS_CHECK(!rejected.valid);
    CONTROLS_CHECK(rejected.reason ==
                   algo::StructuredGrowthReason::AdjacentLimitExceeded);
}

void gradedHexAndPoly() {
    CONTROLS_TEST("hex/poly.v6 graded direct generation preserves volume");
    algo::MappedHexBlockRequest hexRequest;
    hexRequest.geometry = block(12, 10, 8);
    hexRequest.growth.u.mode = algo::StructuredAxisGrowthMode::
        EndToStartCellWidthRatio;
    hexRequest.growth.u.endToStartCellWidthRatio = 3.0;
    hexRequest.growth.v.mode = algo::StructuredAxisGrowthMode::
        EndToStartCellWidthRatio;
    hexRequest.growth.v.endToStartCellWidthRatio = 0.5;
    hexRequest.validation = algo::MappedHexValidationMode::FullTopologyAudit;
    mesh::Mesh hex;
    const auto ho = algo::generateMappedHexBlock(hexRequest, hex);
    CONTROLS_CHECK(ho.published());
    CONTROLS_CHECK(ho.metrics.audit.valid);
    CONTROLS_CHECK(ho.metrics.maximumAdjacentGrowth <= 1.3);
    CONTROLS_CHECK(ho.metrics.minimumCellFraction > 0.0);
    CONTROLS_CHECK(std::abs(ho.metrics.audit.totalVolume - 24.0) < 1.0e-10);
    const auto firstU = (hex.node(1) - hex.node(0)).norm();
    const auto lastBase = static_cast<Index>(hexRequest.geometry.cellsU - 1U);
    const auto lastU = (hex.node(lastBase + 1) - hex.node(lastBase)).norm();
    CONTROLS_CHECK(std::abs(lastU / firstU - 3.0) < 1.0e-11);

    algo::StructuredPolyBlockRequest polyRequest;
    polyRequest.geometry = hexRequest.geometry;
    polyRequest.growth = hexRequest.growth;
    polyRequest.validation = algo::StructuredPolyValidationMode::FullAudit;
    mesh::PolyhedralMesh poly;
    const auto po = algo::generateStructuredPolyhedralBlock(polyRequest, poly);
    CONTROLS_CHECK(po.published());
    CONTROLS_CHECK(po.metrics.audit.valid);
    CONTROLS_CHECK(poly.cellCount() == 960);
    CONTROLS_CHECK(std::abs(po.metrics.audit.totalVolume - 24.0) < 1.0e-10);
    CONTROLS_CHECK(std::abs(poly.node(1).x - hex.node(1).x) < 1.0e-15);
}

void straightPrismLayers() {
    CONTROLS_TEST("prism.v6 straight graded layers publish audited Prism6");
    const mesh::Mesh source = squareSurface(4);
    algo::StraightPrismLayerRequest request;
    request.layers = 5;
    request.firstLayerThickness = 0.01;
    request.cellToCellGrowth = 1.2;
    request.maximumAllowedGrowth = 1.3;
    request.maximumTotalThickness = 0.1;
    request.volumeClassification = {3, 77};
    request.targetClassification = {2, 21};
    request.sideClassification = {2, 22};
    request.validation = algo::PrismLayerValidationMode::FullTopologyAudit;
    mesh::Mesh layers;
    const auto out = algo::generateStraightPrismLayers(
        source, request, layers);
    CONTROLS_CHECK(out.published());
    CONTROLS_CHECK(out.metrics.audit.valid);
    CONTROLS_CHECK(out.metrics.sourceTriangles == 32);
    CONTROLS_CHECK(out.metrics.sourceBoundaryEdges == 16);
    CONTROLS_CHECK(out.metrics.prisms == 160);
    CONTROLS_CHECK(out.metrics.boundaryTriangles == 64);
    CONTROLS_CHECK(out.metrics.boundaryQuads == 80);
    const double expectedThickness = 0.01 *
        (1.0 - std::pow(1.2, 5.0)) / (1.0 - 1.2);
    CONTROLS_CHECK(std::abs(out.metrics.totalThickness -
                            expectedThickness) < 1.0e-14);
    CONTROLS_CHECK(std::abs(out.metrics.audit.totalVolume -
                            expectedThickness) < 1.0e-12);
    CONTROLS_CHECK(out.metrics.audit.minimumPrismVolume > 0.0);
    CONTROLS_CHECK(layers.countCells(mesh::CellType::Prism6) == 160);

    mesh::PolyhedralMesh poly;
    const auto converted = algo::convertSupportedVolumeMeshToPolyhedral(
        layers, poly);
    CONTROLS_CHECK(converted.published());
    CONTROLS_CHECK(converted.metrics.audit.valid);
    CONTROLS_CHECK(poly.cellCount() == 160);
    CONTROLS_CHECK(std::abs(converted.metrics.audit.totalVolume -
                            expectedThickness) < 1.0e-12);
}

void prismContracts() {
    CONTROLS_TEST("prism.v6 invalid growth and orientation are atomic");
    const mesh::Mesh source = squareSurface(1);
    mesh::Mesh destination;
    (void)destination.addNode({9, 8, 7}, {0, 99});
    const auto before = mesh::meshHash(destination);
    algo::StraightPrismLayerRequest request;
    request.layers = 4;
    request.firstLayerThickness = 0.1;
    request.cellToCellGrowth = 1.5;
    request.maximumAllowedGrowth = 1.3;
    const auto growth = algo::generateStraightPrismLayers(
        source, request, destination);
    CONTROLS_CHECK(!growth.published());
    CONTROLS_CHECK(growth.reason == algo::PrismLayerReason::
                   GrowthLimitExceeded);
    CONTROLS_CHECK(growth.destinationUnchanged);
    CONTROLS_CHECK(mesh::meshHash(destination) == before);

    request.cellToCellGrowth = 1.1;
    request.direction = {0, 0, -1};
    const auto orientation = algo::generateStraightPrismLayers(
        source, request, destination);
    CONTROLS_CHECK(!orientation.published());
    CONTROLS_CHECK(orientation.reason == algo::PrismLayerReason::
                   InconsistentSourceOrientation);
    CONTROLS_CHECK(mesh::meshHash(destination) == before);
}

}  // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* meshControlLayerTests(std::size_t& count) {
    static const Entry tests[]{
        {"controls.v6.compile", compiledControls},
        {"controls.v6.structured_plan", structuredBodyPlanner},
        {"growth.v6.axis", geometricAxisGrowth},
        {"hex.v6.graded", gradedHexAndPoly},
        {"prism.v6.layers", straightPrismLayers},
        {"prism.v6.contracts", prismContracts}};
    count = std::size(tests);
    return tests;
}
}  // namespace paacs_tests
