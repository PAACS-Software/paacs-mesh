// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "geometry_checks.hpp"
#include "geometry_fixtures.hpp"
#include "geometry_model_views.hpp"
#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/PreparedSweepMesher.hpp"
#include <iostream>

using namespace geometry_test;
namespace {
int checks = 0;
void verify(const Case &c, const geo::Model &model, algo::GeometryMeshingOptions o,
            unsigned expectedLayers, double expectedRatio, int expectedSource = -1) {
    mesh::Mesh m;
    const auto r = algo::meshGeometryModel(model, o, m);
    require(r.published && r.audit.valid,
            std::string("controlled mesh: ") + algo::geometryReasonName(r.reason));
    const auto &body = r.analysis.bodies.at(0);
    require(body.layers == expectedLayers, "reported layer count");
    require(expectedSource < 0 || body.sourceFaceTag == expectedSource, "selected source face");
    const auto *source = model.face(body.sourceFaceTag);
    Vec2 lo, hi;
    source->paramBounds(lo, hi);
    const auto origin = source->pointAt(lo), normal = source->normalAt(lo);
    const auto denom = body.sweep.dot(normal);
    std::vector<double> levels;
    for (Index i = 0; i < m.nodeCount(); ++i)
        levels.push_back((m.node(i) - origin).dot(normal) / denom);
    std::sort(levels.begin(), levels.end());
    levels.erase(std::unique(levels.begin(), levels.end(),
                             [](double a, double b) { return std::abs(a - b) < 1e-10; }),
                 levels.end());
    require(levels.size() == expectedLayers + 1, "actual node layer count");
    require(std::abs(levels.front()) < 1e-9 && std::abs(levels.back() - 1) < 1e-9,
            "layer endpoints");
    const double ratio = (levels.back() - levels[levels.size() - 2]) / (levels[1] - levels[0]);
    require(std::abs(ratio - expectedRatio) < 1e-7 * expectedRatio, "measured layer width ratio");
    require(std::abs(positiveVolume(m) - c.volume) < 1e-8 * c.volume, "controlled analytic volume");
    Vec3 a, b;
    model.solids()[0]->bbox(a, b);
    checkClassification(model, m, (b - a).norm() * 1e-7);
    const auto quality = algo::analyzeMixedMeshQuality(m);
    require(quality.valid && quality.invertedCells == 0, "controlled cell quality");
    o.threads = 4;
    mesh::Mesh again;
    require(algo::meshGeometryModel(model, o, again).published && digest(m) == digest(again),
            "controlled thread determinism");
    ++checks;
}
void controlledShapes() {
    const auto all = cases();
    for (const auto name : {"block_thin", "block_rotated", "block_skew", "concave_l_oblique",
                            "hole_oblique", "holes_two"}) {
        const auto &c =
            *std::find_if(all.begin(), all.end(), [&](const auto &v) { return v.name == name; });
        auto built = geo::makePolyhedralModel(c.input);
        require(built.built(), "controlled geometry");
        for (double ratio : {1., 4., .25}) {
            algo::GeometryMeshingOptions o;
            o.targetSize = c.targetSize;
            algo::GeometryBodyControl control;
            control.solidTag = 1;
            control.layers = 12;
            control.endToStartLayerWidthRatio = ratio;
            o.bodyControls = {control};
            verify(c, *built.model, o, 12, ratio);
            const auto analysis = algo::analyzeGeometryStrategies(*built.model, o);
            require(analysis.complete, "source selection analysis");
            const auto other = analysis.bodies[0].targetFaceTag;
            o.bodyControls[0].sourceFaceTag = other;
            verify(c, *built.model, o, 12, ratio, other);
        }
    }
}
void bodySpecific() {
    auto all = cases();
    auto a = geo::makePolyhedralModel(all[0].input);
    auto input = all[9].input;
    for (auto &p : input.vertices)
        p.x += 4;
    auto b = geo::makePolyhedralModel(input);
    CombinedModel combined;
    combined.add(*a.model, 0);
    combined.add(*b.model, 1000);
    algo::GeometryMeshingOptions o;
    o.targetSize = .5;
    algo::GeometryBodyControl control;
    control.solidTag = 1001;
    control.targetSize = .15;
    control.layers = 7;
    o.bodyControls = {control};
    mesh::Mesh m;
    const auto r = algo::meshGeometryModel(combined, o, m);
    require(r.published && r.audit.valid && r.analysis.bodies.size() == 2, "per-body mesh");
    require(r.analysis.bodies[0].targetSize == .5 && r.analysis.bodies[1].targetSize == .15 &&
                r.analysis.bodies[1].layers == 7,
            "per-body settings");
    require(std::abs(positiveVolume(m) - 2) < 1e-8, "per-body total volume");
    algo::GeometryBodyControl box;
    box.solidTag = 1;
    box.strategy = algo::GeometryStrategy::LinearPrismSweep;
    o.bodyControls.push_back(box);
    mesh::Mesh swept;
    const auto s = algo::meshGeometryModel(combined, o, swept);
    require(s.published && s.audit.valid &&
                s.analysis.bodies[0].strategy == algo::GeometryStrategy::LinearPrismSweep,
            "per-body strategy override");
    checks += 2;
}
void refusals() {
    auto model = geo::makePolyhedralModel(cases()[0].input);
    const auto reject = [&](algo::GeometryMeshingOptions o, algo::GeometryReason expected) {
        mesh::Mesh m;
        (void)m.addNode({9, 0, 0}, {0, 9});
        const auto before = digest(m);
        const auto r = algo::meshGeometryModel(*model.model, o, m);
        require(!r.published && r.destinationUnchanged && r.reason == expected &&
                    digest(m) == before,
                std::string("control refusal: ") + algo::geometryReasonName(r.reason));
        ++checks;
    };
    using R = algo::GeometryReason;
    for (int test = 0; test < 13; ++test) {
        algo::GeometryMeshingOptions o;
        algo::GeometryBodyControl c;
        c.solidTag = 1;
        auto reason = R::InvalidRequest;
        switch (test) {
        case 0:
            c.solidTag = 999;
            break;
        case 1:
            c.targetSize = -1;
            break;
        case 2:
            c.targetSize = NAN;
            break;
        case 3:
            c.endToStartLayerWidthRatio = 0;
            break;
        case 4:
            c.endToStartLayerWidthRatio = INFINITY;
            break;
        case 5:
            c.layers = 2;
            c.endToStartLayerWidthRatio = 100;
            break;
        case 6:
            c.layers = 1;
            c.endToStartLayerWidthRatio = 2;
            break;
        case 7:
            c.sourceFaceTag = 999;
            reason = R::NoMatchingStrategy;
            break;
        case 8:
            c.layers = 3;
            c.strategy = algo::GeometryStrategy::BoundaryConformingTet4;
            reason = R::NoMatchingStrategy;
            break;
        case 9:
            c.layers = UINT32_MAX;
            o.maximumNodes = 100;
            reason = R::OutputBudget;
            break;
        case 10:
            c.maximumAdjacentLayerRatio = .5;
            break;
        case 11:
            c.sourceFaceTag = -2;
            break;
        case 12:
            c.strategy = static_cast<algo::GeometryStrategy>(255);
            break;
        }
        o.bodyControls = {c};
        reject(o, reason);
    }
    algo::GeometryMeshingOptions o;
    algo::GeometryBodyControl c;
    c.solidTag = 1;
    o.bodyControls = {c, c};
    reject(o, R::InvalidRequest);
    o = {};
    o.maximumAuditTemporaryBytes = 1;
    reject(o, R::OutputBudget);
    // Direct sweep resource checks must precede allocation of the layer table.
    mesh::Mesh source;
    for (Vec3 p : {Vec3{0, 0, 0}, {1, 0, 0}, {0, 1, 0}})
        (void)source.addNode(p, {2, 1});
    (void)source.addCell(mesh::CellType::Tri3, {0, 1, 2}, {2, 1});
    algo::PreparedMixedSweepRequest request;
    request.layers = UINT32_MAX;
    request.maximumNodes = 100;
    mesh::Mesh held;
    (void)held.addNode({9, 0, 0}, {0, 9});
    const auto before = digest(held);
    const auto result = algo::generatePreparedMixedSweep(source, request, held);
    require(!result.published() && result.reason == algo::PreparedSweepReason::BudgetExceeded &&
                digest(held) == before,
            "sweep pre-allocation budget");
    ++checks;
}
void coordinateResolution() {
    // Each mathematical cell is positive, but double coordinates collapse
    // adjacent layers at this translation. An analytic frame is insufficient.
    algo::MappedHexBlockRequest block;
    block.geometry.origin = {0, 0, 1e15};
    block.geometry.cellsW = 12;
    block.executionThreads = 1;
    mesh::Mesh held;
    (void)held.addNode({9, 0, 0}, {0, 9});
    const auto before = digest(held);
    const auto hex = algo::generateMappedHexBlock(block, held);
    require(!hex.published() && digest(held) == before, "rounded mapped cells accepted");
    ++checks;
    mesh::Mesh cap;
    for (Vec3 p : {Vec3{0, 0, 1e15}, {1, 0, 1e15}, {0, 1, 1e15}})
        (void)cap.addNode(p, {2, 1});
    (void)cap.addCell(mesh::CellType::Tri3, {0, 1, 2}, {2, 1});
    algo::PreparedMixedSweepRequest sweep;
    sweep.layers = 12;
    sweep.fullAudit = false;
    sweep.executionThreads = 1;
    const auto prism = algo::generatePreparedMixedSweep(cap, sweep, held);
    require(!prism.published() && digest(held) == before, "rounded swept cells accepted");
    ++checks;
    require(hex.metrics.usedFullGeometricAudit && prism.metrics.usedFullGeometricAudit,
            "coordinate guard did not inspect stored cells");
    block = {};
    block.geometry.cellsW = UINT32_MAX;
    block.geometry.maxNodes = 100;
    const auto bounded = algo::generateMappedHexBlock(block, held);
    require(!bounded.published() && bounded.reason == algo::MappedHexReason::BudgetExceeded &&
                digest(held) == before,
            "mapped pre-allocation budget");
    ++checks;
    mesh::Mesh folded;
    for (Vec3 p : {Vec3{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, -.1}})
        (void)folded.addNode(p, {3, 1});
    (void)folded.addCell(mesh::CellType::Prism6, {0, 1, 2, 3, 4, 5}, {3, 1});
    const auto bad = algo::auditPreparedMixedSweep(folded);
    require(!bad.valid && bad.nonPositiveCells == 1, "folded positive-total-volume prism accepted");
    ++checks;
}

} // namespace
int main() {
    try {
        controlledShapes();
        bodySpecific();
        refusals();
        coordinateResolution();
        std::cout << checks << " graded-layer, body-control and rollback cases passed\n";
    } catch (const std::exception &e) {
        std::cerr << "case " << checks << ": " << e.what() << '\n';
        return 1;
    }
}
