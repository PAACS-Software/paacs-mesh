// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "geometry_checks.hpp"
#include "geometry_fixtures.hpp"
#include "geometry_model_views.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/algo/PreparedSweepMesher.hpp"
#include <iostream>
#include <limits>
#include <random>

using namespace geometry_test;
namespace {
int checks = 0;
void check(bool value, const std::string &name) {
    require(value, name);
    ++checks;
}
mesh::Mesh held() {
    mesh::Mesh m;
    (void)m.addNode({17, 19, 23}, {0, 42});
    return m;
}
void polyhedral() {
    const auto box = cases().front().input;
    const auto rejects = [&](const geo::PolyhedralInput &p, const std::string &name) {
        check(!geo::makePolyhedralModel(p).built(), name);
    };
    auto p = box;
    p.vertices[0].x = std::numeric_limits<double>::quiet_NaN();
    rejects(p, "NaN geometry");
    p = box;
    p.vertices[0].z = std::numeric_limits<double>::infinity();
    rejects(p, "infinite geometry");
    p = box;
    p.faces.pop_back();
    rejects(p, "open shell");
    p = box;
    p.faces.push_back(p.faces[0]);
    rejects(p, "nonmanifold shell");
    p = box;
    std::reverse(p.faces[0].loops[0].begin(), p.faces[0].loops[0].end());
    rejects(p, "one reversed face");
    p = box;
    p.vertices[6].z += .1;
    rejects(p, "nonplanar face");
    p = box;
    p.faces[0].loops[0][0] = 999;
    rejects(p, "invalid vertex index");
    p = box;
    p.faces[0].loops[0][1] = p.faces[0].loops[0][0];
    rejects(p, "repeated loop vertex");
    p = box;
    p.vertices.push_back({5, 5, 5});
    rejects(p, "unused vertex");
    p = box;
    for (auto &f : p.faces)
        for (auto &l : f.loops)
            std::reverse(l.begin(), l.end());
    check(geo::makePolyhedralModel(p).built(), "global inward shell normalization");
    p = box;
    for (auto &v : p.vertices)
        v.z = 0;
    rejects(p, "zero volume");
    std::vector<std::vector<Vec2>> loops{{{0, 0}, {2, 2}, {2, 0}, {0, 2}}};
    rejects(geo::extrudedPolygon(loops, {0, 0, 1}), "crossed outer loop");
    loops = {{{0, 0}, {2, 0}, {2, 2}, {0, 2}}, {{3, 3}, {3, 4}, {4, 4}, {4, 3}}};
    rejects(geo::extrudedPolygon(loops, {0, 0, 1}), "outside hole");
    loops = {{{0, 0}, {2, 0}, {2, 2}, {0, 2}},
             {{.5, .5}, {.5, 1.5}, {1.5, 1.5}, {1.5, .5}},
             {{1, 1}, {1, 1.7}, {1.7, 1.7}, {1.7, 1}}};
    rejects(geo::extrudedPolygon(loops, {0, 0, 1}), "crossing holes");
    loops.back() = {{.7, .7}, {.7, 1}, {1, 1}, {1, .7}};
    rejects(geo::extrudedPolygon(loops, {0, 0, 1}), "nested holes");
    p = box;
    const auto n = p.vertices.size();
    for (auto v : box.vertices)
        p.vertices.push_back(v + Vec3{3, 0, 0});
    for (auto f : box.faces) {
        for (auto &l : f.loops)
            for (auto &v : l)
                v += static_cast<std::uint32_t>(n);
        p.faces.push_back(std::move(f));
    }
    check(geo::makePolyhedralModel(p).reason == geo::PolyhedralReason::DisconnectedShell,
          "disconnected shell");
    geo::PolyhedralLimits limit;
    limit.maximumVertices = 7;
    check(geo::makePolyhedralModel(box, limit).reason == geo::PolyhedralReason::BudgetExceeded,
          "vertex budget");
    limit = {};
    limit.maximumFaces = 5;
    check(!geo::makePolyhedralModel(box, limit).built(), "face budget");
    limit = {};
    limit.maximumWork = 2;
    check(!geo::makePolyhedralModel(box, limit).built(), "geometry work budget");
    limit = {};
    limit.relativeTolerance = 0;
    check(!geo::makePolyhedralModel(box, limit).built(), "bad geometry tolerance");
}
void options() {
    auto model = geo::makePolyhedralModel(cases().front().input);
    const auto rejects = [&](algo::GeometryMeshingOptions o, algo::GeometryReason expected,
                             const std::string &name) {
        auto m = held();
        const auto hash = digest(m);
        const auto r = algo::meshGeometryModel(*model.model, o, m);
        check(!r.published && r.destinationUnchanged && r.reason == expected && digest(m) == hash,
              name);
    };
    using R = algo::GeometryReason;
    algo::GeometryMeshingOptions o;
    o.targetSize = 0;
    rejects(o, R::InvalidRequest, "zero size");
    o = {};
    o.targetSize = std::numeric_limits<double>::quiet_NaN();
    rejects(o, R::InvalidRequest, "NaN size");
    o = {};
    o.threads = -1;
    rejects(o, R::InvalidRequest, "negative threads");
    o = {};
    o.relativeRecognitionTolerance = std::numeric_limits<double>::infinity();
    rejects(o, R::InvalidRequest, "bad recognition tolerance");
    o = {};
    o.strategy = static_cast<algo::GeometryStrategy>(255);
    rejects(o, R::InvalidRequest, "bad strategy enum");
    o = {};
    o.maximumAnalysisWork = 1;
    rejects(o, R::AnalysisBudget, "analysis work bound");
    o = {};
    o.maximumNodes = 1;
    rejects(o, R::OutputBudget, "node bound");
    o = {};
    o.maximumCells = 1;
    rejects(o, R::OutputBudget, "cell bound");
    o = {};
    o.maximumOutputBytes = 1;
    rejects(o, R::OutputBudget, "byte bound");
    o = {};
    o.maximumAuditFaceUses = 1;
    rejects(o, R::OutputBudget, "audit bound");
    o = {};
    o.targetSize = 1e-300;
    rejects(o, R::OutputBudget, "unrepresentable density");
    std::atomic<bool> stop{true};
    o = {};
    o.cancellation = &stop;
    rejects(o, R::Cancelled, "pre-cancelled");
    o = {};
    o.strategy = algo::GeometryStrategy::LinearPrismSweep;
    o.maximumSourcePoints = 1;
    rejects(o, R::OutputBudget, "source bound");
    o = {};
    o.strategy = algo::GeometryStrategy::LinearPrismSweep;
    o.maximumMeshingWork = 1;
    rejects(o, R::MeshingBudget, "constraint work bound");
    for (const auto strategy : {algo::GeometryStrategy::MappedHex8,
                               algo::GeometryStrategy::LinearPrismSweep}) {
        for (const auto work : {0ULL, 1ULL}) {
            o = {};
            o.strategy = strategy;
            o.maximumMeshingWork = work;
            rejects(o, R::MeshingBudget, "structured construction work floor");
        }
    }
    o = {};
    o.targetSize = .5;
    o.strategy = algo::GeometryStrategy::MappedHex8;
    mesh::Mesh counted;
    const auto countedResult = algo::meshGeometryModel(*model.model, o, counted);
    check(countedResult.published, "construction work floor reference");
    const auto constructionCount = static_cast<std::uint64_t>(counted.nodeCount()) +
                                   static_cast<std::uint64_t>(counted.cellCount());
    o.maximumMeshingWork = constructionCount - 1;
    rejects(o, R::MeshingBudget, "exact mapped construction floor rollback");
    o.maximumMeshingWork = constructionCount;
    auto atFloor = held();
    check(algo::meshGeometryModel(*model.model, o, atFloor).published &&
              digest(atFloor) == digest(counted),
          "exact mapped construction floor accepted");
    o = {};
    o.minimumElementQuality = 1;
    rejects(o, R::QualityRequirementNotMet, "quality threshold rollback");
    o = {};
    o.minimumElementQuality = std::numeric_limits<double>::quiet_NaN();
    rejects(o, R::InvalidRequest, "NaN quality threshold");
    auto triangle = geo::makePolyhedralModel(cases()[9].input);
    o = {};
    o.strategy = algo::GeometryStrategy::MappedHex8;
    auto m = held();
    check(algo::meshGeometryModel(*triangle.model, o, m).reason == R::NoMatchingStrategy &&
              m.node(0).x == 17,
          "unproved forced strategy");
    auto input = cases()[9].input;
    for (auto &v : input.vertices)
        v.x += 4;
    auto separated = geo::makePolyhedralModel(input);
    CombinedModel combined;
    combined.add(*model.model, 0);
    combined.add(*separated.model, 1000);
    o = {};
    o.targetSize = .5;
    mesh::Mesh mixed;
    auto r = algo::meshGeometryModel(combined, o, mixed);
    check(r.published && r.audit.valid && r.analysis.bodies.size() == 2, "mixed separated bodies");
    check(r.analysis.bodies[0].strategy == algo::GeometryStrategy::MappedHex8 &&
              r.analysis.bodies[1].strategy == algo::GeometryStrategy::LinearPrismSweep,
          "per-body selection");
    check(std::abs(positiveVolume(mixed) - 2) < 1e-9, "mixed volume");
    checkClassification(combined, mixed, 1e-7);
    combined.reversed = true;
    mesh::Mesh reversed;
    check(algo::meshGeometryModel(combined, o, reversed).published &&
              digest(mixed) == digest(reversed),
          "enumeration determinism");
    o.maximumOutputBytes = r.outputBytes - 1;
    m = held();
    check(!algo::meshGeometryModel(combined, o, m).published && m.nodeCount() == 1 &&
              m.node(0).x == 17,
          "aggregate bound rollback");
    CombinedModel touching;
    touching.add(*model.model, 0);
    touching.add(*triangle.model, 1000);
    o = {};
    m = held();
    check(algo::meshGeometryModel(touching, o, m).reason == R::UnresolvedBodyInterface &&
              m.node(0).x == 17,
          "shared interface refusal");
}
void audits() {
    auto built = geo::makePolyhedralModel(cases().front().input);
    algo::GeometryMeshingOptions o;
    o.targetSize = 2;
    mesh::Mesh valid;
    check(algo::meshGeometryModel(*built.model, o, valid).published, "audit fixture");
    check(mesh::auditLinearVolumeBoundary(valid).valid, "valid hex boundary");
    mesh::Mesh noBoundary;
    for (Index i = 0; i < valid.nodeCount(); ++i)
        (void)noBoundary.addNode(valid.node(i), valid.nodeClass(i));
    for (Index i = 0; i < valid.cellCount(); ++i)
        if (valid.cell(i).type == mesh::CellType::Hex8)
            (void)noBoundary.addCell(valid.cell(i).type, {valid.cellNodes(i), 8},
                                     valid.cell(i).onEntity);
    check(mesh::auditLinearVolumeBoundary(noBoundary).missingBoundary == 6, "missing exterior");
    auto duplicate = valid;
    for (Index i = 0; i < valid.cellCount(); ++i)
        if (valid.cell(i).type == mesh::CellType::Quad4) {
            (void)duplicate.addCell(mesh::CellType::Quad4, {valid.cellNodes(i), 4},
                                    valid.cell(i).onEntity);
            break;
        }
    check(mesh::auditLinearVolumeBoundary(duplicate).duplicateBoundary == 1, "duplicate exterior");
    auto nonfinite = valid;
    nonfinite.setNode(0, {std::numeric_limits<double>::quiet_NaN(), 0, 0});
    check(mesh::auditLinearVolumeBoundary(nonfinite).nonfiniteNodes == 1, "nonfinite audit");
    auto invalid = valid;
    (void)invalid.addCell(mesh::CellType::Tet4, {0, 1, 2, 2}, {3, 1});
    check(mesh::auditLinearVolumeBoundary(invalid).invalidConnectivity == 1, "connectivity audit");
    check(mesh::auditLinearVolumeBoundary(valid, {1}).budgetExceeded, "face storage bound");
    mesh::Mesh tet;
    for (const auto p : std::array<Vec3, 4>{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}})
        (void)tet.addNode(p, {3, 1});
    (void)tet.addCell(mesh::CellType::Tet4, {0, 2, 1, 3}, {3, 1});
    const auto quality = algo::analyzeMixedMeshQuality(tet);
    check(!quality.valid && quality.invertedCells == 1,
          "inverted tet must not become positive by absolute value");
}
void sweep() {
    mesh::Mesh source;
    for (const auto p : std::array<Vec3, 3>{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}})
        (void)source.addNode(p, {2, 1});
    (void)source.addCell(mesh::CellType::Tri3, {0, 1, 2}, {2, 1});
    algo::PreparedMixedSweepRequest r;
    r.targetNodeClassifications = {{2, 2}};
    auto m = held();
    const auto hash = digest(m);
    check(!algo::generatePreparedMixedSweep(source, r, m).published() && digest(m) == hash,
          "implicit target array bounds");
    r = {};
    r.intermediateNodeClassifications = {{3, 1}};
    check(!algo::generatePreparedMixedSweep(source, r, m).published() && digest(m) == hash,
          "middle array bounds");
    r = {};
    r.targetNodeClassifications = {{9, 1}, {9, 1}, {9, 1}};
    check(!algo::generatePreparedMixedSweep(source, r, m).published(), "invalid target classes");
    r = {};
    r.sideEdgeClassifications = {{0, 1, {2, 3}}, {0, 1, {2, 4}}, {2, 0, {2, 5}}};
    check(!algo::generatePreparedMixedSweep(source, r, m).published() && digest(m) == hash,
          "duplicate side map");
    r = {};
    check(algo::generatePreparedMixedSweep(source, r, m).published(), "ordinary sweep retained");
}
} // namespace
int main() {
    try {
        polyhedral();
        options();
        audits();
        sweep();
        std::cout << checks << " geometry hardening checks passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
