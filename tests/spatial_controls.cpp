// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "audit_reference.hpp"
#include "geometry_checks.hpp"
#include "geometry_model_views.hpp"
#include "paacs/algo/AdaptiveBoxMesher.hpp"
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/geometry/PolyhedralModel.hpp"
#include "spatial_control_fixture.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
using namespace paacs;
using geometry_test::require;
namespace {
unsigned configurations = 0, refusals = 0;
std::filesystem::path gallery;
void save(const std::string &name, const mesh::Mesh &m, const sizing::BoundedSizeField &f) {
    if (gallery.empty())
        return;
    std::filesystem::create_directories(gallery);
    require(sdk::writeMeshArchive(m, (gallery / (name + ".paacsmesh")).string()).completed,
            "save archive");
    std::ofstream csv(gallery / (name + "-field.csv"));
    csv << "x,z,size\n" << std::setprecision(17);
    for (int j = 0; j <= 100; ++j)
        for (int i = 0; i <= 100; ++i)
            csv << i / 100. << ',' << j / 100. << ',' << f.eval({i / 100., .5, j / 100.}) << '\n';
    require(bool(csv), "save exact field samples");
}
void boundaryCheck(const mesh::Mesh &m, double scale) {
    std::set<Index> vertices;
    std::set<std::pair<Index, Index>> edges;
    std::uint64_t faces = 0;
    for (Index i = 0; i < m.cellCount(); ++i)
        if (m.cell(i).type == mesh::CellType::Tri3) {
            const auto &cell = m.cell(i);
            const auto *ids = m.cellNodes(i);
            require(cell.onEntity.dim == 2 && cell.onEntity.tag >= 1 && cell.onEntity.tag <= 6,
                    "box face tag");
            const int axis = (cell.onEntity.tag - 1) / 2;
            const double position = cell.onEntity.tag % 2 ? 0 : scale;
            for (int j = 0; j < 3; ++j) {
                const auto p = m.node(ids[j]);
                const double value = axis == 0 ? p.x : axis == 1 ? p.y : p.z;
                require(std::abs(value - position) <= scale * 1e-12, "face location");
                vertices.insert(ids[j]);
                edges.insert(std::minmax(ids[j], ids[(j + 1) % 3]));
            }
            ++faces;
        }
    require(vertices.size() + faces == edges.size() + 2, "boundary Euler characteristic");
}
// An independent dense check within each Tet4, not only generator leaf queries.
void interiorSizes(const mesh::Mesh &m, const sizing::BoundedSizeField &f) {
    for (Index i = 0; i < m.cellCount(); ++i)
        if (m.cell(i).type == mesh::CellType::Tet4) {
            const auto *ids = m.cellNodes(i);
            std::array<Vec3, 4> p{};
            double longest = 0;
            for (int j = 0; j < 4; ++j)
                p[j] = m.node(ids[j]);
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b) {
                    const auto d = p[a] - p[b];
                    longest = std::max(longest, std::hypot(d.x, d.y, d.z));
                }
            for (int a = 0; a <= 3; ++a)
                for (int b = 0; b <= 3 - a; ++b)
                    for (int c = 0; c <= 3 - a - b; ++c) {
                        const auto q =
                            (p[0] * a + p[1] * b + p[2] * c + p[3] * (3 - a - b - c)) * (1. / 3);
                        require(longest <= f.eval(q) * (1 + 1e-10), "interior target respected");
                    }
        }
}
void adaptive(const std::string &name, double scale = 1) {
    auto f = spatial_test::make(name, scale);
    algo::AdaptiveBoxRequest o;
    o.upper = {scale, scale, scale};
    o.field = &f;
    o.targetSize = .4 * scale;
    o.threads = 1;
    o.maximumSizingWork = 50'000'000;
    mesh::Mesh m;
    const auto start = std::chrono::steady_clock::now();
    const auto r = algo::meshAdaptiveBox(o, m);
    require(r.published,
            std::string("adaptive ") + name + ": " + algo::adaptiveBoxReasonName(r.reason));
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    require(mesh::referenceLinearVolumeBoundary(m).valid, "independent full boundary audit");
    require(std::abs(geometry_test::positiveVolume(m) / (scale * scale * scale) - 1) < 1e-10,
            "independent volume");
    boundaryCheck(m, scale);
    interiorSizes(m, f);
    algo::MeshSizeAuditOptions limits;
    limits.maximumEdgeToTargetRatio = 1.2;
    const auto audit = algo::auditMeshSize(m, f, limits);
    require(audit.complete, "size audit complete");
    require(audit.maximumSampledRatio <= 1 + 1e-10, "size diagnostic oracle agreement");
    o.threads = 4;
    mesh::Mesh parallel;
    const auto again = algo::meshAdaptiveBox(o, parallel);
    require(again.published && geometry_test::digest(m) == geometry_test::digest(parallel),
            "thread-identical mesh");
    if (scale == 1)
        save(name, m, f);
    std::cout << std::setprecision(17) << "{\"case\":\"" << name << "\",\"scale\":" << scale
              << ",\"nodes\":" << m.nodeCount() << ",\"tetrahedra\":" << r.tetrahedra
              << ",\"output_bytes\":" << r.outputBytes << ",\"working_bytes\":" << r.workingBytes
              << ",\"sizing_queries\":" << r.sizingQueries << ",\"sizing_work\":" << r.sizingWork
              << ",\"field_nodes\":" << f.nodeCount() << ",\"field_bytes\":" << f.logicalBytes()
              << ",\"size_ratio_bound\":" << audit.maximumConservativeRatio
              << ",\"size_ratio_sampled\":" << audit.maximumSampledRatio
              << ",\"minimum_quality\":" << r.minimumQuality << ",\"seconds\":" << seconds
              << ",\"hash\":\"" << geometry_test::digest(m) << "\"}\n"
              << std::flush;
    ++configurations;
    if (name == "point" && scale == 1) {
        auto reject = [&](algo::AdaptiveBoxRequest request, algo::AdaptiveBoxReason reason) {
            mesh::Mesh held;
            held.addNode({9, 8, 7}, {0, 9});
            const auto before = geometry_test::digest(held);
            const auto bad = algo::meshAdaptiveBox(request, held);
            require(!bad.published && bad.destinationUnchanged && bad.reason == reason &&
                        geometry_test::digest(held) == before,
                    "adaptive bounded refusal/rollback");
            ++refusals;
        };
        auto exact = o;
        exact.maximumSizingWork = r.sizingWork;
        mesh::Mesh dst;
        require(algo::meshAdaptiveBox(exact, dst).published, "exact sizing work admitted");
        exact.maximumSizingWork = r.sizingWork - 1;
        reject(exact, algo::AdaptiveBoxReason::SizingBudget);
        exact = o;
        exact.maximumLeaves = 1;
        reject(exact, algo::AdaptiveBoxReason::LeafBudget);
        exact = o;
        exact.maximumDepth = 1;
        reject(exact, algo::AdaptiveBoxReason::DepthLimit);
        exact = o;
        exact.maximumOutputBytes = 1;
        reject(exact, algo::AdaptiveBoxReason::OutputBudget);
        exact = o;
        exact.maximumWorkingBytes = 1;
        reject(exact, algo::AdaptiveBoxReason::WorkingBudget);
        exact = o;
        std::atomic<bool> stop = true;
        exact.cancellation = &stop;
        reject(exact, algo::AdaptiveBoxReason::Cancelled);
        auto bound = limits;
        bound.maximumSizingWork = audit.sizingWork;
        require(algo::auditMeshSize(m, f, bound).complete, "exact audit work admitted");
        bound.maximumSizingWork = audit.sizingWork - 1;
        require(algo::auditMeshSize(m, f, bound).reason ==
                    algo::MeshSizeAuditReason::BudgetExceeded,
                "audit one unit below");
        ++refusals;
        bound = limits;
        bound.cancellation = &stop;
        require(algo::auditMeshSize(m, f, bound).reason == algo::MeshSizeAuditReason::Cancelled,
                "audit cancelled");
        ++refusals;
        bound = limits;
        bound.maximumEdgeToTargetRatio = .9;
        require(algo::auditMeshSize(m, f, bound).reason ==
                    algo::MeshSizeAuditReason::InvalidRequest,
                "audit invalid ratio");
        ++refusals;
    }
}
void geometry() {
    const std::vector<std::vector<Vec2>> loops{{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
    auto model = geo::makePolyhedralModel(geo::extrudedPolygon(loops, {0, 0, 1}));
    require(model.built(), "model");
    sizing::ControlFieldSpecification s;
    s.maximumSize = .6;
    s.nodes = {sizing::field::constant(.6)};
    s.output = 0;
    auto constant = sizing::compileControlField(s);
    require(constant.compiled, "geometry constant");
    algo::GeometryMeshingOptions o;
    o.targetSize = .6;
    o.strategy = algo::GeometryStrategy::BoundaryConformingTet4;
    o.threads = 1;
    mesh::Mesh original;
    require(algo::meshGeometryModel(*model.model, o, original).published, "constant original Tet4");
    o.strategy = algo::GeometryStrategy::Auto;
    o.spatialField = &constant.field;
    mesh::Mesh same;
    auto r = algo::meshGeometryModel(*model.model, o, same);
    require(r.published && geometry_test::digest(same) == geometry_test::digest(original),
            "constant field preserves Tet4 result");
    ++configurations;
    s.nodes = {sizing::field::distancePoint({.5, .5, .5}),
               sizing::field::threshold(0, 0, .7, .3, .6)};
    s.output = 1;
    auto local = sizing::compileControlField(s);
    require(local.compiled, "geometry local");
    o.spatialField = &local.field;
    o.warningElementSizeRatio = 1;
    mesh::Mesh m;
    r = algo::meshGeometryModel(*model.model, o, m);
    require(r.published && r.analysis.bodies[0].usedSpatialField, "geometry field consumed");
    require(r.analysis.bodies[0].strategy == algo::GeometryStrategy::BoundaryConformingTet4,
            "field routes to compatible algorithm");
    require(geometry_test::digest(m) != geometry_test::digest(original),
            "field changes emitted geometry mesh");
    require(mesh::referenceLinearVolumeBoundary(m).valid &&
                std::abs(geometry_test::positiveVolume(m) - 1) < 1e-10,
            "geometry independent checks");
    geometry_test::checkClassification(*model.model, m, 1e-8);
    const auto audit = r.analysis.bodies[0].sizeAudit;
    require(audit.complete, "geometry reports achieved sizing");
    o.threads = 4;
    mesh::Mesh parallel;
    require(algo::meshGeometryModel(*model.model, o, parallel).published &&
                geometry_test::digest(m) == geometry_test::digest(parallel),
            "geometry field thread identity");
    save("geometry-point", m, local.field);
    ++configurations;
    std::cout << "{\"case\":\"geometry-point\",\"nodes\":" << m.nodeCount()
              << ",\"cells\":" << m.cellCount() << ",\"ratio\":" << audit.maximumConservativeRatio
              << ",\"warning_cells\":" << audit.cellsAboveWarning << ",\"hash\":\""
              << geometry_test::digest(m) << "\"}\n"
              << std::flush;
    auto reject = [&](algo::GeometryMeshingOptions request, algo::GeometryReason reason) {
        mesh::Mesh held;
        held.addNode({9, 8, 7}, {0, 9});
        const auto before = geometry_test::digest(held);
        auto bad = algo::meshGeometryModel(*model.model, request, held);
        require(!bad.published && bad.destinationUnchanged && bad.reason == reason &&
                    before == geometry_test::digest(held),
                "geometry field refusal rollback");
        ++refusals;
    };
    auto exact = o;
    exact.maximumElementSizeRatio = audit.maximumConservativeRatio;
    mesh::Mesh accepted;
    require(algo::meshGeometryModel(*model.model, exact, accepted).published,
            "size ratio at reported bound admitted");
    exact.maximumElementSizeRatio = 1;
    require(audit.maximumConservativeRatio > 1, "fixture has size warning");
    reject(exact, algo::GeometryReason::SizeRequirementNotMet);
    exact = o;
    exact.maximumSizingAuditWork = 1;
    reject(exact, algo::GeometryReason::MeshingBudget);
    exact = o;
    exact.strategy = algo::GeometryStrategy::MappedHex8;
    reject(exact, algo::GeometryReason::NoMatchingStrategy);
    exact = o;
    exact.strategy = algo::GeometryStrategy::LinearPrismSweep;
    reject(exact, algo::GeometryReason::NoMatchingStrategy);
    exact = o;
    algo::GeometryBodyControl body;
    body.solidTag = 1;
    body.layers = 5;
    exact.bodyControls = {body};
    reject(exact, algo::GeometryReason::NoMatchingStrategy);
    exact = o;
    exact.maximumElementSizeRatio = .5;
    reject(exact, algo::GeometryReason::InvalidRequest);
    exact = o;
    exact.spatialField = nullptr;
    mesh::Mesh uniform;
    const auto uniformResult = algo::meshGeometryModel(*model.model, exact, uniform);
    require(uniformResult.published && uniformResult.analysis.bodies[0].sizeAudit.complete,
            "uniform mapped size audit");
    ++configurations;
    exact = o;
    exact.spatialField = nullptr;
    body = {};
    body.solidTag = 1;
    body.spatialField = &local.field;
    exact.bodyControls = {body};
    mesh::Mesh scoped;
    require(algo::meshGeometryModel(*model.model, exact, scoped).published &&
                geometry_test::digest(scoped) == geometry_test::digest(m),
            "per-body field selection");
    ++configurations;
    // Opt-in refinement has a separate sampled engine target and conservative
    // publication gate. A protected-boundary blocker must remain a refusal.
    s.nodes[1] = sizing::field::threshold(0, 0, .7, .4, .6);
    auto moderate = sizing::compileControlField(s);
    require(moderate.compiled, "refinement field");
    exact = o;
    exact.spatialField = &moderate.field;
    exact.threads = 1;
    mesh::Mesh beforeRefinement;
    const auto before = algo::meshGeometryModel(*model.model, exact, beforeRefinement);
    require(before.published, "refinement baseline");
    ++configurations;
    body = {};
    body.solidTag = 1;
    body.tetRefinement.enabled = true;
    body.tetRefinement.maximumEdgeToTarget = 2;
    exact.bodyControls = {body};
    exact.maximumElementSizeRatio = 2.1;
    mesh::Mesh refined;
    const auto done = algo::meshGeometryModel(*model.model, exact, refined);
    require(done.published && done.analysis.bodies[0].hasRefinementDiagnostic &&
                done.analysis.bodies[0].refinementInsertions > 0,
            "interior refinement applied");
    require(done.analysis.bodies[0].sizeAudit.maximumConservativeRatio <
                before.analysis.bodies[0].sizeAudit.maximumConservativeRatio,
            "refinement improves measured size");
    require(mesh::referenceLinearVolumeBoundary(refined).valid &&
                std::abs(geometry_test::positiveVolume(refined) - 1) < 1e-10,
            "refined independent geometry");
    geometry_test::checkClassification(*model.model, refined, 1e-8);
    exact.threads = 4;
    mesh::Mesh refinedAgain;
    require(algo::meshGeometryModel(*model.model, exact, refinedAgain).published &&
                geometry_test::digest(refined) == geometry_test::digest(refinedAgain),
            "refinement thread identity");
    save("geometry-refined", refined, moderate.field);
    ++configurations;
    std::cout << "{\"case\":\"geometry-refined\",\"nodes\":" << refined.nodeCount()
              << ",\"cells\":" << refined.cellCount()
              << ",\"ratio\":" << done.analysis.bodies[0].sizeAudit.maximumConservativeRatio
              << ",\"insertions\":" << done.analysis.bodies[0].refinementInsertions
              << ",\"hash\":\"" << geometry_test::digest(refined) << "\"}\n";
    auto blocked = exact;
    blocked.bodyControls[0].tetRefinement.maximumEdgeToTarget = 1.5;
    reject(blocked, algo::GeometryReason::VolumeMeshingFailed);
    mesh::Mesh blockedMesh;
    const auto blockedResult = algo::meshGeometryModel(*model.model, blocked, blockedMesh);
    require(blockedResult.analysis.bodies[0].hasRefinementDiagnostic &&
                blockedResult.analysis.bodies[0].refinementReason != 0,
            "refinement blocker exposed");
    blocked = exact;
    blocked.bodyControls[0].tetRefinement.maximumWork = 1;
    reject(blocked, algo::GeometryReason::VolumeMeshingFailed);
    blocked = exact;
    blocked.bodyControls[0].tetRefinement.maximumInsertions = 0;
    reject(blocked, algo::GeometryReason::InvalidRequest);
    blocked = exact;
    blocked.bodyControls[0].tetRefinement.maximumEdgeToTarget = NAN;
    reject(blocked, algo::GeometryReason::InvalidRequest);
    // Two separated bodies exercise aggregate audit work and scoped fields.
    auto input = geo::extrudedPolygon(loops, {0, 0, 1});
    for (auto &p : input.vertices)
        p.x += 3;
    auto other = geo::makePolyhedralModel(input);
    require(other.built(), "second model");
    geometry_test::CombinedModel combined;
    combined.add(*model.model, 0);
    combined.add(*other.model, 1000);
    algo::GeometryMeshingOptions multiple;
    multiple.targetSize = .6;
    multiple.threads = 1;
    body = {};
    body.solidTag = 1001;
    body.spatialField = &constant.field;
    multiple.bodyControls = {body};
    multiple.warningElementSizeRatio = 1.5;
    mesh::Mesh joined;
    const auto joinedResult = algo::meshGeometryModel(combined, multiple, joined);
    require(joinedResult.published && joinedResult.analysis.bodies.size() == 2,
            "scoped multi-body controls");
    require(joinedResult.analysis.bodies[0].strategy == algo::GeometryStrategy::MappedHex8 &&
                joinedResult.analysis.bodies[1].usedSpatialField,
            "scoped controls preserve other strategy");
    require(std::abs(geometry_test::positiveVolume(joined) - 2) < 1e-10, "multi-body volume");
    ++configurations;
    const auto work = joinedResult.analysis.bodies[0].sizeAudit.sizingWork +
                      joinedResult.analysis.bodies[1].sizeAudit.sizingWork;
    multiple.maximumSizingAuditWork = work;
    mesh::Mesh atBudget;
    require(algo::meshGeometryModel(combined, multiple, atBudget).published,
            "aggregate exact audit budget");
    multiple.maximumSizingAuditWork = work - 1;
    const auto held = geometry_test::digest(atBudget);
    const auto rejected = algo::meshGeometryModel(combined, multiple, atBudget);
    require(!rejected.published && rejected.reason == algo::GeometryReason::MeshingBudget &&
                geometry_test::digest(atBudget) == held,
            "aggregate audit refusal rollback");
    ++refusals;
}
void auditTypes() {
    auto field = spatial_test::make("coarse");
    const std::array<mesh::CellType, 4> types{mesh::CellType::Tet4, mesh::CellType::Hex8,
                                              mesh::CellType::Prism6, mesh::CellType::Pyramid5};
    const std::array<std::vector<Vec3>, 4> points{
        {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
         {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
         {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}},
         {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {.5, .5, 1}}}};
    const std::array<double, 4> longest{std::sqrt(2.), 1., std::sqrt(2.), std::sqrt(1.5)};
    for (std::size_t i = 0; i < types.size(); ++i) {
        mesh::Mesh mesh;
        std::vector<Index> ids;
        for (auto p : points[i])
            ids.push_back(mesh.addNode(p, {3, 1}));
        mesh.addCell(types[i], ids, {3, 1});
        algo::MeshSizeAuditOptions o;
        o.warningEdgeToTargetRatio = 1;
        o.maximumEdgeToTargetRatio = 1;
        auto r = algo::auditMeshSize(mesh, field, o);
        require(r.complete && !r.requirementMet && r.cellsAboveLimit == 1 &&
                    r.cellsAboveWarning == 1,
                "linear family strict size gate");
        require(std::abs(r.maximumConservativeRatio - longest[i] / .4) < 1e-12,
                "known linear cell edge ratio");
        o.maximumEdgeToTargetRatio = r.maximumConservativeRatio;
        require(algo::auditMeshSize(mesh, field, o).requirementMet, "linear family exact ratio");
    }
    struct InvalidField final : sizing::BoundedSizeField {
        bool valid() const noexcept override { return true; }
        double eval(const Vec3 &) const override { return NAN; }
        double minimumOnBox(const Vec3 &, const Vec3 &) const noexcept override { return .4; }
    } invalid;
    mesh::Mesh mesh;
    for (auto p : points[0])
        mesh.addNode(p, {3, 1});
    mesh.addCell(mesh::CellType::Tet4, {0, 1, 2, 3}, {3, 1});
    require(algo::auditMeshSize(mesh, invalid).reason ==
                algo::MeshSizeAuditReason::FieldEvaluationFailed,
            "nonfinite field refusal");
    ++refusals;
    struct ThrowingField final : sizing::BoundedSizeField {
        bool valid() const noexcept override { return true; }
        double eval(const Vec3 &) const override { throw std::runtime_error("user field"); }
        double minimumOnBox(const Vec3 &, const Vec3 &) const noexcept override { return .4; }
    } throwing;
    require(algo::auditMeshSize(mesh, throwing).reason ==
                algo::MeshSizeAuditReason::FieldEvaluationFailed,
            "throwing field converted to diagnostic");
    ++refusals;
    mesh::Mesh quadratic;
    for (int i = 0; i < 10; ++i)
        quadratic.addNode({double(i), double(i % 3), 0}, {3, 1});
    quadratic.addCell(mesh::CellType::Tet10, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {3, 1});
    require(algo::auditMeshSize(quadratic, field).reason ==
                algo::MeshSizeAuditReason::UnsupportedCell,
            "curved sizing not falsely certified");
    ++refusals;
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--gallery")
            gallery = argv[2];
        else if (argc != 1)
            return 2;
        for (const auto &name : spatial_test::cases())
            adaptive(name);
        for (double scale : {1e-6, 1e6})
            for (const auto &name : {"point", "segment"})
                adaptive(name, scale);
        geometry();
        auditTypes();
        std::cout << "spatial controls: PASS (" << configurations << " mesh configurations, "
                  << refusals << " refusals)\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
