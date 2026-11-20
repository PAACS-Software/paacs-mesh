// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "audit_reference.hpp"
#include "paacs/algo/AdaptiveBoxMesher.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "primitive_size_checks.hpp"
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>

using namespace paacs;
using primitive_test::require;
namespace {
struct Stats {
    double adjacent = 1, sizeRatio = 0, quality = 1, meshVolume = 0;
    std::uint64_t interfaces = 0;
};
Stats inspect(const mesh::Mesh &m, const sizing::RegionSizeField &field) {
    Stats out;
    struct Face {
        double volume = 0;
        unsigned uses = 0;
    };
    std::map<std::array<Index, 3>, Face> faces;
    long double total = 0, remainder = 0;
    for (Index i = 0; i < m.cellCount(); ++i)
        if (m.cell(i).type == mesh::CellType::Tet4) {
            const auto *v = m.cellNodes(i);
            const auto a = m.node(v[0]), b = m.node(v[1]), c = m.node(v[2]), d = m.node(v[3]);
            const auto x = b - a, y = c - a, z = d - a;
            const long double det =
                static_cast<long double>(x.x) *
                    (static_cast<long double>(y.y) * z.z - static_cast<long double>(y.z) * z.y) -
                static_cast<long double>(x.y) *
                    (static_cast<long double>(y.x) * z.z - static_cast<long double>(y.z) * z.x) +
                static_cast<long double>(x.z) *
                    (static_cast<long double>(y.x) * z.y - static_cast<long double>(y.y) * z.x);
            require(det > 0, "independent positive tet");
            const long double measuredVolume = det / 6;
            const long double next = total + measuredVolume;
            // Independent raw-coordinate determinants and a Neumaier sum;
            // long double alone does not provide extra precision on MSVC.
            remainder += std::abs(total) >= std::abs(measuredVolume)
                             ? (total - next) + measuredVolume
                             : (measuredVolume - next) + total;
            total = next;
            const double volume = static_cast<double>(measuredVolume);
            double edge = 0;
            for (unsigned j = 0; j < 4; ++j)
                for (unsigned k = 0; k < j; ++k)
                    edge = std::max(edge, (m.node(v[j]) - m.node(v[k])).norm());
            const double target = field.eval((a + b + c + d) / 4);
            out.sizeRatio = std::max(out.sizeRatio, edge / target);
            for (unsigned omitted = 0; omitted < 4; ++omitted) {
                std::array<Index, 3> f{};
                unsigned at = 0;
                for (unsigned j = 0; j < 4; ++j)
                    if (j != omitted)
                        f[at++] = v[j];
                std::sort(f.begin(), f.end());
                auto &row = faces[f];
                if (row.uses) {
                    out.adjacent =
                        std::max(out.adjacent,
                                 std::cbrt(std::max(volume / row.volume, row.volume / volume)));
                    ++out.interfaces;
                }
                row.volume = volume;
                require(++row.uses <= 2, "independent manifold faces");
            }
        }
    out.meshVolume = static_cast<double>(total + remainder);
    out.quality = primitive_test::sizes(m).minimumQuality;
    require(out.sizeRatio <= 1 + 1e-10, "actual edge-to-field bound");
    const auto audit = mesh::referenceLinearVolumeBoundary(m);
    require(audit.complete && audit.valid, "independent boundary oracle");
    return out;
}
sizing::RegionSizeSpecification specification(std::string name, double growth = 1.3) {
    sizing::RegionSizeSpecification s;
    s.baseSize = .25;
    if (name == "uniform")
        return s;
    sizing::RegionSizeControl c;
    c.targetSize = .0625;
    c.growthRatio = growth;
    if (name == "cube") {
        c.kind = sizing::SizeRegionKind::Box;
        c.lower = {.375, .375, .375};
        c.upper = {.625, .625, .625};
    } else if (name == "sphere") {
        c.kind = sizing::SizeRegionKind::Ball;
        c.lower = {.5, .5, .5};
        c.radius = .16;
    } else if (name == "face")
        c = sizing::boxFaceSizeControl({0, 0, 0}, {1, 1, 1}, sizing::BoxFace::XMin, .0625, growth);
    else
        throw std::runtime_error("unknown case");
    s.controls.push_back(c);
    return s;
}
} // namespace
int main(int argc, char **argv) {
    try {
        std::filesystem::path gallery;
        if (argc == 3 && std::string(argv[1]) == "--gallery")
            gallery = argv[2];
        else
            require(argc == 1, "usage: adaptive-sizing [--gallery directory]");
        if (!gallery.empty())
            std::filesystem::create_directories(gallery);
        unsigned cases = 0, refusals = 0;
        auto rollback = [&](const algo::AdaptiveBoxRequest &q, algo::AdaptiveBoxReason expected) {
            mesh::Mesh old;
            old.addNode({9, 8, 7}, {0, 9});
            const auto hash = mesh::meshHash(old);
            const auto result = algo::meshAdaptiveBox(q, old);
            require(!result.published && result.destinationUnchanged && result.reason == expected &&
                        mesh::meshHash(old) == hash,
                    std::string("adaptive rollback: ") +
                        algo::adaptiveBoxReasonName(result.reason));
            ++refusals;
        };
        for (const auto *name : {"uniform", "cube", "sphere", "face"})
            for (double g : {1.15, 1.3, 1.7}) {
                if (std::string(name) == "uniform" && g != 1.3)
                    continue;
                const auto f = sizing::compileRegionSizeField(specification(name, g));
                require(f.compiled, "compile field");
                algo::AdaptiveBoxRequest q;
                q.field = &f.field;
                mesh::Mesh m;
                const auto r = algo::meshAdaptiveBox(q, m);
                require(r.published,
                        std::string("mesh ") + name + " " + algo::adaptiveBoxReasonName(r.reason));
                const auto stats = inspect(m, f.field);
                require(r.audit.volumeGrowthMeasured &&
                            std::abs(stats.adjacent - r.audit.maximumAdjacentSizeRatio) < 1e-10,
                        "independent growth audit");
                require(std::abs(stats.quality - r.minimumQuality) < 1e-10,
                        "independent mean ratio");
                require(std::abs(r.meshVolume - 1) < 1e-12, "cube volume");
                require(std::abs(stats.meshVolume - 1) < 1e-12 &&
                            std::abs(stats.meshVolume - r.meshVolume) < 1e-12,
                        "independent cube volume and reported volume");
                q.threads = 4;
                mesh::Mesh repeat;
                require(algo::meshAdaptiveBox(q, repeat).published &&
                            mesh::meshHash(m) == mesh::meshHash(repeat),
                        "thread determinism");
                std::cout << std::setprecision(17) << "{\"case\":\"" << name
                          << "\",\"growth\":" << g << ",\"leaves\":" << r.leaves
                          << ",\"nodes\":" << r.nodes << ",\"tets\":" << r.tetrahedra
                          << ",\"quality\":" << r.minimumQuality
                          << ",\"meshVolume\":" << r.meshVolume
                          << ",\"adjacentSizeRatio\":" << stats.adjacent
                          << ",\"maximumEdgeToField\":" << stats.sizeRatio
                          << ",\"workingBytes\":" << r.workingBytes
                          << ",\"outputBytes\":" << r.outputBytes << ",\"seconds\":"
                          << r.planningSeconds + r.generationSeconds + r.validationSeconds +
                                 r.auditSeconds
                          << ",\"hash\":" << mesh::meshHash(m) << "}\n";
                if (!gallery.empty()) {
                    const auto file =
                        gallery /
                        (std::string(name) + "-" +
                         std::to_string(static_cast<int>(std::round(g * 100))) + ".paacsmesh");
                    require(sdk::writeMeshArchive(m, file.string()).completed, "export");
                }
                ++cases;
            }
        // Analytic field checks, including control-order independence and a
        // sub-grid ball whose interior would be missed by sparse sampling.
        unsigned fieldChecks = 0;
        for (const auto *name : {"cube", "sphere", "face"}) {
            auto spec = specification(name);
            const auto compiled = sizing::compileRegionSizeField(spec);
            require(compiled.compiled, "analytic compile");
            const auto &f = compiled.field;
            const Vec3 inner = std::string(name) == "face" ? Vec3{0, .5, .5} : Vec3{.5, .5, .5};
            require(f.eval(inner) == .0625, "inner target plateau");
            std::mt19937_64 rng(77);
            std::uniform_real_distribution<double> u(0, 1);
            for (unsigned i = 0; i < 500; ++i) {
                const Vec3 a{u(rng), u(rng), u(rng)}, b{u(rng), u(rng), u(rng)};
                require(std::abs(f.eval(a) - f.eval(b)) <=
                            f.gradientBound() * (a - b).norm() + 1e-14,
                        "physical Lipschitz bound");
            }
            auto overlap = specification("sphere");
            overlap.controls.push_back(specification("face").controls.front());
            const auto a = sizing::compileRegionSizeField(overlap);
            std::reverse(overlap.controls.begin(), overlap.controls.end());
            const auto b = sizing::compileRegionSizeField(overlap);
            for (unsigned i = 0; i < 100; ++i) {
                const Vec3 p{u(rng), u(rng), u(rng)};
                require(a.field.eval(p) == b.field.eval(p), "permutation invariance");
            }
            ++fieldChecks;
        }
        {
            auto spec = specification("sphere");
            spec.controls[0].lower = {.431, .527, .619};
            spec.controls[0].radius = 1e-8;
            auto f = sizing::compileRegionSizeField(spec);
            require(f.compiled && f.field.minimumOnBox({0, 0, 0}, {1, 1, 1}) == .0625,
                    "tiny off-grid region");
            spec.minimumSize = .1;
            f = sizing::compileRegionSizeField(spec);
            require(f.field.eval(spec.controls[0].lower) == .1, "minimum-size floor");
            spec.minimumSize = 0;
            spec.controls[0].growthRatio = 1;
            f = sizing::compileRegionSizeField(spec);
            require(f.field.eval({1e100, 1e100, 1e100}) == .0625,
                    "growth one has global influence");
            ++fieldChecks;
        }
        auto fieldReject = [&](const sizing::RegionSizeSpecification &spec,
                               sizing::RegionSizeReason why) {
            const auto r = sizing::compileRegionSizeField(spec);
            require(!r.compiled && r.reason == why, "field refusal");
            ++refusals;
        };
        auto s = specification("cube");
        s.baseSize = 0;
        fieldReject(s, sizing::RegionSizeReason::InvalidSize);
        s = specification("cube");
        s.minimumSize = 1;
        fieldReject(s, sizing::RegionSizeReason::InvalidSize);
        s = specification("cube");
        s.controls[0].growthRatio = .9;
        fieldReject(s, sizing::RegionSizeReason::InvalidGrowth);
        s = specification("cube");
        s.controls[0].lower = s.controls[0].upper;
        fieldReject(s, sizing::RegionSizeReason::InvalidRegion);
        s = specification("sphere");
        s.controls[0].radius = -1;
        fieldReject(s, sizing::RegionSizeReason::InvalidRegion);
        s = specification("face");
        s.controls[0].upper = s.controls[0].lower;
        fieldReject(s, sizing::RegionSizeReason::InvalidRegion);
        s = specification("face");
        s.controls[0].kind = static_cast<sizing::SizeRegionKind>(255);
        fieldReject(s, sizing::RegionSizeReason::InvalidRegion);
        s = specification("cube");
        s.maximumControls = 0;
        fieldReject(s, sizing::RegionSizeReason::BudgetExceeded);
        s = specification("cube");
        s.maximumBytes = sizeof(sizing::RegionSizeControl) - 1;
        fieldReject(s, sizing::RegionSizeReason::BudgetExceeded);
        s = specification("cube");
        s.controls[0].targetSize = std::numeric_limits<double>::quiet_NaN();
        fieldReject(s, sizing::RegionSizeReason::InvalidSize);
        s = specification("cube");
        s.controls[0].growthRatio = std::numeric_limits<double>::infinity();
        fieldReject(s, sizing::RegionSizeReason::InvalidGrowth);
        {
            const auto f = sizing::compileRegionSizeField(specification("cube"));
            algo::AdaptiveBoxRequest q;
            q.field = &f.field;
            mesh::Mesh m;
            const auto r = algo::meshAdaptiveBox(q, m);
            require(r.published, "budget baseline");
            q.maximumNodes = r.nodes;
            q.maximumCells = r.tetrahedra + r.boundaryTriangles;
            q.maximumLeaves = r.leaves;
            q.maximumOutputBytes = r.outputBytes;
            q.maximumWorkingBytes = r.workingBytes;
            q.auditLimits.maximumFaceUses = r.audit.faceUses;
            q.auditLimits.maximumTemporaryBytes = r.audit.temporaryBytes;
            q.maximumAdjacentSizeRatio = 1.3;
            mesh::Mesh exact;
            require(algo::meshAdaptiveBox(q, exact).published &&
                        mesh::meshHash(m) == mesh::meshHash(exact),
                    "exact memory/count/audit budgets");
            using R = algo::AdaptiveBoxReason;
            auto bad = q;
            --bad.maximumNodes;
            rollback(bad, R::OutputBudget);
            bad = q;
            --bad.maximumCells;
            rollback(bad, R::OutputBudget);
            bad = q;
            --bad.maximumLeaves;
            rollback(bad, R::LeafBudget);
            bad = q;
            --bad.maximumOutputBytes;
            rollback(bad, R::OutputBudget);
            bad = q;
            --bad.maximumWorkingBytes;
            rollback(bad, R::WorkingBudget);
            bad = q;
            --bad.auditLimits.maximumFaceUses;
            rollback(bad, R::WorkingBudget);
            bad = q;
            --bad.auditLimits.maximumTemporaryBytes;
            rollback(bad, R::WorkingBudget);
            bad = q;
            bad.maximumAdjacentSizeRatio = 1.15;
            rollback(bad, R::GrowthRequirement);
            bad = q;
            bad.minimumQuality = .7;
            rollback(bad, R::QualityRequirement);
            bad = q;
            bad.maximumDepth = 1;
            rollback(bad, R::DepthLimit);
            bad = q;
            bad.maximumDepth = 21;
            rollback(bad, R::InvalidRequest);
            bad = q;
            bad.lower = bad.upper;
            rollback(bad, R::InvalidRequest);
            bad = q;
            bad.lower.x = std::numeric_limits<double>::infinity();
            rollback(bad, R::InvalidRequest);
            bad = q;
            bad.targetSize = 0;
            rollback(bad, R::InvalidRequest);
            bad = q;
            bad.threads = -1;
            rollback(bad, R::InvalidRequest);
            bad = q;
            bad.maximumAdjacentSizeRatio = .99;
            rollback(bad, R::InvalidRequest);
            bad = q;
            bad.minimumQuality = std::numeric_limits<double>::quiet_NaN();
            rollback(bad, R::InvalidRequest);
            std::atomic<bool> stop{true};
            bad = q;
            bad.cancellation = &stop;
            rollback(bad, R::Cancelled);
            bad = q;
            bad.auditLimits.cancellation = &stop;
            rollback(bad, R::Cancelled);
            sizing::RegionSizeField invalid;
            bad = q;
            bad.field = &invalid;
            rollback(bad, R::InvalidRequest);
            // The new growth audit rejects missing/nonpositive volume data;
            // it never publishes a meaningless ratio as a successful check.
            std::vector<double> volumes(1, 1);
            mesh::LinearVolumeAuditLimits limits;
            limits.cellVolumes = volumes;
            auto audited = mesh::auditLinearVolumeBoundary(m, limits);
            require(audited.complete && !audited.valid && audited.invalidVolumeData,
                    "volume span mismatch");
            volumes.assign(static_cast<std::size_t>(m.cellCount()), 1);
            volumes[0] = 0;
            limits.cellVolumes = volumes;
            audited = mesh::auditLinearVolumeBoundary(m, limits);
            require(!audited.valid && audited.invalidVolumeData, "zero supplied volume");
            ++cases;
        }
        for (double scale : {1e-6, 1., 1e6})
            for (unsigned face = 0; face < 6; ++face) {
                sizing::RegionSizeSpecification spec;
                spec.baseSize = .5 * scale;
                const Vec3 lo = Vec3{2, -3, 4} * scale, hi = lo + Vec3{scale, scale, scale};
                spec.controls.push_back(sizing::boxFaceSizeControl(
                    lo, hi, static_cast<sizing::BoxFace>(face), .125 * scale, 1.7));
                const auto f = sizing::compileRegionSizeField(spec);
                require(f.compiled, "scaled face field");
                algo::AdaptiveBoxRequest q;
                q.lower = lo;
                q.upper = hi;
                q.field = &f.field;
                mesh::Mesh m;
                const auto r = algo::meshAdaptiveBox(q, m);
                require(r.published, "all six scaled faces");
                const auto stats = inspect(m, f.field);
                require(stats.quality > .55 && stats.adjacent < 1.2601,
                        "scaled quality and growth");
                ++cases;
            }
        std::mt19937_64 rng(0x73697a696e67);
        std::uniform_real_distribution<double> u(0, 1);
        for (unsigned i = 0; i < 12; ++i) {
            sizing::RegionSizeSpecification spec;
            spec.baseSize = .3;
            sizing::RegionSizeControl c;
            c.kind = i % 2 ? sizing::SizeRegionKind::Box : sizing::SizeRegionKind::Ball;
            c.lower = {.1 + .5 * u(rng), .1 + .5 * u(rng), .1 + .5 * u(rng)};
            c.upper = c.lower + Vec3{.05, .09, .12};
            c.radius = .03 + .1 * u(rng);
            c.targetSize = .055 + .08 * u(rng);
            c.growthRatio = 1.4 + u(rng);
            spec.controls.push_back(c);
            const auto f = sizing::compileRegionSizeField(spec);
            algo::AdaptiveBoxRequest q;
            q.field = &f.field;
            mesh::Mesh m;
            const auto r = algo::meshAdaptiveBox(q, m);
            require(r.published, "seeded local controls");
            const auto stats = inspect(m, f.field);
            require(stats.quality > .55 && stats.adjacent < 1.2601, "seeded quality and growth");
            ++cases;
        }
        {
            auto spec = specification("sphere", 100);
            spec.baseSize = 1;
            spec.controls[0].lower = {.213, .437, .591};
            spec.controls[0].radius = .005;
            spec.controls[0].targetSize = 1. / 64;
            const auto f = sizing::compileRegionSizeField(spec);
            algo::AdaptiveBoxRequest q;
            q.field = &f.field;
            mesh::Mesh m;
            const auto r = algo::meshAdaptiveBox(q, m);
            require(r.published, "steep field balancing");
            const auto stats = inspect(m, f.field);
            require(r.balanceSplits > 0 && stats.quality > .55 && stats.adjacent < 1.2601,
                    "strong balance exercised");
            ++cases;
            std::cout << "{\"balanceSplits\":" << r.balanceSplits
                      << ",\"balancePasses\":" << r.balancePasses << "}\n";
        }
        std::cout << "{\"passed\":" << cases << ",\"fieldChecks\":" << fieldChecks
                  << ",\"refusals\":" << refusals << "}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
