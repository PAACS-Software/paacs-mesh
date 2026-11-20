// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "audit_reference.hpp"
#include "paacs/algo/Delaunay2DScalable.hpp"
#include "paacs/algo/LayeredPatchMesher.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

using namespace paacs;
namespace {
void require(bool v, const char *message) {
    if (!v)
        throw std::runtime_error(message);
}
mesh::Mesh patch(std::string kind, unsigned n) {
    mesh::Mesh m;
    std::vector<Vec2> points;
    std::mt19937_64 rng(0x707269736d);
    std::uniform_real_distribution<double> jitter(-.28, .28);
    for (unsigned j = 0; j <= n; ++j)
        for (unsigned i = 0; i <= n; ++i) {
            const bool interior = i && j && i < n && j < n;
            const double dx = kind == "unstructured" && interior ? jitter(rng) : 0;
            const double dy = kind == "unstructured" && interior ? jitter(rng) : 0;
            points.push_back({(i + dx) / n, (j + dy) / n});
            m.addNode({points.back().u, points.back().v, 0}, {2, 1});
        }
    if (kind == "unstructured") {
        for (const auto &t : algo::delaunay2DScalable(points))
            m.addCell(mesh::CellType::Tri3, std::span<const Index>(t.v), {2, 1});
    } else
        for (unsigned j = 0; j < n; ++j)
            for (unsigned i = 0; i < n; ++i) {
                const Index a = static_cast<Index>(j) * (n + 1) + i, b = a + 1, d = a + n + 1,
                            c = d + 1;
                if (kind == "hex" || (kind == "mixed" && i < n / 2))
                    m.addCell(mesh::CellType::Quad4, {a, b, c, d}, {2, 1});
                else {
                    m.addCell(mesh::CellType::Tri3, {a, b, c}, {2, 1});
                    m.addCell(mesh::CellType::Tri3, {a, c, d}, {2, 1});
                }
            }
    return m;
}
void inspect(const mesh::Mesh &source, const mesh::Mesh &m, const algo::LayeredPatchRequest &q,
             const algo::LayeredPatchResult &r) {
    require(r.schedule.planned && r.sweep.metrics.usedExplicitLayerFractions, "explicit schedule");
    const auto audit = mesh::referenceLinearVolumeBoundary(m);
    require(audit.valid && audit.complete, "independent mixed face audit");
    const auto normal = q.direction / std::hypot(q.direction.x, q.direction.y, q.direction.z);
    double previous = 0;
    for (std::size_t k = 1; k < r.schedule.offsets.size(); ++k) {
        double minimum = 1e300, maximum = 0;
        for (Index j = 0; j < source.nodeCount(); ++j) {
            const auto a = m.node(static_cast<Index>(k - 1) * source.nodeCount() + j),
                       b = m.node(static_cast<Index>(k) * source.nodeCount() + j);
            const double h = (b - a).dot(normal);
            require(h > 0, "positive actual layer");
            minimum = std::min(minimum, h);
            maximum = std::max(maximum, h);
        }
        require(maximum - minimum < 1e-11 * r.schedule.totalDepth, "constant translated layers");
        const double expected =
            k <= q.schedule.wallLayers
                ? q.schedule.firstLayerThickness *
                      std::pow(q.schedule.growthRatio, static_cast<double>(k - 1))
                : r.schedule.offsets[k] - r.schedule.offsets[k - 1];
        require(std::abs(minimum - expected) < 1e-10 * expected,
                "actual thickness matches prescribed wall/core");
        if (previous)
            require(std::max(minimum / previous, previous / minimum) <=
                        q.schedule.maximumAllowedGrowth * (1 + 1e-10),
                    "actual adjacent height ratio");
        previous = minimum;
    }
    require(std::abs(r.sweep.metrics.totalVolume - r.schedule.totalDepth) <
                1e-10 * r.schedule.totalDepth,
            "independent unit patch volume");
}
} // namespace
int main(int argc, char **argv) {
    try {
        std::filesystem::path gallery;
        if (argc == 3 && std::string(argv[1]) == "--gallery")
            gallery = argv[2];
        else
            require(argc == 1, "usage: layered-patch [--gallery directory]");
        if (!gallery.empty())
            std::filesystem::create_directories(gallery);
        unsigned cases = 0, refusals = 0;
        auto rollback = [&](const mesh::Mesh &source, const algo::LayeredPatchRequest &q,
                            algo::LayeredPatchReason expected) {
            mesh::Mesh m;
            m.addNode({9, 8, 7}, {0, 9});
            const auto hash = mesh::meshHash(m), sourceHash = mesh::meshHash(source);
            const auto r = algo::meshLayeredPatch(source, q, m);
            if (r.reason != expected)
                std::cerr << "rollback expected " << int(expected) << " got " << int(r.reason)
                          << '\n';
            require(!r.published && r.destinationUnchanged && r.reason == expected &&
                        mesh::meshHash(m) == hash && mesh::meshHash(source) == sourceHash,
                    "layer rollback");
            ++refusals;
        };
        for (const auto *kind : {"structured", "unstructured", "hex", "mixed"}) {
            const auto source = patch(kind, 12);
            algo::LayeredPatchRequest q;
            mesh::Mesh m;
            auto r = algo::meshLayeredPatch(source, q, m);
            if (!r.published) {
                std::cerr << kind << ": " << static_cast<int>(r.reason) << " schedule "
                          << static_cast<int>(r.schedule.reason) << " sweep "
                          << static_cast<int>(r.sweep.reason) << '\n';
            }
            require(r.published, "layered mesh published");
            inspect(source, m, q, r);
            q.threads = 4;
            mesh::Mesh repeated;
            require(algo::meshLayeredPatch(source, q, repeated).published &&
                        mesh::meshHash(m) == mesh::meshHash(repeated),
                    "layer thread determinism");
            std::cout << std::setprecision(17) << "{\"case\":\"" << kind
                      << "\",\"sourceCells\":" << source.cellCount()
                      << ",\"nodes\":" << m.nodeCount() << ",\"wedges\":" << r.sweep.metrics.wedges
                      << ",\"hexes\":" << r.sweep.metrics.hexes
                      << ",\"wallLayers\":" << r.schedule.wallLayers
                      << ",\"coreLayers\":" << r.schedule.coreLayers
                      << ",\"firstThickness\":" << r.schedule.offsets[1]
                      << ",\"wallThickness\":" << r.schedule.wallThickness
                      << ",\"maximumGrowth\":" << r.schedule.maximumAdjacentGrowth
                      << ",\"sourcePairChecks\":" << r.sourcePairChecks
                      << ",\"outputBytes\":" << r.outputBytes
                      << ",\"workingBytes\":" << r.plannedWorkingBytes
                      << ",\"hash\":" << mesh::meshHash(m) << "}\n";
            if (!gallery.empty())
                require(sdk::writeMeshArchive(
                            m, (gallery / (std::string(kind) + ".paacsmesh")).string())
                            .completed,
                        "layer export");
            ++cases;
        }
        const auto small = patch("structured", 2);
        for (unsigned wall : {1U, 4U, 8U})
            for (double growth : {1., 1.05, 1.2, 1.5})
                for (bool core : {false, true}) {
                    algo::LayeredPatchRequest q;
                    q.schedule.wallLayers = wall;
                    q.schedule.growthRatio = growth;
                    q.schedule.maximumAllowedGrowth = growth + .05;
                    q.schedule.totalDepth = core ? 1 : 0;
                    mesh::Mesh m;
                    const auto r = algo::meshLayeredPatch(small, q, m);
                    require(r.published, "layer parameter matrix");
                    inspect(small, m, q, r);
                    ++cases;
                }
        {
            algo::LayeredPatchRequest q;
            mesh::Mesh m;
            const auto r = algo::meshLayeredPatch(small, q, m);
            require(r.published, "layer budget baseline");
            q.maximumNodes = static_cast<std::uint64_t>(m.nodeCount());
            q.maximumCells = static_cast<std::uint64_t>(m.cellCount());
            q.maximumOutputBytes = r.plannedOutputBytes;
            q.maximumWorkingBytes = r.plannedWorkingBytes;
            q.maximumSourcePairChecks = r.sourcePairChecks;
            q.maximumSourceCells = static_cast<std::uint64_t>(small.cellCount());
            q.maximumSourceNodes = static_cast<std::uint64_t>(small.nodeCount());
            q.auditLimits.maximumTemporaryBytes = r.audit.temporaryBytes;
            q.auditLimits.maximumFaceUses = r.audit.faceUses;
            mesh::Mesh repeat;
            require(algo::meshLayeredPatch(small, q, repeat).published &&
                        mesh::meshHash(m) == mesh::meshHash(repeat),
                    "exact layer budgets");
            using R = algo::LayeredPatchReason;
            auto bad = q;
            --bad.maximumNodes;
            rollback(small, bad, R::OutputBudget);
            bad = q;
            --bad.maximumCells;
            rollback(small, bad, R::OutputBudget);
            bad = q;
            --bad.maximumOutputBytes;
            rollback(small, bad, R::OutputBudget);
            bad = q;
            --bad.maximumWorkingBytes;
            rollback(small, bad, R::WorkingBudget);
            bad = q;
            --bad.maximumSourcePairChecks;
            rollback(small, bad, R::SourceBudget);
            bad = q;
            --bad.maximumSourceCells;
            rollback(small, bad, R::SourceBudget);
            bad = q;
            --bad.maximumSourceNodes;
            rollback(small, bad, R::SourceBudget);
            bad = q;
            --bad.auditLimits.maximumTemporaryBytes;
            rollback(small, bad, R::WorkingBudget);
            bad = q;
            --bad.auditLimits.maximumFaceUses;
            rollback(small, bad, R::WorkingBudget);
            bad = q;
            bad.direction = {0, 0, 0};
            rollback(small, bad, R::InvalidRequest);
            bad = q;
            bad.direction.x = std::numeric_limits<double>::quiet_NaN();
            rollback(small, bad, R::InvalidRequest);
            bad = q;
            bad.threads = -1;
            rollback(small, bad, R::InvalidRequest);
            bad = q;
            bad.schedule.wallLayers = 0;
            rollback(small, bad, R::ScheduleRejected);
            std::atomic<bool> stop{true};
            bad = q;
            bad.cancellation = &stop;
            rollback(small, bad, R::Cancelled);
            bad = q;
            bad.auditLimits.cancellation = &stop;
            rollback(small, bad, R::Cancelled);
            auto corrupt = small;
            corrupt.setNode(0, {0, 0, .1});
            rollback(corrupt, {}, R::NonplanarSource);
            corrupt = small;
            corrupt.setNode(0, {std::numeric_limits<double>::infinity(), 0, 0});
            rollback(corrupt, {}, R::InvalidSource);
            corrupt = small;
            corrupt.addNode({.4, .4, 0}, {2, 1});
            rollback(corrupt, {}, R::InvalidSource);
            corrupt = small;
            corrupt.addCell(mesh::CellType::Tri3, {0, 1, 4}, {2, 1});
            rollback(corrupt, {}, R::IntersectingSource);
            mesh::Mesh overlap;
            for (Vec3 p : std::array<Vec3, 6>{
                     {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {.1, .1, 0}, {.5, .1, 0}, {.1, .5, 0}}})
                overlap.addNode(p, {2, 1});
            overlap.addCell(mesh::CellType::Tri3, {0, 1, 2}, {2, 1});
            overlap.addCell(mesh::CellType::Tri3, {3, 4, 5}, {2, 1});
            rollback(overlap, {}, R::IntersectingSource);
            mesh::Mesh junction;
            for (Vec3 p :
                 std::array<Vec3, 5>{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {.5, 0, 0}, {1, -1, 0}}})
                junction.addNode(p, {2, 1});
            junction.addCell(mesh::CellType::Tri3, {0, 1, 2}, {2, 1});
            junction.addCell(mesh::CellType::Tri3, {3, 4, 1}, {2, 1});
            rollback(junction, {}, R::IntersectingSource);
            mesh::Mesh bowtie;
            for (Vec3 p : std::array<Vec3, 4>{{{0, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 0, 0}}})
                bowtie.addNode(p, {2, 1});
            bowtie.addCell(mesh::CellType::Quad4, {0, 1, 2, 3}, {2, 1});
            rollback(bowtie, {}, R::InvalidSource);
            ++cases;
        }
        unsigned scheduleChecks = 0;
        auto rejectSchedule = [&](const algo::LayerScheduleRequest &q,
                                  algo::LayerScheduleReason expected) {
            const auto r = algo::planLayerSchedule(q);
            require(!r.planned && r.reason == expected, "schedule refusal");
            ++refusals;
        };
        using S = algo::LayerScheduleReason;
        algo::LayerScheduleRequest q;
        q.wallLayers = 0;
        rejectSchedule(q, S::InvalidRequest);
        q = {};
        q.firstLayerThickness = 0;
        rejectSchedule(q, S::InvalidRequest);
        q = {};
        q.growthRatio = .9;
        rejectSchedule(q, S::InvalidRequest);
        q = {};
        q.maximumAllowedGrowth = 1.1;
        rejectSchedule(q, S::InvalidRequest);
        q = {};
        q.totalDepth = -1;
        rejectSchedule(q, S::InvalidRequest);
        q = {};
        q.maximumLayers = 4;
        rejectSchedule(q, S::LayerBudget);
        q = {};
        q.maximumBytes = 0;
        rejectSchedule(q, S::ByteBudget);
        q = {};
        q.totalDepth = .05;
        rejectSchedule(q, S::WallLayersExceedDepth);
        q = {};
        q.totalDepth = .0745;
        rejectSchedule(q, S::TransitionDoesNotFit);
        q = {};
        q.growthRatio = std::numeric_limits<double>::infinity();
        rejectSchedule(q, S::InvalidRequest);
        q = {};
        q.totalDepth = std::numeric_limits<double>::quiet_NaN();
        rejectSchedule(q, S::InvalidRequest);
        q = {};
        q.firstLayerThickness = 1e308;
        q.wallLayers = 10;
        rejectSchedule(q, S::NumericallyUnresolved);
        q = {};
        const auto schedule = algo::planLayerSchedule(q);
        require(schedule.planned, "schedule baseline");
        q.maximumLayers = static_cast<std::uint32_t>(schedule.offsets.size() - 1);
        q.maximumBytes = schedule.offsets.size() * sizeof(double);
        require(algo::planLayerSchedule(q).planned, "schedule exact budgets");
        ++scheduleChecks;
        --q.maximumLayers;
        rejectSchedule(q, S::LayerBudget);
        ++q.maximumLayers;
        --q.maximumBytes;
        rejectSchedule(q, S::ByteBudget);
        // Explicit fractions are also a standalone addition to the existing
        // mixed-sweep API; malformed schedules must preserve its destination.
        for (auto fractions :
             std::vector<std::vector<double>>{{0, .5},
                                              {0, .5, .4},
                                              {0, 0, 1},
                                              {0, .01, 1},
                                              {0, std::numeric_limits<double>::quiet_NaN(), 1}}) {
            algo::PreparedMixedSweepRequest p;
            p.layers = 2;
            p.layerFractions = fractions;
            mesh::Mesh old;
            old.addNode({9, 0, 0}, {0, 9});
            const auto h = mesh::meshHash(old);
            const auto r = algo::generatePreparedMixedSweep(small, p, old);
            require(!r.published() && r.reason == algo::PreparedSweepReason::InvalidGrowth &&
                        mesh::meshHash(old) == h,
                    "explicit schedule rollback");
            ++refusals;
        }
        std::cout << "{\"scheduleChecks\":" << scheduleChecks << "}\n";
        std::cout << "{\"passed\":" << cases << ",\"refusals\":" << refusals << "}\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
