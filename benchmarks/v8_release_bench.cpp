// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MeshingPlan.hpp"
#include "paacs/algo/PreparedSweepMesher.hpp"
#include "paacs/algo/TransitionTemplateCatalog.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sizing/CompiledSpatialField.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;

double median(std::vector<double> values) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2U;
    return values.size() & 1U ? values[middle]
                              : 0.5 * (values[middle - 1U] + values[middle]);
}

template <class Function>
double measure(Function&& function) {
    const auto begin = Clock::now();
    function();
    return std::chrono::duration<double>(Clock::now() - begin).count();
}

std::uint32_t argument(char** argv, int argc, int index,
                       std::uint32_t fallback) {
    if (index >= argc) return fallback;
    const unsigned long value = std::stoul(argv[index]);
    if (value == 0 || value > 1'000'000UL)
        throw std::runtime_error("numeric argument outside supported range");
    return static_cast<std::uint32_t>(value);
}

paacs::mesh::Mesh mixedSurface(std::uint32_t n) {
    paacs::mesh::Mesh mesh;
    const std::uint64_t nodes =
        static_cast<std::uint64_t>(n + 1U) * (n + 1U);
    const std::uint64_t leftColumns = n / 2U;
    const std::uint64_t quads = leftColumns * n;
    const std::uint64_t triangles =
        2U * static_cast<std::uint64_t>(n - leftColumns) * n;
    mesh.reserveAppend(static_cast<std::size_t>(nodes),
                       static_cast<std::size_t>(quads + triangles),
                       static_cast<std::size_t>(4U * quads + 3U * triangles));
    const auto node = [n](std::uint32_t i, std::uint32_t j) {
        return static_cast<paacs::Index>(
            static_cast<std::uint64_t>(j) * (n + 1U) + i);
    };
    for (std::uint32_t j = 0; j <= n; ++j)
        for (std::uint32_t i = 0; i <= n; ++i)
            (void)mesh.addNode({static_cast<double>(i) / n,
                                static_cast<double>(j) / n, 0.0},
                               {2, 100});
    for (std::uint32_t j = 0; j < n; ++j) {
        for (std::uint32_t i = 0; i < n; ++i) {
            const auto a = node(i, j);
            const auto b = node(i + 1U, j);
            const auto c = node(i + 1U, j + 1U);
            const auto d = node(i, j + 1U);
            if (i < leftColumns) {
                (void)mesh.addCell(paacs::mesh::CellType::Quad4,
                                   {a, b, c, d}, {2, 100});
            } else {
                (void)mesh.addCell(paacs::mesh::CellType::Tri3,
                                   {a, b, c}, {2, 100});
                (void)mesh.addCell(paacs::mesh::CellType::Tri3,
                                   {a, c, d}, {2, 100});
            }
        }
    }
    return mesh;
}

paacs::mesh::Mesh quadSurface(std::uint32_t n) {
    paacs::mesh::Mesh mesh;
    const std::uint64_t nodes =
        static_cast<std::uint64_t>(n + 1U) * (n + 1U);
    const std::uint64_t cells = static_cast<std::uint64_t>(n) * n;
    mesh.reserveAppend(static_cast<std::size_t>(nodes),
                       static_cast<std::size_t>(cells),
                       static_cast<std::size_t>(4U * cells));
    const auto node = [n](std::uint32_t i, std::uint32_t j) {
        return static_cast<paacs::Index>(
            static_cast<std::uint64_t>(j) * (n + 1U) + i);
    };
    for (std::uint32_t j = 0; j <= n; ++j)
        for (std::uint32_t i = 0; i <= n; ++i)
            (void)mesh.addNode({static_cast<double>(i) / n,
                                static_cast<double>(j) / n, 0.0},
                               {2, 200});
    for (std::uint32_t j = 0; j < n; ++j)
        for (std::uint32_t i = 0; i < n; ++i)
            (void)mesh.addCell(paacs::mesh::CellType::Quad4,
                               {node(i, j), node(i + 1U, j),
                                node(i + 1U, j + 1U), node(i, j + 1U)},
                               {2, 200});
    return mesh;
}

int runSweep(int argc, char** argv) {
    const std::uint32_t n = argument(argv, argc, 2, 256);
    const std::uint32_t layers = argument(argv, argc, 3, 8);
    const std::uint32_t runs = argument(argv, argc, 4, 9);
    const std::uint32_t threads = argument(argv, argc, 5, 4);
    const bool fullAudit = argc > 6 && std::string(argv[6]) == "full";
    const paacs::mesh::Mesh source = mixedSurface(n);
    paacs::algo::PreparedMixedSweepRequest request;
    request.layers = layers;
    request.direction = {0.0, 0.0, 1.0};
    request.totalThickness = 1.0;
    request.cellToCellGrowth = 1.15;
    request.maximumAllowedGrowth = 1.2;
    request.executionThreads = static_cast<int>(threads);
    request.fullAudit = fullAudit;

    paacs::mesh::Mesh warm;
    if (!paacs::algo::generatePreparedMixedSweep(source, request, warm)
             .published())
        return 2;
    std::vector<double> times;
    paacs::mesh::Mesh retained;
    paacs::algo::PreparedSweepOutcome last;
    for (std::uint32_t run = 0; run < runs; ++run) {
        paacs::mesh::Mesh current;
        times.push_back(measure([&] {
            last = paacs::algo::generatePreparedMixedSweep(
                source, request, current);
        }));
        if (!last.published()) return 3;
        if (run + 1U == runs) retained.swap(current);
    }
    const double seconds = median(times);
    const std::uint64_t volumeCells = last.metrics.wedges + last.metrics.hexes;
    std::cout << "mode=mixed-sweep\n"
              << "audit=" << (fullAudit ? "full" : "certified") << '\n'
              << "surface_n=" << n << '\n'
              << "layers=" << layers << '\n'
              << "source_triangles=" << last.metrics.sourceTriangles << '\n'
              << "source_quads=" << last.metrics.sourceQuads << '\n'
              << "volume_cells=" << volumeCells << '\n'
              << "median_s=" << std::setprecision(12) << seconds << '\n'
              << "million_cells_per_s="
              << static_cast<double>(volumeCells) / seconds / 1.0e6 << '\n'
              << "minimum_cell_volume="
              << last.metrics.minimumCellVolume << '\n'
              << "total_volume=" << last.metrics.totalVolume << '\n'
              << "mesh_hash=" << paacs::mesh::meshHash(retained) << '\n';
    return 0;
}

int runSpatial(int argc, char** argv) {
    const std::uint32_t n = argument(argv, argc, 2, 96);
    const std::uint32_t runs = argument(argv, argc, 3, 9);
    const std::uint32_t threads = argument(argv, argc, 4, 4);
    paacs::sizing::CompiledSpatialFieldRequest request;
    request.lower = {-1.0, -1.0, -1.0};
    request.upper = {1.0, 1.0, 1.0};
    request.samplesX = request.samplesY = request.samplesZ = n;
    request.baseSize = 0.4;
    request.minimumSize = 0.01;
    request.maximumSize = 0.4;
    request.maximumAdjacentGrowth = 1.2;
    request.maximumGrowthSweeps = 32;
    request.executionThreads = static_cast<int>(threads);
    request.controls.push_back({paacs::sizing::SpatialControlKind::Ball,
                                {0.0, 0.0, 0.0}, {}, 0.35, 0.0, 0.025, 0.0});
    request.controls.push_back({paacs::sizing::SpatialControlKind::Cylinder,
                                {-0.8, 0.0, -0.8}, {0.8, 0.0, 0.8},
                                0.18, 0.0, 0.04, 0.0});
    auto warm = paacs::sizing::compileSpatialField(request);
    if (!warm.compiled) return 2;
    std::vector<double> times;
    paacs::sizing::CompiledSpatialFieldOutcome last;
    for (std::uint32_t run = 0; run < runs; ++run) {
        times.push_back(measure([&] {
            last = paacs::sizing::compileSpatialField(request);
        }));
        if (!last.compiled) return 3;
    }
    constexpr std::uint64_t queryCount = 10'000'000ULL;
    volatile double sink = 0.0;
    const double querySeconds = measure([&] {
        double value = 0.0;
        for (std::uint64_t q = 0; q < queryCount; ++q) {
            const double x = static_cast<double>((q * 48271U) % 1'000'003U) /
                             500'001.5 - 1.0;
            const double y = static_cast<double>((q * 69621U) % 1'000'033U) /
                             500'016.5 - 1.0;
            const double z = static_cast<double>((q * 40699U) % 1'000'037U) /
                             500'018.5 - 1.0;
            value += last.field.eval({x, y, z});
        }
        sink = value;
    });
    const double seconds = median(times);
    std::cout << "mode=spatial-field\n"
              << "samples=" << last.metrics.samples << '\n'
              << "compile_median_s=" << std::setprecision(12) << seconds << '\n'
              << "million_samples_per_s="
              << static_cast<double>(last.metrics.samples) / seconds / 1.0e6
              << '\n'
              << "million_queries_per_s="
              << static_cast<double>(queryCount) / querySeconds / 1.0e6 << '\n'
              << "maximum_growth="
              << last.metrics.measuredMaximumAdjacentGrowth << '\n'
              << "field_hash=" << last.field.deterministicHash() << '\n'
              << "sink=" << sink << '\n';
    return 0;
}

int runPlan(int argc, char** argv) {
    const std::uint32_t bodies = argument(argv, argc, 2, 100000);
    const std::uint32_t runs = argument(argv, argc, 3, 9);
    paacs::algo::MeshingPlanRequest request;
    request.bodies.reserve(bodies);
    request.portals.reserve(bodies > 0 ? bodies - 1U : 0U);
    for (std::uint32_t i = 0; i < bodies; ++i) {
        const bool mapped = (i & 1U) == 0U;
        request.bodies.push_back({
            static_cast<int>(i + 1U),
            paacs::sizing::VolumeEngineKind::Auto,
            mapped
                ? paacs::algo::BodyMeshingCapability::GeneralTet4 |
                      paacs::algo::BodyMeshingCapability::MappedBlock |
                      paacs::algo::BodyMeshingCapability::PyramidTransition
                : paacs::algo::BodyMeshingCapability::GeneralTet4 |
                      paacs::algo::BodyMeshingCapability::PyramidTransition,
            mapped ? static_cast<std::uint64_t>(0x1000U + i) : 0U,
            64, 32, true});
        if (i != 0U) {
            request.portals.push_back({
                static_cast<int>(1'000'000U + i),
                static_cast<int>(i), static_cast<int>(i + 1U),
                mapped ? paacs::algo::PortalTopology::Quad4
                       : paacs::algo::PortalTopology::Tri3,
                16, 8, static_cast<std::uint64_t>(0xABC000U + i), true});
        }
    }
    auto warm = paacs::algo::compileMeshingPlan(request);
    if (!warm.compiled) return 2;
    std::vector<double> times;
    paacs::algo::MeshingPlanOutcome last;
    for (std::uint32_t run = 0; run < runs; ++run) {
        times.push_back(measure([&] {
            last = paacs::algo::compileMeshingPlan(request);
        }));
        if (!last.compiled) return 3;
    }
    const double seconds = median(times);
    const std::uint64_t rows =
        static_cast<std::uint64_t>(last.metrics.bodies) + last.metrics.portals;
    std::cout << "mode=meshing-plan\n"
              << "bodies=" << last.metrics.bodies << '\n'
              << "portals=" << last.metrics.portals << '\n'
              << "median_s=" << std::setprecision(12) << seconds << '\n'
              << "million_rows_per_s="
              << static_cast<double>(rows) / seconds / 1.0e6 << '\n'
              << "plan_hash=" << last.plan.deterministicHash() << '\n';
    return 0;
}

int runPyramid(int argc, char** argv) {
    const std::uint32_t n = argument(argv, argc, 2, 512);
    const std::uint32_t runs = argument(argv, argc, 3, 9);
    const paacs::mesh::Mesh source = quadSurface(n);
    paacs::algo::PyramidTransitionRequest request;
    request.apex = {0.5, 0.5, 1.0};
    request.emitBoundaryFaces = false;
    paacs::mesh::Mesh warm;
    if (!paacs::algo::generatePyramidTransitionFan(source, request, warm)
             .published())
        return 2;
    std::vector<double> times;
    paacs::mesh::Mesh retained;
    paacs::algo::TransitionOutcome last;
    for (std::uint32_t run = 0; run < runs; ++run) {
        paacs::mesh::Mesh current;
        times.push_back(measure([&] {
            last = paacs::algo::generatePyramidTransitionFan(
                source, request, current);
        }));
        if (!last.published()) return 3;
        if (run + 1U == runs) retained.swap(current);
    }
    const double seconds = median(times);
    std::cout << "mode=pyramid-transition\n"
              << "cells=" << last.metrics.outputPyramid5 << '\n'
              << "median_s=" << std::setprecision(12) << seconds << '\n'
              << "million_cells_per_s="
              << static_cast<double>(last.metrics.outputPyramid5) /
                     seconds / 1.0e6 << '\n'
              << "total_volume=" << last.metrics.totalVolume << '\n'
              << "mesh_hash=" << paacs::mesh::meshHash(retained) << '\n';
    return 0;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cerr << "usage: paacs-v8-bench MODE [arguments]\n";
            return 1;
        }
        const std::string mode = argv[1];
        if (mode == "mixed-sweep") return runSweep(argc, argv);
        if (mode == "spatial") return runSpatial(argc, argv);
        if (mode == "plan") return runPlan(argc, argv);
        if (mode == "pyramid") return runPyramid(argc, argv);
        throw std::runtime_error("unknown benchmark mode");
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
