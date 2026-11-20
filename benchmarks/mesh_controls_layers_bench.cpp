// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/StraightPrismLayerMesher.hpp"
#include "paacs/algo/StructuredPolyhedralMesher.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/mesh/PolyhedralMesh.hpp"
#include "paacs/sizing/MeshControls.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;

struct Options {
    std::string mode;
    std::uint32_t n = 0;
    std::uint32_t layers = 8;
    int warmups = 1;
    int runs = 7;
};

Options parse(int argc, char** argv) {
    if (argc < 3 || argc > 6)
        throw std::runtime_error(
            "usage: paacs_v6_bench MODE N [LAYERS] [WARMUPS] [RUNS]");
    Options out;
    out.mode = argv[1];
    const auto parsedN = std::stoul(argv[2]);
    const unsigned long maximumN = out.mode == "controls" ? 1'000'000UL : 4096UL;
    if (parsedN == 0 || parsedN > maximumN)
        throw std::runtime_error("N outside supported benchmark range");
    out.n = static_cast<std::uint32_t>(parsedN);
    if (argc >= 4) out.layers = static_cast<std::uint32_t>(std::stoul(argv[3]));
    if (argc >= 5) out.warmups = std::stoi(argv[4]);
    if (argc >= 6) out.runs = std::stoi(argv[5]);
    if (out.layers == 0 || out.warmups < 0 || out.runs <= 0)
        throw std::runtime_error("invalid layer/warmup/run count");
    return out;
}

template <class Fn>
double measure(Fn&& fn) {
    const auto start = Clock::now();
    fn();
    return std::chrono::duration<double>(Clock::now() - start).count();
}

double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    const std::size_t m = values.size() / 2;
    return (values.size() & 1U) != 0U
        ? values[m] : 0.5 * (values[m - 1] + values[m]);
}

void printTimes(const std::vector<double>& values) {
    std::cout << "times_s=";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << std::setprecision(9) << values[i];
    }
    std::cout << '\n';
}

paacs::algo::StructuredBlockGeometry block(std::uint32_t n) {
    paacs::algo::StructuredBlockGeometry out;
    out.cellsU = n;
    out.cellsV = n;
    out.cellsW = n;
    out.axisU = {2.0, 0.0, 0.0};
    out.axisV = {0.2, 3.0, 0.0};
    out.axisW = {0.1, 0.15, 4.0};
    out.volumeClassification = {3, 100};
    for (std::size_t i = 0; i < out.boundaryClassifications.size(); ++i)
        out.boundaryClassifications[i] = {2, static_cast<int>(201 + i)};
    return out;
}

paacs::algo::StructuredBlockGrowth gradedGrowth() {
    paacs::algo::StructuredBlockGrowth out;
    out.u.mode = paacs::algo::StructuredAxisGrowthMode::
        EndToStartCellWidthRatio;
    out.u.endToStartCellWidthRatio = 8.0;
    out.v.mode = paacs::algo::StructuredAxisGrowthMode::
        EndToStartCellWidthRatio;
    out.v.endToStartCellWidthRatio = 0.25;
    out.w.mode = paacs::algo::StructuredAxisGrowthMode::
        EndToStartCellWidthRatio;
    out.w.endToStartCellWidthRatio = 4.0;
    return out;
}

paacs::mesh::Mesh squareSurface(std::uint32_t n) {
    paacs::mesh::Mesh out;
    const std::uint64_t side = static_cast<std::uint64_t>(n) + 1U;
    const std::uint64_t triangles = 2ULL * n * n;
    out.reserveAppend(static_cast<std::size_t>(side * side),
                      static_cast<std::size_t>(triangles),
                      static_cast<std::size_t>(3ULL * triangles));
    for (std::uint32_t j = 0; j <= n; ++j)
        for (std::uint32_t i = 0; i <= n; ++i) {
            const bool boundary = i == 0 || j == 0 || i == n || j == n;
            (void)out.addNode(
                {static_cast<double>(i) / n,
                 static_cast<double>(j) / n, 0.0},
                boundary ? paacs::mesh::Classification{1, 10}
                         : paacs::mesh::Classification{2, 20});
        }
    const auto node = [n](std::uint32_t i, std::uint32_t j) {
        return static_cast<paacs::Index>(
            static_cast<std::uint64_t>(j) * (n + 1U) + i);
    };
    for (std::uint32_t j = 0; j < n; ++j)
        for (std::uint32_t i = 0; i < n; ++i) {
            const auto a = node(i, j);
            const auto b = node(i + 1U, j);
            const auto c = node(i + 1U, j + 1U);
            const auto d = node(i, j + 1U);
            (void)out.addCell(paacs::mesh::CellType::Tri3,
                              {a, b, c}, {2, 20});
            (void)out.addCell(paacs::mesh::CellType::Tri3,
                              {a, c, d}, {2, 20});
        }
    return out;
}

int runHex(const Options& options, bool graded) {
    paacs::algo::MappedHexBlockRequest request;
    request.geometry = block(options.n);
    if (graded) request.growth = gradedGrowth();
    request.validation = paacs::algo::MappedHexValidationMode::CertifiedAffine;

    for (int i = 0; i < options.warmups; ++i) {
        paacs::mesh::Mesh warm;
        if (!paacs::algo::generateMappedHexBlock(request, warm).published())
            return 2;
    }
    std::vector<double> times;
    paacs::algo::MappedHexOutcome last;
    paacs::mesh::Mesh retained;
    for (int i = 0; i < options.runs; ++i) {
        paacs::mesh::Mesh current;
        times.push_back(measure([&] {
            last = paacs::algo::generateMappedHexBlock(request, current);
        }));
        if (!last.published()) return 3;
        if (i + 1 == options.runs) retained.swap(current);
    }
    const double med = median(times);
    std::cout << "mode=hex-" << (graded ? "graded" : "uniform") << '\n'
              << "n=" << options.n << '\n'
              << "cells=" << last.metrics.hexes << '\n'
              << "median_s=" << std::setprecision(9) << med << '\n'
              << "million_cells_per_s=" << std::setprecision(9)
              << static_cast<double>(last.metrics.hexes) / med / 1.0e6 << '\n'
              << "growth_work=" << last.metrics.growthPreparationWorkUnits << '\n'
              << "max_adjacent_growth=" << std::setprecision(17)
              << last.metrics.maximumAdjacentGrowth << '\n'
              << "minimum_scaled_jacobian="
              << last.metrics.audit.minimumScaledJacobian << '\n'
              << "mesh_hash=" << paacs::mesh::meshHash(retained) << '\n';
    printTimes(times);
    return 0;
}

int runPoly(const Options& options, bool graded) {
    paacs::algo::StructuredPolyBlockRequest request;
    request.geometry = block(options.n);
    if (graded) request.growth = gradedGrowth();
    request.validation = paacs::algo::StructuredPolyValidationMode::
        CertifiedStructured;
    for (int i = 0; i < options.warmups; ++i) {
        paacs::mesh::PolyhedralMesh warm;
        if (!paacs::algo::generateStructuredPolyhedralBlock(request, warm)
                 .published()) return 2;
    }
    std::vector<double> times;
    paacs::algo::StructuredPolyOutcome last;
    paacs::mesh::PolyhedralMesh retained;
    for (int i = 0; i < options.runs; ++i) {
        paacs::mesh::PolyhedralMesh current;
        times.push_back(measure([&] {
            last = paacs::algo::generateStructuredPolyhedralBlock(
                request, current);
        }));
        if (!last.published()) return 3;
        if (i + 1 == options.runs) retained.swap(current);
    }
    const double med = median(times);
    std::cout << "mode=poly-" << (graded ? "graded" : "uniform") << '\n'
              << "n=" << options.n << '\n'
              << "cells=" << last.metrics.cells << '\n'
              << "faces=" << retained.faceCount() << '\n'
              << "median_s=" << std::setprecision(9) << med << '\n'
              << "million_cells_per_s=" << std::setprecision(9)
              << static_cast<double>(last.metrics.cells) / med / 1.0e6 << '\n'
              << "growth_work=" << last.metrics.growthPreparationWorkUnits << '\n'
              << "max_adjacent_growth=" << std::setprecision(17)
              << last.metrics.maximumAdjacentGrowth << '\n'
              << "poly_hash="
              << paacs::mesh::polyhedralStrictHash(retained) << '\n';
    printTimes(times);
    return 0;
}

int runPrism(const Options& options) {
    const paacs::mesh::Mesh source = squareSurface(options.n);
    paacs::algo::StraightPrismLayerRequest request;
    request.layers = options.layers;
    request.firstLayerThickness = 1.0e-4;
    request.cellToCellGrowth = 1.2;
    request.maximumAllowedGrowth = 1.3;
    request.maximumTotalThickness = 1.0;
    request.validation = paacs::algo::PrismLayerValidationMode::
        CertifiedStraightSweep;
    for (int i = 0; i < options.warmups; ++i) {
        paacs::mesh::Mesh warm;
        if (!paacs::algo::generateStraightPrismLayers(source, request, warm)
                 .published()) return 2;
    }
    std::vector<double> times;
    paacs::algo::PrismLayerOutcome last;
    paacs::mesh::Mesh retained;
    for (int i = 0; i < options.runs; ++i) {
        paacs::mesh::Mesh current;
        times.push_back(measure([&] {
            last = paacs::algo::generateStraightPrismLayers(
                source, request, current);
        }));
        if (!last.published()) return 3;
        if (i + 1 == options.runs) retained.swap(current);
    }
    const double med = median(times);
    std::cout << "mode=prism-straight" << '\n'
              << "surface_n=" << options.n << '\n'
              << "layers=" << options.layers << '\n'
              << "source_triangles=" << last.metrics.sourceTriangles << '\n'
              << "prisms=" << last.metrics.prisms << '\n'
              << "median_s=" << std::setprecision(9) << med << '\n'
              << "million_cells_per_s=" << std::setprecision(9)
              << static_cast<double>(last.metrics.prisms) / med / 1.0e6 << '\n'
              << "total_thickness=" << std::setprecision(17)
              << last.metrics.totalThickness << '\n'
              << "mesh_hash=" << paacs::mesh::meshHash(retained) << '\n';
    printTimes(times);
    return 0;
}

int runControls(const Options& options) {
    const std::uint32_t count = options.n;
    paacs::sizing::MeshControlSpecification spec;
    spec.minimumSize = 1.0e-4;
    spec.maximumSize = 10.0;
    spec.entitySizes.reserve(count);
    spec.bodyMethods.reserve(count);
    for (std::uint32_t i = 0; i < count; ++i) {
        spec.entitySizes.push_back({3, static_cast<int>(i),
                                    0.01 + 1.0e-7 * i});
        spec.bodyMethods.push_back({static_cast<int>(i),
            (i & 1U) ? paacs::sizing::VolumeEngineKind::Tet4
                     : paacs::sizing::VolumeEngineKind::MappedHex8});
    }
    std::vector<double> compileTimes;
    paacs::sizing::MeshControlCompileOutcome compiled;
    for (int i = 0; i < options.warmups; ++i) {
        const auto warm = paacs::sizing::compileMeshControls(spec);
        if (!warm.compiled) return 2;
    }
    for (int i = 0; i < options.runs; ++i)
        compileTimes.push_back(measure([&] {
            compiled = paacs::sizing::compileMeshControls(spec);
        }));
    if (!compiled.compiled) return 3;

    constexpr std::uint64_t kQueries = 5'000'000ULL;
    volatile double sink = 0.0;
    const double querySeconds = measure([&] {
        double local = 0.0;
        for (std::uint64_t i = 0; i < kQueries; ++i) {
            const int tag = static_cast<int>((i * 48271ULL) % count);
            local += compiled.controls.effectiveSize(3, tag, 1.0);
            local += static_cast<unsigned>(compiled.controls.engineForBody(tag));
        }
        sink = local;
    });
    std::cout << "mode=controls" << '\n'
              << "rows=" << count << '\n'
              << "compile_median_s=" << std::setprecision(9)
              << median(compileTimes) << '\n'
              << "queries=" << kQueries * 2ULL << '\n'
              << "query_s=" << querySeconds << '\n'
              << "million_queries_per_s="
              << (static_cast<double>(kQueries) * 2.0 / querySeconds / 1.0e6)
              << '\n'
              << "sink=" << sink << '\n';
    printTimes(compileTimes);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse(argc, argv);
        if (options.mode == "hex-uniform") return runHex(options, false);
        if (options.mode == "hex-graded") return runHex(options, true);
        if (options.mode == "poly-uniform") return runPoly(options, false);
        if (options.mode == "poly-graded") return runPoly(options, true);
        if (options.mode == "prism") return runPrism(options);
        if (options.mode == "controls") return runControls(options);
        throw std::runtime_error("unknown mode");
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
