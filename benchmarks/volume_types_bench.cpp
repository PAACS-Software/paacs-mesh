// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/StructuredPolyhedralMesher.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/mesh/PolyhedralMesh.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Options {
    std::string mode;
    std::uint32_t n = 0;
    int warmups = 1;
    int runs = 7;
};

Options parse(int argc, char** argv) {
    if (argc < 3 || argc > 5) {
        throw std::runtime_error(
            "usage: paacs_v5_volume_bench MODE N [WARMUPS] [RUNS]");
    }
    Options out;
    out.mode = argv[1];
    const unsigned long parsedN = std::stoul(argv[2]);
    if (parsedN == 0 || parsedN > 4096)
        throw std::runtime_error("N must be in [1,4096]");
    out.n = static_cast<std::uint32_t>(parsedN);
    if (argc >= 4) out.warmups = std::stoi(argv[3]);
    if (argc >= 5) out.runs = std::stoi(argv[4]);
    if (out.warmups < 0 || out.runs <= 0)
        throw std::runtime_error("invalid warmup/run count");
    return out;
}

paacs::algo::StructuredBlockGeometry geometry(std::uint32_t n) {
    paacs::algo::StructuredBlockGeometry out;
    out.cellsU = n;
    out.cellsV = n;
    out.cellsW = n;
    out.axisU = {2.0, 0.0, 0.0};
    out.axisV = {0.2, 3.0, 0.0};
    out.axisW = {0.1, 0.15, 4.0};
    out.volumeClassification = {3, 100};
    for (std::size_t i = 0; i < out.boundaryClassifications.size(); ++i)
        out.boundaryClassifications[i] =
            {2, static_cast<int>(201 + i)};
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
    const std::size_t middle = values.size() / 2;
    if ((values.size() & 1U) != 0U) return values[middle];
    return 0.5 * (values[middle - 1] + values[middle]);
}

void printTimes(const std::vector<double>& values) {
    std::cout << "times_s=";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << std::setprecision(9) << values[i];
    }
    std::cout << '\n';
}

int runHex(const Options& options, bool fullAudit) {
    paacs::algo::MappedHexBlockRequest request;
    request.geometry = geometry(options.n);
    request.validation = fullAudit
        ? paacs::algo::MappedHexValidationMode::FullTopologyAudit
        : paacs::algo::MappedHexValidationMode::CertifiedAffine;

    for (int i = 0; i < options.warmups; ++i) {
        paacs::mesh::Mesh warm;
        const auto out = paacs::algo::generateMappedHexBlock(request, warm);
        if (!out.published()) return 2;
    }

    std::vector<double> times;
    times.reserve(static_cast<std::size_t>(options.runs));
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
    const double rate = static_cast<double>(last.metrics.hexes) / med / 1.0e6;
    std::cout << "mode=hex-" << (fullAudit ? "full" : "certified") << '\n'
              << "n=" << options.n << '\n'
              << "nodes=" << retained.nodeCount() << '\n'
              << "hexes=" << last.metrics.hexes << '\n'
              << "boundary_quads=" << last.metrics.boundaryQuads << '\n'
              << "median_s=" << std::setprecision(9) << med << '\n'
              << "million_cells_per_s=" << std::setprecision(9) << rate << '\n'
              << "minimum_scaled_jacobian="
              << std::setprecision(17)
              << last.metrics.audit.minimumScaledJacobian << '\n'
              << "total_volume=" << std::setprecision(17)
              << last.metrics.audit.totalVolume << '\n'
              << "mesh_hash=" << paacs::mesh::meshHash(retained) << '\n';
    printTimes(times);
    return 0;
}

int runPoly(const Options& options, bool fullAudit) {
    const auto g = geometry(options.n);
    const auto validation = fullAudit
        ? paacs::algo::StructuredPolyValidationMode::FullAudit
        : paacs::algo::StructuredPolyValidationMode::CertifiedStructured;

    for (int i = 0; i < options.warmups; ++i) {
        paacs::mesh::PolyhedralMesh warm;
        const auto out = paacs::algo::generateStructuredPolyhedralBlock(
            g, warm, validation);
        if (!out.published()) return 2;
    }

    std::vector<double> times;
    times.reserve(static_cast<std::size_t>(options.runs));
    paacs::algo::StructuredPolyOutcome last;
    paacs::mesh::PolyhedralMesh retained;
    for (int i = 0; i < options.runs; ++i) {
        paacs::mesh::PolyhedralMesh current;
        times.push_back(measure([&] {
            last = paacs::algo::generateStructuredPolyhedralBlock(
                g, current, validation);
        }));
        if (!last.published()) return 3;
        if (i + 1 == options.runs) retained.swap(current);
    }

    const double med = median(times);
    const double rate = static_cast<double>(last.metrics.cells) / med / 1.0e6;
    std::cout << "mode=poly-" << (fullAudit ? "full" : "certified") << '\n'
              << "n=" << options.n << '\n'
              << "nodes=" << retained.nodeCount() << '\n'
              << "faces=" << retained.faceCount() << '\n'
              << "cells=" << retained.cellCount() << '\n'
              << "median_s=" << std::setprecision(9) << med << '\n'
              << "million_cells_per_s=" << std::setprecision(9) << rate << '\n'
              << "minimum_cell_volume=" << std::setprecision(17)
              << last.metrics.audit.minimumCellVolume << '\n'
              << "total_volume=" << std::setprecision(17)
              << last.metrics.audit.totalVolume << '\n'
              << "poly_hash="
              << paacs::mesh::polyhedralStrictHash(retained) << '\n';
    printTimes(times);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse(argc, argv);
        if (options.mode == "hex-certified")
            return runHex(options, false);
        if (options.mode == "hex-full")
            return runHex(options, true);
        if (options.mode == "poly-certified")
            return runPoly(options, false);
        if (options.mode == "poly-full")
            return runPoly(options, true);
        throw std::runtime_error("unknown mode");
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
