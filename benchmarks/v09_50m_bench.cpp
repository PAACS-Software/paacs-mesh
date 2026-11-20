// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sdk/Version.hpp"

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

struct Arguments {
    std::uint32_t cells = 128;
    int threads = 0;
    std::uint32_t runs = 11;
    double thresholdMillionCells = 50.0;
    bool cold = false;
    bool emitBoundary = true;
};

std::uint32_t parseUnsigned(const char* text, const char* name) {
    const unsigned long value = std::stoul(text);
    if (value == 0 || value > 4096UL)
        throw std::invalid_argument(std::string(name) +
                                    " outside [1,4096]");
    return static_cast<std::uint32_t>(value);
}

Arguments parseArguments(int argc, char** argv) {
    Arguments out;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--cells" && i + 1 < argc) {
            out.cells = parseUnsigned(argv[++i], "cells");
        } else if (option == "--threads" && i + 1 < argc) {
            out.threads = std::stoi(argv[++i]);
            if (out.threads < 0)
                throw std::invalid_argument("threads must be nonnegative");
        } else if (option == "--runs" && i + 1 < argc) {
            out.runs = parseUnsigned(argv[++i], "runs");
            if (out.runs < 3)
                throw std::invalid_argument("runs must be at least three");
        } else if (option == "--threshold" && i + 1 < argc) {
            out.thresholdMillionCells = std::stod(argv[++i]);
            if (!std::isfinite(out.thresholdMillionCells) ||
                out.thresholdMillionCells < 0.0) {
                throw std::invalid_argument(
                    "threshold must be finite and nonnegative");
            }
        } else if (option == "--cold") {
            out.cold = true;
        } else if (option == "--no-boundary") {
            out.emitBoundary = false;
        } else if (option == "--help") {
            std::cout
                << "paacs-v09-50m-bench [options]\n"
                << "  --cells N       N x N x N volume cells (default 128)\n"
                << "  --threads N     request-local threads; 0=auto\n"
                << "  --runs N        measured runs after warmup (default 11)\n"
                << "  --threshold M   pass gate in million Hex8/s (default 50)\n"
                << "  --cold          benchmark fresh one-shot allocation\n"
                << "  --no-boundary   omit explicit exterior Quad4 cells\n";
            std::exit(EXIT_SUCCESS);
        } else {
            throw std::invalid_argument("unknown option: " + option);
        }
    }
    return out;
}

paacs::algo::MappedHexBlockRequest makeRequest(const Arguments& args) {
    paacs::algo::MappedHexBlockRequest request;
    request.geometry.origin = {1.0, -2.0, 0.25};
    request.geometry.axisU = {2.0, 0.0, 0.0};
    request.geometry.axisV = {0.2, 3.0, 0.0};
    request.geometry.axisW = {0.1, 0.15, 4.0};
    request.geometry.cellsU = args.cells;
    request.geometry.cellsV = args.cells;
    request.geometry.cellsW = args.cells;
    request.geometry.volumeClassification = {3, 900};
    for (std::size_t i = 0;
         i < request.geometry.boundaryClassifications.size(); ++i) {
        request.geometry.boundaryClassifications[i] = {
            2, 910 + static_cast<int>(i)};
    }
    request.executionThreads = args.threads;
    request.emitBoundaryQuads = args.emitBoundary;
    request.validation =
        paacs::algo::MappedHexValidationMode::CertifiedAffine;
    return request;
}

struct Distribution {
    double minimum = 0.0;
    double median = 0.0;
    double p95 = 0.0;
    double maximum = 0.0;
};

Distribution summarize(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    Distribution out;
    out.minimum = values.front();
    out.median = values[values.size() / 2U];
    const auto p95Index = static_cast<std::size_t>(
        std::floor(0.95 * static_cast<double>(values.size() - 1U)));
    out.p95 = values[p95Index];
    out.maximum = values.back();
    return out;
}

template <class Function>
double measure(Function&& function) {
    const auto begin = Clock::now();
    function();
    return std::chrono::duration<double>(Clock::now() - begin).count();
}

int run(const Arguments& args) {
    const auto request = makeRequest(args);
    paacs::algo::MappedHexWorkspace workspace;
    paacs::mesh::Mesh retained;
    paacs::algo::MappedHexOutcome last;

    // Full-size warmup is mandatory. The reusable mode needs two ownership
    // exchanges before the workspace and destination form a stable double
    // buffer; the third call proves the fast path is active.
    const int warmupCalls = args.cold ? 1 : 3;
    for (int i = 0; i < warmupCalls; ++i) {
        if (args.cold) {
            paacs::mesh::Mesh warm;
            last = paacs::algo::generateMappedHexBlock(request, warm);
        } else {
            last = paacs::algo::generateMappedHexBlock(
                request, workspace, retained);
        }
        if (!last.published()) {
            std::cerr << "warmup generation failed\n";
            return 2;
        }
    }
    if (!args.cold && !last.metrics.reusedConstructionStorage) {
        std::cerr << "workspace did not enter the reusable fast path\n";
        return 3;
    }

    std::vector<double> seconds;
    seconds.reserve(args.runs);
    for (std::uint32_t run = 0; run < args.runs; ++run) {
        if (args.cold) {
            paacs::mesh::Mesh current;
            const double elapsed = measure([&] {
                last = paacs::algo::generateMappedHexBlock(
                    request, current);
            });
            if (!last.published()) return 4;
            retained.swap(current);
            seconds.push_back(elapsed);
        } else {
            const double elapsed = measure([&] {
                last = paacs::algo::generateMappedHexBlock(
                    request, workspace, retained);
            });
            if (!last.published() ||
                !last.metrics.reusedConstructionStorage) {
                return 5;
            }
            seconds.push_back(elapsed);
        }
    }

    const Distribution timing = summarize(seconds);
    const double rate = static_cast<double>(last.metrics.hexes) /
        timing.median / 1.0e6;
    const double logicalBandwidth =
        static_cast<double>(last.metrics.logicalBytes) /
        timing.median / 1.0e9;
    const bool passed = rate >= args.thresholdMillionCells;
    const std::uint64_t hash = paacs::mesh::meshHash(retained);

    std::cout << std::setprecision(12)
              << "{\"api\":\"" << paacs::sdk::kApiVersion
              << "\",\"mode\":\"" << (args.cold ? "cold" : "reused")
              << "\",\"cells_per_axis\":" << args.cells
              << ",\"threads\":" << args.threads
              << ",\"runs\":" << args.runs
              << ",\"hex8\":" << last.metrics.hexes
              << ",\"boundary_quad4\":" << last.metrics.boundaryQuads
              << ",\"nodes\":" << last.metrics.nodes
              << ",\"median_seconds\":" << timing.median
              << ",\"minimum_seconds\":" << timing.minimum
              << ",\"p95_seconds\":" << timing.p95
              << ",\"maximum_seconds\":" << timing.maximum
              << ",\"million_hex8_per_second\":" << rate
              << ",\"logical_gb_per_second\":" << logicalBandwidth
              << ",\"workspace_retained_bytes\":"
              << (args.cold ? 0U : workspace.retainedLogicalBytes())
              << ",\"mesh_hash\":" << hash
              << ",\"threshold_million_hex8_per_second\":"
              << args.thresholdMillionCells
              << ",\"gate_passed\":" << (passed ? "true" : "false")
              << "}\n";
    return passed ? EXIT_SUCCESS : 10;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(parseArguments(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "invalid benchmark request: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
