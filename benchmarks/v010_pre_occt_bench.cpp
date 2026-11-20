// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/CanonicalCurveMesher.hpp"
#include "paacs/algo/ModelGeometryProviders.hpp"
#include "paacs/algo/PreOcctMeshingCore.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/sizing/SizeField.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string_view>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;

double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    return values[values.size() / 2U];
}

struct CaseResult {
    double seconds = 0.0;
    paacs::Index nodes = 0;
    paacs::Index cells = 0;
    std::uint64_t hash = 0;
    double minimumRadiusRatio = 0.0;
};

CaseResult runCase(const paacs::geo::Model& model,
                   double targetSize,
                   paacs::algo::PreOcctTetOutput output,
                   int runs) {
    auto ownership =
        paacs::algo::buildModelBoundaryOwnershipProvider(model);
    if (!ownership.built) throw std::runtime_error("ownership preparation failed");
    paacs::algo::ModelTet10PlacementProvider placement(model);
    paacs::sizing::MeshSpec spec;
    spec.field = std::make_shared<paacs::sizing::ConstantField>(targetSize);
    spec.threads = 1;
    auto request = paacs::algo::makePreOcctTetRequest(
        output, ownership.provider.stamp(), placement.stamp());
    request.pipeline.source.source.refinement.maxInsertions = 0;
    request.pipeline.source.source.refinement.maxEdgeToTarget = 1.0e9;
    request.pipeline.source.source.refinement.maxRadiusEdge = 1.0e9;

    auto once = [&]() {
        paacs::mesh::Mesh mesh;
        const auto begin = Clock::now();
        auto outcome = paacs::algo::meshPreparedTetModel(
            model, spec, request, mesh, ownership.provider, &placement);
        const double seconds = std::chrono::duration<double>(
            Clock::now() - begin).count();
        if (!outcome.published())
            throw std::runtime_error("production Tet route stopped");
        CaseResult result;
        result.seconds = seconds;
        result.nodes = mesh.nodeCount();
        result.cells = mesh.cellCount();
        result.hash = outcome.strictHash;
        if (outcome.tetQuality.radiusRatio.n)
            result.minimumRadiusRatio = outcome.tetQuality.radiusRatio.min;
        return result;
    };

    (void)once();
    std::vector<double> times;
    CaseResult last;
    times.reserve(static_cast<std::size_t>(runs));
    for (int i = 0; i < runs; ++i) {
        last = once();
        times.push_back(last.seconds);
    }
    last.seconds = median(std::move(times));
    return last;
}

void printCase(std::string_view name, const CaseResult& result) {
    std::cout << std::setprecision(17)
              << "{\"case\":\"" << name << "\",\"seconds\":"
              << result.seconds << ",\"nodes\":" << result.nodes
              << ",\"cells\":" << result.cells
              << ",\"minimum_radius_ratio\":"
              << result.minimumRadiusRatio << ",\"hash\":"
              << result.hash << "}\n";
}
}  // namespace

int main() {
    try {
        auto kernel = paacs::geo::makeBuiltinKernel();
        auto box = kernel->makeBox({0,0,0}, {1,1,1});
        auto cylinder = kernel->makeCylinder({0,0,0}, 1.0, 1.5);
        printCase("box_quality_tet4", runCase(
            *box, 0.45, paacs::algo::PreOcctTetOutput::QualityTet4, 5));
        printCase("box_certified_tet10", runCase(
            *box, 0.45, paacs::algo::PreOcctTetOutput::CertifiedTet10, 5));
        printCase("cylinder_quality_tet4", runCase(
            *cylinder, 0.45, paacs::algo::PreOcctTetOutput::QualityTet4, 5));
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
