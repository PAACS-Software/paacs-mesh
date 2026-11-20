// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "geometry_checks.hpp"
#include "geometry_fixtures.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>

using namespace geometry_test;
int main(int argc, char **argv) {
    std::string selected;
    int threads = 1;
    double factor = 1;
    bool compare = true, list = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--case" && i + 1 < argc)
            selected = argv[++i];
        else if (a == "--threads" && i + 1 < argc)
            threads = std::stoi(argv[++i]);
        else if (a == "--size-factor" && i + 1 < argc)
            factor = std::stod(argv[++i]);
        else if (a == "--single")
            compare = false;
        else if (a == "--list")
            list = true;
        else {
            std::cerr << "invalid argument\n";
            return 2;
        }
    }
    int failures = 0, count = 0;
    for (const auto &c : cases()) {
        if (!selected.empty() && c.name != selected)
            continue;
        ++count;
        if (list) {
            std::cout << c.name << '\n';
            continue;
        }
        try {
            const auto begin = std::chrono::steady_clock::now();
            auto built = geo::makePolyhedralModel(c.input);
            require(built.built(),
                    std::string("geometry: ") + geo::polyhedralReasonName(built.reason));
            require(std::abs(built.volume - c.volume) <= 1e-8 * c.volume,
                    "geometry analytic volume");
            algo::GeometryMeshingOptions options;
            options.targetSize = c.targetSize * factor;
            options.threads = threads;
            options.strategy = c.strategy;
            mesh::Mesh m;
            const auto result = algo::meshGeometryModel(*built.model, options, m);
            if (!result.published) {
                std::string error = algo::geometryReasonName(result.reason);
                for (const auto &b : result.analysis.bodies)
                    error += " engine=" + std::to_string(b.engineStatus) + "/" +
                             std::to_string(b.engineReason) + "/" + std::to_string(b.enginePhase);
                for (const auto &b : result.analysis.bodies)
                    error += " recovery=" + std::to_string(b.recoveryKind) + "/" +
                             std::to_string(b.recoveryDetail);
                throw std::runtime_error(error);
            }
            require(result.audit.valid && result.analysis.bodies.size() == 1, "boundary audit");
            require(result.analysis.bodies[0].strategy == c.expected, "wrong strategy");
            const double volume = positiveVolume(m);
            require(std::abs(volume - c.volume) <= 2e-8 * c.volume, "mesh analytic volume");
            Vec3 lo, hi;
            built.model->solids()[0]->bbox(lo, hi);
            checkClassification(*built.model, m, (hi - lo).norm() * 1e-7);
            const auto quality = algo::analyzeMixedMeshQuality(m);
            require(quality.valid && quality.invertedCells == 0, "quality audit");
            const auto hash = digest(m);
            if (compare) {
                auto alternate = options;
                alternate.threads = threads == 1 ? 4 : 1;
                mesh::Mesh repeated;
                require(algo::meshGeometryModel(*built.model, alternate, repeated).published,
                        "parallel repeat failed");
                require(digest(repeated) == hash, "thread-count nondeterminism");
            }
            const double wall =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
            std::cout << std::setprecision(17) << "{\"case\":\"" << c.name
                      << "\",\"passed\":true,\"strategy\":\""
                      << algo::geometryStrategyName(c.expected) << "\",\"threads\":" << threads
                      << ",\"nodes\":" << m.nodeCount() << ",\"cells\":" << m.cellCount()
                      << ",\"volume_cells\":" << quality.volumeCells
                      << ",\"output_bytes\":" << result.outputBytes << ",\"volume\":" << volume
                      << ",\"expected_volume\":" << c.volume
                      << ",\"minimum_quality\":" << quality.minimumQuality
                      << ",\"mean_quality\":" << quality.meanQuality
                      << ",\"analysis_seconds\":" << result.analysisSeconds
                      << ",\"generation_seconds\":" << result.generationSeconds
                      << ",\"audit_seconds\":" << result.auditSeconds
                      << ",\"audit_temporary_bytes\":" << result.audit.temporaryBytes
                      << ",\"wall_seconds\":" << wall << ",\"hash\":\"" << hash
                      << "\",\"thread_comparison\":" << (compare ? "true" : "false") << "}\n";
        } catch (const std::exception &e) {
            ++failures;
            std::cerr << c.name << ": " << e.what() << '\n';
            std::cout << "{\"case\":\"" << c.name << "\",\"passed\":false,\"error\":\"" << e.what()
                      << "\"}\n";
        }
    }
    if (!count) {
        std::cerr << "unknown case\n";
        return 2;
    }
    return failures ? 1 : 0;
}
