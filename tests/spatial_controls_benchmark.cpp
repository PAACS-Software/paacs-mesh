// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "geometry_checks.hpp"
#include "paacs/algo/AdaptiveBoxMesher.hpp"
#include "paacs/algo/MeshSizeAudit.hpp"
#include "spatial_control_fixture.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
int main(int argc, char **argv) {
    try {
        using namespace paacs;
        if (argc != 3)
            return 2;
        const std::string name = argv[1];
        const int threads = std::stoi(argv[2]);
        auto field = spatial_test::make(name);
        algo::AdaptiveBoxRequest o;
        o.field = &field;
        o.threads = threads;
        o.maximumSizingWork = 1'000'000;
        mesh::Mesh mesh;
        const auto start = std::chrono::steady_clock::now();
        auto result = algo::meshAdaptiveBox(o, mesh);
        const double seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        if (!result.published) {
            std::cerr << algo::adaptiveBoxReasonName(result.reason) << '\n';
            return 3;
        }
        auto audit = algo::auditMeshSize(mesh, field);
        if (!audit.complete)
            return 4;
        std::cout << std::setprecision(17) << "{\"case\":\"" << name << "\",\"threads\":" << threads
                  << ",\"nodes\":" << mesh.nodeCount() << ",\"tetrahedra\":" << result.tetrahedra
                  << ",\"output_bytes\":" << result.outputBytes
                  << ",\"working_bytes\":" << result.workingBytes
                  << ",\"field_bytes\":" << field.logicalBytes()
                  << ",\"field_nodes\":" << field.nodeCount()
                  << ",\"sizing_queries\":" << result.sizingQueries
                  << ",\"sizing_work\":" << result.sizingWork
                  << ",\"minimum_quality\":" << result.minimumQuality
                  << ",\"size_ratio\":" << audit.maximumConservativeRatio
                  << ",\"seconds\":" << seconds << ",\"hash\":\"" << geometry_test::digest(mesh)
                  << "\"}\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
