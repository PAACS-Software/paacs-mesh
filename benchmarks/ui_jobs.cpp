// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sdk/JobService.hpp"
#include <chrono>
#include <iostream>
#include <stdexcept>
using namespace paacs;
using Clock = std::chrono::steady_clock;
int main() {
    try {
        for (bool computeHash : {false, true}) {
            sdk::JobServiceLimits limits;
            limits.computeResultHash = computeHash;
            sdk::JobService service(limits);
            for (double size : {.4, .2, .1}) {
                algo::PrimitiveShape shape;
                algo::PrimitiveMeshingOptions options;
                options.elements = algo::PrimitiveElements::Tet4;
                options.targetSize = size;
                options.maximumChordError = .01;
                for (int sample = -2; sample < 15; ++sample) {
                    mesh::Mesh direct;
                    auto begin = Clock::now();
                    auto native = algo::meshPrimitive(shape, options, direct);
                    auto end = Clock::now();
                    if (!native.published)
                        throw std::runtime_error("native meshing failed");
                    const auto hash = mesh::meshHash(direct);
                    auto startJob = Clock::now();
                    auto job = service.submit(sdk::primitiveTask(shape, options)).job;
                    if (!job || !job.waitFor(std::chrono::seconds(30)) || !job.result())
                        throw std::runtime_error("job failed");
                    auto endJob = Clock::now();
                    auto result = job.result();
                    if (hash != mesh::meshHash(result->mesh) ||
                        (computeHash && hash != result->summary.hash) ||
                        (computeHash != result->summary.hashAvailable))
                        throw std::runtime_error("mesh identity changed");
                    if (sample >= 0)
                        std::cout
                            << "{\"target\":" << size
                            << ",\"hash_enabled\":" << (computeHash ? "true" : "false")
                            << ",\"sample\":" << sample << ",\"native_ms\":"
                            << std::chrono::duration<double, std::milli>(end - begin).count()
                            << ",\"job_ms\":"
                            << std::chrono::duration<double, std::milli>(endJob - startJob).count()
                            << ",\"nodes\":" << result->summary.nodes
                            << ",\"cells\":" << result->summary.cells
                            << ",\"logical_bytes\":" << result->summary.logicalBytes
                            << ",\"hash\":\"" << hash << "\"}\n";
                }
            }
        }
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
