// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sdk/JobService.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <chrono>
#include <iostream>
int main() {
    using namespace paacs;
    sdk::JobService service({4, 64 * 1024 * 1024, 128 * 1024 * 1024});
    algo::PrimitiveMeshingOptions options;
    options.elements = algo::PrimitiveElements::Tet4;
    options.targetSize = .4;
    options.maximumChordError = .03;
    auto submitted = service.submit(sdk::primitiveTask({}, options));
    if (!submitted.job)
        return 1;
    // In a real UI use its timer/event loop to poll; do not wait on that thread.
    while (!submitted.job.waitFor(std::chrono::milliseconds(50)))
        std::cout << sdk::jobStageName(submitted.job.snapshot().stage) << '\n';
    auto result = submitted.job.result();
    if (!result) {
        std::cerr << submitted.job.snapshot().diagnostic.data() << '\n';
        return 1;
    }
    std::cout << result->summary.nodes << " nodes; " << result->summary.cells << " cells\n";
    const auto &report = std::get<algo::PrimitiveMeshingResult>(result->report);
    std::cout << "Minimum quality " << report.minimumQuality << "; warnings "
              << report.qualityWarning << '\n';
    // Keep this shared const result for rendering. It survives the service and
    // does not copy the mesh. Persist an archive from a host I/O thread if needed.
    return 0;
}
