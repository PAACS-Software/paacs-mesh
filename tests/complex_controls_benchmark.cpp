// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "complex_geometry_fixture.hpp"
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
int main(int argc, char **argv) {
    using namespace paacs;
    if (argc != 5)
        return 2;
    const std::string name = argv[1], mode = argv[2], family = argv[3];
    auto shapes = complex_test::shapes();
    auto it =
        std::find_if(shapes.begin(), shapes.end(), [&](const auto &s) { return s.name == name; });
    if (it == shapes.end() || (mode != "fine" && mode != "start") ||
        (family != "native" && family != "tet4"))
        return 2;
    auto model = geo::makePolyhedralModel(complex_test::input(*it));
    if (!model.built())
        return 3;
    sizing::CappedSizeSpecification spec;
    spec.length = std::hypot(it->sweep.x, it->sweep.y, it->sweep.z);
    spec.globalSize = .45;
    spec.controls = {{0, mode == "fine" ? spec.length : 0, .12, 1.25}};
    auto field = sizing::compileCappedSizeProfile(spec);
    if (!field.compiled)
        return 4;
    algo::GeometryMeshingOptions o;
    o.targetSize = it->target;
    o.threads = std::stoi(argv[4]);
    o.warningElementQuality = .15;
    algo::GeometryBodyControl c;
    c.solidTag = 1;
    c.sourceFaceTag = 1;
    c.sizeProfile = &field.profile;
    c.elements = family == "tet4" ? algo::GeometryElements::Tet4 : algo::GeometryElements::Native;
    o.bodyControls = {c};
    mesh::Mesh m;
    const auto begin = std::chrono::steady_clock::now();
    const auto r = algo::meshGeometryModel(*model.model, o, m);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    if (!r.published) {
        std::cerr << algo::geometryReasonName(r.reason) << '\n';
        return 5;
    }
    std::vector<std::uint8_t> bytes;
    if (!sdk::serializeMeshArchive(m, bytes).completed)
        return 6;
    std::uint64_t hash = 1469598103934665603ULL;
    for (auto b : bytes) {
        hash ^= b;
        hash *= 1099511628211ULL;
    }
    const auto &body = r.analysis.bodies[0];
    std::cout << std::setprecision(17) << "{\"case\":\"" << name << "\",\"mode\":\"" << mode
              << "\",\"family\":\"" << family << "\",\"threads\":" << o.threads
              << ",\"nodes\":" << m.nodeCount() << ",\"volume_cells\":" << r.quality.volumeCells
              << ",\"layers\":" << body.layers << ",\"output_bytes\":" << r.outputBytes
              << ",\"profile_working_bytes\":" << body.profileWorkingBytes
              << ",\"conversion_working_bytes\":" << body.conversion.workingBytes
              << ",\"audit_temporary_bytes\":" << r.audit.temporaryBytes
              << ",\"minimum_quality\":" << r.quality.minimumQuality
              << ",\"quality_warning\":" << (r.qualityWarning ? "true" : "false")
              << ",\"seconds\":" << seconds << ",\"hash\":\"" << hash << "\"}\n";
}
