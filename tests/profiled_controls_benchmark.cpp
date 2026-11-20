// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "profile_controls_fixture.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace profile_test;
int main(int argc, char **argv) {
    try {
        if (argc < 7 || argc > 8)
            throw std::runtime_error("usage: benchmark mode tet4|tet10|hex sweep|volume local "
                                     "threads reuse[0|1] [archive]");
        const std::string family = argv[2], strategy = argv[3];
        if (family != "tet4" && family != "tet10" && family != "hex")
            throw std::runtime_error("family");
        if (strategy != "sweep" && strategy != "volume")
            throw std::runtime_error("strategy");
        auto s = shape();
        auto spec = specification(argv[1], 1, std::stod(argv[4]));
        const auto field = sizing::compileCappedSizeProfile(spec);
        if (!field.compiled)
            throw std::runtime_error("field refused");
        auto o = options(family == "tet4"    ? algo::PrimitiveElements::Tet4
                         : family == "tet10" ? algo::PrimitiveElements::Tet10
                                             : algo::PrimitiveElements::HexDominant,
                         strategy == "volume" ? algo::PrimitiveProfileStrategy::Volume
                                              : algo::PrimitiveProfileStrategy::Sweep);
        o.threads = std::stoi(argv[5]);
        o.reuseProfileSlabs = std::stoi(argv[6]) != 0;
        o.sizeProfile = &field.profile;
        mesh::Mesh m;
        const auto begin = std::chrono::steady_clock::now();
        const auto r = algo::meshPrimitive(s, o, m);
        const double seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        if (!r.published)
            throw std::runtime_error(algo::primitiveReasonName(r.reason));
        if (argc == 8 && !sdk::writeMeshArchive(m, argv[7]).completed)
            throw std::runtime_error("archive write");
        std::cout << std::setprecision(17) << "{\"mode\":\"" << argv[1] << "\",\"family\":\""
                  << family << "\",\"strategy\":\"" << strategy
                  << "\",\"local\":" << std::stod(argv[4]) << ",\"threads\":" << o.threads
                  << ",\"reuse\":" << o.reuseProfileSlabs << ",\"nodes\":" << r.nodes
                  << ",\"cells\":" << r.tet4s + r.tet10s + r.hexes << ",\"seconds\":" << seconds
                  << ",\"quality_min\":" << r.minimumQuality << ",\"volume\":" << r.meshVolume
                  << ",\"output_bytes\":" << r.outputBytes
                  << ",\"working_bytes\":" << r.plannedWorkingBytes
                  << ",\"delaunay_bytes\":" << r.profileDelaunayBytes
                  << ",\"built_slabs\":" << r.profileDelaunaySlabs
                  << ",\"reused_slabs\":" << r.profileReusedSlabs << ",\"hash\":\"" << std::hex
                  << mesh::meshHash(m) << "\"}\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
