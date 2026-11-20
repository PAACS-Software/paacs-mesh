// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sdk/MeshArchive.hpp"
#include "primitive_checks.hpp"
#include <charconv>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
using namespace primitive_test;
int main(int argc, char **argv) {
    try {
        std::string name;
        algo::PrimitiveElements family = algo::PrimitiveElements::HexDominant;
        algo::PrimitiveConeStrategy strategy = algo::PrimitiveConeStrategy::Automatic;
        algo::PrimitiveCylinderStrategy cylinderStrategy =
            algo::PrimitiveCylinderStrategy::Automatic;
        algo::PrimitiveCircularStrategy circularStrategy = algo::PrimitiveCircularStrategy::Automatic;
        int threads = 1;
        double factor = 1;
        bool single = false, list = false;
        std::filesystem::path archives;
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--elements" && i + 1 < argc) {
                const std::string value = argv[++i];
                bool found = false;
                for (auto e : {algo::PrimitiveElements::HexDominant, algo::PrimitiveElements::Tet4,
                               algo::PrimitiveElements::Tet10})
                    if (value == algo::primitiveElementsName(e)) {
                        family = e;
                        found = true;
                    }
                require(found, "invalid elements");
            } else if (a == "--swept")
                strategy = algo::PrimitiveConeStrategy::Swept;
            else if (a == "--cylinder-swept")
                cylinderStrategy = algo::PrimitiveCylinderStrategy::Swept;
            else if (a == "--circular-swept")
                circularStrategy = algo::PrimitiveCircularStrategy::Swept;
            else if (a == "--case" && i + 1 < argc)
                name = argv[++i];
            else if (a == "--threads" && i + 1 < argc)
                threads = std::stoi(argv[++i]);
            else if (a == "--size-factor" && i + 1 < argc)
                factor = std::stod(argv[++i]);
            else if (a == "--archives" && i + 1 < argc)
                archives = argv[++i];
            else if (a == "--single")
                single = true;
            else if (a == "--list")
                list = true;
            else
                throw std::runtime_error("invalid argument");
        }
        if (!archives.empty())
            std::filesystem::create_directories(archives);
        int failures = 0, count = 0;
        for (const auto &c : cases()) {
            if (!name.empty() && c.name != name)
                continue;
            ++count;
            if (list) {
                std::cout << c.name << '\n';
                continue;
            }
            try {
                const auto begin = std::chrono::steady_clock::now();
                algo::PrimitiveMeshingOptions o;
                o.targetSize = c.size * factor;
                o.threads = threads;
                o.elements = family;
                o.coneStrategy = strategy;
                o.cylinderStrategy = cylinderStrategy;
                o.circularStrategy = circularStrategy;
                mesh::Mesh m;
                const auto r = algo::meshPrimitive(c.shape, o, m);
                require(r.published, algo::primitiveReasonName(r.reason));
                require(r.audit.valid, "incidence audit");
                const auto proof = check(c.shape, m, r);
                const auto digest = hash(m);
                if (!single) {
                    auto other = o;
                    other.threads = threads == 1 ? 4 : 1;
                    mesh::Mesh repeated;
                    require(algo::meshPrimitive(c.shape, other, repeated).published,
                            "repeat failed");
                    require(hash(repeated) == digest, "thread nondeterminism");
                }
                if (!archives.empty())
                    require(sdk::writeMeshArchive(m, (archives / (c.name + ".paacsmesh")).string())
                                .completed,
                            "write archive");
                const double wall =
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
                std::cout << std::setprecision(17) << "{\"case\":\"" << c.name
                          << "\",\"passed\":true,\"threads\":" << threads
                          << ",\"nodes\":" << r.nodes
                          << ",\"volume_cells\":" << r.hexes + r.pyramids + r.tet4s + r.tet10s
                          << ",\"element_family\":\"" << algo::primitiveElementsName(family) << "\""
                          << ",\"tet4s\":" << r.tet4s << ",\"tet10s\":" << r.tet10s
                          << ",\"low_quality_cells\":" << r.quality.lowQualityCells
                          << ",\"graded_cone\":" << (r.gradedCone ? "true" : "false")
                          << ",\"uniform_cylinder\":" << (r.uniformCylinder ? "true" : "false")
                          << ",\"uniform_circular\":" << (r.uniformCircular ? "true" : "false")
                          << ",\"hexes\":" << r.hexes << ",\"pyramids\":" << r.pyramids
                          << ",\"output_bytes\":" << r.outputBytes
                          << ",\"planned_working_bytes\":" << r.plannedWorkingBytes
                          << ",\"volume\":" << proof.volume
                          << ",\"expected_volume\":" << r.analyticVolume
                          << ",\"relative_volume_error\":" << r.relativeVolumeError
                          << ",\"minimum_quality\":" << r.minimumQuality
                          << ",\"chord_error_bound\":" << r.chordErrorBound
                          << ",\"measured_surface_error\":" << proof.maximumSurfaceError
                          << ",\"boundary_euler\":" << proof.euler
                          << ",\"generation_seconds\":" << r.generationSeconds
                          << ",\"validation_seconds\":" << r.validationSeconds
                          << ",\"audit_seconds\":" << r.auditSeconds
                          << ",\"audit_temporary_bytes\":" << r.audit.temporaryBytes
                          << ",\"wall_seconds\":" << wall << ",\"hash\":\"" << digest
                          << "\",\"thread_comparison\":" << (single ? "false" : "true") << "}\n";
            } catch (const std::exception &e) {
                ++failures;
                std::cerr << c.name << ": " << e.what() << '\n';
            }
        }
        require(count > 0, "unknown case");
        return failures ? 1 : 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 2;
    }
}
