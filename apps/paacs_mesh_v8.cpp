// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/sdk/MeshingSession.hpp"
#include "paacs/sdk/Version.hpp"

#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
constexpr int kUsage = 2;
constexpr int kInput = 3;
constexpr int kMeshing = 4;
constexpr int kAudit = 5;
constexpr int kOutput = 6;
constexpr int kDifferent = 7;

void usage() {
    std::cerr
        << "PAACS Mesh V8 " << paacs::sdk::kApiVersion << "\n"
        << "usage:\n"
        << "  paacs-mesh-v8 --version\n"
        << "  paacs-mesh-v8 hex NX NY NZ OUTPUT [--threads N] [--full-audit]\n"
        << "  paacs-mesh-v8 audit ARCHIVE\n"
        << "  paacs-mesh-v8 validate ARCHIVE\n"
        << "  paacs-mesh-v8 compare FIRST_ARCHIVE SECOND_ARCHIVE\n";
}

std::uint32_t parseCount(const char* text) {
    const unsigned long value = std::stoul(text);
    if (value == 0 || value > 1'000'000UL)
        throw std::invalid_argument("cell count outside [1,1000000]");
    return static_cast<std::uint32_t>(value);
}

int runHex(int argc, char** argv) {
    if (argc < 6) { usage(); return kUsage; }
    paacs::algo::MappedHexBlockRequest request;
    request.geometry.cellsU = parseCount(argv[2]);
    request.geometry.cellsV = parseCount(argv[3]);
    request.geometry.cellsW = parseCount(argv[4]);
    request.geometry.axisU = {1,0,0};
    request.geometry.axisV = {0,1,0};
    request.geometry.axisW = {0,0,1};
    request.geometry.volumeClassification = {3,1};
    for (std::size_t i=0;i<request.geometry.boundaryClassifications.size();++i)
        request.geometry.boundaryClassifications[i] =
            {2,static_cast<int>(101+i)};
    int threads = 0;
    for (int i=6;i<argc;++i) {
        const std::string option = argv[i];
        if (option == "--full-audit") {
            request.validation =
                paacs::algo::MappedHexValidationMode::FullTopologyAudit;
        } else if (option == "--threads" && i + 1 < argc) {
            threads = std::stoi(argv[++i]);
            if (threads <= 0) throw std::invalid_argument("threads must be positive");
        } else {
            throw std::invalid_argument("unknown hex option: " + option);
        }
    }
    paacs::sdk::SessionOptions sessionOptions;
    sessionOptions.threads = threads;
    paacs::sdk::MeshingSession session(sessionOptions);
    paacs::mesh::Mesh mesh;
    const auto result = session.generateMappedHex(request, mesh);
    if (!result.completed) {
        std::cerr << result.diagnostic << '\n';
        return kMeshing;
    }
    const auto quality = paacs::algo::analyzeMixedMeshQuality(mesh);
    if (!quality.valid) {
        std::cerr << "generated mesh failed mixed-cell quality audit\n";
        return kAudit;
    }
    const auto archived = paacs::sdk::writeMeshArchive(mesh, argv[5]);
    if (!archived.completed) {
        std::cerr << "archive write failed\n";
        return kOutput;
    }
    std::cout << std::setprecision(17)
              << "{\"status\":\"complete\",\"api\":\""
              << paacs::sdk::kApiVersion << "\",\"nodes\":"
              << mesh.nodeCount() << ",\"cells\":" << mesh.cellCount()
              << ",\"minimum_quality\":" << quality.minimumQuality
              << ",\"mesh_hash\":" << paacs::mesh::meshHash(mesh)
              << ",\"archive_bytes\":" << archived.bytes << "}\n";
    return 0;
}

bool sameMesh(const paacs::mesh::Mesh& first,
              const paacs::mesh::Mesh& second) {
    if (first.nodeCount() != second.nodeCount() ||
        first.cellCount() != second.cellCount() ||
        first.connectivityCount() != second.connectivityCount())
        return false;
    for (paacs::Index node = 0; node < first.nodeCount(); ++node) {
        const auto& a = first.node(node);
        const auto& b = second.node(node);
        if (a.x != b.x || a.y != b.y || a.z != b.z ||
            first.nodeClass(node) != second.nodeClass(node))
            return false;
    }
    for (paacs::Index cell = 0; cell < first.cellCount(); ++cell) {
        const auto& a = first.cell(cell);
        const auto& b = second.cell(cell);
        if (a.type != b.type || a.onEntity != b.onEntity)
            return false;
        const int count = paacs::mesh::nodesPerCell(a.type);
        const paacs::Index* an = first.cellNodes(cell);
        const paacs::Index* bn = second.cellNodes(cell);
        for (int local = 0; local < count; ++local)
            if (an[local] != bn[local]) return false;
    }
    return true;
}

int runAudit(int argc, char** argv) {
    if (argc != 3) { usage(); return kUsage; }
    paacs::mesh::Mesh mesh;
    const auto loaded = paacs::sdk::readMeshArchive(argv[2], mesh);
    if (!loaded.completed) {
        std::cerr << "archive read failed\n";
        return kInput;
    }
    const auto quality = paacs::algo::analyzeMixedMeshQuality(mesh);
    std::cout << std::setprecision(17)
              << "{\"valid\":" << (quality.valid ? "true" : "false")
              << ",\"nodes\":" << mesh.nodeCount()
              << ",\"cells\":" << mesh.cellCount()
              << ",\"tet4\":" << quality.tet4
              << ",\"pyramid5\":" << quality.pyramid5
              << ",\"prism6\":" << quality.prism6
              << ",\"hex8\":" << quality.hex8
              << ",\"minimum_quality\":" << quality.minimumQuality
              << ",\"mesh_hash\":" << paacs::mesh::meshHash(mesh)
              << "}\n";
    return quality.valid ? 0 : kAudit;
}
int runCompare(int argc, char** argv) {
    if (argc != 4) { usage(); return kUsage; }
    paacs::mesh::Mesh first;
    paacs::mesh::Mesh second;
    const auto firstLoaded = paacs::sdk::readMeshArchive(argv[2], first);
    const auto secondLoaded = paacs::sdk::readMeshArchive(argv[3], second);
    if (!firstLoaded.completed || !secondLoaded.completed) {
        std::cerr << "archive read failed\n";
        return kInput;
    }
    const bool same = sameMesh(first, second);
    std::cout << "{\"same\":" << (same ? "true" : "false")
              << ",\"first_hash\":" << paacs::mesh::meshHash(first)
              << ",\"second_hash\":" << paacs::mesh::meshHash(second)
              << ",\"first_nodes\":" << first.nodeCount()
              << ",\"second_nodes\":" << second.nodeCount()
              << ",\"first_cells\":" << first.cellCount()
              << ",\"second_cells\":" << second.cellCount()
              << "}\n";
    return same ? 0 : kDifferent;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--version") {
            std::cout << paacs::sdk::kApiVersion << '\n';
            return 0;
        }
        if (argc >= 2 && std::string(argv[1]) == "hex")
            return runHex(argc, argv);
        if (argc >= 2 && (std::string(argv[1]) == "audit" ||
                          std::string(argv[1]) == "validate"))
            return runAudit(argc, argv);
        if (argc >= 2 && std::string(argv[1]) == "compare")
            return runCompare(argc, argv);
        usage();
        return kUsage;
    } catch (const std::exception& error) {
        std::cerr << "invalid request: " << error.what() << '\n';
        return kInput;
    }
}
