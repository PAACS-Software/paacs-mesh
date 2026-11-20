// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/algo/ModelGeometryProviders.hpp"
#include "paacs/algo/PreOcctMeshingCore.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/sdk/MeshingSession.hpp"
#include "paacs/sdk/Version.hpp"
#include "paacs/sizing/SizeField.hpp"

#include "NumberArgument.hpp"
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
constexpr int kUsage = 2;
constexpr int kInput = 3;
constexpr int kMeshing = 4;
constexpr int kAudit = 5;
constexpr int kOutput = 6;
constexpr int kDifferent = 7;

void usage() {
    std::cerr << "PAACS Mesh " << paacs::sdk::kApiVersion << "\n"
              << "usage:\n"
              << "  paacs-mesh --version\n"
              << "  paacs-mesh hex NX NY NZ OUTPUT [--threads N] [--full-audit]\n"
              << "  paacs-mesh tet-box SIZE OUTPUT [--tet10] [--threads N] [--refine N]\n"
              << "  paacs-mesh tet-cylinder RADIUS HEIGHT SIZE OUTPUT [--tet10] [--threads N] "
                 "[--refine N]\n"
              << "  paacs-mesh audit ARCHIVE\n"
              << "  paacs-mesh validate ARCHIVE\n"
              << "  paacs-mesh compare FIRST_ARCHIVE SECOND_ARCHIVE\n";
}

template <class T> T parseNumber(const char *text) {
    const std::string_view input(text);
    T value{};
    if (!paacs::cli::parseNumber(input, value))
        throw std::invalid_argument("invalid number: " + std::string(input));
    return value;
}
std::uint32_t parseCount(const char *text) {
    const auto value = parseNumber<std::uint32_t>(text);
    if (value == 0 || value > 1'000'000U)
        throw std::invalid_argument("cell count outside [1,1000000]");
    return value;
}
double parsePositive(const char *text) {
    const auto value = parseNumber<double>(text);
    if (!(value > 0.0) || !std::isfinite(value))
        throw std::invalid_argument("expected a finite positive number");
    return value;
}

int runHex(int argc, char **argv) {
    if (argc < 6) {
        usage();
        return kUsage;
    }
    paacs::algo::MappedHexBlockRequest request;
    request.geometry.cellsU = parseCount(argv[2]);
    request.geometry.cellsV = parseCount(argv[3]);
    request.geometry.cellsW = parseCount(argv[4]);
    request.geometry.axisU = {1, 0, 0};
    request.geometry.axisV = {0, 1, 0};
    request.geometry.axisW = {0, 0, 1};
    request.geometry.volumeClassification = {3, 1};
    for (std::size_t i = 0; i < request.geometry.boundaryClassifications.size(); ++i)
        request.geometry.boundaryClassifications[i] = {2, static_cast<int>(101 + i)};
    int threads = 0;
    for (int i = 6; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--full-audit") {
            request.validation = paacs::algo::MappedHexValidationMode::FullTopologyAudit;
        } else if (option == "--threads" && i + 1 < argc) {
            threads = parseNumber<int>(argv[++i]);
            if (threads <= 0)
                throw std::invalid_argument("threads must be positive");
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
    std::cout << std::setprecision(17) << "{\"status\":\"complete\",\"api\":\""
              << paacs::sdk::kApiVersion << "\",\"nodes\":" << mesh.nodeCount()
              << ",\"cells\":" << mesh.cellCount()
              << ",\"minimum_quality\":" << quality.minimumQuality
              << ",\"mesh_hash\":" << paacs::mesh::meshHash(mesh)
              << ",\"archive_bytes\":" << archived.bytes << "}\n";
    return 0;
}

int runTetPrimitive(int argc, char **argv, bool cylinder) {
    const int required = cylinder ? 6 : 4;
    if (argc < required) {
        usage();
        return kUsage;
    }
    const double radius = cylinder ? parsePositive(argv[2]) : 0.0;
    const double height = cylinder ? parsePositive(argv[3]) : 0.0;
    const double size = parsePositive(argv[cylinder ? 4 : 2]);
    const char *output = argv[cylinder ? 5 : 3];
    if (!(size > 0.0) || (cylinder && (!(radius > 0.0) || !(height > 0.0))))
        throw std::invalid_argument("positive dimensions and size required");
    bool tet10 = false;
    int threads = 0;
    std::uint64_t insertions = 0;
    for (int i = required; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--tet10")
            tet10 = true;
        else if (option == "--threads" && i + 1 < argc) {
            threads = parseNumber<int>(argv[++i]);
            if (threads <= 0)
                throw std::invalid_argument("threads must be positive");
        } else if (option == "--refine" && i + 1 < argc) {
            insertions = parseNumber<std::uint64_t>(argv[++i]);
        } else
            throw std::invalid_argument("unknown Tet option: " + option);
    }

    auto kernel = paacs::geo::makeBuiltinKernel();
    auto model = cylinder ? kernel->makeCylinder({0, 0, 0}, radius, height)
                          : kernel->makeBox({0, 0, 0}, {1, 1, 1});
    auto ownership = paacs::algo::buildModelBoundaryOwnershipProvider(*model);
    if (!ownership.built) {
        std::cerr << "model ownership preparation failed\n";
        return kMeshing;
    }
    paacs::algo::ModelTet10PlacementProvider placement(*model);
    paacs::sizing::MeshSpec spec;
    spec.field = std::make_shared<paacs::sizing::ConstantField>(size);
    spec.threads = threads;
    auto request =
        paacs::algo::makePreOcctTetRequest(tet10 ? paacs::algo::PreOcctTetOutput::CertifiedTet10
                                                 : paacs::algo::PreOcctTetOutput::QualityTet4,
                                           ownership.provider.stamp(), placement.stamp());
    request.pipeline.source.source.refinement.maxInsertions = insertions;
    if (insertions == 0) {
        request.pipeline.source.source.refinement.maxEdgeToTarget = 1.0e9;
        request.pipeline.source.source.refinement.maxRadiusEdge = 1.0e9;
    }
    paacs::mesh::Mesh mesh;
    const auto result = paacs::algo::meshPreparedTetModel(*model, spec, request, mesh,
                                                          ownership.provider, &placement);
    if (!result.published()) {
        std::cerr << "Tet production route stopped at phase " << static_cast<unsigned>(result.phase)
                  << " reason " << static_cast<unsigned>(result.reason) << '\n';
        return kMeshing;
    }
    const auto archived = paacs::sdk::writeMeshArchive(mesh, output);
    if (!archived.completed)
        return kOutput;
    std::cout << std::setprecision(17) << "{\"status\":\"complete\",\"api\":\""
              << paacs::sdk::kApiVersion << "\",\"nodes\":" << mesh.nodeCount()
              << ",\"cells\":" << mesh.cellCount()
              << ",\"tet4\":" << mesh.countCells(paacs::mesh::CellType::Tet4)
              << ",\"tet10\":" << mesh.countCells(paacs::mesh::CellType::Tet10)
              << ",\"minimum_radius_ratio\":"
              << (result.tetQuality.radiusRatio.n ? result.tetQuality.radiusRatio.min : 0.0)
              << ",\"mesh_hash\":" << result.strictHash << "}\n";
    return 0;
}

bool sameMesh(const paacs::mesh::Mesh &first, const paacs::mesh::Mesh &second) {
    if (first.nodeCount() != second.nodeCount() || first.cellCount() != second.cellCount() ||
        first.connectivityCount() != second.connectivityCount())
        return false;
    for (paacs::Index node = 0; node < first.nodeCount(); ++node) {
        const auto &a = first.node(node);
        const auto &b = second.node(node);
        if (a.x != b.x || a.y != b.y || a.z != b.z ||
            first.nodeClass(node) != second.nodeClass(node))
            return false;
    }
    for (paacs::Index cell = 0; cell < first.cellCount(); ++cell) {
        const auto &a = first.cell(cell);
        const auto &b = second.cell(cell);
        if (a.type != b.type || a.onEntity != b.onEntity)
            return false;
        const int count = paacs::mesh::nodesPerCell(a.type);
        const paacs::Index *an = first.cellNodes(cell);
        const paacs::Index *bn = second.cellNodes(cell);
        for (int local = 0; local < count; ++local)
            if (an[local] != bn[local])
                return false;
    }
    return true;
}

int runAudit(int argc, char **argv) {
    if (argc != 3) {
        usage();
        return kUsage;
    }
    paacs::mesh::Mesh mesh;
    const auto loaded = paacs::sdk::readMeshArchive(argv[2], mesh);
    if (!loaded.completed) {
        std::cerr << "archive read failed\n";
        return kInput;
    }
    const auto mixed = paacs::algo::analyzeMixedMeshQuality(mesh);
    const auto tet = paacs::mesh::analyzeTets(mesh);
    paacs::mesh::AuditOptions options;
    options.boundaryMode = paacs::mesh::BoundaryComparisonMode::WhenSurfacePresent;
    options.boundaryOrphanBudget = 0;
    options.boundaryUnexpectedBudget = 0;
    options.duplicateFacetsAreErrors = true;
    options.unusedNodesAreErrors = true;
    const auto audit = paacs::mesh::auditMesh(mesh, options);
    const auto tet10 = mesh.countCells(paacs::mesh::CellType::Tet10);
    const bool tet10VolumeOnly =
        tet10 > 0 && mixed.tet4 == 0 && mixed.pyramid5 == 0 && mixed.prism6 == 0 && mixed.hex8 == 0;
    const bool qualityValid =
        mixed.valid || (tet10VolumeOnly && tet.nTet10Invalid == 0 && tet.nInverted == 0);
    const bool valid = audit.passed() && qualityValid;
    const double minimumQuality =
        mixed.valid ? mixed.minimumQuality : (tet.radiusRatio.n ? tet.radiusRatio.min : 0.0);
    std::cout << std::setprecision(17) << "{\"valid\":" << (valid ? "true" : "false")
              << ",\"nodes\":" << mesh.nodeCount() << ",\"cells\":" << mesh.cellCount()
              << ",\"tet4\":" << mixed.tet4 << ",\"tet10\":" << tet10
              << ",\"pyramid5\":" << mixed.pyramid5 << ",\"prism6\":" << mixed.prism6
              << ",\"hex8\":" << mixed.hex8 << ",\"minimum_quality\":" << minimumQuality
              << ",\"invalid_tet10\":" << tet.nTet10Invalid
              << ",\"audit_passed\":" << (audit.passed() ? "true" : "false")
              << ",\"inverted\":" << audit.nInverted << ",\"nonfinite\":" << audit.nNonFinite
              << ",\"malformed\":" << audit.malformedCells
              << ",\"deferred\":" << audit.deferredAuditCells
              << ",\"orphans\":" << audit.boundaryOrphans
              << ",\"unexpected\":" << audit.boundaryUnexpected
              << ",\"duplicates\":" << audit.boundaryDuplicates
              << ",\"nonmanifold\":" << audit.nonManifoldFaces
              << ",\"mesh_hash\":" << paacs::mesh::meshHash(mesh) << "}\n";
    return valid ? 0 : kAudit;
}
int runCompare(int argc, char **argv) {
    if (argc != 4) {
        usage();
        return kUsage;
    }
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
              << ",\"second_cells\":" << second.cellCount() << "}\n";
    return same ? 0 : kDifferent;
}

} // namespace

int main(int argc, char **argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--version") {
            std::cout << paacs::sdk::kApiVersion << '\n';
            return 0;
        }
        if (argc >= 2 && std::string(argv[1]) == "hex")
            return runHex(argc, argv);
        if (argc >= 2 && std::string(argv[1]) == "tet-box")
            return runTetPrimitive(argc, argv, false);
        if (argc >= 2 && std::string(argv[1]) == "tet-cylinder")
            return runTetPrimitive(argc, argv, true);
        if (argc >= 2 && (std::string(argv[1]) == "audit" || std::string(argv[1]) == "validate"))
            return runAudit(argc, argv);
        if (argc >= 2 && std::string(argv[1]) == "compare")
            return runCompare(argc, argv);
        usage();
        return kUsage;
    } catch (const std::exception &error) {
        std::cerr << "invalid request: " << error.what() << '\n';
        return kInput;
    }
}
