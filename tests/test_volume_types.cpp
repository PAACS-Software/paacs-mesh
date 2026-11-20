// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/Meshers.hpp"
#include "paacs/algo/StructuredPolyhedralMesher.hpp"
#include "paacs/mesh/PolyhedralMesh.hpp"
#include "paacs/mesh/MeshAudit.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

extern int gFailures;
static const char* gVolumeTypeCurrent = "";
#define VOLUME_TYPE_TEST(name)             \
    do {                                   \
        gVolumeTypeCurrent = name;         \
        std::printf("[test] %s\n", name); \
    } while (0)
#define VOLUME_TYPE_CHECK(cond)                                      \
    do {                                                             \
        if (!(cond)) {                                               \
            std::printf("  FAIL %s:%d  %s  (in %s)\n", __FILE__,   \
                        __LINE__, #cond, gVolumeTypeCurrent);         \
            ++gFailures;                                             \
        }                                                            \
    } while (0)

namespace {
namespace algo = paacs::algo;
namespace mesh = paacs::mesh;
using paacs::Index;
using paacs::Vec3;

algo::StructuredBlockGeometry block(std::uint32_t u, std::uint32_t v,
                                    std::uint32_t w) {
    algo::StructuredBlockGeometry out;
    out.cellsU = u;
    out.cellsV = v;
    out.cellsW = w;
    out.volumeClassification = {3, 41};
    for (std::size_t i = 0; i < out.boundaryClassifications.size(); ++i)
        out.boundaryClassifications[i] =
            {2, static_cast<int>(101 + i)};
    return out;
}

bool sameClass(const mesh::Classification& a,
               const mesh::Classification& b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}

std::string readFile(const std::filesystem::path& path) {
    std::ifstream input(path);
    return {std::istreambuf_iterator<char>(input),
            std::istreambuf_iterator<char>()};
}

void mappedHexUnitCube() {
    VOLUME_TYPE_TEST("hex.v5 mapped unit block publishes Hex8 and Quad4");
    algo::MappedHexBlockRequest request;
    request.geometry = block(1, 1, 1);
    mesh::Mesh destination;
    const auto outcome = algo::generateMappedHexBlock(request, destination);
    VOLUME_TYPE_CHECK(outcome.published());
    VOLUME_TYPE_CHECK(!outcome.destinationUnchanged);
    VOLUME_TYPE_CHECK(outcome.metrics.nodes == 8);
    VOLUME_TYPE_CHECK(outcome.metrics.hexes == 1);
    VOLUME_TYPE_CHECK(outcome.metrics.boundaryQuads == 6);
    VOLUME_TYPE_CHECK(destination.nodeCount() == 8);
    VOLUME_TYPE_CHECK(destination.countCells(mesh::CellType::Hex8) == 1);
    VOLUME_TYPE_CHECK(destination.countCells(mesh::CellType::Quad4) == 6);
    const std::array<Index, 8> expected{0, 1, 3, 2, 4, 5, 7, 6};
    const Index* conn = destination.cellNodes(0);
    for (std::size_t i = 0; i < expected.size(); ++i)
        VOLUME_TYPE_CHECK(conn[i] == expected[i]);
    VOLUME_TYPE_CHECK(outcome.metrics.audit.valid);
    VOLUME_TYPE_CHECK(outcome.metrics.audit.boundaryFaces == 6);
    VOLUME_TYPE_CHECK(outcome.metrics.audit.interiorFaces == 0);
    VOLUME_TYPE_CHECK(std::abs(outcome.metrics.audit.totalVolume - 1.0) <
                      1.0e-12);
    VOLUME_TYPE_CHECK(outcome.metrics.audit.minimumScaledJacobian > 0.999999);
    for (Index i = 0; i < destination.nodeCount(); ++i)
        VOLUME_TYPE_CHECK(destination.nodeClass(i).dim == 0);

    const auto base = std::filesystem::temp_directory_path() /
        ("paacs_v5_hex_" + std::to_string(
            static_cast<unsigned long long>(mesh::meshHash(destination))));
    const auto vtk = base.string() + ".vtk";
    const auto msh = base.string() + ".msh";
    paacs::io::writeVtk(destination, vtk, false);
    paacs::io::writeMsh(destination, msh);
    const auto vtkText = readFile(vtk);
    const auto mshText = readFile(msh);
    VOLUME_TYPE_CHECK(vtkText.find("CELL_TYPES 7") != std::string::npos);
    VOLUME_TYPE_CHECK(vtkText.find("\n12\n") != std::string::npos);
    VOLUME_TYPE_CHECK(mshText.find(" 5 1\n") != std::string::npos);
    std::filesystem::remove(vtk);
    std::filesystem::remove(msh);
}

void mappedHexSkewedBlock() {
    VOLUME_TYPE_TEST("hex.v5 skewed affine block has positive sampled Jacobians");
    algo::MappedHexBlockRequest request;
    request.geometry = block(2, 3, 4);
    request.geometry.origin = {-3.0, 5.0, 7.0};
    request.geometry.axisU = {2.0, 0.0, 0.0};
    request.geometry.axisV = {0.5, 3.0, 0.0};
    request.geometry.axisW = {0.2, 0.1, 4.0};
    request.validation = algo::MappedHexValidationMode::FullTopologyAudit;
    mesh::Mesh destination;
    const auto outcome = algo::generateMappedHexBlock(request, destination);
    VOLUME_TYPE_CHECK(outcome.published());
    VOLUME_TYPE_CHECK(destination.nodeCount() == 60);
    VOLUME_TYPE_CHECK(destination.countCells(mesh::CellType::Hex8) == 24);
    VOLUME_TYPE_CHECK(destination.countCells(mesh::CellType::Quad4) == 52);
    VOLUME_TYPE_CHECK(outcome.metrics.audit.valid);
    VOLUME_TYPE_CHECK(outcome.metrics.audit.interiorFaces == 46);
    VOLUME_TYPE_CHECK(outcome.metrics.audit.boundaryFaces == 52);
    VOLUME_TYPE_CHECK(outcome.metrics.audit.minimumJacobian > 0.0);
    VOLUME_TYPE_CHECK(outcome.metrics.audit.minimumScaledJacobian > 0.0);
    VOLUME_TYPE_CHECK(std::abs(outcome.metrics.audit.totalVolume - 24.0) <
                      1.0e-10);
    VOLUME_TYPE_CHECK(sameClass(destination.nodeClass(0),
                                request.geometry.vertexClassifications[0]));
    VOLUME_TYPE_CHECK(sameClass(destination.nodeClass(1),
                                request.geometry.edgeClassifications[0]));
    VOLUME_TYPE_CHECK(sameClass(destination.nodeClass(4),
                                request.geometry.boundaryClassifications[
                                    static_cast<std::size_t>(
                                        algo::StructuredBlockFace::WMin)]));
    VOLUME_TYPE_CHECK(sameClass(destination.nodeClass(16),
                                request.geometry.volumeClassification));
}

void mappedHexContractsAndAudit() {
    VOLUME_TYPE_TEST("hex.v5 preflight and audit reject unsafe blocks atomically");
    mesh::Mesh destination;
    (void)destination.addNode({-1, -2, -3}, {0, 9});
    const auto before = mesh::meshHash(destination);
    {
        algo::MappedHexBlockRequest request;
        request.geometry = block(0, 1, 1);
        const auto out = algo::generateMappedHexBlock(request, destination);
        VOLUME_TYPE_CHECK(!out.published());
        VOLUME_TYPE_CHECK(out.reason ==
                          algo::MappedHexReason::InvalidCellCounts);
        VOLUME_TYPE_CHECK(out.destinationUnchanged);
        VOLUME_TYPE_CHECK(mesh::meshHash(destination) == before);
    }
    {
        algo::MappedHexBlockRequest request;
        request.geometry = block(1, 1, 1);
        request.geometry.axisU = {0, 1, 0};
        request.geometry.axisV = {1, 0, 0};
        const auto out = algo::generateMappedHexBlock(request, destination);
        VOLUME_TYPE_CHECK(!out.published());
        VOLUME_TYPE_CHECK(out.reason ==
                          algo::MappedHexReason::NonPositiveFrame);
        VOLUME_TYPE_CHECK(mesh::meshHash(destination) == before);
    }
    {
        algo::MappedHexBlockRequest request;
        request.geometry = block(2, 2, 2);
        request.geometry.maxCells = 7;
        const auto out = algo::generateMappedHexBlock(request, destination);
        VOLUME_TYPE_CHECK(!out.published());
        VOLUME_TYPE_CHECK(out.reason ==
                          algo::MappedHexReason::BudgetExceeded);
        VOLUME_TYPE_CHECK(mesh::meshHash(destination) == before);
    }
    {
        algo::MappedHexBlockRequest request;
        request.geometry = block(1, 1, 1);
        request.geometry.edgeClassifications[3] = {2, 99};
        const auto out = algo::generateMappedHexBlock(request, destination);
        VOLUME_TYPE_CHECK(!out.published());
        VOLUME_TYPE_CHECK(out.reason ==
                          algo::MappedHexReason::InvalidClassification);
        VOLUME_TYPE_CHECK(mesh::meshHash(destination) == before);
    }
    {
        algo::MappedHexBlockRequest request;
        request.geometry = block(1, 1, 1);
        request.geometry.axisW.z =
            (std::numeric_limits<double>::quiet_NaN)();
        const auto out = algo::generateMappedHexBlock(request, destination);
        VOLUME_TYPE_CHECK(!out.published());
        VOLUME_TYPE_CHECK(out.reason ==
                          algo::MappedHexReason::NonfiniteFrame);
        VOLUME_TYPE_CHECK(mesh::meshHash(destination) == before);
    }

    mesh::Mesh inverted;
    for (const Vec3 p : std::array<Vec3, 8>{
             Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{1, 1, 0}, Vec3{0, 1, 0},
             Vec3{0, 0, 1}, Vec3{1, 0, 1}, Vec3{1, 1, 1}, Vec3{0, 1, 1}})
        (void)inverted.addNode(p, {3, 1});
    (void)inverted.addCell(mesh::CellType::Hex8,
                           {0, 3, 2, 1, 4, 7, 6, 5}, {3, 1});
    const auto audit = algo::auditHex8Mesh(inverted, false);
    VOLUME_TYPE_CHECK(!audit.valid);
    VOLUME_TYPE_CHECK(audit.nonPositiveJacobians != 0);
}

void mappedHexDeterminism() {
    VOLUME_TYPE_TEST("hex.v5 direct indexing is bitwise deterministic");
    algo::MappedHexBlockRequest request;
    request.geometry = block(7, 5, 3);
    request.geometry.axisU = {2.0, 0.1, 0.0};
    request.geometry.axisV = {0.0, 3.0, 0.2};
    request.geometry.axisW = {0.1, 0.0, 4.0};
    mesh::Mesh a, b;
    const auto oa = algo::generateMappedHexBlock(request, a);
    const auto ob = algo::generateMappedHexBlock(request, b);
    VOLUME_TYPE_CHECK(oa.published() && ob.published());
    VOLUME_TYPE_CHECK(mesh::meshHash(a) == mesh::meshHash(b));
    VOLUME_TYPE_CHECK(oa.metrics.constructionWorkUnits ==
                      ob.metrics.constructionWorkUnits);
    VOLUME_TYPE_CHECK(oa.metrics.audit.totalVolume ==
                      ob.metrics.audit.totalVolume);
}

void structuredPolyUnitAndSharedFace() {
    VOLUME_TYPE_TEST("poly.v5 face-based block has exact owner/neighbour senses");
    auto geometry = block(2, 1, 1);
    mesh::PolyhedralMesh destination;
    const auto outcome = algo::generateStructuredPolyhedralBlock(
        geometry, destination);
    VOLUME_TYPE_CHECK(outcome.published());
    VOLUME_TYPE_CHECK(destination.nodeCount() == 12);
    VOLUME_TYPE_CHECK(destination.faceCount() == 11);
    VOLUME_TYPE_CHECK(destination.cellCount() == 2);
    VOLUME_TYPE_CHECK(outcome.metrics.audit.valid);
    VOLUME_TYPE_CHECK(outcome.metrics.audit.boundaryFaces == 10);
    VOLUME_TYPE_CHECK(outcome.metrics.audit.interiorFaces == 1);
    VOLUME_TYPE_CHECK(std::abs(outcome.metrics.audit.totalVolume - 1.0) <
                      1.0e-12);
    VOLUME_TYPE_CHECK(sameClass(destination.nodeClass(0),
                                geometry.vertexClassifications[0]));
    VOLUME_TYPE_CHECK(sameClass(destination.nodeClass(1),
                                geometry.edgeClassifications[0]));
    std::uint64_t shared = 0;
    for (Index face = 0; face < destination.faceCount(); ++face)
        shared += static_cast<std::uint64_t>(
            destination.face(face).neighbour != paacs::kInvalidIndex);
    VOLUME_TYPE_CHECK(shared == 1);
    const auto hash = mesh::polyhedralStrictHash(destination);
    mesh::PolyhedralMesh again;
    const auto repeated = algo::generateStructuredPolyhedralBlock(
        geometry, again);
    VOLUME_TYPE_CHECK(repeated.published());
    VOLUME_TYPE_CHECK(hash == mesh::polyhedralStrictHash(again));
}

void structuredPolySkewedAndConversion() {
    VOLUME_TYPE_TEST("poly.v5 direct and Hex8 conversion agree geometrically");
    auto geometry = block(2, 3, 2);
    geometry.axisU = {2.0, 0.0, 0.0};
    geometry.axisV = {0.25, 3.0, 0.0};
    geometry.axisW = {0.1, 0.2, 2.0};
    mesh::PolyhedralMesh direct;
    const auto directOut = algo::generateStructuredPolyhedralBlock(
        geometry, direct, algo::StructuredPolyValidationMode::FullAudit);
    VOLUME_TYPE_CHECK(directOut.published());

    algo::MappedHexBlockRequest request;
    request.geometry = geometry;
    mesh::Mesh hex;
    const auto hexOut = algo::generateMappedHexBlock(request, hex);
    VOLUME_TYPE_CHECK(hexOut.published());
    mesh::PolyhedralMesh converted;
    const auto convertedOut = algo::convertHex8ToPolyhedral(hex, converted);
    VOLUME_TYPE_CHECK(convertedOut.published());
    VOLUME_TYPE_CHECK(converted.nodeCount() == direct.nodeCount());
    VOLUME_TYPE_CHECK(converted.faceCount() == direct.faceCount());
    VOLUME_TYPE_CHECK(converted.cellCount() == direct.cellCount());
    VOLUME_TYPE_CHECK(convertedOut.metrics.audit.valid);
    VOLUME_TYPE_CHECK(std::abs(convertedOut.metrics.audit.totalVolume -
                               directOut.metrics.audit.totalVolume) < 1.0e-10);
    mesh::PolyhedralMesh convertedAgain;
    const auto convertedAgainOut = algo::convertHex8ToPolyhedral(
        hex, convertedAgain);
    VOLUME_TYPE_CHECK(convertedAgainOut.published());
    VOLUME_TYPE_CHECK(mesh::polyhedralStrictHash(converted) ==
                      mesh::polyhedralStrictHash(convertedAgain));
}


void polyhedralMixedFaceArities() {
    VOLUME_TYPE_TEST("poly.v5 mixed triangle/quad pyramid audits as one cell");
    mesh::PolyhedralMesh value;
    for (const Vec3 p : std::array<Vec3, 5>{
             Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{1, 1, 0},
             Vec3{0, 1, 0}, Vec3{0.5, 0.5, 1}})
        (void)value.addNode(p, {0, static_cast<int>(value.nodeCount() + 1)});
    const std::array<Index, 4> base{0, 3, 2, 1};
    const std::array<Index, 3> side0{0, 1, 4};
    const std::array<Index, 3> side1{1, 2, 4};
    const std::array<Index, 3> side2{2, 3, 4};
    const std::array<Index, 3> side3{3, 0, 4};
    const std::array<Index, 5> faces{
        value.addFace(base, 0, paacs::kInvalidIndex, {2, 1}),
        value.addFace(side0, 0, paacs::kInvalidIndex, {2, 2}),
        value.addFace(side1, 0, paacs::kInvalidIndex, {2, 3}),
        value.addFace(side2, 0, paacs::kInvalidIndex, {2, 4}),
        value.addFace(side3, 0, paacs::kInvalidIndex, {2, 5})};
    const std::array<mesh::OrientedPolyFace, 5> refs{{
        {faces[0], 1}, {faces[1], 1}, {faces[2], 1},
        {faces[3], 1}, {faces[4], 1}}};
    (void)value.addCell(refs, {3, 90});
    const auto audit = mesh::auditPolyhedralMesh(value);
    VOLUME_TYPE_CHECK(audit.valid);
    VOLUME_TYPE_CHECK(audit.cells == 1);
    VOLUME_TYPE_CHECK(audit.faces == 5);
    VOLUME_TYPE_CHECK(std::abs(audit.totalVolume - (1.0 / 3.0)) < 1.0e-12);
}

void polyhedralContractsAndAudit() {
    VOLUME_TYPE_TEST("poly.v5 rejects invalid generation and ownership records");
    mesh::PolyhedralMesh destination;
    (void)destination.addNode({-7, -11, -13}, {0, 8});
    const auto before = mesh::polyhedralStrictHash(destination);
    auto geometry = block(1, 1, 1);
    geometry.maxFaces = 5;
    const auto out = algo::generateStructuredPolyhedralBlock(
        geometry, destination);
    VOLUME_TYPE_CHECK(!out.published());
    VOLUME_TYPE_CHECK(out.reason ==
                      algo::StructuredPolyReason::BudgetExceeded);
    VOLUME_TYPE_CHECK(out.destinationUnchanged);
    VOLUME_TYPE_CHECK(mesh::polyhedralStrictHash(destination) == before);

    mesh::PolyhedralMesh broken;
    for (const Vec3 p : std::array<Vec3, 4>{
             Vec3{0, 0, 0}, Vec3{1, 0, 0},
             Vec3{0, 1, 0}, Vec3{0, 0, 1}})
        (void)broken.addNode(p, {3, 1});
    const std::array<Index, 3> f0{1, 2, 3};
    const std::array<Index, 3> f1{0, 3, 2};
    const std::array<Index, 3> f2{0, 1, 3};
    const std::array<Index, 3> f3{0, 2, 1};
    const Index a = broken.addFace(f0, 1, paacs::kInvalidIndex, {2, 1});
    const Index b = broken.addFace(f1, 0, paacs::kInvalidIndex, {2, 2});
    const Index c = broken.addFace(f2, 0, paacs::kInvalidIndex, {2, 3});
    const Index d = broken.addFace(f3, 0, paacs::kInvalidIndex, {2, 4});
    const std::array<mesh::OrientedPolyFace, 4> refs{{
        {a, 1}, {b, 1}, {c, 1}, {d, 1}}};
    (void)broken.addCell(refs, {3, 1});
    const auto audit = mesh::auditPolyhedralMesh(broken);
    VOLUME_TYPE_CHECK(!audit.valid);
    VOLUME_TYPE_CHECK(audit.invalidOwnerNeighbour != 0 ||
                      audit.ownerReferenceMismatch != 0);

    mesh::Mesh mixed;
    for (const Vec3 p : std::array<Vec3, 4>{
             Vec3{0, 0, 0}, Vec3{1, 0, 0},
             Vec3{0, 1, 0}, Vec3{0, 0, 1}})
        (void)mixed.addNode(p, {3, 1});
    (void)mixed.addCell(mesh::CellType::Tet4, {0, 1, 2, 3}, {3, 1});
    mesh::PolyhedralMesh unchanged;
    (void)unchanged.addNode({5, 6, 7}, {0, 3});
    const auto unchangedHash = mesh::polyhedralStrictHash(unchanged);
    const auto converted = algo::convertHex8ToPolyhedral(mixed, unchanged);
    VOLUME_TYPE_CHECK(!converted.published());
    VOLUME_TYPE_CHECK(converted.reason ==
                      algo::StructuredPolyReason::UnsupportedCellType);
    VOLUME_TYPE_CHECK(mesh::polyhedralStrictHash(unchanged) == unchangedHash);
}

}  // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* volumeTypeTests(std::size_t& count) {
    static const Entry tests[]{
        {"hex.v5.unit", mappedHexUnitCube},
        {"hex.v5.skewed", mappedHexSkewedBlock},
        {"hex.v5.contracts", mappedHexContractsAndAudit},
        {"hex.v5.determinism", mappedHexDeterminism},
        {"poly.v5.unit", structuredPolyUnitAndSharedFace},
        {"poly.v5.convert", structuredPolySkewedAndConversion},
        {"poly.v5.mixed_faces", polyhedralMixedFaceArities},
        {"poly.v5.contracts", polyhedralContractsAndAudit}};
    count = std::size(tests);
    return tests;
}
}  // namespace paacs_tests
