// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/CurvedPrismLayerMesher.hpp"
#include "paacs/algo/HybridMeshingExecutor.hpp"
#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/MeshingPlan.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/algo/PreparedSweepMesher.hpp"
#include "paacs/algo/SurfaceSizingPlanner.hpp"
#include "paacs/algo/TransitionTemplateCatalog.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/sdk/CAPI.h"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/sdk/Version.hpp"
#include "paacs/sizing/CompiledSpatialField.hpp"
#include "paacs/sizing/SizeField.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using paacs::Index;
using paacs::Vec3;
namespace algo = paacs::algo;
namespace mesh = paacs::mesh;
namespace sizing = paacs::sizing;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

bool close(double a, double b, double tol = 1.0e-10) {
    return std::abs(a - b) <= tol * std::max({1.0, std::abs(a), std::abs(b)});
}

void addPositiveTet(mesh::Mesh& value, std::array<Index, 4> nodes,
                    mesh::Classification classification) {
    const auto sign = paacs::predicates::orient3d(
        value.node(nodes[0]), value.node(nodes[1]),
        value.node(nodes[2]), value.node(nodes[3]));
    if (sign == paacs::predicates::Sign::Negative)
        std::swap(nodes[2], nodes[3]);
    require(paacs::predicates::orient3d(
                value.node(nodes[0]), value.node(nodes[1]),
                value.node(nodes[2]), value.node(nodes[3])) ==
                paacs::predicates::Sign::Positive,
            "tetrahedron is not positively oriented");
    (void)value.addCell(mesh::CellType::Tet4,
                        {nodes[0], nodes[1], nodes[2], nodes[3]},
                        classification);
}

mesh::Mesh mixedSweepSource() {
    mesh::Mesh source;
    const mesh::Classification surface{2, 20};
    // Unit Quad4 at x=[0,1].
    for (const Vec3 p : std::array<Vec3, 4>{{
             {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}}})
        (void)source.addNode(p, surface);
    (void)source.addCell(mesh::CellType::Quad4, {0,1,2,3}, surface);
    // Disconnected right triangle of area 0.5 at x=[2,3].
    for (const Vec3 p : std::array<Vec3, 3>{{
             {2,0,0}, {3,0,0}, {2,1,0}}})
        (void)source.addNode(p, surface);
    (void)source.addCell(mesh::CellType::Tri3, {4,5,6}, surface);
    return source;
}

mesh::Mesh quadSurface() {
    mesh::Mesh surface;
    for (const Vec3 p : std::array<Vec3,4>{{
             {0,0,0},{1,0,0},{1,1,0},{0,1,0}}})
        (void)surface.addNode(p, {2, 30});
    (void)surface.addCell(mesh::CellType::Quad4, {0,1,2,3}, {2,30});
    return surface;
}

mesh::Mesh triangulatedSquare(double z, double halfWidth,
                              int faceTag) {
    mesh::Mesh surface;
    const mesh::Classification cls{2, faceTag};
    const std::array<Vec3,5> points{{
        {-halfWidth,-halfWidth,z}, {halfWidth,-halfWidth,z},
        {halfWidth,halfWidth,z}, {-halfWidth,halfWidth,z}, {0,0,z}}};
    for (const auto& p : points) (void)surface.addNode(p, cls);
    (void)surface.addCell(mesh::CellType::Tri3, {0,1,4}, cls);
    (void)surface.addCell(mesh::CellType::Tri3, {1,2,4}, cls);
    (void)surface.addCell(mesh::CellType::Tri3, {2,3,4}, cls);
    (void)surface.addCell(mesh::CellType::Tri3, {3,0,4}, cls);
    return surface;
}

void testSpatialField() {
    sizing::CompiledSpatialFieldRequest request;
    request.lower = {-1,-1,-1};
    request.upper = {1,1,1};
    request.samplesX = request.samplesY = request.samplesZ = 28;
    request.baseSize = 0.4;
    request.minimumSize = 0.05;
    request.maximumSize = 0.4;
    request.maximumAdjacentGrowth = 1.2;
    request.maximumGrowthSweeps = 32;
    request.controls.push_back({sizing::SpatialControlKind::Ball,
                                {0,0,0}, {}, 0.35, 0.0, 0.08, 0.0});
    request.controls.push_back({sizing::SpatialControlKind::DistanceRamp,
                                {0.65,0,0}, {}, 0.0, 0.45, 0.12, 0.35});
    request.executionThreads = 1;
    auto serial = sizing::compileSpatialField(request);
    require(serial.compiled, "serial spatial field compilation failed");
    request.executionThreads = 4;
    auto parallel = sizing::compileSpatialField(request);
    require(parallel.compiled, "parallel spatial field compilation failed");
    require(serial.field.deterministicHash() ==
                parallel.field.deterministicHash(),
            "spatial field changed with thread count");
    require(serial.field.eval({0,0,0}) < serial.field.eval({-1,-1,-1}),
            "spatial refinement was not retained");
    require(serial.metrics.measuredMaximumAdjacentGrowth <= 1.20001,
            "spatial-field growth bound was exceeded");
}

algo::MeshingPlan makeHybridPlan() {
    algo::MeshingPlanRequest request;
    request.compiledControlHash = 0x1234;
    request.compiledSpatialFieldHash = 0x5678;
    request.bodies.push_back({1, sizing::VolumeEngineKind::Tet4,
                              algo::BodyMeshingCapability::GeneralTet4,
                              0, 4, 1, true});
    request.bodies.push_back({2, sizing::VolumeEngineKind::Tet4,
                              algo::BodyMeshingCapability::GeneralTet4,
                              0, 4, 1, true});
    request.portals.push_back({10,1,2,algo::PortalTopology::Tri3,
                               3,1,0xAABBCCDD,true});
    const auto result = algo::compileMeshingPlan(request);
    require(result.compiled, "hybrid Tet plan compilation failed");
    require(result.metrics.tetBodies == 2, "unexpected Tet body count");
    require(result.plan.portals().front().strategy ==
                algo::PortalStrategy::SharedTri3,
            "unexpected Tri3 portal strategy");

    // Separately prove automatic mapped selection and Pyramid5 bridging.
    algo::MeshingPlanRequest autoRequest;
    autoRequest.bodies.push_back({3, sizing::VolumeEngineKind::Auto,
        algo::BodyMeshingCapability::GeneralTet4 |
        algo::BodyMeshingCapability::MappedBlock |
        algo::BodyMeshingCapability::PyramidTransition,
        0x99, 8, 1, true});
    autoRequest.bodies.push_back({4, sizing::VolumeEngineKind::Tet4,
        algo::BodyMeshingCapability::GeneralTet4 |
        algo::BodyMeshingCapability::PyramidTransition,
        0, 4, 1, true});
    autoRequest.portals.push_back({11,3,4,algo::PortalTopology::Quad4,
                                   4,1,0x77,true});
    const auto automatic = algo::compileMeshingPlan(autoRequest);
    require(automatic.compiled, "automatic plan compilation failed");
    require(automatic.plan.body(3)->selectedEngine ==
                sizing::VolumeEngineKind::MappedHex8,
            "Auto did not select mapped Hex8");
    require(automatic.plan.portals().front().strategy ==
                algo::PortalStrategy::Pyramid5Bridge,
            "Hex/Tet interface did not select a Pyramid5 bridge");
    return result.plan;
}

void testPreparedSweep() {
    const mesh::Mesh source = mixedSweepSource();
    algo::PreparedMixedSweepRequest request;
    request.layers = 3;
    request.cellToCellGrowth = 1.15;
    request.maximumAllowedGrowth = 1.2;
    request.direction = {0,0,1};
    request.totalThickness = 2.0;
    request.executionThreads = 4;
    mesh::Mesh output;
    const auto result = algo::generatePreparedMixedSweep(source, request, output);
    require(result.published(), "mixed prepared sweep failed");
    require(result.metrics.hexes == 3 && result.metrics.wedges == 3,
            "mixed sweep cell counts are wrong");
    require(close(result.metrics.totalVolume, 3.0, 1.0e-9),
            "mixed sweep volume is wrong");
    require(result.metrics.audit.valid, "mixed sweep audit failed");

    request.fullAudit = false;
    mesh::Mesh certified;
    const auto certifiedResult =
        algo::generatePreparedMixedSweep(source, request, certified);
    require(certifiedResult.published(),
            "certified mixed prepared sweep failed");
    require(certifiedResult.metrics.audit.valid,
            "certified mixed sweep audit failed");
    require(mesh::meshHash(output) == mesh::meshHash(certified),
            "certified sweep changed explicit mesh output");
    require(close(certifiedResult.metrics.totalVolume,
                  result.metrics.totalVolume, 1.0e-12),
            "certified sweep volume differs from full audit");
}

mesh::Mesh testTransitionsAndQuality() {
    const mesh::Mesh base = quadSurface();
    algo::PyramidTransitionRequest request;
    request.apex = {0.5,0.5,1.0};
    request.emitBoundaryFaces = true;
    mesh::Mesh pyramid;
    const auto generated = algo::generatePyramidTransitionFan(base, request, pyramid);
    require(generated.published(), "Pyramid5 transition generation failed");
    require(generated.metrics.outputPyramid5 == 1,
            "Pyramid5 transition count is wrong");
    require(close(generated.metrics.totalVolume, 1.0 / 3.0),
            "Pyramid5 transition volume is wrong");
    const auto quality = algo::analyzeMixedMeshQuality(pyramid);
    require(quality.valid && quality.pyramid5 == 1,
            "mixed-cell quality did not validate Pyramid5");
    mesh::PolyhedralMesh polyhedral;
    const auto polyResult = algo::convertSupportedVolumeMeshToPolyhedral(
        pyramid, polyhedral);
    require(polyResult.published() && polyhedral.cellCount() == 1,
            "Pyramid5-to-polyhedral conversion failed");
    require(mesh::auditPolyhedralMesh(polyhedral).valid,
            "Pyramid5 polyhedral audit failed");

    algo::MixedToTet4Request tetRequest;
    mesh::Mesh tet;
    const auto converted = algo::convertSupportedMixedMeshToTet4(
        pyramid, tetRequest, tet);
    require(converted.published(), "Pyramid5-to-Tet4 conversion failed");
    require(converted.metrics.outputTet4 == 2,
            "Pyramid5 should convert to two Tet4s");
    require(close(converted.metrics.totalVolume, 1.0 / 3.0),
            "converted Tet4 volume is wrong");
    return pyramid;
}

void testSurfaceSizing() {
    mesh::Mesh surface;
    for (const auto p : std::array<Vec3,3>{{{0,0,0},{2,0,0},{0,2,0}}})
        (void)surface.addNode(p, {2,40});
    (void)surface.addCell(mesh::CellType::Tri3, {0,1,2}, {2,40});
    sizing::ConstantField field(0.5);
    algo::SurfaceSizingPlanRequest request;
    request.field = &field;
    request.maximumEdgeToTargetRatio = 1.25;
    request.executionThreads = 4;
    request.project = [](int, const Vec3& in, Vec3& out) {
        out = in; return true;
    };
    const auto result = algo::planSurfaceSizingPass(surface, request);
    require(result.planned, "surface sizing pass failed");
    require(result.placements.size() == 3,
            "surface sizing did not identify all three long edges");
    require(std::is_sorted(result.placements.begin(), result.placements.end(),
        [](const auto& a, const auto& b) {
            if (a.violationRatio != b.violationRatio)
                return a.violationRatio > b.violationRatio;
            if (a.faceTag != b.faceTag) return a.faceTag < b.faceTag;
            if (a.firstNode != b.firstNode) return a.firstNode < b.firstNode;
            return a.secondNode < b.secondNode;
        }), "surface sizing placements are not canonical");
}

void testCurvedLayers() {
    const mesh::Mesh source = triangulatedSquare(0.0, 0.5, 50);
    const mesh::Mesh obstacle = triangulatedSquare(0.2, 0.8, 51);
    algo::CurvedPrismLayerRequest request;
    request.layers = 4;
    request.firstLayerThickness = 0.08;
    request.cellToCellGrowth = 1.15;
    request.maximumAllowedGrowth = 1.2;
    request.maximumTotalThickness = 0.5;
    request.minimumTotalThickness = 1.0e-4;
    request.collisionSurface = &obstacle;
    request.collisionSafetyFactor = 0.75;
    request.executionThreads = 4;
    mesh::Mesh layers;
    const auto result = algo::generateCurvedPrismLayers(source, request, layers);
    require(result.published, "collision-aware prism layer generation failed");
    require(result.metrics.clippedNodes ==
                static_cast<std::uint64_t>(source.nodeCount()),
            "not every layer node was clipped by the ceiling");
    require(result.metrics.maximumThickness <= 0.150001,
            "collision-safe layer exceeded the expected ceiling");
    require(result.metrics.sweep.audit.valid,
            "curved layer sweep failed its volume audit");
}

void testMixedOptimizer() {
    mesh::Mesh value;
    const std::array<Vec3,5> points{{
        {0,0,0},{1,0,0},{0,1,0},{0,0,1},{0.55,0.12,0.12}}};
    for (std::size_t i=0;i<points.size();++i)
        (void)value.addNode(points[i], i == 4 ? mesh::Classification{3,60}
                                               : mesh::Classification{2,61});
    addPositiveTet(value, {0,1,2,4}, {3,60});
    addPositiveTet(value, {0,1,4,3}, {3,60});
    addPositiveTet(value, {1,2,4,3}, {3,60});
    addPositiveTet(value, {2,0,4,3}, {3,60});
    const auto before = algo::analyzeMixedMeshQuality(value);
    require(before.valid, "optimizer fixture is invalid before smoothing");
    algo::MixedMeshOptimizationRequest request;
    request.maximumIterations = 8;
    request.relaxation = 0.5;
    const auto result = algo::optimizeMixedMeshQuality(value, request);
    require(result.completed, "mixed-cell smart smoothing failed");
    require(result.after.minimumQuality + 1.0e-14 >= before.minimumQuality,
            "mixed-cell optimization reduced the worst quality");
}

mesh::Mesh sharedTetBody(bool upper, int solidTag) {
    mesh::Mesh body;
    const std::array<Vec3,4> points{{
        {0,0,0},{1,0,0},{0,1,0}, upper ? Vec3{0,0,1} : Vec3{0,0,-1}}};
    for (std::size_t i=0;i<points.size();++i)
        (void)body.addNode(points[i], i < 3 ? mesh::Classification{2,10}
                                            : mesh::Classification{3,solidTag});
    addPositiveTet(body, {0,1,2,3}, {3,solidTag});
    return body;
}

void testHybridExecutor(const algo::MeshingPlan& plan) {
    const mesh::Mesh first = sharedTetBody(true, 1);
    const mesh::Mesh second = sharedTetBody(false, 2);
    algo::HybridExecutionRequest request;
    request.plan = &plan;
    request.bodies.push_back({1,algo::HybridBodySourceKind::PreparedFixedMesh,
                              &first});
    request.bodies.push_back({2,algo::HybridBodySourceKind::PreparedFixedMesh,
                              &second});
    for (Index n=0;n<3;++n)
        request.portalNodeBindings.push_back({10,1,n,2,n});
    std::uint64_t progressCalls = 0;
    request.progress = [&](algo::HybridExecutionStage,
                           std::uint64_t, std::uint64_t) {
        ++progressCalls; return true;
    };
    algo::HybridMesh output;
    const auto result = algo::executeHybridMeshingPlan(request, output);
    require(result.published, "hybrid model executor failed");
    require(output.fixed.nodeCount() == 5 && output.fixed.cellCount() == 2,
            "hybrid executor did not merge the conformal interface");
    require(result.metrics.mergedPortalNodes == 3,
            "hybrid portal merge count is wrong");
    require(output.portals.size() == 1 && output.portals[0].sharedNodeCount == 3,
            "hybrid portal receipt is incomplete");
    require(progressCalls >= 4, "hybrid progress callback was not exercised");
}

void testArchiveAndCapi(const mesh::Mesh& source) {
    std::vector<std::uint8_t> bytes;
    const auto written = paacs::sdk::serializeMeshArchive(source, bytes);
    require(written.completed && !bytes.empty(), "mesh archive serialization failed");
    mesh::Mesh restored;
    const auto read = paacs::sdk::deserializeMeshArchive(bytes, restored);
    require(read.completed, "mesh archive deserialization failed");
    require(mesh::meshHash(source) == mesh::meshHash(restored),
            "mesh archive round-trip changed the strict mesh hash");
    require(std::string(paacs_v8_api_version()) == paacs::sdk::kApiVersion,
            "C API and C++ API versions disagree");
    paacs_v8_session* session = paacs_v8_session_create();
    require(session != nullptr, "C API session allocation failed");
    paacs_v8_hex_request request{};
    request.axis_u = {1,0,0}; request.axis_v = {0,1,0};
    request.axis_w = {0,0,1}; request.cells_u = 2;
    request.cells_v = 2; request.cells_w = 2;
    request.volume_tag = 70; request.boundary_tag_base = 100;
    require(paacs_v8_generate_mapped_hex(session, &request) == PAACS_V8_OK,
            std::string("C API mapped Hex8 failed: ") +
                paacs_v8_session_last_error(session));
    require(paacs_v8_session_node_count(session) == 27,
            "C API Hex8 node count is wrong");
    require(paacs_v8_session_cell_count(session) == 32,
            "C API Hex8 cell count including boundary faces is wrong");
    const auto path = std::filesystem::temp_directory_path() /
                      "paacs_v8_capi_smoke.paacsmesh";
    require(paacs_v8_session_write_archive(session, path.string().c_str()) ==
                PAACS_V8_OK,
            "C API archive write failed");
    require(std::filesystem::file_size(path) > 0,
            "C API archive file is empty");
    require(paacs_v8_session_mesh_hash(session) != 0,
            "C API mesh hash was not exposed");
    paacs_v8_vec3 position{};
    paacs_v8_classification nodeClass{};
    require(paacs_v8_session_node_at(session, 0, &position, &nodeClass) ==
                PAACS_V8_OK,
            "C API node view failed");
    require(nodeClass.dimension == 0 && nodeClass.tag >= 0,
            "C API node classification is wrong");
    paacs_v8_cell_view cell{};
    require(paacs_v8_session_cell_at(session, 0, &cell) == PAACS_V8_OK,
            "C API cell view failed");
    require(cell.type == PAACS_V8_CELL_HEX8 && cell.node_count == 8,
            "C API Hex8 cell view is wrong");
    require(paacs_v8_session_cell_at(session, 999999, &cell) ==
                PAACS_V8_OUT_OF_RANGE,
            "C API did not report an out-of-range cell");

    paacs_v8_session* loadedSession = paacs_v8_session_create();
    require(loadedSession != nullptr, "second C API session allocation failed");
    require(paacs_v8_session_read_archive(
                loadedSession, path.string().c_str()) == PAACS_V8_OK,
            "C API archive read failed");
    require(paacs_v8_session_mesh_hash(loadedSession) ==
                paacs_v8_session_mesh_hash(session),
            "C API archive round-trip changed the mesh hash");
    paacs_v8_session_destroy(loadedSession);
    std::filesystem::remove(path);
    paacs_v8_session_destroy(session);
}

void testRequestLocalParallel() {
    std::vector<std::uint64_t> one(20000), four(20000);
    paacs::core::parallelForWithThreads(0, one.size(), 1,
        [&](std::size_t i) { one[i] = i * i + 17; });
    paacs::core::parallelForWithThreads(0, four.size(), 4,
        [&](std::size_t i) { four[i] = i * i + 17; });
    require(one == four, "request-local parallel execution changed output");
}

}  // namespace

int main() {
    try {
        testSpatialField();
        const auto plan = makeHybridPlan();
        testPreparedSweep();
        const auto pyramid = testTransitionsAndQuality();
        testSurfaceSizing();
        testCurvedLayers();
        testMixedOptimizer();
        testHybridExecutor(plan);
        testArchiveAndCapi(pyramid);
        testRequestLocalParallel();
        std::cout << "V8_RELEASE_PROBE_PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "V8_RELEASE_PROBE_FAIL: " << error.what() << '\n';
        return 1;
    }
}
