// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#ifndef PAACS_V09_CORE_ONLY
#include "paacs/sdk/CAPI.h"
#endif
#include "paacs/sdk/Version.hpp"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using paacs::Index;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

paacs::algo::MappedHexBlockRequest makeRequest(
        std::uint32_t n, int threads) {
    paacs::algo::MappedHexBlockRequest request;
    request.geometry.origin = {1.25, -2.0, 0.5};
    request.geometry.axisU = {2.0, 0.0, 0.0};
    request.geometry.axisV = {0.2, 3.0, 0.0};
    request.geometry.axisW = {0.1, 0.15, 4.0};
    request.geometry.cellsU = n;
    request.geometry.cellsV = n;
    request.geometry.cellsW = n;
    request.geometry.volumeClassification = {3, 901};
    for (std::size_t i = 0;
         i < request.geometry.boundaryClassifications.size(); ++i) {
        request.geometry.boundaryClassifications[i] = {
            2, 910 + static_cast<int>(i)};
    }
    request.executionThreads = threads;
    request.validation =
        paacs::algo::MappedHexValidationMode::CertifiedAffine;
    return request;
}


bool exactlyEqual(const paacs::mesh::Mesh& first,
                  const paacs::mesh::Mesh& second) {
    if (first.nodeCount() != second.nodeCount() ||
        first.cellCount() != second.cellCount() ||
        first.connectivityCount() != second.connectivityCount()) {
        return false;
    }
    for (Index i = 0; i < first.nodeCount(); ++i) {
        const auto& a = first.node(i);
        const auto& b = second.node(i);
        if (a.x != b.x || a.y != b.y || a.z != b.z ||
            first.nodeClass(i) != second.nodeClass(i)) {
            return false;
        }
    }
    for (Index i = 0; i < first.cellCount(); ++i) {
        const auto& a = first.cell(i);
        const auto& b = second.cell(i);
        if (a.type != b.type || a.first != b.first ||
            a.onEntity != b.onEntity) {
            return false;
        }
        const int count = paacs::mesh::nodesPerCell(a.type);
        const Index* ac = first.cellNodes(i);
        const Index* bc = second.cellNodes(i);
        for (int local = 0; local < count; ++local) {
            if (ac[local] != bc[local]) return false;
        }
    }
    return true;
}

void testWorkspacePrimingAndIdentity() {
    const auto request = makeRequest(24, 4);
    paacs::mesh::Mesh baseline;
    const auto baselineOutcome =
        paacs::algo::generateMappedHexBlock(request, baseline);
    require(baselineOutcome.published(), "V8-compatible one-shot failed");

    paacs::algo::MappedHexWorkspace workspace;
    paacs::mesh::Mesh generated;
    const auto first = paacs::algo::generateMappedHexBlock(
        request, workspace, generated);
    require(first.published(), "workspace first call failed");
    require(!first.metrics.reusedConstructionStorage,
            "first call unexpectedly reported reused storage");
    const auto second = paacs::algo::generateMappedHexBlock(
        request, workspace, generated);
    require(second.published(), "workspace second call failed");
    require(!second.metrics.reusedConstructionStorage,
            "second call unexpectedly reported reused storage");
    const auto third = paacs::algo::generateMappedHexBlock(
        request, workspace, generated);
    require(third.published(), "workspace third call failed");
    require(third.metrics.reusedConstructionStorage,
            "third call did not reuse the double buffer");
    require(exactlyEqual(baseline, generated),
            "workspace path changed the explicit mesh");
    // Hash fields explicitly: Cell alignment padding is not mesh state.
    require(paacs::mesh::meshHash(baseline) == paacs::mesh::meshHash(generated),
            "workspace path changed the canonical mesh hash");
    require(workspace.retainedLogicalBytes() > 0U,
            "workspace did not retain the spare output buffer");
    workspace.release();
    require(workspace.retainedLogicalBytes() == 0U,
            "workspace release did not return spare memory");
    const auto afterRelease = paacs::algo::generateMappedHexBlock(
        request, workspace, generated);
    require(afterRelease.published() &&
                !afterRelease.metrics.reusedConstructionStorage,
            "generation after release did not re-prime storage");
}

void testFailureAtomicityAndResize() {
    auto request = makeRequest(18, 3);
    paacs::algo::MappedHexWorkspace workspace;
    paacs::mesh::Mesh destination;
    for (int i = 0; i < 3; ++i) {
        require(paacs::algo::generateMappedHexBlock(
                    request, workspace, destination).published(),
                "workspace priming failed");
    }
    const auto beforeHash = paacs::mesh::meshHash(destination);
    const auto beforeNodes = destination.nodeCount();
    const auto beforeCells = destination.cellCount();
    auto invalid = request;
    invalid.geometry.cellsU = 0;
    const auto rejected = paacs::algo::generateMappedHexBlock(
        invalid, workspace, destination);
    require(!rejected.published(), "invalid request was published");
    require(rejected.destinationUnchanged,
            "invalid request did not report unchanged destination");
    require(destination.nodeCount() == beforeNodes &&
                destination.cellCount() == beforeCells &&
                paacs::mesh::meshHash(destination) == beforeHash,
            "invalid request changed the destination");

    request.geometry.cellsU = 20;
    const auto resized1 = paacs::algo::generateMappedHexBlock(
        request, workspace, destination);
    const auto resized2 = paacs::algo::generateMappedHexBlock(
        request, workspace, destination);
    const auto resized3 = paacs::algo::generateMappedHexBlock(
        request, workspace, destination);
    require(resized1.published() && resized2.published() &&
                resized3.published(),
            "resized requests failed");
    require(!resized1.metrics.reusedConstructionStorage &&
                !resized2.metrics.reusedConstructionStorage &&
                resized3.metrics.reusedConstructionStorage,
            "resized double buffer did not re-prime deterministically");
}

void testThreadDeterminismAndFullAudit() {
    auto oneThread = makeRequest(16, 1);
    auto fourThreads = oneThread;
    fourThreads.executionThreads = 4;
    paacs::algo::MappedHexWorkspace serialWorkspace;
    paacs::algo::MappedHexWorkspace parallelWorkspace;
    paacs::mesh::Mesh serial;
    paacs::mesh::Mesh parallel;
    for (int i = 0; i < 3; ++i) {
        require(paacs::algo::generateMappedHexBlock(
                    oneThread, serialWorkspace, serial).published(),
                "serial workspace generation failed");
        require(paacs::algo::generateMappedHexBlock(
                    fourThreads, parallelWorkspace, parallel).published(),
                "parallel workspace generation failed");
    }
    require(exactlyEqual(serial, parallel),
            "workspace result changed with thread count");

    auto auditedRequest = makeRequest(8, 4);
    auditedRequest.validation =
        paacs::algo::MappedHexValidationMode::FullTopologyAudit;
    paacs::algo::MappedHexWorkspace auditedWorkspace;
    paacs::mesh::Mesh audited;
    for (int i = 0; i < 3; ++i) {
        const auto outcome = paacs::algo::generateMappedHexBlock(
            auditedRequest, auditedWorkspace, audited);
        require(outcome.published(), "full-audit workspace call failed");
        require(outcome.metrics.audit.complete &&
                    outcome.metrics.audit.valid,
                "full-audit workspace result was invalid");
    }
}

void testMovedFromWorkspaceRefuses() {
    const auto request = makeRequest(4, 1);
    paacs::algo::MappedHexWorkspace original;
    paacs::algo::MappedHexWorkspace moved(std::move(original));
    paacs::mesh::Mesh result;
    require(paacs::algo::generateMappedHexBlock(
                request, moved, result).published(),
            "moved-to workspace failed");
    const auto refused = paacs::algo::generateMappedHexBlock(
        request, original, result);
    require(!refused.published() && refused.destinationUnchanged,
            "moved-from workspace was not safely refused");
}

#ifndef PAACS_V09_CORE_ONLY
void testCapiSessionReuseCompatibility() {
    paacs_v8_session* session = paacs_v8_session_create();
    require(session != nullptr, "C API session allocation failed");
    paacs_v8_hex_request request{};
    request.axis_u = {1.0, 0.0, 0.0};
    request.axis_v = {0.0, 1.0, 0.0};
    request.axis_w = {0.0, 0.0, 1.0};
    request.cells_u = 10;
    request.cells_v = 9;
    request.cells_w = 8;
    request.volume_tag = 1;
    request.boundary_tag_base = 100;
    request.threads = 4;
    std::uint64_t retainedHash = 0;
    for (int i = 0; i < 4; ++i) {
        require(paacs_v8_generate_mapped_hex(session, &request) ==
                    PAACS_V8_OK,
                "C API repeated mapped-Hex generation failed");
        const auto hash = paacs_v8_session_mesh_hash(session);
        if (i == 0) retainedHash = hash;
        require(hash == retainedHash,
                "C API repeated generation changed the mesh");
    }
    require(std::string(paacs_v8_api_version()) ==
                paacs::sdk::kApiVersion,
            "C API version does not match the SDK version");
    paacs_v8_session_destroy(session);
}
#endif

}  // namespace

int main() {
    try {
        testWorkspacePrimingAndIdentity();
        testFailureAtomicityAndResize();
        testThreadDeterminismAndFullAudit();
        testMovedFromWorkspaceRefuses();
#ifndef PAACS_V09_CORE_ONLY
        testCapiSessionReuseCompatibility();
#endif
        std::cout << "PAACS V0.9 reusable mapped-Hex probe: PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "PAACS V0.9 reusable mapped-Hex probe: FAIL: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
