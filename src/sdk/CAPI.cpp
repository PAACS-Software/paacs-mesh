// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sdk/CAPI.h"

#include "paacs/algo/ModelGeometryProviders.hpp"
#include "paacs/algo/PreOcctMeshingCore.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sizing/SizeField.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/sdk/MeshingSession.hpp"

#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <string>

struct paacs_v8_session {
    paacs::mesh::Mesh mesh;
    // Retained across calls so UI/Python sessions automatically use V0.9's
    // allocation-free mapped-Hex double buffer after priming.
    paacs::algo::MappedHexWorkspace mappedHexWorkspace;
    std::string error;
};

namespace {
paacs_v8_status failed(paacs_v8_session* session,
                       const char* message,
                       paacs_v8_status status = PAACS_V8_FAILED) noexcept {
    if (session) {
        try {
            session->error = message;
        } catch (...) {
            session->error.clear();
        }
    }
    return status;
}

paacs_v8_classification classification(
        const paacs::mesh::Classification value) noexcept {
    return {value.dim, value.tag};
}

paacs::algo::PreOcctTetOutput tetOutput(
        paacs_v010_tet_output value) noexcept {
    using Output = paacs::algo::PreOcctTetOutput;
    switch (value) {
        case PAACS_V010_TET4_CONFORMING:
            return Output::BoundaryConformingTet4;
        case PAACS_V010_TET4_REFINED:
            return Output::RefinedTet4;
        case PAACS_V010_TET4_QUALITY:
            return Output::QualityTet4;
        case PAACS_V010_TET10_CERTIFIED:
            return Output::CertifiedTet10;
    }
    return Output::QualityTet4;
}

bool finitePoint(const paacs_v8_vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

bool validTetOptions(const paacs_v010_tet_options& value) noexcept {
    return std::isfinite(value.target_size) && value.target_size > 0.0 &&
           value.maximum_nodes > 0 && value.maximum_tets > 0 &&
           value.threads >= 0 &&
           value.output >= PAACS_V010_TET4_CONFORMING &&
           value.output <= PAACS_V010_TET10_CERTIFIED;
}

paacs_v8_status generatePreparedTet(
        paacs_v8_session* session,
        std::unique_ptr<paacs::geo::Model> model,
        const paacs_v010_tet_options& options) noexcept {
    if (!session || !model || !validTetOptions(options))
        return PAACS_V8_INVALID_ARGUMENT;
    try {
        auto ownership =
            paacs::algo::buildModelBoundaryOwnershipProvider(*model);
        if (!ownership.built)
            return failed(session, "model ownership preparation failed");
        paacs::algo::ModelTet10PlacementProvider placement(*model);
        paacs::sizing::MeshSpec spec;
        spec.field = std::make_shared<paacs::sizing::ConstantField>(
            options.target_size);
        spec.threads = options.threads;

        paacs::algo::PreOcctResourceLimits limits;
        limits.maximumNodes = options.maximum_nodes;
        limits.maximumTets = options.maximum_tets;
        limits.maximumWorkUnits = 2000000;
        limits.maximumReportRows = 200000;
        auto request = paacs::algo::makePreOcctTetRequest(
            tetOutput(options.output), ownership.provider.stamp(),
            placement.stamp(), limits);
        request.pipeline.source.source.refinement.maxInsertions =
            options.maximum_insertions;
        if (options.maximum_insertions == 0) {
            request.pipeline.source.source.refinement.maxEdgeToTarget = 1.0e9;
            request.pipeline.source.source.refinement.maxRadiusEdge = 1.0e9;
        }
        const auto result = paacs::algo::meshPreparedTetModel(
            *model, spec, request, session->mesh, ownership.provider,
            &placement);
        if (!result.published())
            return failed(session, "prepared Tet production route stopped");
        session->error.clear();
        return PAACS_V8_OK;
    } catch (const std::bad_alloc&) {
        return failed(session, "allocation failed",
                      PAACS_V8_ALLOCATION_FAILED);
    } catch (...) {
        return failed(session, "unexpected Tet generation failure");
    }
}
}  // namespace

extern "C" {

const char* paacs_v8_api_version(void) {
    return paacs::sdk::kApiVersion;
}

paacs_v8_session* paacs_v8_session_create(void) {
    try {
        return new paacs_v8_session;
    } catch (...) {
        return nullptr;
    }
}

void paacs_v8_session_destroy(paacs_v8_session* session) {
    delete session;
}

paacs_v8_status paacs_v8_generate_mapped_hex(
        paacs_v8_session* session,
        const paacs_v8_hex_request* request) {
    if (!session || !request) return PAACS_V8_INVALID_ARGUMENT;
    if (request->boundary_tag_base < 0 ||
        request->boundary_tag_base > std::numeric_limits<int32_t>::max() - 5 ||
        request->volume_tag < 0 || request->threads < 0 ||
        !finitePoint(request->origin) || !finitePoint(request->axis_u) ||
        !finitePoint(request->axis_v) || !finitePoint(request->axis_w))
        return failed(session, "invalid mapped Hex8 request", PAACS_V8_INVALID_ARGUMENT);
    try {
        paacs::algo::MappedHexBlockRequest internal;
        internal.geometry.origin = {
            request->origin.x, request->origin.y, request->origin.z};
        internal.geometry.axisU = {
            request->axis_u.x, request->axis_u.y, request->axis_u.z};
        internal.geometry.axisV = {
            request->axis_v.x, request->axis_v.y, request->axis_v.z};
        internal.geometry.axisW = {
            request->axis_w.x, request->axis_w.y, request->axis_w.z};
        internal.geometry.cellsU = request->cells_u;
        internal.geometry.cellsV = request->cells_v;
        internal.geometry.cellsW = request->cells_w;
        internal.geometry.volumeClassification = {3, request->volume_tag};
        internal.executionThreads = request->threads;
        for (std::size_t i = 0;
             i < internal.geometry.boundaryClassifications.size(); ++i) {
            internal.geometry.boundaryClassifications[i] = {
                2, request->boundary_tag_base + static_cast<int32_t>(i)};
        }
        const auto outcome = paacs::algo::generateMappedHexBlock(
            internal, session->mappedHexWorkspace, session->mesh);
        if (!outcome.published())
            return failed(session, "mapped Hex8 generation failed");
        session->error.clear();
        return PAACS_V8_OK;
    } catch (const std::bad_alloc&) {
        return failed(session, "allocation failed",
                      PAACS_V8_ALLOCATION_FAILED);
    } catch (...) {
        return failed(session, "unexpected generation failure");
    }
}

paacs_v8_status paacs_v010_generate_box_tet(
        paacs_v8_session* session,
        const paacs_v010_box_request* request) {
    if (!session || !request) return PAACS_V8_INVALID_ARGUMENT;
    if (!validTetOptions(request->options) ||
        !finitePoint(request->lower) || !finitePoint(request->upper) ||
        !(request->upper.x > request->lower.x) ||
        !(request->upper.y > request->lower.y) ||
        !(request->upper.z > request->lower.z))
        return PAACS_V8_INVALID_ARGUMENT;
    try {
        auto kernel = paacs::geo::makeBuiltinKernel();
        return generatePreparedTet(
            session,
            kernel->makeBox(
                {request->lower.x, request->lower.y, request->lower.z},
                {request->upper.x, request->upper.y, request->upper.z}),
            request->options);
    } catch (...) {
        return failed(session, "box model construction failed");
    }
}

paacs_v8_status paacs_v010_generate_cylinder_tet(
        paacs_v8_session* session,
        const paacs_v010_cylinder_request* request) {
    if (!session || !request) return PAACS_V8_INVALID_ARGUMENT;
    if (!validTetOptions(request->options) ||
        !finitePoint(request->base_center) ||
        !std::isfinite(request->radius) || request->radius <= 0.0 ||
        !std::isfinite(request->height) || request->height <= 0.0)
        return PAACS_V8_INVALID_ARGUMENT;
    try {
        auto kernel = paacs::geo::makeBuiltinKernel();
        return generatePreparedTet(
            session,
            kernel->makeCylinder(
                {request->base_center.x, request->base_center.y,
                 request->base_center.z},
                request->radius, request->height),
            request->options);
    } catch (...) {
        return failed(session, "cylinder model construction failed");
    }
}

paacs_v8_status paacs_v8_session_read_archive(
        paacs_v8_session* session,
        const char* path) {
    if (!session || !path) return PAACS_V8_INVALID_ARGUMENT;
    paacs::mesh::Mesh loaded;
    const auto outcome = paacs::sdk::readMeshArchive(path, loaded);
    if (!outcome.completed)
        return failed(session, "archive read failed");
    loaded.swap(session->mesh);
    session->error.clear();
    return PAACS_V8_OK;
}

paacs_v8_status paacs_v8_session_write_archive(
        const paacs_v8_session* session,
        const char* path) {
    if (!session || !path) return PAACS_V8_INVALID_ARGUMENT;
    const auto outcome = paacs::sdk::writeMeshArchive(session->mesh, path);
    return outcome.completed ? PAACS_V8_OK : PAACS_V8_FAILED;
}

uint64_t paacs_v8_session_node_count(const paacs_v8_session* session) {
    return session
        ? static_cast<uint64_t>(session->mesh.nodeCount())
        : 0U;
}

uint64_t paacs_v8_session_cell_count(const paacs_v8_session* session) {
    return session
        ? static_cast<uint64_t>(session->mesh.cellCount())
        : 0U;
}

uint64_t paacs_v8_session_mesh_hash(const paacs_v8_session* session) {
    return session ? paacs::mesh::meshHash(session->mesh) : 0U;
}

paacs_v8_status paacs_v8_session_node_at(
        const paacs_v8_session* session,
        uint64_t index,
        paacs_v8_vec3* position,
        paacs_v8_classification* nodeClassification) {
    if (!session || !position || !nodeClassification)
        return PAACS_V8_INVALID_ARGUMENT;
    if (index >= static_cast<uint64_t>(session->mesh.nodeCount()))
        return PAACS_V8_OUT_OF_RANGE;
    const auto& point = session->mesh.node(static_cast<paacs::Index>(index));
    *position = {point.x, point.y, point.z};
    *nodeClassification = classification(
        session->mesh.nodeClass(static_cast<paacs::Index>(index)));
    return PAACS_V8_OK;
}

paacs_v8_status paacs_v8_session_cell_at(
        const paacs_v8_session* session,
        uint64_t index,
        paacs_v8_cell_view* view) {
    if (!session || !view) return PAACS_V8_INVALID_ARGUMENT;
    if (index >= static_cast<uint64_t>(session->mesh.cellCount()))
        return PAACS_V8_OUT_OF_RANGE;
    const auto cellIndex = static_cast<paacs::Index>(index);
    const auto& cell = session->mesh.cell(cellIndex);
    const int count = paacs::mesh::nodesPerCell(cell.type);
    if (count <= 0 || count > 10)
        return PAACS_V8_FAILED;
    view->type = static_cast<paacs_v8_cell_type>(cell.type);
    view->node_count = static_cast<uint32_t>(count);
    view->classification = classification(cell.onEntity);
    const paacs::Index* nodes = session->mesh.cellNodes(cellIndex);
    for (int i = 0; i < count; ++i)
        view->nodes[i] = static_cast<uint64_t>(nodes[i]);
    for (int i = count; i < 10; ++i) view->nodes[i] = 0U;
    return PAACS_V8_OK;
}

const char* paacs_v8_session_last_error(
        const paacs_v8_session* session) {
    return session ? session->error.c_str() : "invalid session";
}

}  // extern "C"
