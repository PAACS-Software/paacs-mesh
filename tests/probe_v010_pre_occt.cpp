// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/CanonicalCurveMesher.hpp"
#include "paacs/algo/ModelGeometryProviders.hpp"
#include "paacs/algo/PreOcctMeshingCore.hpp"
#include "paacs/algo/PreparedSweepMesher.hpp"
#include "paacs/algo/TransitionTemplateCatalog.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sizing/MetricField.hpp"
#include "paacs/sdk/CAPI.h"
#include "geometry_model_views.hpp"
#include "shared_interface_fixture.hpp"

#include <cmath>
#include <cstdio>
#include <memory>

namespace {
int failures = 0;
#define CHECK(expr) do { if (!(expr)) { \
    std::printf("FAIL line %d: %s\n", __LINE__, #expr); ++failures; } \
} while (false)

void metricAndCurves(const paacs::geo::Model& model,
                     const paacs::sizing::MeshSpec& spec) {
    using namespace paacs;
    const auto curves = algo::buildCanonicalCurves(model, spec);
    CHECK(curves.built());
    CHECK(!curves.result.empty());
    for (const auto& record : curves.result.curves()) {
        const auto forward = curves.result.orientedNodes(record.edgeTag, false);
        const auto reverse = curves.result.orientedNodes(record.edgeTag, true);
        CHECK(!forward.empty());
        CHECK(forward.size() == reverse.size());
        if (!forward.empty()) {
            CHECK(forward.front() == reverse.back());
            CHECK(forward.back() == reverse.front());
        }
    }

    sizing::DirectionalMetricControl directional;
    directional.axis = {1.0, 0.0, 0.0};
    directional.alongSize = 0.1;
    directional.transverseSize = 0.25;
    const auto metric = sizing::directionalMetric(directional);
    CHECK(metric.positiveDefinite());
    const auto bounds = sizing::metricEigenBounds(metric);
    CHECK(bounds.valid);
    CHECK(bounds.maximum > bounds.minimum);
    CHECK(std::abs(metric.metricLength({0.1, 0, 0}) - 1.0) < 1e-10);

    auto field = std::make_shared<sizing::ConstantMetricField>(metric);
    sizing::MetricSizeFieldAdapter scalar(field, 0.05, 1.0);
    CHECK(scalar.valid());
    CHECK(std::abs(scalar.eval({0,0,0}) - 0.1) < 1e-10);
}

void productionTet(const paacs::geo::Model& model,
                   const paacs::sizing::MeshSpec& spec) {
    using namespace paacs;
    auto ownership = algo::buildModelBoundaryOwnershipProvider(model);
    CHECK(ownership.built);
    if (!ownership.built) return;
    algo::ModelTet10PlacementProvider placement(model);

    algo::PreOcctResourceLimits limits;
    limits.maximumNodes = 100000;
    limits.maximumTets = 500000;
    limits.maximumWorkUnits = 10000000;
    limits.maximumReportRows = 100000;

    // Report-only quality path is the robust default production route.
    auto qualityRequest = algo::makePreOcctTetRequest(
        algo::PreOcctTetOutput::QualityTet4,
        ownership.provider.stamp(), placement.stamp(), limits);
    qualityRequest.pipeline.source.source.refinement.maxEdgeToTarget = 1.0e9;
    qualityRequest.pipeline.source.source.refinement.maxRadiusEdge = 1.0e9;
    qualityRequest.pipeline.source.source.refinement.maxInsertions = 0;
    mesh::Mesh qualityMesh;
    const auto quality = algo::meshPreparedTetModel(
        model, spec, qualityRequest, qualityMesh, ownership.provider,
        &placement);
    CHECK(quality.published());
    CHECK(qualityMesh.countCells(mesh::CellType::Tet4) > 0);
    CHECK(quality.audit.passed());
    CHECK(quality.tetQuality.nInverted == 0);

    auto tet10Request = algo::makePreOcctTetRequest(
        algo::PreOcctTetOutput::CertifiedTet10,
        ownership.provider.stamp(), placement.stamp(), limits);
    tet10Request.pipeline.source.source.refinement.maxEdgeToTarget = 1.0e9;
    tet10Request.pipeline.source.source.refinement.maxRadiusEdge = 1.0e9;
    tet10Request.pipeline.source.source.refinement.maxInsertions = 0;
    mesh::Mesh tet10Mesh;
    const auto tet10 = algo::meshPreparedTetModel(
        model, spec, tet10Request, tet10Mesh, ownership.provider,
        &placement);
    CHECK(tet10.published());
    CHECK(tet10Mesh.countCells(mesh::CellType::Tet10) > 0);
    CHECK(tet10.audit.nInvalidTet10 == 0);
    CHECK(tet10.tetQuality.nTet10Invalid == 0);

    // Missing placement provider must be a structured, atomic rejection.
    mesh::Mesh untouched;
    untouched.addNode({7,8,9}, {0,77});
    const auto before = mesh::meshHash(untouched);
    const auto rejected = algo::meshPreparedTetModel(
        model, spec, tet10Request, untouched, ownership.provider, nullptr);
    CHECK(!rejected.published());
    CHECK(rejected.reason ==
          algo::PreOcctMeshingReason::MissingTet10PlacementProvider);
    CHECK(mesh::meshHash(untouched) == before);
}

void componentAccounting() {
    using namespace paacs;
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(10.0);
    auto run = [&](const geo::Model& model, mesh::Mesh& destination) {
        auto ownership = algo::buildModelBoundaryOwnershipProvider(model);
        CHECK(ownership.built);
        auto request = algo::makePreOcctTetRequest(
            algo::PreOcctTetOutput::BoundaryConformingTet4,
            ownership.provider.stamp());
        return algo::meshPreparedTetModel(
            model, spec, request, destination, ownership.provider);
    };
    const auto componentMismatch = [](const mesh::AuditResult& audit) {
        return std::any_of(audit.report.items.begin(), audit.report.items.end(),
            [](const diag::Diagnostic& row) {
                return row.code == diag::Code::ComponentCountMismatch;
            });
    };

    geometry_test::SharedFaceTetrahedraModel interfaced;
    mesh::Mesh unchanged;
    unchanged.addNode({7,8,9}, {0,77});
    const auto before = mesh::meshHash(unchanged);
    const auto shared = run(interfaced, unchanged);
    // The lower engine supports and proves this exact interface. The wrapper
    // still has an exterior-only final boundary audit, an independent limit.
    CHECK(shared.conforming && shared.conforming->published());
    CHECK(!shared.published());
    CHECK(shared.reason == algo::PreOcctMeshingReason::FinalAuditFailed);
    CHECK(shared.phase == algo::PreOcctMeshingPhase::FinalAudit);
    CHECK(shared.destinationUnchanged && mesh::meshHash(unchanged) == before);
    CHECK(shared.audit.connectedComponents == 1);
    CHECK(!componentMismatch(shared.audit));
    if (shared.conforming && shared.conforming->retainedResult) {
        const auto& exact = shared.conforming->retainedResult->boundaryAudit;
        CHECK(exact.passed && exact.prescribedInterface > 0);
        CHECK(exact.matchedInterface == exact.prescribedInterface);
        CHECK(shared.audit.boundaryOrphans ==
              static_cast<Index>(exact.prescribedInterface));
    }
    for (const auto& row : shared.audit.report.items)
        CHECK(row.code == diag::Code::BoundaryOrphanFacet);

    auto kernel = geo::makeBuiltinKernel();
    auto first = kernel->makeBox({-3,-1,-1}, {-1,1,1});
    auto second = kernel->makeBox({1,-1,-1}, {3,1,1});
    geometry_test::CombinedModel disconnected;
    disconnected.add(*first, 0);
    disconnected.add(*second, 1000);
    mesh::Mesh separate;
    const auto two = run(disconnected, separate);
    CHECK(two.published() && two.audit.passed());
    CHECK(two.audit.connectedComponents == 2);
    CHECK(two.audit.boundaryOrphans == 0 && two.audit.boundaryUnexpected == 0);
    CHECK(!componentMismatch(two.audit));
    std::printf("  shared solid count=%zu measured components=%d; "
                "disconnected solid count=%zu measured components=%d\n",
        interfaced.solids().size(), shared.audit.connectedComponents,
        disconnected.solids().size(), two.audit.connectedComponents);
}

void cApiFlow() {
    paacs_v8_session* session = paacs_v8_session_create();
    CHECK(session != nullptr);
    if (!session) return;
    paacs_v010_box_request request{};
    request.lower = {0,0,0};
    request.upper = {1,1,1};
    request.options.target_size = 0.45;
    request.options.output = PAACS_V010_TET4_QUALITY;
    request.options.maximum_insertions = 0;
    request.options.maximum_nodes = 100000;
    request.options.maximum_tets = 500000;
    request.options.threads = 1;
    CHECK(paacs_v010_generate_box_tet(session, &request) == PAACS_V8_OK);
    CHECK(paacs_v8_session_node_count(session) > 0);
    CHECK(paacs_v8_session_cell_count(session) > 0);
    paacs_v8_session_destroy(session);
}

void layersAndTransitions() {
    using namespace paacs;
    mesh::Mesh surface;
    const Index a = surface.addNode({0,0,0}, {2,1});
    const Index b = surface.addNode({1,0,0}, {2,1});
    const Index c = surface.addNode({0,1,0}, {2,1});
    surface.addCell(mesh::CellType::Tri3, {a,b,c}, {2,1});

    algo::PreparedMixedSweepRequest sweep;
    sweep.direction = {0,0,1};
    sweep.layers = 2;
    sweep.totalThickness = 0.22;
    sweep.cellToCellGrowth = 1.2;
    sweep.volumeClassification = {3,11};
    sweep.targetClassification = {2,2};
    sweep.sideClassification = {2,3};
    mesh::Mesh layered;
    const auto swept = algo::generatePreparedMixedSweep(surface, sweep, layered);
    CHECK(swept.published());
    CHECK(layered.countCells(mesh::CellType::Prism6) == 2);

    algo::MixedToTet4Request toTet;
    mesh::Mesh tet;
    const auto converted = algo::convertSupportedMixedMeshToTet4(
        layered, toTet, tet);
    CHECK(converted.published());
    CHECK(tet.countCells(mesh::CellType::Tet4) > 0);
    CHECK(mesh::analyzeTets(tet).nInverted == 0);
}
}

int main() {
    using namespace paacs;
    auto kernel = geo::makeBuiltinKernel();
    auto model = kernel->makeBox({0,0,0}, {1,1,1});
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.45);
    spec.threads = 1;
    metricAndCurves(*model, spec);
    productionTet(*model, spec);
    componentAccounting();
    layersAndTransitions();
    cApiFlow();
    if (failures == 0) {
        std::puts("V0.10 pre-OCCT integration probe passed");
        return 0;
    }
    std::printf("V0.10 probe failures: %d\n", failures);
    return 1;
}
