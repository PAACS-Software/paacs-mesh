// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "paacs/bnd/TetFacetRecovery.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "test_harness.hpp"
#include "workspace_digest.hpp"

namespace {
using namespace paacs;

bnd::TetFacetCell positiveCell(
    const std::vector<bnd::TetFacetVertex>& vertices,
    std::array<bnd::TetFacetVertexId, 4> ids, std::uint64_t sequence,
    mesh::Classification classification = {3, 9}) {
    const auto orientation = predicates::orient3d(
        vertices[ids[0]].position, vertices[ids[1]].position,
        vertices[ids[2]].position, vertices[ids[3]].position);
    if (orientation == predicates::Sign::Negative)
        std::swap(ids[0], ids[1]);
    return {ids, sequence, classification};
}

struct BipyramidFixture {
    std::vector<bnd::TetFacetVertex> vertices{
        {{0, 0, 0}, 10}, {{4, 0, 0}, 20}, {{0, 4, 0}, 30},
        {{1, 1, 3}, 40}, {{1, 1, -3}, 50}};
    std::vector<bnd::TetFacetCell> cells;
    std::vector<bnd::TetFacetEdgeKey> protectedEdges;
    std::vector<bnd::TetFacetProtectedFace> protectedFaces;
    bnd::TetFacetTriangle target{{0, 1, 2}};

    BipyramidFixture() {
        cells.push_back(positiveCell(vertices, {3, 4, 0, 1}, 100));
        cells.push_back(positiveCell(vertices, {3, 4, 1, 2}, 110));
        cells.push_back(positiveCell(vertices, {3, 4, 2, 0}, 120));
        protectedEdges = {
            bnd::canonicalTetFacetEdge(0, 1),
            bnd::canonicalTetFacetEdge(1, 2),
            bnd::canonicalTetFacetEdge(2, 0)};
    }

    bnd::TetFacetComplexView view() const {
        return {vertices, cells, target, protectedEdges, protectedFaces,
                {3, 9}};
    }
};

struct TetWorkspaceFixture {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::vector<topo::NodeHandle> nodes;
    std::vector<topo::CellHandle> cells;

    template <class Source>
    explicit TetWorkspaceFixture(const Source& source) {
        topo::NodePool nodePool;
        for (const auto& vertex : source.vertices)
            nodes.push_back(nodePool.add(vertex.position, {3, 9}));
        topo::CellPool cellPool;
        for (const auto& cell : source.cells) {
            std::array<topo::NodeHandle, 4> connectivity{};
            for (std::size_t i = 0; i < 4; ++i)
                connectivity[i] = nodes[cell.vertices[i]];
            cells.push_back(cellPool.add(topo::ElementType::Tet4,
                                         connectivity, cell.classification,
                                         nodePool));
        }
        topo::EdgePool edgePool;
        for (const auto cell : cells)
            (void)edgePool.insertCellEdges(nodePool, cellPool, cell);
        for (const auto& edge : source.protectedEdges) {
            const auto handle = edgePool.find(
                nodePool, nodes[edge.lo], nodes[edge.hi]);
            if (!handle) throw std::logic_error("missing fixture target edge");
            edgePool.setRequiredConstraint(*handle, true);
            edgePool.setClassification(*handle, {1, 41});
        }
        topo::FacePool facePool;
        for (const auto cell : cells)
            (void)facePool.insertCellFaces(nodePool, cellPool, edgePool,
                                           cell);
        for (const auto& protectedFace : source.protectedFaces) {
            std::array<topo::NodeHandle, 3> corners{};
            for (std::size_t i = 0; i < 3; ++i)
                corners[i] = nodes[protectedFace.key.vertices[i]];
            const auto handle = facePool.find(nodePool, corners);
            if (!handle)
                throw std::logic_error("missing fixture protected face");
            facePool.setClassification(*handle,
                                       protectedFace.classification);
        }
        workspace = std::make_unique<topo::TopologyWorkspace>(
            std::move(nodePool), std::move(cellPool), std::move(edgePool),
            std::move(facePool));
    }
};

topo::TetCavitySpec lowerCavitySpec(
    const BipyramidFixture& source, const TetWorkspaceFixture& destination,
    const bnd::TetFacetCavityPlan& plan,
    const bnd::TetFacetCandidate& candidate) {
    topo::TetCavitySpec spec;
    spec.stamp = topo::WorkspaceStamp::of(*destination.workspace);
    for (const auto index : plan.selectedCells)
        spec.removed.push_back(destination.cells[index]);
    for (const auto& sourceCell : candidate.cells) {
        topo::TetCavityCell cell;
        for (std::size_t i = 0; i < 4; ++i)
            cell.nodes[i] = destination.nodes[sourceCell.vertices[i]];
        cell.cls = sourceCell.classification;
        spec.added.push_back(cell);
    }
    for (std::size_t i = 0; i < 3; ++i)
        spec.target[i] = destination.nodes[source.target.vertices[i]];
    spec.targetClassification = {2, 77};
    spec.region = {3, 9};
    return spec;
}

// Scope-valid valence-four edge star. X is strictly inside the original
// D-E-A-B Tet; its stellar subdivision increases the D-E star valence while
// retaining all three target edges A-B, B-C and C-A as real protected edges.
struct ScopeValidHighValenceFixture {
    std::vector<bnd::TetFacetVertex> vertices{
        {{0, 0, 0}, 10},   {{4, 0, 0}, 20}, {{0, 4, 0}, 30},
        {{1, 1, 3}, 40},   {{1, 1, -3}, 50},
        {{1.5, 0.5, 0.6}, 60}};
    std::vector<bnd::TetFacetCell> cells;
    std::vector<bnd::TetFacetEdgeKey> protectedEdges{
        bnd::canonicalTetFacetEdge(0, 1),
        bnd::canonicalTetFacetEdge(1, 2),
        bnd::canonicalTetFacetEdge(2, 0)};
    std::vector<bnd::TetFacetProtectedFace> protectedFaces;
    bnd::TetFacetTriangle target{{0, 1, 2}};

    ScopeValidHighValenceFixture() {
        cells.push_back(positiveCell(vertices, {5, 4, 0, 1}, 100));
        cells.push_back(positiveCell(vertices, {3, 5, 0, 1}, 110));
        cells.push_back(positiveCell(vertices, {3, 4, 5, 1}, 120));
        cells.push_back(positiveCell(vertices, {3, 4, 0, 5}, 130));
        cells.push_back(positiveCell(vertices, {3, 4, 1, 2}, 140));
        cells.push_back(positiveCell(vertices, {3, 4, 2, 0}, 150));
    }

    bnd::TetFacetComplexView view() const {
        return {vertices, cells, target, protectedEdges, protectedFaces,
                {3, 9}};
    }
};

struct HighValenceFixture {
    std::vector<bnd::TetFacetVertex> vertices{
        {{0, 0, -2}, 10}, {{0, 0, 2}, 20},
        {{4, 0, 0}, 30},   {{2, 3, 0}, 40},
        {{-2, 3, 0}, 50},  {{-4, 0, 0}, 60},
        {{-2, -3, 0}, 70}, {{2, -3, 0}, 80}};
    std::vector<bnd::TetFacetCell> cells;
    std::vector<bnd::TetFacetEdgeKey> protectedEdges;
    std::vector<bnd::TetFacetProtectedFace> protectedFaces;
    bnd::TetFacetTriangle target{{2, 4, 6}};

    HighValenceFixture() {
        for (std::uint32_t i = 0; i < 6; ++i) {
            const auto a = static_cast<bnd::TetFacetVertexId>(2 + i);
            const auto b = static_cast<bnd::TetFacetVertexId>(2 + (i + 1) % 6);
            cells.push_back(positiveCell(vertices, {0, 1, a, b},
                                         100 + i * 10));
        }
    }

    bnd::TetFacetComplexView view() const {
        return {vertices, cells, target, protectedEdges, protectedFaces,
                {3, 9}};
    }
};

bool samePlainComplex(const HighValenceFixture& a,
                      const HighValenceFixture& b) {
    if (a.vertices.size() != b.vertices.size() ||
        a.cells.size() != b.cells.size() ||
        a.protectedEdges != b.protectedEdges ||
        a.target.vertices != b.target.vertices)
        return false;
    for (std::size_t i = 0; i < a.vertices.size(); ++i) {
        const auto& av = a.vertices[i];
        const auto& bv = b.vertices[i];
        if (av.position.x != bv.position.x ||
            av.position.y != bv.position.y ||
            av.position.z != bv.position.z || av.sequence != bv.sequence)
            return false;
    }
    for (std::size_t i = 0; i < a.cells.size(); ++i) {
        const auto& ac = a.cells[i];
        const auto& bc = b.cells[i];
        if (ac.vertices != bc.vertices || ac.sequence != bc.sequence ||
            ac.classification.dim != bc.classification.dim ||
            ac.classification.tag != bc.classification.tag)
            return false;
    }
    return true;
}

std::vector<std::uint64_t> selectedSequences(
    const bnd::TetFacetComplexView& view,
    const bnd::TetFacetCavityPlan& plan) {
    std::vector<std::uint64_t> out;
    for (const auto ci : plan.selectedCells)
        out.push_back(view.cells[ci].sequence);
    return out;
}

void exactContactTaxonomy() {
    TEST("t2a.geometry.exact_contact_taxonomy");
    BipyramidFixture fixture;
    for (const auto& cell : fixture.cells)
        CHECK(bnd::classifyTetFacetContact(fixture.view(), cell) ==
              bnd::TetFacetContact::ProperInteriorCut);

    auto exact = positiveCell(fixture.vertices, {3, 0, 1, 2}, 200);
    CHECK(bnd::classifyTetFacetContact(fixture.view(), exact) ==
          bnd::TetFacetContact::CoplanarExactFace);

    std::vector<bnd::TetFacetVertex> touchVertices = fixture.vertices;
    touchVertices.push_back({{0, 0, 1}, 60});
    touchVertices.push_back({{2, 0, 1}, 70});
    touchVertices.push_back({{0, 2, 1}, 80});
    auto declared = positiveCell(touchVertices, {0, 5, 6, 7}, 210);
    std::vector<bnd::TetFacetCell> declaredCells{declared};
    const bnd::TetFacetComplexView declaredView{
        touchVertices, declaredCells, fixture.target, fixture.protectedEdges,
        fixture.protectedFaces, {3, 9}};
    CHECK(bnd::classifyTetFacetContact(declaredView, declared) ==
          bnd::TetFacetContact::DeclaredBoundaryTouch);

    touchVertices.push_back({{1, 1, 0}, 90});
    auto undeclared = positiveCell(touchVertices, {8, 5, 6, 7}, 220);
    std::vector<bnd::TetFacetCell> undeclaredCells{undeclared};
    const bnd::TetFacetComplexView undeclaredView{
        touchVertices, undeclaredCells, fixture.target,
        fixture.protectedEdges, fixture.protectedFaces, {3, 9}};
    CHECK(bnd::classifyTetFacetContact(undeclaredView, undeclared) ==
          bnd::TetFacetContact::UndeclaredVertexTouch);

    touchVertices.push_back({{0.5, 0.5, 0}, 100});
    touchVertices.push_back({{2.0, 0.5, 0}, 110});
    touchVertices.push_back({{0.5, 2.0, 0}, 120});
    auto coplanar = positiveCell(touchVertices, {5, 9, 10, 11}, 230);
    std::vector<bnd::TetFacetCell> coplanarCells{coplanar};
    const bnd::TetFacetComplexView coplanarView{
        touchVertices, coplanarCells, fixture.target,
        fixture.protectedEdges, fixture.protectedFaces, {3, 9}};
    CHECK(bnd::classifyTetFacetContact(coplanarView, coplanar) ==
          bnd::TetFacetContact::CoplanarPartialOverlap);
}

void highValenceEdgeStarClosure() {
    TEST("t2a.geometry.high_valence_edge_star_closure");
    ScopeValidHighValenceFixture fixture;
    for (const auto& targetEdge : fixture.protectedEdges) {
        bool found = false;
        for (const auto& cell : fixture.cells)
            for (std::size_t i = 0; i < cell.vertices.size(); ++i)
                for (std::size_t j = i + 1; j < cell.vertices.size(); ++j)
                    found = found ||
                            bnd::canonicalTetFacetEdge(cell.vertices[i],
                                                       cell.vertices[j]) ==
                                targetEdge;
        CHECK(found);
    }
    const auto plan = bnd::planTetFacetCavity(fixture.view());
    CHECK(plan.planned());
    CHECK(plan.reason == bnd::TetFacetGeometryReason::None);
    CHECK(plan.plan.counters.cellsScanned == fixture.cells.size());
    CHECK(plan.plan.counters.crossedEdgeStars >= 1);
    CHECK(plan.plan.selectedCells.size() >= 4);
    std::size_t centralStar = 0;
    for (std::size_t ci = 0; ci < fixture.cells.size(); ++ci) {
        const auto& vertices = fixture.cells[ci].vertices;
        const bool hasD =
            std::find(vertices.begin(), vertices.end(),
                      bnd::TetFacetVertexId{3}) != vertices.end();
        const bool hasE =
            std::find(vertices.begin(), vertices.end(),
                      bnd::TetFacetVertexId{4}) != vertices.end();
        if (!hasD || !hasE) continue;
        ++centralStar;
        CHECK(std::find(plan.plan.selectedCells.begin(),
                        plan.plan.selectedCells.end(), ci) !=
              plan.plan.selectedCells.end());
    }
    CHECK(centralStar == 4);
    CHECK(std::find(fixture.protectedEdges.begin(),
                    fixture.protectedEdges.end(),
                    bnd::canonicalTetFacetEdge(3, 4)) ==
          fixture.protectedEdges.end());
    const auto built =
        bnd::reconstructTetFacetCandidate(fixture.view(), plan.plan);
    CHECK(built.built());
    CHECK(built.reason == bnd::TetFacetGeometryReason::None);
}

void testingOnlySelectiveFixedPointClosure() {
    TEST("t2a.geometry.testing_only_selective_fixed_point");
    HighValenceFixture fixture;
    const auto before = fixture;
    const std::array<std::size_t, 1> seed{0};
    const auto plan =
        bnd::planTetFacetCavityWithSeedsForTesting(fixture.view(), seed);
    CHECK(plan.planned());
    CHECK(plan.reason == bnd::TetFacetGeometryReason::None);
    CHECK(plan.plan.seedCells.size() == 1);
    CHECK(plan.plan.selectedCells.size() == 6);
    CHECK(plan.plan.counters.cellsScanned == 6);
    CHECK(plan.plan.counters.crossedEdgeStars == 1);
    CHECK(plan.plan.counters.closureInsertions == 5);
    CHECK(plan.plan.counters.fixedPointPasses >= 1);
    CHECK(bnd::auditTetFacetBall(fixture.view(), plan.plan.selectedCells).valid);
    CHECK(samePlainComplex(fixture, before));

    fixture.protectedEdges.push_back(
        bnd::canonicalTetFacetEdge(0, 1));
    const auto protectedBefore = fixture;
    const auto refused =
        bnd::planTetFacetCavityWithSeedsForTesting(fixture.view(), seed);
    CHECK(!refused.planned());
    CHECK(refused.status == bnd::TetFacetPlanStatus::Refused);
    CHECK(refused.reason ==
          bnd::TetFacetGeometryReason::CrossedProtectedEdge);
    CHECK(refused.plan.counters.crossedEdgeStars == 1);
    CHECK(samePlainComplex(fixture, protectedBefore));
}

void deterministicClosureUnderChurn() {
    TEST("t2a.geometry.fixed_point_determinism");
    HighValenceFixture reference;
    const std::array<std::size_t, 1> referenceSeed{0};
    const auto a = bnd::planTetFacetCavityWithSeedsForTesting(
        reference.view(), referenceSeed);
    CHECK(a.planned());

    HighValenceFixture churned;
    std::rotate(churned.cells.begin(), churned.cells.begin() + 2,
                churned.cells.end());
    std::size_t churnedSeed = churned.cells.size();
    for (std::size_t i = 0; i < churned.cells.size(); ++i) {
        churned.cells[i].sequence += 10000;
        if (churned.cells[i].vertices == reference.cells[0].vertices)
            churnedSeed = i;
    }
    // Keep the rotated storage order: stable logical order must come from
    // sequence identity, not container position. Sequence gaps model churn.
    const std::array<std::size_t, 1> seed{churnedSeed};
    const auto b = bnd::planTetFacetCavityWithSeedsForTesting(
        churned.view(), seed);
    CHECK(b.planned());
    CHECK(a.plan.selectedCells.size() == b.plan.selectedCells.size());
    CHECK(a.plan.counters.cellsScanned == b.plan.counters.cellsScanned);
    CHECK(a.plan.counters.crossedEdgeStars ==
          b.plan.counters.crossedEdgeStars);
    CHECK(a.plan.counters.closureInsertions ==
          b.plan.counters.closureInsertions);

    auto aSeq = selectedSequences(reference.view(), a.plan);
    auto bSeq = selectedSequences(churned.view(), b.plan);
    for (auto& sequence : bSeq) sequence -= 10000;
    CHECK(aSeq == bSeq);
}

void ballAuditAndCandidateReconstruction() {
    TEST("t2a.geometry.ball_and_candidate");
    BipyramidFixture fixture;
    const auto plan = bnd::planTetFacetCavity(fixture.view());
    CHECK(plan.planned());
    CHECK(plan.plan.seedCells.size() == 3);
    CHECK(plan.plan.selectedCells.size() == 3);
    CHECK(plan.plan.counters.crossedEdgeStars == 1);
    CHECK(bnd::auditTetFacetBall(fixture.view(), plan.plan.selectedCells).valid);

    const auto built = bnd::reconstructTetFacetCandidate(fixture.view(),
                                                         plan.plan);
    CHECK(built.built());
    CHECK(built.reason == bnd::TetFacetGeometryReason::None);
    CHECK(built.candidate.cells.size() == 2);
    CHECK(built.candidate.targetKey ==
          bnd::canonicalTetFacetFace(fixture.target.vertices));
    CHECK(built.candidate.counters.apexAttempts == 2);
    CHECK(built.candidate.counters.tetTetIntersectionTests >= 1);
    bnd::TetFacetGeometryCounters auditCounters{};
    const auto audit = bnd::validateTetFacetCandidate(
        fixture.view(), plan.plan.selectedCells, built.candidate.cells,
        auditCounters);
    CHECK(audit.valid);
    CHECK(audit.reason == bnd::TetFacetGeometryReason::None);

    auto reversed = fixture;
    std::swap(reversed.target.vertices[0], reversed.target.vertices[1]);
    const auto reversedPlan = bnd::planTetFacetCavity(reversed.view());
    const auto reversedBuilt = bnd::reconstructTetFacetCandidate(
        reversed.view(), reversedPlan.plan);
    CHECK(reversedPlan.planned());
    CHECK(reversedBuilt.built());
    if (reversedBuilt.built()) {
        std::vector<std::array<bnd::TetFacetVertexId, 4>> directKeys;
        std::vector<std::array<bnd::TetFacetVertexId, 4>> reverseKeys;
        for (const auto& cell : built.candidate.cells) {
            auto key = cell.vertices;
            std::sort(key.begin(), key.end());
            directKeys.push_back(key);
        }
        for (const auto& cell : reversedBuilt.candidate.cells) {
            auto key = cell.vertices;
            std::sort(key.begin(), key.end());
            reverseKeys.push_back(key);
        }
        std::sort(directKeys.begin(), directKeys.end());
        std::sort(reverseKeys.begin(), reverseKeys.end());
        CHECK(directKeys == reverseKeys);
    }
}

void candidateIntersectionAndProtectionAudits() {
    TEST("t2a.geometry.candidate_intersection_protection");
    BipyramidFixture fixture;
    const auto plan = bnd::planTetFacetCavity(fixture.view());
    const auto built = bnd::reconstructTetFacetCandidate(fixture.view(),
                                                         plan.plan);
    CHECK(plan.planned());
    CHECK(built.built());
    if (!built.built()) return;

    const auto boundaryFace = bnd::canonicalTetFacetFace({3, 0, 1});
    fixture.protectedFaces.push_back({boundaryFace, 700, {2, 70}});
    const auto retainedPlan = bnd::planTetFacetCavity(fixture.view());
    const auto retained = bnd::reconstructTetFacetCandidate(
        fixture.view(), retainedPlan.plan);
    CHECK(retainedPlan.planned());
    CHECK(retained.built());

    auto protectedInternal = fixture;
    protectedInternal.protectedFaces.push_back(
        {bnd::canonicalTetFacetFace({3, 4, 0}), 701, {2, 71}});
    const auto internalPlan =
        bnd::planTetFacetCavity(protectedInternal.view());
    const auto internal = bnd::reconstructTetFacetCandidate(
        protectedInternal.view(), internalPlan.plan);
    CHECK(internalPlan.planned());
    CHECK(!internal.built());
    CHECK(internal.reason == bnd::TetFacetGeometryReason::ProtectedFaceLost);

    auto withSurvivor = fixture;
    withSurvivor.cells.push_back(built.candidate.cells[0]);
    withSurvivor.cells.back().sequence = 900;
    bnd::TetFacetGeometryCounters counters{};
    const auto intersection = bnd::validateTetFacetCandidate(
        withSurvivor.view(), plan.plan.selectedCells, built.candidate.cells,
        counters);
    CHECK(!intersection.valid);
    CHECK(intersection.reason ==
          bnd::TetFacetGeometryReason::CandidateIntersection);
    CHECK(counters.tetTetIntersectionTests > 0);

    std::vector<std::size_t> duplicateSelection = plan.plan.selectedCells;
    duplicateSelection.push_back(duplicateSelection.front());
    const auto duplicate =
        bnd::auditTetFacetBall(fixture.view(), duplicateSelection);
    CHECK(!duplicate.valid);
    CHECK(duplicate.reason == bnd::TetFacetGeometryReason::InvalidInput);
}

void atomicTetCavityReplacement() {
    TEST("t2a.topology.atomic_tet_cavity");
    BipyramidFixture source;
    const auto plan = bnd::planTetFacetCavity(source.view());
    const auto candidate =
        bnd::reconstructTetFacetCandidate(source.view(), plan.plan);
    CHECK(plan.planned());
    CHECK(candidate.built());
    if (!candidate.built()) return;

    TetWorkspaceFixture fixture(source);
    auto& workspace = *fixture.workspace;
    const std::array<topo::NodeHandle, 3> retainedFaceNodes{
        fixture.nodes[3], fixture.nodes[0], fixture.nodes[1]};
    const auto retainedBefore =
        workspace.faces().find(workspace.nodes(), retainedFaceNodes);
    CHECK(retainedBefore.has_value());
    const auto revisionBefore = workspace.topologyRevision();

    const auto spec = lowerCavitySpec(source, fixture, plan.plan,
                                      candidate.candidate);
    auto transaction = workspace.beginTetCavityTransaction();
    const auto& preview = transaction.stage(spec);
    CHECK(preview.cellsToRetire == 3);
    CHECK(preview.cellsToCreate == 2);
    CHECK(preview.facesToCreate >= 1);
    transaction.prepare();
    CHECK(transaction.prepared());
    const auto predicted = transaction.predictedTargetFace();
    CHECK(predicted.isSpecified());
    const auto committed = transaction.commit();
    CHECK(committed.targetFace == predicted);
    CHECK(committed.resultCells.size() == 2);
    CHECK(committed.revision == revisionBefore + 1);
    CHECK(workspace.state() == topo::WorkspaceState::CommittedAdjacencyStale);
    CHECK(topo::workspaceTopologyReady(workspace.nodes(), workspace.cells(),
                                       workspace.edges(), workspace.faces()));
    const auto targetData = workspace.faces().get(committed.targetFace);
    CHECK(targetData.cls.dim == 2);
    CHECK(targetData.cls.tag == 77);
    CHECK(targetData.incident.size() == 2);
    CHECK(targetData.incident[0].sense * targetData.incident[1].sense == -1);
    const auto retainedAfter =
        workspace.faces().find(workspace.nodes(), retainedFaceNodes);
    CHECK(retainedAfter == retainedBefore);
}

template <class Source>
bnd::TetFacetRequest makeRecoveryRequest(
    const Source& source, const TetWorkspaceFixture& fixture,
    std::uint64_t entityId, mesh::Classification classification = {2, 77}) {
    bnd::TetFacetRequest request;
    request.stamp = topo::WorkspaceStamp::of(*fixture.workspace);
    request.entityId = entityId;
    for (std::size_t i = 0; i < 3; ++i)
        request.vertices[i] = fixture.nodes[source.target.vertices[i]];
    request.classification = classification;
    return request;
}

void structuredTetFacetRecoveryRoutes() {
    TEST("t2a.recovery.structured_routes");
    BipyramidFixture absentSource;
    TetWorkspaceFixture absent(absentSource);
    auto state = bnd::TetFacetConstraintState::forWorkspace(*absent.workspace);
    const auto request = makeRecoveryRequest(absentSource, absent, 9001);
    const auto recovered =
        bnd::recoverTetFacet(*absent.workspace, state, request);
    CHECK(recovered.status == bnd::TetFacetRecoveryStatus::Recovered);
    CHECK(recovered.reason == bnd::TetFacetRecoveryReason::None);
    CHECK(recovered.face.isSpecified());
    CHECK(state.stamp.revision == request.stamp.revision + 1);
    CHECK(state.facets.size() == 1);
    CHECK(state.facets.at(9001).recoveredFace == recovered.face);

    auto retryRequest = request;
    retryRequest.stamp = topo::WorkspaceStamp::of(*absent.workspace);
    const auto already =
        bnd::recoverTetFacet(*absent.workspace, state, retryRequest);
    CHECK(already.status ==
          bnd::TetFacetRecoveryStatus::AlreadyPresent);
    CHECK(already.face == recovered.face);
    CHECK(absent.workspace->topologyRevision() == state.stamp.revision);

    BipyramidFixture basis;
    const auto plan = bnd::planTetFacetCavity(basis.view());
    const auto built = bnd::reconstructTetFacetCandidate(basis.view(),
                                                         plan.plan);
    CHECK(built.built());
    if (!built.built()) return;

    BipyramidFixture presentSource = basis;
    presentSource.cells = built.candidate.cells;
    TetWorkspaceFixture present(presentSource);
    auto presentState =
        bnd::TetFacetConstraintState::forWorkspace(*present.workspace);
    const auto markRequest = makeRecoveryRequest(presentSource, present, 9002,
                                                 {2, 78});
    const auto revisionBefore = present.workspace->topologyRevision();
    const auto marked =
        bnd::recoverTetFacet(*present.workspace, presentState, markRequest);
    CHECK(marked.status == bnd::TetFacetRecoveryStatus::Marked);
    CHECK(marked.face.isSpecified());
    CHECK(present.workspace->topologyRevision() == revisionBefore + 1);
    CHECK(present.workspace->faces().get(marked.face).cls.tag == 78);

    BipyramidFixture oneSidedSource = basis;
    oneSidedSource.cells = {built.candidate.cells.front()};
    TetWorkspaceFixture oneSided(oneSidedSource);
    auto oneSidedState =
        bnd::TetFacetConstraintState::forWorkspace(*oneSided.workspace);
    const auto oneSidedRequest =
        makeRecoveryRequest(oneSidedSource, oneSided, 9003, {2, 79});
    const auto oneRevision = oneSided.workspace->topologyRevision();
    const auto refused = bnd::recoverTetFacet(
        *oneSided.workspace, oneSidedState, oneSidedRequest);
    CHECK(refused.status == bnd::TetFacetRecoveryStatus::Refused);
    CHECK(refused.reason ==
          bnd::TetFacetRecoveryReason::OneSidedTargetUnsupported);
    CHECK(oneSided.workspace->topologyRevision() == oneRevision);
    CHECK(oneSidedState.facets.empty());
}

void highValenceRecoveryAndProtectedRefusal() {
    TEST("t2a.recovery.high_valence_protected_pair");
    ScopeValidHighValenceFixture source;
    TetWorkspaceFixture recoverable(source);
    auto recoveredState =
        bnd::TetFacetConstraintState::forWorkspace(*recoverable.workspace);
    const auto recoveredRequest =
        makeRecoveryRequest(source, recoverable, 9100, {2, 88});
    const auto recovered = bnd::recoverTetFacet(
        *recoverable.workspace, recoveredState, recoveredRequest);
    CHECK(recovered.status == bnd::TetFacetRecoveryStatus::Recovered);
    CHECK(recovered.counters.crossedEdgeStars >= 1);
    CHECK(recoverable.workspace->faces().get(recovered.face).incident.size() ==
          2);

    auto protectedSource = source;
    protectedSource.protectedEdges.push_back(
        bnd::canonicalTetFacetEdge(3, 4));
    TetWorkspaceFixture blocked(protectedSource);
    auto blockedState =
        bnd::TetFacetConstraintState::forWorkspace(*blocked.workspace);
    const auto blockedRequest =
        makeRecoveryRequest(protectedSource, blocked, 9101, {2, 89});
    const auto fullBefore = paacs_tests::observableDigest(*blocked.workspace);
    const auto refused = bnd::recoverTetFacet(
        *blocked.workspace, blockedState, blockedRequest);
    CHECK(refused.status == bnd::TetFacetRecoveryStatus::Refused);
    CHECK(refused.reason ==
          bnd::TetFacetRecoveryReason::CrossedProtectedEdge);
    CHECK(paacs_tests::observableDigest(*blocked.workspace) == fullBefore);
    CHECK(blockedState.facets.empty());
    CHECK(blockedState.stamp == blockedRequest.stamp);
}

void protectedFaceIdentityAndReversal() {
    TEST("t2a.recovery.protected_face_identity_and_reversal");
    BipyramidFixture protectedSource;
    protectedSource.protectedFaces.push_back(
        {bnd::canonicalTetFacetFace({3, 0, 1}), 700, {2, 70}});
    TetWorkspaceFixture protectedFixture(protectedSource);
    const std::array<topo::NodeHandle, 3> protectedNodes{
        protectedFixture.nodes[3], protectedFixture.nodes[0],
        protectedFixture.nodes[1]};
    const auto protectedBefore = protectedFixture.workspace->faces().find(
        protectedFixture.workspace->nodes(), protectedNodes);
    CHECK(protectedBefore.has_value());
    auto protectedState = bnd::TetFacetConstraintState::forWorkspace(
        *protectedFixture.workspace);
    const auto protectedRequest =
        makeRecoveryRequest(protectedSource, protectedFixture, 9150, {2, 91});
    const auto protectedResult = bnd::recoverTetFacet(
        *protectedFixture.workspace, protectedState, protectedRequest);
    CHECK(protectedResult.status == bnd::TetFacetRecoveryStatus::Recovered);
    const auto protectedAfter = protectedFixture.workspace->faces().find(
        protectedFixture.workspace->nodes(), protectedNodes);
    CHECK(protectedAfter == protectedBefore);
    if (protectedAfter) {
        const auto data = protectedFixture.workspace->faces().get(
            *protectedAfter);
        CHECK(data.cls.dim == 2);
        CHECK(data.cls.tag == 70);
        CHECK(data.incident.size() == 1);
    }

    auto conflictSource = protectedSource;
    conflictSource.protectedFaces.push_back(
        {bnd::canonicalTetFacetFace({3, 4, 0}), 701, {2, 71}});
    TetWorkspaceFixture conflict(conflictSource);
    auto conflictState =
        bnd::TetFacetConstraintState::forWorkspace(*conflict.workspace);
    const auto conflictRequest =
        makeRecoveryRequest(conflictSource, conflict, 9151, {2, 92});
    const auto conflictDigest =
        paacs_tests::observableDigest(*conflict.workspace);
    const auto conflictResult = bnd::recoverTetFacet(
        *conflict.workspace, conflictState, conflictRequest);
    CHECK(conflictResult.status == bnd::TetFacetRecoveryStatus::Refused);
    CHECK(conflictResult.reason ==
          bnd::TetFacetRecoveryReason::CrossedProtectedFacet);
    CHECK(paacs_tests::observableDigest(*conflict.workspace) ==
          conflictDigest);

    BipyramidFixture forwardSource;
    TetWorkspaceFixture forward(forwardSource);
    TetWorkspaceFixture reversed(forwardSource);
    auto forwardState =
        bnd::TetFacetConstraintState::forWorkspace(*forward.workspace);
    auto reversedState =
        bnd::TetFacetConstraintState::forWorkspace(*reversed.workspace);
    auto forwardRequest =
        makeRecoveryRequest(forwardSource, forward, 9152, {2, 93});
    auto reversedRequest =
        makeRecoveryRequest(forwardSource, reversed, 9152, {2, 93});
    std::swap(reversedRequest.vertices[0], reversedRequest.vertices[1]);
    const auto forwardResult = bnd::recoverTetFacet(
        *forward.workspace, forwardState, forwardRequest);
    const auto reversedResult = bnd::recoverTetFacet(
        *reversed.workspace, reversedState, reversedRequest);
    CHECK(forwardResult.status == bnd::TetFacetRecoveryStatus::Recovered);
    CHECK(reversedResult.status == bnd::TetFacetRecoveryStatus::Recovered);
    CHECK(forwardResult.counters == reversedResult.counters);
    CHECK(forward.workspace->semanticTopologyHash() ==
          reversed.workspace->semanticTopologyHash());
}

void outsideTetIdentityAndPayloadInvariance() {
    TEST("t2a.recovery.outside_identity_payload_invariance");
    BipyramidFixture source;
    source.vertices.push_back({{10, 0, 10}, 60});
    source.vertices.push_back({{11, 0, 10}, 70});
    source.vertices.push_back({{10, 1, 10}, 80});
    source.vertices.push_back({{10, 0, 11}, 90});
    source.cells.push_back(positiveCell(source.vertices, {5, 6, 7, 8},
                                        130, {3, 9}));
    TetWorkspaceFixture fixture(source);
    auto state = bnd::TetFacetConstraintState::forWorkspace(*fixture.workspace);
    const auto request = makeRecoveryRequest(source, fixture, 9160, {2, 94});

    const auto outsideCell = fixture.cells.back();
    const auto outsideCellBefore = fixture.workspace->cells().get(outsideCell);
    const auto outsideCellSequence =
        fixture.workspace->cells().seqOf(outsideCell);
    const auto outsideEdges = fixture.workspace->edges().cellEdges(
        fixture.workspace->cells(), outsideCell);
    std::vector<topo::EdgePool::EdgeData> edgeDataBefore;
    for (const auto edge : outsideEdges)
        edgeDataBefore.push_back(fixture.workspace->edges().get(edge));
    const auto outsideFaces = fixture.workspace->faces().cellFaces(
        fixture.workspace->cells(), outsideCell);
    std::vector<topo::FacePool::FaceData> faceDataBefore;
    for (const auto face : outsideFaces)
        faceDataBefore.push_back(fixture.workspace->faces().get(face));
    std::vector<topo::NodePool::NodeData> nodeDataBefore;
    std::vector<std::uint64_t> nodeSequenceBefore;
    for (std::size_t i = 5; i < 9; ++i) {
        nodeDataBefore.push_back(
            fixture.workspace->nodes().get(fixture.nodes[i]));
        nodeSequenceBefore.push_back(
            fixture.workspace->nodes().seqOf(fixture.nodes[i]));
    }

    const auto recovered =
        bnd::recoverTetFacet(*fixture.workspace, state, request);
    CHECK(recovered.status == bnd::TetFacetRecoveryStatus::Recovered);
    CHECK(fixture.workspace->cells().isCurrent(outsideCell));
    CHECK(fixture.workspace->cells().seqOf(outsideCell) ==
          outsideCellSequence);
    const auto outsideCellAfter = fixture.workspace->cells().get(outsideCell);
    CHECK(outsideCellAfter.type == outsideCellBefore.type);
    CHECK(outsideCellAfter.count == outsideCellBefore.count);
    CHECK(outsideCellAfter.nodes == outsideCellBefore.nodes);
    CHECK(outsideCellAfter.cls.dim == outsideCellBefore.cls.dim);
    CHECK(outsideCellAfter.cls.tag == outsideCellBefore.cls.tag);
    const auto outsideEdgesAfter = fixture.workspace->edges().cellEdges(
        fixture.workspace->cells(), outsideCell);
    CHECK(outsideEdgesAfter.size() == outsideEdges.size());
    CHECK(std::equal(outsideEdgesAfter.begin(), outsideEdgesAfter.end(),
                     outsideEdges.begin(), outsideEdges.end()));
    for (std::size_t i = 0; i < outsideEdges.size(); ++i) {
        const auto after = fixture.workspace->edges().get(outsideEdges[i]);
        const auto& before = edgeDataBefore[i];
        CHECK(after.a == before.a);
        CHECK(after.b == before.b);
        CHECK(after.cls.dim == before.cls.dim);
        CHECK(after.cls.tag == before.cls.tag);
        CHECK(after.requiredConstraint == before.requiredConstraint);
        CHECK(after.useCount == before.useCount);
    }
    const auto outsideFacesAfter = fixture.workspace->faces().cellFaces(
        fixture.workspace->cells(), outsideCell);
    CHECK(outsideFacesAfter.size() == outsideFaces.size());
    CHECK(std::equal(outsideFacesAfter.begin(), outsideFacesAfter.end(),
                     outsideFaces.begin(), outsideFaces.end()));
    for (std::size_t i = 0; i < outsideFaces.size(); ++i) {
        const auto after = fixture.workspace->faces().get(outsideFaces[i]);
        const auto& before = faceDataBefore[i];
        CHECK(after.shape == before.shape);
        CHECK(after.cornerCount == before.cornerCount);
        CHECK(after.corners == before.corners);
        CHECK(after.cls.dim == before.cls.dim);
        CHECK(after.cls.tag == before.cls.tag);
        CHECK(after.incident.size() == before.incident.size());
        for (std::size_t j = 0; j < after.incident.size(); ++j) {
            CHECK(after.incident[j].cellSlot == before.incident[j].cellSlot);
            CHECK(after.incident[j].cellGen == before.incident[j].cellGen);
            CHECK(after.incident[j].sense == before.incident[j].sense);
        }
    }
    for (std::size_t i = 5; i < 9; ++i) {
        const auto after = fixture.workspace->nodes().get(fixture.nodes[i]);
        const auto& before = nodeDataBefore[i - 5];
        CHECK(fixture.workspace->nodes().isCurrent(fixture.nodes[i]));
        CHECK(fixture.workspace->nodes().seqOf(fixture.nodes[i]) ==
              nodeSequenceBefore[i - 5]);
        CHECK(after.coord.x == before.coord.x);
        CHECK(after.coord.y == before.coord.y);
        CHECK(after.coord.z == before.coord.z);
        CHECK(after.cls.dim == before.cls.dim);
        CHECK(after.cls.tag == before.cls.tag);
    }
}

void persistentStatePreflightAdversaries() {
    TEST("t2a.recovery.persistent_state_preflight");
    BipyramidFixture source;
    const auto plan = bnd::planTetFacetCavity(source.view());
    const auto candidate =
        bnd::reconstructTetFacetCandidate(source.view(), plan.plan);
    CHECK(candidate.built());
    if (!candidate.built()) return;
    source.cells = candidate.candidate.cells;
    TetWorkspaceFixture fixture(source);
    auto request = makeRecoveryRequest(source, fixture, 9200, {2, 90});
    auto state = bnd::TetFacetConstraintState::forWorkspace(*fixture.workspace);
    const auto marked =
        bnd::recoverTetFacet(*fixture.workspace, state, request);
    CHECK(marked.status == bnd::TetFacetRecoveryStatus::Marked);
    if (!marked.succeeded()) return;

    const auto stableDigest = paacs_tests::observableDigest(*fixture.workspace);
    const auto stableRevision = fixture.workspace->topologyRevision();

    auto staleState = state;
    --staleState.stamp.revision;
    request.stamp = topo::WorkspaceStamp::of(*fixture.workspace);
    const auto stale =
        bnd::recoverTetFacet(*fixture.workspace, staleState, request);
    CHECK(stale.status == bnd::TetFacetRecoveryStatus::Rejected);
    CHECK(stale.reason ==
          bnd::TetFacetRecoveryReason::StaleWorkspaceStamp);

    // T2A design section 4 makes request/state stamp validation the first
    // preflight gate.  A stale caller must not be allowed to probe handle
    // currency (or any later request property) through the returned reason.
    auto staleMalformedRequest = request;
    staleMalformedRequest.vertices[0] = topo::NodeHandle{};
    const auto staleMalformed = bnd::recoverTetFacet(
        *fixture.workspace, staleState, staleMalformedRequest);
    CHECK(staleMalformed.status == bnd::TetFacetRecoveryStatus::Rejected);
    CHECK(staleMalformed.reason ==
          bnd::TetFacetRecoveryReason::StaleWorkspaceStamp);

    auto wrongClass = state;
    wrongClass.facets.at(9200).classification = {2, 901};
    const auto wrong =
        bnd::recoverTetFacet(*fixture.workspace, wrongClass, request);
    CHECK(wrong.status == bnd::TetFacetRecoveryStatus::Rejected);
    CHECK(wrong.reason ==
          bnd::TetFacetRecoveryReason::WrongBindingClassification);

    auto keyMismatch = state;
    auto binding = keyMismatch.facets.extract(9200);
    binding.key() = 9201;
    keyMismatch.facets.insert(std::move(binding));
    const auto malformed =
        bnd::recoverTetFacet(*fixture.workspace, keyMismatch, request);
    CHECK(malformed.status == bnd::TetFacetRecoveryStatus::Rejected);
    CHECK(malformed.reason == bnd::TetFacetRecoveryReason::StateKeyMismatch);

    auto faceAlias = state;
    auto duplicate = faceAlias.facets.at(9200);
    duplicate.entityId = 9201;
    faceAlias.facets.emplace(9201, duplicate);
    const auto aliasedFace =
        bnd::recoverTetFacet(*fixture.workspace, faceAlias, request);
    CHECK(aliasedFace.status == bnd::TetFacetRecoveryStatus::Rejected);
    CHECK(aliasedFace.reason == bnd::TetFacetRecoveryReason::FaceAlias);

    auto targetAliasRequest = request;
    targetAliasRequest.entityId = 9201;
    const auto aliasedTarget =
        bnd::recoverTetFacet(*fixture.workspace, state, targetAliasRequest);
    CHECK(aliasedTarget.status == bnd::TetFacetRecoveryStatus::Rejected);
    CHECK(aliasedTarget.reason == bnd::TetFacetRecoveryReason::TargetAlias);

    auto wrongFace = state;
    const std::array<topo::NodeHandle, 3> otherFaceNodes{
        fixture.nodes[3], fixture.nodes[0], fixture.nodes[1]};
    const auto otherFace = fixture.workspace->faces().find(
        fixture.workspace->nodes(), otherFaceNodes);
    CHECK(otherFace.has_value());
    if (otherFace)
        wrongFace.facets.at(9200).recoveredFace = *otherFace;
    const auto wrongFaceResult =
        bnd::recoverTetFacet(*fixture.workspace, wrongFace, request);
    CHECK(wrongFaceResult.status == bnd::TetFacetRecoveryStatus::Rejected);
    CHECK(wrongFaceResult.reason ==
          bnd::TetFacetRecoveryReason::WrongBindingFace);

    auto rebindingRequest = request;
    rebindingRequest.vertices = {fixture.nodes[0], fixture.nodes[1],
                                 fixture.nodes[3]};
    const auto rebinding =
        bnd::recoverTetFacet(*fixture.workspace, state, rebindingRequest);
    CHECK(rebinding.status == bnd::TetFacetRecoveryStatus::Rejected);
    CHECK(rebinding.reason == bnd::TetFacetRecoveryReason::IdRebinding);
    CHECK(fixture.workspace->topologyRevision() == stableRevision);
    CHECK(paacs_tests::observableDigest(*fixture.workspace) == stableDigest);
}

void closedIdentityContactContracts() {
    TEST("t2a.geometry.closed_identity_contacts");
    const std::vector<bnd::TetFacetEdgeKey> noProtectedEdges;
    const std::vector<bnd::TetFacetProtectedFace> noProtectedFaces;

    // Identity-distinct vertex/edge contact: the lower Tet reaches z=0 only
    // at vertex 0, which lies in the relative interior of upper edge {4,5}.
    std::vector<bnd::TetFacetVertex> vertexEdgeVertices{
        {{0, 0, 0}, 10},   {{0, 0, -2}, 20}, {{1, 0, -2}, 30},
        {{0, 1, -2}, 40},  {{-2, 0, 0}, 50}, {{2, 0, 0}, 60},
        {{0, -1, 2}, 70},  {{0, 1, 2}, 80}};
    const auto lower = positiveCell(vertexEdgeVertices, {0, 1, 2, 3}, 100);
    const auto upper = positiveCell(vertexEdgeVertices, {4, 5, 6, 7}, 110);
    const bnd::TetFacetComplexView vertexEdgeView{
        vertexEdgeVertices, std::span<const bnd::TetFacetCell>{},
        {{4, 5, 6}}, noProtectedEdges, noProtectedFaces, {3, 9}};
    CHECK(bnd::tet_facet_detail::tetImproperlyIntersects(
        vertexEdgeView, lower, upper));

    // A shared vertex backed by the same identity is a legal closed contact.
    std::vector<bnd::TetFacetVertex> sharedVertexVertices{
        {{0, 0, 0}, 10}, {{1, 0, 0}, 20},  {{0, 1, 0}, 30},
        {{0, 0, 1}, 40}, {{-1, 0, 0}, 50}, {{0, -1, 0}, 60},
        {{0, 0, -1}, 70}};
    const auto positiveOctant =
        positiveCell(sharedVertexVertices, {0, 1, 2, 3}, 100);
    const auto negativeOctant =
        positiveCell(sharedVertexVertices, {0, 4, 5, 6}, 110);
    const bnd::TetFacetComplexView sharedVertexView{
        sharedVertexVertices, std::span<const bnd::TetFacetCell>{},
        {{0, 1, 2}}, noProtectedEdges, noProtectedFaces, {3, 9}};
    CHECK(!bnd::tet_facet_detail::tetImproperlyIntersects(
        sharedVertexView, positiveOctant, negativeOctant));
    CHECK(bnd::classifyTetFacetContact(sharedVertexView, negativeOctant) ==
          bnd::TetFacetContact::DeclaredBoundaryTouch);

    // A complete shared edge is legal; an unsplit transverse edge contact is
    // covered by the unchanged reviewer regression below and is illegal.
    std::vector<bnd::TetFacetVertex> sharedEdgeVertices{
        {{-1, 0, 0}, 10}, {{1, 0, 0}, 20},  {{0, 1, 0}, 30},
        {{0, 0, 1}, 40},  {{0, -1, 0}, 50}, {{0, 0, -1}, 60}};
    const auto edgeUpper =
        positiveCell(sharedEdgeVertices, {0, 1, 2, 3}, 100);
    const auto edgeLower =
        positiveCell(sharedEdgeVertices, {0, 1, 4, 5}, 110);
    const bnd::TetFacetComplexView sharedEdgeView{
        sharedEdgeVertices, std::span<const bnd::TetFacetCell>{},
        {{0, 1, 2}}, noProtectedEdges, noProtectedFaces, {3, 9}};
    CHECK(!bnd::tet_facet_detail::tetImproperlyIntersects(
        sharedEdgeView, edgeUpper, edgeLower));

    // A complete shared face with opposite owner senses remains legal.
    std::vector<bnd::TetFacetVertex> sharedFaceVertices{
        {{0, 0, 0}, 10}, {{2, 0, 0}, 20}, {{0, 2, 0}, 30},
        {{0, 0, 2}, 40}, {{0, 0, -2}, 50}};
    const auto faceUpper =
        positiveCell(sharedFaceVertices, {0, 1, 2, 3}, 100);
    const auto faceLower =
        positiveCell(sharedFaceVertices, {0, 1, 2, 4}, 110);
    const bnd::TetFacetComplexView sharedFaceView{
        sharedFaceVertices, std::span<const bnd::TetFacetCell>{},
        {{0, 1, 2}}, noProtectedEdges, noProtectedFaces, {3, 9}};
    CHECK(!bnd::tet_facet_detail::tetImproperlyIntersects(
        sharedFaceView, faceUpper, faceLower));
}

// External-review discriminators for revision 2. These pin three contracts
// that the production recovery route will rely on: candidate classifications
// are authoritative, every removed cavity node remains represented, and a
// geometric contact is legal only when the shared subentity has identical
// topology identity.
void reviewerCandidateContractAdversaries() {
    TEST("t2a.geometry.reviewer_candidate_contracts");

    BipyramidFixture fixture;
    const auto plan = bnd::planTetFacetCavity(fixture.view());
    const auto built = bnd::reconstructTetFacetCandidate(fixture.view(),
                                                         plan.plan);
    CHECK(plan.planned());
    CHECK(built.built());
    if (!built.built()) return;

    auto wrongClass = built.candidate.cells;
    wrongClass.front().classification = {3, 700};
    bnd::TetFacetGeometryCounters classCounters{};
    const auto classAudit = bnd::validateTetFacetCandidate(
        fixture.view(), plan.plan.selectedCells, wrongClass, classCounters);
    CHECK(!classAudit.valid);

    // Subdivide one removed Tet around a genuine interior node. The original
    // two-Tet candidate has the same cavity boundary but orphans that node.
    auto subdivided = fixture;
    subdivided.vertices.push_back({{1.5, 0.5, 0.6}, 60});
    const auto interior = static_cast<bnd::TetFacetVertexId>(
        subdivided.vertices.size() - 1);
    std::vector<bnd::TetFacetCell> refined;
    refined.push_back(positiveCell(subdivided.vertices,
                                   {interior, 4, 0, 1}, 101));
    refined.push_back(positiveCell(subdivided.vertices,
                                   {interior, 3, 1, 0}, 102));
    refined.push_back(positiveCell(subdivided.vertices,
                                   {interior, 3, 4, 1}, 103));
    refined.push_back(positiveCell(subdivided.vertices,
                                   {interior, 3, 0, 4}, 104));
    refined.push_back(subdivided.cells[1]);
    refined.push_back(subdivided.cells[2]);
    subdivided.cells = std::move(refined);
    std::vector<std::size_t> removed(subdivided.cells.size());
    for (std::size_t i = 0; i < removed.size(); ++i) removed[i] = i;
    bnd::TetFacetGeometryCounters orphanCounters{};
    const auto orphanAudit = bnd::validateTetFacetCandidate(
        subdivided.view(), removed, built.candidate.cells, orphanCounters);
    CHECK(!orphanAudit.valid);

    // These positive Tets meet only at one coordinate, but the coincident
    // vertices have different identities. That is not an identical shared
    // subentity and therefore must be reported as an improper intersection.
    std::vector<bnd::TetFacetVertex> contactVertices{
        {{0, 0, 0}, 10}, {{1, 0, 0}, 20}, {{0, 1, 0}, 30},
        {{0, 0, 1}, 40}, {{0, 0, 0}, 50}, {{-1, 0, 0}, 60},
        {{0, -1, 0}, 70}, {{0, 0, -1}, 80}};
    std::vector<bnd::TetFacetCell> contactCells{
        positiveCell(contactVertices, {0, 1, 2, 3}, 100),
        positiveCell(contactVertices, {4, 5, 6, 7}, 110)};
    const std::vector<bnd::TetFacetEdgeKey> noProtectedEdges;
    const std::vector<bnd::TetFacetProtectedFace> noProtectedFaces;
    const bnd::TetFacetComplexView contactView{
        contactVertices, contactCells, {{0, 1, 2}}, noProtectedEdges,
        noProtectedFaces, {3, 9}};
    CHECK(bnd::tet_facet_detail::tetImproperlyIntersects(
        contactView, contactCells[0], contactCells[1]));

    // A geometrically coincident face backed by distinct identities is not
    // the target FaceKey. It is an area-overlap degeneracy, not merely a
    // vertex touch.
    std::vector<bnd::TetFacetVertex> coincidentVertices{
        {{0, 0, 0}, 10}, {{2, 0, 0}, 20}, {{0, 2, 0}, 30},
        {{0, 0, 0}, 40}, {{2, 0, 0}, 50}, {{0, 2, 0}, 60},
        {{0, 0, 1}, 70}};
    std::vector<bnd::TetFacetCell> coincidentCells{
        positiveCell(coincidentVertices, {3, 4, 5, 6}, 100)};
    const bnd::TetFacetComplexView coincidentView{
        coincidentVertices, coincidentCells, {{0, 1, 2}},
        noProtectedEdges, noProtectedFaces, {3, 9}};
    CHECK(bnd::classifyTetFacetContact(coincidentView,
                                       coincidentCells.front()) ==
          bnd::TetFacetContact::CoplanarPartialOverlap);

    // A vertex-on-face touch is also illegal unless it is backed by the
    // identical topology vertex.  The two positive Tets below lie on
    // opposite sides of z=0 and meet only where vertex 0 lies strictly in
    // the interior of the other Tet's face {4,5,6}.  No vertex identity is
    // shared, so accepting this contact would create a hidden T-junction in
    // the candidate complex.
    std::vector<bnd::TetFacetVertex> vertexFaceVertices{
        {{0, 0, 0}, 10},   {{0, 0, -2}, 20}, {{1, 0, -2}, 30},
        {{0, 1, -2}, 40},  {{-2, -2, 0}, 50}, {{2, -2, 0}, 60},
        {{0, 2, 0}, 70},   {{0, 0, 2}, 80}};
    const auto lower = positiveCell(vertexFaceVertices, {0, 1, 2, 3}, 100);
    const auto upper = positiveCell(vertexFaceVertices, {4, 5, 6, 7}, 110);
    const bnd::TetFacetComplexView vertexFaceView{
        vertexFaceVertices, std::span<const bnd::TetFacetCell>{},
        {{4, 5, 6}}, noProtectedEdges, noProtectedFaces, {3, 9}};
    CHECK(bnd::tet_facet_detail::tetImproperlyIntersects(
        vertexFaceView, lower, upper));

    // The routing taxonomy must see the same closed-boundary contact.  This
    // target lies strictly below the Tet except for its first vertex, which
    // touches the interior of the Tet's z=0 face through a different node
    // identity.  It is tangential, not disjoint.
    std::vector<bnd::TetFacetVertex> targetTouchVertices{
        {{0, 0, 0}, 10},     {{1, 0, -100}, 20}, {{0, 1, -100}, 30},
        {{-2, -2, 0}, 40},   {{3, -2, 0}, 50},   {{0, 3, 0}, 60},
        {{0, 0, 2}, 70}};
    std::vector<bnd::TetFacetCell> targetTouchCells{
        positiveCell(targetTouchVertices, {3, 4, 5, 6}, 100)};
    const bnd::TetFacetComplexView targetTouchView{
        targetTouchVertices, targetTouchCells, {{0, 1, 2}},
        noProtectedEdges, noProtectedFaces, {3, 9}};
    CHECK(bnd::classifyTetFacetContact(targetTouchView,
                                       targetTouchCells.front()) ==
          bnd::TetFacetContact::TangentialTouch);

    // Likewise, an isolated transverse edge/edge touch is not a legal shared
    // edge when the endpoint identities differ.  Both Tets stay in opposite
    // closed half-spaces of z=0; their boundary edges cross at the origin.
    std::vector<bnd::TetFacetVertex> edgeEdgeVertices{
        {{-2, 0, 0}, 10}, {{2, 0, 0}, 20},  {{0, -1, -2}, 30},
        {{0, 1, -2}, 40}, {{0, -2, 0}, 50}, {{0, 2, 0}, 60},
        {{-1, 0, 2}, 70}, {{1, 0, 2}, 80}};
    const auto below = positiveCell(edgeEdgeVertices, {0, 1, 2, 3}, 100);
    const auto above = positiveCell(edgeEdgeVertices, {4, 5, 6, 7}, 110);
    const bnd::TetFacetComplexView edgeEdgeView{
        edgeEdgeVertices, std::span<const bnd::TetFacetCell>{},
        {{4, 5, 6}}, noProtectedEdges, noProtectedFaces, {3, 9}};
    CHECK(bnd::tet_facet_detail::tetImproperlyIntersects(
        edgeEdgeView, below, above));
}

}  // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };

const Entry* tetFacetTests(std::size_t& count) {
    static const Entry tests[] = {
        {"t2a.contact", exactContactTaxonomy},
        {"t2a.edge_star", highValenceEdgeStarClosure},
        {"t2a.fixed_point_seam", testingOnlySelectiveFixedPointClosure},
        {"t2a.determinism", deterministicClosureUnderChurn},
        {"t2a.candidate", ballAuditAndCandidateReconstruction},
        {"t2a.audits", candidateIntersectionAndProtectionAudits},
        {"t2a.tet_cavity", atomicTetCavityReplacement},
        {"t2a.recovery", structuredTetFacetRecoveryRoutes},
        {"t2a.recovery_high_valence",
         highValenceRecoveryAndProtectedRefusal},
        {"t2a.recovery_protected", protectedFaceIdentityAndReversal},
        {"t2a.recovery_outside", outsideTetIdentityAndPayloadInvariance},
        {"t2a.recovery_state", persistentStatePreflightAdversaries},
        {"t2a.closed_identity", closedIdentityContactContracts},
        {"t2a.reviewer_contracts", reviewerCandidateContractAdversaries},
    };
    count = sizeof(tests) / sizeof(tests[0]);
    return tests;
}
}  // namespace paacs_tests
