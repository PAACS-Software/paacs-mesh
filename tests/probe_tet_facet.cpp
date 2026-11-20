// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// T2-A standalone allocation qualification.
// This target replaces global operator new and is intentionally excluded from
// sanitizer builds.  It proves complete structured allocation attribution,
// phase-tier digest guarantees, same-workspace retry, all Tet/edge/face
// placeholder seams, and zero allocation after Phase 3 begins.
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <new>
#include <vector>

#include "paacs/bnd/TetFacetRecovery.hpp"
#include "workspace_digest.hpp"

static long gPhase3Allocations = 0;
static long gFailCountdown = 0;

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t bytes) {
    if (paacs::topo::detail::phase3Active()) ++gPhase3Allocations;
    if (gFailCountdown > 0 && --gFailCountdown == 0) throw std::bad_alloc{};
    if (void* memory = std::malloc(bytes ? bytes : 1)) return memory;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete(void* memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void* memory, std::size_t) noexcept {
    std::free(memory);
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
using namespace paacs;

int failures = 0;
#define REQUIRE(condition)                                                   \
    do {                                                                     \
        if (!(condition)) {                                                  \
            std::printf("FAIL %d: %s\n", __LINE__, #condition);            \
            ++failures;                                                      \
        }                                                                    \
    } while (0)

struct Fixture {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::vector<topo::NodeHandle> nodes;
    bnd::TetFacetConstraintState state;
    bnd::TetFacetRequest request;
};

std::array<topo::NodeHandle, 4> positive(
    const topo::NodePool& nodes,
    std::array<topo::NodeHandle, 4> connectivity) {
    const auto sign = predicates::orient3d(
        nodes.get(connectivity[0]).coord, nodes.get(connectivity[1]).coord,
        nodes.get(connectivity[2]).coord, nodes.get(connectivity[3]).coord);
    if (sign == predicates::Sign::Negative)
        std::swap(connectivity[0], connectivity[1]);
    return connectivity;
}

Fixture buildFixture() {
    Fixture fixture;
    topo::NodePool nodes;
    fixture.nodes.push_back(nodes.add({0, 0, 0}, {3, 9}));
    fixture.nodes.push_back(nodes.add({4, 0, 0}, {3, 9}));
    fixture.nodes.push_back(nodes.add({0, 4, 0}, {3, 9}));
    fixture.nodes.push_back(nodes.add({1, 1, 3}, {3, 9}));
    fixture.nodes.push_back(nodes.add({1, 1, -3}, {3, 9}));

    topo::CellPool cells;
    std::vector<topo::CellHandle> cellHandles;
    for (auto connectivity : {
             std::array<topo::NodeHandle, 4>{
                 fixture.nodes[3], fixture.nodes[4], fixture.nodes[0],
                 fixture.nodes[1]},
             std::array<topo::NodeHandle, 4>{
                 fixture.nodes[3], fixture.nodes[4], fixture.nodes[1],
                 fixture.nodes[2]},
             std::array<topo::NodeHandle, 4>{
                 fixture.nodes[3], fixture.nodes[4], fixture.nodes[2],
                 fixture.nodes[0]}}) {
        connectivity = positive(nodes, connectivity);
        cellHandles.push_back(cells.add(topo::ElementType::Tet4,
                                        connectivity, {3, 9}, nodes));
    }
    topo::EdgePool edges;
    for (const auto cell : cellHandles)
        (void)edges.insertCellEdges(nodes, cells, cell);
    for (std::size_t i = 0; i < 3; ++i) {
        const auto edge = edges.find(nodes, fixture.nodes[i],
                                     fixture.nodes[(i + 1) % 3]);
        if (!edge) throw std::logic_error("probe target edge missing");
        edges.setRequiredConstraint(*edge, true);
        edges.setClassification(*edge, {1, 41});
    }
    topo::FacePool faces;
    for (const auto cell : cellHandles)
        (void)faces.insertCellFaces(nodes, cells, edges, cell);
    fixture.workspace = std::make_unique<topo::TopologyWorkspace>(
        std::move(nodes), std::move(cells), std::move(edges),
        std::move(faces));
    fixture.state =
        bnd::TetFacetConstraintState::forWorkspace(*fixture.workspace);
    fixture.request.stamp = topo::WorkspaceStamp::of(*fixture.workspace);
    fixture.request.entityId = 7001;
    fixture.request.vertices = {fixture.nodes[0], fixture.nodes[1],
                                fixture.nodes[2]};
    fixture.request.classification = {2, 77};
    return fixture;
}

bool sameState(const bnd::TetFacetConstraintState& a,
               const bnd::TetFacetConstraintState& b) {
    if (!(a.stamp == b.stamp) || a.facets.size() != b.facets.size())
        return false;
    auto ia = a.facets.begin();
    auto ib = b.facets.begin();
    for (; ia != a.facets.end(); ++ia, ++ib) {
        const auto& x = ia->second;
        const auto& y = ib->second;
        if (ia->first != ib->first || x.entityId != y.entityId ||
            x.canonicalVertices != y.canonicalVertices ||
            x.classification.dim != y.classification.dim ||
            x.classification.tag != y.classification.tag ||
            x.recoveredFace != y.recoveredFace)
            return false;
    }
    return true;
}

void requireRetry(Fixture& fixture) {
    const auto retry = bnd::recoverTetFacet(*fixture.workspace, fixture.state,
                                            fixture.request);
    REQUIRE(retry.status == bnd::TetFacetRecoveryStatus::Recovered);
    REQUIRE(retry.face.isSpecified());
    REQUIRE(fixture.workspace->state() ==
            topo::WorkspaceState::CommittedAdjacencyStale);
    REQUIRE(topo::workspaceTopologyReady(
        fixture.workspace->nodes(), fixture.workspace->cells(),
        fixture.workspace->edges(), fixture.workspace->faces()));
}

}  // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    static_assert(noexcept(bnd::recoverTetFacet(
        std::declval<topo::TopologyWorkspace&>(),
        std::declval<bnd::TetFacetConstraintState&>(),
        std::declval<const bnd::TetFacetRequest&>())));

    // A normal run instruments the complete commit boundary.
    {
        auto fixture = buildFixture();
        gPhase3Allocations = 0;
        const auto result = bnd::recoverTetFacet(
            *fixture.workspace, fixture.state, fixture.request);
        REQUIRE(result.status == bnd::TetFacetRecoveryStatus::Recovered);
        REQUIRE(gPhase3Allocations == 0);
    }

    long allocationFaults = 0;
    long fullTier = 0;
    long logicalTier = 0;
    for (long ordinal = 1;; ++ordinal) {
        auto fixture = buildFixture();
        const auto stateBefore = fixture.state;
        const auto fullBefore =
            paacs_tests::observableDigest(*fixture.workspace);
        const auto logicalBefore =
            paacs_tests::observableLogicalDigest(*fixture.workspace);
        gFailCountdown = ordinal;
        const auto result = bnd::recoverTetFacet(
            *fixture.workspace, fixture.state, fixture.request);
        gFailCountdown = 0;
        if (result.status == bnd::TetFacetRecoveryStatus::Recovered) break;
        ++allocationFaults;
        REQUIRE(result.status == bnd::TetFacetRecoveryStatus::Failed);
        REQUIRE(result.reason == bnd::TetFacetRecoveryReason::AllocationFailed);
        REQUIRE(sameState(fixture.state, stateBefore));
        if (result.capacityResiduePermitted) {
            ++logicalTier;
            REQUIRE(paacs_tests::observableLogicalDigest(*fixture.workspace) ==
                    logicalBefore);
        } else {
            ++fullTier;
            REQUIRE(paacs_tests::observableDigest(*fixture.workspace) ==
                    fullBefore);
        }
        requireRetry(fixture);
    }

    // The explicit preparation seam fires immediately after each table and
    // per-cell incidence placeholder, even when reserve made the associated
    // insertion allocation-free.
    long placeholderFaults = 0;
    for (int ordinal = 1;; ++ordinal) {
        auto fixture = buildFixture();
        const auto stateBefore = fixture.state;
        const auto logicalBefore =
            paacs_tests::observableLogicalDigest(*fixture.workspace);
        topo::faultCountdown() = ordinal;
        const auto result = bnd::recoverTetFacet(
            *fixture.workspace, fixture.state, fixture.request);
        topo::faultCountdown() = 0;
        if (result.status == bnd::TetFacetRecoveryStatus::Recovered) break;
        ++placeholderFaults;
        REQUIRE(result.status == bnd::TetFacetRecoveryStatus::Failed);
        REQUIRE(result.reason == bnd::TetFacetRecoveryReason::AllocationFailed);
        REQUIRE(sameState(fixture.state, stateBefore));
        REQUIRE(paacs_tests::observableLogicalDigest(*fixture.workspace) ==
                logicalBefore);
        requireRetry(fixture);
    }

    std::printf(
        "tet-facet sweep: %ld allocation faults (%ld full/%ld logical), "
        "%ld placeholder faults, %ld Phase-3 allocations\n",
        allocationFaults, fullTier, logicalTier, placeholderFaults,
        gPhase3Allocations);
    REQUIRE(allocationFaults >= 100);
    REQUIRE(fullTier > 0);
    REQUIRE(logicalTier > 0);
    REQUIRE(placeholderFaults >= 5);
    if (failures != 0) {
        std::printf("TET FACET PROBE FAILED (%d)\n", failures);
        return 1;
    }
    std::printf("TET FACET PROBE PASSED\n");
    return 0;
}
