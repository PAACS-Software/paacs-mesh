// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// T2-B indexed mutation allocation qualification.  The end-to-end v3 box
// probe exercises every pipeline phase but all of its convex boundary
// constraints are already present.  These two small fixtures force the new
// target-edge transaction and the indexed facet recovery/update paths.
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <new>
#include <stdexcept>
#include <vector>

#include "paacs/bnd/TetFacetSetRecovery.hpp"
#include "workspace_digest.hpp"

static std::size_t gCalls = 0;
static std::size_t gFailAt = 0;
static std::size_t gPhase3Allocations = 0;

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t bytes) {
    ++gCalls;
    if (paacs::topo::detail::phase3Active()) ++gPhase3Allocations;
    if (gFailAt != 0 && gCalls == gFailAt) throw std::bad_alloc{};
    if (void* memory = std::malloc(bytes == 0 ? 1 : bytes)) return memory;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete(void* memory,std::size_t) noexcept { std::free(memory); }
void operator delete[](void* memory,std::size_t) noexcept { std::free(memory); }
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
using namespace paacs;
int failures = 0;
#define REQUIRE(condition)                                                   \
    do {                                                                     \
        if (!(condition)) {                                                  \
            std::printf("FAIL %d: %s\n",__LINE__,#condition);             \
            ++failures;                                                      \
        }                                                                    \
    } while (0)

std::array<topo::NodeHandle,4> positive(
    const topo::NodePool& nodes,std::array<topo::NodeHandle,4> value) {
    if (predicates::orient3d(
            nodes.get(value[0]).coord,nodes.get(value[1]).coord,
            nodes.get(value[2]).coord,nodes.get(value[3]).coord) ==
        predicates::Sign::Negative)
        std::swap(value[0],value[1]);
    return value;
}

bool sameState(const bnd::TetFacetSetState& a,
               const bnd::TetFacetSetState& b) {
    if (!(a.workspaceStamp == b.workspaceStamp) ||
        !(a.artifactStamp == b.artifactStamp) ||
        a.canonicalEdgeSchedule != b.canonicalEdgeSchedule ||
        a.canonicalSchedule != b.canonicalSchedule ||
        a.edgeBindings.size() != b.edgeBindings.size() ||
        a.bindings.size() != b.bindings.size())
        return false;
    auto ea = a.edgeBindings.begin(),eb = b.edgeBindings.begin();
    for (;ea != a.edgeBindings.end();++ea,++eb) {
        const auto& x = ea->second; const auto& y = eb->second;
        if (!(ea->first == eb->first) || !(x.key == y.key) ||
            x.classification.dim != y.classification.dim ||
            x.classification.tag != y.classification.tag ||
            x.required != y.required || x.edge != y.edge) return false;
    }
    auto fa = a.bindings.begin(),fb = b.bindings.begin();
    for (;fa != a.bindings.end();++fa,++fb) {
        const auto& x = fa->second; const auto& y = fb->second;
        if (!(fa->first == fb->first) || !(x.key == y.key) ||
            x.canonicalVertices != y.canonicalVertices ||
            x.classification.dim != y.classification.dim ||
            x.classification.tag != y.classification.tag ||
            x.recoveredFace != y.recoveredFace) return false;
    }
    return true;
}

struct EdgeFixture {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::unique_ptr<bnd::TetRecoveryContext> context;
    std::vector<topo::NodeHandle> nodes;
    bnd::TetFacetSetState state;
    bnd::TetEdgeRecoveryRequest request;
};

EdgeFixture edgeFixture() {
    EdgeFixture fixture;
    topo::NodePool nodes;
    for (const Vec3& p : std::array<Vec3,5>{
             Vec3{0,0,0},Vec3{4,0,0},Vec3{0,4,0},
             Vec3{1,1,3},Vec3{1,1,-3}})
        fixture.nodes.push_back(nodes.add(p,{3,9}));
    topo::CellPool cells;
    for (auto value : std::array<std::array<std::size_t,4>,2>{
             std::array<std::size_t,4>{3,0,1,2},
             std::array<std::size_t,4>{4,0,2,1}})
        (void)cells.add(topo::ElementType::Tet4,
            positive(nodes,{fixture.nodes[value[0]],fixture.nodes[value[1]],
                            fixture.nodes[value[2]],fixture.nodes[value[3]]}),
            {3,9},nodes);
    topo::EdgePool edges;
    for (const auto cell : cells.activeBySequence())
        (void)edges.insertCellEdges(nodes,cells,cell);
    topo::FacePool faces;
    for (const auto cell : cells.activeBySequence())
        (void)faces.insertCellFaces(nodes,cells,edges,cell);
    fixture.workspace = std::make_unique<topo::TopologyWorkspace>(
        std::move(nodes),std::move(cells),std::move(edges),std::move(faces));
    auto indexed = bnd::buildTetRecoveryContext(*fixture.workspace,{3,9});
    if (!indexed.built()) throw std::logic_error("edge index fixture");
    fixture.context = std::move(indexed.context);
    fixture.state.workspaceStamp = topo::WorkspaceStamp::of(*fixture.workspace);
    fixture.state.artifactStamp = {700,1};
    fixture.request.stamp = fixture.state.workspaceStamp;
    fixture.request.artifactStamp = fixture.state.artifactStamp;
    fixture.request.key = bnd::VolumeArtifactEdgeKey::of(30,40);
    fixture.request.vertices = {fixture.nodes[3],fixture.nodes[4]};
    fixture.request.classification = {1,55};
    fixture.request.region = {3,9};
    return fixture;
}

topo::TetEdgeCavitySpec dim2CavitySpec(const EdgeFixture& fixture) {
    topo::TetEdgeCavitySpec spec;
    spec.stamp=topo::WorkspaceStamp::of(*fixture.workspace);
    spec.removed=fixture.workspace->cells().activeBySequence();
    spec.target={fixture.nodes[3],fixture.nodes[4]};
    spec.targetClassification={2,155};
    spec.region={3,9};
    constexpr std::array<std::array<std::size_t,4>,3> added{{
        {{3,4,0,1}},{{3,4,1,2}},{{3,4,2,0}}}};
    for (const auto ids : added) {
        topo::TetCavityCell cell;
        cell.nodes=positive(fixture.workspace->nodes(),{
            fixture.nodes[ids[0]],fixture.nodes[ids[1]],
            fixture.nodes[ids[2]],fixture.nodes[ids[3]]});
        cell.cls={3,9};
        spec.added.push_back(cell);
    }
    return spec;
}

enum class Dim2DirectStatus { Committed,AllocationFailed,Unexpected };

Dim2DirectStatus runDim2Direct(
    EdgeFixture& fixture,std::size_t allocationOrdinal,
    int placeholderOrdinal) {
    auto spec=dim2CavitySpec(fixture);
    auto transaction=fixture.workspace->beginTetCavityTransaction();
    try {
        (void)transaction.stage(spec);
        gCalls=0;
        gFailAt=allocationOrdinal;
        topo::faultCountdown()=placeholderOrdinal;
        transaction.prepare();
    } catch (const std::bad_alloc&) {
        gFailAt=0;
        topo::faultCountdown()=0;
        transaction.abandon();
        return Dim2DirectStatus::AllocationFailed;
    } catch (...) {
        gFailAt=0;
        topo::faultCountdown()=0;
        transaction.abandon();
        return Dim2DirectStatus::Unexpected;
    }
    gFailAt=0;
    topo::faultCountdown()=0;
    const auto predicted=transaction.predictedTargetEdge();
    const auto committed=transaction.commitEdge();
    const auto actual=fixture.workspace->edges().find(
        fixture.workspace->nodes(),fixture.nodes[3],fixture.nodes[4]);
    REQUIRE(predicted.isSpecified());
    REQUIRE(committed.targetEdge==predicted);
    REQUIRE(actual&&*actual==predicted);
    if (actual) {
        const auto data=fixture.workspace->edges().get(*actual);
        REQUIRE(data.cls.dim==2&&data.cls.tag==155);
        REQUIRE(data.requiredConstraint);
    }
    return Dim2DirectStatus::Committed;
}

std::size_t sweepDim2Allocations() {
    std::size_t faults=0;
    for (std::size_t ordinal=1;;++ordinal) {
        auto fixture=edgeFixture();
        const auto before=paacs_tests::observableLogicalDigest(
            *fixture.workspace);
        const auto result=runDim2Direct(fixture,ordinal,0);
        if (result==Dim2DirectStatus::Committed) break;
        ++faults;
        REQUIRE(result==Dim2DirectStatus::AllocationFailed);
        REQUIRE(paacs_tests::observableLogicalDigest(*fixture.workspace)==before);
        REQUIRE(runDim2Direct(fixture,0,0)==Dim2DirectStatus::Committed);
        if (failures!=0) break;
    }
    return faults;
}

std::size_t sweepDim2Placeholders() {
    std::size_t faults=0;
    for (int ordinal=1;;++ordinal) {
        auto fixture=edgeFixture();
        const auto before=paacs_tests::observableLogicalDigest(
            *fixture.workspace);
        const auto result=runDim2Direct(fixture,0,ordinal);
        if (result==Dim2DirectStatus::Committed) break;
        ++faults;
        REQUIRE(result==Dim2DirectStatus::AllocationFailed);
        REQUIRE(paacs_tests::observableLogicalDigest(*fixture.workspace)==before);
        REQUIRE(runDim2Direct(fixture,0,0)==Dim2DirectStatus::Committed);
        if (failures!=0) break;
    }
    return faults;
}

struct FacetFixture {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::unique_ptr<bnd::TetRecoveryContext> context;
    std::vector<topo::NodeHandle> nodes;
    bnd::TetFacetSetState state;
    bnd::TetFacetIndexedRequest request;
};

FacetFixture facetFixture() {
    FacetFixture fixture;
    topo::NodePool nodes;
    for (const Vec3& p : std::array<Vec3,5>{
             Vec3{0,0,0},Vec3{4,0,0},Vec3{0,4,0},
             Vec3{1,1,3},Vec3{1,1,-3}})
        fixture.nodes.push_back(nodes.add(p,{3,9}));
    topo::CellPool cells;
    for (const auto ids : std::array<std::array<std::size_t,4>,3>{
             std::array<std::size_t,4>{3,4,0,1},
             std::array<std::size_t,4>{3,4,1,2},
             std::array<std::size_t,4>{3,4,2,0}})
        (void)cells.add(topo::ElementType::Tet4,
            positive(nodes,{fixture.nodes[ids[0]],fixture.nodes[ids[1]],
                            fixture.nodes[ids[2]],fixture.nodes[ids[3]]}),
            {3,9},nodes);
    topo::EdgePool edges;
    for (const auto cell : cells.activeBySequence())
        (void)edges.insertCellEdges(nodes,cells,cell);
    for (std::size_t i = 0;i < 3;++i) {
        const auto edge = edges.find(
            nodes,fixture.nodes[i],fixture.nodes[(i+1)%3]);
        if (!edge) throw std::logic_error("facet edge fixture");
        edges.setClassification(*edge,{1,41});
        edges.setRequiredConstraint(*edge,true);
    }
    topo::FacePool faces;
    for (const auto cell : cells.activeBySequence())
        (void)faces.insertCellFaces(nodes,cells,edges,cell);
    fixture.workspace = std::make_unique<topo::TopologyWorkspace>(
        std::move(nodes),std::move(cells),std::move(edges),std::move(faces));
    auto indexed = bnd::buildTetRecoveryContext(*fixture.workspace,{3,9});
    if (!indexed.built()) throw std::logic_error("facet index fixture");
    fixture.context = std::move(indexed.context);
    fixture.state.workspaceStamp = topo::WorkspaceStamp::of(*fixture.workspace);
    fixture.state.artifactStamp = {701,1};
    fixture.request.stamp = fixture.state.workspaceStamp;
    fixture.request.artifactStamp = fixture.state.artifactStamp;
    fixture.request.key = {800,0};
    fixture.request.vertices = {
        fixture.nodes[0],fixture.nodes[1],fixture.nodes[2]};
    fixture.request.classification = {2,77};
    fixture.request.region = {3,9};
    return fixture;
}

template<class Fixture,class Recover,class Succeeded,class Allocation>
std::array<std::size_t,4> sweep(
    Fixture (*makeFixture)(),Recover recover,Succeeded succeeded,
    Allocation allocation) {
    std::size_t faults = 0,full = 0,logical = 0,postCommit = 0;
    for (std::size_t ordinal = 1;;++ordinal) {
        auto fixture = makeFixture();
        const auto stateBefore = fixture.state;
        const auto fullBefore =
            paacs_tests::observableDigest(*fixture.workspace);
        const auto logicalBefore =
            paacs_tests::observableLogicalDigest(*fixture.workspace);
        const auto revision = fixture.workspace->topologyRevision();
        gCalls = 0; gFailAt = ordinal;
        const auto outcome = recover(fixture);
        gFailAt = 0;
        if (succeeded(outcome)) break;
        ++faults;
        REQUIRE(allocation(outcome));
        REQUIRE(sameState(fixture.state,stateBefore));
        if (fixture.workspace->topologyRevision() == revision) {
            REQUIRE(paacs_tests::observableLogicalDigest(*fixture.workspace) ==
                    logicalBefore);
            if (paacs_tests::observableDigest(*fixture.workspace) == fullBefore)
                ++full;
            else
                ++logical;
            const auto retry = recover(fixture);
            REQUIRE(succeeded(retry));
        } else {
            ++postCommit;
            REQUIRE(fixture.context->poisoned);
            REQUIRE(topo::workspaceTopologyReady(
                fixture.workspace->nodes(),fixture.workspace->cells(),
                fixture.workspace->edges(),fixture.workspace->faces()));
            auto retryFixture = makeFixture();
            const auto retry = recover(retryFixture);
            REQUIRE(succeeded(retry));
        }
        if (failures != 0) break;
    }
    return {faults,full,logical,postCommit};
}

template<class Fixture,class Recover,class Succeeded,class Allocation>
std::size_t placeholderSweep(
    Fixture (*makeFixture)(),Recover recover,Succeeded succeeded,
    Allocation allocation) {
    std::size_t faults = 0;
    for (int ordinal = 1;;++ordinal) {
        auto fixture = makeFixture();
        const auto stateBefore = fixture.state;
        const auto logicalBefore =
            paacs_tests::observableLogicalDigest(*fixture.workspace);
        topo::faultCountdown() = ordinal;
        const auto outcome = recover(fixture);
        topo::faultCountdown() = 0;
        if (succeeded(outcome)) break;
        ++faults;
        REQUIRE(allocation(outcome));
        REQUIRE(sameState(fixture.state,stateBefore));
        REQUIRE(paacs_tests::observableLogicalDigest(*fixture.workspace) ==
                logicalBefore);
        const auto retry = recover(fixture);
        REQUIRE(succeeded(retry));
        if (failures != 0) break;
    }
    return faults;
}
}  // namespace

int main() {
    std::setvbuf(stdout,nullptr,_IONBF,0);
    auto edgeRecover = [](EdgeFixture& f) {
        return bnd::recoverTetEdgeIndexed(
            f.request,f.state,*f.workspace,*f.context);
    };
    auto edgeSucceeded = [](const auto& out) {
        return out.status == bnd::TetEdgeRecoveryStatus::Recovered;
    };
    auto edgeAllocation = [](const auto& out) {
        return out.status == bnd::TetEdgeRecoveryStatus::Failed &&
               out.reason == bnd::TetEdgeRecoveryReason::AllocationFailed;
    };
    auto facetRecover = [](FacetFixture& f) {
        return bnd::recoverTetFacetIndexed(
            f.request,f.state,*f.workspace,*f.context);
    };
    auto facetSucceeded = [](const auto& out) {
        return out.recovery.status == bnd::TetFacetRecoveryStatus::Recovered;
    };
    auto facetAllocation = [](const auto& out) {
        return out.recovery.status == bnd::TetFacetRecoveryStatus::Failed &&
               out.recovery.reason ==
                   bnd::TetFacetRecoveryReason::AllocationFailed;
    };

    gPhase3Allocations = 0;
    {
        auto fixture = edgeFixture();
        const auto result = edgeRecover(fixture);
        REQUIRE(edgeSucceeded(result));
        const auto actual = fixture.workspace->edges().find(
            fixture.workspace->nodes(),fixture.nodes[3],fixture.nodes[4]);
        REQUIRE(actual && *actual == result.edge);
    }
    {
        auto fixture = facetFixture();
        const auto result = facetRecover(fixture);
        REQUIRE(facetSucceeded(result));
        REQUIRE(result.recovery.face.isSpecified());
    }
    REQUIRE(gPhase3Allocations == 0);

    const auto edge = sweep(
        edgeFixture,edgeRecover,edgeSucceeded,edgeAllocation);
    const auto facet = sweep(
        facetFixture,facetRecover,facetSucceeded,facetAllocation);
    const auto edgePlaceholders = placeholderSweep(
        edgeFixture,edgeRecover,edgeSucceeded,edgeAllocation);
    const auto facetPlaceholders = placeholderSweep(
        facetFixture,facetRecover,facetSucceeded,facetAllocation);
    const auto dim2Allocations=sweepDim2Allocations();
    const auto dim2Placeholders=sweepDim2Placeholders();
    REQUIRE(edge[0] > 0 && facet[0] > 0);
    REQUIRE(edge[1]+edge[2] > 0 && facet[1]+facet[2] > 0);
    REQUIRE(edgePlaceholders > 0 && facetPlaceholders > 0);
    REQUIRE(dim2Allocations > 0 && dim2Placeholders > 0);
    REQUIRE(gPhase3Allocations == 0);
    std::printf("t2b indexed faults: edge=%zu (%zu full/%zu logical/%zu "
                "postcommit) facet=%zu (%zu/%zu/%zu) placeholders=%zu/%zu "
                "dim2=%zu/%zu phase3=%zu\n",
        edge[0],edge[1],edge[2],edge[3],facet[0],facet[1],facet[2],facet[3],
        edgePlaceholders,facetPlaceholders,dim2Allocations,dim2Placeholders,
        gPhase3Allocations);
    if (failures != 0) return 1;
    std::printf("probe_t2b_indexed PASSED\n");
    return 0;
}
