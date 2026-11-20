// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Focused allocation qualification for the v4-private mixed exact-present
// scheduler.  Edge and facet fixtures each contain exactly one binding batch.
// Sweeping every allocation on their clean paths proves rollback and retry;
// the dedicated mixed-binding commit-zone marker separately proves that the
// report/state/reverse-index publication tail performs no allocation.

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <map>
#include <memory>
#include <new>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#if defined(_WIN32)
#include <malloc.h>
#endif

#include "paacs/bnd/TetFacetSetRecovery.hpp"
#include "workspace_digest.hpp"

namespace {
std::size_t gCalls=0;
std::size_t gBytes=0;
std::size_t gFailAt=0;
std::size_t gThrown=0;
std::size_t gPhase3Allocations=0;
std::size_t gMixedCommitAllocations=0;
std::size_t gCurrentOrdinal=0;
int gCurrentKind=-1;

void noteAllocation(std::size_t bytes) {
    ++gCalls;
    gBytes+=bytes;
    if(paacs::topo::detail::phase3Active())++gPhase3Allocations;
    if(paacs::bnd::mixedPresentCommitZoneActiveForTesting())
        ++gMixedCommitAllocations;
    if(gFailAt!=0&&gCalls==gFailAt) {
        ++gThrown;
        throw std::bad_alloc{};
    }
}

void* allocateAligned(std::size_t bytes,std::size_t alignment) {
    noteAllocation(bytes);
    const auto count=bytes==0 ? std::size_t{1}:bytes;
#if defined(_WIN32)
    if(void* memory=_aligned_malloc(count,alignment))return memory;
#else
    void* memory=nullptr;
    if(posix_memalign(&memory,alignment,count)==0)return memory;
#endif
    throw std::bad_alloc{};
}
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t bytes) {
    noteAllocation(bytes);
    if(void* memory=std::malloc(bytes==0 ? 1:bytes))return memory;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void* operator new(std::size_t bytes,std::align_val_t alignment) {
    return allocateAligned(bytes,static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t bytes,std::align_val_t alignment) {
    return ::operator new(bytes,alignment);
}
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete(void* memory,std::size_t) noexcept { std::free(memory); }
void operator delete[](void* memory,std::size_t) noexcept {
    std::free(memory);
}
void operator delete(void* memory,std::align_val_t) noexcept {
#if defined(_WIN32)
    _aligned_free(memory);
#else
    std::free(memory);
#endif
}
void operator delete[](void* memory,std::align_val_t alignment) noexcept {
    ::operator delete(memory,alignment);
}
void operator delete(void* memory,std::size_t,
                     std::align_val_t alignment) noexcept {
    ::operator delete(memory,alignment);
}
void operator delete[](void* memory,std::size_t,
                       std::align_val_t alignment) noexcept {
    ::operator delete(memory,alignment);
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
using namespace paacs;
int failures=0;
#define REQUIRE(condition)                                                   \
    do {                                                                     \
        if(!(condition)) {                                                   \
            std::printf("FAIL %d: %s\n",__LINE__,#condition);              \
            ++failures;                                                      \
        }                                                                    \
    } while(false)

enum class BatchKind { Edge,Facet };

std::array<topo::NodeHandle,4> positive(
    const topo::NodePool& nodes,std::array<topo::NodeHandle,4> value) {
    if(predicates::orient3d(
           nodes.get(value[0]).coord,nodes.get(value[1]).coord,
           nodes.get(value[2]).coord,nodes.get(value[3]).coord)==
       predicates::Sign::Negative)
        std::swap(value[0],value[1]);
    return value;
}

struct Fixture {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::unique_ptr<bnd::TetRecoveryContext> context;
    bnd::TetFacetSetState state;
    bnd::ClosedTri3Artifact artifact;
    std::map<Index,topo::NodeHandle> localByGlobal;
};

Fixture makeFixture(BatchKind kind) {
    Fixture fixture;
    topo::NodePool nodes;
    std::array<topo::NodeHandle,5> handles{
        nodes.add({0,0,0},{3,9}),nodes.add({2,0,0},{3,9}),
        nodes.add({0,2,0},{3,9}),nodes.add({0,0,2},{3,9}),
        nodes.add({0,0,-2},{3,9})};
    topo::CellPool cells;
    (void)cells.add(topo::ElementType::Tet4,
        positive(nodes,{handles[0],handles[1],handles[2],handles[3]}),
        {3,9},nodes);
    (void)cells.add(topo::ElementType::Tet4,
        positive(nodes,{handles[0],handles[2],handles[1],handles[4]}),
        {3,9},nodes);
    auto edges=topo::EdgePool::buildComplete(nodes,cells);
    const std::array<std::array<std::size_t,2>,3> faceEdges{{
        {{0,1}},{{1,2}},{{2,0}}}};
    for(const auto ids:faceEdges) {
        const auto edge=edges.find(nodes,handles[ids[0]],handles[ids[1]]);
        if(!edge)throw std::logic_error("missing exact-present edge");
        edges.setClassification(*edge,kind==BatchKind::Edge
            ? mesh::Classification{1,41}:mesh::Classification{2,79});
        edges.setRequiredConstraint(*edge,true);
    }
    auto faces=topo::FacePool::buildComplete(nodes,cells,edges);
    if(kind==BatchKind::Facet) {
        const auto face=faces.find(nodes,
            std::array<topo::NodeHandle,3>{handles[0],handles[1],handles[2]});
        if(!face)throw std::logic_error("missing exact-present facet");
        faces.setClassification(*face,{2,79});
    }
    fixture.workspace=std::make_unique<topo::TopologyWorkspace>(
        std::move(nodes),std::move(cells),std::move(edges),std::move(faces));
    auto indexed=bnd::buildTetRecoveryContext(*fixture.workspace,{3,9});
    if(!indexed.built()||!indexed.context)
        throw std::logic_error("mixed-present index fixture");
    fixture.context=std::move(indexed.context);
    fixture.state.workspaceStamp=topo::WorkspaceStamp::of(*fixture.workspace);
    fixture.state.artifactStamp={900+(kind==BatchKind::Edge ? 0U:1U),1};
    fixture.artifact.surfaceSetStamp=fixture.state.artifactStamp;
    fixture.artifact.regionTags={9};
    if(kind==BatchKind::Edge) {
        for(const auto ids:faceEdges) {
            const auto key=bnd::VolumeArtifactEdgeKey::of(
                static_cast<Index>(ids[0]),static_cast<Index>(ids[1]));
            fixture.artifact.edges.push_back({key,{1,41},{}});
        }
        std::sort(fixture.artifact.edges.begin(),fixture.artifact.edges.end(),
            [](const auto& a,const auto& b) { return a.key<b.key; });
    } else {
        fixture.artifact.facets.push_back(
            {{910,0},{0,1,2},{2,79},{{9,1}}});
        // Seed one valid retained binding so the mixed entry receipt has a
        // real reserve/copy surface.  The sweep below injects at its first
        // allocation, before either canonical schedule is prepared.
        std::array<topo::NodeHandle,3> triple{
            handles[0],handles[1],handles[2]};
        std::sort(triple.begin(),triple.end(),[&](const auto a,const auto b) {
            return fixture.workspace->nodes().seqOf(a)<
                   fixture.workspace->nodes().seqOf(b);
        });
        const auto face=fixture.workspace->faces().find(
            fixture.workspace->nodes(),triple);
        if(!face)throw std::logic_error("missing retained entry facet");
        const bnd::VolumeFacetKey key{910,0};
        fixture.state.bindings.emplace(key,
            bnd::PersistentTetFacetBinding{
                key,triple,{2,79},*face});
    }
    for(std::size_t i=0;i<handles.size();++i)
        fixture.localByGlobal.emplace(static_cast<Index>(i),handles[i]);
    return fixture;
}

std::vector<bnd::VolumeArtifactEdgeKey> expectedEdges(
    const Fixture& fixture) {
    std::vector<bnd::VolumeArtifactEdgeKey> result;
    result.reserve(fixture.artifact.edges.size());
    for(const auto& row:fixture.artifact.edges)result.push_back(row.key);
    std::sort(result.begin(),result.end());
    return result;
}

std::vector<bnd::VolumeFacetKey> expectedFacets(const Fixture& fixture) {
    std::vector<const bnd::VolumeFacet*> sorted;
    sorted.reserve(fixture.artifact.facets.size());
    for(const auto& row:fixture.artifact.facets)sorted.push_back(&row);
    std::sort(sorted.begin(),sorted.end(),[](const auto* a,const auto* b) {
        const int amin=a->owners.empty() ? -1:a->owners.front().regionTag;
        const int bmin=b->owners.empty() ? -1:b->owners.front().regionTag;
        const int amax=a->owners.empty() ? -1:a->owners.back().regionTag;
        const int bmax=b->owners.empty() ? -1:b->owners.back().regionTag;
        const bool ai=a->owners.size()==2;
        const bool bi=b->owners.size()==2;
        return std::tie(amin,ai,amax,a->key,a->directedGlobalNodes)<
               std::tie(bmin,bi,bmax,b->key,b->directedGlobalNodes);
    });
    std::vector<bnd::VolumeFacetKey> result;
    result.reserve(sorted.size());
    for(const auto* row:sorted)result.push_back(row->key);
    return result;
}

bool sameHandle(topo::EdgeHandle a,topo::EdgeHandle b) {
    return a==b;
}
bool sameHandle(topo::FaceHandle a,topo::FaceHandle b) {
    return a==b;
}

bool sameState(const bnd::TetFacetSetState& a,
               const bnd::TetFacetSetState& b) {
    if(a.workspaceStamp!=b.workspaceStamp||a.artifactStamp!=b.artifactStamp||
       a.canonicalEdgeSchedule!=b.canonicalEdgeSchedule||
       a.canonicalSchedule!=b.canonicalSchedule||
       a.edgeBindings.size()!=b.edgeBindings.size()||
       a.bindings.size()!=b.bindings.size())return false;
    auto ae=a.edgeBindings.begin();
    auto be=b.edgeBindings.begin();
    for(;ae!=a.edgeBindings.end();++ae,++be)
        if(ae->first!=be->first||ae->second.key!=be->second.key||
           ae->second.classification.dim!=be->second.classification.dim||
           ae->second.classification.tag!=be->second.classification.tag||
           ae->second.required!=be->second.required||
           !sameHandle(ae->second.edge,be->second.edge))return false;
    auto af=a.bindings.begin();
    auto bf=b.bindings.begin();
    for(;af!=a.bindings.end();++af,++bf)
        if(af->first!=bf->first||af->second.key!=bf->second.key||
           af->second.canonicalVertices!=bf->second.canonicalVertices||
           af->second.classification.dim!=bf->second.classification.dim||
           af->second.classification.tag!=bf->second.classification.tag||
           !sameHandle(af->second.recoveredFace,
                       bf->second.recoveredFace))return false;
    return true;
}

struct CleanReference {
    std::size_t allocationCalls=0;
    std::size_t requestedBytes=0;
    bnd::TetFacetSetMixedPresentMetrics metrics{};
    std::size_t processedRows=0;
};

bnd::TetFacetSetRecoveryOutcome run(
    Fixture& fixture,bnd::TetFacetSetMixedPresentBatchSink& sink) {
    return bnd::recoverTetFacetSetIndexed(
        fixture.artifact,fixture.localByGlobal,{3,9},fixture.state,
        *fixture.workspace,*fixture.context,nullptr,nullptr,&sink);
}

void requireClean(BatchKind kind,const Fixture& fixture,
                  const bnd::TetFacetSetRecoveryOutcome& outcome,
                  const bnd::TetFacetSetMixedPresentBatchSink& sink,
                  const CleanReference* reference=nullptr) {
    if(outcome.status!=bnd::TetFacetSetRecoveryStatus::Complete)
        std::printf("  clean %s status=%d reason=%d processed=%zu "
                    "unattempted=%zu/%zu\n",
            kind==BatchKind::Edge ? "edge":"facet",
            static_cast<int>(outcome.status),
            static_cast<int>(outcome.reason),outcome.processed.size(),
            outcome.unattemptedEdges.size(),
            outcome.unattemptedFacets.size());
    REQUIRE(outcome.status==bnd::TetFacetSetRecoveryStatus::Complete);
    REQUIRE(outcome.reason==bnd::TetFacetSetRecoveryReason::None);
    REQUIRE(outcome.reportComplete);
    REQUIRE(outcome.unattemptedEdges.empty());
    REQUIRE(outcome.unattemptedFacets.empty());
    const auto expected=kind==BatchKind::Edge
        ? fixture.artifact.edges.size():fixture.artifact.facets.size();
    REQUIRE(outcome.processed.size()==expected);
    REQUIRE(fixture.state.canonicalEdgeSchedule==expectedEdges(fixture));
    REQUIRE(fixture.state.canonicalSchedule==expectedFacets(fixture));
    REQUIRE(sink.metrics.complete);
    REQUIRE(sink.metrics.presentBindingBatches==1);
    REQUIRE(sink.metrics.presentBindingEntities==expected);
    REQUIRE(sink.metrics.stateCopies==1);
    REQUIRE(sink.metrics.nontrivialRecoveries==0);
    REQUIRE(sink.metrics.modeledBatchPreparationAllocationRequests>0);
    REQUIRE(sink.metrics.modeledBatchPreparationLogicalBytes>0);
    const auto& assessment=sink.metrics.assessmentBindingWork;
    const auto& recovery=sink.metrics.recoveryBindingWork;
    REQUIRE(assessment.complete);
    REQUIRE(recovery.complete);
    REQUIRE(assessment.edgeDataCopyRows>=
            assessment.edgeBindingRowsVisited);
    REQUIRE(assessment.faceDataCopyRows>=
            assessment.facetBindingRowsVisited);
    REQUIRE(assessment.faceIncidenceRowsCopied>=
            assessment.faceDataCopyRows);
    REQUIRE(recovery.stateCopies==0);
    REQUIRE(recovery.stateCopyRows==0);
    if(kind==BatchKind::Edge) {
        REQUIRE(fixture.state.edgeBindings.size()==expected);
        REQUIRE(fixture.state.bindings.empty());
        for(const auto& row:outcome.processed) {
            REQUIRE(row.kind==bnd::TetFacetSetWorkKind::Edge);
            REQUIRE(row.edgeOutcome.status==
                    bnd::TetEdgeRecoveryStatus::AlreadyPresent);
        }
    } else {
        // The retained entry binding is validated before schedule setup, the
        // one facet assessment reuses that receipt, and the prepared prefix
        // validates its scratch state as the second receipt build.
        REQUIRE(fixture.artifact.facets.size()==1);
        REQUIRE(assessment.facetBindingReceiptBuilds==2);
        REQUIRE(assessment.facetBindingReceiptReuses==1);
        REQUIRE(assessment.facetBindingReceiptReserveCalls==1);
        REQUIRE(assessment.facetBindingReceiptRowsCopied==2);
        REQUIRE(assessment.facetBindingReceiptLookups>0);
        REQUIRE(assessment.facetBindingRequestRowsVisited==0);
        REQUIRE(fixture.state.bindings.size()==expected);
        REQUIRE(fixture.state.edgeBindings.empty());
        for(const auto& row:outcome.processed) {
            REQUIRE(row.kind==bnd::TetFacetSetWorkKind::Facet);
            REQUIRE(row.facetOutcome.status==
                    bnd::TetFacetRecoveryStatus::AlreadyPresent);
        }
    }
    REQUIRE(fixture.context->currentFor(*fixture.workspace));
    if(reference) {
        REQUIRE(outcome.processed.size()==reference->processedRows);
        REQUIRE(sink.metrics==reference->metrics);
    }
}

CleanReference cleanReference(BatchKind kind) {
    auto fixture=makeFixture(kind);
    const auto before=paacs_tests::observableDigest(*fixture.workspace);
    bnd::TetFacetSetMixedPresentBatchSink sink;
    gCalls=0; gBytes=0; gFailAt=0; gThrown=0; gPhase3Allocations=0;
    gMixedCommitAllocations=0;
    const auto outcome=run(fixture,sink);
    CleanReference reference;
    reference.allocationCalls=gCalls;
    reference.requestedBytes=gBytes;
    reference.metrics=sink.metrics;
    reference.processedRows=outcome.processed.size();
    requireClean(kind,fixture,outcome,sink);
    REQUIRE(gPhase3Allocations==0);
    REQUIRE(gMixedCommitAllocations==0);
    REQUIRE(paacs_tests::observableDigest(*fixture.workspace)==before);
    return reference;
}

struct SweepResult {
    std::size_t faults=0;
    std::size_t envelopeFaults=0;
    std::size_t earlyReceiptFaults=0;
};

SweepResult sweep(BatchKind kind,const CleanReference& reference) {
    SweepResult swept;
    for(std::size_t ordinal=1;ordinal<=reference.allocationCalls;++ordinal) {
        gCurrentOrdinal=ordinal;
        gCurrentKind=kind==BatchKind::Edge ? 0:1;
        auto fixture=makeFixture(kind);
        const auto workspaceBefore=
            paacs_tests::observableDigest(*fixture.workspace);
        const bnd::TetFacetSetState stateBefore=fixture.state;
        bnd::TetFacetSetMixedPresentBatchSink sink;
        gCalls=0; gBytes=0; gFailAt=ordinal; gThrown=0;
        gPhase3Allocations=0;
        gMixedCommitAllocations=0;
        const auto outcome=run(fixture,sink);
        gFailAt=0;
        REQUIRE(gThrown==1);
        REQUIRE(gCalls==ordinal);
        REQUIRE(outcome.status==bnd::TetFacetSetRecoveryStatus::Failed);
        REQUIRE(outcome.reason==
                bnd::TetFacetSetRecoveryReason::AllocationFailed);
        if(kind==BatchKind::Facet&&ordinal==1) {
            ++swept.earlyReceiptFaults;
            REQUIRE(!outcome.reportComplete);
            REQUIRE(sink.metrics.assessmentBindingWork.
                    facetBindingAuditCalls==1);
            REQUIRE(sink.metrics.assessmentBindingWork.
                    facetBindingReceiptBuilds==1);
            REQUIRE(sink.metrics.assessmentBindingWork.
                    facetBindingReceiptReserveCalls==1);
            REQUIRE(sink.metrics.assessmentBindingWork.
                    facetBindingReceiptRowsCopied==0);
            REQUIRE(fixture.state.canonicalEdgeSchedule.empty());
            REQUIRE(fixture.state.canonicalSchedule.empty());
        }
        if(outcome.reportComplete) {
            REQUIRE(outcome.processed.empty());
            if(kind==BatchKind::Edge) {
                REQUIRE(outcome.unattemptedEdges==
                        expectedEdges(fixture));
                REQUIRE(outcome.unattemptedFacets.empty());
            } else {
                REQUIRE(outcome.unattemptedEdges.empty());
                REQUIRE(outcome.unattemptedFacets==
                        expectedFacets(fixture));
            }
        } else {
            REQUIRE(outcome.processed.empty());
            REQUIRE(outcome.unattemptedEdges.empty());
            REQUIRE(outcome.unattemptedFacets.empty());
            ++swept.envelopeFaults;
        }
        REQUIRE(sameState(fixture.state,stateBefore));
        REQUIRE(paacs_tests::observableDigest(*fixture.workspace)==
                workspaceBefore);
        REQUIRE(fixture.context->currentFor(*fixture.workspace));
        REQUIRE(gPhase3Allocations==0);
        REQUIRE(gMixedCommitAllocations==0);

        bnd::TetFacetSetMixedPresentBatchSink retrySink;
        gCalls=0; gBytes=0; gThrown=0; gPhase3Allocations=0;
        gMixedCommitAllocations=0;
        const auto retry=run(fixture,retrySink);
        const auto retryCalls=gCalls;
        const auto retryBytes=gBytes;
        requireClean(kind,fixture,retry,retrySink,&reference);
        REQUIRE(retryCalls==reference.allocationCalls);
        REQUIRE(retryBytes==reference.requestedBytes);
        REQUIRE(gPhase3Allocations==0);
        REQUIRE(gMixedCommitAllocations==0);
        REQUIRE(paacs_tests::observableDigest(*fixture.workspace)==
                workspaceBefore);
        ++swept.faults;
        if(failures!=0)break;
    }
    return swept;
}

Fixture makeMixedFixture() {
    Fixture fixture;
    topo::NodePool nodes;
    std::vector<topo::NodeHandle> handles;
    for(const Vec3& p:std::array<Vec3,11>{
            Vec3{-1,-1,0},Vec3{1,-1,0},Vec3{1,1,0},Vec3{-1,1,0},
            Vec3{0,0,2},Vec3{0,0,-2},
            Vec3{10,0,0},Vec3{11,0,0},Vec3{10,1,0},
            Vec3{10,0,1},Vec3{10,0,-1}})
        handles.push_back(nodes.add(p,{3,9}));
    topo::CellPool cells;
    for(std::size_t i=0;i<4;++i)
        (void)cells.add(topo::ElementType::Tet4,
            positive(nodes,{handles[4],handles[5],handles[i],
                            handles[(i+1)%4]}),{3,9},nodes);
    (void)cells.add(topo::ElementType::Tet4,
        positive(nodes,{handles[6],handles[7],handles[8],handles[9]}),
        {3,9},nodes);
    (void)cells.add(topo::ElementType::Tet4,
        positive(nodes,{handles[6],handles[8],handles[7],handles[10]}),
        {3,9},nodes);
    auto edges=topo::EdgePool::buildComplete(nodes,cells);
    for(std::size_t i=0;i<4;++i) {
        const auto edge=edges.find(nodes,handles[i],handles[(i+1)%4]);
        if(!edge)throw std::logic_error("missing mixed ring edge");
        edges.setClassification(*edge,{1,41});
        edges.setRequiredConstraint(*edge,true);
    }
    for(const auto ids:std::array<std::array<std::size_t,2>,3>{
            {{{6,7}},{{7,8}},{{8,6}}}}) {
        const auto edge=edges.find(nodes,handles[ids[0]],handles[ids[1]]);
        if(!edge)throw std::logic_error("missing mixed face edge");
        edges.setClassification(*edge,{2,79});
        edges.setRequiredConstraint(*edge,true);
    }
    auto faces=topo::FacePool::buildComplete(nodes,cells,edges);
    const auto protectedFace=faces.find(nodes,
        std::array<topo::NodeHandle,3>{handles[6],handles[7],handles[8]});
    if(!protectedFace)throw std::logic_error("missing mixed face");
    faces.setClassification(*protectedFace,{2,79});
    fixture.workspace=std::make_unique<topo::TopologyWorkspace>(
        std::move(nodes),std::move(cells),std::move(edges),std::move(faces));
    auto indexed=bnd::buildTetRecoveryContext(*fixture.workspace,{3,9});
    if(!indexed.built()||!indexed.context)
        throw std::logic_error("mixed-present recovery context");
    fixture.context=std::move(indexed.context);
    fixture.state.workspaceStamp=topo::WorkspaceStamp::of(*fixture.workspace);
    fixture.state.artifactStamp={704,1};
    fixture.artifact.surfaceSetStamp=fixture.state.artifactStamp;
    fixture.artifact.regionTags={9};
    fixture.artifact.facets={
        {{809,0},{6,7,8},{2,79},{{9,1}}},
        {{810,0},{0,1,2},{2,77},{{9,1}}},
        {{811,0},{0,2,3},{2,78},{{9,1}}}};
    std::map<bnd::VolumeArtifactEdgeKey,bnd::VolumeArtifactEdge> edgeRows;
    for(const auto& facet:fixture.artifact.facets)
        for(std::size_t i=0;i<3;++i) {
            const auto key=bnd::VolumeArtifactEdgeKey::of(
                facet.directedGlobalNodes[i],
                facet.directedGlobalNodes[(i+1)%3]);
            const auto classification=key.canonicalGlobalNodes[0]>=6
                ? mesh::Classification{2,79}:mesh::Classification{1,41};
            auto [it,inserted]=edgeRows.emplace(
                key,bnd::VolumeArtifactEdge{key,classification,{}});
            if(!inserted&&
               (it->second.classification.dim!=classification.dim||
                it->second.classification.tag!=classification.tag))
                throw std::logic_error("mixed edge classification");
            it->second.incidentFacets.push_back(facet.key);
        }
    for(auto& [key,row]:edgeRows) {
        (void)key;
        std::sort(row.incidentFacets.begin(),row.incidentFacets.end());
        fixture.artifact.edges.push_back(std::move(row));
    }
    for(std::size_t i=0;i<handles.size();++i)
        fixture.localByGlobal.emplace(static_cast<Index>(i),handles[i]);
    return fixture;
}

bool edgeSucceeded(bnd::TetEdgeRecoveryStatus status) {
    return status==bnd::TetEdgeRecoveryStatus::AlreadyPresent||
           status==bnd::TetEdgeRecoveryStatus::Marked||
           status==bnd::TetEdgeRecoveryStatus::Recovered;
}

bool facetSucceeded(bnd::TetFacetRecoveryStatus status) {
    return status==bnd::TetFacetRecoveryStatus::AlreadyPresent||
           status==bnd::TetFacetRecoveryStatus::Marked||
           status==bnd::TetFacetRecoveryStatus::Recovered;
}

void requireCompleteXor(const Fixture& fixture,
                        const bnd::TetFacetSetRecoveryOutcome& outcome) {
    const auto edgeKeys=expectedEdges(fixture);
    const auto facetKeys=expectedFacets(fixture);
    for(const auto& key:edgeKeys) {
        std::size_t processed=0;
        std::size_t unattempted=0;
        for(const auto& row:outcome.processed)
            if(row.kind==bnd::TetFacetSetWorkKind::Edge&&row.edgeKey==key)
                ++processed;
        for(const auto& row:outcome.unattemptedEdges)
            if(row==key)++unattempted;
        REQUIRE(processed+unattempted==1);
        REQUIRE(processed==0||unattempted==0);
    }
    for(const auto& key:facetKeys) {
        std::size_t processed=0;
        std::size_t unattempted=0;
        for(const auto& row:outcome.processed)
            if(row.kind==bnd::TetFacetSetWorkKind::Facet&&row.facetKey==key)
                ++processed;
        for(const auto& row:outcome.unattemptedFacets)
            if(row==key)++unattempted;
        REQUIRE(processed+unattempted==1);
        REQUIRE(processed==0||unattempted==0);
    }
    REQUIRE(outcome.processed.size()+outcome.unattemptedEdges.size()+
            outcome.unattemptedFacets.size()==edgeKeys.size()+facetKeys.size());
}

void requireStateMatchesProcessed(
    const Fixture& fixture,const bnd::TetFacetSetRecoveryOutcome& outcome,
    bool allowTerminalCommittedBinding=false) {
    std::size_t successfulEdges=0;
    std::size_t successfulFacets=0;
    for(const auto& row:outcome.processed) {
        if(row.kind==bnd::TetFacetSetWorkKind::Edge&&
           edgeSucceeded(row.edgeOutcome.status)) {
            ++successfulEdges;
            REQUIRE(fixture.state.edgeBindings.find(row.edgeKey)!=
                    fixture.state.edgeBindings.end());
        }
        if(row.kind==bnd::TetFacetSetWorkKind::Facet&&
           facetSucceeded(row.facetOutcome.status)) {
            ++successfulFacets;
            REQUIRE(fixture.state.bindings.find(row.facetKey)!=
                    fixture.state.bindings.end());
        }
    }
    if(allowTerminalCommittedBinding&&!outcome.processed.empty()) {
        const auto& row=outcome.processed.back();
        if(row.kind==bnd::TetFacetSetWorkKind::Edge&&
           !edgeSucceeded(row.edgeOutcome.status)&&
           fixture.state.edgeBindings.find(row.edgeKey)!=
               fixture.state.edgeBindings.end())
            ++successfulEdges;
        if(row.kind==bnd::TetFacetSetWorkKind::Facet&&
           !facetSucceeded(row.facetOutcome.status)&&
           fixture.state.bindings.find(row.facetKey)!=
               fixture.state.bindings.end())
            ++successfulFacets;
    }
    if(fixture.state.edgeBindings.size()!=successfulEdges||
       fixture.state.bindings.size()!=successfulFacets) {
        std::printf("  mixed state/report mismatch ordinal=%zu "
                    "edge=%zu/%zu facet=%zu/%zu processed=%zu\n",
            gCurrentOrdinal,fixture.state.edgeBindings.size(),successfulEdges,
            fixture.state.bindings.size(),successfulFacets,
            outcome.processed.size());
        for(const auto& row:outcome.processed)
            if(row.kind==bnd::TetFacetSetWorkKind::Facet)
                std::printf("    facet key=%llu/%zu status=%d reason=%d "
                            "phase=%d\n",
                    static_cast<unsigned long long>(
                        row.facetKey.facePersistentId),
                    row.facetKey.normalizedLocalCellOrdinal,
                    static_cast<int>(row.facetOutcome.status),
                    static_cast<int>(row.facetOutcome.reason),
                    static_cast<int>(row.facetOutcome.phase));
        for(const auto& [key,binding]:fixture.state.bindings) {
            (void)binding;
            std::printf("    state facet key=%llu/%zu\n",
                static_cast<unsigned long long>(key.facePersistentId),
                key.normalizedLocalCellOrdinal);
        }
    }
    REQUIRE(fixture.state.edgeBindings.size()==successfulEdges);
    REQUIRE(fixture.state.bindings.size()==successfulFacets);
}

bool hasTerminalAllocationRow(
    const bnd::TetFacetSetRecoveryOutcome& outcome) {
    if(outcome.processed.empty())return false;
    const auto& row=outcome.processed.back();
    if(row.kind==bnd::TetFacetSetWorkKind::Edge)
        return row.edgeOutcome.status==bnd::TetEdgeRecoveryStatus::Failed&&
               row.edgeOutcome.reason==
                   bnd::TetEdgeRecoveryReason::AllocationFailed;
    return row.facetOutcome.status==bnd::TetFacetRecoveryStatus::Failed&&
           row.facetOutcome.reason==
               bnd::TetFacetRecoveryReason::AllocationFailed;
}

std::uint64_t recoveredRows(
    const bnd::TetFacetSetRecoveryOutcome& outcome) {
    std::uint64_t count=0;
    for(const auto& row:outcome.processed) {
        if(row.kind==bnd::TetFacetSetWorkKind::Edge) {
            if(row.edgeOutcome.status==bnd::TetEdgeRecoveryStatus::Recovered)
                ++count;
        } else if(row.facetOutcome.status==
                  bnd::TetFacetRecoveryStatus::Recovered) {
            ++count;
        }
    }
    return count;
}

struct MixedReference {
    std::size_t allocationCalls=0;
    std::size_t requestedBytes=0;
    std::uint64_t semanticHash=0;
    bnd::TetFacetSetMixedPresentMetrics metrics{};
};

MixedReference mixedReference() {
    auto fixture=makeMixedFixture();
    bnd::TetFacetSetMixedPresentBatchSink sink;
    gCalls=0; gBytes=0; gFailAt=0; gThrown=0;
    gPhase3Allocations=0; gMixedCommitAllocations=0;
    const auto outcome=run(fixture,sink);
    const auto routeCalls=gCalls;
    const auto routeBytes=gBytes;
    REQUIRE(outcome.status==bnd::TetFacetSetRecoveryStatus::Complete);
    REQUIRE(outcome.reason==bnd::TetFacetSetRecoveryReason::None);
    REQUIRE(outcome.reportComplete);
    requireCompleteXor(fixture,outcome);
    requireStateMatchesProcessed(fixture,outcome);
    REQUIRE(outcome.unattemptedEdges.empty());
    REQUIRE(outcome.unattemptedFacets.empty());
    REQUIRE(sink.metrics.presentBindingBatches>=2);
    REQUIRE(sink.metrics.presentBindingEntities>0);
    REQUIRE(sink.metrics.nontrivialRecoveries>0);
    REQUIRE(sink.metrics.nontrivialRecoveries==recoveredRows(outcome));
    REQUIRE(sink.metrics.assessmentBindingWork.complete);
    REQUIRE(sink.metrics.recoveryBindingWork.complete);
    REQUIRE(sink.metrics.assessmentBindingWork.edgeBindingRowsVisited>0);
    REQUIRE(sink.metrics.assessmentBindingWork.facetBindingRowsVisited>0);
    REQUIRE(sink.metrics.assessmentBindingWork.faceIncidenceRowsCopied>=
            sink.metrics.assessmentBindingWork.faceDataCopyRows);
    REQUIRE(sink.metrics.assessmentBindingWork.
            facetBindingReceiptBuilds>0);
    REQUIRE(sink.metrics.assessmentBindingWork.
            facetBindingReceiptReuses>0);
    REQUIRE(sink.metrics.assessmentBindingWork.
            facetBindingRequestRowsVisited==0);
    REQUIRE(sink.metrics.recoveryBindingWork.stateCopies>0);
    REQUIRE(sink.metrics.recoveryBindingWork.stateCopyRows>=
            sink.metrics.recoveryBindingWork.stateCopies);
    REQUIRE(gMixedCommitAllocations==0);
    REQUIRE(gPhase3Allocations==0);
    MixedReference reference;
    reference.allocationCalls=routeCalls;
    reference.requestedBytes=routeBytes;
    reference.semanticHash=fixture.workspace->semanticTopologyHashV2();
    reference.metrics=sink.metrics;
    return reference;
}

struct MixedSweepResult {
    std::size_t faults=0;
    std::size_t postCommitFaults=0;
    std::size_t staleUnpoisonedFaults=0;
    std::size_t postFaultAllocationRuns=0;
    std::size_t postFaultAllocations=0;
};

MixedSweepResult sweepMixed(const MixedReference& reference) {
    MixedSweepResult swept;
    for(std::size_t ordinal=1;ordinal<=reference.allocationCalls;++ordinal) {
        gCurrentOrdinal=ordinal;
        gCurrentKind=2;
        auto fixture=makeMixedFixture();
        const auto before=paacs_tests::observableDigest(*fixture.workspace);
        const auto revision=fixture.workspace->topologyRevision();
        const auto initialState=fixture.state;
        bnd::TetFacetSetMixedPresentBatchSink sink;
        gCalls=0; gBytes=0; gFailAt=ordinal; gThrown=0;
        gPhase3Allocations=0; gMixedCommitAllocations=0;
        const auto outcome=run(fixture,sink);
        gFailAt=0;
        REQUIRE(gThrown==1);
        if(gCalls!=ordinal&&swept.postFaultAllocationRuns==0)
            std::printf("  mixed post-fault allocation ordinal=%zu "
                        "calls=%zu processed=%zu outer=%d/%d\n",
                ordinal,gCalls,outcome.processed.size(),
                static_cast<int>(outcome.status),
                static_cast<int>(outcome.reason));
        if(gCalls>ordinal) {
            ++swept.postFaultAllocationRuns;
            swept.postFaultAllocations+=gCalls-ordinal;
        }
        REQUIRE(outcome.status==bnd::TetFacetSetRecoveryStatus::Failed);
        const bool allocationFailure=
            outcome.reason==bnd::TetFacetSetRecoveryReason::AllocationFailed||
            (outcome.reason==bnd::TetFacetSetRecoveryReason::None&&
             hasTerminalAllocationRow(outcome));
        if(!allocationFailure)
            std::printf("  mixed ordinal=%zu status=%d reason=%d "
                        "processed=%zu unattempted=%zu/%zu calls=%zu\n",
                ordinal,static_cast<int>(outcome.status),
                static_cast<int>(outcome.reason),outcome.processed.size(),
                outcome.unattemptedEdges.size(),
                outcome.unattemptedFacets.size(),gCalls);
        REQUIRE(allocationFailure);
        REQUIRE(sink.metrics.nontrivialRecoveries==recoveredRows(outcome));
        REQUIRE(gMixedCommitAllocations==0);
        REQUIRE(gPhase3Allocations==0);
        const bool contextCurrent=
            fixture.context->currentFor(*fixture.workspace);
        if(outcome.reportComplete) {
            requireCompleteXor(fixture,outcome);
            requireStateMatchesProcessed(
                fixture,outcome,!contextCurrent);
            if(outcome.processed.empty()) {
                REQUIRE(sameState(fixture.state,initialState));
                REQUIRE(paacs_tests::observableDigest(*fixture.workspace)==
                        before);
            }
        } else {
            REQUIRE(outcome.processed.empty());
            REQUIRE(outcome.unattemptedEdges.empty());
            REQUIRE(outcome.unattemptedFacets.empty());
            REQUIRE(sameState(fixture.state,initialState));
            REQUIRE(paacs_tests::observableDigest(*fixture.workspace)==before);
        }
        if(!contextCurrent&&swept.postCommitFaults==0) {
            std::printf("  mixed poisoned ordinal=%zu processed=%zu "
                        "outer=%d/%d phase3=%zu commit=%zu",
                ordinal,outcome.processed.size(),
                static_cast<int>(outcome.status),
                static_cast<int>(outcome.reason),gPhase3Allocations,
                gMixedCommitAllocations);
            if(!outcome.processed.empty()) {
                const auto& row=outcome.processed.back();
                if(row.kind==bnd::TetFacetSetWorkKind::Edge)
                    std::printf(" edge=%d/%d/%d",
                        static_cast<int>(row.edgeOutcome.status),
                        static_cast<int>(row.edgeOutcome.reason),
                        static_cast<int>(row.edgeOutcome.phase));
                else
                    std::printf(" facet=%d/%d/%d",
                        static_cast<int>(row.facetOutcome.status),
                        static_cast<int>(row.facetOutcome.reason),
                        static_cast<int>(row.facetOutcome.phase));
            }
            std::printf("\n");
        }
        if(contextCurrent) {
            REQUIRE(fixture.workspace->topologyRevision()==revision||
                    !outcome.processed.empty());
        } else {
            ++swept.postCommitFaults;
            if(!fixture.context->poisoned) {
                ++swept.staleUnpoisonedFaults;
                if(swept.staleUnpoisonedFaults==1)
                  std::printf("  mixed stale-unpoisoned ordinal=%zu "
                            "processed=%zu revisions=%llu/%llu/%llu "
                            "state=%llu",ordinal,outcome.processed.size(),
                            static_cast<unsigned long long>(revision),
                            static_cast<unsigned long long>(
                                fixture.context->contextRevision),
                            static_cast<unsigned long long>(
                                fixture.workspace->topologyRevision()),
                            static_cast<unsigned long long>(
                                fixture.state.workspaceStamp.revision));
                if(swept.staleUnpoisonedFaults==1&&
                   !outcome.processed.empty()) {
                    const auto& row=outcome.processed.back();
                    if(row.kind==bnd::TetFacetSetWorkKind::Edge)
                        std::printf(" edge=%d/%d/%d",
                            static_cast<int>(row.edgeOutcome.status),
                            static_cast<int>(row.edgeOutcome.reason),
                            static_cast<int>(row.edgeOutcome.phase));
                    else
                        std::printf(" facet=%d/%d/%d",
                            static_cast<int>(row.facetOutcome.status),
                            static_cast<int>(row.facetOutcome.reason),
                            static_cast<int>(row.facetOutcome.phase));
                }
                if(swept.staleUnpoisonedFaults==1)std::printf("\n");
            }
            REQUIRE(fixture.workspace->topologyRevision()!=revision);
            REQUIRE(topo::workspaceTopologyReady(
                fixture.workspace->nodes(),fixture.workspace->cells(),
                fixture.workspace->edges(),fixture.workspace->faces()));
        }

        bnd::TetFacetSetMixedPresentBatchSink retrySink;
        if(!contextCurrent)fixture=makeMixedFixture();
        gCalls=0; gBytes=0; gThrown=0;
        gPhase3Allocations=0; gMixedCommitAllocations=0;
        const auto retry=run(fixture,retrySink);
        REQUIRE(retry.status==bnd::TetFacetSetRecoveryStatus::Complete);
        REQUIRE(retry.reason==bnd::TetFacetSetRecoveryReason::None);
        REQUIRE(retry.reportComplete);
        requireCompleteXor(fixture,retry);
        requireStateMatchesProcessed(fixture,retry);
        REQUIRE(fixture.state.canonicalEdgeSchedule==expectedEdges(fixture));
        REQUIRE(fixture.state.canonicalSchedule==expectedFacets(fixture));
        REQUIRE(fixture.workspace->semanticTopologyHashV2()==
                reference.semanticHash);
        // A same-fixture retry may begin after a successfully published
        // canonical-prefix binding batch. Its remaining work is therefore
        // intentionally smaller than the fresh-fixture reference, while the
        // completed report/state/hash must still reconcile exactly.
        REQUIRE(retrySink.metrics.complete);
        REQUIRE(retrySink.metrics.assessmentBindingWork.complete);
        REQUIRE(retrySink.metrics.recoveryBindingWork.complete);
        REQUIRE(retrySink.metrics.nontrivialRecoveries==
                recoveredRows(retry));
        REQUIRE(retrySink.metrics.exactPresentSelectedRebuilds==0);
        REQUIRE(retrySink.metrics.exactPresentDependencyTests==0);
        REQUIRE(gMixedCommitAllocations==0);
        REQUIRE(gPhase3Allocations==0);
        ++swept.faults;
        if(failures!=0)break;
    }
    return swept;
}

} // namespace

int main() {
    std::set_terminate([] {
        std::fprintf(stderr,"TERMINATE kind=%d ordinal=%zu calls=%zu "
                    "fail=%zu thrown=%zu uncaught=%d current=%d\n",
                    gCurrentKind,gCurrentOrdinal,gCalls,gFailAt,gThrown,
                    std::uncaught_exceptions(),std::current_exception()?1:0);
        std::fflush(stderr);
        std::_Exit(99);
    });
    try {
    const auto edge=cleanReference(BatchKind::Edge);
    const auto facet=cleanReference(BatchKind::Facet);
    // Exclude one-time CRT/STL initialization in the first genuinely
    // mutating mixed route, matching the canonical warmed fault probes.
    (void)mixedReference();
    const auto mixed=mixedReference();
    std::printf("  clean edge=%zu facet=%zu mixed=%zu\n",
                edge.allocationCalls,facet.allocationCalls,
                mixed.allocationCalls);
    std::fflush(stdout);
    const auto edgeFaults=sweep(BatchKind::Edge,edge);
    const auto facetFaults=sweep(BatchKind::Facet,facet);
    const auto mixedFaults=sweepMixed(mixed);
    std::printf("  swept edge=%zu/%zu prep_model=%llu facet=%zu/%zu "
                "prep_model=%llu early_receipt=%zu mixed=%zu/%zu "
                "post_commit=%zu "
                "stale_unpoisoned=%zu "
                "post_fault_alloc=%zu/%zu\n",
        edgeFaults.faults,edge.allocationCalls,
        static_cast<unsigned long long>(
            edge.metrics.modeledBatchPreparationAllocationRequests),
        facetFaults.faults,facet.allocationCalls,
        static_cast<unsigned long long>(
            facet.metrics.modeledBatchPreparationAllocationRequests),
        facetFaults.earlyReceiptFaults,
        mixedFaults.faults,mixed.allocationCalls,
        mixedFaults.postCommitFaults,mixedFaults.staleUnpoisonedFaults,
        mixedFaults.postFaultAllocationRuns,
        mixedFaults.postFaultAllocations);
    REQUIRE(edgeFaults.faults==edge.allocationCalls);
    REQUIRE(facetFaults.faults==facet.allocationCalls);
    REQUIRE(facetFaults.earlyReceiptFaults==1);
    REQUIRE(mixedFaults.faults==mixed.allocationCalls);
    REQUIRE(mixedFaults.staleUnpoisonedFaults==0);
    if(failures==0) {
        std::printf(
            "probe_t2b_mixed_present PASSED: edge=%zu/%zu prep=%zu "
            "envelope=%zu (%zu bytes), facet=%zu/%zu prep=%zu "
            "envelope=%zu early_receipt=%zu (%zu bytes), "
            "mixed=%zu/%zu post_commit=%zu "
            "stale_unpoisoned=%zu "
            "post_fault_alloc=%zu/%zu "
            "(%zu bytes), "
            "zero mixed commit-tail allocations\n",
            edgeFaults.faults,edge.allocationCalls,
            static_cast<std::size_t>(
                edge.metrics.modeledBatchPreparationAllocationRequests),
            edgeFaults.envelopeFaults,
            edge.requestedBytes,facetFaults.faults,facet.allocationCalls,
            static_cast<std::size_t>(
            facet.metrics.modeledBatchPreparationAllocationRequests),
            facetFaults.envelopeFaults,facetFaults.earlyReceiptFaults,
            facet.requestedBytes,
            mixedFaults.faults,mixed.allocationCalls,
            mixedFaults.postCommitFaults,
            mixedFaults.staleUnpoisonedFaults,
            mixedFaults.postFaultAllocationRuns,
            mixedFaults.postFaultAllocations,mixed.requestedBytes);
    }
    return failures==0 ? 0:1;
    } catch(const std::exception& error) {
        std::printf("probe_t2b_mixed_present exception: %s\n",error.what());
    } catch(...) {
        std::printf("probe_t2b_mixed_present unknown exception\n");
    }
    return 2;
}
