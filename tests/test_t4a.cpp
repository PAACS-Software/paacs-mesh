// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "test_harness.hpp"

#include "paacs/algo/TetQualityMesher.hpp"
#include "paacs/algo/Meshers.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/sizing/SizeField.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#if defined(_WIN32)
#include <process.h>
#define getpid _getpid
#else
#include <unistd.h>
#endif

namespace {
using namespace paacs;

class Provider final : public algo::VolumeBoundaryOwnershipProvider {
public:
    algo::VolumeBoundaryProviderStamp current{1701,1};
    std::map<algo::FacePersistentId,std::vector<bnd::VolumeRegionUse>> rows;
    mutable std::uint64_t ownershipCalls=0,envelopeCalls=0;

    algo::VolumeBoundaryProviderStamp stamp() const noexcept override {
        return current;
    }
    algo::VolumeFaceOwnership ownership(
            algo::FacePersistentId face) const override {
        algo::VolumeFaceOwnership result;
        ++ownershipCalls;
        result.providerStamp=current;
        result.facePersistentId=face;
        const auto found=rows.find(face);
        if (found!=rows.end()) result.owners=found->second;
        return result;
    }
    algo::ProviderRecoveryEnvelope recoveryEnvelope(
            std::span<const Vec3>) const override {
        algo::ProviderRecoveryEnvelope result;
        ++envelopeCalls;
        result.providerStamp=current;
        result.availability=
            algo::VolumeProviderAvailability::RecoveryEnvelopeUnavailable;
        return result;
    }
};

class CountingField final : public sizing::SizeField {
public:
    explicit CountingField(double value) : value_(value) {}
    double eval(const Vec3&) const override {
        ++calls;
        return value_;
    }
    mutable std::uint64_t calls=0;
private:
    double value_=1.0;
};

class ReversedFaceModel final : public geo::Model {
public:
    explicit ReversedFaceModel(const geo::Model& source) : source_(source) {}
    std::vector<const geo::GVertex*> vertices() const override {
        return source_.vertices();
    }
    std::vector<const geo::GEdge*> edges() const override {
        return source_.edges();
    }
    std::vector<const geo::GFace*> faces() const override {
        auto result=source_.faces();
        std::reverse(result.begin(),result.end());
        return result;
    }
    std::vector<const geo::GSolid*> solids() const override {
        return source_.solids();
    }
    const geo::GVertex* vertex(int tag) const override {
        return source_.vertex(tag);
    }
    const geo::GEdge* edge(int tag) const override { return source_.edge(tag); }
    const geo::GFace* face(int tag) const override { return source_.face(tag); }
    const geo::GSolid* solid(int tag) const override {
        return source_.solid(tag);
    }
private:
    const geo::Model& source_;
};

void configureProvider(const geo::Model& model,Provider& provider) {
    const auto solids=model.solids();
    CHECK(solids.size()==1);
    if (solids.size()!=1) return;
    Vec3 lo{},hi{};
    solids.front()->bbox(lo,hi);
    const Vec3 centre=(lo+hi)*0.5;
    for (const auto* face:model.faces()) {
        Vec2 uvLo{},uvHi{};
        face->paramBounds(uvLo,uvHi);
        const Vec2 uv=(uvLo+uvHi)*0.5;
        const int sense=face->normalAt(uv).dot(
            face->pointAt(uv)-centre)>0.0?1:-1;
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            {{0,sense}};
    }
}

algo::T4ARequest reportRequest(const Provider& provider) {
    algo::T4ARequest request;
    request.source.conformingRequest.expectedSurfaceSetStamp={17001,1};
    request.source.conformingRequest.expectedOwnershipProviderStamp=
        provider.current;
    request.source.conformingRequest.envelopeMode=
        algo::RecoveryEnvelopeMode::DerivedFiniteBox;
    request.source.refinement.maxEdgeToTarget=1.0e9;
    request.source.refinement.maxRadiusEdge=1.0e9;
    request.source.refinement.maxInsertions=0;
    request.source.refinement.maxWorkUnits=200000;
    request.source.refinement.maxOutputNodes=100000;
    request.source.refinement.maxOutputTets=100000;
    request.source.refinement.maxReportRows=5;
    request.source.refinement.maxReportWitnesses=100000;
    request.source.refinement.maxReportTopologyHandles=300000;
    request.quality.mode=algo::T4AQualityMode::ReportOnly;
    request.quality.minimumRadiusRatio=0.2;
    request.quality.maxResultNodes=100000;
    request.quality.maxResultTets=100000;
    request.quality.maxReportRows=1000;
    request.quality.maxWorstTailRows=64;
    return request;
}

algo::T4ARequest boundedRequest(const Provider& provider,
                                std::uint64_t operations=1) {
    auto request=reportRequest(provider);
    request.quality.mode=algo::T4AQualityMode::BoundedImprove;
    request.quality.minimumRadiusRatio=0.2;
    request.quality.minimumImprovement=1.0e-8;
    request.quality.maxOperations=operations;
    request.quality.maxSelectedGenerations=64;
    request.quality.maxCandidateAttemptsPerCell=1;
    request.quality.maxCreatedTetsPerOperation=128;
    request.quality.maxWorkUnits=1000000;
    request.quality.maxNewNodes=operations;
    request.quality.maxResultNodes=100000;
    request.quality.maxResultTets=100000;
    request.quality.maxReportRows=200000;
    request.quality.maxWorstTailRows=64;
    request.quality.maxDependencyRowsPerAttempt=1000;
    request.quality.maxProtectedWitnessesPerAttempt=1000;
    request.quality.maxTopologyRowsPerOperation=10000;
    return request;
}

void requireGenerationAccounting(const algo::T4AOutcome& outcome) {
    CHECK(outcome.report!=nullptr&&outcome.retainedResult!=nullptr);
    if (!outcome.report||!outcome.retainedResult) return;
    const auto& report=*outcome.report;
    std::set<algo::T4AQualityKey> owners;
    for (std::uint64_t i=0;i<report.commitsUsed;++i)
        CHECK(owners.insert(report.commitStorage[
            static_cast<std::size_t>(i)].selected).second);
    for (std::uint64_t i=0;i<report.blockedUsed;++i)
        CHECK(owners.insert(report.blockedStorage[
            static_cast<std::size_t>(i)].selected).second);
    if (report.hasTerminal) CHECK(owners.insert(report.terminal.selected).second);
    CHECK(owners.size()==outcome.retainedResult->metrics.selectedGenerations);
    for (std::uint64_t i=0;i<report.unattemptedUsed;++i) {
        const auto& row=report.unattemptedStorage[static_cast<std::size_t>(i)];
        CHECK(!owners.count(row));
        if (i) CHECK(report.unattemptedStorage[static_cast<std::size_t>(i-1)]<
                     row);
    }
    CHECK(report.attemptsUsed==
          outcome.retainedResult->metrics.selectedGenerations);
    const auto& trace=algo::t4a_detail::finalLiveTraceForTesting();
    if (trace.enabled) {
        CHECK(trace.complete);
        CHECK(trace.used==report.unattemptedUsed);
        if (trace.complete&&trace.used==report.unattemptedUsed)
            for (std::uint64_t i=0;i<report.unattemptedUsed;++i)
                CHECK(trace.rows[static_cast<std::size_t>(i)]==
                      report.unattemptedStorage[static_cast<std::size_t>(i)]);
    }
}

void appendVertex(std::vector<std::uint64_t>& words,
                  algo::RefinementVertexKey value) {
    words.push_back(static_cast<std::uint64_t>(value.kind));
    words.push_back(value.value);
}

void appendQualityKey(std::vector<std::uint64_t>& words,
                      const algo::T4AQualityKey& value) {
    words.push_back(value.radiusRatioBits);
    for (const auto vertex:value.tet.vertices) appendVertex(words,vertex);
    words.push_back(value.generation);
}

template<class Tag>
void appendHandle(std::vector<std::uint64_t>& words,
                  topo::Handle<Tag> value) {
    words.push_back(static_cast<std::uint32_t>(value.slot));
    words.push_back(value.gen);
}

void appendSummary(std::vector<std::uint64_t>& words,
                   const algo::T4AQualitySummary& value) {
    words.insert(words.end(),{value.tetCount,value.exactPositive,
        value.nonfinite,value.belowThreshold,
        std::bit_cast<std::uint64_t>(value.minimum),
        std::bit_cast<std::uint64_t>(value.p1),
        std::bit_cast<std::uint64_t>(value.p5),
        std::bit_cast<std::uint64_t>(value.p50)});
    words.insert(words.end(),value.histogram.begin(),value.histogram.end());
    words.insert(words.end(),{value.tailOffset,value.tailCount,
                              value.tailTruncated?1U:0U});
}

std::vector<std::uint64_t> normalizedOutcomeWords(
    const algo::T4AOutcome& outcome) {
    std::vector<std::uint64_t> words;
    words.reserve(1024+(outcome.report?
        static_cast<std::size_t>(outcome.report->attemptsUsed*50+
            outcome.report->blockedUsed*30+outcome.report->topologyUsed*12):0));
    words.insert(words.end(),{
        static_cast<std::uint64_t>(outcome.status),
        static_cast<std::uint64_t>(outcome.failure.reason),
        static_cast<std::uint64_t>(outcome.failure.phase),
        outcome.failure.hasCell?1U:0U,
        outcome.destinationUnchanged?1U:0U,
        outcome.reportComplete?1U:0U});
    if (outcome.failure.hasCell)
        for (const auto vertex:outcome.failure.cell.vertices)
            appendVertex(words,vertex);
    if (outcome.sourceReceipt) {
        words.insert(words.end(),{1,
            static_cast<std::uint64_t>(outcome.sourceReceipt->status),
            outcome.sourceReceipt->reportComplete?1U:0U,
            outcome.sourceReceipt->strictHash,
            outcome.sourceReceipt->semanticHash.value,
            outcome.sourceReceipt->semanticHash.canonical?1U:0U,
            outcome.sourceReceipt->nodeCount,outcome.sourceReceipt->cellCount,
            outcome.sourceReceipt->nodeMappingDigest,
            outcome.sourceReceipt->cellMappingDigest});
    } else words.push_back(0);
    if (outcome.report) {
        const auto& report=*outcome.report;
        words.insert(words.end(),{1,report.shellPrepared?1U:0U,
            report.reportComplete?1U:0U,report.attemptsUsed,
            report.blockedUsed,report.commitsUsed,report.dependenciesUsed,
            report.witnessesUsed,report.topologyUsed,
            report.unattemptedUsed,report.tailUsed,
            report.hasTerminal?1U:0U});
        appendSummary(words,report.before);
        appendSummary(words,report.after);
        for (const auto& row:report.attempts()) {
            words.insert(words.end(),{row.attemptOrdinal,row.operationOrdinal});
            appendQualityKey(words,row.selected);
            words.insert(words.end(),{
                std::bit_cast<std::uint64_t>(row.proposedCoordinate.x),
                std::bit_cast<std::uint64_t>(row.proposedCoordinate.y),
                std::bit_cast<std::uint64_t>(row.proposedCoordinate.z),
                row.removedCells,row.createdCells,row.dependencyOffset,
                row.dependencyCount,row.witnessOffset,row.witnessCount,
                std::bit_cast<std::uint64_t>(row.selectedQuality),
                std::bit_cast<std::uint64_t>(row.createdMinimumQuality),
                static_cast<std::uint64_t>(row.disposition),
                static_cast<std::uint64_t>(row.reason),
                row.workBefore,row.workAfter});
        }
        for (const auto& row:report.blocked()) {
            appendQualityKey(words,row.selected);
            words.insert(words.end(),{static_cast<std::uint64_t>(row.reason),
                row.attemptOrdinal,row.dependencyOffset,row.dependencyCount,
                row.witnessOffset,row.witnessCount,
                row.invalidatedByOperation});
        }
        for (std::uint64_t i=0;i<report.commitsUsed;++i) {
            const auto& row=report.commitStorage[static_cast<std::size_t>(i)];
            words.insert(words.end(),{row.operationOrdinal,row.attemptOrdinal});
            appendQualityKey(words,row.selected);
            appendHandle(words,row.predictedNode);
            appendHandle(words,row.actualNode);
            words.insert(words.end(),{row.topologyOffset,row.topologyCount,
                std::bit_cast<std::uint64_t>(row.selectedQuality),
                std::bit_cast<std::uint64_t>(row.createdMinimumQuality),
                row.postCommitStamp.revision});
        }
        for (std::uint64_t i=0;i<report.dependenciesUsed;++i) {
            const auto& row=report.dependencyStorage[static_cast<std::size_t>(i)];
            words.push_back(static_cast<std::uint64_t>(row.kind));
            words.push_back(row.vertexCount);
            for (const auto vertex:row.vertices) appendVertex(words,vertex);
        }
        for (std::uint64_t i=0;i<report.witnessesUsed;++i) {
            const auto& row=report.witnessStorage[static_cast<std::size_t>(i)];
            words.push_back(static_cast<std::uint64_t>(row.kind));
            words.push_back(row.edge.canonicalGlobalNodes[0]);
            words.push_back(row.edge.canonicalGlobalNodes[1]);
            words.push_back(row.facet.facePersistentId);
            words.push_back(row.facet.normalizedLocalCellOrdinal);
        }
        for (std::uint64_t i=0;i<report.topologyUsed;++i) {
            const auto& row=report.topologyStorage[static_cast<std::size_t>(i)];
            words.push_back(static_cast<std::uint64_t>(row.kind));
            words.push_back(static_cast<std::uint64_t>(row.transition));
            appendHandle(words,row.predictedNode); appendHandle(words,row.actualNode);
            appendHandle(words,row.predictedCell); appendHandle(words,row.actualCell);
            appendHandle(words,row.predictedEdge); appendHandle(words,row.actualEdge);
            appendHandle(words,row.predictedFace); appendHandle(words,row.actualFace);
        }
        for (std::uint64_t i=0;i<report.unattemptedUsed;++i)
            appendQualityKey(words,report.unattemptedStorage[
                static_cast<std::size_t>(i)]);
        for (std::uint64_t i=0;i<report.tailUsed;++i)
            appendQualityKey(words,report.tailStorage[static_cast<std::size_t>(i)]);
        if (report.hasTerminal) {
            appendQualityKey(words,report.terminal.selected);
            words.insert(words.end(),{
                static_cast<std::uint64_t>(report.terminal.reason),
                static_cast<std::uint64_t>(report.terminal.phase),
                report.terminal.attemptOrdinal});
        }
    } else words.push_back(0);
    if (outcome.retainedResult) {
        const auto& result=*outcome.retainedResult;
        const auto& m=result.metrics;
        words.push_back(1);
        appendSummary(words,result.before); appendSummary(words,result.after);
        words.insert(words.end(),{
            m.initialCellsScanned,m.finalCellsAudited,m.qualityEvaluations,
            m.queuePushes,m.queuePops,m.staleRows,m.queueHighWater,
            m.selectedGenerations,m.candidateAttempts,m.formulaRejects,
            m.duplicateRejects,m.protectedIndexVisits,m.locationIndexVisits,
            m.locationTetTests,m.cavityCellsVisited,m.cavityFaceVisits,
            m.survivorIndexVisits,m.candidateSurvivorTests,
            m.localRowsPrepared,m.localRowsApplied,m.blockedRows,
            m.blockedInvalidations,m.preparedTransactions,
            m.committedOperations,m.perOperationGlobalTetScans,
            m.contextRebuilds,m.phase3Allocations,
            m.publicationTailAllocations,m.workUnitsConsumed,
            m.transientAccountedHighWaterBytes,m.runComplete?1U:0U,
            m.qualityComplete?1U:0U,
            result.boundaryAudit.prescribedExterior,
            result.boundaryAudit.matchedExterior,
            result.boundaryAudit.orphanExterior,
            result.boundaryAudit.unexpectedExterior,
            result.boundaryAudit.invertedTets,
            result.finalStrictHash,result.finalSemanticHash.value,
            result.finalSemanticHash.canonical?1U:0U,
            result.publicationReady?1U:0U});
    } else words.push_back(0);
    return words;
}

void exactQualityMetric() {
    TEST("T4-A canonical exact-positive radius-ratio metric");
    std::array<Vec3,4> points{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}};
    std::array<algo::RefinementVertexKey,4> keys{};
    for (std::size_t i=0;i<keys.size();++i) {
        keys[i].kind=algo::RefinementVertexKind::ExistingGlobal;
        keys[i].value=static_cast<std::uint64_t>(i);
    }
    double first=0.0;
    CHECK(algo::t4a_detail::checkedTetQualityForTesting(points,keys,first));
    CHECK(first>0.0&&first<=1.0);
    std::reverse(points.begin(),points.end());
    std::reverse(keys.begin(),keys.end());
    double reversed=0.0;
    CHECK(algo::t4a_detail::checkedTetQualityForTesting(
        points,keys,reversed));
    CHECK(reversed==first);
    points[3]=points[2];
    CHECK(!algo::t4a_detail::checkedTetQualityForTesting(
        points,keys,reversed));

    const std::array<Vec3,4> source{{
        {0,0,0},{1,0,0},{0,1,0},{0,0,1}}};
    std::array<std::size_t,4> permutation{{0,1,2,3}};
    do {
        std::array<Vec3,4> p{};
        std::array<algo::RefinementVertexKey,4> k{};
        for (std::size_t i=0;i<4;++i) {
            p[i]=source[permutation[i]];
            k[i]={algo::RefinementVertexKind::ExistingGlobal,
                  permutation[i]};
        }
        if (predicates::orient3d(p[0],p[1],p[2],p[3])==
                predicates::Sign::Negative) {
            std::swap(p[0],p[1]);
            std::swap(k[0],k[1]);
        }
        double value=0.0;
        CHECK(algo::t4a_detail::checkedTetQualityForTesting(p,k,value));
        CHECK(std::bit_cast<std::uint64_t>(value)==
              std::bit_cast<std::uint64_t>(first));
    } while (std::next_permutation(permutation.begin(),permutation.end()));

    std::array<Vec3,4> poor{{
        {0,0,0},{1,0,0},{0,1,0},{1.0e-12,1.0e-12,1.0e-16}}};
    double poorQuality=0.0;
    CHECK(algo::t4a_detail::checkedTetQualityForTesting(
        poor,keys,poorQuality));
    CHECK(poorQuality>=0.0&&poorQuality<first);
    std::array<Vec3,4> tiny{{
        {0,0,0},{1.0e-100,0,0},{0,1.0e-100,0},{0,0,1.0e-100}}};
    double tinyQuality=-1.0;
    CHECK(algo::t4a_detail::checkedTetQualityForTesting(
        tiny,keys,tinyQuality));
    CHECK(tinyQuality==0.0&&!std::signbit(tinyQuality));
    tiny[0].x=(std::numeric_limits<double>::infinity)();
    CHECK(!algo::t4a_detail::checkedTetQualityForTesting(
        tiny,keys,tinyQuality));
}

void preflightAndSizingBoundary() {
    TEST("T4-A preflight precedence and zero post-T3 sizing calls");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Provider invalidProvider;
    configureProvider(*model,invalidProvider);
    auto invalid=reportRequest(invalidProvider);
    invalid.quality.maxResultNodes=0;
    sizing::MeshSpec constant;
    constant.field=std::make_shared<sizing::ConstantField>(0.3);
    mesh::Mesh destination;
    const auto rejected=algo::meshModelQualityTet4(
        *model,constant,invalid,destination,invalidProvider);
    CHECK(rejected.status==algo::T4AStatus::Rejected);
    CHECK(rejected.failure.reason==algo::T4AReason::InvalidRequest);
    CHECK(invalidProvider.ownershipCalls==0&&invalidProvider.envelopeCalls==0);

    Provider sourceProvider;
    configureProvider(*model,sourceProvider);
    auto sourceInvalid=reportRequest(sourceProvider);
    sourceInvalid.source.refinement.maxEdgeToTarget=0.0;
    const auto sourceRejected=algo::meshModelQualityTet4(
        *model,constant,sourceInvalid,destination,sourceProvider);
    CHECK(sourceRejected.status==algo::T4AStatus::Rejected);
    CHECK(sourceRejected.failure.reason==algo::T4AReason::SourceT3AStopped);
    CHECK(sourceRejected.sourceOutcome&&
          sourceRejected.sourceOutcome->status==algo::T3AStatus::Rejected);
    CHECK(sourceRejected.sourceOutcome&&
          sourceRejected.sourceOutcome->failure.reason==
              algo::T3AReason::InvalidRefinementPolicy);
    CHECK(sourceProvider.ownershipCalls==0&&sourceProvider.envelopeCalls==0);

    Provider directProvider;
    configureProvider(*model,directProvider);
    auto directField=std::make_shared<CountingField>(0.3);
    sizing::MeshSpec directSpec;
    directSpec.field=directField;
    const auto request=reportRequest(directProvider);
    mesh::Mesh directMesh;
    const auto direct=algo::meshModelRefinedTet4(
        *model,directSpec,request.source,directMesh,directProvider);
    CHECK(direct.published());
    const auto directCalls=directField->calls;
    CHECK(directCalls>0);

    Provider t4Provider;
    configureProvider(*model,t4Provider);
    auto t4Field=std::make_shared<CountingField>(0.3);
    sizing::MeshSpec t4Spec;
    t4Spec.field=t4Field;
    auto t4Request=reportRequest(t4Provider);
    mesh::Mesh t4Mesh;
    const auto t4=algo::meshModelQualityTet4(
        *model,t4Spec,t4Request,t4Mesh,t4Provider);
    CHECK(t4.published());
    CHECK(t4Field->calls==directCalls);
    CHECK(mesh::meshHash(t4Mesh)==mesh::meshHash(directMesh));
}

void reportOnlyAndCapacity() {
    TEST("T4-A consuming report-only handoff and exact result capacity gate");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureProvider(*model,provider);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.3);

    auto request=reportRequest(provider);
    mesh::Mesh published;
    const auto outcome=algo::meshModelQualityTet4(
        *model,spec,request,published,provider);
    if (!outcome.published())
        std::printf("  t4a report diagnostic status=%d reason=%d phase=%d "
                    "source=%d/%d\n",static_cast<int>(outcome.status),
            static_cast<int>(outcome.failure.reason),
            static_cast<int>(outcome.failure.phase),
            outcome.sourceOutcome?static_cast<int>(outcome.sourceOutcome->status):-1,
            outcome.sourceOutcome?static_cast<int>(outcome.sourceOutcome->failure.reason):-1);
    CHECK(outcome.published());
    CHECK(!outcome.destinationUnchanged);
    CHECK(outcome.report&&outcome.report->reportComplete);
    CHECK(outcome.sourceReceipt&&outcome.sourceReceipt->reportComplete);
    CHECK(outcome.retainedResult&&outcome.retainedResult->publicationReady);
    CHECK(outcome.retainedResult&&
          std::string_view(outcome.retainedResult->routeId)==
              algo::kT4ARouteId);
    CHECK(outcome.retainedResult&&outcome.retainedResult->before.tetCount>0);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->before.tetCount==
              outcome.retainedResult->after.tetCount);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.committedOperations==0);
    CHECK(outcome.retainedResult&&
          !outcome.retainedResult->metrics.qualityComplete);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->boundaryAudit.orphanExterior==0);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->boundaryAudit.unexpectedExterior==0);
    CHECK(outcome.sourceReceipt&&
          outcome.sourceReceipt->strictHash==mesh::meshHash(published));

    auto insufficient=request;
    insufficient.quality.maxResultNodes=1;
    insufficient.quality.maxResultTets=1;
    mesh::Mesh unchanged;
    unchanged.addNode({-7,-8,-9},{0,77});
    const auto digest=mesh::meshHash(unchanged);
    const auto stopped=algo::meshModelQualityTet4(
        *model,spec,insufficient,unchanged,provider);
    if (stopped.status!=algo::T4AStatus::Incomplete)
        std::printf("  t4a capacity diagnostic status=%d reason=%d phase=%d "
                    "source=%d/%d\n",static_cast<int>(stopped.status),
            static_cast<int>(stopped.failure.reason),
            static_cast<int>(stopped.failure.phase),
            stopped.sourceOutcome?static_cast<int>(stopped.sourceOutcome->status):-1,
            stopped.sourceOutcome?static_cast<int>(stopped.sourceOutcome->failure.reason):-1);
    CHECK(stopped.status==algo::T4AStatus::Incomplete);
    CHECK(stopped.failure.reason==algo::T4AReason::ResultBudgetReached);
    CHECK(stopped.failure.phase==algo::T4APhase::ReportEnvelope);
    CHECK(stopped.destinationUnchanged);
    CHECK(mesh::meshHash(unchanged)==digest);
    CHECK(stopped.report&&stopped.report->reportComplete);
    CHECK(stopped.sourceOutcome&&stopped.sourceOutcome->published());
    CHECK(stopped.sourceOutcome&&stopped.sourceOutcome->reportComplete);
    CHECK(!stopped.sourceReceipt&&!stopped.retainedResult);

    auto bounded=boundedRequest(provider);
    mesh::Mesh improved;
    const auto quality=algo::meshModelQualityTet4(
        *model,spec,bounded,improved,provider);
    if (!quality.retainedResult||
        quality.retainedResult->metrics.committedOperations==0)
        std::printf("  t4a bounded diagnostic status=%d reason=%d phase=%d "
                    "before=%.17g commits=%llu blocked=%llu\n",
            static_cast<int>(quality.status),
            static_cast<int>(quality.failure.reason),
            static_cast<int>(quality.failure.phase),
            outcome.retainedResult?outcome.retainedResult->before.minimum:-1.0,
            static_cast<unsigned long long>(quality.retainedResult?
                quality.retainedResult->metrics.committedOperations:0),
            static_cast<unsigned long long>(quality.report?
                quality.report->blockedUsed:0));
    if (quality.report&&quality.retainedResult&&
        quality.retainedResult->metrics.committedOperations==0)
        for (const auto& row:quality.report->attempts())
            std::printf("    attempt=%llu disposition=%d reason=%d q=%.17g "
                        "created=%.17g removed=%llu made=%llu\n",
                static_cast<unsigned long long>(row.attemptOrdinal),
                static_cast<int>(row.disposition),static_cast<int>(row.reason),
                row.selectedQuality,row.createdMinimumQuality,
                static_cast<unsigned long long>(row.removedCells),
                static_cast<unsigned long long>(row.createdCells));
    CHECK(quality.status==algo::T4AStatus::Published||
          quality.status==algo::T4AStatus::Incomplete);
    CHECK(quality.retainedResult!=nullptr);
    CHECK(quality.retainedResult&&
          quality.retainedResult->metrics.committedOperations==1);
    CHECK(quality.status==algo::T4AStatus::Incomplete);
    CHECK(quality.failure.reason==algo::T4AReason::OperationBudgetReached);
    CHECK(quality.destinationUnchanged);
    CHECK(quality.retainedResult&&
          quality.retainedResult->metrics.perOperationGlobalTetScans==0);
    CHECK(quality.retainedResult&&
          quality.retainedResult->metrics.contextRebuilds==0);
    CHECK(quality.retainedResult&&
          quality.retainedResult->metrics.phase3Allocations==0);
    CHECK(quality.retainedResult&&
          quality.retainedResult->metrics.publicationTailAllocations==0);
    CHECK(quality.retainedResult&&quality.retainedResult->metrics.runComplete);
    CHECK(quality.retainedResult&&
          quality.retainedResult->after.minimum>=
              quality.retainedResult->before.minimum);
    CHECK(quality.retainedResult&&
          quality.retainedResult->boundaryAudit.orphanExterior==0);
    CHECK(quality.retainedResult&&
          quality.retainedResult->boundaryAudit.unexpectedExterior==0);
    CHECK(quality.sourceReceipt&&outcome.sourceReceipt&&
          quality.sourceReceipt->strictHash==outcome.sourceReceipt->strictHash);
    CHECK(quality.sourceReceipt&&outcome.sourceReceipt&&
          quality.sourceReceipt->semanticHash.value==
              outcome.sourceReceipt->semanticHash.value&&
          quality.sourceReceipt->semanticHash.canonical==
              outcome.sourceReceipt->semanticHash.canonical);
    CHECK(quality.sourceReceipt&&outcome.sourceReceipt&&
          quality.sourceReceipt->nodeMappingDigest==
              outcome.sourceReceipt->nodeMappingDigest);
    CHECK(quality.sourceReceipt&&outcome.sourceReceipt&&
          quality.sourceReceipt->cellMappingDigest==
              outcome.sourceReceipt->cellMappingDigest);
    CHECK(quality.report&&quality.report->commitsUsed==1);
    CHECK(quality.report&&quality.report->blockedUsed==0);
    CHECK(quality.report&&quality.report->attemptsUsed==1);
    CHECK(quality.retainedResult&&
          quality.retainedResult->metrics.staleRows>0);
    if (quality.report&&quality.retainedResult) {
        const auto& qm=quality.retainedResult->metrics;
        // The selected generation owns the commit. Queue rows retired only
        // because they share its cavity are represented by staleRows and are
        // deliberately absent from both processed and unattempted owners.
        CHECK(qm.queuePushes==quality.report->commitsUsed+
              quality.report->blockedUsed+
              (quality.report->hasTerminal?1U:0U)+qm.staleRows+
              quality.report->unattemptedUsed);
        CHECK(quality.report->commitStorage[0].predictedNode==
              quality.report->commitStorage[0].actualNode);
        const auto& commit=quality.report->commitStorage[0];
        bool retainedProtectedFace=false;
        for (std::uint64_t i=0;i<commit.topologyCount;++i) {
            const auto& row=quality.report->topologyStorage[
                static_cast<std::size_t>(commit.topologyOffset+i)];
            switch (row.kind) {
            case algo::T3ATopologyKind::Node:
                CHECK(row.predictedNode==row.actualNode); break;
            case algo::T3ATopologyKind::Cell:
                CHECK(row.predictedCell==row.actualCell); break;
            case algo::T3ATopologyKind::Edge:
                CHECK(row.predictedEdge==row.actualEdge); break;
            case algo::T3ATopologyKind::Face:
                CHECK(row.predictedFace==row.actualFace);
                if (row.transition==algo::T3ATopologyTransition::Retained&&
                    quality.retainedResult->artifact&&
                    std::any_of(
                        quality.retainedResult->artifact->immutableBoundary.
                            begin(),
                        quality.retainedResult->artifact->immutableBoundary.
                            end(),
                        [&](const auto& boundary) {
                            return boundary.authoritativeFace==
                                   row.predictedFace;
                        }))
                    retainedProtectedFace=true;
                break;
            }
        }
        CHECK(retainedProtectedFace);
        CHECK(quality.report->commitStorage[0].createdMinimumQuality>=
              quality.report->commitStorage[0].selectedQuality+
                  bounded.quality.minimumImprovement);
        CHECK(quality.retainedResult->artifact&&
              quality.retainedResult->artifact->authoritative&&
              quality.retainedResult->artifact->authoritative->workspace);
        if (quality.retainedResult->artifact&&
            quality.retainedResult->artifact->authoritative&&
            quality.retainedResult->artifact->authoritative->workspace) {
            const auto& artifact=*quality.retainedResult->artifact;
            const auto& workspace=*artifact.authoritative->workspace;
            const auto identity=artifact.vertexKeyByHandle.find(
                commit.actualNode);
            CHECK(identity!=artifact.vertexKeyByHandle.end());
            if (identity!=artifact.vertexKeyByHandle.end()) {
                CHECK(identity->second.kind==
                      algo::RefinementVertexKind::InteriorInsertion);
                CHECK(identity->second.value==
                      artifact.t3CommittedInsertionCount);
            }
            const auto newSlot=static_cast<std::size_t>(commit.actualNode.slot);
            CHECK(newSlot<artifact.publishedGlobalByAuthoritativeNodeSlot.size());
            CHECK(newSlot>=artifact.sourceGlobalByAuthoritativeNodeSlot.size()||
                  artifact.sourceGlobalByAuthoritativeNodeSlot[newSlot]==
                      kInvalidIndex);
            bool selectedSurvived=false;
            for (const auto cell:workspace.cells().activeBySequence()) {
                const auto data=workspace.cells().get(cell);
                algo::T3ATetKey key;
                for (std::size_t i=0;i<4;++i)
                    key.vertices[i]=artifact.vertexKeyByHandle.at(data.nodes[i]);
                std::sort(key.vertices.begin(),key.vertices.end());
                selectedSurvived=selectedSurvived||
                    key==commit.selected.tet;
            }
            CHECK(!selectedSurvived);
        }
    }
    requireGenerationAccounting(quality);
}

void blockerContinuationAndInvalidation() {
    TEST("T4-A blocker continuation and dependency invalidation");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeCylinder({0,0,0},0.5,1.0);
    Provider provider;
    configureProvider(*model,provider);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.15);
    auto request=boundedRequest(provider,1);
    request.quality.maxCreatedTetsPerOperation=2048;
    request.quality.maxTopologyRowsPerOperation=100000;
    request.quality.maxDependencyRowsPerAttempt=10000;
    request.quality.maxProtectedWitnessesPerAttempt=10000;
    request.quality.maxReportRows=2000000;
    request.quality.maxResultNodes=1000000;
    request.quality.maxResultTets=1000000;

    algo::t4a_detail::finalLiveTraceForTesting().reset(true);
    mesh::Mesh destination;
    const auto outcome=algo::meshModelQualityTet4(
        *model,spec,request,destination,provider);
    CHECK(outcome.status==algo::T4AStatus::Incomplete);
    CHECK(outcome.failure.reason==algo::T4AReason::OperationBudgetReached);
    CHECK(outcome.destinationUnchanged);
    CHECK(outcome.report&&outcome.report->reportComplete);
    CHECK(outcome.retainedResult&&outcome.retainedResult->metrics.runComplete);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.committedOperations==1);
    CHECK(outcome.report&&outcome.report->blockedUsed==5);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.blockedInvalidations==5);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.perOperationGlobalTetScans==0);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.contextRebuilds==0);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.phase3Allocations==0);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.publicationTailAllocations==0);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->boundaryAudit.orphanExterior==0);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->boundaryAudit.unexpectedExterior==0);

    bool sawBlock=false,sawCommitAfterBlock=false,sawInvalidation=false;
    if (outcome.report) {
        for (const auto& row:outcome.report->attempts()) {
            if (row.disposition==algo::T4AAttemptDisposition::
                    NoQualityProgress||
                row.disposition==algo::T4AAttemptDisposition::Encroached||
                row.disposition==algo::T4AAttemptDisposition::Blocked)
                sawBlock=true;
            if (sawBlock&&row.disposition==
                    algo::T4AAttemptDisposition::Committed)
                sawCommitAfterBlock=true;
        }
        for (const auto& row:outcome.report->blocked()) {
            CHECK(row.attemptOrdinal<outcome.report->attemptsUsed);
            CHECK(row.dependencyCount>0);
            if (row.invalidatedByOperation!=algo::kCurrentT4Operation) {
                sawInvalidation=true;
                CHECK(row.invalidatedByOperation<
                      outcome.retainedResult->metrics.committedOperations);
            }
        }
    }
    CHECK(sawBlock&&sawCommitAfterBlock&&sawInvalidation);
    requireGenerationAccounting(outcome);
    algo::t4a_detail::finalLiveTraceForTesting().reset(false);
}

void protectedBlockerContinuation() {
    TEST("T4-A protected blocker continues to an independent improvement");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureProvider(*model,provider);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.3);
    auto request=boundedRequest(provider,1);
    algo::t4a_detail::forcedCandidateForTesting()={0.471234,0.4321,0.0};
    algo::t4a_detail::forceCandidateForTesting()=true;
    algo::t4a_detail::finalLiveTraceForTesting().reset(true);
    mesh::Mesh destination;
    const auto outcome=algo::meshModelQualityTet4(
        *model,spec,request,destination,provider);
    algo::t4a_detail::forceCandidateForTesting()=false;
    CHECK(outcome.status==algo::T4AStatus::Incomplete);
    CHECK(outcome.failure.reason==algo::T4AReason::OperationBudgetReached);
    CHECK(outcome.report&&outcome.report->attemptsUsed>=2);
    CHECK(outcome.report&&outcome.report->blockedUsed>=1);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.committedOperations==1);
    if (outcome.report&&outcome.report->attemptsUsed>=2) {
        const auto& first=outcome.report->attemptStorage[0];
        CHECK(first.disposition==algo::T4AAttemptDisposition::Encroached);
        CHECK(first.reason==algo::T4AReason::ProtectedFeatureBlocked);
        CHECK(first.witnessCount>0);
        CHECK(std::any_of(outcome.report->attempts().begin()+1,
                          outcome.report->attempts().end(),
            [](const auto& row) {
                return row.disposition==
                    algo::T4AAttemptDisposition::Committed;
            }));
    }
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->boundaryAudit.orphanExterior==0);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->boundaryAudit.unexpectedExterior==0);
    requireGenerationAccounting(outcome);
    algo::t4a_detail::finalLiveTraceForTesting().reset(false);
}

void budgetsAndResidualTaxonomy() {
    TEST("T4-A exact budget precedence and honest residual outcome");
    auto kernel=geo::makeBuiltinKernel();
    auto box=kernel->makeBox({0,0,0},{1,1,1});
    Provider boxProvider;
    configureProvider(*box,boxProvider);
    sizing::MeshSpec boxSpec;
    boxSpec.field=std::make_shared<sizing::ConstantField>(0.3);

    auto runBox=[&](algo::T4ARequest request) {
        mesh::Mesh destination;
        destination.addNode({-19,-23,-29},{0,880});
        const auto before=mesh::meshHash(destination);
        algo::t4a_detail::finalLiveTraceForTesting().reset(true);
        auto outcome=algo::meshModelQualityTet4(
            *box,boxSpec,request,destination,boxProvider);
        CHECK(outcome.destinationUnchanged);
        CHECK(mesh::meshHash(destination)==before);
        requireGenerationAccounting(outcome);
        algo::t4a_detail::finalLiveTraceForTesting().reset(false);
        return outcome;
    };

    auto residualRequest=boundedRequest(boxProvider,8);
    residualRequest.quality.minimumImprovement=0.02;
    residualRequest.quality.maxReportRows=1000000;
    auto residual=runBox(residualRequest);
    CHECK(residual.status==algo::T4AStatus::Incomplete);
    CHECK(residual.failure.reason==algo::T4AReason::ResidualQuality);
    CHECK(residual.retainedResult&&
          residual.retainedResult->metrics.committedOperations==1);
    CHECK(residual.report&&residual.report->blockedUsed==1);
    CHECK(residual.retainedResult&&
          !residual.retainedResult->metrics.qualityComplete);

    auto workRequest=boundedRequest(boxProvider,1);
    workRequest.quality.maxWorkUnits=1;
    auto work=runBox(workRequest);
    CHECK(work.status==algo::T4AStatus::Incomplete);
    CHECK(work.failure.reason==algo::T4AReason::WorkBudgetReached);
    CHECK(work.report&&work.report->attemptsUsed==0);
    CHECK(work.report&&work.report->unattemptedUsed>0);

    auto selectedRequest=boundedRequest(boxProvider,1);
    selectedRequest.quality.maxSelectedGenerations=1;
    selectedRequest.quality.maxReportRows=20000;
    algo::t4a_detail::forcedCandidateForTesting()={0.471234,0.4321,0.0};
    algo::t4a_detail::forceCandidateForTesting()=true;
    auto selected=runBox(selectedRequest);
    algo::t4a_detail::forceCandidateForTesting()=false;
    CHECK(selected.status==algo::T4AStatus::Incomplete);
    // With the normative one-attempt-per-generation policy, the attempt
    // envelope and selected-generation cap become exhausted together.
    // Canonical terminal precedence assigns that collision to report budget.
    CHECK(selected.failure.reason==algo::T4AReason::ReportBudgetReached);
    CHECK(selected.report&&selected.report->attemptsUsed==1);
    CHECK(selected.report&&selected.report->blockedUsed==1);

    auto reportRequestValue=boundedRequest(boxProvider,1);
    reportRequestValue.quality.minimumImprovement=0.24;
    reportRequestValue.quality.maxDependencyRowsPerAttempt=1;
    auto reportBudget=runBox(reportRequestValue);
    CHECK(reportBudget.status==algo::T4AStatus::Incomplete);
    CHECK(reportBudget.failure.reason==algo::T4AReason::ReportBudgetReached);
    CHECK(reportBudget.report&&reportBudget.report->hasTerminal);
    CHECK(reportBudget.report&&reportBudget.report->attemptsUsed==1);
    CHECK(reportBudget.report&&reportBudget.report->blockedUsed==0);
    if (reportBudget.report) {
        CHECK(reportBudget.report->terminal.reason==
              algo::T4AReason::ReportBudgetReached);
        CHECK(reportBudget.report->terminal.selected==
              reportBudget.report->attemptStorage[0].selected);
    }

    // Candidate-created-Tet capacity is a pre-transaction gate in the
    // accepted contract.  A deliberately tiny cap must therefore stop
    // without preparing topology or consuming placeholders.
    auto createdTetCapRequest=boundedRequest(boxProvider,1);
    createdTetCapRequest.quality.maxCreatedTetsPerOperation=1;
    auto createdTetCap=runBox(createdTetCapRequest);
    CHECK(createdTetCap.status==algo::T4AStatus::Incomplete);
    CHECK(createdTetCap.failure.reason==algo::T4AReason::ResultBudgetReached);
    CHECK(createdTetCap.retainedResult&&
          createdTetCap.retainedResult->metrics.preparedTransactions==0);

    auto cylinder=kernel->makeCylinder({0,0,0},0.5,1.0);
    Provider cylinderProvider;
    configureProvider(*cylinder,cylinderProvider);
    sizing::MeshSpec cylinderSpec;
    cylinderSpec.field=std::make_shared<sizing::ConstantField>(0.15);
    auto nodeRequest=boundedRequest(cylinderProvider,2);
    nodeRequest.quality.maxNewNodes=1;
    nodeRequest.quality.maxCreatedTetsPerOperation=2048;
    nodeRequest.quality.maxTopologyRowsPerOperation=100000;
    nodeRequest.quality.maxDependencyRowsPerAttempt=10000;
    nodeRequest.quality.maxProtectedWitnessesPerAttempt=10000;
    nodeRequest.quality.maxReportRows=2000000;
    nodeRequest.quality.maxResultNodes=1000000;
    nodeRequest.quality.maxResultTets=1000000;
    mesh::Mesh nodeDestination;
    const auto nodeBefore=mesh::meshHash(nodeDestination);
    algo::t4a_detail::finalLiveTraceForTesting().reset(true);
    const auto node=algo::meshModelQualityTet4(
        *cylinder,cylinderSpec,nodeRequest,nodeDestination,cylinderProvider);
    CHECK(node.status==algo::T4AStatus::Incomplete);
    CHECK(node.failure.reason==algo::T4AReason::NodeBudgetReached);
    CHECK(node.destinationUnchanged);
    CHECK(mesh::meshHash(nodeDestination)==nodeBefore);
    CHECK(node.retainedResult&&
          node.retainedResult->metrics.committedOperations==1);
    requireGenerationAccounting(node);
    algo::t4a_detail::finalLiveTraceForTesting().reset(false);
}

void postCommitQuarantine() {
    TEST("T4-A post-commit invariant defect quarantines private result");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureProvider(*model,provider);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.3);
    auto request=boundedRequest(provider,1);
    mesh::Mesh destination;
    destination.addNode({-83,-89,-97},{0,996});
    const auto before=mesh::meshHash(destination);
    algo::t4a_detail::injectPostCommitAuditFailureForTesting()=true;
    const auto outcome=algo::meshModelQualityTet4(
        *model,spec,request,destination,provider);
    algo::t4a_detail::injectPostCommitAuditFailureForTesting()=false;
    CHECK(outcome.status==algo::T4AStatus::Failed);
    CHECK(outcome.failure.reason==algo::T4AReason::InternalError);
    CHECK(outcome.failure.phase==algo::T4APhase::LocalAudit);
    CHECK(outcome.destinationUnchanged);
    CHECK(mesh::meshHash(destination)==before);
    CHECK(outcome.sourceReceipt!=nullptr);
    CHECK(outcome.retainedResult==nullptr);

    mesh::Mesh retryDestination;
    const auto retry=algo::meshModelQualityTet4(
        *model,spec,request,retryDestination,provider);
    CHECK(retry.status==algo::T4AStatus::Incomplete);
    CHECK(retry.failure.reason==algo::T4AReason::OperationBudgetReached);
    CHECK(retry.retainedResult&&
          retry.retainedResult->metrics.committedOperations==1);
}

void candidateTaxonomy() {
    TEST("T4-A duplicate nonfinite and required-over-one candidate taxonomy");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureProvider(*model,provider);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.3);

    auto runForced=[&](Vec3 point) {
        auto request=boundedRequest(provider,1);
        algo::t4a_detail::forcedCandidateForTesting()=point;
        algo::t4a_detail::forceCandidateForTesting()=true;
        mesh::Mesh destination;
        destination.addNode({-101,-103,-107},{0,998});
        const auto before=mesh::meshHash(destination);
        auto outcome=algo::meshModelQualityTet4(
            *model,spec,request,destination,provider);
        algo::t4a_detail::forceCandidateForTesting()=false;
        CHECK(outcome.destinationUnchanged);
        CHECK(mesh::meshHash(destination)==before);
        return outcome;
    };
    const auto duplicate=runForced({0,0,0});
    CHECK(duplicate.status==algo::T4AStatus::Incomplete);
    CHECK(duplicate.report&&duplicate.report->attemptsUsed>=2);
    if (duplicate.report&&duplicate.report->attemptsUsed>=2) {
        CHECK(duplicate.report->attemptStorage[0].disposition==
              algo::T4AAttemptDisposition::DuplicateProposal);
        CHECK(duplicate.report->attemptStorage[0].reason==
              algo::T4AReason::NoImprovingCavity);
        CHECK(std::any_of(duplicate.report->attempts().begin()+1,
                          duplicate.report->attempts().end(),
            [](const auto& row) {
                return row.disposition==
                    algo::T4AAttemptDisposition::Committed;
            }));
    }

    const auto nonfinite=runForced({
        (std::numeric_limits<double>::infinity)(),0,0});
    CHECK(nonfinite.status==algo::T4AStatus::Failed);
    CHECK(nonfinite.failure.reason==algo::T4AReason::InvalidCandidate);
    CHECK(nonfinite.failure.phase==algo::T4APhase::Encroachment);
    CHECK(nonfinite.retainedResult==nullptr);
    CHECK(nonfinite.sourceReceipt!=nullptr);

    auto impossible=boundedRequest(provider,1);
    impossible.quality.minimumImprovement=0.95;
    impossible.quality.maxSelectedGenerations=1;
    impossible.quality.maxReportRows=20000;
    mesh::Mesh impossibleDestination;
    const auto impossibleOutcome=algo::meshModelQualityTet4(
        *model,spec,impossible,impossibleDestination,provider);
    CHECK(impossibleOutcome.status==algo::T4AStatus::Incomplete);
    CHECK(impossibleOutcome.report&&
          impossibleOutcome.report->attemptsUsed==1);
    CHECK(impossibleOutcome.report&&
          impossibleOutcome.report->attemptStorage[0].disposition==
              algo::T4AAttemptDisposition::NoQualityProgress);
    CHECK(impossibleOutcome.report&&
          impossibleOutcome.report->attemptStorage[0].reason==
              algo::T4AReason::NoImprovingCavity);
}

const algo::t4a_detail::T4AStateDigest* findT4Digest(
    const algo::t4a_detail::T4ADigestTrace& trace,
    algo::t4a_detail::T4ADigestCheckpoint checkpoint,bool postFailure) {
    for (std::uint32_t i=0;i<trace.used;++i)
        if (trace.rows[i].checkpoint==checkpoint&&
            trace.rows[i].postFailure==postFailure)
            return &trace.rows[i];
    return nullptr;
}

void checkpointDigests() {
    TEST("T4-A four checkpoint digests and allocation rollback tiers");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureProvider(*model,provider);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.3);
    const std::array checkpoints{
        algo::t4a_detail::T4ADigestCheckpoint::PreSelection,
        algo::t4a_detail::T4ADigestCheckpoint::SelectedRowSecured,
        algo::t4a_detail::T4ADigestCheckpoint::PrePreparation,
        algo::t4a_detail::T4ADigestCheckpoint::Prepared};
    for (const auto checkpoint:checkpoints) {
        auto request=boundedRequest(provider,1);
        auto& trace=algo::t4a_detail::digestTraceForTesting();
        trace.arm(checkpoint);
        mesh::Mesh destination;
        destination.addNode({-109,-113,-127},{0,995});
        const auto before=mesh::meshHash(destination);
        auto outcome=algo::meshModelQualityTet4(
            *model,spec,request,destination,provider);
        CHECK(outcome.status==algo::T4AStatus::Failed);
        CHECK(outcome.failure.reason==algo::T4AReason::AllocationFailed);
        CHECK(outcome.destinationUnchanged);
        CHECK(mesh::meshHash(destination)==before);
        CHECK(outcome.sourceReceipt!=nullptr);
        CHECK(outcome.retainedResult==nullptr);
        CHECK(!trace.faultPending);
        const auto* pre=findT4Digest(trace,checkpoint,false);
        const auto* post=findT4Digest(trace,checkpoint,true);
        CHECK(pre!=nullptr&&post!=nullptr);
        if (!pre||!post) continue;
        if (checkpoint==
                algo::t4a_detail::T4ADigestCheckpoint::Prepared) {
            CHECK(pre->prepared!=0);
            CHECK(post->prepared==0);
            const auto* logical=findT4Digest(
                trace,algo::t4a_detail::T4ADigestCheckpoint::PrePreparation,
                false);
            CHECK(logical!=nullptr);
            if (logical) {
                CHECK(post->structural==logical->structural);
                CHECK(post->revision==logical->revision);
                CHECK(post->pools==logical->pools);
            }
        } else {
            CHECK(post->structural==pre->structural);
            CHECK(post->revision==pre->revision);
            CHECK(post->pools==pre->pools);
        }
        if (checkpoint==
                algo::t4a_detail::T4ADigestCheckpoint::PreSelection) {
            CHECK(pre->report==post->report);
            CHECK(post->attemptsUsed==0);
            CHECK(!outcome.report->hasTerminal);
        } else {
            CHECK(pre->report!=post->report);
            CHECK(post->attemptsUsed==1);
            CHECK(outcome.report&&outcome.report->hasTerminal);
            CHECK(outcome.report&&outcome.report->terminal.reason==
                  algo::T4AReason::AllocationFailed);
        }
        trace.reset(false);
        mesh::Mesh retryDestination;
        auto retry=algo::meshModelQualityTet4(
            *model,spec,request,retryDestination,provider);
        CHECK(retry.published()||
              retry.status==algo::T4AStatus::Incomplete);
        CHECK(retry.retainedResult&&
              retry.retainedResult->metrics.committedOperations==1);
    }
    algo::t4a_detail::digestTraceForTesting().reset(false);
}

void completeDeterminism() {
    TEST("T4-A complete report deterministic for 1/4/auto and churn");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeCylinder({0,0,0},0.5,1.0);
    ReversedFaceModel reversed(*model);
    const std::array<int,3> threads{{1,4,0}};
    std::vector<std::uint64_t> reference;
    for (std::size_t run=0;run<threads.size();++run) {
        topo::NodePool churn;
        for (std::size_t i=0;i<run*17;++i)
            (void)churn.add({static_cast<double>(i),-61,71},{0,903});
        const geo::Model& presented=run==1
            ? static_cast<const geo::Model&>(reversed)
            : static_cast<const geo::Model&>(*model);
        Provider provider;
        configureProvider(presented,provider);
        sizing::MeshSpec spec;
        spec.field=std::make_shared<sizing::ConstantField>(0.15);
        spec.threads=threads[run];
        auto request=boundedRequest(provider,1);
        request.quality.maxCreatedTetsPerOperation=2048;
        request.quality.maxTopologyRowsPerOperation=100000;
        request.quality.maxDependencyRowsPerAttempt=10000;
        request.quality.maxProtectedWitnessesPerAttempt=10000;
        request.quality.maxReportRows=2000000;
        request.quality.maxResultNodes=1000000;
        request.quality.maxResultTets=1000000;
        algo::t4a_detail::finalLiveTraceForTesting().reset(true);
        mesh::Mesh destination;
        const auto outcome=algo::meshModelQualityTet4(
            presented,spec,request,destination,provider);
        CHECK(outcome.status==algo::T4AStatus::Incomplete);
        CHECK(outcome.failure.reason==algo::T4AReason::OperationBudgetReached);
        CHECK(outcome.report&&outcome.report->blockedUsed==5);
        CHECK(outcome.retainedResult&&
              outcome.retainedResult->metrics.blockedInvalidations==5);
        requireGenerationAccounting(outcome);
        const auto words=normalizedOutcomeWords(outcome);
        if (run==0) reference=words;
        else CHECK(words==reference);
        algo::t4a_detail::finalLiveTraceForTesting().reset(false);
    }
    core::setThreadCount(0);
}

void publicOperationScaling() {
    TEST("T4-A four-row public indexed operation scaling");
    auto kernel=geo::makeBuiltinKernel();
    const std::array<double,4> spacings{{0.3,0.24,0.18,0.12}};
    std::array<std::uint64_t,4> seed{},commits{},work{};
    for (std::size_t row=0;row<spacings.size();++row) {
        auto model=kernel->makeBox({0,0,0},{1,1,1});
        Provider provider;
        configureProvider(*model,provider);
        sizing::MeshSpec spec;
        spec.field=std::make_shared<sizing::ConstantField>(spacings[row]);
        spec.threads=1;
        auto request=boundedRequest(provider,8);
        request.source.conformingRequest.expectedSurfaceSetStamp={18000+row,1};
        request.quality.minimumRadiusRatio=0.4;
        request.quality.minimumImprovement=0.001;
        request.quality.maxCreatedTetsPerOperation=2048;
        request.quality.maxTopologyRowsPerOperation=100000;
        request.quality.maxDependencyRowsPerAttempt=10000;
        request.quality.maxProtectedWitnessesPerAttempt=10000;
        request.quality.maxReportRows=3000000;
        request.quality.maxWorkUnits=5000000;
        request.quality.maxResultNodes=1000000;
        request.quality.maxResultTets=1000000;
        request.source.refinement.maxOutputNodes=1000000;
        request.source.refinement.maxOutputTets=1000000;
        mesh::Mesh destination;
        const auto outcome=algo::meshModelQualityTet4(
            *model,spec,request,destination,provider);
        CHECK(outcome.status==algo::T4AStatus::Incomplete||outcome.published());
        CHECK(outcome.retainedResult!=nullptr);
        if (!outcome.retainedResult) continue;
        seed[row]=outcome.retainedResult->before.tetCount;
        commits[row]=outcome.retainedResult->metrics.committedOperations;
        work[row]=outcome.retainedResult->metrics.workUnitsConsumed;
        CHECK(outcome.retainedResult->metrics.perOperationGlobalTetScans==0);
        CHECK(outcome.retainedResult->metrics.contextRebuilds==0);
        CHECK(outcome.retainedResult->boundaryAudit.orphanExterior==0);
        CHECK(outcome.retainedResult->boundaryAudit.unexpectedExterior==0);
        std::printf("  t4a scale row=%zu spacing=%.3f seed=%llu "
                    "commits=%llu work=%llu\n",row,spacings[row],
            static_cast<unsigned long long>(seed[row]),
            static_cast<unsigned long long>(commits[row]),
            static_cast<unsigned long long>(work[row]));
        CHECK(commits[row]>=2);
    }
    core::setThreadCount(0);
    for (std::size_t row=1;row<spacings.size();++row) {
        CHECK(seed[row]>seed[row-1]);
        CHECK(work[row]<12*work[row-1]);
    }
    CHECK(seed.back()>=4*seed.front());
    CHECK(work.back()<12*work.front());
}

void shapeBoundaryBijection() {
    TEST("T4-A report-only box cylinder and thin-box boundary bijection");
    auto kernel=geo::makeBuiltinKernel();
    struct Fixture { std::unique_ptr<geo::Model> model; double size; };
    std::array<Fixture,3> fixtures{{
        {kernel->makeBox({0,0,0},{1,1,1}),0.3},
        {kernel->makeCylinder({0,0,0},0.5,1.0),0.15},
        {kernel->makeBox({0,0,0},{1,1,1.0e-3}),0.2}}};
    for (std::size_t i=0;i<fixtures.size();++i) {
        Provider provider;
        configureProvider(*fixtures[i].model,provider);
        sizing::MeshSpec spec;
        spec.field=std::make_shared<sizing::ConstantField>(fixtures[i].size);
        auto request=reportRequest(provider);
        request.source.conformingRequest.expectedSurfaceSetStamp={17100+i,1};
        request.quality.minimumRadiusRatio=0.0;
        mesh::Mesh output;
        const auto outcome=algo::meshModelQualityTet4(
            *fixtures[i].model,spec,request,output,provider);
        CHECK(outcome.published());
        CHECK(outcome.retainedResult&&outcome.retainedResult->metrics.qualityComplete);
        CHECK(outcome.retainedResult&&
              outcome.retainedResult->boundaryAudit.orphanExterior==0);
        CHECK(outcome.retainedResult&&
              outcome.retainedResult->boundaryAudit.unexpectedExterior==0);
        CHECK(outcome.retainedResult&&
              outcome.retainedResult->finalStrictHash==mesh::meshHash(output));
        CHECK(outcome.retainedResult&&
              outcome.retainedResult->finalSemanticHash.canonical);
    }

}

struct ReloadedCell {
    mesh::CellType type=mesh::CellType::Tet4;
    std::vector<Index> nodes;
    mesh::Classification classification{};

    friend bool operator<(const ReloadedCell& a,const ReloadedCell& b) {
        return std::tie(a.classification.dim,a.classification.tag,a.type,a.nodes)<
               std::tie(b.classification.dim,b.classification.tag,b.type,b.nodes);
    }
    friend bool operator==(const ReloadedCell& a,const ReloadedCell& b) {
        return a.type==b.type&&
               a.classification.dim==b.classification.dim&&
               a.classification.tag==b.classification.tag&&a.nodes==b.nodes;
    }
};

struct ReloadedMesh {
    std::vector<Vec3> nodes;
    std::vector<mesh::Classification> nodeClasses;
    std::vector<ReloadedCell> cells;
    bool hasNodeClasses=false;
};

double exportedDouble(double value) {
    std::ostringstream encoded;
    encoded<<value;
    std::istringstream decoded(encoded.str());
    double result=0.0;
    decoded>>result;
    return result;
}

std::optional<mesh::CellType> vtkCellType(int value) {
    if (value==10) return mesh::CellType::Tet4;
    return std::nullopt;
}

std::optional<mesh::CellType> mshCellType(int value) {
    if (value==1) return mesh::CellType::Line2;
    if (value==2) return mesh::CellType::Tri3;
    if (value==4) return mesh::CellType::Tet4;
    return std::nullopt;
}

std::optional<ReloadedMesh> reloadVtk(const std::string& path) {
    std::ifstream input(path);
    if (!input) return std::nullopt;
    std::string token;
    while (input>>token&&token!="POINTS") {}
    std::size_t nodeCount=0;
    std::string scalarType;
    if (token!="POINTS"||!(input>>nodeCount>>scalarType)) return std::nullopt;
    ReloadedMesh result;
    result.nodes.resize(nodeCount);
    result.nodeClasses.assign(nodeCount,mesh::Classification{});
    for (auto& point:result.nodes)
        if (!(input>>point.x>>point.y>>point.z)) return std::nullopt;
    std::size_t cellCount=0,listLength=0;
    if (!(input>>token>>cellCount>>listLength)||token!="CELLS")
        return std::nullopt;
    result.cells.resize(cellCount);
    std::size_t observedLength=0;
    for (auto& cell:result.cells) {
        std::size_t count=0;
        if (!(input>>count)) return std::nullopt;
        observedLength+=count+1;
        cell.nodes.resize(count);
        for (auto& node:cell.nodes)
            if (!(input>>node)) return std::nullopt;
    }
    if (observedLength!=listLength) return std::nullopt;
    std::size_t typeCount=0;
    if (!(input>>token>>typeCount)||token!="CELL_TYPES"||
        typeCount!=cellCount)
        return std::nullopt;
    for (auto& cell:result.cells) {
        int type=0;
        if (!(input>>type)) return std::nullopt;
        const auto mapped=vtkCellType(type);
        if (!mapped) return std::nullopt;
        cell.type=*mapped;
        cell.classification.dim=3;
    }
    std::size_t dataCount=0;
    std::string name,integerType,components,lookup,table;
    if (!(input>>token>>dataCount)||token!="CELL_DATA"||
        dataCount!=cellCount||
        !(input>>token>>name>>integerType>>components)||token!="SCALARS"||
        name!="entity_tag"||
        !(input>>lookup>>table)||lookup!="LOOKUP_TABLE")
        return std::nullopt;
    for (auto& cell:result.cells)
        if (!(input>>cell.classification.tag)) return std::nullopt;
    return result;
}

std::optional<ReloadedMesh> reloadMsh(const std::string& path) {
    std::ifstream input(path);
    if (!input) return std::nullopt;
    ReloadedMesh result;
    std::string token;
    while (input>>token&&token!="$Nodes") {}
    std::size_t blockCount=0,nodeCount=0;
    std::uint64_t minimumTag=0,maximumTag=0;
    if (token!="$Nodes"||
        !(input>>blockCount>>nodeCount>>minimumTag>>maximumTag)||
        minimumTag!=1||maximumTag!=nodeCount)
        return std::nullopt;
    result.nodes.resize(nodeCount);
    result.nodeClasses.resize(nodeCount);
    result.hasNodeClasses=true;
    std::vector<bool> nodeSeen(nodeCount,false);
    for (std::size_t block=0;block<blockCount;++block) {
        int dimension=0,tag=0,parametric=0;
        std::size_t count=0;
        if (!(input>>dimension>>tag>>parametric>>count)||parametric!=0)
            return std::nullopt;
        std::vector<std::uint64_t> tags(count);
        for (auto& value:tags)
            if (!(input>>value)||value==0||value>nodeCount)
                return std::nullopt;
        for (const auto value:tags) {
            const auto index=static_cast<std::size_t>(value-1);
            if (nodeSeen[index]) return std::nullopt;
            nodeSeen[index]=true;
            result.nodeClasses[index]={dimension,tag};
            auto& point=result.nodes[index];
            if (!(input>>point.x>>point.y>>point.z)) return std::nullopt;
        }
    }
    if (!std::all_of(nodeSeen.begin(),nodeSeen.end(),[](bool v){return v;})||
        !(input>>token)||token!="$EndNodes")
        return std::nullopt;
    while (input>>token&&token!="$Elements") {}
    std::size_t elementBlocks=0,elementCount=0;
    std::uint64_t minimumElement=0,maximumElement=0;
    if (token!="$Elements"||
        !(input>>elementBlocks>>elementCount>>minimumElement>>maximumElement)||
        (elementCount>0&&(minimumElement!=1||maximumElement!=elementCount)))
        return std::nullopt;
    result.cells.reserve(elementCount);
    std::set<std::uint64_t> elementTags;
    for (std::size_t block=0;block<elementBlocks;++block) {
        int dimension=0,tag=0,typeCode=0;
        std::size_t count=0;
        if (!(input>>dimension>>tag>>typeCode>>count)) return std::nullopt;
        const auto type=mshCellType(typeCode);
        if (!type) return std::nullopt;
        const auto nodePerCell=static_cast<std::size_t>(mesh::nodesPerCell(*type));
        for (std::size_t row=0;row<count;++row) {
            std::uint64_t elementTag=0;
            if (!(input>>elementTag)||!elementTags.insert(elementTag).second)
                return std::nullopt;
            ReloadedCell cell;
            cell.type=*type;
            cell.classification={dimension,tag};
            cell.nodes.resize(nodePerCell);
            for (auto& node:cell.nodes) {
                std::uint64_t nodeTag=0;
                if (!(input>>nodeTag)||nodeTag==0||nodeTag>nodeCount)
                    return std::nullopt;
                node=static_cast<Index>(nodeTag-1);
            }
            result.cells.push_back(std::move(cell));
        }
    }
    if (result.cells.size()!=elementCount||
        !(input>>token)||token!="$EndElements")
        return std::nullopt;
    return result;
}

ReloadedMesh exportProjection(const mesh::Mesh& source,bool volumeOnly,
                              bool nodeClasses,bool fullPrecision=false) {
    ReloadedMesh result;
    result.hasNodeClasses=nodeClasses;
    result.nodes.reserve(static_cast<std::size_t>(source.nodeCount()));
    result.nodeClasses.reserve(static_cast<std::size_t>(source.nodeCount()));
    for (Index i=0;i<source.nodeCount();++i) {
        const auto point=source.node(i);
        result.nodes.push_back(fullPrecision?point:Vec3{
            exportedDouble(point.x),exportedDouble(point.y),exportedDouble(point.z)});
        result.nodeClasses.push_back(
            nodeClasses?source.nodeClass(i):mesh::Classification{});
    }
    for (Index i=0;i<source.cellCount();++i) {
        const auto cell=source.cell(i);
        if (volumeOnly&&cell.type!=mesh::CellType::Tet4) continue;
        ReloadedCell row;
        row.type=cell.type;
        row.classification=cell.onEntity;
        const auto count=mesh::nodesPerCell(cell.type);
        row.nodes.assign(source.cellNodes(i),source.cellNodes(i)+count);
        result.cells.push_back(std::move(row));
    }
    return result;
}

mesh::Mesh normalizedReloadMesh(ReloadedMesh value) {
    mesh::Mesh result;
    result.reserveAppend(value.nodes.size(),value.cells.size(),
        std::accumulate(value.cells.begin(),value.cells.end(),std::size_t{0},
            [](std::size_t total,const ReloadedCell& cell) {
                return total+cell.nodes.size();
            }));
    for (std::size_t i=0;i<value.nodes.size();++i)
        result.addNode(value.nodes[i],value.nodeClasses[i]);
    std::sort(value.cells.begin(),value.cells.end());
    for (const auto& cell:value.cells)
        result.addCell(cell.type,cell.nodes,cell.classification);
    return result;
}

bool sameReloadNodes(const ReloadedMesh& a,const ReloadedMesh& b) noexcept {
    if (a.nodes.size()!=b.nodes.size()) return false;
    for (std::size_t i=0;i<a.nodes.size();++i) {
        const auto& x=a.nodes[i];
        const auto& y=b.nodes[i];
        if (x.x!=y.x||x.y!=y.y||x.z!=y.z) return false;
    }
    return true;
}

bool sameReloadNodeClasses(const ReloadedMesh& a,
                           const ReloadedMesh& b) noexcept {
    if (a.nodeClasses.size()!=b.nodeClasses.size()) return false;
    for (std::size_t i=0;i<a.nodeClasses.size();++i)
        if (a.nodeClasses[i].dim!=b.nodeClasses[i].dim||
            a.nodeClasses[i].tag!=b.nodeClasses[i].tag)
            return false;
    return true;
}

std::set<std::array<Index,3>> tetBoundary(const ReloadedMesh& value) {
    std::map<std::array<Index,3>,std::size_t> incidence;
    static constexpr std::array<std::array<std::size_t,3>,4> faces{{
        {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
    for (const auto& cell:value.cells) {
        if (cell.type!=mesh::CellType::Tet4||cell.nodes.size()!=4) continue;
        for (const auto& local:faces) {
            std::array<Index,3> key{
                cell.nodes[local[0]],cell.nodes[local[1]],cell.nodes[local[2]]};
            std::sort(key.begin(),key.end());
            ++incidence[key];
        }
    }
    std::set<std::array<Index,3>> result;
    for (const auto& [key,count]:incidence)
        if (count==1) result.insert(key);
    return result;
}

std::set<std::array<Index,3>> prescribedTriangles(const ReloadedMesh& value) {
    std::set<std::array<Index,3>> result;
    for (const auto& cell:value.cells) {
        if (cell.type!=mesh::CellType::Tri3||cell.nodes.size()!=3) continue;
        std::array<Index,3> key{cell.nodes[0],cell.nodes[1],cell.nodes[2]};
        std::sort(key.begin(),key.end());
        result.insert(key);
    }
    return result;
}

void auditedExport() {
    TEST("T4-A audited improved Tet4 exports VTK and MSH");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureProvider(*model,provider);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.4);
    auto request=boundedRequest(provider,1);
    mesh::Mesh output;
    const auto outcome=algo::meshModelQualityTet4(
        *model,spec,request,output,provider);
    CHECK(outcome.published());
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.committedOperations==1);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->after.minimum>
              outcome.retainedResult->before.minimum);
    const auto base=std::filesystem::temp_directory_path()/
        ("paacs_t4a_audited_test_"+std::to_string(::getpid()));
    const auto vtk=base.string()+".vtk";
    const auto msh=base.string()+".msh";
    io::writeVtk(output,vtk,true);
    io::writeMsh(output,msh);
    const auto vtkReload=reloadVtk(vtk);
    const auto mshReload=reloadMsh(msh);
    CHECK(vtkReload.has_value());
    CHECK(mshReload.has_value());
    const auto vtkExpected=exportProjection(output,true,false,true);
    const auto mshExpected=exportProjection(output,false,true);
    if (vtkReload) {
        CHECK(sameReloadNodes(*vtkReload,vtkExpected));
        auto actualCells=vtkReload->cells;
        auto expectedCells=vtkExpected.cells;
        std::sort(actualCells.begin(),actualCells.end());
        std::sort(expectedCells.begin(),expectedCells.end());
        CHECK(actualCells==expectedCells);
        CHECK(tetBoundary(*vtkReload)==prescribedTriangles(mshExpected));
        const auto actualMesh=normalizedReloadMesh(*vtkReload);
        const auto expectedMesh=normalizedReloadMesh(vtkExpected);
        CHECK(mesh::meshHash(actualMesh)==mesh::meshHash(expectedMesh));
        const auto actualSemantic=mesh::semanticHashEx(actualMesh);
        const auto expectedSemantic=mesh::semanticHashEx(expectedMesh);
        CHECK(actualSemantic.canonical&&expectedSemantic.canonical);
        CHECK(actualSemantic.value==expectedSemantic.value);
    }
    if (mshReload) {
        CHECK(sameReloadNodes(*mshReload,mshExpected));
        CHECK(sameReloadNodeClasses(*mshReload,mshExpected));
        auto actualCells=mshReload->cells;
        auto expectedCells=mshExpected.cells;
        std::sort(actualCells.begin(),actualCells.end());
        std::sort(expectedCells.begin(),expectedCells.end());
        CHECK(actualCells==expectedCells);
        CHECK(tetBoundary(*mshReload)==prescribedTriangles(*mshReload));
        const auto actualMesh=normalizedReloadMesh(*mshReload);
        const auto expectedMesh=normalizedReloadMesh(mshExpected);
        CHECK(mesh::meshHash(actualMesh)==mesh::meshHash(expectedMesh));
        const auto actualSemantic=mesh::semanticHashEx(actualMesh);
        const auto expectedSemantic=mesh::semanticHashEx(expectedMesh);
        CHECK(actualSemantic.canonical&&expectedSemantic.canonical);
        CHECK(actualSemantic.value==expectedSemantic.value);
    }
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->finalStrictHash==mesh::meshHash(output));
    std::filesystem::remove(vtk);
    std::filesystem::remove(msh);
}

} // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* t4aTests(std::size_t& count) {
    static const Entry tests[] = {
        {"t4a.quality.exact",exactQualityMetric},
        {"t4a.route.preflight_sizing",preflightAndSizingBoundary},
        {"t4a.route.report_capacity",reportOnlyAndCapacity},
        {"t4a.route.blocker_invalidation",blockerContinuationAndInvalidation},
        {"t4a.route.protected_blocker",protectedBlockerContinuation},
        {"t4a.route.budget_taxonomy",budgetsAndResidualTaxonomy},
        {"t4a.route.post_commit_quarantine",postCommitQuarantine},
        {"t4a.route.candidate_taxonomy",candidateTaxonomy},
        {"t4a.route.checkpoint_digests",checkpointDigests},
        {"t4a.route.determinism",completeDeterminism},
        {"t4a.metrics.public_scaling",publicOperationScaling},
        {"t4a.route.boundary_shapes",shapeBoundaryBijection},
        {"t4a.export.audited",auditedExport},
    };
    count=sizeof(tests)/sizeof(tests[0]);
    return tests;
}
} // namespace paacs_tests
