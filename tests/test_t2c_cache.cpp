// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/algo/T2CSurfaceCache.hpp"
#include "paacs/algo/FacetBarycentre.hpp"
#include "paacs/algo/FacetConflictResolver.hpp"
#include "test_harness.hpp"

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <vector>

namespace {
namespace algo=paacs::algo;
namespace cache=paacs::algo::t2c_detail;

class FixedSizing final : public cache::SurfaceCacheSizingEvaluator {
public:
    explicit FixedSizing(double value) : value_(value) {}
    double evaluate(const paacs::Vec3&,double) const override {
        ++calls; return value_;
    }
    mutable std::uint64_t calls=0;
private:
    double value_=0.0;
};

class OrderedSizing final : public cache::SurfaceCacheSizingEvaluator {
public:
    double evaluate(const paacs::Vec3& point,double) const override {
        if(calls<seen.size())seen[calls]=point.x;
        ++calls;
        return point.x*11.0;
    }
    mutable std::array<double,4> seen{};
    mutable std::size_t calls=0;
};

std::shared_ptr<const cache::CachedSurfaceFaceArtifact> entry(
    std::uint64_t face,std::uint8_t key) {
    auto result=std::make_shared<cache::CachedSurfaceFaceArtifact>();
    result->facePersistentId=face;
    result->faceStamp={face,1};
    result->authorityBytes={key};
    result->complete=true;
    result->audit.workspaceReady=true;
    result->metrics.facePersistentId=face;
    result->metrics.construction.adjacencyAuditPassed=true;
    result->metrics.construction.runComplete=true;
    return result;
}

algo::ChartCertificate barycentreCertificate() {
    algo::ChartCertificate result;
    result.providerStamp={41,3};
    result.chart=7;
    result.certificateStamp={42,5};
    result.certifiedRegionBoundary={
        {-2.0,-2.0},{2.0,-2.0},{2.0,2.0},{-2.0,2.0}};
    result.chartPositiveToNormalSign=1;
    return result;
}

bool sameChartCertificate(const algo::ChartCertificate& a,
                          const algo::ChartCertificate& b) {
    if(a.providerStamp!=b.providerStamp||a.chart!=b.chart||
       a.certificateStamp!=b.certificateStamp||
       a.chartPositiveToNormalSign!=b.chartPositiveToNormalSign||
       a.regular!=b.regular||a.periodic!=b.periodic||a.seam!=b.seam||
       a.certifiedRegionBoundary.size()!=b.certifiedRegionBoundary.size())
        return false;
    for(std::size_t i=0;i<a.certifiedRegionBoundary.size();++i)
        if(a.certifiedRegionBoundary[i].u!=b.certifiedRegionBoundary[i].u||
           a.certifiedRegionBoundary[i].v!=b.certifiedRegionBoundary[i].v)
            return false;
    return true;
}

bool sameChartCertificate(
    const algo::FacetBarycentreCertificateView& a,
    const algo::ChartCertificate& b) {
    if(!(a.providerStamp==b.providerStamp)||a.chart!=b.chart||
       a.certificateStamp!=b.certificateStamp||
       a.chartPositiveToNormalSign!=b.chartPositiveToNormalSign||
       a.regular!=b.regular||a.periodic!=b.periodic||a.seam!=b.seam||
       a.certifiedRegionBoundary.size()!=b.certifiedRegionBoundary.size())
        return false;
    for(std::size_t i=0;i<a.certifiedRegionBoundary.size();++i)
        if(a.certifiedRegionBoundary[i].u!=b.certifiedRegionBoundary[i].u||
           a.certifiedRegionBoundary[i].v!=b.certifiedRegionBoundary[i].v)
            return false;
    return true;
}

bool sameChartCertificate(
    const algo::FacetBarycentreCertificateView& a,
    const algo::FacetBarycentreCertificateView& b) {
    if(!(a.providerStamp==b.providerStamp)||a.chart!=b.chart||
       a.certificateStamp!=b.certificateStamp||
       a.chartPositiveToNormalSign!=b.chartPositiveToNormalSign||
       a.regular!=b.regular||a.periodic!=b.periodic||a.seam!=b.seam||
       a.certifiedRegionBoundary.size()!=b.certifiedRegionBoundary.size())
        return false;
    for(std::size_t i=0;i<a.certifiedRegionBoundary.size();++i)
        if(a.certifiedRegionBoundary[i].u!=b.certifiedRegionBoundary[i].u||
           a.certifiedRegionBoundary[i].v!=b.certifiedRegionBoundary[i].v)
            return false;
    return true;
}

algo::FacetBarycentreRequest barycentreRequest() {
    algo::FacetBarycentreRequest result;
    result.faceStamp={40,9};
    result.chart=7;
    result.expectedProviderStamp={41,3};
    result.expectedCertificateStamp={42,5};
    result.vertices={{{30,7,{0.25,0.75}},
                      {10,7,{0.0,0.0}},
                      {20,7,{0.5,0.25}}}};
    return result;
}

void facetBarycentreCanonicalBits() {
    TEST("t2c facet barycentre has canonical exact bits");
    const auto certificate=barycentreCertificate();
    const auto request=barycentreRequest();
    const auto reference=algo::computeFacetBarycentre(request,certificate);
    CHECK(reference.ready());
    CHECK(reference.workUnits>0);
    CHECK(reference.canonicalVertices[0].publishedGlobalNode==10);
    CHECK(reference.canonicalVertices[1].publishedGlobalNode==20);
    CHECK(reference.canonicalVertices[2].publishedGlobalNode==30);
    CHECK(std::bit_cast<std::uint64_t>(reference.parameter.u)==
          UINT64_C(0x3fd0000000000000));
    CHECK(std::bit_cast<std::uint64_t>(reference.parameter.v)==
          UINT64_C(0x3fd5555555555555));
    CHECK(reference.identity.domain==
          algo::SurfaceNodeIdentityDomain::RefinedFacetInterior);
    CHECK(reference.identity.parentArtifact==request.faceStamp);
    CHECK(reference.identity.parentFacetGlobalNodes==
          (std::array<paacs::Index,3>{{10,20,30}}));
    CHECK(reference.identity.parentOrdinal==0);
    CHECK(reference.identity.parentSubordinal==0);
    CHECK(reference.identity.parameterUBits==
          UINT64_C(0x3fd0000000000000));
    CHECK(reference.identity.parameterVBits==
          UINT64_C(0x3fd5555555555555));
    CHECK(reference.identity.dyadicNumerator==0);
    CHECK(reference.identity.dyadicLevel==0);
    CHECK(algo::validSurfaceNodeIdentity(reference.identity));

    constexpr std::array<std::array<std::size_t,3>,6> orders{{
        {{0,1,2}},{{0,2,1}},{{1,0,2}},{{1,2,0}},{{2,0,1}},{{2,1,0}}}};
    for(const auto order:orders) {
        auto reversed=request;
        for(std::size_t i=0;i<3;++i)
            reversed.vertices[i]=request.vertices[order[i]];
        const auto replay=algo::computeFacetBarycentre(reversed,certificate);
        CHECK(replay.ready());
        CHECK(replay.canonicalVertices==reference.canonicalVertices);
        CHECK(std::bit_cast<std::uint64_t>(replay.parameter.u)==
              reference.identity.parameterUBits);
        CHECK(std::bit_cast<std::uint64_t>(replay.parameter.v)==
              reference.identity.parameterVBits);
        CHECK(replay.identity==reference.identity);
        CHECK(replay.workUnits==reference.workUnits);
    }

    for(std::uint64_t budget=0;budget<reference.workUnits;++budget) {
        const auto stopped=algo::computeFacetBarycentre(
            request,certificate,budget);
        CHECK(stopped.status==algo::FacetBarycentreStatus::WorkBudgetReached);
        CHECK(stopped.workUnits<=budget);
    }
    const auto exact=algo::computeFacetBarycentre(
        request,certificate,reference.workUnits);
    CHECK(exact.ready());
    CHECK(exact.workUnits==reference.workUnits);
    CHECK(exact.identity==reference.identity);
}

void facetBarycentrePreflight() {
    TEST("t2c facet barycentre rejects invalid plain inputs");
    const auto requireStatus=[](algo::FacetBarycentreRequest request,
                                algo::ChartCertificate certificate,
                                algo::FacetBarycentreStatus status) {
        CHECK(algo::computeFacetBarycentre(request,certificate).status==status);
    };
    const auto request=barycentreRequest();
    const auto certificate=barycentreCertificate();
    auto changed=request;
    changed.faceStamp={};
    requireStatus(changed,certificate,
                  algo::FacetBarycentreStatus::InvalidFaceIdentity);
    changed=request;changed.vertices[0].publishedGlobalNode=-1;
    requireStatus(changed,certificate,
                  algo::FacetBarycentreStatus::InvalidVertexIdentity);
    changed=request;changed.vertices[0].publishedGlobalNode=20;
    requireStatus(changed,certificate,
                  algo::FacetBarycentreStatus::DuplicateVertexIdentity);
    changed=request;changed.vertices[0].parameter.u=
        std::numeric_limits<double>::infinity();
    requireStatus(changed,certificate,
                  algo::FacetBarycentreStatus::NonfiniteParameter);
    changed=request;changed.vertices[0].chart=8;
    requireStatus(changed,certificate,
                  algo::FacetBarycentreStatus::ChartMismatch);
    auto changedCertificate=certificate;changedCertificate.chart=8;
    requireStatus(request,changedCertificate,
                  algo::FacetBarycentreStatus::ChartMismatch);
    changedCertificate=certificate;changedCertificate.periodic=true;
    requireStatus(request,changedCertificate,
                  algo::FacetBarycentreStatus::UnsupportedChart);
    changedCertificate=certificate;changedCertificate.seam=true;
    requireStatus(request,changedCertificate,
                  algo::FacetBarycentreStatus::UnsupportedChart);
    changedCertificate=certificate;changedCertificate.regular=false;
    requireStatus(request,changedCertificate,
                  algo::FacetBarycentreStatus::UnsupportedChart);
    changed=request;changed.expectedProviderStamp={};
    requireStatus(changed,certificate,
                  algo::FacetBarycentreStatus::InvalidCertificateAuthority);
    changed=request;changed.expectedCertificateStamp={};
    requireStatus(changed,certificate,
                  algo::FacetBarycentreStatus::InvalidCertificateAuthority);
    changedCertificate=certificate;changedCertificate.providerStamp={99,1};
    requireStatus(request,changedCertificate,
                  algo::FacetBarycentreStatus::InvalidCertificateAuthority);
    changedCertificate=certificate;changedCertificate.certificateStamp={};
    requireStatus(request,changedCertificate,
                  algo::FacetBarycentreStatus::InvalidCertificateAuthority);
    changedCertificate=certificate;
    changedCertificate.certifiedRegionBoundary.resize(2);
    requireStatus(request,changedCertificate,
                  algo::FacetBarycentreStatus::InvalidChartBoundary);
    changedCertificate=certificate;
    changedCertificate.certifiedRegionBoundary[0].u=
        std::numeric_limits<double>::quiet_NaN();
    requireStatus(request,changedCertificate,
                  algo::FacetBarycentreStatus::InvalidChartBoundary);
    changedCertificate=certificate;
    changedCertificate.certifiedRegionBoundary={
        {0.0,0.0},{0.1,0.0},{0.1,0.1},{0.0,0.1}};
    requireStatus(request,changedCertificate,
                  algo::FacetBarycentreStatus::OutsideCertifiedRegion);
    changed=request;
    changed.vertices={{{10,7,{0.0,0.0}},
                       {20,7,{0.5,0.5}},
                       {30,7,{1.0,1.0}}}};
    requireStatus(changed,certificate,
                  algo::FacetBarycentreStatus::DegenerateFacet);

    changed=request;
    changed.vertices={{{10,7,{-1.0,0.0}},
                       {20,7,{-0.0,1.0}},
                       {30,7,{1.0,0.0}}}};
    const auto zero=algo::computeFacetBarycentre(changed,certificate);
    CHECK(zero.ready());
    CHECK(std::bit_cast<std::uint64_t>(zero.parameter.u)==0);
    CHECK(zero.identity.parameterUBits==0);
}

void facetBarycentreRealLayerGoldens() {
    TEST("t2c facet barycentre pins both real layer SCC facets");
    struct Golden {
        std::uint64_t face;
        std::uint64_t chart;
        std::array<paacs::Index,3> globals;
    };
    constexpr std::array<Golden,2> goldens{{
        {40,40,{{77,169,171}}},
        {41,41,{{89,180,183}}}}};
    const std::array<paacs::Vec2,3> parameters{{
        {0.0,0.33333333333333331},
        {0.21458333333333329,0.60416666666666652},
        {0.29947916666666663,0.31556488160479113}}};
    constexpr std::array<std::array<std::size_t,3>,6> orders{{
        {{0,1,2}},{{0,2,1}},{{1,0,2}},{{1,2,0}},{{2,0,1}},{{2,1,0}}}};
    for(const auto& golden:goldens) {
        algo::FacetBarycentreRequest request;
        request.faceStamp={golden.face,1};
        request.chart=golden.chart;
        for(std::size_t i=0;i<3;++i)
            request.vertices[i]={golden.globals[i],golden.chart,parameters[i]};
        algo::ChartCertificate certificate;
        certificate.chart=golden.chart;
        certificate.providerStamp={900+golden.face,1};
        certificate.certificateStamp={1000+golden.face,1};
        certificate.certifiedRegionBoundary={
            {-1.0,-1.0},{2.0,-1.0},{2.0,2.0},{-1.0,2.0}};
        certificate.chartPositiveToNormalSign=1;
        request.expectedProviderStamp=certificate.providerStamp;
        request.expectedCertificateStamp=certificate.certificateStamp;
        const auto reference=algo::computeFacetBarycentre(request,certificate);
        CHECK(reference.ready());
        CHECK(reference.identity.parentArtifact==request.faceStamp);
        CHECK(reference.identity.parentFacetGlobalNodes==golden.globals);
        CHECK(reference.identity.parentOrdinal==0);
        CHECK(reference.identity.parentSubordinal==0);
        CHECK(reference.identity.dyadicNumerator==0);
        CHECK(reference.identity.dyadicLevel==0);
        CHECK(reference.identity.parameterUBits==
              UINT64_C(0x3fc5eeeeeeeeeeee));
        CHECK(reference.identity.parameterVBits==
              UINT64_C(0x3fdabb67ae8584ca));
        for(const auto order:orders) {
            auto reversed=request;
            for(std::size_t i=0;i<3;++i)
                reversed.vertices[i]=request.vertices[order[i]];
            const auto replay=algo::computeFacetBarycentre(
                reversed,certificate);
            CHECK(replay.ready());
            CHECK(replay.canonicalVertices==reference.canonicalVertices);
            CHECK(replay.identity==reference.identity);
        }
    }
}

void immutableSlotReplacement() {
    TEST("t2c cache immutable slots and checked budgets");
    cache::SurfaceFaceCacheLimits limits;
    limits.maxFaces=2;
    limits.maxCachedNodeRows=2;
    limits.maxCachedCellRows=2;
    limits.maxCachedConstraintRows=2;
    limits.maxCachedMappingRows=4;
    cache::SurfaceFaceCacheSession session(limits);
    const auto first=entry(10,1);
    const auto second=entry(20,2);
    const std::array<algo::FacePersistentId,2> faces{10,20};
    std::vector<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>> next{
        first,second};
    CHECK(session.prepareReplacement(std::move(next),faces));
    session.commitPrepared();
    const std::array<std::uint8_t,1> key1{1},key2{2},bad{3};
    CHECK(session.find(10,key1)==first);
    CHECK(session.find(20,key2)==second);
    CHECK(!session.find(10,bad));
    CHECK(session.metrics().cacheCommits==1);

    auto tooLarge=std::make_shared<cache::CachedSurfaceFaceArtifact>();
    *tooLarge=*first;
    tooLarge->localNodes.resize(3);
    std::vector<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>> badNext{
        tooLarge,second};
    CHECK(!session.prepareReplacement(std::move(badNext),faces));
    CHECK(session.diagnostic().status==
          cache::SurfaceFaceCacheStatus::Incomplete);
    CHECK(session.find(10,key1)==first);
}

void sizingPrefixReplay() {
    TEST("t2c cache sizing prefix replay and authority change");
    cache::SurfaceFaceCacheSession session;
    cache::CachedSurfaceFaceArtifact stored;
    stored.facePersistentId=77;
    stored.complete=true;
    stored.sizingSamples={
        {77,0,{1.0,2.0,3.0},0.25,2.0,0},
        {77,1,{4.0,5.0,6.0},0.5,3.0,1}};
    session.beginAttempt();
    CHECK(session.recordSizing(77,{1.0,2.0,3.0},0.25,2.0));
    FixedSizing same(3.0);
    CHECK(session.verifyAndReplaySizing(stored,same));
    CHECK(same.calls==1);
    CHECK(session.metrics().sizingReplayCalls==1);
    const std::array<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>,1>
        complete{std::make_shared<cache::CachedSurfaceFaceArtifact>(stored)};
    CHECK(session.verifyFinalSizing(complete,same));
    CHECK(same.calls==1); // the final checkpoint is comparison-only

    session.beginAttempt();
    CHECK(session.recordSizing(77,{1.0,2.0,3.0},0.25,2.0));
    FixedSizing changed(4.0);
    CHECK(!session.verifyAndReplaySizing(stored,changed));
    CHECK(session.diagnostic().reason==
          cache::SurfaceFaceCacheReason::AuthorityChanged);
    CHECK(session.diagnostic().phase==
          cache::SurfaceFaceCachePhase::SizingReplay);

    session.beginAttempt();
    auto wrongStatus=stored;
    wrongStatus.sizingSamples[0].status=
        algo::ProductionSurfaceReason::SizingEvaluationFailed;
    CHECK(session.recordSizing(77,{1.0,2.0,3.0},0.25,2.0));
    FixedSizing statusEvaluator(3.0);
    CHECK(!session.verifyAndReplaySizing(wrongStatus,statusEvaluator));
    CHECK(session.diagnostic().reason==
          cache::SurfaceFaceCacheReason::AuthorityChanged);

    cache::SurfaceFaceCacheSession orderedSession;
    auto orderedEntry=std::make_shared<cache::CachedSurfaceFaceArtifact>();
    orderedEntry->facePersistentId=88;
    orderedEntry->complete=true;
    orderedEntry->sizingSamples={
        {88,0,{1.0,0.0,0.0},0.0,11.0,0},
        {88,1,{2.0,0.0,0.0},0.0,22.0,1},
        {88,2,{3.0,0.0,0.0},0.0,33.0,2}};
    orderedSession.beginAttempt();
    CHECK(orderedSession.recordSizing(88,{1.0,0.0,0.0},0.0,11.0));
    OrderedSizing ordered;
    CHECK(orderedSession.verifyAndReplaySizing(*orderedEntry,ordered));
    CHECK(ordered.calls==2);
    CHECK(ordered.seen[0]==2.0);
    CHECK(ordered.seen[1]==3.0);
    const std::array<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>,1>
        orderedEntries{orderedEntry};
    CHECK(orderedSession.verifyFinalSizing(orderedEntries,ordered));
    CHECK(ordered.calls==2);
}

void structuralKeyCannotLaunderSizingOrFlatIndices() {
    TEST("t2c cache structural key excludes sizing outputs and flat indices");
    algo::SurfaceFaceRequest request;
    request.faceStamp={77,3};
    request.chart=77;
    request.faceClass={2,7};
    request.desiredNormalSign=1;
    request.expectedGeometryStamp={9,4};
    request.expectedChartCertificateStamp={10,4};
    request.expectedCandidateProviderStamp={11,4};
    request.physicalPolicy.minPhysicalDoubleArea=1.0;
    request.physicalPolicy.maxVertexMappingError=2.0;
    request.physicalPolicy.maxMappingDeviation=3.0;
    request.physicalPolicy.minNormalAlignment=0.8;
    request.physicalPolicy.maxArtifactPositionError=4.0;
    algo::SharedCurveNodeKey key;
    key.curvePersistentId=91;
    key.curveRevision=5;
    key.nodeOrdinal=2;
    algo::SurfaceChain chain;
    chain.chainPersistentId=100;
    chain.role=algo::SurfaceChainRole::OuterBoundary;
    chain.classification={1,3};
    chain.uses.push_back({1,key});
    request.outer.pieces.push_back(chain);
    paacs::mesh::Mesh destination;
    (void)destination.addNode({2.0,3.0,4.0},{0,1});
    (void)destination.addNode({2.0,3.0,4.0},{0,1});
    algo::SharedCurveNodeBinding binding{key,0,{0,1}};
    std::map<algo::SharedCurveNodeKey,algo::SharedCurveNodeBinding> ledger;
    ledger.emplace(key,binding);
    algo::ChartCertificate certificate;
    certificate.providerStamp={9,4};
    certificate.chart=77;
    certificate.certificateStamp={10,4};
    certificate.certifiedRegionBoundary={{0.0,0.0},{1.0,0.0},{0.0,1.0}};
    const auto first=cache::encodeSurfaceFaceCacheAuthority(
        {70,2},request,certificate,ledger,destination,
        {0.0,0.0},{1.0,1.0},{1.0,1.0});
    request.physicalPolicy.minPhysicalDoubleArea=101.0;
    request.physicalPolicy.maxVertexMappingError=102.0;
    request.physicalPolicy.maxMappingDeviation=103.0;
    request.physicalPolicy.minNormalAlignment=0.7;
    request.physicalPolicy.maxArtifactPositionError=104.0;
    ledger.begin()->second.existingGlobalIndex=1;
    const auto second=cache::encodeSurfaceFaceCacheAuthority(
        {70,2},request,certificate,ledger,destination,
        {0.0,0.0},{1.0,1.0},{1.0,1.0});
    CHECK(first==second);

    cache::SurfaceFaceCacheSession session;
    auto stored=entry(77,0);
    auto mutableStored=std::make_shared<cache::CachedSurfaceFaceArtifact>(
        *stored);
    mutableStored->authorityBytes=first;
    mutableStored->sizingSamples={{77,0,{0.5,0.5,0.0},0.0,2.0,0}};
    std::vector<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>> slots{
        mutableStored};
    const std::array<algo::FacePersistentId,1> faces{77};
    CHECK(session.prepareReplacement(std::move(slots),faces));
    session.commitPrepared();
    session.beginAttempt();
    CHECK(session.recordSizing(77,{0.5,0.5,0.0},0.0,3.0));
    CHECK(session.find(77,second)==mutableStored);
    FixedSizing stable(2.0);
    CHECK(!session.verifyAndReplaySizing(*mutableStored,stable));
    CHECK(session.diagnostic().reason==
          cache::SurfaceFaceCacheReason::AuthorityChanged);
}

std::vector<std::uint8_t> invariantBytes(
    const std::vector<std::uint8_t>& authority) {
    CHECK(authority.size()>=sizeof(std::uint32_t)+sizeof(std::uint64_t));
    std::uint64_t count=0;
    std::memcpy(&count,authority.data()+sizeof(std::uint32_t),sizeof(count));
    const auto first=sizeof(std::uint32_t)+sizeof(std::uint64_t);
    CHECK(count<=authority.size()-first);
    return {authority.begin()+static_cast<std::ptrdiff_t>(first),
            authority.begin()+static_cast<std::ptrdiff_t>(first+count)};
}

void invariantAuthorityCannotBecomeDirty() {
    TEST("t2c cache invariant authority mismatch is terminal");
    algo::SurfaceFaceRequest request;
    request.faceStamp={51,2}; request.chart=51; request.faceClass={2,4};
    request.desiredNormalSign=1;
    request.expectedGeometryStamp={700,3};
    request.expectedChartCertificateStamp={701,3};
    request.expectedCandidateProviderStamp={702,3};
    algo::ChartCertificate certificate;
    certificate.providerStamp={700,3}; certificate.chart=51;
    certificate.certificateStamp={701,3};
    certificate.regular=true;
    certificate.certifiedRegionBoundary={
        {0.0,0.0},{1.0,0.0},{0.0,1.0}};
    paacs::mesh::Mesh destination;
    const std::map<algo::SharedCurveNodeKey,
                   algo::SharedCurveNodeBinding> ledger;
    const auto accepted=cache::encodeSurfaceFaceCacheAuthority(
        {900,8},request,certificate,ledger,destination,
        {0.0,0.0},{1.0,1.0},{1.0,1.0});
    auto stored=std::make_shared<cache::CachedSurfaceFaceArtifact>();
    stored->facePersistentId=51; stored->faceStamp={51,2};
    stored->authorityBytes=accepted;
    stored->invariantAuthorityBytes=invariantBytes(accepted);
    stored->complete=true; stored->audit.workspaceReady=true;
    stored->metrics.facePersistentId=51;
    stored->metrics.construction.adjacencyAuditPassed=true;
    stored->metrics.construction.runComplete=true;
    cache::SurfaceFaceCacheSession session;
    const std::array<algo::FacePersistentId,1> faces{51};
    std::vector<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>> slots{
        stored};
    CHECK(session.prepareReplacement(std::move(slots),faces));
    session.commitPrepared(); session.beginAttempt();
    const std::array<algo::FacePersistentId,0> clean{};
    CHECK(session.configureExpectedDirty(faces,clean));

    ++certificate.providerStamp.revision;
    const auto changed=cache::encodeSurfaceFaceCacheAuthority(
        {900,8},request,certificate,ledger,destination,
        {0.0,0.0},{1.0,1.0},{1.0,1.0});
    CHECK(!session.find(51,changed));
    CHECK(session.diagnostic().status==cache::SurfaceFaceCacheStatus::Failed);
    CHECK(session.diagnostic().reason==
          cache::SurfaceFaceCacheReason::AuthorityChanged);
    CHECK(!session.noteLookup(51,true));
    CHECK(session.metrics().dirtyPropagationVisits==1);
}

void exactPredictedDirtyDependencies() {
    TEST("t2c cache rejects omitted and extra dirty dependencies");
    const std::array<algo::FacePersistentId,2> faces{10,20};
    const std::array<algo::FacePersistentId,1> first{10};

    cache::SurfaceFaceCacheSession omitted;
    omitted.beginAttempt();
    CHECK(omitted.configureExpectedDirty(faces,first));
    CHECK(!omitted.noteLookup(10,false));
    CHECK(omitted.diagnostic().reason==cache::SurfaceFaceCacheReason::Invalid);
    CHECK(omitted.metrics().dirtyPropagationVisits==2);

    const std::array<algo::FacePersistentId,0> none{};
    cache::SurfaceFaceCacheSession extra;
    extra.beginAttempt();
    CHECK(extra.configureExpectedDirty(faces,none));
    CHECK(!extra.noteLookup(10,true));
    CHECK(extra.diagnostic().reason==cache::SurfaceFaceCacheReason::Invalid);
    CHECK(extra.metrics().dirtyPropagationVisits==2);

    cache::SurfaceFaceCacheLimits limits;
    limits.maxDirtyFacesPerBatch=1;
    cache::SurfaceFaceCacheSession bounded(limits);
    const std::array<algo::FacePersistentId,2> both{10,20};
    bounded.beginAttempt();
    CHECK(!bounded.configureExpectedDirty(faces,both));
    CHECK(bounded.diagnostic().status==
          cache::SurfaceFaceCacheStatus::Incomplete);
}

void globalSizingOrderIsCanonical() {
    TEST("t2c cache pins cross-face sizing evaluation order");
    cache::SurfaceFaceCacheSession session;
    session.beginAttempt();
    CHECK(session.recordSizing(10,{1.0,0.0,0.0},0.0,1.0));
    CHECK(session.recordSizing(20,{2.0,0.0,0.0},0.0,2.0));
    CHECK(session.recordSizing(10,{3.0,0.0,0.0},0.0,3.0));
    auto a=entry(10,1),b=entry(20,2);
    auto ma=std::make_shared<cache::CachedSurfaceFaceArtifact>(*a);
    auto mb=std::make_shared<cache::CachedSurfaceFaceArtifact>(*b);
    ma->sizingSamples.assign(session.observedSizing(10).begin(),
                             session.observedSizing(10).end());
    mb->sizingSamples.assign(session.observedSizing(20).begin(),
                             session.observedSizing(20).end());
    const std::array<algo::FacePersistentId,2> faces{10,20};
    const std::array<algo::FacePersistentId,2> dirty{10,20};
    CHECK(session.configureExpectedDirty(faces,dirty));
    CHECK(session.noteLookup(10,true)); CHECK(session.noteLookup(20,true));
    std::vector<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>> slots{
        ma,mb};
    CHECK(session.prepareReplacement(std::move(slots),faces));
    session.commitPrepared(); session.beginAttempt();
    CHECK(session.recordSizing(20,{2.0,0.0,0.0},0.0,2.0)==false);
    CHECK(session.diagnostic().reason==
          cache::SurfaceFaceCacheReason::AuthorityChanged);
}

void completeFaceXorDiagnostic() {
    TEST("t2c cache diagnostic pins complete face XOR counts");
    cache::SurfaceFaceCacheSession session;
    std::array<algo::SurfaceCacheFaceReportRow,5> processed{};
    std::array<algo::FacePersistentId,5> unattempted{};
    cache::SurfaceCacheDiagnosticChannel channel{processed,unattempted};
    session.attachDiagnosticChannel(&channel);
    session.setDiagnostic(cache::SurfaceFaceCacheStatus::Failed,
        cache::SurfaceFaceCacheReason::FaceBuildStopped,
        cache::SurfaceFaceCachePhase::FaceBuild,true,20,2,3);
    const auto& row=session.diagnostic();
    CHECK(row.reportComplete);
    CHECK(row.failedFace==20);
    CHECK(row.processedFaces==2);
    CHECK(row.unattemptedFaces==3);
    CHECK(row.processedFaces+row.unattemptedFaces==5);
    CHECK(session.metrics().facesRejected==1);
    CHECK(channel.status==cache::SurfaceFaceCacheStatus::Failed);
    CHECK(channel.reason==cache::SurfaceFaceCacheReason::FaceBuildStopped);
    CHECK(channel.phase==cache::SurfaceFaceCachePhase::FaceBuild);
    CHECK(channel.cacheWorkUnits==session.attemptWorkUnits());
    CHECK(channel.cacheAccountedBytes==
          session.metrics().accountedByteHighWater);
    session.setDiagnostic(cache::SurfaceFaceCacheStatus::Failed,
        cache::SurfaceFaceCacheReason::FaceBuildStopped,
        cache::SurfaceFaceCachePhase::SetAudit,true,20,2,3);
    CHECK(session.metrics().facesRejected==1);
}

void exactCoreMetricsAndDirtyLimit() {
    TEST("t2c cache exact work classification and dirty limit");
    auto rich=std::make_shared<cache::CachedSurfaceFaceArtifact>();
    rich->facePersistentId=41;
    rich->faceStamp={41,1};
    rich->authorityBytes={7};
    rich->localNodes.resize(1);
    rich->localParams.resize(1);
    rich->localMappingParams.resize(1);
    rich->localCells.resize(1);
    rich->persistentEdges.resize(1);
    rich->complete=true;
    rich->audit.workspaceReady=true;
    rich->metrics.facePersistentId=41;
    rich->metrics.construction.adjacencyAuditPassed=true;
    rich->metrics.construction.runComplete=true;

    cache::SurfaceFaceCacheSession session;
    session.beginAttempt();
    CHECK(session.noteBuilt(*rich));
    CHECK(session.noteReused(*rich));
    CHECK(session.noteSetAudit(3));
    // Seven logical rows: entry, authority, node, param, mapping param,
    // cell and persistent constraint.  Dirty snapshot work is preparation;
    // only rehydration is replay.
    CHECK(session.metrics().preparationWorkUnits==7);
    CHECK(session.metrics().replayWorkUnits==7);
    CHECK(session.metrics().setAuditWorkUnits==3);
    CHECK(session.attemptWorkUnits()==17);
    CHECK(session.metrics().dirtyNodeRows==1);
    CHECK(session.metrics().cleanNodeRows==1);
    CHECK(session.metrics().dirtyCellRows==1);
    CHECK(session.metrics().cleanCellRows==1);
    CHECK(session.metrics().constraintRows==2);
    CHECK(session.metrics().copiedBytes>=sizeof(*rich));

    const std::array<algo::FacePersistentId,1> faces{41};
    std::vector<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>> next{
        rich};
    CHECK(session.prepareReplacement(std::move(next),faces));
    // Replacement validation and conservative overlap accounting are two
    // distinct seven-row passes.
    CHECK(session.metrics().preparationWorkUnits==21);
    CHECK(session.attemptWorkUnits()==31);
    session.commitPrepared();
    CHECK(session.metrics().highWaterSlots==1);
    CHECK(session.metrics().highWaterRows==7);
    CHECK(session.metrics().accountedByteHighWater>sizeof(session));

    cache::SurfaceFaceCacheSession bounded;
    bounded.beginAttempt();
    const std::array<algo::FacePersistentId,2> dependencyFaces{10,20};
    CHECK(bounded.configureExpectedDirty(dependencyFaces,dependencyFaces));
    CHECK(bounded.noteLookup(10,true));
    CHECK(bounded.noteLookup(20,true));
    CHECK(bounded.metrics().keyRowsVisited==2);
    CHECK(bounded.metrics().dirtyPropagationVisits==4);
    CHECK(bounded.metrics().facesRejected==0);
    CHECK(bounded.metrics().preparationWorkUnits==6);
}

void fourBatchReplacementIsLocalAndDeterministic() {
    TEST("t2c cache four prepared batches replace only one immutable owner");
    constexpr std::array<algo::FacePersistentId,4> faces{10,20,30,40};
    cache::SurfaceFaceCacheSession session;
    std::array<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>,4>
        owners{entry(10,1),entry(20,2),entry(30,3),entry(40,4)};
    std::array<std::uint8_t,4> keys{1,2,3,4};
    std::vector<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>>
        initial(owners.begin(),owners.end());
    CHECK(session.prepareReplacement(std::move(initial),faces));
    session.commitPrepared();

    for(std::size_t round=0;round<faces.size();++round) {
        session.beginAttempt();
        const std::array<algo::FacePersistentId,1> dirty{faces[round]};
        CHECK(session.configureExpectedDirty(faces,dirty));
        auto replacement=entry(faces[round],
            static_cast<std::uint8_t>(20+round));
        std::vector<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>>
            next;
        next.reserve(faces.size());
        for(std::size_t i=0;i<faces.size();++i) {
            const std::uint8_t query=i==round
                ? static_cast<std::uint8_t>(20+round):keys[i];
            const std::array<std::uint8_t,1> authority{query};
            const auto found=session.find(faces[i],authority);
            CHECK(session.noteLookup(faces[i],!found));
            if(i==round) {
                CHECK(!found);
                CHECK(session.noteBuilt(*replacement));
                next.push_back(replacement);
            } else {
                CHECK(found==owners[i]);
                CHECK(session.noteReused(*found));
                next.push_back(found);
            }
        }
        CHECK(session.prepareReplacement(std::move(next),faces));
        session.commitPrepared();
        owners[round]=replacement;
        keys[round]=static_cast<std::uint8_t>(20+round);
    }
    CHECK(session.metrics().cacheCommits==5);
    CHECK(session.metrics().facesBuilt==4);
    CHECK(session.metrics().facesReused==12);
    CHECK(session.metrics().keyRowsVisited==16);
    // Four canonical dependency visits per batch plus the one exact dirty
    // transition selected in each of the four batches.
    CHECK(session.metrics().dirtyPropagationVisits==20);
    CHECK(session.metrics().highWaterSlots==8);
}

void facetProposalSecuresCacheBeforePlacement() {
    TEST("t2c facet proposal secures dirty report overlay and slots");
    cache::SurfaceFaceCacheLimits limits;
    limits.maxFaces=2;
    limits.maxDirtyFacesPerBatch=2;
    cache::SurfaceFaceCacheSession session(limits);
    const auto first=entry(10,1);
    const auto second=entry(20,2);
    constexpr std::array<algo::FacePersistentId,2> faces{10,20};
    std::vector<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>> slots{
        first,second};
    CHECK(session.prepareReplacement(std::move(slots),faces));
    session.commitPrepared();
    const auto beforeWork=session.metrics().preparationWorkUnits;
    const auto beforeOwners=session.slotsForTesting();
    CHECK(beforeOwners.size()==2);
    CHECK(beforeOwners[0]==first);
    CHECK(beforeOwners[1]==second);

    std::array<algo::SurfaceCacheFaceReportRow,2> processed{};
    std::array<algo::FacePersistentId,2> unattempted{};
    cache::SurfaceCacheDiagnosticChannel channel{processed,unattempted};
    channel.processedUsed=2;
    channel.reportComplete=true;
    session.attachDiagnosticChannel(&channel);
    constexpr std::array<paacs::bnd::VolumeFacetKey,3> planned{{
        {10,0},{10,1},{20,0}}};
    cache::SurfaceFacetCacheProposalCapacities capacity;
    capacity.committedFacePointCapacity=3;
    capacity.preparedFacePointCapacity=3;
    capacity.committedLedgerCapacity=3;
    capacity.preparedLedgerCapacity=3;
    capacity.committedSemanticStateCapacity=3;
    capacity.preparedSemanticStateCapacity=3;
    cache::SurfaceFacetCacheProposal proposal;
    CHECK(session.prepareFacetProposal(faces,planned,capacity,proposal));
    CHECK(proposal.ready);
    CHECK(proposal.plannedRows==3);
    CHECK(proposal.dirtyFaceRows==2);
    CHECK(proposal.replacementSlotRows==2);
    CHECK(proposal.replacementSlotCapacity>=2);
    CHECK(proposal.workUnits==15);
    CHECK(session.metrics().preparationWorkUnits-beforeWork==15);
    CHECK(proposal.accountedBytes<=limits.maxSurfaceCacheAccountedBytes);
    const auto dirty=session.preparedFacetDirtyFaces();
    CHECK(dirty.size()==2);
    CHECK(dirty[0]==10);
    CHECK(dirty[1]==20);
    CHECK(session.slotsForTesting()[0]==first);
    CHECK(session.slotsForTesting()[1]==second);

    std::vector<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>>
        replacement;
    CHECK(session.acquireReplacementScratch(faces.size(),replacement));
    CHECK(replacement.empty());
    CHECK(replacement.capacity()>=faces.size());
    session.recycleReplacementScratch(std::move(replacement));
    session.abandonFacetProposal();
    CHECK(session.preparedFacetDirtyFaces().empty());
    CHECK(session.slotsForTesting()[0]==first);
    CHECK(session.slotsForTesting()[1]==second);
}

void facetProposalStopsAtEveryCapacityBeforePlacement() {
    TEST("t2c facet proposal stops at one-over capacity without owner change");
    constexpr std::array<algo::FacePersistentId,2> faces{10,20};
    constexpr std::array<paacs::bnd::VolumeFacetKey,3> planned{{
        {10,0},{10,1},{20,0}}};
    const auto first=entry(10,1);
    const auto second=entry(20,2);
    const auto install=[&](cache::SurfaceFaceCacheSession& session) {
        std::vector<std::shared_ptr<const cache::CachedSurfaceFaceArtifact>>
            slots{first,second};
        CHECK(session.prepareReplacement(std::move(slots),faces));
        session.commitPrepared();
    };

    cache::SurfaceFaceCacheLimits dirtyLimits;
    dirtyLimits.maxFaces=2;
    dirtyLimits.maxDirtyFacesPerBatch=1;
    cache::SurfaceFaceCacheSession dirty(dirtyLimits);
    install(dirty);
    std::array<algo::SurfaceCacheFaceReportRow,2> dirtyProcessed{};
    std::array<algo::FacePersistentId,2> dirtyUnattempted{};
    cache::SurfaceCacheDiagnosticChannel dirtyChannel{
        dirtyProcessed,dirtyUnattempted};
    dirtyChannel.processedUsed=2;
    dirtyChannel.reportComplete=true;
    dirty.attachDiagnosticChannel(&dirtyChannel);
    cache::SurfaceFacetCacheProposalCapacities exact;
    exact.committedFacePointCapacity=3;
    exact.preparedFacePointCapacity=3;
    exact.committedLedgerCapacity=3;
    exact.preparedLedgerCapacity=3;
    exact.committedSemanticStateCapacity=3;
    exact.preparedSemanticStateCapacity=3;
    cache::SurfaceFacetCacheProposal result;
    CHECK(!dirty.prepareFacetProposal(faces,planned,exact,result));
    CHECK(dirty.diagnostic().status==
          cache::SurfaceFaceCacheStatus::Incomplete);
    CHECK(dirty.diagnostic().reason==
          cache::SurfaceFaceCacheReason::BudgetReached);
    CHECK(dirty.slotsForTesting()[0]==first);
    CHECK(dirty.slotsForTesting()[1]==second);

    cache::SurfaceFaceCacheLimits ordinaryLimits;
    ordinaryLimits.maxFaces=2;
    std::array<cache::SurfaceFacetCacheProposalCapacities,6> oneUnder{};
    oneUnder.fill(exact);
    oneUnder[0].committedFacePointRows=1;
    oneUnder[1].preparedFacePointRows=1;
    oneUnder[2].committedLedgerRows=1;
    oneUnder[3].preparedLedgerRows=1;
    oneUnder[4].committedSemanticStateRows=1;
    oneUnder[5].preparedSemanticStateRows=1;
    for(const auto& capacity:oneUnder) {
        cache::SurfaceFaceCacheSession overlay(ordinaryLimits);
        install(overlay);
        std::array<algo::SurfaceCacheFaceReportRow,2> overlayProcessed{};
        std::array<algo::FacePersistentId,2> overlayUnattempted{};
        cache::SurfaceCacheDiagnosticChannel overlayChannel{
            overlayProcessed,overlayUnattempted};
        overlayChannel.processedUsed=2;
        overlayChannel.reportComplete=true;
        overlay.attachDiagnosticChannel(&overlayChannel);
        CHECK(!overlay.prepareFacetProposal(
            faces,planned,capacity,result));
        CHECK(overlay.diagnostic().status==
              cache::SurfaceFaceCacheStatus::Incomplete);
        CHECK(overlay.preparedFacetDirtyFaces().empty());
        CHECK(overlay.slotsForTesting()[0]==first);
        CHECK(overlay.slotsForTesting()[1]==second);
    }

    cache::SurfaceFaceCacheSession shortReport(ordinaryLimits);
    install(shortReport);
    std::array<algo::SurfaceCacheFaceReportRow,1> shortProcessed{};
    std::array<algo::FacePersistentId,2> shortUnattempted{};
    cache::SurfaceCacheDiagnosticChannel shortChannel{
        shortProcessed,shortUnattempted};
    shortReport.attachDiagnosticChannel(&shortChannel);
    CHECK(!shortReport.prepareFacetProposal(faces,planned,exact,result));
    CHECK(shortReport.diagnostic().status==
          cache::SurfaceFaceCacheStatus::Incomplete);
    CHECK(shortReport.slotsForTesting()[0]==first);
    CHECK(shortReport.slotsForTesting()[1]==second);

    cache::SurfaceFaceCacheSession invalid(ordinaryLimits);
    install(invalid);
    std::array<algo::SurfaceCacheFaceReportRow,2> invalidProcessed{};
    std::array<algo::FacePersistentId,2> invalidUnattempted{};
    cache::SurfaceCacheDiagnosticChannel invalidChannel{
        invalidProcessed,invalidUnattempted};
    invalid.attachDiagnosticChannel(&invalidChannel);
    constexpr std::array<paacs::bnd::VolumeFacetKey,2> reversed{{
        {20,0},{10,0}}};
    CHECK(!invalid.prepareFacetProposal(faces,reversed,exact,result));
    CHECK(invalid.diagnostic().status==
          cache::SurfaceFaceCacheStatus::Failed);
    CHECK(invalid.diagnostic().reason==
          cache::SurfaceFaceCacheReason::Invalid);
    CHECK(invalid.slotsForTesting()[0]==first);
    CHECK(invalid.slotsForTesting()[1]==second);
}

std::unique_ptr<algo::SurfaceFaceResult> retainedFacetFace(
    algo::FacePersistentId faceId) {
    auto face=std::make_unique<algo::SurfaceFaceResult>();
    face->faceStamp={faceId,9};
    face->chart=7;
    face->faceClass={2,static_cast<int>(faceId-1)};
    face->chartCertificate=barycentreCertificate();
    face->expectedGeometryStamp=face->chartCertificate.providerStamp;
    face->localNodes.resize(4);
    face->localMappingParams={
        {0.25,0.75},{0.0,0.0},{0.75,0.5},{0.5,0.25}};
    face->publishedGlobalByLocalNode={30,10,40,20};
    // Deliberately opposite normalized order: local cell zero has canonical
    // triple (20,30,40), while local cell one has (10,20,30).
    face->localCells={
        {{{0,2,3}},{2,static_cast<int>(faceId-1)}},
        {{{1,3,0}},{2,static_cast<int>(faceId-1)}}};
    face->publishedGlobalCellByLocalCell={700,701};
    return face;
}

algo::FacetBarycentreResult childFromResolvedFacet(
    const algo::FacetConflictResolveResult& resolved) {
    algo::FacetBarycentreRequest request;
    request.faceStamp=resolved.faceStamp;
    request.chart=resolved.chart;
    request.expectedProviderStamp=resolved.chartCertificate.providerStamp;
    request.expectedCertificateStamp=
        resolved.chartCertificate.certificateStamp;
    request.vertices=resolved.vertices;
    return algo::computeFacetBarycentre(
        request,resolved.chartCertificate,
        std::numeric_limits<std::uint64_t>::max());
}

void facetConflictResolverUsesNormalizedRank() {
    TEST("t2c facet resolver uses canonical rank not local cell index");
    algo::SurfaceSetResult retained;
    retained.faces.push_back(retainedFacetFace(41));
    retained.faces.push_back(retainedFacetFace(40));
    const paacs::bnd::VolumeFacetKey key{40,0};
    std::array<algo::FacetConflictRankedFacetScratch,2> scratch{};
    const auto reference=algo::resolveFacetConflict(
        retained,key,scratch,std::numeric_limits<std::uint64_t>::max());
    CHECK(reference.resolved());
    CHECK(reference.reason==algo::FacetConflictResolveReason::None);
    CHECK(reference.workUnits>0);
    CHECK(reference.faceStamp==(algo::ArtifactStamp{40,9}));
    CHECK(reference.faceTag==39);
    CHECK(reference.chart==7);
    CHECK(sameChartCertificate(reference.chartCertificate,
                               barycentreCertificate()));
    CHECK(reference.vertices[0].publishedGlobalNode==10);
    CHECK(reference.vertices[0].parameter.u==0.0);
    CHECK(reference.vertices[0].parameter.v==0.0);
    CHECK(reference.vertices[1].publishedGlobalNode==20);
    CHECK(reference.vertices[1].parameter.u==0.5);
    CHECK(reference.vertices[1].parameter.v==0.25);
    CHECK(reference.vertices[2].publishedGlobalNode==30);
    CHECK(reference.vertices[2].parameter.u==0.25);
    CHECK(reference.vertices[2].parameter.v==0.75);

    std::swap(retained.faces[1]->localCells[0],
              retained.faces[1]->localCells[1]);
    std::swap(retained.faces[1]->publishedGlobalCellByLocalCell[0],
              retained.faces[1]->publishedGlobalCellByLocalCell[1]);
    std::reverse(retained.faces.begin(),retained.faces.end());
    const auto replay=algo::resolveFacetConflict(
        retained,key,scratch,std::numeric_limits<std::uint64_t>::max());
    CHECK(replay.resolved());
    CHECK(replay.workUnits==reference.workUnits);
    CHECK(replay.faceStamp==reference.faceStamp);
    CHECK(replay.faceTag==reference.faceTag);
    CHECK(replay.chart==reference.chart);
    CHECK(sameChartCertificate(replay.chartCertificate,
                               reference.chartCertificate));
    CHECK(replay.vertices==reference.vertices);

    algo::SurfaceSetOutcome outcome;
    outcome.retainedArtifact=std::make_unique<algo::SurfaceSetResult>(
        std::move(retained));
    const auto throughOutcome=algo::resolveFacetConflict(
        outcome,key,scratch,std::numeric_limits<std::uint64_t>::max());
    CHECK(throughOutcome.resolved());
    CHECK(throughOutcome.workUnits==reference.workUnits);
    CHECK(throughOutcome.vertices==reference.vertices);

    for(std::uint64_t budget=0;budget<reference.workUnits;++budget) {
        const auto stopped=algo::resolveFacetConflict(
            *outcome.retainedArtifact,key,scratch,budget);
        CHECK(stopped.status==algo::FacetConflictResolveStatus::Unsupported);
        CHECK(stopped.reason==
              algo::FacetConflictResolveReason::WorkBudgetReached);
        CHECK(stopped.workUnits<=budget);
    }
    const auto exact=algo::resolveFacetConflict(
        *outcome.retainedArtifact,key,scratch,reference.workUnits);
    CHECK(exact.resolved());
    CHECK(exact.workUnits==reference.workUnits);
    const auto noScratch=algo::resolveFacetConflict(
        *outcome.retainedArtifact,key,
        std::span<algo::FacetConflictRankedFacetScratch>{},
        std::numeric_limits<std::uint64_t>::max());
    CHECK(noScratch.status==algo::FacetConflictResolveStatus::Unsupported);
    CHECK(noScratch.reason==
          algo::FacetConflictResolveReason::ScratchCapacityReached);
}

void facetConflictResolverRejectsMissingAmbiguousAndCorrupt() {
    TEST("t2c facet resolver rejects nonunique and corrupt retained mapping");
    algo::SurfaceSetOutcome empty;
    const auto missingArtifact=algo::resolveFacetConflict(
        empty,paacs::bnd::VolumeFacetKey{40,0});
    CHECK(missingArtifact.status==algo::FacetConflictResolveStatus::Invalid);
    CHECK(missingArtifact.reason==
          algo::FacetConflictResolveReason::MissingRetainedArtifact);
    CHECK(missingArtifact.workUnits==0);

    algo::SurfaceSetResult retained;
    retained.faces.push_back(retainedFacetFace(40));
    auto result=algo::resolveFacetConflict(
        retained,paacs::bnd::VolumeFacetKey{41,0});
    CHECK(result.status==algo::FacetConflictResolveStatus::Unsupported);
    CHECK(result.reason==algo::FacetConflictResolveReason::FaceNotFound);
    CHECK(result.workUnits==1);
    result=algo::resolveFacetConflict(
        retained,paacs::bnd::VolumeFacetKey{40,2});
    CHECK(result.status==algo::FacetConflictResolveStatus::Unsupported);
    CHECK(result.reason==
          algo::FacetConflictResolveReason::FacetOrdinalOutOfRange);
    CHECK(result.workUnits>1);

    {
        algo::SurfaceSetResult compatible;
        compatible.faces.push_back(retainedFacetFace(40));
        compatible.faces.push_back(retainedFacetFace(40));
        compatible.faces.back()->chart=8;
        compatible.faces.back()->chartCertificate.chart=8;
        compatible.faces.back()->localMappingParams[1]={-0.0,-0.0};

        const auto reconciled=algo::resolveFacetConflict(
            compatible,paacs::bnd::VolumeFacetKey{40,0});
        CHECK(reconciled.resolved());
        CHECK(reconciled.chart==7);
        const auto reconciledChild=childFromResolvedFacet(reconciled);
        CHECK(reconciledChild.ready());

        algo::SurfaceSetResult primaryOnly;
        primaryOnly.faces.push_back(retainedFacetFace(40));
        const auto primary=algo::resolveFacetConflict(
            primaryOnly,paacs::bnd::VolumeFacetKey{40,0});
        CHECK(primary.resolved());
        const auto primaryChild=childFromResolvedFacet(primary);
        CHECK(primaryChild.ready());

        algo::SurfaceSetResult signedZeroOnly;
        signedZeroOnly.faces.push_back(retainedFacetFace(40));
        signedZeroOnly.faces.back()->chart=8;
        signedZeroOnly.faces.back()->chartCertificate.chart=8;
        signedZeroOnly.faces.back()->localMappingParams[1]={-0.0,-0.0};
        const auto signedZero=algo::resolveFacetConflict(
            signedZeroOnly,paacs::bnd::VolumeFacetKey{40,0});
        CHECK(signedZero.resolved());
        const auto signedZeroChild=childFromResolvedFacet(signedZero);
        CHECK(signedZeroChild.ready());

        // Compatible duplicate authority produces one semantic child identity,
        // even when one chart presents both zero components with the other
        // IEEE-754 sign bit.
        CHECK(reconciledChild.identity==primaryChild.identity);
        CHECK(reconciledChild.identity==signedZeroChild.identity);
        CHECK(reconciledChild.identity.parameterUBits==
              primaryChild.identity.parameterUBits);
        CHECK(reconciledChild.identity.parameterVBits==
              primaryChild.identity.parameterVBits);

        compatible.faces.back()->localMappingParams[1].u=0.125;
        CHECK(compatible.faces.back()->localMappingParams[1].u!=0.0);
        const auto differentNonzero=algo::resolveFacetConflict(
            compatible,paacs::bnd::VolumeFacetKey{40,0});
        CHECK(differentNonzero.status==
              algo::FacetConflictResolveStatus::Unsupported);
        CHECK(differentNonzero.reason==
              algo::FacetConflictResolveReason::AmbiguousFace);
    }

    {
        algo::SurfaceSetResult coincidentAuthorities;
        coincidentAuthorities.faces.push_back(retainedFacetFace(40));
        coincidentAuthorities.faces.push_back(retainedFacetFace(41));
        const auto face40=algo::resolveFacetConflict(
            coincidentAuthorities,paacs::bnd::VolumeFacetKey{40,0});
        const auto face41=algo::resolveFacetConflict(
            coincidentAuthorities,paacs::bnd::VolumeFacetKey{41,0});
        CHECK(face40.resolved());
        CHECK(face41.resolved());
        const auto child40=childFromResolvedFacet(face40);
        const auto child41=childFromResolvedFacet(face41);
        CHECK(child40.ready());
        CHECK(child41.ready());
        CHECK(child40.parameter.u==child41.parameter.u);
        CHECK(child40.parameter.v==child41.parameter.v);
        CHECK(child40.identity.parentFacetGlobalNodes==
              child41.identity.parentFacetGlobalNodes);
        CHECK(!(child40.identity.parentArtifact==
                child41.identity.parentArtifact));
        CHECK(!(child40.identity==child41.identity));
    }

    retained.faces[0]->localCells[1]=retained.faces[0]->localCells[0];
    result=algo::resolveFacetConflict(
        retained,paacs::bnd::VolumeFacetKey{40,0});
    CHECK(result.status==algo::FacetConflictResolveStatus::Unsupported);
    CHECK(result.reason==
          algo::FacetConflictResolveReason::DuplicateCanonicalFacet);
    CHECK(result.workUnits>1);

    retained.faces[0]=retainedFacetFace(40);
    retained.faces[0]->publishedGlobalByLocalNode.pop_back();
    result=algo::resolveFacetConflict(
        retained,paacs::bnd::VolumeFacetKey{40,0});
    CHECK(result.status==algo::FacetConflictResolveStatus::Invalid);
    CHECK(result.reason==
          algo::FacetConflictResolveReason::IncompleteFacePayload);

    retained.faces[0]=retainedFacetFace(40);
    retained.faces[0]->chartCertificate.certificateStamp={};
    result=algo::resolveFacetConflict(
        retained,paacs::bnd::VolumeFacetKey{40,0});
    CHECK(result.status==algo::FacetConflictResolveStatus::Invalid);
    CHECK(result.reason==
          algo::FacetConflictResolveReason::InvalidChartCertificate);
    retained.faces[0]=retainedFacetFace(40);
    retained.faces[0]->expectedGeometryStamp={99,1};
    result=algo::resolveFacetConflict(
        retained,paacs::bnd::VolumeFacetKey{40,0});
    CHECK(result.status==algo::FacetConflictResolveStatus::Invalid);
    CHECK(result.reason==
          algo::FacetConflictResolveReason::InvalidChartCertificate);
}

} // namespace

namespace paacs_tests {
struct Entry { const char* name;void(*fn)(); };
const Entry* t2cCacheTests(std::size_t& count) {
    static const Entry tests[]{
        {"t2c_cache.immutable_slots",immutableSlotReplacement},
        {"t2c_cache.sizing_replay",sizingPrefixReplay},
        {"t2c_cache.structural_key",structuralKeyCannotLaunderSizingOrFlatIndices},
        {"t2c_cache.invariant_authority",invariantAuthorityCannotBecomeDirty},
        {"t2c_cache.predicted_dependencies",exactPredictedDirtyDependencies},
        {"t2c_cache.global_sizing_order",globalSizingOrderIsCanonical},
        {"t2c_cache.report_xor",completeFaceXorDiagnostic},
        {"t2c_cache.exact_core_metrics",exactCoreMetricsAndDirtyLimit},
        {"t2c_cache.four_batch",fourBatchReplacementIsLocalAndDeterministic},
        {"t2c_cache.facet_proposal",facetProposalSecuresCacheBeforePlacement},
        {"t2c_cache.facet_proposal_capacity",facetProposalStopsAtEveryCapacityBeforePlacement},
        {"t2c_cache.facet_barycentre_bits",facetBarycentreCanonicalBits},
        {"t2c_cache.facet_barycentre_preflight",facetBarycentrePreflight},
        {"t2c_cache.facet_barycentre_real",facetBarycentreRealLayerGoldens},
        {"t2c_cache.facet_resolver_rank",facetConflictResolverUsesNormalizedRank},
        {"t2c_cache.facet_resolver_reject",facetConflictResolverRejectsMissingAmbiguousAndCorrupt}};
    count=std::size(tests); return tests;
}
} // namespace paacs_tests
