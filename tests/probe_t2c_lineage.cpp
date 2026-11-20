// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Testing-only structural oracle for T2-C outer-attempt edge closure.
// This executable deliberately does not qualify work or memory accounting.

#include "t2c_facet_probe_fixture.hpp"

#include "paacs/bnd/TetFacetSetRecovery.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/predicates/Predicates.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <set>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

namespace fixture=paacs::tests::t2c_facet_probe;
namespace {

using namespace paacs;
using namespace paacs::bnd;

int failures=0;
#define REQUIRE(condition) do { \
    if(!(condition)) { \
        std::printf("FAIL %d: %s\n",__LINE__,#condition); \
        ++failures; \
    } \
} while(false)

constexpr std::array<std::array<int,2>,6> kPlainEdges{{
    {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
constexpr std::array<std::array<int,3>,4> kPlainFaces{{
    {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};

bool classEqual(const mesh::Classification& a,
                const mesh::Classification& b) noexcept {
    return a.dim==b.dim&&a.tag==b.tag;
}

bool handleLess(topo::CellHandle a,topo::CellHandle b) noexcept {
    return std::tie(a.pool,a.slot,a.gen)<std::tie(b.pool,b.slot,b.gen);
}

const T2CEdgeLineageVertexForTesting* vertex(
    const T2CEdgeLineageSinkForTesting& sink,
    const T2CEdgeLineageAttemptForTesting& attempt,
    DenseTetNodeId dense) noexcept {
    if(dense>=attempt.vertexCount)return nullptr;
    const auto index=attempt.vertexOffset+static_cast<std::size_t>(dense);
    if(index>=sink.vertexRows.size())return nullptr;
    const auto& row=sink.vertexRows[index];
    return row.dense==dense?&row:nullptr;
}

bool endpointWedgePlain(
    const Vec3& endpoint,const Vec3& other,
    const std::array<Vec3,4>& tet,
    const std::array<DenseTetNodeId,4>& nodes,
    DenseTetNodeId endpointId) {
    if(std::find(nodes.begin(),nodes.end(),endpointId)==nodes.end())
        return false;
    bool constrained=false;
    for(const auto& face:kPlainFaces) {
        const auto atEndpoint=predicates::orient3d(
            tet[static_cast<std::size_t>(face[0])],
            tet[static_cast<std::size_t>(face[1])],
            tet[static_cast<std::size_t>(face[2])],endpoint);
        if(atEndpoint!=predicates::Sign::Zero)continue;
        constrained=true;
        if(predicates::orient3d(
               tet[static_cast<std::size_t>(face[0])],
               tet[static_cast<std::size_t>(face[1])],
               tet[static_cast<std::size_t>(face[2])],other)==
           predicates::Sign::Positive)return false;
    }
    return constrained;
}

template<class Cell>
bool cellTouchesTargetPlain(
    const T2CEdgeLineageSinkForTesting& sink,
    const T2CEdgeLineageAttemptForTesting& attempt,const Cell& cell,
    const std::array<DenseTetNodeId,2>& targetDense) {
    const auto* a=vertex(sink,attempt,targetDense[0]);
    const auto* b=vertex(sink,attempt,targetDense[1]);
    if(!a||!b)return true;
    std::array<Vec3,4> tet{};
    for(std::size_t i=0;i<tet.size();++i) {
        const auto* row=vertex(sink,attempt,cell.denseNodes[i]);
        if(!row)return true;
        tet[i]=row->position;
    }
    int crossings=0;
    for(const auto& face:kPlainFaces)
        crossings+=tet_facet_detail::segmentCrossesOpenTriangle(
            a->position,b->position,
            tet[static_cast<std::size_t>(face[0])],
            tet[static_cast<std::size_t>(face[1])],
            tet[static_cast<std::size_t>(face[2])])?1:0;
    if(crossings>=2)return true;
    return endpointWedgePlain(a->position,b->position,tet,
                              cell.denseNodes,targetDense[0])||
           endpointWedgePlain(b->position,a->position,tet,
                              cell.denseNodes,targetDense[1]);
}

bool rangeFits(std::size_t offset,std::size_t count,
               std::size_t size) noexcept {
    return offset<=size&&count<=size-offset;
}

bool cellTouchesBlockerPlain(
    const T2CEdgeLineageSinkForTesting& sink,
    const T2CEdgeLineageAttemptForTesting& attempt,
    const T2CEdgeLineageSupportCellForTesting& cell,
    const T2CEdgeLineageBlockerForTesting& blocker) {
    if(cellTouchesTargetPlain(sink,attempt,cell,blocker.targetDense))
        return true;
    if(!rangeFits(blocker.supportEdgeOffset,blocker.supportEdgeCount,
                  sink.edgeRows.size())||
       !rangeFits(blocker.supportFaceOffset,blocker.supportFaceCount,
                  sink.faceRows.size()))return true;
    for(const auto& edge:kPlainEdges) {
        const auto key=DenseTetEdgeKey::of(
            cell.denseNodes[static_cast<std::size_t>(edge[0])],
            cell.denseNodes[static_cast<std::size_t>(edge[1])]);
        for(std::size_t i=0;i<blocker.supportEdgeCount;++i)
            if(sink.edgeRows[blocker.supportEdgeOffset+i].key==key)
                return true;
    }
    for(const auto& face:kPlainFaces) {
        const auto key=DenseTetFaceKey::of({
            cell.denseNodes[static_cast<std::size_t>(face[0])],
            cell.denseNodes[static_cast<std::size_t>(face[1])],
            cell.denseNodes[static_cast<std::size_t>(face[2])]});
        for(std::size_t i=0;i<blocker.supportFaceCount;++i)
            if(sink.faceRows[blocker.supportFaceOffset+i].key==key)
                return true;
    }
    return false;
}

bool containsHandle(const T2CEdgeLineageSinkForTesting& sink,
                    const T2CEdgeLineageBlockerForTesting& blocker,
                    topo::CellHandle handle) {
    if(!rangeFits(blocker.supportCellOffset,blocker.supportCellCount,
                  sink.cellRows.size()))return true;
    for(std::size_t i=0;i<blocker.supportCellCount;++i)
        if(sink.cellRows[blocker.supportCellOffset+i]==handle)return true;
    return false;
}

bool blockersIndependentPlain(
    const T2CEdgeLineageSinkForTesting& sink,
    const T2CEdgeLineageAttemptForTesting& attempt,
    const T2CEdgeLineageBlockerForTesting& a,
    const T2CEdgeLineageBlockerForTesting& b) {
    if(!rangeFits(a.supportCellOffset,a.supportCellCount,sink.cellRows.size())||
       !rangeFits(b.supportCellOffset,b.supportCellCount,sink.cellRows.size())||
       !rangeFits(a.supportPayloadOffset,a.supportPayloadCount,
                  sink.supportCellRows.size())||
       !rangeFits(b.supportPayloadOffset,b.supportPayloadCount,
                  sink.supportCellRows.size())||
       !rangeFits(a.supportEdgeOffset,a.supportEdgeCount,sink.edgeRows.size())||
       !rangeFits(b.supportEdgeOffset,b.supportEdgeCount,sink.edgeRows.size())||
       !rangeFits(a.supportFaceOffset,a.supportFaceCount,sink.faceRows.size())||
       !rangeFits(b.supportFaceOffset,b.supportFaceCount,sink.faceRows.size()))
        return false;
    for(std::size_t i=0;i<a.supportCellCount;++i)
        for(std::size_t j=0;j<b.supportCellCount;++j)
            if(sink.cellRows[a.supportCellOffset+i]==
               sink.cellRows[b.supportCellOffset+j])return false;
    for(std::size_t i=0;i<a.supportEdgeCount;++i)
        for(std::size_t j=0;j<b.supportEdgeCount;++j)
            if(sink.edgeRows[a.supportEdgeOffset+i].key==
               sink.edgeRows[b.supportEdgeOffset+j].key)return false;
    for(std::size_t i=0;i<a.supportFaceCount;++i)
        for(std::size_t j=0;j<b.supportFaceCount;++j)
            if(sink.faceRows[a.supportFaceOffset+i].key==
               sink.faceRows[b.supportFaceOffset+j].key)return false;
    for(std::size_t i=0;i<a.supportPayloadCount;++i)
        if(cellTouchesBlockerPlain(
               sink,attempt,sink.supportCellRows[a.supportPayloadOffset+i],b))
            return false;
    for(std::size_t i=0;i<b.supportPayloadCount;++i)
        if(cellTouchesBlockerPlain(
               sink,attempt,sink.supportCellRows[b.supportPayloadOffset+i],a))
            return false;
    return true;
}

bool selectedTouchesTargetPlain(
    const T2CEdgeLineageSinkForTesting& sink,
    const T2CEdgeLineageAttemptForTesting& attempt,
    const T2CEdgeLineageAssessmentForTesting& source,
    const std::array<DenseTetNodeId,2>& targetDense) {
    if(!rangeFits(source.selectedPayloadOffset,source.selectedPayloadCount,
                  sink.supportCellRows.size()))return true;
    for(std::size_t i=0;i<source.selectedPayloadCount;++i)
        if(cellTouchesTargetPlain(
               sink,attempt,
               sink.supportCellRows[source.selectedPayloadOffset+i],
               targetDense))return true;
    return false;
}

void mix(std::uint64_t& digest,std::uint64_t value) noexcept {
    digest^=value+UINT64_C(0x9e3779b97f4a7c15)+(digest<<6)+(digest>>2);
}

template<class Handle>
void mixHandle(std::uint64_t& digest,Handle handle) noexcept {
    mix(digest,handle.slot);mix(digest,handle.gen);mix(digest,handle.pool);
}

struct PlainEdge {
    T2CEdgeLineageSupportEdgeForTesting row{};
    std::vector<topo::CellHandle> star;
};
struct PlainFace {
    T2CEdgeLineageSupportFaceForTesting row{};
    std::vector<T2CEdgeLineageSupportFaceOwnerForTesting> owners;
};
struct PlainCertificate {
    T2CEdgeLineageBlockerForTesting row{};
    std::array<T2CEdgeLineageVertexForTesting,2> targets{};
    std::vector<topo::CellHandle> selected,support;
    std::vector<T2CEdgeLineageSupportCellForTesting> cells;
    std::vector<PlainEdge> edges;
    std::vector<PlainFace> faces;
};

bool extractCertificate(const T2CEdgeLineageSinkForTesting& sink,
                        const T2CEdgeLineageAttemptForTesting& attempt,
                        const T2CEdgeLineageBlockerForTesting& row,
                        PlainCertificate& out) {
    out={};out.row=row;
    for(std::size_t i=0;i<2;++i) {
        const auto* found=vertex(sink,attempt,row.targetDense[i]);
        if(!found)return false;
        out.targets[i]=*found;
    }
    if(!rangeFits(row.selectedCellOffset,row.selectedCellCount,
                  sink.cellRows.size())||
       !rangeFits(row.supportCellOffset,row.supportCellCount,
                  sink.cellRows.size())||
       !rangeFits(row.supportPayloadOffset,row.supportPayloadCount,
                  sink.supportCellRows.size())||
       !rangeFits(row.supportEdgeOffset,row.supportEdgeCount,
                  sink.edgeRows.size())||
       !rangeFits(row.supportFaceOffset,row.supportFaceCount,
                  sink.faceRows.size()))return false;
    out.selected.assign(sink.cellRows.begin()+
                            static_cast<std::ptrdiff_t>(row.selectedCellOffset),
                        sink.cellRows.begin()+static_cast<std::ptrdiff_t>(
                            row.selectedCellOffset+row.selectedCellCount));
    out.support.assign(sink.cellRows.begin()+
                           static_cast<std::ptrdiff_t>(row.supportCellOffset),
                       sink.cellRows.begin()+static_cast<std::ptrdiff_t>(
                           row.supportCellOffset+row.supportCellCount));
    out.cells.assign(sink.supportCellRows.begin()+
                         static_cast<std::ptrdiff_t>(row.supportPayloadOffset),
                     sink.supportCellRows.begin()+static_cast<std::ptrdiff_t>(
                         row.supportPayloadOffset+row.supportPayloadCount));
    for(std::size_t i=0;i<row.supportEdgeCount;++i) {
        PlainEdge edge;edge.row=sink.edgeRows[row.supportEdgeOffset+i];
        if(!rangeFits(edge.row.completeStarOffset,
                      edge.row.completeStarCount,sink.cellRows.size()))
            return false;
        edge.star.assign(sink.cellRows.begin()+static_cast<std::ptrdiff_t>(
                             edge.row.completeStarOffset),
                         sink.cellRows.begin()+static_cast<std::ptrdiff_t>(
                             edge.row.completeStarOffset+
                             edge.row.completeStarCount));
        out.edges.push_back(std::move(edge));
    }
    for(std::size_t i=0;i<row.supportFaceCount;++i) {
        PlainFace face;face.row=sink.faceRows[row.supportFaceOffset+i];
        if(!rangeFits(face.row.ownerOffset,face.row.ownerCount,
                      sink.faceOwnerRows.size()))return false;
        face.owners.assign(sink.faceOwnerRows.begin()+
                               static_cast<std::ptrdiff_t>(face.row.ownerOffset),
                           sink.faceOwnerRows.begin()+static_cast<std::ptrdiff_t>(
                               face.row.ownerOffset+face.row.ownerCount));
        out.faces.push_back(std::move(face));
    }
    return true;
}

std::uint64_t digestCertificate(const PlainCertificate& value) noexcept {
    std::uint64_t digest=UINT64_C(0xcbf29ce484222325);
    mix(digest,UINT64_C(0x4c494e4541474533));
    mix(digest,value.row.artifactStamp.persistentId);
    mix(digest,value.row.artifactStamp.revision);
    mix(digest,static_cast<std::uint64_t>(
        value.row.key.canonicalGlobalNodes[0]));
    mix(digest,static_cast<std::uint64_t>(
        value.row.key.canonicalGlobalNodes[1]));
    mix(digest,static_cast<std::uint64_t>(value.row.classification.dim));
    mix(digest,static_cast<std::uint64_t>(value.row.classification.tag));
    mix(digest,static_cast<std::uint64_t>(value.row.region.dim));
    mix(digest,static_cast<std::uint64_t>(value.row.region.tag));
    for(std::size_t i=0;i<2;++i) {
        mixHandle(digest,value.targets[i].handle);
        mix(digest,value.row.targetDense[i]);
        mix(digest,std::bit_cast<std::uint64_t>(value.targets[i].position.x));
        mix(digest,std::bit_cast<std::uint64_t>(value.targets[i].position.y));
        mix(digest,std::bit_cast<std::uint64_t>(value.targets[i].position.z));
        mix(digest,static_cast<std::uint64_t>(
            value.targets[i].classification.dim));
        mix(digest,static_cast<std::uint64_t>(
            value.targets[i].classification.tag));
    }
    mix(digest,UINT64_C(0x53454c4543544544));mix(digest,value.selected.size());
    for(const auto handle:value.selected)mixHandle(digest,handle);
    mix(digest,UINT64_C(0x535550504f525443));mix(digest,value.support.size());
    for(const auto handle:value.support)mixHandle(digest,handle);
    mix(digest,UINT64_C(0x43454c4c5041594c));mix(digest,value.cells.size());
    for(const auto& cell:value.cells) {
        mixHandle(digest,cell.handle);mix(digest,cell.sequence);
        for(const auto node:cell.denseNodes)mix(digest,node);
        mix(digest,static_cast<std::uint64_t>(cell.classification.dim));
        mix(digest,static_cast<std::uint64_t>(cell.classification.tag));
    }
    mix(digest,UINT64_C(0x45444745524f5753));mix(digest,value.edges.size());
    for(const auto& edge:value.edges) {
        mix(digest,edge.row.key.nodes[0]);mix(digest,edge.row.key.nodes[1]);
        mixHandle(digest,edge.row.handle);mix(digest,edge.row.required?1U:0U);
        mix(digest,static_cast<std::uint64_t>(edge.row.classification.dim));
        mix(digest,static_cast<std::uint64_t>(edge.row.classification.tag));
        mix(digest,edge.star.size());
        for(const auto handle:edge.star)mixHandle(digest,handle);
    }
    mix(digest,UINT64_C(0x46414345524f5753));mix(digest,value.faces.size());
    for(const auto& face:value.faces) {
        for(const auto node:face.row.key.nodes)mix(digest,node);
        mixHandle(digest,face.row.handle);
        mix(digest,static_cast<std::uint64_t>(face.row.classification.dim));
        mix(digest,static_cast<std::uint64_t>(face.row.classification.tag));
        mix(digest,face.owners.size());
        for(const auto& owner:face.owners) {
            mixHandle(digest,owner.cell);
            mix(digest,static_cast<std::uint64_t>(
                static_cast<std::int64_t>(owner.sense)));
        }
    }
    return digest;
}

const T2CEdgeLineageAttemptForTesting* attemptForBlocker(
    const T2CEdgeLineageSinkForTesting& sink,std::size_t blocker) noexcept {
    for(const auto& attempt:sink.attempts)
        if(blocker>=attempt.blockerOffset&&
           blocker<attempt.blockerOffset+attempt.blockerCount)return &attempt;
    return nullptr;
}

const T2CEdgeLineageAttemptForTesting* attemptForRecovery(
    const T2CEdgeLineageSinkForTesting& sink,std::size_t recovery) noexcept {
    for(const auto& attempt:sink.attempts)
        if(recovery>=attempt.recoveryOffset&&
           recovery<attempt.recoveryOffset+attempt.recoveryCount)return &attempt;
    return nullptr;
}

void proveEndpointWedgeDiscriminator() {
    const std::array<Vec3,4> tet{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}};
    const std::array<DenseTetNodeId,4> nodes{{0,1,2,3}};
    const Vec3 endpoint{0,0,0},other{0.125,0.125,0.125};
    int crossings=0;
    for(const auto& face:kPlainFaces)
        crossings+=tet_facet_detail::segmentCrossesOpenTriangle(
            endpoint,other,tet[static_cast<std::size_t>(face[0])],
            tet[static_cast<std::size_t>(face[1])],
            tet[static_cast<std::size_t>(face[2])])?1:0;
    REQUIRE(crossings<2);
    REQUIRE(endpointWedgePlain(endpoint,other,tet,nodes,0));
    REQUIRE(!endpointWedgePlain(other,endpoint,tet,nodes,4));
}

void proveBooleanExhaustion() {
    for(unsigned mask=0;mask<32;++mask) {
        const bool cell=(mask&1U)!=0,edge=(mask&2U)!=0;
        const bool face=(mask&4U)!=0,ab=(mask&8U)!=0,ba=(mask&16U)!=0;
        const bool independent=!(cell||edge||face||ab||ba);
        REQUIRE(independent==(mask==0));
    }
    for(unsigned mask=0;mask<32;++mask) {
        const bool retired=(mask&1U)!=0,candidate=(mask&2U)!=0;
        const bool alias=(mask&4U)!=0,stale=(mask&8U)!=0;
        const bool preparation=(mask&16U)!=0;
        const bool allow=!retired&&!candidate&&!alias&&!stale&&preparation;
        REQUIRE(allow==(mask==16U));
    }
}

void replayCertificatesAndCorruptions(
    const T2CEdgeLineageSinkForTesting& sink) {
    bool discriminated=false;
    for(std::size_t i=0;i<sink.blockers.size();++i) {
        const auto* attempt=attemptForBlocker(sink,i);
        REQUIRE(attempt!=nullptr);if(!attempt)continue;
        PlainCertificate certificate;
        REQUIRE(extractCertificate(sink,*attempt,sink.blockers[i],certificate));
        if(!extractCertificate(sink,*attempt,sink.blockers[i],certificate))
            continue;
        REQUIRE(digestCertificate(certificate)==sink.blockers[i].supportDigest);
        if(discriminated||certificate.edges.empty()||certificate.faces.empty()||
           certificate.cells.empty())continue;
        const auto baseline=digestCertificate(certificate);
        ++certificate.row.artifactStamp.revision;
        REQUIRE(digestCertificate(certificate)!=baseline);
        --certificate.row.artifactStamp.revision;
        ++certificate.row.region.tag;
        REQUIRE(digestCertificate(certificate)!=baseline);
        --certificate.row.region.tag;
        ++certificate.targets[0].classification.tag;
        REQUIRE(digestCertificate(certificate)!=baseline);
        --certificate.targets[0].classification.tag;
        ++certificate.edges[0].row.classification.tag;
        REQUIRE(digestCertificate(certificate)!=baseline);
        --certificate.edges[0].row.classification.tag;
        ++certificate.faces[0].row.classification.tag;
        REQUIRE(digestCertificate(certificate)!=baseline);
        --certificate.faces[0].row.classification.tag;
        ++certificate.cells[0].sequence;
        REQUIRE(digestCertificate(certificate)!=baseline);
        --certificate.cells[0].sequence;
        discriminated=true;
    }
    if(!sink.blockers.empty())REQUIRE(discriminated);
}

void replayGraph(const T2CEdgeLineageSinkForTesting& sink) {
    std::size_t blockerDecisions=0,candidateDecisions=0,preparedDecisions=0;
    std::vector<unsigned> passUse(sink.passRows.size(),0);
    std::vector<unsigned> assessmentUse(sink.assessmentRows.size(),0);
    std::vector<unsigned> decisionUse(sink.decisionRows.size(),0);
    std::vector<unsigned> blockerIndexUse(sink.blockerIndexRows.size(),0);
    for(const auto& attempt:sink.attempts) {
        REQUIRE(rangeFits(attempt.passOffset,attempt.passCount,
                          sink.passRows.size()));
        for(std::size_t pi=0;pi<attempt.passCount;++pi) {
            const auto passIndex=attempt.passOffset+pi;
            ++passUse[passIndex];
            const auto& pass=sink.passRows[passIndex];
            REQUIRE(pass.canonicalScanComplete);
            REQUIRE(rangeFits(pass.assessmentOffset,pass.assessmentCount,
                              sink.assessmentRows.size()));
            REQUIRE(rangeFits(pass.decisionOffset,pass.decisionCount,
                              sink.decisionRows.size()));
            REQUIRE(rangeFits(pass.frozenBlockerOffset,
                              pass.frozenBlockerCount,
                              sink.blockerIndexRows.size()));
            for(std::size_t i=0;i<pass.assessmentCount;++i)
                ++assessmentUse[pass.assessmentOffset+i];
            for(std::size_t i=0;i<pass.decisionCount;++i)
                ++decisionUse[pass.decisionOffset+i];
            for(std::size_t i=0;i<pass.frozenBlockerCount;++i)
                ++blockerIndexUse[pass.frozenBlockerOffset+i];
            for(std::size_t ai=0;ai<pass.assessmentCount;++ai) {
                const auto& row=sink.assessmentRows[pass.assessmentOffset+ai];
                REQUIRE(row.canonicalOrdinal==ai);
                if(ai!=0)REQUIRE(sink.assessmentRows[
                    pass.assessmentOffset+ai-1].key<row.key);
                REQUIRE(row.candidate==
                        (row.outcome.status==TetEdgeRecoveryStatus::Recovered));
            }
            std::size_t expectedPermanent=(std::numeric_limits<std::size_t>::max)();
            for(std::size_t ai=0;ai<pass.assessmentCount;++ai) {
                const auto& row=sink.assessmentRows[pass.assessmentOffset+ai];
                if(row.excludedFrozen)continue;
                if(!row.mapped||
                   (row.outcome.status!=TetEdgeRecoveryStatus::AlreadyPresent&&
                    row.outcome.status!=TetEdgeRecoveryStatus::Recovered)) {
                    expectedPermanent=pass.assessmentOffset+ai;break;
                }
            }
            REQUIRE(pass.permanentAssessment==expectedPermanent);
            std::vector<T2CEdgeLineageDependencyForTesting> expected;
            std::set<std::size_t> blocked;
            const auto appendExpected=[&](auto kind,std::size_t source,
                                          std::size_t target) {
                T2CEdgeLineageDependencyForTesting row;
                row.kind=kind;row.sourceRow=source;row.targetRow=target;
                expected.push_back(row);
            };
            if(expectedPermanent!=(std::numeric_limits<std::size_t>::max)()) {
                const auto& boundary=sink.assessmentRows[expectedPermanent];
                const bool blocker=boundary.outcome.status==
                        TetEdgeRecoveryStatus::Refused&&
                    boundary.outcome.reason==TetEdgeRecoveryReason::NonCycleLink&&
                    boundary.outcome.phase==TetEdgeRecoveryPhase::Discovery;
                if(blocker) {
                    REQUIRE(pass.discoveredBlocker<sink.blockers.size());
                    bool discoveredTouches=false;
                    if(pass.discoveredBlocker<sink.blockers.size()) {
                        const auto& discovered=
                            sink.blockers[pass.discoveredBlocker];
                        REQUIRE(discovered.key==boundary.key);
                        REQUIRE(discovered.observedStamp==pass.stamp);
                        REQUIRE(discovered.artifactStamp==
                                boundary.artifactStamp);
                        REQUIRE(classEqual(discovered.classification,
                                           boundary.classification));
                        REQUIRE(classEqual(discovered.region,boundary.region));
                        REQUIRE(discovered.outcome.status==
                                boundary.outcome.status);
                        REQUIRE(discovered.outcome.reason==
                                boundary.outcome.reason);
                        REQUIRE(discovered.outcome.phase==
                                boundary.outcome.phase);
                    }
                    for(std::size_t i=0;i<pass.frozenBlockerCount;++i) {
                        const auto frozenIndex=
                            sink.blockerIndexRows[pass.frozenBlockerOffset+i];
                        appendExpected(
                            T2CEdgeLineageDependencyKindForTesting::
                                BlockerToBlocker,
                            frozenIndex,
                            pass.discoveredBlocker);
                        if(frozenIndex<sink.blockers.size()&&
                           pass.discoveredBlocker<sink.blockers.size())
                            discoveredTouches=discoveredTouches||
                                !blockersIndependentPlain(
                                    sink,attempt,sink.blockers[frozenIndex],
                                    sink.blockers[pass.discoveredBlocker]);
                    }
                    if(pass.discoveredBlocker<sink.blockers.size()) {
                        const auto& discovered=
                            sink.blockers[pass.discoveredBlocker];
                        REQUIRE(discovered.deferred==discoveredTouches);
                        REQUIRE(discovered.pairwiseIndependent==
                                !discoveredTouches);
                        if(discoveredTouches)REQUIRE(!attempt.frontierReplaced);
                    }
                } else REQUIRE(pass.discoveredBlocker==
                    (std::numeric_limits<std::size_t>::max)());
            } else {
                std::vector<bool> eligible(pass.assessmentCount,false);
                for(std::size_t ai=0;ai<pass.assessmentCount;++ai) {
                    const auto sourceIndex=pass.assessmentOffset+ai;
                    const auto& source=sink.assessmentRows[sourceIndex];
                    if(!source.candidate||source.deferredAtEntry)continue;
                    bool frozenTouch=false;
                    for(std::size_t bi=0;bi<pass.frozenBlockerCount;++bi) {
                        const auto blockerIndex=
                            sink.blockerIndexRows[pass.frozenBlockerOffset+bi];
                        REQUIRE(blockerIndex<sink.blockers.size());
                        appendExpected(
                            T2CEdgeLineageDependencyKindForTesting::
                                CandidateToFrozenBlocker,
                            sourceIndex,blockerIndex);
                        if(blockerIndex<sink.blockers.size())
                            frozenTouch=frozenTouch||selectedTouchesTargetPlain(
                                sink,attempt,source,
                                sink.blockers[blockerIndex].targetDense);
                    }
                    REQUIRE(source.deferred==
                            (source.deferredAtEntry||frozenTouch||
                             source.preparedDeferred));
                    eligible[ai]=!frozenTouch;
                }
                for(std::size_t ai=0;ai<pass.assessmentCount;++ai) {
                    if(!eligible[ai])continue;
                    for(std::size_t aj=0;aj<pass.assessmentCount;++aj) {
                        if(ai==aj||!eligible[aj])continue;
                        const auto source=pass.assessmentOffset+ai;
                        const auto target=pass.assessmentOffset+aj;
                        appendExpected(
                            T2CEdgeLineageDependencyKindForTesting::
                                CandidateToCandidate,source,target);
                        if(selectedTouchesTargetPlain(
                               sink,attempt,sink.assessmentRows[source],
                               sink.assessmentRows[target].targetDense))
                            blocked.insert(source);
                    }
                }
                std::size_t selected=(std::numeric_limits<std::size_t>::max)();
                for(std::size_t ai=0;ai<pass.assessmentCount;++ai)
                    if(eligible[ai]&&
                       blocked.count(pass.assessmentOffset+ai)==0) {
                        selected=pass.assessmentOffset+ai;break;
                    }
                REQUIRE(pass.selectedAssessment==selected);
                REQUIRE(pass.edgeFixedPoint==
                        (selected==(std::numeric_limits<std::size_t>::max)()));
                if(selected!=(std::numeric_limits<std::size_t>::max)()) {
                    REQUIRE(pass.preparedRecovery<sink.recoveries.size());
                    for(std::size_t bi=0;bi<pass.frozenBlockerCount;++bi)
                        appendExpected(
                            T2CEdgeLineageDependencyKindForTesting::
                                PreparedDeltaToFrozenBlocker,
                            pass.preparedRecovery,
                            sink.blockerIndexRows[pass.frozenBlockerOffset+bi]);
                } else REQUIRE(pass.preparedRecovery==
                    (std::numeric_limits<std::size_t>::max)());
            }
            REQUIRE(pass.decisionCount==expected.size());
            for(std::size_t i=0;i<(std::min)(pass.decisionCount,
                                             expected.size());++i) {
                const auto& actual=sink.decisionRows[pass.decisionOffset+i];
                REQUIRE(actual.kind==expected[i].kind);
                REQUIRE(actual.sourceRow==expected[i].sourceRow);
                REQUIRE(actual.targetRow==expected[i].targetRow);
            }
            for(std::size_t di=0;di<pass.decisionCount;++di) {
                const auto& decision=
                    sink.decisionRows[pass.decisionOffset+di];
                if(decision.kind==
                   T2CEdgeLineageDependencyKindForTesting::BlockerToBlocker) {
                    ++blockerDecisions;
                    REQUIRE(decision.sourceRow<sink.blockers.size());
                    REQUIRE(decision.targetRow<sink.blockers.size());
                    if(decision.sourceRow>=sink.blockers.size()||
                       decision.targetRow>=sink.blockers.size())continue;
                    const bool independent=blockersIndependentPlain(
                        sink,attempt,sink.blockers[decision.sourceRow],
                        sink.blockers[decision.targetRow]);
                    REQUIRE(decision.touches==!independent);
                } else if(decision.kind==
                   T2CEdgeLineageDependencyKindForTesting::
                       CandidateToFrozenBlocker) {
                    ++candidateDecisions;
                    REQUIRE(decision.sourceRow<sink.assessmentRows.size());
                    REQUIRE(decision.targetRow<sink.blockers.size());
                    if(decision.sourceRow>=sink.assessmentRows.size()||
                       decision.targetRow>=sink.blockers.size())continue;
                    const bool touches=selectedTouchesTargetPlain(
                        sink,attempt,sink.assessmentRows[decision.sourceRow],
                        sink.blockers[decision.targetRow].targetDense);
                    REQUIRE(touches==decision.touches);
                    if(touches)blocked.insert(decision.sourceRow);
                } else if(decision.kind==
                   T2CEdgeLineageDependencyKindForTesting::
                       CandidateToCandidate) {
                    ++candidateDecisions;
                    REQUIRE(decision.sourceRow<sink.assessmentRows.size());
                    REQUIRE(decision.targetRow<sink.assessmentRows.size());
                    if(decision.sourceRow>=sink.assessmentRows.size()||
                       decision.targetRow>=sink.assessmentRows.size())continue;
                    const bool touches=selectedTouchesTargetPlain(
                        sink,attempt,sink.assessmentRows[decision.sourceRow],
                        sink.assessmentRows[decision.targetRow].targetDense);
                    REQUIRE(touches==decision.touches);
                    if(touches)blocked.insert(decision.sourceRow);
                } else {
                    ++preparedDecisions;
                    REQUIRE(decision.sourceRow<sink.recoveries.size());
                    REQUIRE(decision.targetRow<sink.blockers.size());
                    if(decision.sourceRow>=sink.recoveries.size()||
                       decision.targetRow>=sink.blockers.size())continue;
                    const auto& recovery=sink.recoveries[decision.sourceRow];
                    const auto& blocker=sink.blockers[decision.targetRow];
                    bool retired=false,candidate=false;
                    REQUIRE(rangeFits(recovery.retiredPayloadOffset,
                                      recovery.retiredPayloadCount,
                                      sink.supportCellRows.size()));
                    REQUIRE(rangeFits(recovery.candidateOffset,
                                      recovery.candidateCount,
                                      sink.candidateRows.size()));
                    for(std::size_t i=0;i<recovery.retiredPayloadCount;++i) {
                        const auto& cell=sink.supportCellRows[
                            recovery.retiredPayloadOffset+i];
                        retired=retired||containsHandle(sink,blocker,cell.handle)||
                            cellTouchesBlockerPlain(sink,attempt,cell,blocker);
                    }
                    for(std::size_t i=0;i<recovery.candidateCount;++i) {
                        const auto& source=
                            sink.candidateRows[recovery.candidateOffset+i];
                        T2CEdgeLineageSupportCellForTesting cell;
                        cell.denseNodes=source.denseNodes;
                        cell.classification=source.classification;
                        candidate=candidate||
                            cellTouchesBlockerPlain(sink,attempt,cell,blocker);
                    }
                    REQUIRE(retired==decision.retiredTouches);
                    REQUIRE(candidate==decision.candidateTouches);
                    REQUIRE((retired||candidate)==decision.touches);
                }
            }
        }
    }
    REQUIRE(!sink.blockers.empty());
    REQUIRE(preparedDecisions!=0||sink.recoveries.empty());
    REQUIRE(blockerDecisions==sink.blockerPairTests);
    REQUIRE(candidateDecisions==sink.candidatePairTests);
    REQUIRE(preparedDecisions==sink.transactionPairTests);
    REQUIRE(blockerDecisions==22);
    REQUIRE(candidateDecisions==330);
    REQUIRE(preparedDecisions==18);
    REQUIRE(std::all_of(passUse.begin(),passUse.end(),
                        [](unsigned used) { return used==1; }));
    REQUIRE(std::all_of(assessmentUse.begin(),assessmentUse.end(),
                        [](unsigned used) { return used==1; }));
    REQUIRE(std::all_of(decisionUse.begin(),decisionUse.end(),
                        [](unsigned used) { return used==1; }));
    REQUIRE(std::all_of(blockerIndexUse.begin(),blockerIndexUse.end(),
                        [](unsigned used) { return used==1; }));
    std::printf("oracle decisions blocker=%zu candidate=%zu prepared=%zu\n",
                blockerDecisions,candidateDecisions,preparedDecisions);
}

void replayPreparedAliases(const T2CEdgeLineageSinkForTesting& sink) {
    for(const auto& recovery:sink.recoveries) {
        REQUIRE(rangeFits(recovery.boundEdgeOffset,recovery.boundEdgeCount,
                          sink.boundEdgeRows.size()));
        bool alias=false;
        for(std::size_t i=0;i<recovery.boundEdgeCount;++i) {
            const auto& row=sink.boundEdgeRows[recovery.boundEdgeOffset+i];
            if(row.handle==recovery.predictedEdge&&!(row.key==recovery.key))
                alias=true;
        }
        REQUIRE(alias==recovery.predictedBindingAlias);
        bool deltaTouches=false;
        for(const auto& decision:sink.decisionRows)
            if(decision.kind==T2CEdgeLineageDependencyKindForTesting::
                    PreparedDeltaToFrozenBlocker&&
               decision.sourceRow==static_cast<std::size_t>(
                   &recovery-sink.recoveries.data()))
                deltaTouches=deltaTouches||decision.touches;
        REQUIRE(recovery.allowed==(!alias&&!deltaTouches));
    }
}

void replayRevalidations(const T2CEdgeLineageSinkForTesting& sink) {
    bool anyDeferred=false;
    for(const auto& attempt:sink.attempts) {
        std::vector<std::size_t> active;
        std::vector<std::uint64_t> expected(attempt.blockerCount,0);
        std::vector<topo::WorkspaceStamp> finalStamp(attempt.blockerCount);
        for(std::size_t i=0;i<attempt.blockerCount;++i)
            finalStamp[i]=sink.blockers[attempt.blockerOffset+i].observedStamp;
        for(std::size_t pi=0;pi<attempt.passCount;++pi) {
            const auto& pass=sink.passRows[attempt.passOffset+pi];
            if(pass.discoveredBlocker<sink.blockers.size()) {
                const auto& blocker=sink.blockers[pass.discoveredBlocker];
                if(blocker.deferred) {
                    anyDeferred=true;
                    REQUIRE(!blocker.pairwiseIndependent);
                } else {
                    REQUIRE(blocker.pairwiseIndependent);
                    active.push_back(pass.discoveredBlocker);
                }
            }
            if(pass.preparedRecovery<sink.recoveries.size()) {
                const auto& recovery=sink.recoveries[pass.preparedRecovery];
                if(recovery.allowed&&recovery.outcome.succeeded())
                    for(const auto blockerIndex:active) {
                        REQUIRE(blockerIndex>=attempt.blockerOffset);
                        const auto local=blockerIndex-attempt.blockerOffset;
                        REQUIRE(local<expected.size());
                        if(local<expected.size()) {
                            ++expected[local];
                            finalStamp[local]=recovery.afterStamp;
                        }
                    }
            }
        }
        for(std::size_t i=0;i<attempt.blockerCount;++i) {
            const auto& blocker=sink.blockers[attempt.blockerOffset+i];
            REQUIRE(blocker.revalidations==expected[i]);
            REQUIRE(blocker.finalValidatedStamp==finalStamp[i]);
        }
    }
    REQUIRE(sink.allBlockersPairwiseIndependent==!anyDeferred);
}

bool sameBits(double a,double b) noexcept {
    return std::bit_cast<std::uint64_t>(a)==std::bit_cast<std::uint64_t>(b);
}

bool ledgerRowEqual(const algo::ProtectedFeatureChildIdentity& a,
                    const algo::ProtectedFeatureChildIdentity& b) noexcept {
    return a.parentEntityPersistentId==b.parentEntityPersistentId&&
        a.childIdentity==b.childIdentity&&a.splitOrdinal==b.splitOrdinal&&
        a.edgeTag==b.edgeTag&&a.faceTag==b.faceTag&&
        sameBits(a.parentParameterLo,b.parentParameterLo)&&
        sameBits(a.parentParameterHi,b.parentParameterHi)&&
        sameBits(a.midpointParameter,b.midpointParameter)&&
        sameBits(a.faceParameter.u,b.faceParameter.u)&&
        sameBits(a.faceParameter.v,b.faceParameter.v)&&
        sameBits(a.coordinate.x,b.coordinate.x)&&
        sameBits(a.coordinate.y,b.coordinate.y)&&
        sameBits(a.coordinate.z,b.coordinate.z)&&
        classEqual(a.classification,b.classification)&&
        a.providerStamp==b.providerStamp;
}

void requirePublishedArtifact(
    const algo::BoundaryConformingV4Outcome& outcome,const mesh::Mesh& mesh,
    const fixture::Placement& placement,
    const fixture::LayeredModel& model) {
    REQUIRE(outcome.published());
    REQUIRE(outcome.status==algo::BoundaryConformingV4Status::Published);
    REQUIRE(outcome.reason==algo::BoundaryConformingV4Reason::None);
    REQUIRE(outcome.phase==algo::BoundaryConformingV4Phase::Publication);
    REQUIRE(!outcome.destinationUnchanged);
    REQUIRE(mesh.nodeCount()==198);
    REQUIRE(mesh.countCells(mesh::CellType::Tet4)==777);
    REQUIRE(mesh::meshHash(mesh)==UINT64_C(2737597274546495434));
    const auto semantic=mesh::semanticHashEx(mesh);
    REQUIRE(semantic.canonical);
    REQUIRE(semantic.value==UINT64_C(16753220727311122020));
    REQUIRE(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    const auto& retained=*outcome.retainedResult;
    REQUIRE(retained.subdivisionLedger.size()==65);
    REQUIRE(retained.surfaceCacheReport.reportComplete);
    REQUIRE(!retained.surfaceCacheReport.hasFailedFace);
    REQUIRE(retained.surfaceCacheReport.unattempted.empty());
    REQUIRE(retained.surfaceCacheReport.processed.size()==model.faces().size());
    std::set<algo::FacePersistentId> cacheFaces;
    for(const auto& row:retained.surfaceCacheReport.processed)
        REQUIRE(cacheFaces.insert(row.facePersistentId).second);
    for(const auto face:retained.surfaceCacheReport.unattempted)
        REQUIRE(cacheFaces.insert(face).second);
    REQUIRE(cacheFaces.size()==model.faces().size());
    std::set<algo::FacePersistentId> expectedFaces;
    for(const auto* face:model.faces())
        expectedFaces.insert(static_cast<algo::FacePersistentId>(
            face->tag())+1U);
    REQUIRE(cacheFaces==expectedFaces);
    REQUIRE(retained.metrics.edgeSplits==63);
    REQUIRE(retained.metrics.facetSplits==2);
    REQUIRE(retained.metrics.uniquePlacements==65);
    REQUIRE(retained.metrics.facetConflictRounds==1);
    REQUIRE(retained.metrics.facetConflictSccMembers==2);
    REQUIRE(retained.metrics.facetConflictDependencyRows==2);
    REQUIRE(retained.metrics.facetConflictProviderCalls==7);
    REQUIRE(placement.callsUsed==196);
    REQUIRE(retained.sourceOutcome!=nullptr);
    if(!retained.sourceOutcome)return;
    REQUIRE(retained.sourceOutcome->reportComplete);
    REQUIRE(retained.sourceOutcome->report!=nullptr);
    if(retained.sourceOutcome->report)
        REQUIRE(retained.sourceOutcome->report->reportComplete);
    REQUIRE(retained.sourceOutcome->retainedResult!=nullptr);
    if(!retained.sourceOutcome->retainedResult)return;
    const auto& audit=retained.sourceOutcome->retainedResult->boundaryAudit;
    REQUIRE(audit.passed);
    REQUIRE(audit.orphanExterior==0);
    REQUIRE(audit.unexpectedExterior==0);
    REQUIRE(audit.missingInterface==0);
    REQUIRE(audit.unexpectedInterface==0);
}

void requirePlacementLedgerOrder(
    const fixture::Placement& placement,
    std::span<const algo::ProtectedFeatureChildIdentity> ledger) {
    std::size_t placementOrdinal=0;
    for(const auto& call:placement.trace()) {
        if(call.kind!=fixture::ProviderCallKind::Placement)continue;
        REQUIRE(placementOrdinal<ledger.size());
        if(placementOrdinal>=ledger.size())continue;
        REQUIRE(call.identity==ledger[placementOrdinal].childIdentity);
        if(placementOrdinal<63)
            REQUIRE(call.entityKind!=
                    algo::ProtectedFeatureEntityKind::FacetInterior);
        else
            REQUIRE(call.entityKind==
                    algo::ProtectedFeatureEntityKind::FacetInterior);
        ++placementOrdinal;
    }
    REQUIRE(placementOrdinal==65);
    REQUIRE(placementOrdinal==ledger.size());
}

} // namespace

int main() {
    proveEndpointWedgeDiscriminator();
    proveBooleanExhaustion();

    fixture::LayeredModel model(8);
    fixture::Ownership ownership;
    fixture::bindOwnership(model,ownership);
    fixture::Placement placement;
    sizing::MeshSpec spec;
    spec.field=std::make_shared<fixture::LayerField>();spec.threads=1;
    fixture::Placement referencePlacement;
    auto referenceRequest=fixture::requestFor(ownership,referencePlacement);
    referenceRequest.recoverySchedulingMode=
        algo::TetRecoverySchedulingMode::MixedExactPresentBatching;
    mesh::Mesh referenceMesh;
    (void)referenceMesh.addNode({-7,-11,-13},{0,91});
    const auto reference=algo::meshModelBoundaryConformingV4(
        model,spec,referenceRequest,referenceMesh,ownership,
        referencePlacement);
    requirePublishedArtifact(reference,referenceMesh,referencePlacement,
                             model);
    REQUIRE(reference.retainedResult!=nullptr);
    if(reference.retainedResult)
        REQUIRE(reference.retainedResult->metrics.rounds==11);

    auto request=fixture::requestFor(ownership,placement);
    request.recoverySchedulingMode=
        algo::TetRecoverySchedulingMode::MixedExactPresentBatching;
    mesh::Mesh destination;
    (void)destination.addNode({-7,-11,-13},{0,91});
    const auto before=mesh::meshHash(destination);

    T2CEdgeLineageSinkForTesting sink;
    sink.limits.maxAttempts=32;
    sink.limits.maxPasses=4096;
    sink.limits.maxBlockers=512;
    sink.limits.maxScratchRecoveries=8192;
    sink.limits.maxCellRows=600000;
    sink.limits.maxSupportCellRows=600000;
    sink.limits.maxVertexRows=65536;
    sink.limits.maxEdgeRows=200000;
    sink.limits.maxFaceRows=200000;
    sink.limits.maxFaceOwnerRows=400000;
    sink.limits.maxCandidateRows=200000;
    sink.limits.maxBoundEdgeRows=400000;
    sink.limits.maxPassRows=8192;
    sink.limits.maxAssessmentRows=600000;
    sink.limits.maxDecisionRows=2000000;
    sink.limits.maxBlockerIndexRows=200000;
    sink.limits.maxAccountedBytes=std::size_t{512}<<20;
    sink.reset();
    setT2CEdgeLineageSinkForTesting(&sink);
    const auto outcome=algo::meshModelBoundaryConformingV4(
        model,spec,request,destination,ownership,placement);
    setT2CEdgeLineageSinkForTesting(nullptr);

    REQUIRE(sink.complete);
    REQUIRE(sink.allPreparedTransactionsGuarded);
    REQUIRE(!sink.workAccountingComplete);
    REQUIRE(!sink.memoryAccountingComplete);
    REQUIRE(!sink.attempts.empty());
    replayCertificatesAndCorruptions(sink);
    replayGraph(sink);
    replayPreparedAliases(sink);
    replayRevalidations(sink);

    std::size_t fixedPoints=0,replaced=0;
    std::size_t deferredBlockers=0,deferredRecoveries=0;
    for(std::size_t i=0;i<sink.attempts.size();++i) {
        const auto& attempt=sink.attempts[i];
        if(attempt.edgeFixedPoint)++fixedPoints;
        deferredBlockers+=attempt.deferredBlockers;
        deferredRecoveries+=attempt.deferredRecoveries;
        if(attempt.frontierReplaced) {
            ++replaced;
            REQUIRE(attempt.edgeFixedPoint);
            REQUIRE(attempt.status==T2CEdgeLineageStatusForTesting::
                    EdgeFixedPoint);
        }
        std::printf("attempt[%zu] status=%u passes=%zu blockers=%zu "
                    "recoveries=%zu deferred_b=%zu deferred_r=%zu "
                    "fixed=%u frontier=%u\n",i,
                    static_cast<unsigned>(attempt.status),attempt.passes,
                    attempt.blockerCount,attempt.recoveryCount,
                    attempt.deferredBlockers,attempt.deferredRecoveries,
                    attempt.edgeFixedPoint?1U:0U,
                    attempt.frontierReplaced?1U:0U);
    }
    REQUIRE(fixedPoints==replaced);
    REQUIRE(sink.attempts.size()==9);
    REQUIRE(sink.recoveries.size()==9);
    REQUIRE(sink.blockers.size()==24);
    REQUIRE(fixedPoints==1);
    REQUIRE(replaced==1);
    REQUIRE(deferredBlockers==8);
    REQUIRE(deferredRecoveries==2);
    for(std::size_t i=0;i<sink.attempts.size();++i) {
        const auto expected=i==7
            ? T2CEdgeLineageStatusForTesting::EdgeFixedPoint
            : T2CEdgeLineageStatusForTesting::BlockerOverlap;
        REQUIRE(sink.attempts[i].status==expected);
    }
    REQUIRE(!sink.allBlockersPairwiseIndependent);
    REQUIRE(sink.status==T2CEdgeLineageStatusForTesting::BlockerOverlap);
    auto firstOverlap=sink.decisionRows.end();
    std::size_t firstOverlapAttempt=(std::numeric_limits<std::size_t>::max)();
    std::size_t firstOverlapPass=(std::numeric_limits<std::size_t>::max)();
    for(std::size_t ai=0;ai<sink.attempts.size()&&
                          firstOverlap==sink.decisionRows.end();++ai) {
        const auto& attempt=sink.attempts[ai];
        for(std::size_t pi=0;pi<attempt.passCount&&
                              firstOverlap==sink.decisionRows.end();++pi) {
            const auto& pass=sink.passRows[attempt.passOffset+pi];
            for(std::size_t di=0;di<pass.decisionCount;++di) {
                auto it=sink.decisionRows.begin()+static_cast<std::ptrdiff_t>(
                    pass.decisionOffset+di);
                if(it->kind==T2CEdgeLineageDependencyKindForTesting::
                                 BlockerToBlocker&&it->touches) {
                    firstOverlap=it;firstOverlapAttempt=ai;
                    firstOverlapPass=pi;break;
                }
            }
        }
    }
    REQUIRE(firstOverlap!=sink.decisionRows.end());
    REQUIRE(firstOverlapAttempt==0);
    REQUIRE(firstOverlapPass==1);
    REQUIRE(!sink.attempts[0].frontierReplaced);
    if(firstOverlap!=sink.decisionRows.end()&&
       firstOverlap->sourceRow<sink.blockers.size()&&
       firstOverlap->targetRow<sink.blockers.size()) {
        const auto& source=sink.blockers[firstOverlap->sourceRow];
        const auto& target=sink.blockers[firstOverlap->targetRow];
        REQUIRE(firstOverlap->sourceRow==0);
        REQUIRE(firstOverlap->targetRow==1);
        REQUIRE(source.key.canonicalGlobalNodes==
                (std::array<Index,2>{{40,127}}));
        REQUIRE(target.key.canonicalGlobalNodes==
                (std::array<Index,2>{{42,50}}));
        REQUIRE(source.supportDigest==UINT64_C(1867039483036265797));
        REQUIRE(target.supportDigest==UINT64_C(5786267217461886265));
        REQUIRE(source.pairwiseIndependent);
        REQUIRE(!source.deferred);
        REQUIRE(!target.pairwiseIndependent);
        REQUIRE(target.deferred);
        std::printf("first_overlap attempt=%zu pass=%zu source=%zu "
                    "key=(%llu,%llu) digest=%llu "
                    "target=%zu key=(%llu,%llu) digest=%llu\n",
                    firstOverlapAttempt,firstOverlapPass,
                    firstOverlap->sourceRow,
                    static_cast<unsigned long long>(
                        source.key.canonicalGlobalNodes[0]),
                    static_cast<unsigned long long>(
                        source.key.canonicalGlobalNodes[1]),
                    static_cast<unsigned long long>(source.supportDigest),
                    firstOverlap->targetRow,
                    static_cast<unsigned long long>(
                        target.key.canonicalGlobalNodes[0]),
                    static_cast<unsigned long long>(
                        target.key.canonicalGlobalNodes[1]),
                    static_cast<unsigned long long>(target.supportDigest));
    }
    requirePublishedArtifact(outcome,destination,placement,model);
    REQUIRE(mesh::meshHash(destination)!=before);
    REQUIRE(referencePlacement.trace().size()==placement.trace().size());
    REQUIRE(std::equal(referencePlacement.trace().begin(),
                       referencePlacement.trace().end(),
                       placement.trace().begin(),placement.trace().end()));
    if(reference.retainedResult&&outcome.retainedResult) {
        const auto& expected=reference.retainedResult->subdivisionLedger;
        const auto& actual=outcome.retainedResult->subdivisionLedger;
        REQUIRE(expected.size()==actual.size());
        for(std::size_t i=0;i<(std::min)(expected.size(),actual.size());++i)
            REQUIRE(ledgerRowEqual(expected[i],actual[i]));
        REQUIRE(actual.size()==65);
        requirePlacementLedgerOrder(referencePlacement,expected);
        requirePlacementLedgerOrder(placement,actual);
        REQUIRE(outcome.retainedResult->metrics.rounds==11);
    }

    const auto rounds=outcome.retainedResult?
        outcome.retainedResult->metrics.rounds:0;
    const bool productionGo=outcome.published()&&
        sink.allBlockersPairwiseIndependent&&rounds<=4&&
        sink.workAccountingComplete&&sink.memoryAccountingComplete&&
        fixedPoints!=0&&fixedPoints==replaced;
    REQUIRE(!productionGo);
    std::printf(
        "lineage structural=%s status=%u attempts=%zu fixed_points=%zu "
        "frontiers=%zu recoveries=%zu blockers=%zu rounds=%llu "
        "provider_events=%zu coordinator_lb=%llu nested_lb=%llu "
        "work_complete=0 memory_complete=0 production_go=%s\n",
        failures==0?"PASS":"FAIL",static_cast<unsigned>(sink.status),
        sink.attempts.size(),fixedPoints,replaced,sink.recoveries.size(),
        sink.blockers.size(),static_cast<unsigned long long>(rounds),
        placement.callsUsed,
        static_cast<unsigned long long>(sink.coordinatorWork),
        static_cast<unsigned long long>(sink.nestedRecoveryWork),
        productionGo?"true":"false");
    std::printf("deferred_blockers=%zu deferred_recoveries=%zu\n",
                deferredBlockers,deferredRecoveries);
    return failures==0?0:1;
}
