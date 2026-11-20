// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// T5 route/publication allocation qualification. Excluded from sanitizer
// builds because it replaces global operator new.
#include "paacs/algo/MultiRegionMesher.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sizing/SizeField.hpp"
#include "paacs/topo/EdgeTopology.hpp"
#include "paacs/topo/FaceTopology.hpp"

#include <algorithm>
#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <map>
#include <memory>
#include <new>
#include <stdexcept>
#include <string_view>

static long gOrdinal=0,gFail=0,gTailAllocations=0;
static bool gTrack=false,gInjected=false,gPublicationOnly=false,
            gSnapshotOnly=false;

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t bytes) {
    if(paacs::algo::t5_detail::publicationTailActive())++gTailAllocations;
    const bool selected=gSnapshotOnly||
        (gPublicationOnly
            ?paacs::algo::t5_detail::publicationPreparationActive()
            :paacs::algo::t5_detail::routeAllocationActive());
    if(gTrack&&selected) {
        ++gOrdinal;
        if(gFail>0&&gOrdinal==gFail) {
            gInjected=true;throw std::bad_alloc{};
        }
    }
    if(void* memory=std::malloc(bytes?bytes:1))return memory;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t bytes){return ::operator new(bytes);}
void operator delete(void* memory) noexcept{std::free(memory);}
void operator delete[](void* memory) noexcept{std::free(memory);}
void operator delete(void* memory,std::size_t) noexcept{std::free(memory);}
void operator delete[](void* memory,std::size_t) noexcept{std::free(memory);}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
using namespace paacs;
int failures=0;
#define REQUIRE(x) do{if(!(x)){std::printf("FAIL %d: %s\n",__LINE__,#x);++failures;}}while(0)

class Ownership final : public algo::VolumeBoundaryOwnershipProvider {
public:
    algo::VolumeBoundaryProviderStamp current{19001,1};
    std::map<algo::FacePersistentId,std::vector<bnd::VolumeRegionUse>> rows;
    algo::VolumeBoundaryProviderStamp stamp() const noexcept override {
        return current;
    }
    algo::VolumeFaceOwnership ownership(
            algo::FacePersistentId face) const override {
        algo::VolumeFaceOwnership result;result.providerStamp=current;
        result.facePersistentId=face;const auto found=rows.find(face);
        if(found!=rows.end())result.owners=found->second;return result;
    }
    algo::ProviderRecoveryEnvelope recoveryEnvelope(
            std::span<const Vec3>) const override {
        algo::ProviderRecoveryEnvelope result;result.providerStamp=current;
        result.availability=
            algo::VolumeProviderAvailability::RecoveryEnvelopeUnavailable;
        return result;
    }
};

void configure(const geo::Model& model,Ownership& provider) {
    Vec3 lo{},hi{};model.solids().front()->bbox(lo,hi);
    const Vec3 centre=(lo+hi)*0.5;
    for(const auto* face:model.faces()) {
        Vec2 uvLo{},uvHi{};face->paramBounds(uvLo,uvHi);
        const Vec2 uv=(uvLo+uvHi)*0.5;
        const int sense=face->normalAt(uv).dot(
            face->pointAt(uv)-centre)>0?1:-1;
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            {{0,sense}};
    }
}

struct Fixture {
    std::unique_ptr<geo::Kernel> kernel;
    std::unique_ptr<geo::Model> model;
    Ownership provider;
    sizing::MeshSpec spec;
    algo::T5RegionBinding region{0,19002,19003,{3,0}};
    algo::T5Request request;
    Fixture() {
        kernel=geo::makeBuiltinKernel();
        model=kernel->makeBox({0,0,0},{1,1,1});configure(*model,provider);
        spec.field=std::make_shared<sizing::ConstantField>(0.55);spec.threads=1;
        request.sourceRequest.expectedSurfaceSetStamp={19004,1};
        request.sourceRequest.expectedOwnershipProviderStamp=provider.current;
        request.regionBindings=std::span<const algo::T5RegionBinding>(&region,1);
    }
};

struct Signature {
    std::uint64_t strictHash=0,semanticHash=0,authority=0,processed=0;
    bool canonical=false;
    friend bool operator==(const Signature&,const Signature&)=default;
};
Signature signature(const algo::T5Outcome& outcome) {
    Signature result;
    if(outcome.retainedResult) {
        result.strictHash=outcome.retainedResult->strictHash;
        result.semanticHash=outcome.retainedResult->semanticHash.value;
        result.canonical=outcome.retainedResult->semanticHash.canonical;
        for(const auto value:outcome.retainedResult->semanticAuthorityCoreDigest)
            result.authority=result.authority*131U+value;
    }
    if(outcome.report)result.processed=outcome.report->processedUsed;
    return result;
}

algo::T5Outcome run(Fixture& fixture,mesh::Mesh& destination,long failure,
                    bool publicationOnly) {
    gOrdinal=0;gFail=failure;gInjected=false;
    gPublicationOnly=publicationOnly;gSnapshotOnly=false;gTrack=true;
    auto outcome=algo::meshModelMultiRegionTet4(
        *fixture.model,fixture.spec,fixture.request,destination,
        fixture.provider);
    gTrack=false;gFail=0;return outcome;
}

long measure(bool publicationOnly,Signature& reference) {
    {Fixture warm;mesh::Mesh mesh;auto outcome=algo::meshModelMultiRegionTet4(
        *warm.model,warm.spec,warm.request,mesh,warm.provider);
     REQUIRE(outcome.published());}
    Fixture fixture;mesh::Mesh mesh;
    auto outcome=run(fixture,mesh,0,publicationOnly);
    REQUIRE(outcome.published());REQUIRE(outcome.retainedResult!=nullptr);
    reference=signature(outcome);REQUIRE(gTailAllocations==0);
    return gOrdinal;
}

void sweep(bool publicationOnly,long first,long last,long expected) {
    Signature reference;const long total=measure(publicationOnly,reference);
    REQUIRE(expected<=0||expected==total);if(last<=0||last>total)last=total;
    long faults=0;
    for(long ordinal=first;ordinal<=last;++ordinal) {
        Fixture fixture;mesh::Mesh destination;
        destination.addNode({-31,-37,-41},{0,991});
        const auto before=mesh::meshHash(destination);
        auto outcome=run(fixture,destination,ordinal,publicationOnly);
        REQUIRE(gInjected);REQUIRE(!outcome.published());
        REQUIRE(outcome.destinationUnchanged);
        REQUIRE(mesh::meshHash(destination)==before);
        REQUIRE(outcome.failure.reason==algo::T5Reason::AllocationFailed||
                outcome.failure.reason==algo::T5Reason::SourceFailed);
        if(publicationOnly) {
            REQUIRE(outcome.sourceReceipt!=nullptr);
            REQUIRE(outcome.retainedResult!=nullptr);
        }
        Fixture retry;auto clean=run(retry,destination,0,publicationOnly);
        REQUIRE(clean.published());REQUIRE(signature(clean)==reference);
        ++faults;
    }
    std::printf("probe_t5 %s range=%ld-%ld total=%ld faults=%ld tail=%ld\n",
        publicationOnly?"publication":"route",first,last,total,faults,
        gTailAllocations);
}

void snapshotProof() {
    mesh::Mesh mesh;mesh.reserveAppend(8,2,8);
    gOrdinal=0;gFail=1;gInjected=false;gSnapshotOnly=true;gTrack=true;
    const auto first=mesh.capacitySnapshot();
    const auto second=mesh.capacitySnapshot();
    gTrack=false;gSnapshotOnly=false;gFail=0;
    REQUIRE(first.complete&&second.complete);REQUIRE(gOrdinal==0);
    REQUIRE(!gInjected);

    Fixture fixture;mesh::Mesh destination;
    auto outcome=run(fixture,destination,0,false);
    REQUIRE(outcome.published());REQUIRE(outcome.retainedResult!=nullptr);
    REQUIRE(outcome.retainedResult->authoritative!=nullptr);
    auto& adjacency=outcome.retainedResult->authoritative->workspace->adjacency();
    const auto expected=adjacency.capacitySnapshot();
    gOrdinal=0;gFail=1;gInjected=false;gSnapshotOnly=true;gTrack=true;
    const auto adjacencyFirst=adjacency.capacitySnapshot();
    const auto adjacencySecond=adjacency.capacitySnapshot();
    gTrack=false;gSnapshotOnly=false;gFail=0;
    const auto same=[](const auto& a,const auto& b) noexcept {
        return a.entryBucketCount==b.entryBucketCount&&
               a.entryCount==b.entryCount&&
               a.entryFaceCapacity==b.entryFaceCapacity&&
               a.entryNeighborCapacity==b.entryNeighborCapacity&&
               a.boundaryCapacity==b.boundaryCapacity&&
               a.interiorCapacity==b.interiorCapacity&&
               a.nonManifoldCapacity==b.nonManifoldCapacity&&
               a.nonManifoldCellCapacity==b.nonManifoldCellCapacity&&
               a.complete==b.complete;
    };
    REQUIRE(expected.complete);REQUIRE(same(expected,adjacencyFirst));
    REQUIRE(same(expected,adjacencySecond));REQUIRE(gOrdinal==0);
    REQUIRE(!gInjected);
}

void completeFactorySweep() {
    topo::NodePool nodes;
    std::vector<topo::NodeHandle> n;
    for(int i=0;i<6;++i)
        n.push_back(nodes.add({double(i),double(i%2),double(i%3)},{3,1}));
    topo::CellPool cells;
    cells.add(topo::ElementType::Tet4,{n[0],n[1],n[2],n[3]},
              {3,1},nodes);
    cells.add(topo::ElementType::Tet4,{n[0],n[2],n[1],n[4]},
              {3,1},nodes);
    cells.add(topo::ElementType::Tet4,{n[0],n[1],n[2],n[5]},
              {3,1},nodes);
    const auto sourceHash=topo::contentHash(nodes,cells).value;

    gOrdinal=0;gFail=0;gInjected=false;gSnapshotOnly=true;gTrack=true;
    topo::EdgePoolCompleteBuildMetrics edgeReference;
    auto completeEdges=topo::EdgePool::buildComplete(
        nodes,cells,&edgeReference);
    gTrack=false;gSnapshotOnly=false;
    const long edgeTotal=gOrdinal;
    REQUIRE(edgeReference.runComplete);
    REQUIRE(edgeReference.validationCellVisits==3);
    REQUIRE(edgeReference.buildCellVisits==3);
    REQUIRE(edgeReference.edgeUses==18);
    REQUIRE(edgeReference.finalEdgeAudits==1);
    REQUIRE(edgeReference.perCellGlobalAudits==0);
    REQUIRE(topo::edgeTopologyReady(nodes,cells,completeEdges));

    long edgeFaults=0;
    for(long ordinal=1;ordinal<=edgeTotal;++ordinal) {
        gOrdinal=0;gFail=ordinal;gInjected=false;
        gSnapshotOnly=true;gTrack=true;bool threw=false;
        try {
            auto scratch=topo::EdgePool::buildComplete(nodes,cells);
            (void)scratch;
        } catch(const std::bad_alloc&) { threw=true; }
        gTrack=false;gSnapshotOnly=false;gFail=0;
        REQUIRE(threw);REQUIRE(gInjected);
        REQUIRE(topo::contentHash(nodes,cells).value==sourceHash);
        REQUIRE(topo::storageEdgeReady(nodes,cells));
        topo::EdgePoolCompleteBuildMetrics retryMetrics;
        auto retry=topo::EdgePool::buildComplete(nodes,cells,&retryMetrics);
        REQUIRE(retryMetrics.runComplete);
        REQUIRE(topo::edgeTopologyReady(nodes,cells,retry));
        ++edgeFaults;
    }

    const auto edgeHandlesBefore=completeEdges.activeBySequence();
    const auto edgeStatsBefore=completeEdges.stats();
    gOrdinal=0;gFail=0;gInjected=false;gSnapshotOnly=true;gTrack=true;
    topo::FacePoolCompleteBuildMetrics faceReference;
    auto completeFaces=topo::FacePool::buildComplete(
        nodes,cells,completeEdges,&faceReference);
    gTrack=false;gSnapshotOnly=false;
    const long faceTotal=gOrdinal;
    REQUIRE(faceReference.runComplete);
    REQUIRE(faceReference.validationCellVisits==3);
    REQUIRE(faceReference.buildCellVisits==3);
    REQUIRE(faceReference.faceUses==12);
    REQUIRE(faceReference.suppliedEdgeAudits==1);
    REQUIRE(faceReference.finalFaceAudits==1);
    REQUIRE(faceReference.perCellGlobalAudits==0);
    REQUIRE(topo::faceTopologyReady(nodes,cells,completeEdges,completeFaces));

    long faceFaults=0;
    for(long ordinal=1;ordinal<=faceTotal;++ordinal) {
        gOrdinal=0;gFail=ordinal;gInjected=false;
        gSnapshotOnly=true;gTrack=true;bool threw=false;
        try {
            auto scratch=topo::FacePool::buildComplete(
                nodes,cells,completeEdges);
            (void)scratch;
        } catch(const std::bad_alloc&) { threw=true; }
        gTrack=false;gSnapshotOnly=false;gFail=0;
        REQUIRE(threw);REQUIRE(gInjected);
        REQUIRE(topo::contentHash(nodes,cells).value==sourceHash);
        REQUIRE(completeEdges.activeBySequence()==edgeHandlesBefore);
        const auto edgeStatsAfter=completeEdges.stats();
        REQUIRE(edgeStatsAfter.live==edgeStatsBefore.live);
        REQUIRE(edgeStatsAfter.totalCreated==edgeStatsBefore.totalCreated);
        REQUIRE(topo::edgeTopologyReady(nodes,cells,completeEdges));
        topo::FacePoolCompleteBuildMetrics retryMetrics;
        auto retry=topo::FacePool::buildComplete(
            nodes,cells,completeEdges,&retryMetrics);
        REQUIRE(retryMetrics.runComplete);
        REQUIRE(topo::faceTopologyReady(nodes,cells,completeEdges,retry));
        ++faceFaults;
    }
    std::printf("t5 complete factories edge=%ld/%ld face=%ld/%ld\n",
                edgeFaults,edgeTotal,faceFaults,faceTotal);
}

long number(const char* text) {
    long value=0;const char* end=text+std::strlen(text);
    const auto parsed=std::from_chars(text,end,value);
    if(parsed.ec!=std::errc{}||parsed.ptr!=end||value<0)
        throw std::invalid_argument("ordinal");
    return value;
}
bool publicationMode(std::string_view value) {
    if(value=="publication")return true;if(value=="route")return false;
    throw std::invalid_argument("mode");
}
} // namespace

int main(int argc,char** argv) {
    try {
        snapshotProof();
        if(argc==2&&std::string_view(argv[1])=="--snapshot") {
            // snapshotProof() already ran before argument dispatch.
        } else if(argc==2&&std::string_view(argv[1])=="--factories") {
            completeFactorySweep();
        } else if(argc==3&&std::string_view(argv[1])=="--count") {
            const bool publication=publicationMode(argv[2]);Signature value;
            std::printf("t5 %s allocation count: %ld\n",argv[2],
                        measure(publication,value));
        } else if(argc==4&&std::string_view(argv[1])=="--expect-count") {
            const bool publication=publicationMode(argv[2]);Signature value;
            const long observed=measure(publication,value);
            REQUIRE(observed==number(argv[3]));
        } else if(argc==5&&std::string_view(argv[1])=="--range") {
            sweep(publicationMode(argv[2]),number(argv[3]),number(argv[4]),0);
        } else if(argc==6&&std::string_view(argv[1])=="--expect-range") {
            sweep(publicationMode(argv[2]),number(argv[3]),number(argv[4]),
                  number(argv[5]));
        } else if(argc==1) {
            completeFactorySweep();
            Signature value;const long total=measure(false,value);
            sweep(false,1,(std::min)(8L,total),total);
            const long publication=measure(true,value);
            sweep(true,1,(std::min)(8L,publication),publication);
        } else return 2;
    } catch(const std::exception& error) {
        gTrack=false;std::printf("probe_t5 exception: %s\n",error.what());
        ++failures;
    }
    REQUIRE(gTailAllocations==0);
    if(failures){std::printf("probe_t5 FAILED failures=%d\n",failures);return 1;}
    std::puts("probe_t5 PASSED");return 0;
}
