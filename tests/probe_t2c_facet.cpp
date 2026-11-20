// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Focused T2-C facet-conflict production-route allocation qualification.
// The real eight-layer route is intentionally used, but allocation tracking
// is restricted to the SCC-specific PAACS_TOPO_TESTING zones, including the
// facet-expanded owner/report storage secured before any provider call.

#include "t2c_facet_probe_fixture.hpp"

#include "../src/algo/T2CSurfaceCache.hpp"
#include "paacs/mesh/MeshAudit.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <new>
#include <set>
#include <span>
#include <string_view>
#include <vector>

namespace cache=paacs::algo::t2c_detail;
namespace fixture=paacs::tests::t2c_facet_probe;

namespace {
long gOrdinal=0;
long gFailOrdinal=0;
long gTailAllocations=0;
bool gTrack=false;
bool gInjected=false;
bool gProviderStarted=false;
cache::SurfaceCacheAllocationZone gInjectedZone=
    cache::SurfaceCacheAllocationZone::None;
constexpr std::size_t kZoneCount=static_cast<std::size_t>(
    cache::SurfaceCacheAllocationZone::Count);
std::array<long,kZoneCount> gZoneAllocations{};

bool selected(cache::SurfaceCacheAllocationZone zone) noexcept {
    return zone==cache::SurfaceCacheAllocationZone::FacetConflictStorage||
           zone==cache::SurfaceCacheAllocationZone::
               FacetConflictPreparedStateStorage||
           zone==cache::SurfaceCacheAllocationZone::FacetConflictPreparation||
           zone==cache::SurfaceCacheAllocationZone::FacetConflictProviderBatch;
}
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t bytes) {
    const auto zone=cache::surfaceCacheAllocationZoneForTesting();
    if(gTrack&&gProviderStarted&&
       (zone==cache::SurfaceCacheAllocationZone::ProposedStateIntegration||
        zone==cache::SurfaceCacheAllocationZone::CommitTail))
        ++gTailAllocations;
    if(gTrack&&selected(zone)) {
        ++gOrdinal;
        ++gZoneAllocations[static_cast<std::size_t>(zone)];
        if(gFailOrdinal>0&&gOrdinal==gFailOrdinal) {
            gInjected=true;
            gInjectedZone=zone;
            throw std::bad_alloc{};
        }
    }
    if(void* memory=std::malloc(bytes?bytes:1))return memory;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete(void* memory,std::size_t) noexcept { std::free(memory); }
void operator delete[](void* memory,std::size_t) noexcept {
    std::free(memory);
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
int failures=0;
#define REQUIRE(condition) do { \
    if(!(condition)) { \
        std::printf("FAIL %d: %s\n",__LINE__,#condition); \
        ++failures; \
    } \
} while(false)

class Digest {
public:
    void scalar(std::uint64_t value) noexcept {
        for(unsigned byte=0;byte<8;++byte) {
            value_^=(value>>(byte*8U))&UINT64_C(0xff);
            value_*=UINT64_C(1099511628211);
        }
    }
    void floating(double value) noexcept {
        scalar(std::bit_cast<std::uint64_t>(value));
    }
    std::uint64_t value() const noexcept { return value_; }
private:
    std::uint64_t value_=UINT64_C(1469598103934665603);
};

void digestIdentity(Digest& digest,
                    const paacs::algo::SurfaceNodeIdentity& identity) noexcept {
    digest.scalar(static_cast<std::uint64_t>(identity.domain));
    digest.scalar(identity.parentArtifact.persistentId);
    digest.scalar(identity.parentArtifact.revision);
    for(const auto node:identity.parentFacetGlobalNodes)
        digest.scalar(static_cast<std::uint64_t>(node));
    digest.scalar(identity.parentOrdinal);
    digest.scalar(identity.parentSubordinal);
    digest.scalar(identity.parameterUBits);
    digest.scalar(identity.parameterVBits);
    digest.scalar(identity.dyadicNumerator);
    digest.scalar(identity.dyadicLevel);
}

std::uint64_t ledgerDigest(
        std::span<const paacs::algo::ProtectedFeatureChildIdentity> rows)
        noexcept {
    Digest digest;
    digest.scalar(rows.size());
    for(const auto& row:rows) {
        digest.scalar(row.parentEntityPersistentId);
        digestIdentity(digest,row.childIdentity);
        digest.scalar(row.splitOrdinal);
        digest.scalar(static_cast<std::uint64_t>(
            static_cast<std::int64_t>(row.edgeTag)));
        digest.scalar(static_cast<std::uint64_t>(
            static_cast<std::int64_t>(row.faceTag)));
        digest.floating(row.parentParameterLo);
        digest.floating(row.parentParameterHi);
        digest.floating(row.midpointParameter);
        digest.floating(row.faceParameter.u);
        digest.floating(row.faceParameter.v);
        digest.floating(row.coordinate.x);
        digest.floating(row.coordinate.y);
        digest.floating(row.coordinate.z);
        digest.scalar(static_cast<std::uint64_t>(
            static_cast<std::int64_t>(row.classification.dim)));
        digest.scalar(static_cast<std::uint64_t>(
            static_cast<std::int64_t>(row.classification.tag)));
        digest.scalar(row.providerStamp.persistentId);
        digest.scalar(row.providerStamp.revision);
    }
    return digest.value();
}

struct RouteFixture {
    fixture::LayeredModel model;
    fixture::Ownership ownership;
    fixture::Placement placement;
    paacs::sizing::MeshSpec spec;
    paacs::algo::BoundaryConformingV4Request request;
    paacs::mesh::Mesh destination;
    std::uint64_t destinationBefore=0;

    RouteFixture() {
        fixture::bindOwnership(model,ownership);
        placement.providerStarted=&gProviderStarted;
        spec.field=std::make_shared<fixture::LayerField>();
        spec.threads=1;
        request=fixture::requestFor(ownership,placement);
        (void)destination.addNode({-7,-11,-13},{0,91});
        destinationBefore=paacs::mesh::meshHash(destination);
    }
    void resetProvider() noexcept {
        placement.callsUsed=0;
        placement.overflow=false;
        gProviderStarted=false;
    }
};

paacs::algo::BoundaryConformingV4Outcome run(RouteFixture& fixture,
                                              long ordinal) {
    gOrdinal=0;
    gFailOrdinal=ordinal;
    gInjected=false;
    gInjectedZone=cache::SurfaceCacheAllocationZone::None;
    gZoneAllocations.fill(0);
    gTailAllocations=0;
    gProviderStarted=false;
    struct TrackGuard {
        TrackGuard() noexcept { gTrack=true; }
        ~TrackGuard() noexcept { gTrack=false;gFailOrdinal=0; }
    } track;
    return paacs::algo::meshModelBoundaryConformingV4(
        fixture.model,fixture.spec,fixture.request,fixture.destination,
        fixture.ownership,fixture.placement);
}

bool providerPrefix(std::span<const fixture::ProviderCall> prefix,
                    std::span<const fixture::ProviderCall> complete) {
    return prefix.size()<=complete.size()&&
           std::equal(prefix.begin(),prefix.end(),complete.begin());
}

void requireCacheXor(const RouteFixture& fixture,
                     const paacs::algo::SurfaceCacheReport& report) {
    REQUIRE(report.reportComplete);
    REQUIRE(report.processed.size()+report.unattempted.size()==
            fixture.model.faces().size());
    std::set<paacs::algo::FacePersistentId> observed;
    for(const auto& row:report.processed)
        REQUIRE(observed.insert(row.facePersistentId).second);
    for(const auto face:report.unattempted)
        REQUIRE(observed.insert(face).second);
    REQUIRE(observed.size()==fixture.model.faces().size());
}

const paacs::algo::BoundaryConformingV3Outcome* retainedSource(
        const paacs::algo::BoundaryConformingV4Outcome& outcome) noexcept {
    if(outcome.lastSourceOutcome)return outcome.lastSourceOutcome.get();
    if(outcome.retainedResult&&outcome.retainedResult->sourceOutcome)
        return outcome.retainedResult->sourceOutcome.get();
    return nullptr;
}

void requireSourceCheckpoint(
        const paacs::algo::BoundaryConformingV3Outcome* source) {
    REQUIRE(source!=nullptr);
    if(!source)return;
    REQUIRE(source->reportComplete);
    REQUIRE(source->report!=nullptr);
    REQUIRE(source->surfaceOutcome!=nullptr);
    REQUIRE(source->surfaceOutcome&&
            source->surfaceOutcome->retainedArtifact!=nullptr);
    REQUIRE(source->failure.reason==
            paacs::algo::BoundaryConformingV3Reason::RecoveryConflictCycle);
    REQUIRE(source->failure.workKind==paacs::algo::VolumeWorkKind::Facet);
    REQUIRE((source->failure.facetKey==paacs::bnd::VolumeFacetKey{40,4}));
    REQUIRE(source->failure.scheduleOrdinal==1109);
    REQUIRE(source->failure.exactWitnessDigest==
            UINT64_C(12033182254251637294));
}

struct Reference {
    std::uint64_t strict=0;
    std::uint64_t semantic=0;
    bool semanticCanonical=false;
    std::uint64_t fullLedger=0;
    std::uint64_t committedPrefix=0;
    std::vector<fixture::ProviderCall> calls;
    std::array<long,kZoneCount> zoneCounts{};
    long total=0;
};

Reference measure() {
    RouteFixture fixture;
    const auto outcome=run(fixture,0);
    REQUIRE(outcome.published());
    REQUIRE(!fixture.placement.overflow);
    REQUIRE(gTailAllocations==0);
    REQUIRE(outcome.retainedResult!=nullptr);
    Reference reference;
    reference.total=gOrdinal;
    reference.zoneCounts=gZoneAllocations;
    reference.strict=paacs::mesh::meshHash(fixture.destination);
    const auto semantic=paacs::mesh::semanticHashEx(fixture.destination);
    reference.semantic=semantic.value;
    reference.semanticCanonical=semantic.canonical;
    reference.calls.assign(fixture.placement.trace().begin(),
                           fixture.placement.trace().end());
    REQUIRE(reference.calls.size()==196);
    if(reference.calls.size()==196) {
        for(std::size_t placement=0;placement<65;++placement) {
            REQUIRE(reference.calls[placement*3].kind==
                    fixture::ProviderCallKind::Stamp);
            REQUIRE(reference.calls[placement*3+1].kind==
                    fixture::ProviderCallKind::Placement);
            REQUIRE(reference.calls[placement*3+2].kind==
                    fixture::ProviderCallKind::Stamp);
        }
        REQUIRE(reference.calls[195].kind==
                fixture::ProviderCallKind::Stamp);
        REQUIRE(reference.calls[190].entityKind==
                paacs::algo::ProtectedFeatureEntityKind::FacetInterior);
        REQUIRE(reference.calls[193].entityKind==
                paacs::algo::ProtectedFeatureEntityKind::FacetInterior);
    }
    if(outcome.retainedResult) {
        const auto& ledger=outcome.retainedResult->subdivisionLedger;
        reference.fullLedger=ledgerDigest(ledger);
        REQUIRE(ledger.size()==65);
        if(ledger.size()>=63)
            reference.committedPrefix=ledgerDigest(
                std::span<const paacs::algo::ProtectedFeatureChildIdentity>(
                    ledger.data(),63));
        requireCacheXor(fixture,outcome.retainedResult->surfaceCacheReport);
    }
    REQUIRE(reference.strict==UINT64_C(2737597274546495434));
    REQUIRE(reference.semantic==UINT64_C(16753220727311122020));
    REQUIRE(reference.semanticCanonical);
    REQUIRE(reference.zoneCounts[static_cast<std::size_t>(
                cache::SurfaceCacheAllocationZone::FacetConflictStorage)]==2);
    REQUIRE(reference.zoneCounts[static_cast<std::size_t>(
                cache::SurfaceCacheAllocationZone::
                    FacetConflictPreparedStateStorage)]>0);
    REQUIRE(reference.zoneCounts[static_cast<std::size_t>(
                cache::SurfaceCacheAllocationZone::FacetConflictPreparation)]==2);
    REQUIRE(reference.zoneCounts[static_cast<std::size_t>(
                cache::SurfaceCacheAllocationZone::FacetConflictProviderBatch)]==0);
    REQUIRE(reference.total==
        reference.zoneCounts[static_cast<std::size_t>(
            cache::SurfaceCacheAllocationZone::FacetConflictStorage)]+
        reference.zoneCounts[static_cast<std::size_t>(
            cache::SurfaceCacheAllocationZone::
                FacetConflictPreparedStateStorage)]+
        reference.zoneCounts[static_cast<std::size_t>(
            cache::SurfaceCacheAllocationZone::FacetConflictPreparation)]+
        reference.zoneCounts[static_cast<std::size_t>(
            cache::SurfaceCacheAllocationZone::FacetConflictProviderBatch)]);
    return reference;
}

void requireCleanRetry(RouteFixture& fixture,const Reference& reference) {
    fixture.resetProvider();
    const auto retry=run(fixture,0);
    REQUIRE(retry.published());
    REQUIRE(paacs::mesh::meshHash(fixture.destination)==reference.strict);
    const auto semantic=paacs::mesh::semanticHashEx(fixture.destination);
    REQUIRE(semantic.value==reference.semantic);
    REQUIRE(semantic.canonical==reference.semanticCanonical);
    REQUIRE(!fixture.placement.overflow);
    REQUIRE(std::equal(fixture.placement.trace().begin(),
                       fixture.placement.trace().end(),
                       reference.calls.begin(),reference.calls.end()));
    REQUIRE(gTailAllocations==0);
    REQUIRE(retry.retainedResult!=nullptr);
    if(retry.retainedResult) {
        REQUIRE(ledgerDigest(retry.retainedResult->subdivisionLedger)==
                reference.fullLedger);
        requireCacheXor(fixture,retry.retainedResult->surfaceCacheReport);
    }
}

void sweep(const Reference& reference) {
    for(long ordinal=1;ordinal<=reference.total;++ordinal) {
        RouteFixture fixture;
        const auto failed=run(fixture,ordinal);
        REQUIRE(gInjected);
        REQUIRE(failed.status==paacs::algo::BoundaryConformingV4Status::Failed);
        REQUIRE(failed.reason==
                paacs::algo::BoundaryConformingV4Reason::AllocationFailed);
        REQUIRE(failed.destinationUnchanged);
        REQUIRE(paacs::mesh::meshHash(fixture.destination)==
                fixture.destinationBefore);
        REQUIRE(failed.retainedResult!=nullptr);
        REQUIRE(!fixture.placement.overflow);
        REQUIRE(providerPrefix(fixture.placement.trace(),reference.calls));
        REQUIRE(gTailAllocations==0);
        if(failed.retainedResult) {
            requireCacheXor(fixture,
                            failed.retainedResult->surfaceCacheReport);
            const auto& ledger=failed.retainedResult->subdivisionLedger;
            if(gInjectedZone==
                   cache::SurfaceCacheAllocationZone::FacetConflictStorage||
               gInjectedZone==cache::SurfaceCacheAllocationZone::
                   FacetConflictPreparedStateStorage) {
                REQUIRE(ledger.empty());
                REQUIRE(fixture.placement.trace().empty());
                REQUIRE(retainedSource(failed)==nullptr);
            } else {
                REQUIRE(ledger.size()==63);
                REQUIRE(ledgerDigest(ledger)==reference.committedPrefix);
                requireSourceCheckpoint(retainedSource(failed));
            }
        }
        requireCleanRetry(fixture,reference);
    }
}

} // namespace

int main(int argc,char** argv) {
    const Reference reference=measure();
    std::printf("probe_t2c_facet allocations=%ld storage=%ld prepared=%ld "
                "preparation=%ld provider=%ld tail=%ld\n",reference.total,
        reference.zoneCounts[static_cast<std::size_t>(
            cache::SurfaceCacheAllocationZone::FacetConflictStorage)],
        reference.zoneCounts[static_cast<std::size_t>(
            cache::SurfaceCacheAllocationZone::
                FacetConflictPreparedStateStorage)],
        reference.zoneCounts[static_cast<std::size_t>(
            cache::SurfaceCacheAllocationZone::FacetConflictPreparation)],
        reference.zoneCounts[static_cast<std::size_t>(
            cache::SurfaceCacheAllocationZone::FacetConflictProviderBatch)],
        gTailAllocations);
    const bool measureOnly=argc==2&&std::string_view(argv[1])=="--measure";
    if(argc>2||(argc==2&&!measureOnly)) {
        std::puts("usage: paacs_probe_t2c_facet [--measure]");
        return 2;
    }
    if(!measureOnly)sweep(reference);
    if(failures) {
        std::printf("probe_t2c_facet FAILED assertions=%d\n",failures);
        return 1;
    }
    std::puts("probe_t2c_facet PASSED");
    return 0;
}
