// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// T4-A standalone allocation qualification. This target replaces global
// operator new and is intentionally excluded from sanitizer builds.
#include <algorithm>
#include <array>
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
#include <vector>

#if defined(_WIN32)
#define NOMINMAX
#include <Windows.h>
#include <intrin.h>
#else
#include <unistd.h>
#endif

#include "paacs/algo/TetQualityMesher.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sizing/SizeField.hpp"

static long gAllocationOrdinal=0;
static long gFailOrdinal=0;
static long gPhase3Allocations=0;
static long gPublicationTailAllocations=0;
static bool gTrackAllocations=false;
static bool gFailureInjected=false;

namespace {
void terminateDiagnostic() noexcept {
    char buffer[192];
    char* cursor=buffer;
    auto literal=[&](const char* value) noexcept {
        const auto count=std::strlen(value);
        std::memcpy(cursor,value,count);
        cursor+=count;
    };
    auto number=[&](long value) noexcept {
        const auto converted=std::to_chars(cursor,buffer+sizeof(buffer),value);
        if (converted.ec==std::errc{}) cursor=converted.ptr;
    };
    literal("probe_t4a terminate ordinal="); number(gAllocationOrdinal);
    literal(" fail="); number(gFailOrdinal);
    literal(" phase3="); number(gPhase3Allocations);
    literal(" publication="); number(gPublicationTailAllocations);
    literal("\n");
#if defined(_WIN32)
    DWORD written=0;
    (void)::WriteFile(::GetStdHandle(STD_ERROR_HANDLE),buffer,
                      static_cast<DWORD>(cursor-buffer),&written,nullptr);
    __fastfail(7);
#else
    (void)::write(STDERR_FILENO,buffer,
                  static_cast<std::size_t>(cursor-buffer));
    std::_Exit(134);
#endif
}
} // namespace

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t bytes) {
    if (paacs::topo::detail::phase3Active()) ++gPhase3Allocations;
    if (paacs::algo::boundary_conforming_detail::publicationTailActive())
        ++gPublicationTailAllocations;
    if (gTrackAllocations&&
        paacs::algo::t4a_detail::allocationScopeForTesting()) {
        ++gAllocationOrdinal;
        if (gFailOrdinal>0&&gAllocationOrdinal==gFailOrdinal) {
            gFailureInjected=true;
            throw std::bad_alloc{};
        }
    }
    if (void* memory=std::malloc(bytes?bytes:1)) return memory;
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
int failures=0;
#define REQUIRE(condition)                                                   \
    do {                                                                     \
        if (!(condition)) {                                                  \
            std::printf("FAIL %d: %s\n",__LINE__,#condition);              \
            ++failures;                                                      \
        }                                                                    \
    } while (0)

class Provider final : public algo::VolumeBoundaryOwnershipProvider {
public:
    algo::VolumeBoundaryProviderStamp current{2701,1};
    std::map<algo::FacePersistentId,std::vector<bnd::VolumeRegionUse>> rows;
    algo::VolumeBoundaryProviderStamp stamp() const noexcept override {
        return current;
    }
    algo::VolumeFaceOwnership ownership(
            algo::FacePersistentId face) const override {
        algo::VolumeFaceOwnership result;
        result.providerStamp=current;
        result.facePersistentId=face;
        const auto found=rows.find(face);
        if (found!=rows.end()) result.owners=found->second;
        return result;
    }
    algo::ProviderRecoveryEnvelope recoveryEnvelope(
            std::span<const Vec3>) const override {
        algo::ProviderRecoveryEnvelope result;
        result.providerStamp=current;
        result.availability=
            algo::VolumeProviderAvailability::RecoveryEnvelopeUnavailable;
        return result;
    }
};

void configureProvider(const geo::Model& model,Provider& provider) {
    const auto solids=model.solids();
    if (solids.size()!=1) throw std::logic_error("probe requires one solid");
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

algo::T4ARequest requestFor(const Provider& provider,bool publication) {
    algo::T4ARequest request;
    request.source.conformingRequest.expectedSurfaceSetStamp={27001,1};
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
    request.quality.minimumRadiusRatio=publication?0.0:0.2;
    request.quality.maxResultNodes=100000;
    request.quality.maxResultTets=100000;
    request.quality.maxReportRows=publication?1000:200000;
    request.quality.maxWorstTailRows=64;
    if (publication) {
        request.quality.mode=algo::T4AQualityMode::ReportOnly;
        return request;
    }
    request.quality.mode=algo::T4AQualityMode::BoundedImprove;
    request.quality.minimumImprovement=1.0e-8;
    request.quality.maxOperations=1;
    request.quality.maxSelectedGenerations=64;
    request.quality.maxCandidateAttemptsPerCell=1;
    request.quality.maxCreatedTetsPerOperation=128;
    request.quality.maxWorkUnits=1000000;
    request.quality.maxNewNodes=1;
    request.quality.maxDependencyRowsPerAttempt=1000;
    request.quality.maxProtectedWitnessesPerAttempt=1000;
    request.quality.maxTopologyRowsPerOperation=10000;
    return request;
}

struct Signature {
    algo::T4AStatus status=algo::T4AStatus::Failed;
    algo::T4AReason reason=algo::T4AReason::None;
    std::uint64_t attempts=0,blocked=0,commits=0,dependencies=0,
        unattempted=0,work=0,strictHash=0,semanticHash=0;
    bool canonical=false,reportComplete=false,runComplete=false;
    friend bool operator==(const Signature&,const Signature&)=default;
};

Signature signatureOf(const algo::T4AOutcome& outcome) {
    Signature result;
    result.status=outcome.status;
    result.reason=outcome.failure.reason;
    result.reportComplete=outcome.reportComplete;
    if (outcome.report) {
        result.attempts=outcome.report->attemptsUsed;
        result.blocked=outcome.report->blockedUsed;
        result.commits=outcome.report->commitsUsed;
        result.dependencies=outcome.report->dependenciesUsed;
        result.unattempted=outcome.report->unattemptedUsed;
    }
    if (outcome.retainedResult) {
        result.work=outcome.retainedResult->metrics.workUnitsConsumed;
        result.strictHash=outcome.retainedResult->finalStrictHash;
        result.semanticHash=outcome.retainedResult->finalSemanticHash.value;
        result.canonical=outcome.retainedResult->finalSemanticHash.canonical;
        result.runComplete=outcome.retainedResult->metrics.runComplete;
    }
    return result;
}

struct RouteFixture {
    std::unique_ptr<geo::Kernel> kernel;
    std::unique_ptr<geo::Model> model;
    Provider provider;
    sizing::MeshSpec spec;
    algo::T4ARequest request;
    explicit RouteFixture(bool publication) {
        kernel=geo::makeBuiltinKernel();
        model=kernel->makeBox({0,0,0},{1,1,1});
        configureProvider(*model,provider);
        spec.field=std::make_shared<sizing::ConstantField>(
            publication?0.5:0.4);
        spec.threads=1;
        request=requestFor(provider,publication);
    }
};

algo::T4AOutcome runTracked(RouteFixture& fixture,mesh::Mesh& destination,
                            long ordinal) {
    gAllocationOrdinal=0;
    gFailOrdinal=ordinal;
    gFailureInjected=false;
    gTrackAllocations=true;
    auto result=algo::meshModelQualityTet4(
        *fixture.model,fixture.spec,fixture.request,destination,
        fixture.provider);
    gTrackAllocations=false;
    gFailOrdinal=0;
    return result;
}

long measure(bool publication,Signature& signature) {
    // Exclude one-time runtime/static initialization from the repeatable
    // route ordinal space.  The injected runs below execute in this same
    // warmed process, so qualification must count that identical state.
    {
        RouteFixture warmFixture(publication);
        mesh::Mesh warmDestination;
        warmDestination.addNode({-41,-43,-47},{0,990});
        auto warm=algo::meshModelQualityTet4(
            *warmFixture.model,warmFixture.spec,warmFixture.request,
            warmDestination,warmFixture.provider);
        REQUIRE(warm.published()||warm.status==algo::T4AStatus::Incomplete);
    }
    RouteFixture fixture(publication);
    mesh::Mesh destination;
    destination.addNode({-41,-43,-47},{0,990});
    auto outcome=runTracked(fixture,destination,0);
    signature=signatureOf(outcome);
    REQUIRE(outcome.published()||
            outcome.status==algo::T4AStatus::Incomplete);
    REQUIRE(outcome.retainedResult!=nullptr);
    REQUIRE(publication||signature.commits==1);
    REQUIRE(outcome.retainedResult&&
            outcome.retainedResult->metrics.phase3Allocations==0);
    REQUIRE(outcome.retainedResult&&
            outcome.retainedResult->metrics.publicationTailAllocations==0);
    REQUIRE(gPhase3Allocations==0);
    REQUIRE(gPublicationTailAllocations==0);
    return gAllocationOrdinal;
}

void requireStructuredFailure(const algo::T4AOutcome& outcome,
                              std::uint64_t before,const mesh::Mesh& destination) {
    REQUIRE(outcome.destinationUnchanged);
    REQUIRE(mesh::meshHash(destination)==before);
    REQUIRE(outcome.status!=algo::T4AStatus::Published);
    if (outcome.failure.reason==algo::T4AReason::AllocationFailed&&
        outcome.failure.phase==algo::T4APhase::PublicationPreparation) {
        // T4A_DESIGN.md section 8 promises that publication-preparation
        // allocation failure retains the already-audited private result for
        // retry evidence.  The immutable source receipt must also survive.
        REQUIRE(outcome.sourceReceipt!=nullptr);
        REQUIRE(outcome.retainedResult!=nullptr);
        REQUIRE(outcome.retainedResult&&
                !outcome.retainedResult->publicationReady);
    }
    if (outcome.failure.reason!=algo::T4AReason::AllocationFailed) {
        REQUIRE(outcome.failure.reason==algo::T4AReason::SourceT3AStopped);
        REQUIRE(outcome.sourceOutcome!=nullptr);
        REQUIRE(outcome.sourceOutcome&&
                outcome.sourceOutcome->status==algo::T3AStatus::Failed);
        REQUIRE(outcome.sourceOutcome&&
                (outcome.sourceOutcome->failure.reason==
                     algo::T3AReason::AllocationFailed||
                 outcome.sourceOutcome->failure.reason==
                     algo::T3AReason::SourceFailed));
    }
}

void sweep(bool publication,long first,long last,long expectedTotal) {
    Signature reference;
    const long total=measure(publication,reference);
    REQUIRE(expectedTotal<=0||total==expectedTotal);
    if (last<=0||last>total) last=total;
    long faults=0;
    for (long ordinal=first;ordinal<=last;++ordinal) {
        RouteFixture fixture(publication);
        mesh::Mesh destination;
        destination.addNode({-41,-43,-47},{0,990});
        const auto before=mesh::meshHash(destination);
        auto outcome=runTracked(fixture,destination,ordinal);
        if (!gFailureInjected)
            std::printf("probe_t4a missed ordinal=%ld observed=%ld\n",
                        ordinal,gAllocationOrdinal);
        REQUIRE(gFailureInjected);
        requireStructuredFailure(outcome,before,destination);
        ++faults;

        RouteFixture retryFixture(publication);
        auto retry=runTracked(retryFixture,destination,0);
        REQUIRE(signatureOf(retry)==reference);
        REQUIRE(retry.retainedResult&&
                retry.retainedResult->metrics.phase3Allocations==0);
        REQUIRE(retry.retainedResult&&
                retry.retainedResult->metrics.publicationTailAllocations==0);
    }
    std::printf("probe_t4a %s range=%ld-%ld total=%ld faults=%ld "
                "phase3=%ld publication=%ld\n",
        publication?"publication":"route",first,last,total,faults,
        gPhase3Allocations,gPublicationTailAllocations);
}

long parseLong(const char* value) {
    long result=0;
    const auto end=value+std::strlen(value);
    const auto parsed=std::from_chars(value,end,result);
    if (parsed.ec!=std::errc{}||parsed.ptr!=end||result<0)
        throw std::invalid_argument("invalid ordinal");
    return result;
}
} // namespace

int main(int argc,char** argv) {
    std::set_terminate(terminateDiagnostic);
    static_assert(std::is_nothrow_move_constructible_v<paacs::algo::T4AOutcome>);
    try {
        gPhase3Allocations=0;
        gPublicationTailAllocations=0;
        if (argc==3&&std::string_view(argv[1])=="--count") {
            const bool publication=std::string_view(argv[2])=="publication";
            REQUIRE(publication||std::string_view(argv[2])=="route");
            Signature signature;
            const auto count=measure(publication,signature);
            std::printf("t4a %s allocation count: %ld\n",
                        publication?"publication":"route",count);
        } else if (argc==4&&
                  std::string_view(argv[1])=="--expect-count") {
            const bool publication=std::string_view(argv[2])=="publication";
            REQUIRE(publication||std::string_view(argv[2])=="route");
            Signature signature;
            const auto count=measure(publication,signature);
            REQUIRE(count==parseLong(argv[3]));
            std::printf("t4a %s allocation count: %ld\n",
                        publication?"publication":"route",count);
        } else if (argc==5&&std::string_view(argv[1])=="--range") {
            const bool publication=std::string_view(argv[2])=="publication";
            REQUIRE(publication||std::string_view(argv[2])=="route");
            sweep(publication,parseLong(argv[3]),parseLong(argv[4]),0);
        } else if (argc==6&&
                  std::string_view(argv[1])=="--expect-range") {
            const bool publication=std::string_view(argv[2])=="publication";
            REQUIRE(publication||std::string_view(argv[2])=="route");
            sweep(publication,parseLong(argv[3]),parseLong(argv[4]),
                  parseLong(argv[5]));
        } else if (argc==1) {
            Signature route;
            const auto routeCount=measure(false,route);
            sweep(false,1,(std::min)(16L,routeCount),routeCount);
        } else {
            std::printf("usage: paacs_probe_t4a [--count route|publication|"
                        "--expect-count route|publication TOTAL|--range "
                        "route|publication FIRST LAST|--expect-range "
                        "route|publication FIRST LAST TOTAL]\n");
            return 2;
        }
    } catch (const std::exception& error) {
        gTrackAllocations=false;
        std::printf("probe_t4a exception: %s\n",error.what());
        ++failures;
    }
    REQUIRE(gPhase3Allocations==0);
    REQUIRE(gPublicationTailAllocations==0);
    if (failures) {
        std::printf("probe_t4a FAILED failures=%d\n",failures);
        return 1;
    }
    std::puts("probe_t4a PASSED");
    return 0;
}
