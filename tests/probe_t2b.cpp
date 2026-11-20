// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// T2-B partitionable end-to-end allocation probe. Every injected allocation
// must be contained by the structured v3 outcome, preserve the exact caller
// mesh, retain report-tier semantics, and permit a disarmed same-fixture
// retry. The global hooks also prove zero Tet Phase-3 and publication-tail
// allocations.
#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <map>
#include <memory>
#include <new>
#include <string_view>

#if defined(_WIN32)
#define NOMINMAX
#include <Windows.h>
#include <intrin.h>
#else
#include <unistd.h>
#endif

#include "paacs/algo/BoundaryConformingMesher.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/sizing/SizeField.hpp"
#include "paacs/topo/TransactionAccess.hpp"

static std::size_t gCalls = 0;
static std::size_t gFailAt = 0;
static std::size_t gPhase3Allocations = 0;
static std::size_t gPublicationAllocations = 0;
static std::uint8_t gFailurePhase = 0;
static bool gFailSessionIndex = false;
static bool gSessionIndexFaultTriggered = false;

namespace {
void terminateDiagnostic() noexcept {
    char buffer[256];
    char* cursor = buffer;
    auto literal = [&](const char* value) noexcept {
        const std::size_t count = std::strlen(value);
        std::memcpy(cursor,value,count);
        cursor += count;
    };
    auto number = [&](std::size_t value) noexcept {
        const auto converted = std::to_chars(cursor,buffer+sizeof(buffer),value);
        if (converted.ec == std::errc{}) cursor = converted.ptr;
    };
    literal("probe_t2b terminate calls="); number(gCalls);
    literal(" fail="); number(gFailAt);
    literal(" phase="); number(gFailurePhase);
    literal(" phase3="); number(gPhase3Allocations);
    literal(" publication="); number(gPublicationAllocations);
    literal(" session_index=");
    number(paacs::bnd::detail::recoverySessionIndexBuildActive() ? 1U : 0U);
    literal("\n");
#if defined(_WIN32)
    DWORD written = 0;
    (void)::WriteFile(::GetStdHandle(STD_ERROR_HANDLE),buffer,
                      static_cast<DWORD>(cursor-buffer),&written,nullptr);
    __fastfail(7);
#else
    (void)::write(STDERR_FILENO,buffer,
                  static_cast<std::size_t>(cursor-buffer));
    std::_Exit(134);
#endif
}
}  // namespace

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t bytes) {
    ++gCalls;
    if (paacs::topo::detail::phase3Active()) ++gPhase3Allocations;
    if (paacs::algo::boundary_conforming_detail::publicationTailActive())
        ++gPublicationAllocations;
    if (gFailSessionIndex && !gSessionIndexFaultTriggered &&
        paacs::bnd::detail::recoverySessionIndexBuildActive()) {
        gSessionIndexFaultTriggered = true;
        gFailurePhase =
            paacs::algo::boundary_conforming_detail::
                testingAllocationPhaseMarker();
        throw std::bad_alloc{};
    }
    if (gFailAt != 0 && gCalls == gFailAt) {
        gFailurePhase =
            paacs::algo::boundary_conforming_detail::
                testingAllocationPhaseMarker();
        throw std::bad_alloc{};
    }
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
            std::printf("FAIL %d: %s\n",__LINE__,#condition);              \
            ++failures;                                                      \
        }                                                                    \
    } while (0)

class Provider final : public algo::VolumeBoundaryOwnershipProvider {
public:
    algo::VolumeBoundaryProviderStamp current{0x54324250524f5649ULL,1};
    std::map<algo::FacePersistentId,std::vector<bnd::VolumeRegionUse>> rows;
    algo::ProviderRecoveryEnvelope envelope;

    algo::VolumeBoundaryProviderStamp stamp() const noexcept override {
        return current;
    }
    algo::VolumeFaceOwnership ownership(
        algo::FacePersistentId face) const override {
        algo::VolumeFaceOwnership result;
        result.providerStamp = current;
        result.facePersistentId = face;
        const auto found = rows.find(face);
        if (found != rows.end()) result.owners = found->second;
        return result;
    }
    algo::ProviderRecoveryEnvelope recoveryEnvelope(
        std::span<const Vec3>) const override {
        auto result = envelope;
        result.providerStamp = current;
        return result;
    }
};

algo::ProviderRecoveryEnvelope boxEnvelope(double lo,double hi) {
    algo::ProviderRecoveryEnvelope result;
    result.availability = algo::VolumeProviderAvailability::Available;
    for (std::size_t i = 0;i < 8;++i)
        result.guardCoordinates.push_back(
            {(i&1U)?hi:lo,(i&2U)?hi:lo,(i&4U)?hi:lo});
    result.outwardHullFaces = {
        {0,4,2},{4,6,2},{1,3,5},{3,7,5},
        {0,1,4},{1,5,4},{2,6,3},{6,7,3},
        {0,2,1},{2,3,1},{4,5,6},{5,7,6}};
    return result;
}

Provider providerFor(const geo::Model& model) {
    Provider provider;
    const auto solid = model.solids().front();
    Vec3 lo{},hi{};
    solid->bbox(lo,hi);
    const Vec3 center = (lo+hi)*0.5;
    const double span = std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z});
    provider.envelope = boxEnvelope(
        std::min({lo.x,lo.y,lo.z})-span,
        std::max({hi.x,hi.y,hi.z})+span);
    for (const geo::GFace* face : model.faces()) {
        Vec2 uvLo{},uvHi{};
        face->paramBounds(uvLo,uvHi);
        const Vec2 uv = (uvLo+uvHi)*0.5;
        const int sense = face->normalAt(uv).dot(
            face->pointAt(uv)-center) > 0.0 ? 1 : -1;
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1] =
            {{0,sense}};
    }
    return provider;
}

std::uint64_t mix(std::uint64_t h,std::uint64_t value) noexcept {
    for (int i = 0;i < 8;++i) {
        h ^= (value>>(8*i))&0xffU;
        h *= 1099511628211ULL;
    }
    return h;
}

std::uint64_t bits(double value) noexcept {
    std::uint64_t result = 0;
    static_assert(sizeof(result) == sizeof(value));
    std::memcpy(&result,&value,sizeof(result));
    return result;
}

std::uint64_t digest(const mesh::Mesh& mesh) noexcept {
    std::uint64_t h = 1469598103934665603ULL;
    h = mix(h,static_cast<std::uint64_t>(mesh.nodeCount()));
    h = mix(h,static_cast<std::uint64_t>(mesh.cellCount()));
    for (Index i = 0;i < mesh.nodeCount();++i) {
        const auto& p = mesh.node(i);
        const auto& c = mesh.nodeClass(i);
        h = mix(mix(mix(h,bits(p.x)),bits(p.y)),bits(p.z));
        h = mix(mix(h,static_cast<std::uint32_t>(c.dim)),
                static_cast<std::uint32_t>(c.tag));
    }
    for (Index i = 0;i < mesh.cellCount();++i) {
        const auto& cell = mesh.cell(i);
        h = mix(h,static_cast<std::uint8_t>(cell.type));
        h = mix(mix(h,static_cast<std::uint32_t>(cell.onEntity.dim)),
                static_cast<std::uint32_t>(cell.onEntity.tag));
        const Index* nodes = mesh.cellNodes(i);
        for (int n = 0;n < mesh::nodesPerCell(cell.type);++n)
            h = mix(h,static_cast<std::uint64_t>(nodes[n]));
    }
    return h;
}

mesh::Mesh sentinelMesh() {
    mesh::Mesh result;
    result.addNode({-7,-8,-9},{0,991});
    return result;
}

struct Fixture {
    std::unique_ptr<geo::Kernel> kernel;
    std::unique_ptr<geo::Model> model;
    sizing::MeshSpec spec;
    algo::BoundaryConformingV3Request request;
    Provider provider;
};

Fixture fixture() {
    Fixture result;
    result.kernel = geo::makeBuiltinKernel();
    result.model = result.kernel->makeBox({0,0,0},{1,1,1});
    // Small but complete vertical slice: all six constrained faces, every
    // protected artifact edge, Tet construction/import, set recovery,
    // topological extraction, exact audit and atomic publication.  Keeping
    // only the authoritative box vertices makes an exhaustive injected-fault
    // campaign practical without weakening any pipeline phase.
    result.spec.field = std::make_shared<sizing::ConstantField>(10.0);
    result.spec.threads = 1;
    result.request.expectedSurfaceSetStamp = {0x5432424152544946ULL,1};
    result.provider = providerFor(*result.model);
    result.request.expectedOwnershipProviderStamp = result.provider.current;
    result.request.envelopeMode =
        algo::RecoveryEnvelopeMode::DerivedFiniteBox;
    return result;
}

algo::BoundaryConformingV3Outcome run(
    Fixture& f,mesh::Mesh& destination,std::size_t failAt) {
    gCalls = 0;
    gFailAt = failAt;
    auto outcome = algo::meshModelBoundaryConforming(
        *f.model,f.spec,f.request,destination,f.provider);
    gFailAt = 0;
    return outcome;
}

std::size_t cleanAllocationCount(std::uint64_t& expectedDigest) {
    {
        auto warm = fixture();
        mesh::Mesh warmDestination = sentinelMesh();
        auto warmOutcome = run(warm,warmDestination,0);
        REQUIRE(warmOutcome.published());
    }
    auto f = fixture();
    mesh::Mesh destination = sentinelMesh();
    auto outcome = run(f,destination,0);
    REQUIRE(outcome.published());
    REQUIRE(outcome.retainedResult != nullptr);
    REQUIRE(outcome.retainedResult && outcome.retainedResult->boundaryAudit.passed);
    REQUIRE(outcome.report != nullptr);
    bool sawEdge = false, sawFacet = false;
    std::size_t recoveredEdges = 0,recoveredFacets = 0;
    if (outcome.report)
        for (const auto& row : outcome.report->processed) {
            sawEdge = sawEdge || row.key.kind == algo::VolumeWorkKind::ArtifactEdge;
            sawFacet = sawFacet || row.key.kind == algo::VolumeWorkKind::Facet;
            recoveredEdges += row.key.kind == algo::VolumeWorkKind::ArtifactEdge &&
                              row.status == algo::BoundaryConformingV3WorkStatus::Recovered;
            recoveredFacets += row.key.kind == algo::VolumeWorkKind::Facet &&
                               row.status == algo::BoundaryConformingV3WorkStatus::Recovered;
        }
    REQUIRE(sawEdge);
    REQUIRE(sawFacet);
    std::printf("probe_t2b route evidence: recovered_edges=%zu "
                "recovered_facets=%zu\n",recoveredEdges,recoveredFacets);
    expectedDigest = digest(destination);
    return gCalls;
}

void sessionIndexAllocationRegression() {
    auto f = fixture();
    mesh::Mesh destination = sentinelMesh();
    const auto before = digest(destination);
    gCalls = 0;
    gFailAt = 0;
    gSessionIndexFaultTriggered = false;
    gFailSessionIndex = true;
    auto failed = algo::meshModelBoundaryConforming(
        *f.model,f.spec,f.request,destination,f.provider);
    gFailSessionIndex = false;
    REQUIRE(gSessionIndexFaultTriggered);
    REQUIRE(!failed.published());
    REQUIRE(failed.status == algo::BoundaryConformingV3Status::Failed);
    REQUIRE(failed.failure.phase ==
            algo::BoundaryConformingV3Phase::SurfaceConstruction);
    REQUIRE(failed.surfaceOutcome != nullptr);
    bool nestedAllocation = failed.surfaceOutcome &&
        failed.surfaceOutcome->setFailure.reason ==
            algo::SurfaceReason::AllocationFailed;
    if (failed.surfaceOutcome)
        for (const auto& row : failed.surfaceOutcome->processed)
            nestedAllocation = nestedAllocation ||
                row.failure.reason == algo::SurfaceReason::AllocationFailed;
    REQUIRE(nestedAllocation);
    REQUIRE(failed.destinationUnchanged);
    REQUIRE(digest(destination) == before);
    auto retry = algo::meshModelBoundaryConforming(
        *f.model,f.spec,f.request,destination,f.provider);
    REQUIRE(retry.published());
}

void sweep(std::size_t begin,std::size_t end,std::uint64_t expectedDigest) {
    for (std::size_t fault = begin;fault <= end;++fault) {
        auto f = fixture();
        mesh::Mesh destination = sentinelMesh();
        const auto before = digest(destination);
        auto failed = run(f,destination,fault);
        REQUIRE(!failed.published());
        REQUIRE(failed.status == algo::BoundaryConformingV3Status::Failed);
        bool nestedAllocation = failed.surfaceOutcome &&
            failed.surfaceOutcome->setFailure.reason ==
                algo::SurfaceReason::AllocationFailed;
        if (failed.surfaceOutcome)
            for (const auto& row : failed.surfaceOutcome->processed) {
                nestedAllocation = nestedAllocation ||
                    row.failure.reason == algo::SurfaceReason::AllocationFailed;
                nestedAllocation = nestedAllocation ||
                    (row.failure.reason == algo::SurfaceReason::RecoveryStopped &&
                     (row.failure.bndFail == bnd::SteinerFail::CandidateAllocFailed ||
                      row.failure.bndFail == bnd::SteinerFail::StagedNodeAllocFailed ||
                      row.failure.bndFail == bnd::SteinerFail::CavityThrew ||
                      row.failure.bndFail == bnd::SteinerFail::IndexRebuildFailed));
            }
        const bool allocationAttributed =
            failed.failure.reason ==
                algo::BoundaryConformingV3Reason::AllocationFailed ||
            (failed.failure.reason ==
                 algo::BoundaryConformingV3Reason::SurfaceConstructionFailed &&
             nestedAllocation);
        REQUIRE(allocationAttributed);
        if (!allocationAttributed)
            std::printf("fault=%zu status=%d reason=%d phase=%d calls=%zu\n",
                fault,static_cast<int>(failed.status),
                static_cast<int>(failed.failure.reason),
                static_cast<int>(failed.failure.phase),gCalls);
        if (!allocationAttributed && failed.surfaceOutcome) {
            std::printf("  surface status=%d setReason=%d rows=%zu "
                        "unattempted=%zu complete=%d\n",
                static_cast<int>(failed.surfaceOutcome->setStatus),
                static_cast<int>(failed.surfaceOutcome->setFailure.reason),
                failed.surfaceOutcome->processed.size(),
                failed.surfaceOutcome->unattempted.size(),
                failed.surfaceOutcome->reportComplete ? 1 : 0);
            for (const auto& row : failed.surfaceOutcome->processed)
                if (row.status == algo::SurfaceFaceStatus::Failed)
                    std::printf("    face=%llu reason=%d bndReject=%d "
                                "bndRefuse=%d bndFail=%d\n",
                        static_cast<unsigned long long>(row.facePersistentId),
                        static_cast<int>(row.failure.reason),
                        static_cast<int>(row.failure.bndReject),
                        static_cast<int>(row.failure.bndRefuse),
                        static_cast<int>(row.failure.bndFail));
        }
        REQUIRE(failed.destinationUnchanged);
        REQUIRE(digest(destination) == before);
        REQUIRE(failed.report == nullptr ||
                failed.report->reportComplete == failed.reportComplete);
        auto retry = run(f,destination,0);
        REQUIRE(retry.published());
        REQUIRE(digest(destination) == expectedDigest);
        if (failures != 0) {
            std::printf("first failing allocation ordinal: %zu\n",fault);
            return;
        }
    }
}

std::size_t parse(const char* text) {
    return static_cast<std::size_t>(std::strtoull(text,nullptr,10));
}

}  // namespace

int main(int argc,char** argv) {
    std::set_terminate(&terminateDiagnostic);
    sessionIndexAllocationRegression();
    std::uint64_t expectedDigest = 0;
    const std::size_t count = cleanAllocationCount(expectedDigest);
    REQUIRE(gPhase3Allocations == 0);
    REQUIRE(gPublicationAllocations == 0);
    if (argc == 3 && std::string_view(argv[1]) == "--expect-count") {
        REQUIRE(count == parse(argv[2]));
    } else if (argc == 4 && std::string_view(argv[1]) == "--range") {
        const std::size_t begin = parse(argv[2]);
        const std::size_t end = parse(argv[3]);
        REQUIRE(begin >= 1 && begin <= end && end <= count);
        if (failures == 0) sweep(begin,end,expectedDigest);
    } else if (argc == 1) {
        std::printf("probe_t2b allocation count: %zu\n",count);
    } else {
        std::printf("usage: probe_t2b [--expect-count N|--range A B]\n");
        return 2;
    }
    if (failures != 0) return 1;
    std::printf("probe_t2b PASSED: count=%zu phase3=0 publication=0\n",count);
    return 0;
}
