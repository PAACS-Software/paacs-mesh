// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// CK-3 candidate/handoff allocation-failure containment. This standalone
// probe replaces global operator new and therefore is excluded from sanitizer
// builds. Every injected failure must be structured, leave source values
// byte/logically unchanged, and permit a clean same-fixture retry.
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>
#include <vector>

#include "paacs/algo/TetCandidate.hpp"
#include "../src/algo/TetCandidateDetail.hpp"

namespace {
bool gArmed = false;
std::size_t gCalls = 0;
std::size_t gFailAt = 0;
int gBad = 0;
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t size) {
    if (gArmed) {
        ++gCalls;
        if (gFailAt != 0 && gCalls == gFailAt) throw std::bad_alloc{};
    }
    if (void* memory = std::malloc(size == 0 ? 1 : size)) return memory;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t size) { return ::operator new(size); }
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

#define REQUIRE(condition)                                                   \
    do {                                                                     \
        if (!(condition)) {                                                  \
            std::printf("FAIL %d: %s\n", __LINE__, #condition);             \
            ++gBad;                                                          \
        }                                                                    \
    } while (0)

std::vector<Vec3> fixture() {
    return {{0, 0, 0}, {3, 0, 0}, {0, 4, 0}, {0, 0, 5},
            {1, 1, 1}, {2, 1, 1}, {1, 2, 1}, {1, 1, 3}, {2, 2, 2}};
}

std::uint64_t digestPoints(const std::vector<Vec3>& points) {
    std::uint64_t hash = 1469598103934665603ULL;
    auto mix = [&](std::uint64_t value) {
        hash = (hash ^ value) * 1099511628211ULL;
    };
    mix(points.size());
    for (const auto& point : points) {
        std::uint64_t x = 0, y = 0, z = 0;
        std::memcpy(&x, &point.x, sizeof(x));
        std::memcpy(&y, &point.y, sizeof(y));
        std::memcpy(&z, &point.z, sizeof(z));
        mix(x); mix(y); mix(z);
    }
    return hash;
}

std::uint64_t digestCandidate(const algo::TetCandidate& candidate) {
    std::uint64_t hash = 1469598103934665603ULL;
    auto mix = [&](std::uint64_t value) {
        hash = (hash ^ value) * 1099511628211ULL;
    };
    mix(candidate.points.size());
    mix(candidate.cells.size());
    mix(candidate.boundaryFaces.size());
    for (const auto& cell : candidate.cells) {
        for (const Index value : cell.vertices)
            mix(static_cast<std::uint64_t>(value));
        for (const Index value : cell.neighbours)
            mix(static_cast<std::uint64_t>(value));
        mix(static_cast<std::uint64_t>(cell.sourceOrdinal));
    }
    return hash;
}

std::uint64_t digestMesh(const mesh::Mesh& value) {
    std::uint64_t hash = 1469598103934665603ULL;
    auto mix = [&](std::uint64_t word) {
        hash = (hash ^ word) * 1099511628211ULL;
    };
    mix(static_cast<std::uint64_t>(value.nodeCount()));
    mix(static_cast<std::uint64_t>(value.cellCount()));
    for (Index i = 0; i < value.nodeCount(); ++i) {
        const auto& point = value.node(i);
        std::uint64_t x = 0, y = 0, z = 0;
        std::memcpy(&x, &point.x, sizeof(x));
        std::memcpy(&y, &point.y, sizeof(y));
        std::memcpy(&z, &point.z, sizeof(z));
        mix(x); mix(y); mix(z);
        const auto cls = value.nodeClass(i);
        mix(static_cast<std::uint64_t>(cls.dim));
        mix(static_cast<std::uint64_t>(cls.tag));
    }
    return hash;
}

}  // namespace

int main() {
    const auto points = fixture();
    const auto pointDigest = digestPoints(points);

    // Warm predicate/thread-local/library state before defining the exact
    // replayable allocation ordinal domain.
    auto warmBuild = algo::buildTetCandidate(points);
    REQUIRE(warmBuild.built());

    gCalls = 0;
    gFailAt = 0;
    gArmed = true;
    auto cleanBuild = algo::buildTetCandidate(points);
    gArmed = false;
    const std::size_t buildCalls = gCalls;
    std::size_t buildAbsorbed = 0;
    REQUIRE(cleanBuild.built());
    for (std::size_t ordinal = 1; ordinal <= buildCalls; ++ordinal) {
        gCalls = 0;
        gFailAt = ordinal;
        gArmed = true;
        auto failed = algo::buildTetCandidate(points);
        gArmed = false;
        // libstdc++'s stable_sort is permitted to absorb a temporary-buffer
        // allocation failure and complete through its allocation-free
        // fallback. That is a successful construction, not an escaped or
        // partially handled failure; pin its complete observable result.
        if (failed.built()) {
            ++buildAbsorbed;
            REQUIRE(cleanBuild.built());
            REQUIRE(digestCandidate(*failed.candidate) ==
                    digestCandidate(*cleanBuild.candidate));
            REQUIRE(algo::validateTetCandidate(*failed.candidate).valid);
        } else {
            REQUIRE(failed.reason ==
                    algo::TetCandidateReason::AllocationFailed);
        }
        REQUIRE(digestPoints(points) == pointDigest);
        auto retry = algo::buildTetCandidate(points);
        REQUIRE(retry.built());
    }

    REQUIRE(cleanBuild.built());
    if (!cleanBuild.built()) return 1;
    const auto& candidate = *cleanBuild.candidate;
    const auto candidateDigest = digestCandidate(candidate);
    std::vector<Index> selected(candidate.cells.size());
    for (std::size_t i = 0; i < selected.size(); ++i)
        selected[i] = static_cast<Index>(i);
    std::vector<Index> globals(candidate.points.size());
    std::vector<mesh::Classification> classes(
        candidate.points.size(), mesh::Classification{3, 8});
    mesh::Mesh global;
    for (std::size_t i = 0; i < globals.size(); ++i)
        globals[i] = global.addNode(candidate.points[i], classes[i]);
    const auto globalDigest = digestMesh(global);

    gCalls = 0;
    gFailAt = 0;
    gArmed = true;
    auto cleanImport = algo::importTetCandidate(
        candidate, selected, globals, classes, global, {3, 8});
    gArmed = false;
    const std::size_t importCalls = gCalls;
    REQUIRE(cleanImport.imported());
    REQUIRE(cleanImport.result->metrics.importCandidateValidationCalls==1);
    REQUIRE(cleanImport.result->metrics.importGlobalMapRows==
            candidate.points.size());
    REQUIRE(cleanImport.result->metrics.importGlobalMapComparisons>0);
    REQUIRE(cleanImport.result->metrics.importGlobalMapInsertions==
            candidate.points.size());
    REQUIRE(cleanImport.result->metrics.importLocalGlobalSlotsInitialized==
            static_cast<std::uint64_t>(global.nodeCount()));
    REQUIRE(cleanImport.result->metrics.authoritativeModeledStorageBytes>0);
    REQUIRE(cleanImport.result->metrics.importScratchModeledPayloadBytes>0);
    for (std::size_t ordinal = 1; ordinal <= importCalls; ++ordinal) {
        gCalls = 0;
        gFailAt = ordinal;
        gArmed = true;
        auto failed = algo::importTetCandidate(
            candidate, selected, globals, classes, global, {3, 8});
        gArmed = false;
        REQUIRE(!failed.imported());
        REQUIRE(failed.reason == algo::TetCandidateReason::AllocationFailed);
        REQUIRE(digestCandidate(candidate) == candidateDigest);
        REQUIRE(digestMesh(global) == globalDigest);
        auto retry = algo::importTetCandidate(
            candidate, selected, globals, classes, global, {3, 8});
        REQUIRE(retry.imported());
    }

    static constexpr std::array<std::array<std::size_t,3>,4> tetFaces{{
        {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
    std::array<Index,3> interiorGlobal{};
    bool foundInterior=false;
    for (const auto& cell:candidate.cells) {
        for (std::size_t localFace=0;localFace<4;++localFace) {
            if (cell.neighbours[localFace]==kInvalidIndex) continue;
            for (std::size_t i=0;i<3;++i)
                interiorGlobal[i]=globals[static_cast<std::size_t>(
                    cell.vertices[tetFaces[localFace][i]])];
            foundInterior=true;
            break;
        }
        if (foundInterior) break;
    }
    REQUIRE(foundInterior);
    REQUIRE(!candidate.boundaryFaces.empty());
    std::array<Index,3> boundaryGlobal{};
    if (!candidate.boundaryFaces.empty())
        for (std::size_t i=0;i<3;++i)
            boundaryGlobal[i]=globals[static_cast<std::size_t>(
                candidate.boundaryFaces[0].vertices[i])];
    std::array<Index,2> missingGlobal{kInvalidIndex,kInvalidIndex};
    for (std::size_t a=0;a<candidate.points.size()&&
                         missingGlobal[0]==kInvalidIndex;++a)
        for (std::size_t b=a+1;b<candidate.points.size();++b) {
            bool present=false;
            for (const auto& cell:candidate.cells) {
                bool hasA=false,hasB=false;
                for (const Index point:cell.vertices) {
                    hasA=hasA||point==static_cast<Index>(a);
                    hasB=hasB||point==static_cast<Index>(b);
                }
                if (hasA&&hasB) { present=true; break; }
            }
            if (!present) {
                missingGlobal={globals[a],globals[b]};
                break;
            }
        }
    REQUIRE(missingGlobal[0]!=kInvalidIndex);
    const std::array<algo::TetImportProtectedEdge,2> partialEdges{{
        {{interiorGlobal[1],interiorGlobal[0]},{2,701}},
        {{missingGlobal[1],missingGlobal[0]},{1,704}}}};
    const std::array<algo::TetImportProtectedFace,2> partialFaces{{
        {{interiorGlobal[2],interiorGlobal[0],interiorGlobal[1]},{2,702}},
        {{boundaryGlobal[1],boundaryGlobal[2],boundaryGlobal[0]},{2,703}}}};

    gCalls=0;
    gFailAt=0;
    gArmed=true;
    auto cleanPartial=
        algo::tet_candidate_detail::
            importTetCandidateWithPartialProtectedBoundary(
            candidate,selected,globals,classes,global,{3,8},
            partialEdges,partialFaces);
    gArmed=false;
    const std::size_t partialCalls=gCalls;
    REQUIRE(cleanPartial.imported());
    REQUIRE(cleanPartial.result->metrics.partialImportFaceDataCopies==2);
    REQUIRE(cleanPartial.result->metrics.
                partialImportFaceIncidenceRowsCopied==3);
    REQUIRE(cleanPartial.result->metrics.partialImportEdgeDataReads==1);
    REQUIRE(cleanPartial.result->metrics.partialImportSortComparisons>0);
    for (std::size_t ordinal=1;ordinal<=partialCalls;++ordinal) {
        gCalls=0;
        gFailAt=ordinal;
        gArmed=true;
        auto failed=algo::tet_candidate_detail::
            importTetCandidateWithPartialProtectedBoundary(
            candidate,selected,globals,classes,global,{3,8},
            partialEdges,partialFaces);
        gArmed=false;
        REQUIRE(!failed.imported());
        REQUIRE(failed.reason==algo::TetCandidateReason::AllocationFailed);
        REQUIRE(digestCandidate(candidate)==candidateDigest);
        REQUIRE(digestMesh(global)==globalDigest);
        auto retry=algo::tet_candidate_detail::
            importTetCandidateWithPartialProtectedBoundary(
            candidate,selected,globals,classes,global,{3,8},
            partialEdges,partialFaces);
        REQUIRE(retry.imported());
    }

    std::printf("probe_tet_handoff: build=%zu (%zu absorbed) import=%zu "
                "partial=%zu zero escapes\n",
                buildCalls, buildAbsorbed, importCalls,partialCalls);
    return gBad == 0 ? 0 : 1;
}
