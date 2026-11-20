// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Increment-5 standalone allocation-failure probe. The bounded builder and
// atomic FaceMeshResult publisher must contain every allocation failure; the
// global mesh remains logically unchanged and a disarmed retry succeeds.
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <new>

#include "paacs/algo/FaceMesherAdapter.hpp"

static long gFailCountdownFace = 0;
static long gCurrentFaceInjection = 0;
static const char* gCurrentFacePhase = "startup";
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t n) {
    if (gFailCountdownFace > 0 && --gFailCountdownFace == 0)
        throw std::bad_alloc();
    if (void* p = std::malloc(n ? n : 1)) return p;
    throw std::bad_alloc();
}
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
namespace algo = paacs::algo;
namespace bnd = paacs::bnd;
namespace mesh = paacs::mesh;
using paacs::Index;
using paacs::Vec2;
using paacs::Vec3;

static int gBad = 0;
#define REQUIRE(c)                                                        \
    do {                                                                  \
        if (!(c)) {                                                       \
            std::printf("FAIL %d: %s\n", __LINE__, #c);                   \
            ++gBad;                                                       \
        }                                                                 \
    } while (0)

struct Provider : bnd::GeometryProvider {
    bnd::GeometryProviderStamp value{700, 3};
    bnd::GeometryProviderStamp stamp() const override { return value; }
    Vec3 pointAt(bnd::ChartId, const Vec2& uv) const override {
        return {uv.u, uv.v, 0};
    }
};

struct Fixture {
    mesh::Mesh global;
    algo::SyntheticFaceInput input;
};

Fixture makeFixture() {
    Fixture f;
    const Vec2 uv[6] = {
        {0, 0}, {2, 0}, {4, 0}, {4, 2}, {2, 2}, {0, 2}};
    f.input.chart = 7;
    f.input.faceClass = {2, 5};
    for (std::size_t i = 0; i < 6; ++i) {
        const Index global =
            f.global.addNode({uv[i].u, uv[i].v, 0},
                             {0, 100 + static_cast<int>(i)});
        algo::SyntheticFaceNode node;
        node.uv = uv[i];
        node.classification = f.global.nodeClass(global);
        node.existingGlobalIndex = global;
        node.constraintVertex = 1000 + i;
        f.input.nodes.push_back(node);
        f.input.outerLoop.push_back(i);
    }
    f.input.triangles = {{{0, 1, 4}},
                         {{0, 4, 5}},
                         {{1, 2, 3}},
                         {{1, 3, 4}}};
    for (std::size_t i = 0; i < 6; ++i) {
        bnd::ConstraintSegment segment;
        segment.id = 10 + i;
        segment.a = 1000 + i;
        segment.b = 1000 + ((i + 1) % 6);
        segment.chart = 7;
        segment.cls = {1, 100 + static_cast<int>(i)};
        f.input.constraints.push_back(segment);
    }
    f.input.constraints.push_back({500, 1001, 1005, 7, {1, 9}});
    f.input.constraints.push_back({501, 1002, 1004, 7, {1, 10}});
    return f;
}

std::uint64_t mix(std::uint64_t h, std::uint64_t v) {
    for (int i = 0; i < 8; ++i) {
        h ^= (v >> (8 * i)) & 0xffU;
        h *= 1099511628211ULL;
    }
    return h;
}
std::uint64_t bits(double d) {
    std::uint64_t value = 0;
    static_assert(sizeof(value) == sizeof(d));
    std::memcpy(&value, &d, sizeof(value));
    return value;
}
std::uint64_t meshDigest(const mesh::Mesh& m) {
    std::uint64_t h = 1469598103934665603ULL;
    h = mix(h, static_cast<std::uint64_t>(m.nodeCount()));
    for (Index i = 0; i < m.nodeCount(); ++i) {
        h = mix(h, bits(m.node(i).x));
        h = mix(h, bits(m.node(i).y));
        h = mix(h, bits(m.node(i).z));
        h = mix(h, static_cast<std::uint64_t>(m.nodeClass(i).dim));
        h = mix(h, static_cast<std::uint64_t>(m.nodeClass(i).tag));
    }
    h = mix(h, static_cast<std::uint64_t>(m.cellCount()));
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto& cell = m.cell(i);
        h = mix(h, static_cast<std::uint64_t>(cell.type));
        h = mix(h, static_cast<std::uint64_t>(cell.onEntity.dim));
        h = mix(h, static_cast<std::uint64_t>(cell.onEntity.tag));
        const Index* conn = m.cellNodes(i);
        for (int j = 0; j < mesh::nodesPerCell(cell.type); ++j)
            h = mix(h, static_cast<std::uint64_t>(conn[j]));
    }
    return h;
}

bnd::SteinerPolicy policy() {
    bnd::SteinerPolicy p;
    p.splitBudget = 8;
    p.minFeatureLevel = 52;
    return p;
}

long builderSweep() {
    long faults = 0;
    long candidateFailures = 0;
    long stagedFailures = 0;
    long cavityFailures = 0;
    long rebuildFailures = 0;
    for (long k = 1;; ++k) {
        gCurrentFacePhase = "builder";
        gCurrentFaceInjection = k;
        Fixture f = makeFixture();
        Provider provider;
        const auto before = meshDigest(f.global);
        bool escaped = false;
        algo::FaceAdapterOutcome built;
        gFailCountdownFace = k;
        try {
            built = algo::buildSyntheticFaceMesh(
                f.input, f.global, provider, provider.value, policy());
        } catch (...) {
            escaped = true;
        }
        gFailCountdownFace = 0;
        REQUIRE(!escaped);
        if (escaped) {
            ++faults;
            continue;
        }
        REQUIRE(meshDigest(f.global) == before);
        if (built.status == algo::FaceAdapterStatus::Ready) {
            REQUIRE(built.result != nullptr);
            break;
        }
        ++faults;
        REQUIRE(built.status == algo::FaceAdapterStatus::Failed);
        REQUIRE(built.reason == algo::FaceAdapterReason::AllocationFailed ||
                built.reason == algo::FaceAdapterReason::RecoveryStopped);
        REQUIRE(built.result == nullptr);
        switch (built.fail) {
            case bnd::SteinerFail::CandidateAllocFailed:
                ++candidateFailures;
                break;
            case bnd::SteinerFail::StagedNodeAllocFailed:
                ++stagedFailures;
                break;
            case bnd::SteinerFail::CavityThrew:
                ++cavityFailures;
                break;
            case bnd::SteinerFail::IndexRebuildFailed:
                ++rebuildFailures;
                break;
            default:
                REQUIRE(false);
                break;
        }

        const auto retry = algo::buildSyntheticFaceMesh(
            f.input, f.global, provider, provider.value, policy());
        REQUIRE(retry.status == algo::FaceAdapterStatus::Ready);
        REQUIRE(retry.result != nullptr);
        REQUIRE(meshDigest(f.global) == before);
    }
    std::printf(
        "face-adapter builder: %ld faults / 0 escapes "
        "(candidate %ld, staged %ld, cavity %ld, rebuild %ld)\n",
        faults, candidateFailures, stagedFailures, cavityFailures,
        rebuildFailures);
    return faults;
}

long publicationSweep() {
    long faults = 0;
    for (long k = 1;; ++k) {
        gCurrentFacePhase = "publication";
        gCurrentFaceInjection = k;
        Fixture f = makeFixture();
        Provider provider;
        auto built = algo::buildSyntheticFaceMesh(
            f.input, f.global, provider, provider.value, policy());
        REQUIRE(built.status == algo::FaceAdapterStatus::Ready);
        REQUIRE(built.result != nullptr);
        if (!built.result) return faults;
        const auto before = meshDigest(f.global);
        bool escaped = false;
        algo::FacePublication publication;
        gFailCountdownFace = k;
        try {
            publication =
                algo::publishFaceMeshResult(f.global, *built.result);
        } catch (...) {
            escaped = true;
        }
        gFailCountdownFace = 0;
        REQUIRE(!escaped);
        if (escaped) {
            ++faults;
            continue;
        }
        if (publication.status == algo::FacePublishStatus::Published)
            break;
        ++faults;
        REQUIRE(publication.status == algo::FacePublishStatus::Failed);
        REQUIRE(publication.reason ==
                algo::FacePublishReason::AllocationFailed);
        REQUIRE(meshDigest(f.global) == before);
        // Capacity residue is permitted; logical state is reusable.
        const auto retry =
            algo::publishFaceMeshResult(f.global, *built.result);
        REQUIRE(retry.status == algo::FacePublishStatus::Published);
    }
    std::printf("face-adapter publication: %ld faults / 0 escapes\n",
                faults);
    return faults;
}

}  // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::set_terminate([] {
        std::printf(
            "TERMINATE in face-adapter %s injection %ld "
            "(active exception %s)\n",
            gCurrentFacePhase, gCurrentFaceInjection,
            std::current_exception() ? "yes" : "no");
        std::_Exit(97);
    });
    const long builderFaults = builderSweep();
    REQUIRE(builderFaults >= 4);
    const long publicationFaults = publicationSweep();
    REQUIRE(publicationFaults >= 4);
    if (gBad) {
        std::printf("FACE ADAPTER PROBE FAILED (%d)\n", gBad);
        return 1;
    }
    std::printf("FACE ADAPTER PROBE PASSED\n");
    return 0;
}
