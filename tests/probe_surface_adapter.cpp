// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Increment-6 exhaustive set construction/publication allocation probe.
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <exception>
#include <new>
#include <set>
#include <type_traits>

#include "paacs/algo/SurfaceMesherAdapter.hpp"

static long gSurfaceFailCountdown = 0;
static long gSurfaceInjection = 0;
static const char* gSurfacePhase = "startup";
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t n) {
    if (gSurfaceFailCountdown > 0 && --gSurfaceFailCountdown == 0)
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

int gBad = 0;
#define REQUIRE(c)                                                        \
    do {                                                                  \
        if (!(c)) {                                                       \
            std::printf("FAIL %d: %s\n", __LINE__, #c);                  \
            ++gBad;                                                       \
        }                                                                 \
    } while (0)

struct Geometry : algo::SurfaceGeometryProvider {
    bnd::GeometryProviderStamp stamp() const override { return {7, 1}; }
    Vec3 pointAt(bnd::ChartId, const Vec2& uv) const override {
        return {uv.u, uv.v, 0};
    }
    Vec3 normalAt(bnd::ChartId, const Vec2&) const override {
        return {0, 0, 1};
    }
    algo::ChartCertificate chartCertificate(bnd::ChartId chart) const override {
        return {{7, 1}, chart, {99, 1},
                {{-2, -2}, {4, -2}, {4, 4}, {-2, 4}}, 1, true, false, false};
    }
};

struct Candidate : algo::SurfaceCandidateProvider {
    algo::CandidateProviderStamp stamp() const override { return {8, 1}; }
    algo::SurfaceCandidate build(const algo::SurfaceFaceRequest& request) const override {
        algo::SurfaceCandidate c;
        c.providerStamp = {8, 1}; c.faceStamp = request.faceStamp;
        c.chart = request.chart;
        const Vec2 uv[6]{{0, 0}, {2, 0}, {0, 2},
                         {-1, -1}, {4, -1}, {-1, 4}};
        for (std::size_t i = 0; i < 6; ++i) {
            algo::SurfaceCandidateNode n;
            n.node = {{uv[i].u, uv[i].v, 0},
                      i < 3 ? request.vertices[i].classification
                            : request.faceClass};
            n.parameter = {uv[i], request.chart}; n.frameOnly = i >= 3;
            if (i < 3) c.declaredVertexMap.emplace(
                request.vertices[i].faceLocalVertexId, i);
            c.nodes.push_back(n);
        }
        const std::array<std::array<std::size_t, 3>, 7> tris{{
            {{0, 1, 2}}, {{3, 4, 1}}, {{3, 1, 0}}, {{4, 5, 2}},
            {{4, 2, 1}}, {{5, 3, 0}}, {{5, 0, 2}}}};
        for (const auto& tri : tris) c.cells.push_back({tri, request.faceClass});
        c.frameBoundary = {3, 4, 5}; c.exteriorSeedCell = 1;
        return c;
    }
};

struct Fixture {
    mesh::Mesh global;
    algo::SurfaceSetRequest request;
};

Fixture fixture() {
    Fixture f;
    f.request.setStamp = {9000, 1};
    algo::SurfaceFaceRequest face;
    face.faceStamp = {10, 1}; face.chart = 3; face.faceClass = {2, 30};
    face.expectedGeometryStamp = {7, 1};
    face.expectedChartCertificateStamp = {99, 1};
    face.expectedCandidateProviderStamp = {8, 1};
    face.outer.role = algo::SurfaceChainRole::OuterBoundary;
    algo::SurfaceChain chain;
    chain.chainPersistentId = 50;
    chain.role = algo::SurfaceChainRole::OuterBoundary;
    chain.classification = {1, 77}; chain.closed = true;
    chain.owningCurve = algo::ArtifactStamp{500, 1};
    const Vec2 uv[3]{{0, 0}, {2, 0}, {0, 2}};
    for (std::uint64_t i = 0; i < 3; ++i) {
        const Index global = f.global.addNode(
            {uv[i].u, uv[i].v, 0}, {0, static_cast<int>(100 + i)});
        const algo::SharedCurveNodeKey key{500, 1, i};
        face.vertices.push_back({1000 + i, {uv[i], 3},
                                 f.global.nodeClass(global)});
        chain.uses.push_back({1000 + i, key});
        chain.segmentEntityIds.push_back(2000 + i);
        f.request.sharedBindings.push_back(
            {key, global, f.global.nodeClass(global)});
    }
    face.outer.pieces.push_back(std::move(chain));
    f.request.faces.push_back(std::move(face));
    return f;
}

std::uint64_t digest(const mesh::Mesh& m) {
    std::uint64_t h = 1469598103934665603ull;
    auto mix = [&](std::uint64_t x) { h = (h ^ x) * 1099511628211ull; };
    mix(static_cast<std::uint64_t>(m.nodeCount()));
    mix(static_cast<std::uint64_t>(m.cellCount()));
    for (Index i = 0; i < m.nodeCount(); ++i) {
        const auto& p = m.node(i); const auto& c = m.nodeClass(i);
        std::uint64_t bx = 0, by = 0, bz = 0;
        std::memcpy(&bx, &p.x, sizeof(double));
        std::memcpy(&by, &p.y, sizeof(double));
        std::memcpy(&bz, &p.z, sizeof(double));
        mix(bx); mix(by); mix(bz); mix(static_cast<std::uint64_t>(c.dim));
        mix(static_cast<std::uint64_t>(c.tag));
    }
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto& c = m.cell(i); const auto* n = m.cellNodes(i);
        mix(static_cast<std::uint64_t>(c.type));
        mix(static_cast<std::uint64_t>(c.onEntity.dim));
        mix(static_cast<std::uint64_t>(c.onEntity.tag));
        for (int j = 0; j < mesh::nodesPerCell(c.type); ++j)
            mix(static_cast<std::uint64_t>(n[j]));
    }
    return h;
}

bool complete(const algo::SurfaceSetOutcome& out) {
    if (!out.reportComplete) return out.processed.empty() && out.unattempted.empty();
    std::set<algo::FacePersistentId> p, u;
    for (const auto& row : out.processed)
        if (!p.insert(row.facePersistentId).second) return false;
    for (const auto id : out.unattempted)
        if (!u.insert(id).second || p.count(id)) return false;
    return p.size() + u.size() == 1;
}

long constructionSweep() {
    long faults = 0, incomplete = 0;
    for (long k = 1;; ++k) {
        gSurfacePhase = "construction"; gSurfaceInjection = k;
        Fixture f = fixture(); Geometry geometry; Candidate candidate;
        const auto before = digest(f.global);
        bool escaped = false;
        algo::SurfaceSetOutcome out;
        gSurfaceFailCountdown = k;
        try {
            out = algo::buildAndPublishSurfaceSet(
                f.request, f.global, geometry, candidate);
        } catch (...) { escaped = true; }
        gSurfaceFailCountdown = 0;
        REQUIRE(!escaped);
        if (escaped) { ++faults; continue; }
        REQUIRE(complete(out));
        if (out.setStatus == algo::SurfaceSetStatus::Published) break;
        ++faults;
        if (!out.reportComplete) ++incomplete;
        REQUIRE(digest(f.global) == before);
        if (out.retainedArtifact) {
            bool published = false;
            const auto retry = algo::publishSurfaceSet(
                f.global, *out.retainedArtifact, geometry, published);
            REQUIRE(retry.reason == algo::SurfaceReason::None);
            REQUIRE(published);
        } else {
            const auto retry = algo::buildAndPublishSurfaceSet(
                f.request, f.global, geometry, candidate);
            REQUIRE(retry.setStatus == algo::SurfaceSetStatus::Published);
        }
    }
    REQUIRE(incomplete == 2);
    std::printf("surface construction: %ld faults / 0 escapes "
                "(%ld incomplete envelope)\n", faults, incomplete);
    return faults;
}

long publicationSweep() {
    long faults = 0;
    for (long k = 1;; ++k) {
        Fixture builtFixture = fixture(); Geometry geometry; Candidate candidate;
        auto built = algo::buildAndPublishSurfaceSet(
            builtFixture.request, builtFixture.global, geometry, candidate);
        REQUIRE(built.setStatus == algo::SurfaceSetStatus::Published);
        REQUIRE(built.retainedArtifact != nullptr);
        Fixture target = fixture();
        const auto before = digest(target.global);
        bool escaped = false, published = false;
        algo::SurfaceFailure result;
        gSurfacePhase = "publication"; gSurfaceInjection = k;
        gSurfaceFailCountdown = k;
        try {
            result = algo::publishSurfaceSet(
                target.global, *built.retainedArtifact, geometry, published);
        } catch (...) { escaped = true; }
        gSurfaceFailCountdown = 0;
        REQUIRE(!escaped);
        if (escaped) { ++faults; continue; }
        if (published) break;
        ++faults;
        REQUIRE(result.reason == algo::SurfaceReason::AllocationFailed);
        REQUIRE(result.allocation ==
                    algo::SurfaceAllocationPhase::PublicationPreparation ||
                result.allocation ==
                    algo::SurfaceAllocationPhase::GeometryProvider);
        REQUIRE(digest(target.global) == before);
        bool retried = false;
        const auto retry = algo::publishSurfaceSet(
            target.global, *built.retainedArtifact, geometry, retried);
        REQUIRE(retry.reason == algo::SurfaceReason::None);
        REQUIRE(retried);
    }
    std::printf("surface publication: %ld faults / 0 escapes\n", faults);
    return faults;
}

}  // namespace

int main() {
    static_assert(std::is_nothrow_move_constructible_v<algo::SurfaceFaceRow>);
    static_assert(std::is_nothrow_move_constructible_v<algo::SurfaceSetOutcome>);
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::set_terminate([] {
        std::printf("TERMINATE in surface %s injection %ld (active %s)\n",
                    gSurfacePhase, gSurfaceInjection,
                    std::current_exception() ? "yes" : "no");
        std::_Exit(97);
    });
    REQUIRE(constructionSweep() >= 4);
    REQUIRE(publicationSweep() >= 4);
    if (gBad) {
        std::printf("SURFACE ADAPTER PROBE FAILED (%d)\n", gBad);
        return 1;
    }
    std::printf("SURFACE ADAPTER PROBE PASSED\n");
    return 0;
}
