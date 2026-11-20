// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// CK-2 end-to-end allocation containment probe. Run without sanitizers: this
// executable replaces global operator new.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paacs/algo/Meshers.hpp"
#include "paacs/algo/ProductionModelMesher.hpp"
#include "paacs/algo/ProductionSurfaceMesher.hpp"

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
        if (gFailAt != 0 && gCalls == gFailAt) throw std::bad_alloc();
    }
    if (void* memory = std::malloc(size == 0 ? 1 : size)) return memory;
    throw std::bad_alloc();
}
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete(void* memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void* memory, std::size_t) noexcept { std::free(memory); }
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
using namespace paacs;

#define REQUIRE(condition)                                                \
    do {                                                                  \
        if (!(condition)) {                                               \
            std::printf("FAIL %d: %s\n", __LINE__, #condition);          \
            ++gBad;                                                       \
        }                                                                 \
    } while (0)

class Edge final : public geo::GEdge {
public:
    Edge(int tag, const geo::GVertex* a, const geo::GVertex* b)
        : GEdge(tag), a_(a), b_(b) {}
    Vec3 pointAt(double t) const override {
        return a_->position() + (b_->position() - a_->position()) * t;
    }
    Vec3 tangentAt(double) const override {
        return (b_->position() - a_->position()).normalized();
    }
    double length() const override {
        return (b_->position() - a_->position()).norm();
    }
    int startVertexTag() const override { return a_->tag(); }
    int endVertexTag() const override { return b_->tag(); }
private:
    const geo::GVertex* a_;
    const geo::GVertex* b_;
};

class Face final : public geo::GFace {
public:
    Face(std::vector<int> loop, Vec2 hi)
        : GFace(1), loop_(std::move(loop)), hi_(hi) {}
    Vec3 pointAt(const Vec2& p) const override { return {p.u, p.v, 0}; }
    Vec2 paramAt(const Vec3& p) const override { return {p.x, p.y}; }
    Vec3 normalAt(const Vec2&) const override { return {0, 0, 1}; }
    void paramBounds(Vec2& lo, Vec2& hi) const override {
        lo = {0, 0}; hi = hi_;
    }
    std::vector<int> boundaryEdgeTags() const override { return loop_; }
    double curvatureAt(const Vec2&) const override { return 0; }
private:
    std::vector<int> loop_;
    Vec2 hi_;
};

class Solid final : public geo::GSolid {
public:
    Solid() : GSolid(1) {}
    bool isInside(const Vec3& p) const override {
        return p.x >= 0.0 && p.y >= 0.0 && p.z >= 0.0 &&
               p.x + p.y + p.z <= 1.0;
    }
    std::vector<int> boundaryFaceTags() const override { return {1}; }
    void bbox(Vec3& lo, Vec3& hi) const override {
        lo = {0, 0, 0}; hi = {1, 1, 1};
    }
};

class Model final : public geo::Model {
public:
    explicit Model(bool minimalVolume = false) {
        const std::vector<Vec2> ring = minimalVolume
            ? std::vector<Vec2>{{0, 0}, {1, 0}, {0, 1}}
            : std::vector<Vec2>{{0, 0}, {4, 0}, {4, 1},
                                {1, 1}, {1, 4}, {0, 4}};
        for (std::size_t i = 0; i < ring.size(); ++i)
            vertices_.push_back(std::make_unique<geo::GVertex>(
                static_cast<int>(i + 1),
                Vec3{ring[i].u, ring[i].v, 0}));
        if (minimalVolume)
            vertices_.push_back(
                std::make_unique<geo::GVertex>(4, Vec3{0, 0, 1}));
        std::vector<int> loop;
        for (std::size_t i = 0; i < ring.size(); ++i) {
            loop.push_back(static_cast<int>(i + 1));
            edges_.push_back(std::make_unique<Edge>(
                static_cast<int>(i + 1), vertices_[i].get(),
                vertices_[(i + 1) % ring.size()].get()));
        }
        if (minimalVolume) {
            edges_.push_back(std::make_unique<Edge>(
                4, vertices_[0].get(), vertices_[3].get()));
            edges_.push_back(std::make_unique<Edge>(
                5, vertices_[1].get(), vertices_[3].get()));
            edges_.push_back(std::make_unique<Edge>(
                6, vertices_[2].get(), vertices_[3].get()));
        }
        Vec2 hi{};
        for (const auto& p : ring) {
            hi.u = std::max(hi.u, p.u);
            hi.v = std::max(hi.v, p.v);
        }
        face_ = std::make_unique<Face>(std::move(loop), hi);
        if (minimalVolume) solid_ = std::make_unique<Solid>();
    }
    std::vector<const geo::GVertex*> vertices() const override {
        std::vector<const geo::GVertex*> out;
        for (const auto& x : vertices_) out.push_back(x.get());
        return out;
    }
    std::vector<const geo::GEdge*> edges() const override {
        std::vector<const geo::GEdge*> out;
        for (const auto& x : edges_) out.push_back(x.get());
        return out;
    }
    std::vector<const geo::GFace*> faces() const override {
        return {face_.get()};
    }
    std::vector<const geo::GSolid*> solids() const override {
        return solid_ ? std::vector<const geo::GSolid*>{solid_.get()}
                      : std::vector<const geo::GSolid*>{};
    }
    const geo::GVertex* vertex(int tag) const override {
        for (const auto& x : vertices_) if (x->tag() == tag) return x.get();
        return nullptr;
    }
    const geo::GEdge* edge(int tag) const override {
        for (const auto& x : edges_) if (x->tag() == tag) return x.get();
        return nullptr;
    }
    const geo::GFace* face(int tag) const override {
        return tag == 1 ? face_.get() : nullptr;
    }
    const geo::GSolid* solid(int tag) const override {
        return solid_ && tag == 1 ? solid_.get() : nullptr;
    }
private:
    std::vector<std::unique_ptr<geo::GVertex>> vertices_;
    std::vector<std::unique_ptr<geo::GEdge>> edges_;
    std::unique_ptr<Face> face_;
    std::unique_ptr<Solid> solid_;
};

struct Fixture {
    Model model;
    sizing::MeshSpec spec;
    mesh::Mesh mesh;
    Fixture() {
        // One segment per polygon side and no interior lattice points keeps
        // the exhaustive sweep fast while the re-entrant corner still forces
        // the qualified transverse boundary-recovery path.
        spec.field = std::make_shared<sizing::ConstantField>(10.0);
        algo::meshEdges(model, spec, mesh);
    }
};

struct ModelFixture {
    std::unique_ptr<geo::Kernel> kernel = geo::makeBuiltinKernel();
    std::unique_ptr<geo::Model> model =
        kernel->makeBox({0, 0, 0}, {1, 1, 1});
    sizing::MeshSpec spec;
    mesh::Mesh mesh;
    ModelFixture() {
        spec.field = std::make_shared<sizing::ConstantField>(10.0);
        spec.threads = 1;
        // The order-2 route is the allocation superset: it exercises the
        // same edge/surface/volume/audit/publication phases plus private
        // quadratic elevation before publication.
        spec.order = 2;
        mesh.setDedupTolerance(0.125);
        const Index a = mesh.addNode({-3, -2, -1}, {0, 91});
        const Index b = mesh.addNode({-2, -2, -1}, {0, 92});
        mesh.addCell(mesh::CellType::Line2, {a, b}, {1, 93});
    }
};

std::uint64_t digest(const mesh::Mesh& value) {
    std::uint64_t h = 1469598103934665603ULL;
    auto mix = [&](std::uint64_t x) { h = (h ^ x) * 1099511628211ULL; };
    std::uint64_t dedup = 0;
    const double dedupTolerance = value.dedupTolerance();
    std::memcpy(&dedup, &dedupTolerance, sizeof(dedup));
    mix(dedup);
    mix(static_cast<std::uint64_t>(value.nodeCount()));
    mix(static_cast<std::uint64_t>(value.cellCount()));
    for (Index i = 0; i < value.nodeCount(); ++i) {
        const auto& p = value.node(i);
        std::uint64_t x = 0, y = 0, z = 0;
        std::memcpy(&x, &p.x, sizeof(x));
        std::memcpy(&y, &p.y, sizeof(y));
        std::memcpy(&z, &p.z, sizeof(z));
        mix(x); mix(y); mix(z);
        mix(static_cast<std::uint64_t>(value.nodeClass(i).dim));
        mix(static_cast<std::uint64_t>(value.nodeClass(i).tag));
    }
    for (Index i = 0; i < value.cellCount(); ++i) {
        const auto& cell = value.cell(i);
        mix(static_cast<std::uint64_t>(cell.type));
        mix(static_cast<std::uint64_t>(cell.onEntity.dim));
        mix(static_cast<std::uint64_t>(cell.onEntity.tag));
        for (int j = 0; j < mesh::nodesPerCell(cell.type); ++j)
            mix(static_cast<std::uint64_t>(value.cellNodes(i)[j]));
    }
    return h;
}

}  // namespace

int main(int argc, char** argv) {
    const bool modelOnly = argc == 2 && std::strcmp(argv[1], "--model-only") == 0;
    const bool surfaceOnly = argc == 2 &&
                             std::strcmp(argv[1], "--surface-only") == 0;
    const bool cleanModelOnly = argc == 2 &&
                                std::strcmp(argv[1], "--clean-model") == 0;
    const bool countModelOnly = argc == 2 &&
                                std::strcmp(argv[1], "--count-model") == 0;
    bool expectModelCount = argc == 3 &&
                            std::strcmp(argv[1],
                                        "--expect-model-count") == 0;
    std::size_t expectedModelCount = 0;
    if (expectModelCount) {
        char* tail = nullptr;
        const auto expected = std::strtoull(argv[2], &tail, 10);
        expectModelCount = tail != nullptr && *tail == '\0' && expected > 0 &&
                           expected <=
                               (std::numeric_limits<std::size_t>::max)();
        if (expectModelCount)
            expectedModelCount = static_cast<std::size_t>(expected);
    }
    bool modelRange = argc == 4 &&
                      std::strcmp(argv[1], "--model-range") == 0;
    std::size_t modelRangeBegin = 0;
    std::size_t modelRangeEnd = 0;
    if (modelRange) {
        char* beginTail = nullptr;
        char* endTail = nullptr;
        const auto begin = std::strtoull(argv[2], &beginTail, 10);
        const auto end = std::strtoull(argv[3], &endTail, 10);
        modelRange = beginTail != nullptr && *beginTail == '\0' &&
                     endTail != nullptr && *endTail == '\0' && begin > 0 &&
                     end >= begin &&
                     end <= (std::numeric_limits<std::size_t>::max)();
        if (modelRange) {
            modelRangeBegin = static_cast<std::size_t>(begin);
            modelRangeEnd = static_cast<std::size_t>(end);
        }
    }
    if ((argc != 1 && argc != 2 && argc != 3 && argc != 4) ||
        (argc == 2 && !modelOnly && !surfaceOnly && !cleanModelOnly &&
         !countModelOnly) || (argc == 3 && !expectModelCount) ||
        (argc == 4 && !modelRange)) {
        std::printf("usage: probe_production_surface "
                    "[--surface-only|--model-only|--clean-model|"
                    "--count-model|--expect-model-count COUNT|"
                    "--model-range BEGIN END]\n");
        return 2;
    }
    if (countModelOnly || expectModelCount) {
        // Match the canonical sweep: its clean reference call warms any
        // process-local runtime/static state before ordinal injection begins.
        ModelFixture warm;
        const auto warmOutcome = algo::meshModelConstrained(
            *warm.model, warm.spec, warm.mesh);
        REQUIRE(warmOutcome.published());
        ModelFixture fixture;
        gCalls = 0;
        gFailAt = 0;
        gArmed = true;
        const auto outcome = algo::meshModelConstrained(
            *fixture.model, fixture.spec, fixture.mesh);
        gArmed = false;
        REQUIRE(outcome.published());
        if (expectModelCount) REQUIRE(gCalls == expectedModelCount);
        if (gBad != 0) return 1;
        std::printf("probe_production_surface model allocation count: %zu\n",
                    gCalls);
        return 0;
    }
    std::size_t faults = 0;
    std::size_t publicationFaults = 0;
    std::size_t toleratedFaults = 0;
    if (!modelOnly && !cleanModelOnly && !modelRange) {
        Fixture cleanFixture;
        const auto clean = algo::meshFacesConstrained(
            cleanFixture.model, cleanFixture.spec, cleanFixture.mesh);
        REQUIRE(clean.published());
        const auto cleanDigest = digest(cleanFixture.mesh);
        for (std::size_t fault = 1;; ++fault) {
            Fixture fixture;
            const auto before = digest(fixture.mesh);
            gCalls = 0;
            gFailAt = fault;
            gArmed = true;
            const auto outcome = algo::meshFacesConstrained(
                fixture.model, fixture.spec, fixture.mesh);
            gArmed = false;
            if (outcome.published()) {
                REQUIRE(digest(fixture.mesh) == cleanDigest);
                if (gCalls >= fault) {
                    // Some standard-library/algorithm internals may catch an
                    // allocation failure and complete through an exact
                    // fallback. That is a valid contained outcome, not the
                    // end-of-sweep sentinel; continue to the next ordinal.
                    ++toleratedFaults;
                    continue;
                }
                break;
            }
            ++faults;
            if (outcome.phase == algo::ProductionSurfacePhase::Publication) {
                ++publicationFaults;
                REQUIRE(outcome.reason ==
                        algo::ProductionSurfaceReason::AllocationFailed);
            }
            REQUIRE(outcome.destinationUnchanged);
            REQUIRE(digest(fixture.mesh) == before);
            const auto retry = algo::meshFacesConstrained(
                fixture.model, fixture.spec, fixture.mesh);
            REQUIRE(retry.published());
            REQUIRE(digest(fixture.mesh) == cleanDigest);
        }
        REQUIRE(faults > 0);
        REQUIRE(publicationFaults > 0);
    }
    std::size_t modelFaults = 0;
    std::size_t modelPublicationFaults = 0;
    std::size_t modelToleratedFaults = 0;
    if (!surfaceOnly) {
        ModelFixture cleanModelFixture;
        const auto cleanModel = algo::meshModelConstrained(
            *cleanModelFixture.model, cleanModelFixture.spec,
            cleanModelFixture.mesh);
        REQUIRE(cleanModel.published());
        const auto cleanModelDigest = digest(cleanModelFixture.mesh);
        if (cleanModelOnly) {
            std::printf("clean model published=%d status=%u reason=%u phase=%u "
                        "surfaceStatus=%u surfaceReason=%u nodes=%lld cells=%lld\n",
                        cleanModel.published() ? 1 : 0,
                        static_cast<unsigned>(cleanModel.status),
                        static_cast<unsigned>(cleanModel.reason),
                        static_cast<unsigned>(cleanModel.phase),
                        static_cast<unsigned>(cleanModel.surface.status),
                        static_cast<unsigned>(cleanModel.surface.reason),
                        (long long)cleanModelFixture.mesh.nodeCount(),
                        (long long)cleanModelFixture.mesh.cellCount());
            return gBad == 0 ? 0 : 1;
        }
        const std::size_t firstFault = modelRange ? modelRangeBegin : 1;
        for (std::size_t fault = firstFault;; ++fault) {
            ModelFixture fixture;
            const auto before = digest(fixture.mesh);
            gCalls = 0;
            gFailAt = fault;
            gArmed = true;
            const auto outcome = algo::meshModelConstrained(
                *fixture.model, fixture.spec, fixture.mesh);
            gArmed = false;
            if (outcome.published()) {
                REQUIRE(digest(fixture.mesh) == cleanModelDigest);
                if (gCalls >= fault) {
                    ++modelToleratedFaults;
                    if (modelRange && fault == modelRangeEnd) break;
                    continue;
                }
                REQUIRE(!modelRange);
                break;
            }
            ++modelFaults;
            if (outcome.phase == algo::ProductionModelPhase::Publication) {
                ++modelPublicationFaults;
                REQUIRE(outcome.reason ==
                        algo::ProductionModelReason::AllocationFailed);
            }
            REQUIRE(outcome.destinationUnchanged);
            REQUIRE(digest(fixture.mesh) == before);
            const auto retry = algo::meshModelConstrained(
                *fixture.model, fixture.spec, fixture.mesh);
            REQUIRE(retry.published());
            REQUIRE(digest(fixture.mesh) == cleanModelDigest);
            if (modelRange && fault == modelRangeEnd) break;
        }
        REQUIRE(modelFaults > 0);
        if (!modelRange) REQUIRE(modelPublicationFaults > 0);
    }
    if (gBad != 0) {
        std::printf("probe_production_surface FAILED: %d assertions\n", gBad);
        return 1;
    }
    std::printf("probe_production_surface PASSED: surface %zu structured/"
                "%zu tolerated/%zu publication; model %zu/%zu/%zu; zero "
                "escapes, exact rollback, clean retry", faults,
                toleratedFaults, publicationFaults, modelFaults,
                modelToleratedFaults, modelPublicationFaults);
    if (modelRange)
        std::printf("; model range %zu-%zu", modelRangeBegin, modelRangeEnd);
    std::printf("\n");
    return 0;
}
