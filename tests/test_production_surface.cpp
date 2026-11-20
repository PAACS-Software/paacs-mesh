// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <algorithm>
#include <array>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paacs/algo/Meshers.hpp"
#include "paacs/algo/ProductionSurfaceMesher.hpp"
#include "paacs/algo/SurfaceMesherAdapter.hpp"

#include "test_harness.hpp"

extern int gFailures;

namespace {
using namespace paacs;

struct FlatSnapshot {
    std::vector<double> coordinates;
    std::vector<std::int64_t> topology;
    friend bool operator==(const FlatSnapshot&, const FlatSnapshot&) = default;
};

FlatSnapshot snapshot(const mesh::Mesh& value) {
    FlatSnapshot out;
    out.topology.push_back(value.nodeCount());
    out.topology.push_back(value.cellCount());
    for (Index i = 0; i < value.nodeCount(); ++i) {
        const auto& p = value.node(i);
        out.coordinates.insert(out.coordinates.end(), {p.x, p.y, p.z});
        out.topology.push_back(value.nodeClass(i).dim);
        out.topology.push_back(value.nodeClass(i).tag);
    }
    for (Index i = 0; i < value.cellCount(); ++i) {
        const auto& c = value.cell(i);
        out.topology.push_back(static_cast<std::int64_t>(c.type));
        out.topology.push_back(c.onEntity.dim);
        out.topology.push_back(c.onEntity.tag);
        for (int n = 0; n < mesh::nodesPerCell(c.type); ++n)
            out.topology.push_back(value.cellNodes(i)[n]);
    }
    return out;
}

bool prefixPreserved(const FlatSnapshot& before, const mesh::Mesh& after) {
    if (after.nodeCount() < before.topology[0] ||
        after.cellCount() < before.topology[1])
        return false;
    FlatSnapshot prefix;
    const Index nodeCount = static_cast<Index>(before.topology[0]);
    const Index cellCount = static_cast<Index>(before.topology[1]);
    prefix.topology.push_back(nodeCount);
    prefix.topology.push_back(cellCount);
    for (Index i = 0; i < nodeCount; ++i) {
        const auto& p = after.node(i);
        prefix.coordinates.insert(prefix.coordinates.end(), {p.x, p.y, p.z});
        prefix.topology.push_back(after.nodeClass(i).dim);
        prefix.topology.push_back(after.nodeClass(i).tag);
    }
    for (Index i = 0; i < cellCount; ++i) {
        const auto& c = after.cell(i);
        prefix.topology.push_back(static_cast<std::int64_t>(c.type));
        prefix.topology.push_back(c.onEntity.dim);
        prefix.topology.push_back(c.onEntity.tag);
        for (int n = 0; n < mesh::nodesPerCell(c.type); ++n)
            prefix.topology.push_back(after.cellNodes(i)[n]);
    }
    return prefix == before;
}

class TestEdge final : public geo::GEdge {
public:
    TestEdge(int tag, const geo::GVertex* a, const geo::GVertex* b)
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

class PolygonFace final : public geo::GFace {
public:
    PolygonFace(int tag, std::vector<int> loop, Vec2 hi,
                bool failGeometry = false)
        : GFace(tag), loop_(std::move(loop)), hi_(hi),
          failGeometry_(failGeometry) {}
    Vec3 pointAt(const Vec2& uv) const override {
        if (failGeometry_) throw std::runtime_error("geometry provider failure");
        return {uv.u, uv.v, 0};
    }
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
    bool failGeometry_;
};

class PolygonModel final : public geo::Model {
public:
    PolygonModel(std::vector<Vec2> ring, bool failGeometry = false) {
        for (std::size_t i = 0; i < ring.size(); ++i)
            vertices_.push_back(std::make_unique<geo::GVertex>(
                static_cast<int>(i + 1), Vec3{ring[i].u, ring[i].v, 0}));
        std::vector<int> loop;
        for (std::size_t i = 0; i < ring.size(); ++i) {
            const int tag = static_cast<int>(i + 1);
            loop.push_back(tag);
            edges_.push_back(std::make_unique<TestEdge>(
                tag, vertices_[i].get(),
                vertices_[(i + 1) % ring.size()].get()));
        }
        Vec2 hi{};
        for (const auto& p : ring) {
            hi.u = std::max(hi.u, p.u);
            hi.v = std::max(hi.v, p.v);
        }
        face_ = std::make_unique<PolygonFace>(1, std::move(loop), hi,
                                              failGeometry);
    }
    std::vector<const geo::GVertex*> vertices() const override {
        std::vector<const geo::GVertex*> out;
        for (const auto& p : vertices_) out.push_back(p.get());
        return out;
    }
    std::vector<const geo::GEdge*> edges() const override {
        std::vector<const geo::GEdge*> out;
        for (const auto& p : edges_) out.push_back(p.get());
        return out;
    }
    std::vector<const geo::GFace*> faces() const override {
        return {face_.get()};
    }
    std::vector<const geo::GSolid*> solids() const override { return {}; }
    const geo::GVertex* vertex(int tag) const override {
        for (const auto& p : vertices_) if (p->tag() == tag) return p.get();
        return nullptr;
    }
    const geo::GEdge* edge(int tag) const override {
        for (const auto& p : edges_) if (p->tag() == tag) return p.get();
        return nullptr;
    }
    const geo::GFace* face(int tag) const override {
        return face_->tag() == tag ? face_.get() : nullptr;
    }
    const geo::GSolid* solid(int) const override { return nullptr; }
private:
    std::vector<std::unique_ptr<geo::GVertex>> vertices_;
    std::vector<std::unique_ptr<geo::GEdge>> edges_;
    std::unique_ptr<PolygonFace> face_;
};

class ReverseFaceView final : public geo::Model {
public:
    explicit ReverseFaceView(const geo::Model& source) : source_(source) {}
    std::vector<const geo::GVertex*> vertices() const override {
        return source_.vertices();
    }
    std::vector<const geo::GEdge*> edges() const override {
        return source_.edges();
    }
    std::vector<const geo::GFace*> faces() const override {
        auto out = source_.faces();
        std::reverse(out.begin(), out.end());
        return out;
    }
    std::vector<const geo::GSolid*> solids() const override {
        return source_.solids();
    }
    const geo::GVertex* vertex(int tag) const override {
        return source_.vertex(tag);
    }
    const geo::GEdge* edge(int tag) const override {
        return source_.edge(tag);
    }
    const geo::GFace* face(int tag) const override {
        return source_.face(tag);
    }
    const geo::GSolid* solid(int tag) const override {
        return source_.solid(tag);
    }
private:
    const geo::Model& source_;
};

sizing::MeshSpec constantSpec(double h) {
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(h);
    return spec;
}

void testConvexAndSharedSet() {
    TEST("production_surface.convex_shared_set");
    auto kernel = geo::makeBuiltinKernel();
    auto model = kernel->makeBox({0, 0, 0}, {1, 1, 1});
    auto spec = constantSpec(0.25);
    mesh::Mesh mesh;
    algo::meshEdges(*model, spec, mesh);
    const auto before = snapshot(mesh);
    const auto result = algo::meshFacesConstrained(*model, spec, mesh);
    CHECK(result.published());
    CHECK(result.reportComplete);
    CHECK(result.processedFaces == 6);
    CHECK(result.unattemptedFaces == 0);
    CHECK(result.faceMetrics.size() == 6);
    CHECK(result.retainedArtifact != nullptr);
    CHECK(result.retainedArtifact->faces.size() == 6);
    CHECK(mesh.countCells(mesh::CellType::Tri3) > 0);
    CHECK(prefixPreserved(before, mesh));
    for (const auto& row : result.faceMetrics) {
        CHECK(row.construction.runComplete);
        CHECK(row.construction.fullScanFallbacks == 0);
        CHECK(row.protectedSegments == row.recoveredProtectedSegments);
        CHECK(row.spacingIndexHighWater >= row.boundaryNodes + 4);
    }
    for (const auto& face : result.retainedArtifact->faces)
        CHECK(face != nullptr && face->audit.passed());
}

void testConcaveRecovery() {
    TEST("production_surface.concave_boundary_recovery");
    PolygonModel model({{0, 0}, {4, 0}, {4, 1},
                        {1, 1}, {1, 4}, {0, 4}});
    auto spec = constantSpec(0.8);
    mesh::Mesh mesh;
    algo::meshEdges(model, spec, mesh);
    const auto before = snapshot(mesh);
    const auto result = algo::meshFacesConstrained(model, spec, mesh);
    CHECK(result.published());
    CHECK(result.faceMetrics.size() == 1);
    CHECK(result.faceMetrics[0].protectedSegments ==
          result.faceMetrics[0].recoveredProtectedSegments);
    CHECK(result.faceMetrics[0].construction.fullScanFallbacks == 0);
    CHECK(prefixPreserved(before, mesh));
    CHECK(mesh.countCells(mesh::CellType::Tri3) > 0);
    CHECK(result.retainedArtifact->faces[0]->audit.passed());
}

void testCurvedOrientation() {
    TEST("production_surface.curved_orientation");
    auto kernel = geo::makeBuiltinKernel();
    auto model = kernel->makeCylinder({0, 0, 0}, 0.5, 1.0);
    auto spec = constantSpec(0.2);
    mesh::Mesh mesh;
    algo::meshEdges(*model, spec, mesh);
    const auto result = algo::meshFacesConstrained(*model, spec, mesh);
    CHECK(result.published());
    CHECK(result.retainedArtifact != nullptr);
    CHECK(result.retainedArtifact->faces.size() == 4);
    for (const auto& face : result.retainedArtifact->faces)
        CHECK(face != nullptr && face->audit.physicalOrientation);
}

void testFailureInvariance() {
    TEST("production_surface.structured_failures_invariant");
    PolygonModel model({{0, 0}, {2, 0}, {2, 2}, {0, 2}});
    auto spec = constantSpec(0.5);
    {
        mesh::Mesh mesh;
        const auto before = snapshot(mesh);
        const auto result = algo::meshFacesConstrained(model, spec, mesh);
        CHECK(result.status == algo::ProductionSurfaceStatus::Rejected);
        CHECK(result.reason ==
              algo::ProductionSurfaceReason::MissingBoundaryDiscretization);
        CHECK(result.destinationUnchanged);
        CHECK(snapshot(mesh) == before);
    }
    {
        PolygonModel failing({{0, 0}, {2, 0}, {2, 2}, {0, 2}}, true);
        mesh::Mesh mesh;
        algo::meshEdges(failing, spec, mesh);
        const auto before = snapshot(mesh);
        const auto result = algo::meshFacesConstrained(failing, spec, mesh);
        CHECK(result.status == algo::ProductionSurfaceStatus::Failed);
        CHECK(result.reason ==
              algo::ProductionSurfaceReason::GeometryProviderFailure);
        CHECK(result.destinationUnchanged);
        CHECK(snapshot(mesh) == before);
    }
}

void testDeterminism() {
    TEST("production_surface.deterministic_ordering");
    auto run = [] {
        PolygonModel model({{0, 0}, {5, 0}, {5, 2},
                            {2, 2}, {2, 5}, {0, 5}});
        auto spec = constantSpec(0.7);
        mesh::Mesh mesh;
        algo::meshEdges(model, spec, mesh);
        auto result = algo::meshFacesConstrained(model, spec, mesh);
        return std::make_pair(snapshot(mesh), std::move(result));
    };
    auto a = run();
    auto b = run();
    CHECK(a.second.published());
    CHECK(b.second.published());
    CHECK(a.first == b.first);
    CHECK(a.second.faceMetrics.size() == b.second.faceMetrics.size());
    CHECK(a.second.faceMetrics[0].spacingQueries ==
          b.second.faceMetrics[0].spacingQueries);
    CHECK(a.second.faceMetrics[0].spacingDistanceTests ==
          b.second.faceMetrics[0].spacingDistanceTests);
    CHECK(a.second.faceMetrics[0].construction.trianglesCommitted ==
          b.second.faceMetrics[0].construction.trianglesCommitted);

    auto kernel = geo::makeBuiltinKernel();
    auto model = kernel->makeBox({0, 0, 0}, {1, 1, 1});
    ReverseFaceView reversed(*model);
    auto spec = constantSpec(0.3);
    mesh::Mesh forwardMesh, reversedMesh;
    algo::meshEdges(*model, spec, forwardMesh);
    algo::meshEdges(*model, spec, reversedMesh);
    const auto forward =
        algo::meshFacesConstrained(*model, spec, forwardMesh);
    const auto backward =
        algo::meshFacesConstrained(reversed, spec, reversedMesh);
    CHECK(forward.published());
    CHECK(backward.published());
    CHECK(snapshot(forwardMesh) == snapshot(reversedMesh));
}

void testSpacingGrowth() {
    TEST("production_surface.spacing_index_growth");
    auto run = [](double h) {
        PolygonModel model({{0, 0}, {8, 0}, {8, 8}, {0, 8}});
        auto spec = constantSpec(h);
        mesh::Mesh mesh;
        algo::meshEdges(model, spec, mesh);
        auto result = algo::meshFacesConstrained(model, spec, mesh);
        CHECK(result.published());
        return result.faceMetrics.front();
    };
    const auto coarse = run(1.0);
    const auto fine = run(0.5);
    CHECK(fine.interiorCandidates > coarse.interiorCandidates * 3);
    CHECK(fine.spacingDistanceTests < coarse.spacingDistanceTests * 7 + 64);
    CHECK(fine.spacingQueries < coarse.spacingQueries * 5);
}

}  // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };

const Entry* productionSurfaceTests(std::size_t& count) {
    static const Entry tests[] = {
        {"production_surface.convex_shared", testConvexAndSharedSet},
        {"production_surface.concave", testConcaveRecovery},
        {"production_surface.curved", testCurvedOrientation},
        {"production_surface.failures", testFailureInvariance},
        {"production_surface.determinism", testDeterminism},
        {"production_surface.spacing", testSpacingGrowth},
    };
    count = sizeof(tests) / sizeof(tests[0]);
    return tests;
}
}  // namespace paacs_tests
