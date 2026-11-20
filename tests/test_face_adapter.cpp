// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// BND-2A Increment 5 bounded synthetic FaceMesher integration proofs.
#include <algorithm>
#include <array>
#include <cstdio>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "paacs/algo/FaceMesherAdapter.hpp"

extern int gFailures;
static const char* gFaceCurrent = "";
#define FACE_TEST(name)                     \
    do {                                    \
        gFaceCurrent = name;                \
        std::printf("[test] %s\n", name);   \
    } while (0)
#define FACE_CHECK(cond)                                                   \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::printf("  FAIL %s:%d  %s  (in %s)\n", __FILE__, __LINE__, \
                        #cond, gFaceCurrent);                              \
            ++gFailures;                                                   \
        }                                                                  \
    } while (0)

namespace {

namespace algo = paacs::algo;
namespace bnd = paacs::bnd;
namespace mesh = paacs::mesh;
using paacs::Index;
using paacs::Vec2;
using paacs::Vec3;

bool sameClass(const mesh::Classification& a,
               const mesh::Classification& b) {
    return a.dim == b.dim && a.tag == b.tag;
}
bool samePoint(const Vec3& a, const Vec3& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

struct MeshSnapshot {
    std::vector<Vec3> nodes;
    std::vector<mesh::Classification> nodeClass;
    struct Cell {
        mesh::CellType type{};
        mesh::Classification cls{};
        std::vector<Index> conn;
    };
    std::vector<Cell> cells;
};

MeshSnapshot snapshot(const mesh::Mesh& m) {
    MeshSnapshot s;
    for (Index i = 0; i < m.nodeCount(); ++i) {
        s.nodes.push_back(m.node(i));
        s.nodeClass.push_back(m.nodeClass(i));
    }
    for (Index i = 0; i < m.cellCount(); ++i) {
        MeshSnapshot::Cell c;
        c.type = m.cell(i).type;
        c.cls = m.cell(i).onEntity;
        const int n = mesh::nodesPerCell(c.type);
        const Index* p = m.cellNodes(i);
        c.conn.assign(p, p + n);
        s.cells.push_back(std::move(c));
    }
    return s;
}

bool matches(const mesh::Mesh& m, const MeshSnapshot& s) {
    if (m.nodeCount() != static_cast<Index>(s.nodes.size()) ||
        m.cellCount() != static_cast<Index>(s.cells.size()))
        return false;
    for (Index i = 0; i < m.nodeCount(); ++i)
        if (!samePoint(m.node(i), s.nodes[std::size_t(i)]) ||
            !sameClass(m.nodeClass(i), s.nodeClass[std::size_t(i)]))
            return false;
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto& want = s.cells[std::size_t(i)];
        if (m.cell(i).type != want.type ||
            !sameClass(m.cell(i).onEntity, want.cls))
            return false;
        const Index* p = m.cellNodes(i);
        for (std::size_t j = 0; j < want.conn.size(); ++j)
            if (p[j] != want.conn[j]) return false;
    }
    return true;
}

struct Provider : bnd::GeometryProvider {
    bnd::GeometryProviderStamp value{700, 3};
    int mode = 0;  // 1 bad_alloc, 2 other, 3 domain, 4 nonfinite, 5 stamp throw
    bnd::GeometryProviderStamp stamp() const override {
        if (mode == 5) throw std::runtime_error("stamp");
        return value;
    }
    Vec3 pointAt(bnd::ChartId, const Vec2& uv) const override {
        if (mode == 1) throw std::bad_alloc();
        if (mode == 2) throw std::runtime_error("mapping");
        if (mode == 3) throw bnd::GeometryProviderDomainError("domain");
        if (mode == 4)
            return {std::numeric_limits<double>::quiet_NaN(), uv.v, 0};
        return {uv.u, uv.v, 0};
    }
};

struct Fixture {
    mesh::Mesh global;
    algo::SyntheticFaceInput input;
    std::vector<Index> boundaryGlobal;
};

Fixture makeFace(const std::vector<Vec2>& polygon,
                 const std::vector<std::array<std::size_t, 3>>& triangles,
                 const std::vector<std::array<std::size_t, 2>>& interiors,
                 bool duplicateInteriorCoordinate = false) {
    Fixture f;
    // An unrelated pre-existing face proves outside/global-prefix invariance.
    const Index o0 = f.global.addNode({-20, -20, 1}, {2, 99});
    const Index o1 = f.global.addNode({-19, -20, 1}, {2, 99});
    const Index o2 = f.global.addNode({-20, -19, 1}, {2, 99});
    f.global.addCell(mesh::CellType::Tri3, {o0, o1, o2}, {2, 99});
    if (duplicateInteriorCoordinate)
        f.global.addNode({1, 1, 0}, {0, 777});

    f.input.chart = 7;
    f.input.faceClass = {2, 5};
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const Index global = f.global.addNode(
            {polygon[i].u, polygon[i].v, 0},
            {0, 100 + static_cast<int>(i)});
        f.boundaryGlobal.push_back(global);
        algo::SyntheticFaceNode n;
        n.uv = polygon[i];
        n.classification = f.global.nodeClass(global);
        n.existingGlobalIndex = global;
        n.constraintVertex = 1000 + i;
        f.input.nodes.push_back(n);
        f.input.outerLoop.push_back(i);
    }
    for (const auto& t : triangles)
        f.input.triangles.push_back({t});
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        bnd::ConstraintSegment s;
        s.id = 10 + i;
        s.a = 1000 + i;
        s.b = 1000 + ((i + 1) % polygon.size());
        s.chart = 7;
        s.cls = {1, 100 + static_cast<int>(i)};
        f.input.constraints.push_back(s);
    }
    for (std::size_t i = 0; i < interiors.size(); ++i) {
        bnd::ConstraintSegment s;
        s.id = 500 + i;
        s.a = 1000 + interiors[i][0];
        s.b = 1000 + interiors[i][1];
        s.chart = 7;
        s.cls = {1, 9 + static_cast<int>(i)};
        f.input.constraints.push_back(s);
    }
    return f;
}

Fixture convexFixture(bool duplicateCenter = false) {
    return makeFace({{0, 0}, {2, 0}, {2, 2}, {0, 2}},
                    {{0, 1, 2}, {0, 2, 3}}, {{1, 3}},
                    duplicateCenter);
}

Fixture concaveFixture() {
    // Strong L-shaped reflex corner at (1,1). Constraint 0->2 remains wholly
    // in the bottom arm and crosses the existing 1->3 diagonal.
    return makeFace(
        {{0, 0}, {4, 0}, {4, 1}, {1, 1}, {1, 4}, {0, 4}},
        {{0, 1, 3}, {1, 2, 3}, {0, 3, 5}, {3, 4, 5}}, {{0, 2}});
}

Fixture multipleFixture() {
    // Two parallel, disjoint prescribed diagonals in adjacent quads.
    return makeFace(
        {{0, 0}, {2, 0}, {4, 0}, {4, 2}, {2, 2}, {0, 2}},
        {{0, 1, 4}, {0, 4, 5}, {1, 2, 3}, {1, 3, 4}},
        {{1, 5}, {2, 4}});
}

bnd::SteinerPolicy policy() {
    bnd::SteinerPolicy p;
    p.splitBudget = 8;
    p.minFeatureLevel = 52;
    return p;
}

void testConvexPublicationIdentity() {
    FACE_TEST("face_adapter.convex (audited local recovery + identity publication)");
    auto f = convexFixture(/*duplicateCenter=*/true);
    Provider provider;
    const MeshSnapshot before = snapshot(f.global);
    const auto built = algo::buildSyntheticFaceMesh(
        f.input, f.global, provider, provider.value, policy());
    FACE_CHECK(built.status == algo::FaceAdapterStatus::Ready);
    FACE_CHECK(built.result != nullptr);
    FACE_CHECK(matches(f.global, before));  // build is strictly local
    if (!built.result) return;
    FACE_CHECK(built.result->audit.passed());
    FACE_CHECK(built.result->localWorkspace != nullptr);
    FACE_CHECK(built.result->existingGlobalIndex.size() == 4);
    FACE_CHECK(built.recovery.processed.size() == 5);
    if (built.recovery.processed.size() == 5) {
        for (std::size_t i = 0; i < 4; ++i)
            FACE_CHECK(built.recovery.processed[i].result.status ==
                       bnd::SteinerStatus::AlreadyPresent);
        FACE_CHECK(built.recovery.processed[4].result.status ==
                   bnd::SteinerStatus::Recovered);
    }

    const Index nodesBefore = f.global.nodeCount();
    const Index cellsBefore = f.global.cellCount();
    const auto publication =
        algo::publishFaceMeshResult(f.global, *built.result);
    FACE_CHECK(publication.status == algo::FacePublishStatus::Published);
    FACE_CHECK(f.global.nodeCount() == nodesBefore + 1);
    FACE_CHECK(f.global.cellCount() ==
               cellsBefore +
                   static_cast<Index>(built.result->localCells.size()));
    // An exact-coordinate node already existed globally, but publication must
    // append the new curve-classified Steiner identity via addNode().
    int centres = 0;
    for (Index i = 0; i < f.global.nodeCount(); ++i)
        if (samePoint(f.global.node(i), {1, 1, 0})) ++centres;
    FACE_CHECK(centres == 2);
    for (const auto& [local, global] :
         built.result->existingGlobalIndex)
        FACE_CHECK(publication.localToGlobal[local] == global);
}

void testConcaveAndMultiple() {
    FACE_TEST("face_adapter.concave_multiple (strong reflex + two constraints)");
    Provider provider;
    for (int which = 0; which < 2; ++which) {
        auto f = which == 0 ? concaveFixture() : multipleFixture();
        const auto before = snapshot(f.global);
        const auto built = algo::buildSyntheticFaceMesh(
            f.input, f.global, provider, provider.value, policy());
        FACE_CHECK(built.status == algo::FaceAdapterStatus::Ready);
        FACE_CHECK(built.result != nullptr);
        FACE_CHECK(matches(f.global, before));
        if (!built.result) continue;
        FACE_CHECK(built.result->audit.oneOuterLoop);
        FACE_CHECK(built.result->audit.constraintsRecovered);
        int recovered = 0;
        for (const auto& row : built.recovery.processed)
            if (row.result.status == bnd::SteinerStatus::Recovered)
                ++recovered;
        FACE_CHECK(recovered == (which == 0 ? 1 : 2));
        const auto published =
            algo::publishFaceMeshResult(f.global, *built.result);
        FACE_CHECK(published.status ==
                   algo::FacePublishStatus::Published);
    }
}

void testDeterministicMappingsAndBoundaryReuse() {
    FACE_TEST("face_adapter.determinism (container order + boundary identity)");
    auto a = multipleFixture();
    auto b = multipleFixture();
    std::reverse(b.input.constraints.begin(), b.input.constraints.end());
    Provider provider;
    auto ra = algo::buildSyntheticFaceMesh(
        a.input, a.global, provider, provider.value, policy());
    auto rb = algo::buildSyntheticFaceMesh(
        b.input, b.global, provider, provider.value, policy());
    FACE_CHECK(ra.status == algo::FaceAdapterStatus::Ready);
    FACE_CHECK(rb.status == algo::FaceAdapterStatus::Ready);
    if (!ra.result || !rb.result) return;
    FACE_CHECK(ra.result->localNodes.size() == rb.result->localNodes.size());
    FACE_CHECK(ra.result->localCells.size() == rb.result->localCells.size());
    FACE_CHECK(ra.result->outerLoop == rb.result->outerLoop);
    FACE_CHECK(ra.result->existingGlobalIndex ==
               rb.result->existingGlobalIndex);
    for (std::size_t i = 0; i < ra.result->localNodes.size(); ++i) {
        FACE_CHECK(samePoint(ra.result->localNodes[i].coord,
                             rb.result->localNodes[i].coord));
        FACE_CHECK(sameClass(ra.result->localNodes[i].cls,
                             rb.result->localNodes[i].cls));
    }
    for (std::size_t i = 0; i < ra.result->localCells.size(); ++i)
        FACE_CHECK(ra.result->localCells[i].nodes ==
                   rb.result->localCells[i].nodes);
    FACE_CHECK(ra.result->localWorkspace->semanticTopologyHash() ==
               rb.result->localWorkspace->semanticTopologyHash());
    for (const std::size_t local : ra.result->outerLoop) {
        const auto it = ra.result->existingGlobalIndex.find(local);
        FACE_CHECK(it != ra.result->existingGlobalIndex.end());
        if (it != ra.result->existingGlobalIndex.end())
            FACE_CHECK(std::find(a.boundaryGlobal.begin(),
                                 a.boundaryGlobal.end(), it->second) !=
                       a.boundaryGlobal.end());
    }
}

void testProviderAndInputFailures() {
    FACE_TEST("face_adapter.failures (provider taxonomy + global invariance)");
    struct Case {
        int mode;
        algo::FaceAdapterStatus status;
        bnd::SteinerReject reject;
        bnd::SteinerFail fail;
    };
    const Case cases[] = {
        {1, algo::FaceAdapterStatus::Failed, bnd::SteinerReject::None,
         bnd::SteinerFail::CandidateAllocFailed},
        {2, algo::FaceAdapterStatus::Failed, bnd::SteinerReject::None,
         bnd::SteinerFail::GeometryProviderError},
        {3, algo::FaceAdapterStatus::Rejected,
         bnd::SteinerReject::ProviderDomainRejected, bnd::SteinerFail::None},
        {4, algo::FaceAdapterStatus::Rejected,
         bnd::SteinerReject::NonfiniteMappedPosition, bnd::SteinerFail::None},
    };
    for (const auto& c : cases) {
        auto f = convexFixture();
        const auto before = snapshot(f.global);
        Provider provider;
        provider.mode = c.mode;
        const auto out = algo::buildSyntheticFaceMesh(
            f.input, f.global, provider, provider.value, policy());
        FACE_CHECK(out.status == c.status);
        FACE_CHECK(out.reject == c.reject);
        FACE_CHECK(out.fail == c.fail);
        FACE_CHECK(out.result == nullptr);
        FACE_CHECK(matches(f.global, before));
    }
    {
        auto f = convexFixture();
        const auto before = snapshot(f.global);
        Provider provider;
        const bnd::GeometryProviderStamp stale{provider.value.identity,
                                                provider.value.revision + 1};
        const auto out = algo::buildSyntheticFaceMesh(
            f.input, f.global, provider, stale, policy());
        FACE_CHECK(out.status == algo::FaceAdapterStatus::Rejected);
        FACE_CHECK(out.reject ==
                   bnd::SteinerReject::StaleGeometryProvider);
        FACE_CHECK(matches(f.global, before));
    }
    {
        auto f = convexFixture();
        f.input.periodic = true;
        const auto before = snapshot(f.global);
        Provider provider;
        const auto out = algo::buildSyntheticFaceMesh(
            f.input, f.global, provider, provider.value, policy());
        FACE_CHECK(out.status == algo::FaceAdapterStatus::Refused);
        FACE_CHECK(out.refuse == bnd::SteinerRefuse::UnsupportedChart);
        FACE_CHECK(matches(f.global, before));
    }
}

void testPublicationFailureAtomicity() {
    FACE_TEST("face_adapter.publication_atomicity (invalid result is no-op)");
    auto f = convexFixture();
    Provider provider;
    auto built = algo::buildSyntheticFaceMesh(
        f.input, f.global, provider, provider.value, policy());
    FACE_CHECK(built.status == algo::FaceAdapterStatus::Ready);
    if (!built.result) return;

    // Existing identity changed after local recovery: publication rejects
    // before reserve/append and preserves the changed pre-publication state.
    const auto mapped = *built.result->existingGlobalIndex.begin();
    f.global.setNode(mapped.second, {99, 99, 99});
    const auto changed = snapshot(f.global);
    const auto rejected =
        algo::publishFaceMeshResult(f.global, *built.result);
    FACE_CHECK(rejected.status == algo::FacePublishStatus::Rejected);
    FACE_CHECK(rejected.reason ==
               algo::FacePublishReason::ExistingIdentityMismatch);
    FACE_CHECK(matches(f.global, changed));

    // Restore identity, then corrupt an exported cell. Re-audit at
    // publication catches it and again leaves the mesh unchanged.
    f.global.setNode(mapped.second,
                     built.result->localNodes[mapped.first].coord);
    built.result->localCells[0].nodes[1] =
        built.result->localCells[0].nodes[0];
    const auto beforeBadResult = snapshot(f.global);
    const auto invalid =
        algo::publishFaceMeshResult(f.global, *built.result);
    FACE_CHECK(invalid.status == algo::FacePublishStatus::Rejected);
    FACE_CHECK(invalid.reason == algo::FacePublishReason::InvalidResult);
    FACE_CHECK(matches(f.global, beforeBadResult));
}

}  // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* faceAdapterTests(std::size_t& count) {
    static const Entry tests[] = {
        {"face_adapter.convex", testConvexPublicationIdentity},
        {"face_adapter.concave_multiple", testConcaveAndMultiple},
        {"face_adapter.determinism", testDeterministicMappingsAndBoundaryReuse},
        {"face_adapter.failures", testProviderAndInputFailures},
        {"face_adapter.publication_atomicity", testPublicationFailureAtomicity},
    };
    count = sizeof(tests) / sizeof(tests[0]);
    return tests;
}
}  // namespace paacs_tests
