// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// BND-2A Increment 6 surface-set contract proofs.
#include <algorithm>
#include <bit>
#include <cstdio>
#include <limits>
#include <set>
#include <stdexcept>
#include <type_traits>

#include "paacs/algo/SurfaceMesherAdapter.hpp"

extern int gFailures;
static const char* gSurfaceCurrent = "";
#define SURFACE_TEST(name)                  \
    do {                                    \
        gSurfaceCurrent = name;             \
        std::printf("[test] %s\n", name);  \
    } while (0)
#define SURFACE_CHECK(cond)                                                \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::printf("  FAIL %s:%d  %s  (in %s)\n", __FILE__, __LINE__, \
                        #cond, gSurfaceCurrent);                           \
            ++gFailures;                                                   \
        }                                                                  \
    } while (0)

namespace {
namespace algo = paacs::algo;
namespace bnd = paacs::bnd;
namespace mesh = paacs::mesh;
using paacs::Vec2;
using paacs::Vec3;

struct FlatMeshSnapshot {
    std::vector<double> coordinates;
    std::vector<std::int64_t> topology;
    friend bool operator==(const FlatMeshSnapshot&, const FlatMeshSnapshot&) =
        default;
};

FlatMeshSnapshot snapshotMesh(const mesh::Mesh& value) {
    FlatMeshSnapshot out;
    out.topology.push_back(value.nodeCount());
    out.topology.push_back(value.cellCount());
    for (paacs::Index i = 0; i < value.nodeCount(); ++i) {
        const auto& p = value.node(i);
        out.coordinates.insert(out.coordinates.end(), {p.x, p.y, p.z});
        out.topology.push_back(value.nodeClass(i).dim);
        out.topology.push_back(value.nodeClass(i).tag);
    }
    for (paacs::Index i = 0; i < value.cellCount(); ++i) {
        const auto& cell = value.cell(i);
        out.topology.push_back(static_cast<std::int64_t>(cell.type));
        out.topology.push_back(cell.onEntity.dim);
        out.topology.push_back(cell.onEntity.tag);
        for (int j = 0; j < mesh::nodesPerCell(cell.type); ++j)
            out.topology.push_back(value.cellNodes(i)[j]);
    }
    return out;
}

bool sameVertexUse(const algo::SurfaceChainVertex& a,
                   const algo::SurfaceChainVertex& b) {
    return a.faceLocalVertexId == b.faceLocalVertexId &&
           a.sharedKey == b.sharedKey;
}

bool sameRetainedFace(const algo::SurfaceFaceResult& a,
                      const algo::SurfaceFaceResult& b) {
    if (a.boundaryLoops != b.boundaryLoops ||
        a.existingGlobalIndex != b.existingGlobalIndex ||
        a.normalizedChains.size() != b.normalizedChains.size() ||
        a.segmentOrigins.size() != b.segmentOrigins.size() ||
        a.persistentEdges.size() != b.persistentEdges.size() ||
        a.refinedIdentityByLocalNode != b.refinedIdentityByLocalNode ||
        a.localMappingParams.size() != b.localMappingParams.size() ||
        a.localCells.size() != b.localCells.size())
        return false;
    for (std::size_t i = 0; i < a.normalizedChains.size(); ++i) {
        const auto& x = a.normalizedChains[i];
        const auto& y = b.normalizedChains[i];
        if (x.chainPersistentId != y.chainPersistentId || x.role != y.role ||
            x.classification.dim != y.classification.dim ||
            x.classification.tag != y.classification.tag ||
            x.closed != y.closed || x.segmentEntityIds != y.segmentEntityIds ||
            x.owningCurve != y.owningCurve || x.uses.size() != y.uses.size())
            return false;
        for (std::size_t j = 0; j < x.uses.size(); ++j)
            if (!sameVertexUse(x.uses[j], y.uses[j])) return false;
    }
    auto ai = a.segmentOrigins.begin();
    auto bi = b.segmentOrigins.begin();
    for (; ai != a.segmentOrigins.end(); ++ai, ++bi) {
        const auto& x = ai->second;
        const auto& y = bi->second;
        if (ai->first != bi->first ||
            x.chainPersistentId != y.chainPersistentId ||
            x.normalizedSegmentOrdinal != y.normalizedSegmentOrdinal ||
            x.classification.dim != y.classification.dim ||
            x.classification.tag != y.classification.tag ||
            !sameVertexUse(x.endpoints[0], y.endpoints[0]) ||
            !sameVertexUse(x.endpoints[1], y.endpoints[1]))
            return false;
    }
    for (std::size_t i = 0; i < a.persistentEdges.size(); ++i) {
        const auto& x = a.persistentEdges[i];
        const auto& y = b.persistentEdges[i];
        if (x.bndEntityId != y.bndEntityId ||
            x.chainPersistentId != y.chainPersistentId ||
            x.normalizedSegmentOrdinal != y.normalizedSegmentOrdinal ||
            x.role != y.role || x.endpointSharedKeys != y.endpointSharedKeys ||
            x.classification.dim != y.classification.dim ||
            x.classification.tag != y.classification.tag ||
            x.publishedGlobalEndpoints != y.publishedGlobalEndpoints)
            return false;
    }
    for (std::size_t i = 0; i < a.localMappingParams.size(); ++i)
        if (a.localMappingParams[i].u != b.localMappingParams[i].u ||
            a.localMappingParams[i].v != b.localMappingParams[i].v)
            return false;
    for (std::size_t i = 0; i < a.localCells.size(); ++i)
        if (a.localCells[i].nodes != b.localCells[i].nodes ||
            a.localCells[i].classification.dim !=
                b.localCells[i].classification.dim ||
            a.localCells[i].classification.tag !=
                b.localCells[i].classification.tag)
            return false;
    return true;
}

const algo::SurfaceChain* requestChainById(
    const algo::SurfaceFaceRequest& request, algo::SurfaceChainId id) {
    for (const auto& chain : request.outer.pieces)
        if (chain.chainPersistentId == id) return &chain;
    for (const auto& hole : request.holes)
        for (const auto& chain : hole.pieces)
            if (chain.chainPersistentId == id) return &chain;
    for (const auto& chain : request.embedded)
        if (chain.chainPersistentId == id) return &chain;
    return nullptr;
}

bool originalProvenanceMatches(const algo::SurfaceFaceRequest& request,
                               const algo::SurfaceFaceResult& result) {
    for (const auto& [entity, origin] : result.segmentOrigins) {
        const auto* chain = requestChainById(request, origin.chainPersistentId);
        if (!chain ||
            origin.originalSegmentOrdinal >= chain->segmentEntityIds.size() ||
            chain->segmentEntityIds[origin.originalSegmentOrdinal] != entity)
            return false;
    }
    for (const auto& edge : result.persistentEdges) {
        const auto origin = result.segmentOrigins.find(edge.bndEntityId);
        if (origin == result.segmentOrigins.end() ||
            edge.originalSegmentOrdinal !=
                origin->second.originalSegmentOrdinal ||
            edge.normalizedSegmentOrdinal !=
                origin->second.normalizedSegmentOrdinal)
            return false;
    }
    return true;
}

struct Geometry : algo::SurfaceGeometryProvider {
    mutable int calls = 0;
    bnd::GeometryProviderStamp stamp() const override {
        ++calls; return {7, 1};
    }
    Vec3 pointAt(bnd::ChartId, const Vec2& uv) const override {
        ++calls; return {uv.u, uv.v, 0};
    }
    Vec3 normalAt(bnd::ChartId, const Vec2&) const override {
        ++calls; return {0, 0, 1};
    }
    algo::ChartCertificate chartCertificate(bnd::ChartId chart) const override {
        ++calls;
        return {{7, 1}, chart, {99, 1},
                {{-2, -2}, {4, -2}, {4, 4}, {-2, 4}}, 1, true, false, false};
    }
};

struct Candidate : algo::SurfaceCandidateProvider {
    mutable int calls = 0;
    algo::CandidateProviderStamp stamp() const override {
        ++calls; return {8, 1};
    }
    algo::SurfaceCandidate build(const algo::SurfaceFaceRequest&) const override {
        ++calls; return {};
    }
};

struct FullFrameCandidate : algo::SurfaceCandidateProvider {
    mutable int calls = 0;
    algo::CandidateProviderStamp stamp() const override {
        ++calls; return {8, 1};
    }
    algo::SurfaceCandidate build(const algo::SurfaceFaceRequest& request) const override {
        ++calls;
        algo::SurfaceCandidate c;
        c.providerStamp = {8, 1};
        c.faceStamp = request.faceStamp;
        c.chart = request.chart;
        const Vec2 uv[6]{{0, 0}, {2, 0}, {0, 2},
                         {-1, -1}, {4, -1}, {-1, 4}};
        for (std::size_t i = 0; i < 6; ++i) {
            algo::SurfaceCandidateNode n;
            n.node.coord = {uv[i].u, uv[i].v, 0};
            n.parameter = {uv[i], request.chart};
            n.frameOnly = i >= 3;
            if (i < 3) {
                n.node.cls = request.vertices[i].classification;
                c.declaredVertexMap.emplace(
                    request.vertices[i].faceLocalVertexId, i);
            } else {
                n.node.cls = request.faceClass;
            }
            c.nodes.push_back(n);
        }
        const std::array<std::array<std::size_t, 3>, 7> tris{{
            {{0, 1, 2}}, {{3, 4, 1}}, {{3, 1, 0}}, {{4, 5, 2}},
            {{4, 2, 1}}, {{5, 3, 0}}, {{5, 0, 2}}}};
        for (const auto& tri : tris)
            c.cells.push_back({tri, request.faceClass});
        c.frameBoundary = {3, 4, 5};
        c.exteriorSeedCell = 1;
        return c;
    }
};

std::uint64_t parameterBits(double value) {
    return std::bit_cast<std::uint64_t>(value == 0.0 ? 0.0 : value);
}

algo::SurfaceNodeIdentity faceInteriorIdentity(
    algo::ArtifactStamp faceStamp,const Vec2& uv) {
    algo::SurfaceNodeIdentity identity;
    identity.domain=algo::SurfaceNodeIdentityDomain::RefinedFacetInterior;
    identity.parentArtifact=faceStamp;
    identity.parentFacetGlobalNodes={0,1,2};
    identity.parameterUBits=parameterBits(uv.u);
    identity.parameterVBits=parameterBits(uv.v);
    return identity;
}

struct RefinedInteriorCandidate : FullFrameCandidate {
    algo::SurfaceCandidate build(
        const algo::SurfaceFaceRequest& request) const override {
        auto candidate=FullFrameCandidate::build(request);
        if(request.refinedNodes.size()!=1)return candidate;
        const auto& refined=request.refinedNodes.front();
        algo::SurfaceCandidateNode node;
        node.parameter=refined.parameter;
        node.node.coord={refined.parameter.p.u,refined.parameter.p.v,0.0};
        node.node.cls=refined.classification;
        node.persistentIdentity=refined.identity;
        const std::size_t middle=candidate.nodes.size();
        candidate.nodes.push_back(node);
        candidate.cells.erase(candidate.cells.begin());
        candidate.cells.insert(candidate.cells.begin(),{
            {{{0,1,middle}},request.faceClass},
            {{{1,2,middle}},request.faceClass},
            {{{2,0,middle}},request.faceClass}});
        candidate.exteriorSeedCell+=2;
        return candidate;
    }
};

struct WideGeometry : Geometry {
    algo::ChartCertificate chartCertificate(bnd::ChartId chart) const override {
        ++calls;
        return {{7, 1}, chart, {99, 1},
                {{-2, -2}, {8, -2}, {8, 6}, {-2, 6}}, 1, true, false, false};
    }
};

struct GridCandidate : algo::SurfaceCandidateProvider {
    mutable int calls = 0;
    algo::CandidateProviderStamp stamp() const override {
        ++calls; return {8, 1};
    }
    algo::SurfaceCandidate build(const algo::SurfaceFaceRequest& request) const override {
        ++calls;
        algo::SurfaceCandidate c;
        c.providerStamp = {8, 1};
        c.faceStamp = request.faceStamp;
        c.chart = request.chart;
        constexpr int nx = 9, ny = 7;
        auto at = [](int x, int y) {
            return static_cast<std::size_t>(y * nx + x);
        };
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x) {
                const Vec2 uv{static_cast<double>(x - 1),
                              static_cast<double>(y - 1)};
                algo::SurfaceCandidateNode n;
                n.node = {{uv.u, uv.v, 0}, request.faceClass};
                n.parameter = {uv, request.chart};
                n.frameOnly = x == 0 || x == nx - 1 || y == 0 || y == ny - 1;
                for (const auto& v : request.vertices)
                    if (v.parameter.p.u == uv.u && v.parameter.p.v == uv.v) {
                        n.node.cls = v.classification;
                        c.declaredVertexMap.emplace(v.faceLocalVertexId,
                                                    at(x, y));
                    }
                c.nodes.push_back(n);
            }
        for (int y = 0; y + 1 < ny; ++y)
            for (int x = 0; x + 1 < nx; ++x) {
                const auto a = at(x, y), b = at(x + 1, y);
                const auto d = at(x, y + 1), cc = at(x + 1, y + 1);
                c.cells.push_back({{a, b, cc}, request.faceClass});
                c.cells.push_back({{a, cc, d}, request.faceClass});
            }
        for (int x = 0; x < nx; ++x) c.frameBoundary.push_back(at(x, 0));
        for (int y = 1; y < ny; ++y) c.frameBoundary.push_back(at(nx - 1, y));
        for (int x = nx - 2; x >= 0; --x) c.frameBoundary.push_back(at(x, ny - 1));
        for (int y = ny - 2; y > 0; --y) c.frameBoundary.push_back(at(0, y));
        c.exteriorSeedCell = 0;
        return c;
    }
};

struct TwoChartGeometry : Geometry {
    Vec3 pointAt(bnd::ChartId chart, const Vec2& uv) const override {
        ++calls;
        return chart == 4 ? Vec3{2.0 - uv.u, -uv.v, 0.0}
                          : Vec3{uv.u, uv.v, 0.0};
    }
};

struct TwoChartCandidate : FullFrameCandidate {
    algo::SurfaceCandidate build(const algo::SurfaceFaceRequest& request) const override {
        auto c = FullFrameCandidate::build(request);
        if (request.chart == 4)
            for (auto& node : c.nodes)
                node.node.coord = {2.0 - node.parameter.p.u,
                                   -node.parameter.p.v, 0.0};
        return c;
    }
};

struct AdaptiveTriangleCandidate : FullFrameCandidate {
    algo::SurfaceCandidate build(const algo::SurfaceFaceRequest& request) const override {
        auto c = FullFrameCandidate::build(request);
        c.declaredVertexMap.clear();
        for (std::size_t i = 0; i < 3; ++i) {
            c.nodes[i].parameter = request.vertices[i].parameter;
            c.nodes[i].node.coord = {request.vertices[i].parameter.p.u,
                                     request.vertices[i].parameter.p.v, 0};
            c.nodes[i].node.cls = request.vertices[i].classification;
            c.declaredVertexMap.emplace(request.vertices[i].faceLocalVertexId,
                                        i);
        }
        const Vec2 frame[3]{{-1, -1}, {3, -1}, {1, 2}};
        for (std::size_t i = 0; i < 3; ++i) {
            c.nodes[3 + i].parameter.p = frame[i];
            c.nodes[3 + i].node.coord = {frame[i].u, frame[i].v, 0};
        }
        c.cells.clear();
        const std::array<std::array<std::size_t, 3>, 7> tris{{
            {{0, 1, 2}}, {{3, 4, 1}}, {{3, 1, 0}}, {{4, 5, 1}},
            {{5, 2, 1}}, {{5, 3, 0}}, {{5, 0, 2}}}};
        for (const auto& tri : tris)
            c.cells.push_back({tri, request.faceClass});
        c.frameBoundary = {3, 4, 5};
        c.exteriorSeedCell = 1;
        return c;
    }
};

struct NegativeFrameCandidate : FullFrameCandidate {
    algo::SurfaceCandidate build(const algo::SurfaceFaceRequest& request) const override {
        auto c = FullFrameCandidate::build(request);
        for (auto& cell : c.cells) std::swap(cell.nodes[1], cell.nodes[2]);
        std::reverse(c.frameBoundary.begin(), c.frameBoundary.end());
        return c;
    }
};

// Reviewer adversary: a provider inserts an otherwise valid, anonymous
// face-classified node exactly on a protected boundary segment.
struct AnonymousBoundaryCandidate : FullFrameCandidate {
    algo::SurfaceCandidate build(const algo::SurfaceFaceRequest& request) const override {
        auto c = FullFrameCandidate::build(request);
        algo::SurfaceCandidateNode midpoint;
        midpoint.node = {{1, 0, 0}, request.faceClass};
        midpoint.parameter = {{1, 0}, request.chart};
        const std::size_t m = c.nodes.size();
        c.nodes.push_back(midpoint);
        c.cells.erase(c.cells.begin() + 2);
        c.cells.erase(c.cells.begin());
        c.cells.push_back({{0, m, 2}, request.faceClass});
        c.cells.push_back({{m, 1, 2}, request.faceClass});
        c.cells.push_back({{3, 1, m}, request.faceClass});
        c.cells.push_back({{3, m, 0}, request.faceClass});
        c.exteriorSeedCell = 0;
        return c;
    }
};

enum class CandidateMutation {
    MalformedMap,
    MissingSeed,
    OverusedEdge,
    Pretrimmed,
    EmptyFrame
};
struct BadCandidate : FullFrameCandidate {
    CandidateMutation mutation;
    explicit BadCandidate(CandidateMutation m) : mutation(m) {}
    algo::SurfaceCandidate build(const algo::SurfaceFaceRequest& request) const override {
        auto c = FullFrameCandidate::build(request);
        switch (mutation) {
            case CandidateMutation::MalformedMap:
                c.declaredVertexMap.erase(c.declaredVertexMap.begin());
                break;
            case CandidateMutation::MissingSeed:
                c.exteriorSeedCell = c.cells.size();
                break;
            case CandidateMutation::OverusedEdge:
                c.cells.push_back({{1, 0, 4}, request.faceClass});
                break;
            case CandidateMutation::Pretrimmed:
                c.cells.erase(c.cells.begin());
                break;
            case CandidateMutation::EmptyFrame:
                c.frameBoundary.clear();
                for (auto& node : c.nodes) node.frameOnly = false;
                break;
        }
        return c;
    }
};

struct SecondFaceBadCandidate : TwoChartCandidate {
    algo::SurfaceCandidate build(const algo::SurfaceFaceRequest& request) const override {
        auto c = TwoChartCandidate::build(request);
        if (request.faceStamp.persistentId == 2)
            c.declaredVertexMap.erase(c.declaredVertexMap.begin());
        return c;
    }
};

struct ThrowingCandidate : algo::SurfaceCandidateProvider {
    bool allocation = false;
    algo::CandidateProviderStamp stamp() const override { return {8, 1}; }
    algo::SurfaceCandidate build(const algo::SurfaceFaceRequest&) const override {
        if (allocation) throw std::bad_alloc();
        throw std::runtime_error("candidate failure");
    }
};

struct ThrowingGeometry : Geometry {
    enum class Mode { Domain, Other } mode = Mode::Domain;
    Vec3 pointAt(bnd::ChartId, const Vec2&) const override {
        ++calls;
        if (mode == Mode::Domain)
            throw bnd::GeometryProviderDomainError("outside");
        throw std::runtime_error("geometry failure");
    }
};

struct CurvedGeometry : Geometry {
    double curvature = 0.1;
    bool reverseCertificate = false;
    bool badNormal = false;
    Vec3 pointAt(bnd::ChartId, const Vec2& uv) const override {
        ++calls; return {uv.u, uv.v, curvature * uv.u * uv.v};
    }
    Vec3 normalAt(bnd::ChartId, const Vec2& uv) const override {
        ++calls;
        if (badNormal) return {0, 0, 0};
        return {-curvature * uv.v, -curvature * uv.u, 1};
    }
    algo::ChartCertificate chartCertificate(bnd::ChartId chart) const override {
        ++calls;
        return {{7, 1}, chart, {99, 1},
                {{-2, -2}, {4, -2}, {4, 4}, {-2, 4}},
                reverseCertificate ? -1 : 1, true, false, false};
    }
};

struct VertexOnlyBadNormalGeometry : Geometry {
    Vec3 normalAt(bnd::ChartId, const Vec2& uv) const override {
        ++calls;
        return uv.u == 0.0 && uv.v == 0.0 ? Vec3{0, 0, 0}
                                           : Vec3{0, 0, 1};
    }
};

struct NonfiniteNormalGeometry : Geometry {
    Vec3 normalAt(bnd::ChartId, const Vec2&) const override {
        ++calls;
        return {0, 0, std::numeric_limits<double>::infinity()};
    }
};

struct PhaseChangingGeometry : Geometry {
    int changeAtCertificateCall = 0;
    mutable int certificateCalls = 0;
    algo::ChartCertificate chartCertificate(bnd::ChartId chart) const override {
        ++calls;
        ++certificateCalls;
        auto value = algo::ChartCertificate{
            {7, 1}, chart, {99, 1},
            {{-2, -2}, {4, -2}, {4, 4}, {-2, 4}},
            1, true, false, false};
        if (changeAtCertificateCall > 0 &&
            certificateCalls >= changeAtCertificateCall)
            value.certifiedRegionBoundary[2].u = 5.0;
        return value;
    }
};

algo::SurfaceSetRequest triangleRequest(mesh::Mesh& meshOut,
                                        mesh::Classification chainClass) {
    algo::SurfaceSetRequest set;
    set.setStamp = {9000, 1};
    algo::SurfaceFaceRequest face;
    face.faceStamp = {10, 1};
    face.chart = 3;
    face.faceClass = {2, 30};
    face.expectedGeometryStamp = {7, 1};
    face.expectedChartCertificateStamp = {99, 1};
    face.expectedCandidateProviderStamp = {8, 1};
    face.outer.role = algo::SurfaceChainRole::OuterBoundary;
    algo::SurfaceChain chain;
    chain.chainPersistentId = 50;
    chain.role = algo::SurfaceChainRole::OuterBoundary;
    chain.classification = chainClass;
    chain.closed = true;
    chain.owningCurve = algo::ArtifactStamp{500, 1};
    const Vec2 uv[3]{{0, 0}, {2, 0}, {0, 2}};
    for (std::uint64_t i = 0; i < 3; ++i) {
        const auto global = meshOut.addNode({uv[i].u, uv[i].v, 0},
                                            {0, static_cast<int>(100 + i)});
        algo::SharedCurveNodeKey key{500, 1, i};
        face.vertices.push_back({1000 + i, {uv[i], 3},
                                 meshOut.nodeClass(global)});
        chain.uses.push_back({1000 + i, key});
        chain.segmentEntityIds.push_back(2000 + i);
        set.sharedBindings.push_back(
            {key, global, meshOut.nodeClass(global)});
    }
    face.outer.pieces.push_back(std::move(chain));
    set.faces.push_back(std::move(face));
    return set;
}

algo::SurfaceSetRequest trianglePiecesRequest(mesh::Mesh& meshOut) {
    algo::SurfaceSetRequest set;
    set.setStamp = {9001, 1};
    algo::SurfaceFaceRequest face;
    face.faceStamp = {12, 1};
    face.chart = 3;
    face.faceClass = {2, 33};
    face.expectedGeometryStamp = {7, 1};
    face.expectedChartCertificateStamp = {99, 1};
    face.expectedCandidateProviderStamp = {8, 1};
    face.outer.role = algo::SurfaceChainRole::OuterBoundary;
    const Vec2 uv[3]{{0, 0}, {2, 0}, {0, 2}};
    paacs::Index globals[3]{};
    for (std::uint64_t i = 0; i < 3; ++i) {
        globals[i] = meshOut.addNode({uv[i].u, uv[i].v, 0},
                                     {0, static_cast<int>(120 + i)});
        face.vertices.push_back({1100 + i, {uv[i], 3},
                                 meshOut.nodeClass(globals[i])});
    }
    for (std::uint64_t i = 0; i < 3; ++i) {
        const auto j = (i + 1) % 3;
        algo::SurfaceChain chain;
        chain.chainPersistentId = 80 + i;
        chain.role = algo::SurfaceChainRole::OuterBoundary;
        chain.classification = {1, static_cast<int>(130 + i)};
        chain.closed = false;
        chain.owningCurve = algo::ArtifactStamp{800 + i, 1};
        const algo::SharedCurveNodeKey a{800 + i, 1, 0};
        const algo::SharedCurveNodeKey b{800 + i, 1, 1};
        chain.uses = {{1100 + i, a}, {1100 + j, b}};
        chain.segmentEntityIds = {2100 + i};
        face.outer.pieces.push_back(std::move(chain));
        set.sharedBindings.push_back(
            {a, globals[i], meshOut.nodeClass(globals[i])});
        set.sharedBindings.push_back(
            {b, globals[j], meshOut.nodeClass(globals[j])});
    }
    set.faces.push_back(std::move(face));
    return set;
}

algo::SurfaceSetRequest reversedSharedFacesRequest(mesh::Mesh& meshOut) {
    algo::SurfaceSetRequest set;
    set.setStamp = {9002, 1};
    const paacs::Index shared0 = meshOut.addNode({0, 0, 0}, {0, 140});
    const paacs::Index shared1 = meshOut.addNode({2, 0, 0}, {0, 141});
    const paacs::Index top = meshOut.addNode({0, 2, 0}, {0, 142});
    const paacs::Index bottom = meshOut.addNode({2, -2, 0}, {0, 143});
    const algo::SharedCurveNodeKey k0{900, 1, 0}, k1{900, 1, 1};
    set.sharedBindings.push_back({k0, shared0, meshOut.nodeClass(shared0)});
    set.sharedBindings.push_back({k1, shared1, meshOut.nodeClass(shared1)});
    auto addFace = [&](std::uint64_t faceId, bnd::ChartId chart,
                       paacs::Index thirdGlobal, bool reverseShared) {
        algo::SurfaceFaceRequest face;
        face.faceStamp = {faceId, 1}; face.chart = chart;
        face.faceClass = {2, static_cast<int>(40 + faceId)};
        face.expectedGeometryStamp = {7, 1};
        face.expectedChartCertificateStamp = {99, 1};
        face.expectedCandidateProviderStamp = {8, 1};
        face.outer.role = algo::SurfaceChainRole::OuterBoundary;
        const std::uint64_t base = faceId * 100;
        const paacs::Index localGlobal[3]{
            reverseShared ? shared1 : shared0,
            reverseShared ? shared0 : shared1, thirdGlobal};
        const Vec2 uv[3]{{0, 0}, {2, 0}, {0, 2}};
        for (std::uint64_t i = 0; i < 3; ++i)
            face.vertices.push_back(
                {base + i, {uv[i], chart}, meshOut.nodeClass(localGlobal[i])});
        for (std::uint64_t edge = 0; edge < 3; ++edge) {
            const auto next = (edge + 1) % 3;
            algo::SurfaceChain chain;
            chain.chainPersistentId = base + 10 + edge;
            chain.role = algo::SurfaceChainRole::OuterBoundary;
            chain.classification =
                edge == 0 ? mesh::Classification{1, 44}
                          : mesh::Classification{1, static_cast<int>(150 + base + edge)};
            chain.closed = false;
            const std::uint64_t curve = edge == 0 ? 900 : 10000 + base + edge;
            chain.owningCurve = algo::ArtifactStamp{curve, 1};
            algo::SharedCurveNodeKey a, b;
            if (edge == 0) {
                a = reverseShared ? k1 : k0;
                b = reverseShared ? k0 : k1;
            } else {
                a = {curve, 1, 0}; b = {curve, 1, 1};
                set.sharedBindings.push_back(
                    {a, localGlobal[edge], meshOut.nodeClass(localGlobal[edge])});
                set.sharedBindings.push_back(
                    {b, localGlobal[next], meshOut.nodeClass(localGlobal[next])});
            }
            chain.uses = {{base + edge, a}, {base + next, b}};
            chain.segmentEntityIds = {base + 20 + edge};
            face.outer.pieces.push_back(std::move(chain));
        }
        set.faces.push_back(std::move(face));
    };
    addFace(1, 3, top, false);
    addFace(2, 4, bottom, true);
    return set;
}

void addClosedLoop(algo::SurfaceSetRequest& set,
                   algo::SurfaceFaceRequest& face,
                   const std::vector<Vec2>& ring,
                   algo::SurfaceChainRole role,
                   std::uint64_t curveId,
                   std::uint64_t chainId,
                   int edgeTag,
                   std::uint64_t& nextVertex,
                   std::uint64_t& nextEntity,
                   mesh::Mesh& global) {
    algo::SurfaceChain chain;
    chain.chainPersistentId = chainId;
    chain.role = role;
    chain.classification = {1, edgeTag};
    chain.closed = true;
    chain.owningCurve = algo::ArtifactStamp{curveId, 1};
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto globalIndex = global.addNode(
            {ring[i].u, ring[i].v, 0},
            {0, static_cast<int>(nextVertex)});
        const algo::SharedCurveNodeKey key{curveId, 1, i};
        face.vertices.push_back(
            {nextVertex, {ring[i], face.chart}, global.nodeClass(globalIndex)});
        chain.uses.push_back({nextVertex, key});
        chain.segmentEntityIds.push_back(nextEntity++);
        set.sharedBindings.push_back(
            {key, globalIndex, global.nodeClass(globalIndex)});
        ++nextVertex;
    }
    if (role == algo::SurfaceChainRole::OuterBoundary) {
        face.outer.role = role;
        face.outer.pieces.push_back(std::move(chain));
    } else {
        algo::SurfaceLoop hole;
        hole.role = role;
        hole.pieces.push_back(std::move(chain));
        face.holes.push_back(std::move(hole));
    }
}

algo::SurfaceSetRequest twoHoleRequest(mesh::Mesh& global) {
    algo::SurfaceSetRequest set;
    set.setStamp = {9003, 1};
    algo::SurfaceFaceRequest face;
    face.faceStamp = {20, 1};
    face.chart = 3;
    face.faceClass = {2, 31};
    face.expectedGeometryStamp = {7, 1};
    face.expectedChartCertificateStamp = {99, 1};
    face.expectedCandidateProviderStamp = {8, 1};
    std::uint64_t vertex = 3000, entity = 4000;
    std::vector<Vec2> outer;
    for (int x = 0; x <= 6; ++x) outer.push_back({double(x), 0});
    for (int y = 1; y <= 4; ++y) outer.push_back({6, double(y)});
    for (int x = 5; x >= 0; --x) outer.push_back({double(x), 4});
    for (int y = 3; y > 0; --y) outer.push_back({0, double(y)});
    addClosedLoop(set, face, outer, algo::SurfaceChainRole::OuterBoundary,
                  600, 60, 80, vertex, entity, global);
    addClosedLoop(set, face, {{1, 1}, {2, 1}, {2, 2}, {1, 2}},
                  algo::SurfaceChainRole::HoleBoundary, 601, 61, 81,
                  vertex, entity, global);
    addClosedLoop(set, face, {{4, 1}, {5, 1}, {5, 2}, {4, 2}},
                  algo::SurfaceChainRole::HoleBoundary, 602, 62, 82,
                  vertex, entity, global);
    set.faces.push_back(std::move(face));
    return set;
}

void reorderClosedChain(algo::SurfaceChain& chain, std::size_t start,
                        bool reverse) {
    const auto sourceUses = chain.uses;
    const auto sourceIds = chain.segmentEntityIds;
    const std::size_t n = sourceUses.size();
    std::vector<algo::SurfaceChainVertex> reordered;
    reordered.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t index = reverse ? (start + n - i) % n
                                          : (start + i) % n;
        reordered.push_back(sourceUses[index]);
    }
    std::vector<bnd::ConstraintEntityId> ids;
    ids.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& a = reordered[i];
        const auto& b = reordered[(i + 1) % n];
        bool found = false;
        for (std::size_t j = 0; j < n; ++j) {
            const auto& x = sourceUses[j];
            const auto& y = sourceUses[(j + 1) % n];
            if ((sameVertexUse(a, x) && sameVertexUse(b, y)) ||
                (sameVertexUse(a, y) && sameVertexUse(b, x))) {
                ids.push_back(sourceIds[j]);
                found = true;
                break;
            }
        }
        if (!found) throw std::logic_error("closed-chain edge disappeared");
    }
    chain.uses = std::move(reordered);
    chain.segmentEntityIds = std::move(ids);
}

algo::SurfaceSetRequest polygonRequest(mesh::Mesh& global,
                                       const std::vector<Vec2>& corners,
                                       std::uint64_t faceId);

algo::SurfaceSetRequest twoEmbeddedRequest(mesh::Mesh& global,
                                           bool reverseInputOrder) {
    auto set = polygonRequest(global, {{0, 0}, {6, 0}, {6, 4}, {0, 4}},
                              56);
    auto& face = set.faces[0];
    const std::array<std::array<Vec2, 2>, 2> endpoints{{
        {{{2, 1}, {3, 1}}}, {{{4, 2}, {4, 3}}}}};
    for (std::size_t c = 0; c < endpoints.size(); ++c) {
        algo::SurfaceChain chain;
        chain.chainPersistentId = 9700 + c;
        chain.role = algo::SurfaceChainRole::EmbeddedConstraint;
        chain.classification = {1, static_cast<int>(97 + c)};
        for (std::size_t endpoint = 0; endpoint < 2; ++endpoint) {
            const auto id = 9710 + c * 2 + endpoint;
            face.vertices.push_back(
                {id, {endpoints[c][endpoint], face.chart},
                 {0, static_cast<int>(id)}});
            chain.uses.push_back({id, std::nullopt});
        }
        chain.segmentEntityIds = {9800 + c};
        face.embedded.push_back(std::move(chain));
    }
    if (reverseInputOrder)
        std::reverse(face.embedded.begin(), face.embedded.end());
    return set;
}

algo::SurfaceSetRequest embeddedRequest(mesh::Mesh& global) {
    algo::SurfaceSetRequest set;
    set.setStamp = {9004, 1};
    algo::SurfaceFaceRequest face;
    face.faceStamp = {30, 1};
    face.chart = 3;
    face.faceClass = {2, 32};
    face.expectedGeometryStamp = {7, 1};
    face.expectedChartCertificateStamp = {99, 1};
    face.expectedCandidateProviderStamp = {8, 1};
    std::uint64_t vertex = 5000, entity = 6000;
    std::vector<Vec2> outer;
    for (int x = 0; x <= 6; ++x) outer.push_back({double(x), 0});
    for (int y = 1; y <= 4; ++y) outer.push_back({6, double(y)});
    for (int x = 5; x >= 0; --x) outer.push_back({double(x), 4});
    for (int y = 3; y > 0; --y) outer.push_back({0, double(y)});
    addClosedLoop(set, face, outer, algo::SurfaceChainRole::OuterBoundary,
                  700, 70, 90, vertex, entity, global);
    algo::SurfaceChain embedded;
    embedded.chainPersistentId = 71;
    embedded.role = algo::SurfaceChainRole::EmbeddedConstraint;
    embedded.classification = {1, 91};
    embedded.closed = false;
    for (const Vec2 p : {Vec2{2, 2}, Vec2{3, 1}}) {
        face.vertices.push_back({vertex, {p, face.chart},
                                 {0, static_cast<int>(vertex)}});
        embedded.uses.push_back({vertex, std::nullopt});
        ++vertex;
    }
    embedded.segmentEntityIds = {entity++};
    face.embedded.push_back(std::move(embedded));
    set.faces.push_back(std::move(face));
    return set;
}

std::vector<Vec2> expandAxisRing(const std::vector<Vec2>& corners) {
    std::vector<Vec2> ring;
    for (std::size_t i = 0; i < corners.size(); ++i) {
        const Vec2 a = corners[i];
        const Vec2 b = corners[(i + 1) % corners.size()];
        int x = static_cast<int>(a.u), y = static_cast<int>(a.v);
        const int bx = static_cast<int>(b.u), by = static_cast<int>(b.v);
        const int dx = bx == x ? 0 : (bx > x ? 1 : -1);
        const int dy = by == y ? 0 : (by > y ? 1 : -1);
        if (ring.empty()) ring.push_back(a);
        while (x != bx || y != by) {
            x += dx; y += dy;
            if (x != static_cast<int>(corners[0].u) ||
                y != static_cast<int>(corners[0].v))
                ring.push_back({double(x), double(y)});
        }
    }
    return ring;
}

algo::SurfaceSetRequest polygonRequest(mesh::Mesh& global,
                                       const std::vector<Vec2>& corners,
                                       std::uint64_t faceId) {
    algo::SurfaceSetRequest set;
    set.setStamp = {10000 + faceId, 1};
    algo::SurfaceFaceRequest face;
    face.faceStamp = {faceId, 1}; face.chart = 3;
    face.faceClass = {2, static_cast<int>(200 + faceId)};
    face.expectedGeometryStamp = {7, 1};
    face.expectedChartCertificateStamp = {99, 1};
    face.expectedCandidateProviderStamp = {8, 1};
    std::uint64_t vertex = 7000 + faceId * 100;
    std::uint64_t entity = 8000 + faceId * 100;
    addClosedLoop(set, face, expandAxisRing(corners),
                  algo::SurfaceChainRole::OuterBoundary, 2000 + faceId,
                  3000 + faceId, static_cast<int>(220 + faceId), vertex,
                  entity, global);
    set.faces.push_back(std::move(face));
    return set;
}

algo::SurfaceSetRequest closedEmbeddedRequest(mesh::Mesh& global, int variant) {
    auto set = polygonRequest(global, {{0, 0}, {6, 0}, {6, 4}, {0, 4}},
                              55);
    auto& face = set.faces[0];
    const Vec2 points[4]{{2, 1}, {3, 1}, {3, 2}, {2, 2}};
    for (std::uint64_t i = 0; i < 4; ++i)
        face.vertices.push_back(
            {9500 + i, {points[i], face.chart},
             {0, static_cast<int>(9500 + i)}});
    const std::array<std::uint64_t, 4> order =
        variant == 0 ? std::array<std::uint64_t, 4>{0, 1, 2, 3}
        : variant == 1 ? std::array<std::uint64_t, 4>{2, 3, 0, 1}
                       : std::array<std::uint64_t, 4>{0, 3, 2, 1};
    const std::array<bnd::ConstraintEntityId, 4> edgeId{
        9600, 9601, 9602, 9603};
    algo::SurfaceChain chain;
    chain.chainPersistentId = 9590;
    chain.role = algo::SurfaceChainRole::EmbeddedConstraint;
    chain.classification = {1, 96};
    chain.closed = true;
    for (std::size_t i = 0; i < order.size(); ++i) {
        const auto a = order[i];
        const auto b = order[(i + 1) % order.size()];
        chain.uses.push_back({9500 + a, std::nullopt});
        bnd::ConstraintEntityId id = 0;
        for (std::size_t e = 0; e < 4; ++e)
            if ((e == a && (e + 1) % 4 == b) ||
                (e == b && (e + 1) % 4 == a))
                id = edgeId[e];
        chain.segmentEntityIds.push_back(id);
    }
    face.embedded.push_back(std::move(chain));
    return set;
}

void testClassificationPreflight() {
    SURFACE_TEST("surface_adapter.classification_preflight");
    for (const mesh::Classification bad :
         {mesh::Classification{0, 4}, mesh::Classification{2, 4},
          mesh::Classification{1, -1}}) {
        mesh::Mesh mesh;
        auto request = triangleRequest(mesh, bad);
        Geometry geometry;
        Candidate candidate;
        const auto result = algo::buildAndPublishSurfaceSet(
            request, mesh, geometry, candidate);
        SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Rejected);
        SURFACE_CHECK(result.setFailure.reason ==
                      algo::SurfaceReason::InvalidClassification);
        SURFACE_CHECK(result.reportComplete);
        SURFACE_CHECK(result.processed.size() == 1);
        SURFACE_CHECK(result.unattempted.empty());
        SURFACE_CHECK(geometry.calls == 0);
        SURFACE_CHECK(candidate.calls == 0);
    }
}

void testBindingAndIdentityPreflight() {
    SURFACE_TEST("surface_adapter.binding_identity_preflight");
    auto check = [](mesh::Mesh& mesh, algo::SurfaceSetRequest request,
                    algo::SurfaceReason reason) {
        Geometry geometry;
        Candidate candidate;
        const auto result = algo::buildAndPublishSurfaceSet(
            request, mesh, geometry, candidate);
        SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Rejected);
        SURFACE_CHECK(result.setFailure.reason == reason);
        SURFACE_CHECK(result.reportComplete);
        SURFACE_CHECK(result.unattempted.empty());
        SURFACE_CHECK(geometry.calls == 0);
        SURFACE_CHECK(candidate.calls == 0);
    };
    {
        mesh::Mesh mesh; auto r = triangleRequest(mesh, {1, 77});
        r.sharedBindings.push_back(r.sharedBindings.front());
        check(mesh, std::move(r), algo::SurfaceReason::DuplicateSharedBinding);
    }
    {
        mesh::Mesh mesh; auto r = triangleRequest(mesh, {1, 77});
        r.sharedBindings.pop_back();
        check(mesh, std::move(r), algo::SurfaceReason::MissingSharedBinding);
    }
    {
        mesh::Mesh mesh; auto r = triangleRequest(mesh, {1, 77});
        const auto extra = mesh.addNode({3, 3, 0}, {0, 199});
        r.sharedBindings.push_back({{999, 1, 0}, extra, mesh.nodeClass(extra)});
        check(mesh, std::move(r), algo::SurfaceReason::UnusedSharedBinding);
    }
    {
        mesh::Mesh mesh; auto r = triangleRequest(mesh, {1, 77});
        r.sharedBindings.front().existingGlobalIndex = paacs::kInvalidIndex;
        check(mesh, std::move(r), algo::SurfaceReason::UnknownSharedGlobalIndex);
    }
    {
        mesh::Mesh mesh; auto r = triangleRequest(mesh, {1, 77});
        r.sharedBindings.front().classification.tag += 1;
        check(mesh, std::move(r),
              algo::SurfaceReason::SharedClassificationConflict);
    }
    {
        mesh::Mesh mesh; auto r = triangleRequest(mesh, {1, 77});
        auto duplicate = r.faces.front();
        r.faces.push_back(std::move(duplicate));
        check(mesh, std::move(r), algo::SurfaceReason::DuplicateFaceIdentity);
    }
    {
        mesh::Mesh mesh; auto r = triangleRequest(mesh, {1, 77});
        r.faces[0].outer.pieces[0].segmentEntityIds.pop_back();
        check(mesh, std::move(r), algo::SurfaceReason::InvalidSegmentIdentity);
    }
    {
        mesh::Mesh mesh; auto r = triangleRequest(mesh, {1, 77});
        r.faces[0].outer.pieces[0].segmentEntityIds[1] =
            r.faces[0].outer.pieces[0].segmentEntityIds[0];
        check(mesh, std::move(r), algo::SurfaceReason::InvalidSegmentIdentity);
    }
    {
        mesh::Mesh mesh; auto r = triangleRequest(mesh, {1, 77});
        r.faces[0].vertices.push_back(r.faces[0].vertices.front());
        check(mesh, std::move(r),
              algo::SurfaceReason::InvalidConstraintEmbedding);
    }
}

void testCurveClassificationConflict() {
    SURFACE_TEST("surface_adapter.shared_curve_classification_conflict");
    mesh::Mesh mesh;
    auto request = triangleRequest(mesh, {1, 9});
    auto second = request.faces.front();
    second.faceStamp = {11, 1};
    second.outer.pieces.front().classification = {1, 10};
    request.faces.push_back(std::move(second));
    Geometry geometry;
    Candidate candidate;
    const auto result = algo::buildAndPublishSurfaceSet(
        request, mesh, geometry, candidate);
    SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Rejected);
    SURFACE_CHECK(result.setFailure.reason ==
                  algo::SurfaceReason::SharedClassificationConflict);
    SURFACE_CHECK(result.processed.size() == 2);
    SURFACE_CHECK(result.unattempted.empty());
    SURFACE_CHECK(geometry.calls == 0);
    SURFACE_CHECK(candidate.calls == 0);
}

void testOwnershipTraits() {
    SURFACE_TEST("surface_adapter.ownership_traits");
    static_assert(std::is_nothrow_move_constructible_v<algo::SurfaceFaceRow>);
    static_assert(std::is_nothrow_move_constructible_v<algo::SurfaceSetOutcome>);
    static_assert(std::is_nothrow_move_constructible_v<
                  std::unique_ptr<algo::SurfaceSetResult>>);
    SURFACE_CHECK(true);
}

void testRefinedIdentityContract() {
    SURFACE_TEST("surface_adapter.refined_identity_contract");
    const Vec2 uv{0.5,0.5};
    const auto identity=faceInteriorIdentity({10,1},uv);
    SURFACE_CHECK(algo::validSurfaceNodeIdentity(identity));
    {
        auto invalid=identity;
        invalid.parameterUBits=std::bit_cast<std::uint64_t>(-0.0);
        SURFACE_CHECK(!algo::validSurfaceNodeIdentity(invalid));
        invalid=identity;
        invalid.parentFacetGlobalNodes={0,2,1};
        SURFACE_CHECK(!algo::validSurfaceNodeIdentity(invalid));
        invalid=identity;
        invalid.parentOrdinal=1;
        SURFACE_CHECK(!algo::validSurfaceNodeIdentity(invalid));
    }
    algo::SurfaceNodeIdentity curveIdentity;
    curveIdentity.domain=algo::SurfaceNodeIdentityDomain::RefinedCurve;
    curveIdentity.parentArtifact={500,1};
    curveIdentity.parentOrdinal=0;
    curveIdentity.parameterUBits=parameterBits(0.5);
    curveIdentity.dyadicNumerator=1;
    curveIdentity.dyadicLevel=1;
    SURFACE_CHECK(algo::validSurfaceNodeIdentity(curveIdentity));

    {
        mesh::Mesh target;
        auto request=triangleRequest(target,{1,77});
        request.faces[0].refinedNodes.push_back(
            {{},{uv,request.faces[0].chart},request.faces[0].faceClass});
        const auto before=snapshotMesh(target);
        Geometry geometry;
        Candidate candidate;
        const auto outcome=algo::buildAndPublishSurfaceSet(
            request,target,geometry,candidate);
        SURFACE_CHECK(outcome.setStatus==algo::SurfaceSetStatus::Rejected);
        SURFACE_CHECK(outcome.setFailure.reason==
                      algo::SurfaceReason::SharedIdentityConflict);
        SURFACE_CHECK(geometry.calls==0);
        SURFACE_CHECK(candidate.calls==0);
        SURFACE_CHECK(snapshotMesh(target)==before);
    }
    {
        mesh::Mesh target;
        auto request=triangleRequest(target,{1,77});
        auto badKey=*request.faces[0].outer.pieces[0].uses[0].sharedKey;
        auto badIdentity=curveIdentity;
        badIdentity.parentArtifact={501,1};
        badKey.nodeOrdinal=0;
        badKey.identity=badIdentity;
        request.faces[0].outer.pieces[0].uses[0].sharedKey=badKey;
        request.sharedBindings[0].key=badKey;
        const auto before=snapshotMesh(target);
        Geometry geometry;
        Candidate candidate;
        const auto outcome=algo::buildAndPublishSurfaceSet(
            request,target,geometry,candidate);
        SURFACE_CHECK(outcome.setStatus==algo::SurfaceSetStatus::Rejected);
        SURFACE_CHECK(outcome.setFailure.reason==
                      algo::SurfaceReason::SharedIdentityConflict);
        SURFACE_CHECK(geometry.calls==0);
        SURFACE_CHECK(candidate.calls==0);
        SURFACE_CHECK(snapshotMesh(target)==before);
    }
    {
        mesh::Mesh target;
        auto request=triangleRequest(target,{1,77});
        request.faces[0].refinedNodes.push_back(
            {identity,{uv,request.faces[0].chart},request.faces[0].faceClass});
        Geometry geometry;
        RefinedInteriorCandidate candidate;
        const auto outcome=algo::buildAndPublishSurfaceSet(
            request,target,geometry,candidate);
        SURFACE_CHECK(outcome.setStatus==algo::SurfaceSetStatus::Published);
        SURFACE_CHECK(outcome.retainedArtifact!=nullptr);
        if(outcome.retainedArtifact) {
            const auto found=outcome.retainedArtifact->
                publishedGlobalByRefinedIdentity.find(identity);
            SURFACE_CHECK(found!=outcome.retainedArtifact->
                                      publishedGlobalByRefinedIdentity.end());
            if(found!=outcome.retainedArtifact->
                          publishedGlobalByRefinedIdentity.end()) {
                SURFACE_CHECK(found->second>=0&&found->second<target.nodeCount());
                SURFACE_CHECK(target.node(found->second).x==uv.u);
                SURFACE_CHECK(target.node(found->second).y==uv.v);
            }
            SURFACE_CHECK(outcome.retainedArtifact->faces.size()==1);
            if(outcome.retainedArtifact->faces.size()==1)
                SURFACE_CHECK(outcome.retainedArtifact->faces[0]->
                    refinedIdentityByLocalNode.size()==1);
        }
    }
}

void testDimensionZeroEndpointsAndPublication() {
    SURFACE_TEST("surface_adapter.dimension0_endpoints_publication");
    mesh::Mesh mesh;
    auto request = triangleRequest(mesh, {1, 77});
    const auto beforeCells = mesh.cellCount();
    Geometry geometry;
    FullFrameCandidate candidate;
    auto result = algo::buildAndPublishSurfaceSet(
        request, mesh, geometry, candidate);
    SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(result.reportComplete);
    SURFACE_CHECK(result.processed.size() == 1);
    SURFACE_CHECK(result.processed[0].status ==
                  algo::SurfaceFaceStatus::Ready);
    SURFACE_CHECK(result.unattempted.empty());
    SURFACE_CHECK(result.retainedArtifact != nullptr);
    SURFACE_CHECK(mesh.nodeCount() == 3);
    SURFACE_CHECK(mesh.cellCount() == beforeCells + 1);
    if (!result.retainedArtifact || result.retainedArtifact->faces.empty())
        return;
    const auto& face = *result.retainedArtifact->faces[0];
    SURFACE_CHECK(face.refinedIdentityByLocalNode.empty());
    SURFACE_CHECK(result.retainedArtifact->
                  publishedGlobalByRefinedIdentity.empty());
    SURFACE_CHECK(face.audit.passed());
    SURFACE_CHECK(face.boundaryLoops.size() == 1);
    SURFACE_CHECK(face.persistentEdges.size() == 3);
    SURFACE_CHECK(face.segmentOrigins.size() == 3);
    for (const auto& edge : face.persistentEdges) {
        SURFACE_CHECK(edge.classification.dim == 1);
        SURFACE_CHECK(edge.classification.tag == 77);
        SURFACE_CHECK(edge.publishedGlobalEndpoints.has_value());
    }
    for (const auto& [entity, origin] : face.segmentOrigins) {
        (void)entity;
        SURFACE_CHECK(origin.classification.dim == 1);
        SURFACE_CHECK(origin.classification.tag == 77);
    }
}

void testAdjacentChainClassifications() {
    SURFACE_TEST("surface_adapter.adjacent_chain_classifications");
    mesh::Mesh mesh;
    auto request = trianglePiecesRequest(mesh);
    Geometry geometry;
    FullFrameCandidate candidate;
    auto result = algo::buildAndPublishSurfaceSet(
        request, mesh, geometry, candidate);
    SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(result.retainedArtifact != nullptr);
    if (!result.retainedArtifact || result.retainedArtifact->faces.empty())
        return;
    const auto& face = *result.retainedArtifact->faces[0];
    SURFACE_CHECK(face.persistentEdges.size() == 3);
    std::set<int> tags;
    for (const auto& edge : face.persistentEdges) {
        tags.insert(edge.classification.tag);
        const auto& origin = face.segmentOrigins.at(edge.bndEntityId);
        SURFACE_CHECK(origin.classification.tag == edge.classification.tag);
        SURFACE_CHECK(origin.chainPersistentId == edge.chainPersistentId);
    }
    SURFACE_CHECK(tags == std::set<int>({130, 131, 132}));
}

void testReversedSharedBoundary() {
    SURFACE_TEST("surface_adapter.reversed_shared_boundary");
    mesh::Mesh mesh;
    auto request = reversedSharedFacesRequest(mesh);
    TwoChartGeometry geometry;
    TwoChartCandidate candidate;
    auto result = algo::buildAndPublishSurfaceSet(
        request, mesh, geometry, candidate);
    SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(result.processed.size() == 2);
    SURFACE_CHECK(result.retainedArtifact != nullptr);
    if (!result.retainedArtifact) return;
    SURFACE_CHECK(result.retainedArtifact->faces.size() == 2);
    for (const auto& face : result.retainedArtifact->faces) {
        SURFACE_CHECK(face->audit.passed());
        bool sawShared = false;
        for (const auto& edge : face->persistentEdges)
            if (edge.chainPersistentId ==
                face->faceStamp.persistentId * 100 + 10) {
                sawShared = true;
                SURFACE_CHECK(edge.classification.dim == 1);
                SURFACE_CHECK(edge.classification.tag == 44);
            }
        SURFACE_CHECK(sawShared);
    }
    SURFACE_CHECK(mesh.cellCount() == 2);

    // Whole-set input order is not authoritative. Reversing the two face
    // requests must preserve both the flat publication and every retained
    // existing-global mapping/provenance row.
    mesh::Mesh reorderedMesh;
    auto reorderedRequest = reversedSharedFacesRequest(reorderedMesh);
    std::reverse(reorderedRequest.faces.begin(), reorderedRequest.faces.end());
    TwoChartGeometry reorderedGeometry;
    TwoChartCandidate reorderedCandidate;
    auto reordered = algo::buildAndPublishSurfaceSet(
        reorderedRequest, reorderedMesh, reorderedGeometry,
        reorderedCandidate);
    SURFACE_CHECK(reordered.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(reordered.reportComplete);
    SURFACE_CHECK(reordered.processed.size() == result.processed.size());
    for (std::size_t i = 0;
         i < reordered.processed.size() && i < result.processed.size(); ++i)
        SURFACE_CHECK(reordered.processed[i].facePersistentId ==
                      result.processed[i].facePersistentId);

    SURFACE_CHECK(reorderedMesh.nodeCount() == mesh.nodeCount());
    for (paacs::Index i = 0;
         i < reorderedMesh.nodeCount() && i < mesh.nodeCount(); ++i) {
        const auto& a = mesh.node(i);
        const auto& b = reorderedMesh.node(i);
        SURFACE_CHECK(a.x == b.x && a.y == b.y && a.z == b.z);
        SURFACE_CHECK(mesh.nodeClass(i).dim == reorderedMesh.nodeClass(i).dim);
        SURFACE_CHECK(mesh.nodeClass(i).tag == reorderedMesh.nodeClass(i).tag);
    }
    SURFACE_CHECK(reorderedMesh.cellCount() == mesh.cellCount());
    for (paacs::Index i = 0;
         i < reorderedMesh.cellCount() && i < mesh.cellCount(); ++i) {
        const auto& a = mesh.cell(i);
        const auto& b = reorderedMesh.cell(i);
        SURFACE_CHECK(a.type == b.type);
        SURFACE_CHECK(a.onEntity.dim == b.onEntity.dim);
        SURFACE_CHECK(a.onEntity.tag == b.onEntity.tag);
        for (int j = 0; j < mesh::nodesPerCell(a.type); ++j)
            SURFACE_CHECK(mesh.cellNodes(i)[j] ==
                          reorderedMesh.cellNodes(i)[j]);
    }

    SURFACE_CHECK(reordered.retainedArtifact != nullptr);
    if (!reordered.retainedArtifact || !result.retainedArtifact) return;
    SURFACE_CHECK(reordered.retainedArtifact->faces.size() ==
                  result.retainedArtifact->faces.size());
    for (std::size_t i = 0;
         i < reordered.retainedArtifact->faces.size() &&
         i < result.retainedArtifact->faces.size(); ++i) {
        const auto& a = *result.retainedArtifact->faces[i];
        const auto& b = *reordered.retainedArtifact->faces[i];
        SURFACE_CHECK(a.faceStamp.persistentId == b.faceStamp.persistentId);
        SURFACE_CHECK(a.existingGlobalIndex == b.existingGlobalIndex);
        SURFACE_CHECK(a.boundaryLoops == b.boundaryLoops);
        SURFACE_CHECK(a.segmentOrigins.size() == b.segmentOrigins.size());
        SURFACE_CHECK(a.persistentEdges.size() == b.persistentEdges.size());
        for (std::size_t j = 0;
             j < a.persistentEdges.size() && j < b.persistentEdges.size();
             ++j) {
            SURFACE_CHECK(a.persistentEdges[j].bndEntityId ==
                          b.persistentEdges[j].bndEntityId);
            SURFACE_CHECK(a.persistentEdges[j].classification.dim ==
                          b.persistentEdges[j].classification.dim);
            SURFACE_CHECK(a.persistentEdges[j].classification.tag ==
                          b.persistentEdges[j].classification.tag);
            SURFACE_CHECK(a.persistentEdges[j].publishedGlobalEndpoints ==
                          b.persistentEdges[j].publishedGlobalEndpoints);
        }
    }
}

void testTwoHoleTopology() {
    SURFACE_TEST("surface_adapter.two_hole_topology");
    mesh::Mesh mesh;
    auto request = twoHoleRequest(mesh);
    const auto boundaryNodeCount = mesh.nodeCount();
    WideGeometry geometry;
    GridCandidate candidate;
    auto result = algo::buildAndPublishSurfaceSet(
        request, mesh, geometry, candidate);
    SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(result.retainedArtifact != nullptr);
    if (!result.retainedArtifact || result.retainedArtifact->faces.empty())
        return;
    const auto& face = *result.retainedArtifact->faces[0];
    SURFACE_CHECK(face.audit.passed());
    SURFACE_CHECK(face.audit.boundaryCycleCount == 3);
    SURFACE_CHECK(face.audit.eulerCharacteristic == -1);
    SURFACE_CHECK(face.boundaryLoops.size() == 3);
    SURFACE_CHECK(face.audit.noFrameSurvivors);
    SURFACE_CHECK(mesh.nodeCount() > boundaryNodeCount);
    SURFACE_CHECK(mesh.cellCount() > 0);
}

void testConcaveAndNarrowChannels() {
    SURFACE_TEST("surface_adapter.concave_narrow_channels");
    const std::vector<std::vector<Vec2>> fixtures{
        {{0, 0}, {6, 0}, {6, 1}, {2, 1}, {2, 4}, {0, 4}},
        {{0, 0}, {6, 0}, {6, 4}, {5, 4}, {5, 1}, {1, 1},
         {1, 4}, {0, 4}}};
    std::uint64_t faceId = 40;
    for (const auto& corners : fixtures) {
        mesh::Mesh mesh;
        auto request = polygonRequest(mesh, corners, faceId++);
        WideGeometry geometry;
        GridCandidate candidate;
        auto result = algo::buildAndPublishSurfaceSet(
            request, mesh, geometry, candidate);
        SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Published);
        SURFACE_CHECK(result.retainedArtifact != nullptr);
        if (!result.retainedArtifact) continue;
        const auto& face = *result.retainedArtifact->faces[0];
        SURFACE_CHECK(face.audit.passed());
        SURFACE_CHECK(face.audit.boundaryCycleCount == 1);
        SURFACE_CHECK(face.audit.eulerCharacteristic == 1);
        SURFACE_CHECK(face.audit.noFrameSurvivors);
    }
}

void testNearCollinearFeature() {
    SURFACE_TEST("surface_adapter.near_collinear_feature");
    mesh::Mesh mesh;
    auto request = triangleRequest(mesh, {1, 77});
    const Vec2 near{1.0, 1.0e-8};
    request.faces[0].vertices[2].parameter.p = near;
    mesh.setNode(2, {near.u, near.v, 0});
    Geometry geometry;
    AdaptiveTriangleCandidate candidate;
    auto result = algo::buildAndPublishSurfaceSet(
        request, mesh, geometry, candidate);
    SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(result.retainedArtifact != nullptr);
    if (result.retainedArtifact) {
        const auto& face = *result.retainedArtifact->faces[0];
        SURFACE_CHECK(face.audit.passed());
        SURFACE_CHECK(face.localCells.size() == 1);
        const auto& tri = face.localCells[0];
        SURFACE_CHECK(paacs::predicates::orient2d(
                          face.localMappingParams[tri.nodes[0]],
                          face.localMappingParams[tri.nodes[1]],
                          face.localMappingParams[tri.nodes[2]]) !=
                      paacs::predicates::Sign::Zero);
    }
}

void testClosedChainProvenance() {
    SURFACE_TEST("surface_adapter.closed_chain_provenance");
    std::map<bnd::ConstraintEntityId, std::size_t> canonicalOrdinals;
    for (int variant = 0; variant < 3; ++variant) {
        mesh::Mesh mesh;
        auto request = closedEmbeddedRequest(mesh, variant);
        WideGeometry geometry;
        GridCandidate candidate;
        auto result = algo::buildAndPublishSurfaceSet(
            request, mesh, geometry, candidate);
        SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Published);
        SURFACE_CHECK(result.retainedArtifact != nullptr);
        if (!result.retainedArtifact) continue;
        const auto& face = *result.retainedArtifact->faces[0];
        for (bnd::ConstraintEntityId id = 9600; id < 9604; ++id) {
            const auto& origin = face.segmentOrigins.at(id);
            SURFACE_CHECK(origin.classification.dim == 1);
            SURFACE_CHECK(origin.classification.tag == 96);
            const auto originalIds =
                request.faces[0].embedded[0].segmentEntityIds;
            const auto expectedOriginal = static_cast<std::size_t>(
                std::find(originalIds.begin(), originalIds.end(), id) -
                originalIds.begin());
            SURFACE_CHECK(origin.originalSegmentOrdinal == expectedOriginal);
            if (variant == 0)
                canonicalOrdinals[id] = origin.normalizedSegmentOrdinal;
            else
                SURFACE_CHECK(origin.normalizedSegmentOrdinal ==
                              canonicalOrdinals.at(id));
        }
        for (const auto& edge : face.persistentEdges)
            if (edge.chainPersistentId == 9590) {
                SURFACE_CHECK(edge.classification.tag == 96);
                const auto& origin = face.segmentOrigins.at(edge.bndEntityId);
                SURFACE_CHECK(edge.originalSegmentOrdinal ==
                              origin.originalSegmentOrdinal);
                SURFACE_CHECK(edge.normalizedSegmentOrdinal ==
                              origin.normalizedSegmentOrdinal);
            }
    }
}

void testFaceLocalEmbeddedClassification() {
    SURFACE_TEST("surface_adapter.face_local_embedded_classification");
    mesh::Mesh mesh;
    auto request = embeddedRequest(mesh);
    WideGeometry geometry;
    GridCandidate candidate;
    bnd::SteinerPolicy policy;
    policy.steinerEnabled = true;
    policy.splitBudget = 1;
    auto result = algo::buildAndPublishSurfaceSet(
        request, mesh, geometry, candidate, policy);
    if (result.setStatus != algo::SurfaceSetStatus::Published)
        std::printf("  embedded stop status=%d reason=%d bnd=(%d,%d,%d)\n",
                    static_cast<int>(result.setStatus),
                    static_cast<int>(result.setFailure.reason),
                    static_cast<int>(result.setFailure.bndReject),
                    static_cast<int>(result.setFailure.bndRefuse),
                    static_cast<int>(result.setFailure.bndFail));
    SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(result.retainedArtifact != nullptr);
    if (!result.retainedArtifact || result.retainedArtifact->faces.empty())
        return;
    const auto& face = *result.retainedArtifact->faces[0];
    SURFACE_CHECK(face.audit.passed());
    SURFACE_CHECK(face.localRecovery.processed.size() == 1);
    SURFACE_CHECK(face.localRecovery.unattempted.empty());
    std::size_t embeddedEdges = 0;
    for (const auto& edge : face.persistentEdges)
        if (edge.role == algo::SurfaceChainRole::EmbeddedConstraint) {
            ++embeddedEdges;
            SURFACE_CHECK(edge.classification.dim == 1);
            SURFACE_CHECK(edge.classification.tag == 91);
            SURFACE_CHECK(edge.publishedGlobalEndpoints.has_value());
        }
    SURFACE_CHECK(embeddedEdges >= 2);
    for (const auto& row : face.localRecovery.processed) {
        SURFACE_CHECK(row.result.status == bnd::SteinerStatus::Recovered ||
                      row.result.status == bnd::SteinerStatus::AlreadyPresent);
        SURFACE_CHECK(row.result.recoveredChain.size() >= 2);
    }
}

void testCandidateFailureDetails() {
    SURFACE_TEST("surface_adapter.candidate_failure_details");
    struct Case { CandidateMutation mutation; algo::SurfaceCandidateDetail detail; };
    const Case cases[] = {
        {CandidateMutation::MalformedMap,
         algo::SurfaceCandidateDetail::MalformedCandidate},
        {CandidateMutation::MissingSeed,
         algo::SurfaceCandidateDetail::MissingExteriorSeed},
        {CandidateMutation::OverusedEdge,
         algo::SurfaceCandidateDetail::DisconnectedFrame},
        {CandidateMutation::Pretrimmed,
         algo::SurfaceCandidateDetail::PretrimmedDomain},
        {CandidateMutation::EmptyFrame,
         algo::SurfaceCandidateDetail::MalformedCandidate},
    };
    for (const auto& c : cases) {
        mesh::Mesh mesh;
        auto request = triangleRequest(mesh, {1, 77});
        const auto beforeNodes = mesh.nodeCount();
        const auto beforeCells = mesh.cellCount();
        Geometry geometry;
        BadCandidate candidate(c.mutation);
        const auto result = algo::buildAndPublishSurfaceSet(
            request, mesh, geometry, candidate);
        SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Failed);
        SURFACE_CHECK(result.setFailure.reason ==
                      algo::SurfaceReason::CandidateValidationFailed);
        if (result.setFailure.candidate != c.detail)
            std::printf("  candidate detail got=%d expected=%d mutation=%d\n",
                        static_cast<int>(result.setFailure.candidate),
                        static_cast<int>(c.detail),
                        static_cast<int>(c.mutation));
        SURFACE_CHECK(result.setFailure.candidate == c.detail);
        SURFACE_CHECK(result.reportComplete);
        SURFACE_CHECK(result.processed.size() == 1);
        SURFACE_CHECK(result.unattempted.empty());
        SURFACE_CHECK(mesh.nodeCount() == beforeNodes);
        SURFACE_CHECK(mesh.cellCount() == beforeCells);
    }
    using S = paacs::predicates::Sign;
    SURFACE_CHECK(algo::testSurfaceBoundaryOwnerGate(
                      2, S::Positive, S::Negative) ==
                  algo::SurfaceDomainDetail::None);
    SURFACE_CHECK(algo::testSurfaceBoundaryOwnerGate(
                      1, S::Positive, S::Negative) ==
                  algo::SurfaceDomainDetail::OneSidedDomainBoundary);
    SURFACE_CHECK(algo::testSurfaceBoundaryOwnerGate(
                      2, S::Positive, S::Positive) ==
                  algo::SurfaceDomainDetail::OneSidedDomainBoundary);
    SURFACE_CHECK(algo::testSurfaceBoundaryOwnerGate(
                      2, S::Zero, S::Negative) ==
                  algo::SurfaceDomainDetail::OneSidedDomainBoundary);
}

void testProviderFailures() {
    SURFACE_TEST("surface_adapter.provider_failures");
    for (const bool alloc : {false, true}) {
        mesh::Mesh mesh;
        auto request = triangleRequest(mesh, {1, 77});
        Geometry geometry;
        ThrowingCandidate candidate;
        candidate.allocation = alloc;
        const auto result = algo::buildAndPublishSurfaceSet(
            request, mesh, geometry, candidate);
        SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Failed);
        SURFACE_CHECK(result.setFailure.reason ==
                      (alloc ? algo::SurfaceReason::AllocationFailed
                             : algo::SurfaceReason::CandidateProviderError));
        if (alloc)
            SURFACE_CHECK(result.setFailure.allocation ==
                          algo::SurfaceAllocationPhase::Candidate);
        SURFACE_CHECK(mesh.cellCount() == 0);
    }
    for (const auto mode : {ThrowingGeometry::Mode::Domain,
                            ThrowingGeometry::Mode::Other}) {
        mesh::Mesh mesh;
        auto request = triangleRequest(mesh, {1, 77});
        ThrowingGeometry geometry;
        geometry.mode = mode;
        FullFrameCandidate candidate;
        const auto result = algo::buildAndPublishSurfaceSet(
            request, mesh, geometry, candidate);
        SURFACE_CHECK(result.setFailure.reason ==
                      (mode == ThrowingGeometry::Mode::Domain
                           ? algo::SurfaceReason::ProviderDomainRejected
                           : algo::SurfaceReason::GeometryProviderError));
        SURFACE_CHECK(mesh.cellCount() == 0);
    }
}

void testSetPrefixAndPublicationValidation() {
    SURFACE_TEST("surface_adapter.set_prefix_publication_validation");
    {
        mesh::Mesh mesh;
        auto request = reversedSharedFacesRequest(mesh);
        const auto nodes = mesh.nodeCount();
        TwoChartGeometry geometry;
        SecondFaceBadCandidate candidate;
        auto result = algo::buildAndPublishSurfaceSet(
            request, mesh, geometry, candidate);
        SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Failed);
        SURFACE_CHECK(result.setFailure.reason ==
                      algo::SurfaceReason::CandidateValidationFailed);
        SURFACE_CHECK(result.processed.size() == 2);
        SURFACE_CHECK(result.processed[0].status ==
                      algo::SurfaceFaceStatus::Ready);
        SURFACE_CHECK(result.processed[0].readyResult != nullptr);
        SURFACE_CHECK(result.processed[1].status ==
                      algo::SurfaceFaceStatus::Failed);
        SURFACE_CHECK(result.unattempted.empty());
        SURFACE_CHECK(mesh.nodeCount() == nodes);
        SURFACE_CHECK(mesh.cellCount() == 0);
    }
    {
        mesh::Mesh builtMesh;
        auto request = triangleRequest(builtMesh, {1, 77});
        Geometry geometry;
        FullFrameCandidate candidate;
        auto built = algo::buildAndPublishSurfaceSet(
            request, builtMesh, geometry, candidate);
        SURFACE_CHECK(built.setStatus == algo::SurfaceSetStatus::Published);
        SURFACE_CHECK(built.retainedArtifact != nullptr);
        if (!built.retainedArtifact) return;
        mesh::Mesh target;
        target.addNode({0, 0, 0}, {0, 100});
        target.addNode({2, 0, 0}, {0, 101});
        target.addNode({0, 2, 0}, {0, 102});
        target.addCell(mesh::CellType::Tri3, {0, 1, 2}, {2, 30});
        const auto beforeNodes = target.nodeCount();
        const auto beforeCells = target.cellCount();
        bool published = false;
        const auto rejected = algo::publishSurfaceSet(
            target, *built.retainedArtifact, geometry, published);
        SURFACE_CHECK(!published);
        SURFACE_CHECK(rejected.reason ==
                      algo::SurfaceReason::PublicationValidationFailed);
        SURFACE_CHECK(target.nodeCount() == beforeNodes);
        SURFACE_CHECK(target.cellCount() == beforeCells);
        SURFACE_CHECK(built.retainedArtifact->faces[0]->audit.passed());
    }
}

void testCurvedPolicyAndOrientation() {
    SURFACE_TEST("surface_adapter.curved_policy_orientation");
    {
        mesh::Mesh mesh;
        auto request = triangleRequest(mesh, {1, 77});
        request.faces[0].physicalPolicy.maxMappingDeviation = 0.11;
        request.faces[0].physicalPolicy.minNormalAlignment = 0.9;
        CurvedGeometry geometry;
        FullFrameCandidate candidate;
        const auto result = algo::buildAndPublishSurfaceSet(
            request, mesh, geometry, candidate);
        SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Published);
        SURFACE_CHECK(result.retainedArtifact != nullptr);
    }
    for (int mode = 0; mode < 3; ++mode) {
        mesh::Mesh mesh;
        auto request = triangleRequest(mesh, {1, 77});
        request.faces[0].physicalPolicy.maxMappingDeviation =
            mode == 0 ? 0.05 : 0.11;
        request.faces[0].physicalPolicy.minNormalAlignment = 0.9;
        CurvedGeometry geometry;
        geometry.reverseCertificate = mode == 1;
        geometry.badNormal = mode == 2;
        FullFrameCandidate candidate;
        const auto result = algo::buildAndPublishSurfaceSet(
            request, mesh, geometry, candidate);
        SURFACE_CHECK(result.setStatus != algo::SurfaceSetStatus::Published);
        SURFACE_CHECK(result.setFailure.reason ==
                      (mode == 2 ? algo::SurfaceReason::NonfiniteMappedNormal
                                 : algo::SurfaceReason::PhysicalOrientationFailed));
        SURFACE_CHECK(mesh.cellCount() == 0);
    }
}

void testNegativeOrientationSuccess() {
    SURFACE_TEST("surface_adapter.negative_orientation_success");
    mesh::Mesh mesh;
    auto request = triangleRequest(mesh, {1, 77});
    auto& face = request.faces[0];
    face.desiredNormalSign = -1;
    auto& chain = face.outer.pieces[0];
    std::reverse(chain.uses.begin(), chain.uses.end());
    chain.segmentEntityIds = {2002, 2001, 2000};
    Geometry geometry;
    NegativeFrameCandidate candidate;
    auto result = algo::buildAndPublishSurfaceSet(
        request, mesh, geometry, candidate);
    SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(result.retainedArtifact != nullptr);
    if (result.retainedArtifact) {
        const auto& local = *result.retainedArtifact->faces[0];
        SURFACE_CHECK(local.audit.passed());
        const auto& tri = local.localCells[0];
        SURFACE_CHECK(paacs::predicates::orient2d(
                          local.localMappingParams[tri.nodes[0]],
                          local.localMappingParams[tri.nodes[1]],
                          local.localMappingParams[tri.nodes[2]]) ==
                      paacs::predicates::Sign::Negative);
    }
}

paacs::predicates::Sign retainedRingWinding(
    const algo::SurfaceFaceResult& face,
    const std::vector<std::size_t>& ring) {
    std::size_t best = 0;
    for (std::size_t i = 1; i < ring.size(); ++i) {
        const auto& p = face.localMappingParams[ring[i]];
        const auto& q = face.localMappingParams[ring[best]];
        if (p.u < q.u || (p.u == q.u && p.v < q.v)) best = i;
    }
    return paacs::predicates::orient2d(
        face.localMappingParams[ring[(best + ring.size() - 1) % ring.size()]],
        face.localMappingParams[ring[best]],
        face.localMappingParams[ring[(best + 1) % ring.size()]]);
}

void testHoleNormalizationDeterminism() {
    SURFACE_TEST("surface_adapter.hole_normalization_determinism");
    mesh::Mesh canonicalMesh;
    auto canonicalRequest = twoHoleRequest(canonicalMesh);
    WideGeometry canonicalGeometry;
    GridCandidate canonicalCandidate;
    auto canonical = algo::buildAndPublishSurfaceSet(
        canonicalRequest, canonicalMesh, canonicalGeometry,
        canonicalCandidate);

    mesh::Mesh permutedMesh;
    auto permutedRequest = twoHoleRequest(permutedMesh);
    std::reverse(permutedRequest.faces[0].holes.begin(),
                 permutedRequest.faces[0].holes.end());
    for (std::size_t i = 0; i < permutedRequest.faces[0].holes.size(); ++i)
        reorderClosedChain(
            permutedRequest.faces[0].holes[i].pieces[0], i + 1, true);
    reorderClosedChain(permutedRequest.faces[0].outer.pieces[0], 5, true);
    WideGeometry permutedGeometry;
    GridCandidate permutedCandidate;
    auto permuted = algo::buildAndPublishSurfaceSet(
        permutedRequest, permutedMesh, permutedGeometry, permutedCandidate);

    SURFACE_CHECK(canonical.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(permuted.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(snapshotMesh(canonicalMesh) == snapshotMesh(permutedMesh));
    SURFACE_CHECK(canonical.retainedArtifact != nullptr);
    SURFACE_CHECK(permuted.retainedArtifact != nullptr);
    if (!canonical.retainedArtifact || !permuted.retainedArtifact) return;
    const auto& a = *canonical.retainedArtifact->faces[0];
    const auto& b = *permuted.retainedArtifact->faces[0];
    SURFACE_CHECK(sameRetainedFace(a, b));
    SURFACE_CHECK(originalProvenanceMatches(canonicalRequest.faces[0], a));
    SURFACE_CHECK(originalProvenanceMatches(permutedRequest.faces[0], b));
    SURFACE_CHECK(a.boundaryLoops.size() == 3);
    if (a.boundaryLoops.size() == 3) {
        SURFACE_CHECK(retainedRingWinding(a, a.boundaryLoops[0]) ==
                      paacs::predicates::Sign::Positive);
        SURFACE_CHECK(retainedRingWinding(a, a.boundaryLoops[1]) ==
                      paacs::predicates::Sign::Negative);
        SURFACE_CHECK(retainedRingWinding(a, a.boundaryLoops[2]) ==
                      paacs::predicates::Sign::Negative);
    }
}

void testEmbeddedInputOrderDeterminism() {
    SURFACE_TEST("surface_adapter.embedded_input_order_determinism");
    mesh::Mesh forwardMesh;
    auto forwardRequest = twoEmbeddedRequest(forwardMesh, false);
    WideGeometry forwardGeometry;
    GridCandidate forwardCandidate;
    auto forward = algo::buildAndPublishSurfaceSet(
        forwardRequest, forwardMesh, forwardGeometry, forwardCandidate);

    mesh::Mesh reversedMesh;
    auto reversedRequest = twoEmbeddedRequest(reversedMesh, true);
    WideGeometry reversedGeometry;
    GridCandidate reversedCandidate;
    auto reversed = algo::buildAndPublishSurfaceSet(
        reversedRequest, reversedMesh, reversedGeometry, reversedCandidate);

    SURFACE_CHECK(forward.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(reversed.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(snapshotMesh(forwardMesh) == snapshotMesh(reversedMesh));
    SURFACE_CHECK(forward.retainedArtifact != nullptr);
    SURFACE_CHECK(reversed.retainedArtifact != nullptr);
    if (!forward.retainedArtifact || !reversed.retainedArtifact) return;
    SURFACE_CHECK(sameRetainedFace(*forward.retainedArtifact->faces[0],
                                   *reversed.retainedArtifact->faces[0]));
}

void testGeometryCertificatePhaseCheckpoints() {
    SURFACE_TEST("surface_adapter.geometry_certificate_phase_checkpoints");
    for (const int changeAt : {2, 3, 4}) {
        mesh::Mesh mesh;
        auto request = triangleRequest(mesh, {1, 77});
        const auto before = snapshotMesh(mesh);
        PhaseChangingGeometry geometry;
        geometry.changeAtCertificateCall = changeAt;
        FullFrameCandidate candidate;
        const auto result = algo::buildAndPublishSurfaceSet(
            request, mesh, geometry, candidate);
        SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Rejected);
        SURFACE_CHECK(result.setFailure.reason ==
                      algo::SurfaceReason::StaleGeometryProvider);
        SURFACE_CHECK(result.reportComplete);
        SURFACE_CHECK(result.processed.size() == 1);
        SURFACE_CHECK(result.unattempted.empty());
        SURFACE_CHECK(result.retainedArtifact == nullptr);
        SURFACE_CHECK(geometry.certificateCalls == changeAt);
        SURFACE_CHECK(snapshotMesh(mesh) == before);
    }
}

void testRepairRejectionInvariance() {
    SURFACE_TEST("surface_adapter.repair_rejection_invariance");
    {
        mesh::Mesh mesh;
        auto request = triangleRequest(mesh, {1, 77});
        const auto before = snapshotMesh(mesh);
        Geometry geometry;
        AnonymousBoundaryCandidate candidate;
        const auto result = algo::buildAndPublishSurfaceSet(
            request, mesh, geometry, candidate);
        SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Failed);
        SURFACE_CHECK(result.setFailure.reason ==
                      algo::SurfaceReason::CandidateValidationFailed);
        SURFACE_CHECK(result.setFailure.candidate ==
                      algo::SurfaceCandidateDetail::MalformedCandidate);
        SURFACE_CHECK(snapshotMesh(mesh) == before);
    }
    for (const bool nonfinite : {false, true}) {
        mesh::Mesh mesh;
        auto request = triangleRequest(mesh, {1, 77});
        const auto before = snapshotMesh(mesh);
        FullFrameCandidate candidate;
        algo::SurfaceSetOutcome result;
        if (nonfinite) {
            NonfiniteNormalGeometry geometry;
            result = algo::buildAndPublishSurfaceSet(
                request, mesh, geometry, candidate);
        } else {
            VertexOnlyBadNormalGeometry geometry;
            result = algo::buildAndPublishSurfaceSet(
                request, mesh, geometry, candidate);
        }
        SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Rejected);
        SURFACE_CHECK(result.setFailure.reason ==
                      algo::SurfaceReason::NonfiniteMappedNormal);
        SURFACE_CHECK(result.reportComplete);
        SURFACE_CHECK(result.processed.size() == 1);
        SURFACE_CHECK(result.unattempted.empty());
        SURFACE_CHECK(snapshotMesh(mesh) == before);
    }
}

void testReviewerCandidateConstraintOwnership() {
    SURFACE_TEST("surface_adapter.reviewer_candidate_constraint_ownership");
    mesh::Mesh mesh;
    auto request = triangleRequest(mesh, {1, 77});
    Geometry geometry;
    AnonymousBoundaryCandidate candidate;
    const auto result = algo::buildAndPublishSurfaceSet(
        request, mesh, geometry, candidate);
    std::printf("  reviewer ownership got status=%d reason=%d detail=%d\n",
                static_cast<int>(result.setStatus),
                static_cast<int>(result.setFailure.reason),
                static_cast<int>(result.setFailure.candidate));
    SURFACE_CHECK(result.setStatus == algo::SurfaceSetStatus::Failed);
    SURFACE_CHECK(result.setFailure.reason ==
                  algo::SurfaceReason::CandidateValidationFailed);
    SURFACE_CHECK(result.setFailure.candidate ==
                  algo::SurfaceCandidateDetail::MalformedCandidate);
    SURFACE_CHECK(mesh.cellCount() == 0);
}

void testReviewerRingReversalDeterminism() {
    SURFACE_TEST("surface_adapter.reviewer_ring_reversal_determinism");
    mesh::Mesh forwardMesh;
    auto forwardRequest = triangleRequest(forwardMesh, {1, 77});
    Geometry forwardGeometry;
    FullFrameCandidate forwardCandidate;
    auto forward = algo::buildAndPublishSurfaceSet(
        forwardRequest, forwardMesh, forwardGeometry, forwardCandidate);

    mesh::Mesh reverseMesh;
    auto reverseRequest = triangleRequest(reverseMesh, {1, 77});
    auto& reversed = reverseRequest.faces[0].outer.pieces[0];
    std::reverse(reversed.uses.begin(), reversed.uses.end());
    reversed.segmentEntityIds = {2001, 2000, 2002};
    Geometry reverseGeometry;
    FullFrameCandidate reverseCandidate;
    auto reverse = algo::buildAndPublishSurfaceSet(
        reverseRequest, reverseMesh, reverseGeometry, reverseCandidate);

    SURFACE_CHECK(forward.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(reverse.setStatus == algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(forward.retainedArtifact != nullptr);
    SURFACE_CHECK(reverse.retainedArtifact != nullptr);
    if (!forward.retainedArtifact || !reverse.retainedArtifact) return;
    const auto& a = *forward.retainedArtifact->faces[0];
    const auto& b = *reverse.retainedArtifact->faces[0];
    SURFACE_CHECK(a.boundaryLoops == b.boundaryLoops);
    SURFACE_CHECK(a.segmentOrigins.size() == b.segmentOrigins.size());
    SURFACE_CHECK(a.persistentEdges.size() == b.persistentEdges.size());
}

void testReviewerVertexNormalSampling() {
    SURFACE_TEST("surface_adapter.reviewer_vertex_normal_sampling");
    mesh::Mesh mesh;
    auto request = triangleRequest(mesh, {1, 77});
    VertexOnlyBadNormalGeometry geometry;
    FullFrameCandidate candidate;
    const auto result = algo::buildAndPublishSurfaceSet(
        request, mesh, geometry, candidate);
    SURFACE_CHECK(result.setStatus != algo::SurfaceSetStatus::Published);
    SURFACE_CHECK(result.setFailure.reason ==
                  algo::SurfaceReason::NonfiniteMappedNormal);
    SURFACE_CHECK(mesh.cellCount() == 0);
}

void testCandidateBoundaryPslgTaxonomy() {
    SURFACE_TEST("surface_adapter.candidate_boundary_pslg_taxonomy");
    mesh::Mesh global;
    auto set = triangleRequest(global,{1,77});
    Geometry geometry;
    FullFrameCandidate provider;
    auto candidate = provider.build(set.faces[0]);

    algo::SurfaceChain overlap;
    overlap.chainPersistentId = 90000;
    overlap.role = algo::SurfaceChainRole::EmbeddedConstraint;
    overlap.classification = {1,900};
    overlap.closed = false;
    overlap.owningCurve = algo::ArtifactStamp{90000,1};
    overlap.uses = {set.faces[0].outer.pieces[0].uses[0],
                    set.faces[0].outer.pieces[0].uses[1]};
    overlap.segmentEntityIds = {90001};
    set.faces[0].embedded.push_back(std::move(overlap));

    const auto result = algo::recoverCandidateBoundaryForTesting(
        set.faces[0],candidate,geometry);
    SURFACE_CHECK(!result.complete);
    SURFACE_CHECK(result.reject == bnd::SteinerReject::PslgRejected);
    SURFACE_CHECK(result.fail == bnd::SteinerFail::None);
}

}  // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* surfaceAdapterTests(std::size_t& count) {
    static const Entry tests[] = {
        {"surface_adapter.classification_preflight", testClassificationPreflight},
        {"surface_adapter.binding_identity_preflight",
         testBindingAndIdentityPreflight},
        {"surface_adapter.shared_curve_classification_conflict",
         testCurveClassificationConflict},
        {"surface_adapter.ownership_traits", testOwnershipTraits},
        {"surface_adapter.refined_identity_contract",
         testRefinedIdentityContract},
        {"surface_adapter.dimension0_endpoints_publication",
         testDimensionZeroEndpointsAndPublication},
        {"surface_adapter.adjacent_chain_classifications",
         testAdjacentChainClassifications},
        {"surface_adapter.reversed_shared_boundary",
         testReversedSharedBoundary},
        {"surface_adapter.two_hole_topology", testTwoHoleTopology},
        {"surface_adapter.concave_narrow_channels",
         testConcaveAndNarrowChannels},
        {"surface_adapter.near_collinear_feature",
         testNearCollinearFeature},
        {"surface_adapter.closed_chain_provenance",
         testClosedChainProvenance},
        {"surface_adapter.face_local_embedded_classification",
         testFaceLocalEmbeddedClassification},
        {"surface_adapter.candidate_failure_details",
         testCandidateFailureDetails},
        {"surface_adapter.provider_failures", testProviderFailures},
        {"surface_adapter.set_prefix_publication_validation",
         testSetPrefixAndPublicationValidation},
        {"surface_adapter.curved_policy_orientation",
         testCurvedPolicyAndOrientation},
        {"surface_adapter.negative_orientation_success",
         testNegativeOrientationSuccess},
        {"surface_adapter.hole_normalization_determinism",
         testHoleNormalizationDeterminism},
        {"surface_adapter.embedded_input_order_determinism",
         testEmbeddedInputOrderDeterminism},
        {"surface_adapter.geometry_certificate_phase_checkpoints",
         testGeometryCertificatePhaseCheckpoints},
        {"surface_adapter.repair_rejection_invariance",
         testRepairRejectionInvariance},
        {"surface_adapter.reviewer_candidate_constraint_ownership",
         testReviewerCandidateConstraintOwnership},
        {"surface_adapter.reviewer_ring_reversal_determinism",
         testReviewerRingReversalDeterminism},
        {"surface_adapter.reviewer_vertex_normal_sampling",
         testReviewerVertexNormalSampling},
        {"surface_adapter.candidate_boundary_pslg_taxonomy",
         testCandidateBoundaryPslgTaxonomy},
    };
    count = sizeof(tests) / sizeof(tests[0]);
    return tests;
}
}  // namespace paacs_tests
