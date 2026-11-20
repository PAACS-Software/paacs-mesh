// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#if defined(_MSC_VER)
// The accepted topology header's non-testing commitNodeBatch declaration has
// a test-only prediction parameter. Instantiating the qualified recovery
// types in this TU makes MSVC diagnose that frozen declaration even though
// CK-2 never calls it directly. Keep the suppression local to this TU; GCC
// and Clang still compile it under the full project warning set.
#pragma warning(push)
#pragma warning(disable : 4100)
#endif
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include "ProductionSurfaceMesherDetail.hpp"
#include "paacs/algo/ProtectedFeatureRefinement.hpp"

#include "CheckedSizing.hpp"
#include "PreparedSurfaceCandidate.hpp"
#include "T2CSurfaceCache.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>

#include "paacs/predicates/Predicates.hpp"
#include "paacs/algo/PlanarConstraintRecovery.hpp"

namespace paacs::algo {
namespace {

constexpr std::uint64_t kGeometryProviderId = 0x5041414353434b32ULL;
constexpr std::uint64_t kCandidateProviderId = 0x50414143434b3201ULL;
constexpr std::uint64_t kProviderRevision = 1;
constexpr std::uint64_t kDirectConstraintCandidateRevision = 2;

CandidateProviderStamp candidateStampFor(
        detail::SurfaceCandidateConstructionMode mode) noexcept {
    if (mode == detail::SurfaceCandidateConstructionMode::ScalablePlanarRecovery)
        return {kCandidateProviderId, 3};
    return {kCandidateProviderId,
            mode==detail::SurfaceCandidateConstructionMode::
                    ScalableDelaunayDirectConstraints
                ? kDirectConstraintCandidateRevision:kProviderRevision};
}

bool finite(const Vec2& p) noexcept {
    return std::isfinite(p.u) && std::isfinite(p.v);
}

bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) &&
           std::isfinite(p.z);
}

bool sameClass(const mesh::Classification& a,
               const mesh::Classification& b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}

std::uint64_t persistentId(int tag) noexcept {
    return static_cast<std::uint64_t>(static_cast<std::uint32_t>(tag)) + 1;
}

bnd::ConstraintEntityId segmentId(int edgeTag,
                                  std::size_t naturalOrdinal) {
    if (edgeTag < 0 || naturalOrdinal >=
                           static_cast<std::size_t>(
                               (std::numeric_limits<std::uint32_t>::max)()))
        throw std::invalid_argument("CK-2 segment identity is out of range");
    return (static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(edgeTag))
            << 32U) |
           (static_cast<std::uint64_t>(naturalOrdinal) + 1U);
}

bool pointOnSegment(const Vec2& a, const Vec2& b, const Vec2& p) noexcept {
    if (predicates::orient2d(a, b, p) != predicates::Sign::Zero) return false;
    return p.u >= std::min(a.u, b.u) && p.u <= std::max(a.u, b.u) &&
           p.v >= std::min(a.v, b.v) && p.v <= std::max(a.v, b.v);
}

struct MeshSnapshotShape {
    Index nodes = 0;
    Index cells = 0;
};

MeshSnapshotShape shapeOf(const mesh::Mesh& mesh) noexcept {
    return {mesh.nodeCount(), mesh.cellCount()};
}

struct EdgeDiscretization {
    const geo::GEdge* edge = nullptr;
    std::vector<Index> naturalNodes;
};

ProductionSurfaceReason extractEdgeDiscretization(
    const geo::GEdge& edge, const mesh::Mesh& destination,
    EdgeDiscretization& out) {
    std::map<Index, std::vector<Index>> graph;
    std::set<std::array<Index, 2>> segments;
    for (Index ci = 0; ci < destination.cellCount(); ++ci) {
        const auto& cell = destination.cell(ci);
        if (cell.type != mesh::CellType::Line2 || cell.onEntity.dim != 1 ||
            cell.onEntity.tag != edge.tag())
            continue;
        const Index* nodes = destination.cellNodes(ci);
        if (nodes[0] == nodes[1])
            return ProductionSurfaceReason::AmbiguousBoundaryDiscretization;
        const auto key = std::minmax(nodes[0], nodes[1]);
        const std::array<Index, 2> stored{key.first, key.second};
        if (!segments.insert(stored).second)
            return ProductionSurfaceReason::AmbiguousBoundaryDiscretization;
        graph[nodes[0]].push_back(nodes[1]);
        graph[nodes[1]].push_back(nodes[0]);
    }
    if (segments.empty())
        return ProductionSurfaceReason::MissingBoundaryDiscretization;
    for (auto& [node, adjacent] : graph) {
        (void)node;
        std::sort(adjacent.begin(), adjacent.end());
        if (adjacent.empty() || adjacent.size() > 2)
            return ProductionSurfaceReason::AmbiguousBoundaryDiscretization;
    }

    std::optional<Index> start;
    std::optional<Index> finish;
    for (const auto& [node, adjacent] : graph) {
        if (adjacent.size() != 1) continue;
        const auto& cls = destination.nodeClass(node);
        if (cls.dim == 0 && cls.tag == edge.startVertexTag()) start = node;
        if (cls.dim == 0 && cls.tag == edge.endVertexTag()) finish = node;
    }
    if (!start || !finish || *start == *finish)
        return ProductionSurfaceReason::AmbiguousBoundaryDiscretization;

    out.edge = &edge;
    out.naturalNodes.reserve(segments.size() + 1);
    Index previous = kInvalidIndex;
    Index current = *start;
    while (true) {
        out.naturalNodes.push_back(current);
        if (current == *finish) break;
        const auto it = graph.find(current);
        if (it == graph.end())
            return ProductionSurfaceReason::AmbiguousBoundaryDiscretization;
        Index next = kInvalidIndex;
        for (const Index candidate : it->second)
            if (candidate != previous) {
                if (next != kInvalidIndex)
                    return ProductionSurfaceReason::AmbiguousBoundaryDiscretization;
                next = candidate;
            }
        if (next == kInvalidIndex ||
            out.naturalNodes.size() > segments.size())
            return ProductionSurfaceReason::AmbiguousBoundaryDiscretization;
        previous = current;
        current = next;
    }
    if (out.naturalNodes.size() != segments.size() + 1)
        return ProductionSurfaceReason::AmbiguousBoundaryDiscretization;
    for (std::size_t i = 0; i < out.naturalNodes.size(); ++i) {
        const auto& cls = destination.nodeClass(out.naturalNodes[i]);
        if (i == 0) {
            if (cls.dim != 0 || cls.tag != edge.startVertexTag())
                return ProductionSurfaceReason::AmbiguousBoundaryDiscretization;
        } else if (i + 1 == out.naturalNodes.size()) {
            if (cls.dim != 0 || cls.tag != edge.endVertexTag())
                return ProductionSurfaceReason::AmbiguousBoundaryDiscretization;
        } else if (cls.dim != 1 || cls.tag != edge.tag()) {
            return ProductionSurfaceReason::AmbiguousBoundaryDiscretization;
        }
    }
    return ProductionSurfaceReason::None;
}

struct FaceContext {
    const geo::GFace* face = nullptr;
    const sizing::MeshSpec* spec = nullptr;
    Vec2 lo{};
    Vec2 hi{};
    Vec2 metricScale{};
    double referenceSize = 0.0;
    int chartPositiveToNormalSign = 0;
    std::map<FaceLocalVertexId, Vec3> declaredCoordinates;
    std::map<FaceLocalVertexId, SurfaceNodeIdentity>
        declaredPersistentIdentities;
    std::vector<t2c_detail::PrivateFacePoint> privateInteriorPoints;
    t2c_detail::SurfaceFaceCacheSession* cache = nullptr;
    FacePersistentId facePersistentId = 0;
};

double evaluatedFaceSize(const FaceContext& context,const Vec3& point,
                         double curvature) {
    const double value=detail::evaluatedGeometrySize(
        *context.spec,point,curvature);
    if(context.cache&&!context.cache->recordSizing(
           context.facePersistentId,point,curvature,value))
        throw std::bad_alloc{};
    return value;
}

bool computeChartSign(const geo::GFace& face, const Vec2& lo, const Vec2& hi,
                      int& sign) {
    const Vec2 center{(lo.u + hi.u) * 0.5, (lo.v + hi.v) * 0.5};
    const Vec3 pu0 = face.pointAt({lo.u, center.v});
    const Vec3 pu1 = face.pointAt({hi.u, center.v});
    const Vec3 pv0 = face.pointAt({center.u, lo.v});
    const Vec3 pv1 = face.pointAt({center.u, hi.v});
    const Vec3 normal = face.normalAt(center);
    if (!finite(pu0) || !finite(pu1) || !finite(pv0) || !finite(pv1) ||
        !finite(normal))
        return false;
    const Vec3 cross = (pu1 - pu0).cross(pv1 - pv0);
    const double d = cross.dot(normal);
    if (!std::isfinite(d) || d == 0.0 || cross.norm() == 0.0 ||
        normal.norm() == 0.0)
        return false;
    sign = d > 0.0 ? 1 : -1;
    return true;
}

struct GridKey {
    long long x = 0;
    long long y = 0;
    friend bool operator<(const GridKey& a, const GridKey& b) noexcept {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    }
};

class DeterministicSpacingIndex {
public:
    DeterministicSpacingIndex(double cellSize,
                              ProductionSurfaceFaceMetrics& metrics)
        : inv_(1.0 / cellSize), metrics_(metrics) {}

    void insert(const Vec2& p) {
        buckets_[key(p)].push_back(p);
        ++size_;
        metrics_.spacingIndexHighWater = std::max<std::uint64_t>(
            metrics_.spacingIndexHighWater, size_);
    }

    bool anyWithin(const Vec2& p, double radius) const {
        ++metrics_.spacingQueries;
        const double r2 = radius * radius;
        const long long span = static_cast<long long>(
            std::ceil(std::max(0.0, radius) * inv_));
        const GridKey center = key(p);
        for (long long y = center.y - span; y <= center.y + span; ++y)
            for (long long x = center.x - span; x <= center.x + span; ++x) {
                ++metrics_.spacingBucketsVisited;
                const auto it = buckets_.find({x, y});
                if (it == buckets_.end()) continue;
                for (const auto& q : it->second) {
                    ++metrics_.spacingDistanceTests;
                    const double du = q.u - p.u;
                    const double dv = q.v - p.v;
                    if (du * du + dv * dv < r2) return true;
                }
            }
        return false;
    }

private:
    long long coordinate(double x) const noexcept {
        const double c = std::floor(x * inv_);
        constexpr double limit = 4.0e18;
        return static_cast<long long>(std::clamp(c, -limit, limit));
    }
    GridKey key(const Vec2& p) const noexcept {
        return {coordinate(p.u), coordinate(p.v)};
    }

    double inv_;
    ProductionSurfaceFaceMetrics& metrics_;
    std::map<GridKey, std::vector<Vec2>> buckets_;
    std::uint64_t size_ = 0;
};

std::vector<std::pair<Vec2, Vec2>> requestSegments(
    const SurfaceFaceRequest& request) {
    std::map<FaceLocalVertexId, Vec2> params;
    for (const auto& vertex : request.vertices)
        params.emplace(vertex.faceLocalVertexId, vertex.parameter.p);
    std::vector<std::pair<Vec2, Vec2>> segments;
    auto add = [&](const SurfaceChain& chain) {
        const std::size_t count = chain.closed ? chain.uses.size()
                                               : chain.uses.size() - 1;
        for (std::size_t i = 0; i < count; ++i)
            segments.push_back({
                params.at(chain.uses[i].faceLocalVertexId),
                params.at(chain.uses[(i + 1) % chain.uses.size()]
                              .faceLocalVertexId)});
    };
    for (const auto& chain : request.outer.pieces) add(chain);
    for (const auto& hole : request.holes)
        for (const auto& chain : hole.pieces) add(chain);
    for (const auto& chain : request.embedded) add(chain);
    return segments;
}

std::vector<bnd::ConstraintSegment> lowerSegments(
    const SurfaceFaceRequest& request) {
    std::vector<bnd::ConstraintSegment> result;
    auto add = [&](const SurfaceChain& chain) {
        const std::size_t count = chain.closed ? chain.uses.size()
                                               : chain.uses.size() - 1;
        for (std::size_t i = 0; i < count; ++i)
            result.push_back({chain.segmentEntityIds[i],
                              chain.uses[i].faceLocalVertexId,
                              chain.uses[(i + 1) % chain.uses.size()]
                                  .faceLocalVertexId,
                              request.chart, chain.classification});
    };
    for (const auto& chain : request.outer.pieces) add(chain);
    for (const auto& hole : request.holes)
        for (const auto& chain : hole.pieces) add(chain);
    for (const auto& chain : request.embedded) add(chain);
    return result;
}

class KernelGeometryProvider final : public SurfaceGeometryProvider {
public:
    explicit KernelGeometryProvider(std::map<bnd::ChartId, FaceContext> faces)
        : faces_(std::move(faces)) {}

    bnd::GeometryProviderStamp stamp() const override {
        return {kGeometryProviderId, kProviderRevision};
    }

    Vec3 pointAt(bnd::ChartId chart, const Vec2& uv) const override {
        return context(chart).face->pointAt(uv);
    }

    Vec3 normalAt(bnd::ChartId chart, const Vec2& uv) const override {
        return context(chart).face->normalAt(uv);
    }

    ChartCertificate chartCertificate(bnd::ChartId chart) const override {
        const auto& c = context(chart);
        const double du = c.hi.u - c.lo.u;
        const double dv = c.hi.v - c.lo.v;
        const double mu = std::max(du * 0.5, c.referenceSize / c.metricScale.u);
        const double mv = std::max(dv * 0.5, c.referenceSize / c.metricScale.v);
        return {{kGeometryProviderId, kProviderRevision}, chart,
                {chart + 1U, kProviderRevision},
                {{c.lo.u - 2.0 * mu, c.lo.v - 2.0 * mv},
                 {c.hi.u + 2.0 * mu, c.lo.v - 2.0 * mv},
                 {c.hi.u + 2.0 * mu, c.hi.v + 2.0 * mv},
                 {c.lo.u - 2.0 * mu, c.hi.v + 2.0 * mv}},
                c.chartPositiveToNormalSign, true, false, false};
    }

private:
    const FaceContext& context(bnd::ChartId chart) const {
        const auto it = faces_.find(chart);
        if (it == faces_.end())
            throw bnd::GeometryProviderDomainError("unknown CK-2 chart");
        return it->second;
    }
    std::map<bnd::ChartId, FaceContext> faces_;
};

struct CandidateRecoveryFailure {
    bnd::SteinerReject reject = bnd::SteinerReject::None;
    bnd::SteinerRefuse refuse = bnd::SteinerRefuse::None;
    bnd::SteinerFail fail = bnd::SteinerFail::None;
};

class ScalableCandidateProvider final : public SurfaceCandidateProvider,
                                        public t2c_detail::SurfaceCacheMetricsSource {
public:
    ScalableCandidateProvider(const std::map<bnd::ChartId, FaceContext>& faces,
                              const KernelGeometryProvider& geometry,
                              detail::SurfaceCandidateConstructionMode mode =
                                  detail::SurfaceCandidateConstructionMode::
                                      ScalableDelaunay,
                              std::uint64_t maximumConstraintWork = 50000000)
        : faces_(faces), geometry_(geometry), mode_(mode),
          maximumConstraintWork_(maximumConstraintWork) {}

    CandidateProviderStamp stamp() const override {
        return candidateStampFor(mode_);
    }

    SurfaceCandidate build(const SurfaceFaceRequest& request) const override {
        const auto fit = faces_.find(request.chart);
        if (fit == faces_.end())
            throw std::invalid_argument("unknown CK-2 candidate chart");
        const auto& context = fit->second;

        ProductionSurfaceFaceMetrics metrics;
        metrics.facePersistentId = request.faceStamp.persistentId;
        metrics.boundaryNodes = request.vertices.size();
        metrics.protectedSegments = lowerSegments(request).size();

        SurfaceCandidate candidate;
        candidate.providerStamp = candidateStampFor(mode_);
        candidate.faceStamp = request.faceStamp;
        candidate.chart = request.chart;

        std::vector<Vec2> metricPoints;
        metricPoints.reserve(request.vertices.size() + 64);
        std::set<std::pair<double, double>> exactCoordinates;
        for (const auto& vertex : request.vertices) {
            const auto coordinate = context.declaredCoordinates.find(
                vertex.faceLocalVertexId);
            if (coordinate == context.declaredCoordinates.end())
                throw std::invalid_argument("missing CK-2 boundary coordinate");
            if (!exactCoordinates.emplace(vertex.parameter.p.u,
                                          vertex.parameter.p.v).second)
                throw std::invalid_argument("duplicate CK-2 candidate coordinate");
            const std::size_t index = candidate.nodes.size();
            const auto persistent=
                context.declaredPersistentIdentities.find(
                    vertex.faceLocalVertexId);
            candidate.nodes.push_back(
                {{coordinate->second, vertex.classification},
                 vertex.parameter, false,
                 persistent==context.declaredPersistentIdentities.end()
                    ? std::optional<SurfaceNodeIdentity>{}
                    : std::optional<SurfaceNodeIdentity>{persistent->second}});
            candidate.declaredVertexMap.emplace(vertex.faceLocalVertexId,
                                                index);
            const Vec2 metric{
                vertex.parameter.p.u * context.metricScale.u,
                vertex.parameter.p.v * context.metricScale.v};
            if (!finite(metric)) detail::failSizing(*context.spec);
            metricPoints.push_back(metric);
        }
        for (const auto& point : context.privateInteriorPoints) {
            const auto uv=point.parameter;
            const auto coordinate=point.coordinate;
            if (!exactCoordinates.emplace(uv.u,uv.v).second)
                throw std::invalid_argument(
                    "duplicate CK-2 private interior coordinate");
            candidate.nodes.push_back(
                {{coordinate,request.faceClass},{uv,request.chart},false,
                 point.identity});
            const Vec2 metric{uv.u*context.metricScale.u,
                              uv.v*context.metricScale.v};
            if (!finite(metric)) detail::failSizing(*context.spec);
            metricPoints.push_back(metric);
        }

        const double du = context.hi.u - context.lo.u;
        const double dv = context.hi.v - context.lo.v;
        const double mu = std::max(du * 0.25,
                                   context.referenceSize /
                                       context.metricScale.u);
        const double mv = std::max(dv * 0.25,
                                   context.referenceSize /
                                       context.metricScale.v);
        const std::array<Vec2, 4> frame{{
            {context.lo.u - mu, context.lo.v - mv},
            {context.hi.u + mu, context.lo.v - mv},
            {context.hi.u + mu, context.hi.v + mv},
            {context.lo.u - mu, context.hi.v + mv}}};
        for (const auto& uv : frame) {
            if (!exactCoordinates.emplace(uv.u, uv.v).second)
                throw std::invalid_argument("duplicate CK-2 frame coordinate");
            candidate.frameBoundary.push_back(candidate.nodes.size());
            candidate.nodes.push_back(
                {{{uv.u, uv.v, 0.0}, request.faceClass},
                 {uv, request.chart}, true, std::nullopt});
            candidate.nodes.back().node.coord =
                geometry_.pointAt(request.chart, uv);
            const Vec2 metric{uv.u * context.metricScale.u,
                              uv.v * context.metricScale.v};
            if (!finite(metric)) detail::failSizing(*context.spec);
            metricPoints.push_back(metric);
        }

        const double spacingCell = context.referenceSize * 0.6;
        if (!detail::positiveFinite(spacingCell) ||
            !detail::positiveFinite(1.0 / spacingCell))
            detail::failSizing(*context.spec);
        DeterministicSpacingIndex spacing(spacingCell, metrics);
        for (const auto& p : metricPoints) spacing.insert(p);
        const auto constraints = requestSegments(request);
        const Vec2 mlo{context.lo.u * context.metricScale.u,
                       context.lo.v * context.metricScale.v};
        const Vec2 mhi{context.hi.u * context.metricScale.u,
                       context.hi.v * context.metricScale.v};
        if (!finite(mlo) || !finite(mhi))
            detail::failSizing(*context.spec);
        const double rowStep =
            context.referenceSize * 0.86602540378443864676;
        if (!detail::positiveFinite(rowStep))
            detail::failSizing(*context.spec);
        const double firstRow = mlo.v + context.referenceSize;
        const double rowStop = mhi.v - 0.4 * context.referenceSize;
        if (!std::isfinite(firstRow) || !std::isfinite(rowStop))
            detail::failSizing(*context.spec);
        if (firstRow < rowStop) {
            detail::requireProgress(*context.spec, firstRow, rowStep);
            const double rowCount =
                std::ceil((rowStop - firstRow) / rowStep);
            if (!detail::positiveFinite(rowCount) ||
                rowCount > static_cast<double>(
                               (std::numeric_limits<std::uint64_t>::max)()))
                detail::failSizing(*context.spec);
        }
        std::uint64_t row = 0;
        for (double v = firstRow; v < rowStop;
             detail::requireProgress(*context.spec, v, rowStep),
                    v += rowStep, ++row) {
            const double shift = (row & 1U) != 0U
                                     ? context.referenceSize * 0.5
                                     : 0.0;
            const double firstColumn =
                mlo.u + context.referenceSize * 0.7 + shift;
            const double columnStop =
                mhi.u - 0.4 * context.referenceSize;
            if (!std::isfinite(firstColumn) || !std::isfinite(columnStop))
                detail::failSizing(*context.spec);
            if (firstColumn < columnStop) {
                detail::requireProgress(*context.spec, firstColumn,
                                        context.referenceSize);
                const double columnCount = std::ceil(
                    (columnStop - firstColumn) / context.referenceSize);
                if (!detail::positiveFinite(columnCount) ||
                    columnCount > static_cast<double>(
                                      (std::numeric_limits<std::size_t>::max)()))
                    detail::failSizing(*context.spec);
            }
            for (double u = firstColumn; u < columnStop;
                 detail::requireProgress(*context.spec, u,
                                         context.referenceSize),
                        u += context.referenceSize) {
                const Vec2 uv{u / context.metricScale.u,
                              v / context.metricScale.v};
                if (!finite(uv)) detail::failSizing(*context.spec);
                bool onConstraint = false;
                for (const auto& segment : constraints) {
                    ++metrics.constraintPointTests;
                    if (pointOnSegment(segment.first, segment.second, uv)) {
                        onConstraint = true;
                        break;
                    }
                }
                if (onConstraint) continue;
                const Vec2 metric{u, v};
                const Vec3 point = geometry_.pointAt(request.chart, uv);
                const double localSize=evaluatedFaceSize(
                    context,point,context.face->curvatureAt(uv));
                const double radius = 0.6 * localSize;
                if (!detail::positiveFinite(radius) ||
                    !detail::positiveFinite(radius * radius) ||
                    !std::isfinite(std::ceil(radius / spacingCell)) ||
                    std::ceil(radius / spacingCell) >
                        static_cast<double>(
                            (std::numeric_limits<long long>::max)() / 4))
                    detail::failSizing(*context.spec);
                if (!spacing.anyWithin(metric, radius) &&
                    exactCoordinates.emplace(uv.u, uv.v).second) {
                    candidate.nodes.push_back(
                        {{point, request.faceClass}, {uv, request.chart}, false,
                         std::nullopt});
                    metricPoints.push_back(metric);
                    spacing.insert(metric);
                    ++metrics.interiorCandidates;
                }
            }
        }

        auto triangles = delaunay2DScalable(
            metricPoints, &metrics.construction);
        if (triangles.empty())
            throw std::invalid_argument("empty CK-2 candidate triangulation");
        if (mode_ == detail::SurfaceCandidateConstructionMode::ScalablePlanarRecovery) {
            std::vector<PlanarSegment> segments;
            for (const auto& s : lowerSegments(request))
                segments.push_back({static_cast<Index>(candidate.declaredVertexMap.at(s.a)),
                                    static_cast<Index>(candidate.declaredVertexMap.at(s.b))});
            const auto recovered = recoverPlanarConstraints(
                metricPoints, segments, triangles, maximumConstraintWork_);
            planarRecovery_ = recovered;
            if (!recovered.completed) {
                if (recovered.reason == PlanarRecoveryReason::AllocationFailed)
                    throw std::bad_alloc{};
                throw std::runtime_error("bounded planar constraint recovery stopped");
            }
        }
        const int requiredWinding = context.chartPositiveToNormalSign *
                                    request.desiredNormalSign;
        candidate.cells.reserve(triangles.size());
        for (const auto& triangle : triangles) {
            std::array<std::size_t, 3> nodes{
                static_cast<std::size_t>(triangle.v[0]),
                static_cast<std::size_t>(triangle.v[1]),
                static_cast<std::size_t>(triangle.v[2])};
            if (requiredWinding < 0) std::swap(nodes[1], nodes[2]);
            candidate.cells.push_back({nodes, request.faceClass});
        }

        recoverBoundary(request, candidate, metrics);
        const auto frameEdges = [&]() {
            std::set<std::array<std::size_t, 2>> edges;
            for (std::size_t i = 0; i < candidate.frameBoundary.size(); ++i) {
                auto a = candidate.frameBoundary[i];
                auto b = candidate.frameBoundary[
                    (i + 1) % candidate.frameBoundary.size()];
                if (b < a) std::swap(a, b);
                edges.insert({a, b});
            }
            return edges;
        }();
        candidate.exteriorSeedCell = candidate.cells.size();
        for (std::size_t i = 0; i < candidate.cells.size(); ++i) {
            const auto& c = candidate.cells[i].nodes;
            for (std::size_t e = 0; e < 3; ++e) {
                auto a = c[e];
                auto b = c[(e + 1) % 3];
                if (b < a) std::swap(a, b);
                if (frameEdges.count({a, b}) != 0) {
                    candidate.exteriorSeedCell = i;
                    break;
                }
            }
            if (candidate.exteriorSeedCell != candidate.cells.size()) break;
        }
        if (candidate.exteriorSeedCell == candidate.cells.size())
            throw std::invalid_argument("missing CK-2 exterior seed");
        if (candidate.preparedRecovery) {
            candidate.preparedRecovery->frameBoundary =
                candidate.frameBoundary;
            candidate.preparedRecovery->exteriorSeedCell =
                candidate.exteriorSeedCell;
            candidate.preparedRecovery->finalized = true;
        }

        rows_.push_back(metrics);
        return candidate;
    }

    std::vector<ProductionSurfaceFaceMetrics> takeRows() noexcept {
        return std::move(rows_);
    }

    bool metricForFace(FacePersistentId face,
                       ProductionSurfaceFaceMetrics& out) const override {
        for(const auto& row:rows_)
            if(row.facePersistentId==face) { out=row; return true; }
        return false;
    }

    const std::optional<CandidateRecoveryFailure>& recoveryFailure() const
        noexcept {
        return recoveryFailure_;
    }
    PlanarRecoveryOutcome planarRecovery() const noexcept { return planarRecovery_; }

private:
    void recoverBoundary(const SurfaceFaceRequest& request,
                         SurfaceCandidate& candidate,
                         ProductionSurfaceFaceMetrics& metrics) const {
        const auto recovered = detail::recoverCandidateBoundaryForCk2(
            request, candidate, geometry_,
            mode_ != detail::SurfaceCandidateConstructionMode::ScalableDelaunay);
        metrics.recoveredProtectedSegments = recovered.recoveredSegments;
        if (!recovered.complete) {
            recoveryFailure_ = CandidateRecoveryFailure{
                recovered.reject, recovered.refuse, recovered.fail};
            if (recovered.fail == bnd::SteinerFail::CandidateAllocFailed ||
                recovered.fail == bnd::SteinerFail::StagedNodeAllocFailed ||
                recovered.fail == bnd::SteinerFail::CavityThrew ||
                recovered.fail == bnd::SteinerFail::IndexRebuildFailed)
                throw std::bad_alloc{};
            throw std::runtime_error("CK-2 boundary recovery stopped");
        }
    }

    const std::map<bnd::ChartId, FaceContext>& faces_;
    const KernelGeometryProvider& geometry_;
    detail::SurfaceCandidateConstructionMode mode_ =
        detail::SurfaceCandidateConstructionMode::ScalableDelaunay;
    std::uint64_t maximumConstraintWork_ = 50000000;
    mutable PlanarRecoveryOutcome planarRecovery_;
    mutable std::vector<ProductionSurfaceFaceMetrics> rows_;
    mutable std::optional<CandidateRecoveryFailure> recoveryFailure_;
};

class CheckedSizingCacheEvaluator final
    : public t2c_detail::SurfaceCacheSizingEvaluator {
public:
    explicit CheckedSizingCacheEvaluator(const sizing::MeshSpec& spec)
        : spec_(spec) {}
    double evaluate(const Vec3& point,double curvature) const override {
        return detail::evaluatedGeometrySize(spec_,point,curvature);
    }
private:
    const sizing::MeshSpec& spec_;
};

struct RequestBuild {
    SurfaceSetRequest request;
    std::map<bnd::ChartId, FaceContext> contexts;
    ProductionSurfaceReason reason = ProductionSurfaceReason::None;
};

RequestBuild buildRequest(const geo::Model& model,
                          const sizing::MeshSpec& spec,
                          const mesh::Mesh& destination,
                          std::span<const t2c_detail::PrivateCurveNodeIdentity>
                              curveNodeIdentities,
                          std::span<const t2c_detail::PrivateFacePoint>
                              facePoints,
                          t2c_detail::SurfaceFaceCacheSession* cache=nullptr,
                          std::span<const t2c_detail::PrivateFacePoint>
                              proposedFacePoints={},
                          detail::SurfaceCandidateConstructionMode mode=
                              detail::SurfaceCandidateConstructionMode::
                                  ScalableDelaunay) {
    RequestBuild built;
    // Constrained-v2 is the internal caller of the accepted surface-set API.
    // Its fixed authority is compatibility metadata only; boundary-
    // conforming-v3 supplies its caller-owned stamp through the v3 request.
    built.request.setStamp = {0x434b325355524641ULL, 1};
    std::vector<const geo::GFace*> faces = model.faces();
    if (faces.empty()) {
        built.reason = ProductionSurfaceReason::EmptyFaceSet;
        return built;
    }
    std::sort(faces.begin(), faces.end(), [](const auto* a, const auto* b) {
        return a->tag() < b->tag();
    });
    std::set<int> faceTags;
    std::map<int, EdgeDiscretization> edgeData;
    std::map<SharedCurveNodeKey, SharedCurveNodeBinding> bindings;
    std::map<Index,SurfaceNodeIdentity> refinedCurveIdentity;
    for(const auto& row:curveNodeIdentities) {
        if(row.globalNode<0||row.globalNode>=destination.nodeCount()||
           row.identity.domain!=SurfaceNodeIdentityDomain::RefinedCurve||
           !refinedCurveIdentity.emplace(row.globalNode,row.identity).second) {
            built.reason=ProductionSurfaceReason::AmbiguousBoundaryDiscretization;
            return built;
        }
    }

    for (const auto* face : faces) {
        if (!face || face->tag() < 0 || !faceTags.insert(face->tag()).second) {
            built.reason = ProductionSurfaceReason::InvalidFaceIdentity;
            return built;
        }
        if (face->periodicU() || face->periodicV()) {
            built.reason = ProductionSurfaceReason::UnsupportedPeriodicChart;
            return built;
        }
        const auto loops = face->boundaryLoops();
        std::vector<int> edgeTags;
        for (const auto& loop : loops) {
            if (loop.empty()) {
                built.reason = ProductionSurfaceReason::InvalidBoundaryContract;
                return built;
            }
            edgeTags.insert(edgeTags.end(), loop.begin(), loop.end());
        }
        if (edgeTags.empty()) {
            built.reason = ProductionSurfaceReason::InvalidBoundaryContract;
            return built;
        }
        for (const int edgeTag : edgeTags) {
            if (edgeTag < 0) {
                built.reason = ProductionSurfaceReason::InvalidBoundaryContract;
                return built;
            }
            if (edgeData.count(edgeTag) != 0) continue;
            const auto* edge = model.edge(edgeTag);
            if (!edge) {
                built.reason = ProductionSurfaceReason::InvalidBoundaryContract;
                return built;
            }
            EdgeDiscretization data;
            built.reason = extractEdgeDiscretization(*edge, destination, data);
            if (built.reason != ProductionSurfaceReason::None) return built;
            edgeData.emplace(edgeTag, std::move(data));
        }
    }

    for (const auto* face : faces) {
        SurfaceFaceRequest request;
        request.faceStamp = {persistentId(face->tag()), kProviderRevision};
        request.chart = persistentId(face->tag());
        request.faceClass = {2, face->tag()};
        request.desiredNormalSign = 1;
        request.expectedGeometryStamp =
            {kGeometryProviderId, kProviderRevision};
        request.expectedChartCertificateStamp =
            {request.chart + 1U, kProviderRevision};
        request.expectedCandidateProviderStamp = candidateStampFor(mode);
        request.outer.role = SurfaceChainRole::OuterBoundary;

        FaceContext context;
        context.face = face;
        context.spec = &spec;
        context.cache=cache;
        context.facePersistentId=request.faceStamp.persistentId;
        face->paramBounds(context.lo, context.hi);
        context.metricScale = face->metricScale();
        if (!finite(context.lo) || !finite(context.hi) ||
            !finite(context.metricScale) ||
            !(context.hi.u > context.lo.u) ||
            !(context.hi.v > context.lo.v) ||
            !(context.metricScale.u > 0.0) ||
            !(context.metricScale.v > 0.0)) {
            built.reason = ProductionSurfaceReason::InvalidBoundaryContract;
            return built;
        }
        const Vec2 center{(context.lo.u + context.hi.u) * 0.5,
                          (context.lo.v + context.hi.v) * 0.5};
        const double du = context.hi.u - context.lo.u;
        const double dv = context.hi.v - context.lo.v;
        if (!finite(center) || !detail::positiveFinite(du) ||
            !detail::positiveFinite(dv) ||
            !detail::positiveFinite(du * context.metricScale.u) ||
            !detail::positiveFinite(dv * context.metricScale.v))
            detail::failSizing(spec);
        const Vec3 centerPoint = face->pointAt(center);
        context.referenceSize=evaluatedFaceSize(
            context,centerPoint,face->curvatureAt(center));
        if (!computeChartSign(*face, context.lo, context.hi,
                              context.chartPositiveToNormalSign)) {
            built.reason = ProductionSurfaceReason::InvalidBoundaryContract;
            return built;
        }

        const double physicalScale = std::max(
            {1.0, std::abs(centerPoint.x), std::abs(centerPoint.y),
             std::abs(centerPoint.z), context.referenceSize});
        const double roundoff = 128.0 *
            std::numeric_limits<double>::epsilon() * physicalScale;
        const double roundoffArea = roundoff * roundoff;
        const double sizeArea =
            context.referenceSize * context.referenceSize * 1.0e-12;
        const double mappingDeviation = context.referenceSize * 0.25;
        if (!detail::positiveFinite(physicalScale) ||
            !detail::positiveFinite(roundoff) ||
            !detail::positiveFinite(roundoffArea) ||
            !detail::positiveFinite(sizeArea) ||
            !detail::positiveFinite(mappingDeviation) ||
            !detail::positiveFinite(
                context.referenceSize / context.metricScale.u) ||
            !detail::positiveFinite(
                context.referenceSize / context.metricScale.v))
            detail::failSizing(spec);
        request.physicalPolicy.minPhysicalDoubleArea =
            std::max(roundoffArea, sizeArea);
        request.physicalPolicy.maxVertexMappingError = roundoff;
        request.physicalPolicy.maxArtifactPositionError = roundoff;
        request.physicalPolicy.maxMappingDeviation =
            std::max(roundoff, mappingDeviation);
        request.physicalPolicy.minNormalAlignment = 0.8;

        std::map<Index, FaceLocalVertexId> localId;
        const auto loops = face->boundaryLoops();
        for (std::size_t loopIndex = 0; loopIndex < loops.size(); ++loopIndex) {
        const auto& edgeTags = loops[loopIndex];
        auto& targetLoop = loopIndex == 0 ? request.outer : request.holes.emplace_back();
        targetLoop.role = loopIndex == 0 ? SurfaceChainRole::OuterBoundary
                                        : SurfaceChainRole::HoleBoundary;
        std::vector<bool> reversed(edgeTags.size(), false);
        auto orient = [&](bool firstReversed, std::vector<bool>& result) {
            result[0] = firstReversed;
            int start = firstReversed
                            ? edgeData.at(edgeTags[0]).edge->endVertexTag()
                            : edgeData.at(edgeTags[0]).edge->startVertexTag();
            int end = firstReversed
                          ? edgeData.at(edgeTags[0]).edge->startVertexTag()
                          : edgeData.at(edgeTags[0]).edge->endVertexTag();
            for (std::size_t i = 1; i < edgeTags.size(); ++i) {
                const auto* edge = edgeData.at(edgeTags[i]).edge;
                if (edge->startVertexTag() == end) {
                    result[i] = false;
                    end = edge->endVertexTag();
                } else if (edge->endVertexTag() == end) {
                    result[i] = true;
                    end = edge->startVertexTag();
                } else {
                    return false;
                }
            }
            return end == start;
        };
        if (!orient(false, reversed) && !orient(true, reversed)) {
            built.reason = ProductionSurfaceReason::InvalidBoundaryContract;
            return built;
        }

        for (std::size_t ei = 0; ei < edgeTags.size(); ++ei) {
            const int edgeTag = edgeTags[ei];
            const auto& data = edgeData.at(edgeTag);
            const auto& natural = data.naturalNodes;
            SurfaceChain chain;
            chain.chainPersistentId = persistentId(edgeTag);
            chain.role = targetLoop.role;
            chain.classification = {1, edgeTag};
            chain.closed = false;
            chain.owningCurve = ArtifactStamp{persistentId(edgeTag),
                                              kProviderRevision};
            const std::size_t count = natural.size();
            for (std::size_t oi = 0; oi < count; ++oi) {
                const std::size_t ni = reversed[ei] ? count - 1 - oi : oi;
                const Index global = natural[ni];
                auto [idIt, inserted] = localId.emplace(
                    global, static_cast<FaceLocalVertexId>(global) + 1U);
                const FaceLocalVertexId id = idIt->second;
                if (inserted) {
                    const Vec2 uv = face->paramAt(destination.node(global));
                    if (!finite(uv)) {
                        built.reason = ProductionSurfaceReason::InvalidBoundaryContract;
                        return built;
                    }
                    request.vertices.push_back(
                        {id, {uv, request.chart}, destination.nodeClass(global)});
                    context.declaredCoordinates.emplace(id,
                                                        destination.node(global));
                }
                SurfaceNodeIdentity nodeIdentity;
                const auto refined=refinedCurveIdentity.find(global);
                if(refined!=refinedCurveIdentity.end()) {
                    nodeIdentity=refined->second;
                    if(nodeIdentity.parentArtifact!=*chain.owningCurve) {
                        built.reason=ProductionSurfaceReason::
                            AmbiguousBoundaryDiscretization;
                        return built;
                    }
                }
                const SharedCurveNodeKey key{
                    persistentId(edgeTag),kProviderRevision,
                    refined==refinedCurveIdentity.end()?ni:0,nodeIdentity};
                if(refined!=refinedCurveIdentity.end()) {
                    const auto [identityIt,identityInserted]=
                        context.declaredPersistentIdentities.emplace(
                            id,nodeIdentity);
                    if(!identityInserted&&identityIt->second!=nodeIdentity) {
                        built.reason=ProductionSurfaceReason::
                            AmbiguousBoundaryDiscretization;
                        return built;
                    }
                }
                chain.uses.push_back({id, key});
                const SharedCurveNodeBinding binding{
                    key, global, destination.nodeClass(global)};
                const auto [bindingIt, bindingInserted] =
                    bindings.emplace(key, binding);
                if (!bindingInserted &&
                    (bindingIt->second.existingGlobalIndex != global ||
                     !sameClass(bindingIt->second.classification,
                                binding.classification))) {
                    built.reason =
                        ProductionSurfaceReason::AmbiguousBoundaryDiscretization;
                    return built;
                }
            }
            for (std::size_t oi = 0; oi + 1 < count; ++oi) {
                const std::size_t naturalOrdinal = reversed[ei]
                    ? count - 2 - oi
                    : oi;
                chain.segmentEntityIds.push_back(
                    segmentId(edgeTag, naturalOrdinal));
            }
            targetLoop.pieces.push_back(std::move(chain));
        }
        }
        const auto addFacePoints=[&](
            std::span<const t2c_detail::PrivateFacePoint> points) {
            for (const auto& point : points) {
                if (point.facePersistentId != request.faceStamp.persistentId)
                    continue;
                if (!finite(point.parameter) || !finite(point.coordinate)) {
                    built.reason =
                        ProductionSurfaceReason::InvalidBoundaryContract;
                    return false;
                }
                context.privateInteriorPoints.push_back(point);
                request.refinedNodes.push_back(
                    {point.identity,{point.parameter,request.chart},
                     request.faceClass});
            }
            return true;
        };
        if(!addFacePoints(facePoints)||!addFacePoints(proposedFacePoints))
            return built;
        std::sort(request.vertices.begin(), request.vertices.end(),
                  [](const auto& a, const auto& b) {
                      return a.faceLocalVertexId < b.faceLocalVertexId;
                  });
        built.contexts.emplace(request.chart, std::move(context));
        built.request.faces.push_back(std::move(request));
    }
    for (const auto& [key, binding] : bindings) {
        (void)key;
        built.request.sharedBindings.push_back(binding);
    }
    return built;
}

ProductionSurfaceStatus statusFor(const SurfaceSetOutcome& surface) noexcept {
    switch (surface.setStatus) {
        case SurfaceSetStatus::Published:
            return ProductionSurfaceStatus::Published;
        case SurfaceSetStatus::Rejected:
            return ProductionSurfaceStatus::Rejected;
        case SurfaceSetStatus::Refused:
            return ProductionSurfaceStatus::Unsupported;
        case SurfaceSetStatus::Failed:
            return ProductionSurfaceStatus::Failed;
    }
    return ProductionSurfaceStatus::Failed;
}

ProductionSurfacePhase phaseFor(const SurfaceFailure& failure) noexcept {
    switch (failure.allocation) {
        case SurfaceAllocationPhase::GeometryProvider:
            return ProductionSurfacePhase::BoundaryMapping;
        case SurfaceAllocationPhase::Candidate:
            return ProductionSurfacePhase::ScalableConstruction;
        case SurfaceAllocationPhase::Recovery:
            return ProductionSurfacePhase::CandidateBoundaryRecovery;
        case SurfaceAllocationPhase::DomainExtraction:
        case SurfaceAllocationPhase::Audit:
        case SurfaceAllocationPhase::SetMapping:
            return ProductionSurfacePhase::SurfaceBuild;
        case SurfaceAllocationPhase::PublicationPreparation:
            return ProductionSurfacePhase::Publication;
        case SurfaceAllocationPhase::None:
            return ProductionSurfacePhase::SurfaceBuild;
    }
    return ProductionSurfacePhase::SurfaceBuild;
}

}  // namespace

ProductionSurfaceOutcome::ProductionSurfaceOutcome() = default;
ProductionSurfaceOutcome::~ProductionSurfaceOutcome() = default;
ProductionSurfaceOutcome::ProductionSurfaceOutcome(
    ProductionSurfaceOutcome&&) noexcept = default;
ProductionSurfaceOutcome& ProductionSurfaceOutcome::operator=(
    ProductionSurfaceOutcome&&) noexcept = default;

detail::ProductionSurfaceSetStage detail::buildProductionSurfaceSetForV3(
    const geo::Model& model, const sizing::MeshSpec& spec,
    ArtifactStamp setStamp, mesh::Mesh& destination,
    std::span<const t2c_detail::PrivateCurveNodeIdentity> curveNodeIdentities,
    std::span<const t2c_detail::PrivateFacePoint> facePoints,
    t2c_detail::SurfaceFaceCacheSession* cacheSession,
    std::span<const t2c_detail::DyadicEdgeSplit> proposedEdgeSplits,
    std::span<const t2c_detail::PrivateFacePoint> proposedFacePoints,
    SurfaceCandidateConstructionMode constructionMode) noexcept {
    ProductionSurfaceSetStage out;
    const auto finalizeCacheException=[&](
            t2c_detail::SurfaceFaceCacheReason reason) noexcept {
        if(!cacheSession)return;
        auto* channel=cacheSession->diagnosticChannel();
        if(!channel)return;
        // The v4 coordinator preallocates both spans to the canonical face
        // count before entering this route.  Whole-cache failures that occur
        // while forming authorities/dirty scratch therefore remain fully
        // itemized without allocating in the catch path.  Never overwrite a
        // more precise diagnostic produced by a deeper cache stage.
        const auto existing=cacheSession->diagnostic();
        if(existing.status==t2c_detail::SurfaceFaceCacheStatus::None)
            cacheSession->setDiagnostic(
                t2c_detail::SurfaceFaceCacheStatus::Failed,reason,
                t2c_detail::SurfaceFaceCachePhase::Preparation,true,0,0,0);
        const auto faceCount=channel->unattemptedStorage.size();
        const bool complete=channel->processedStorage.size()>=faceCount;
        channel->processedUsed=0;
        channel->unattemptedUsed=0;
        SurfaceFailure failure;
        failure.reason=reason==t2c_detail::SurfaceFaceCacheReason::
                AllocationFailed
            ? SurfaceReason::AllocationFailed
            : SurfaceReason::InternalError;
        failure.allocation=reason==t2c_detail::SurfaceFaceCacheReason::
                AllocationFailed
            ? SurfaceAllocationPhase::SetMapping
            : SurfaceAllocationPhase::None;
        if(complete)for(std::size_t i=0;i<faceCount;++i)
            channel->processedStorage[channel->processedUsed++]={
                channel->unattemptedStorage[i],
                SurfaceCacheFaceDisposition::Stopped,
                SurfaceFaceStatus::Failed,failure};
        channel->reportComplete=complete;
        channel->hasFailedFace=false;
        channel->failedFacePersistentId=0;
        const auto exact=cacheSession->diagnostic();
        cacheSession->setDiagnostic(exact.status,exact.reason,exact.phase,
            complete,exact.failedFace,channel->processedUsed,0);
    };
    try {
        if(cacheSession)cacheSession->beginAttempt();
        const sizing::MeshSpec checkedSpec = makeCheckedSpec(spec);
        auto built = buildRequest(model,checkedSpec,destination,
                                  curveNodeIdentities,facePoints,cacheSession,
                                  proposedFacePoints,constructionMode);
        built.request.setStamp = setStamp;
        if (built.reason != ProductionSurfaceReason::None) {
            out.reason = built.reason;
            out.phase = ProductionSurfacePhase::Preflight;
            const bool unsupported =
                built.reason == ProductionSurfaceReason::EmptyFaceSet;
            out.surface.setStatus = unsupported ? SurfaceSetStatus::Refused
                                                : SurfaceSetStatus::Rejected;
            out.surface.setFailure.reason = SurfaceReason::InvalidLoopTopology;
            auto faces = model.faces();
            std::sort(faces.begin(),faces.end(),[](const auto* a,const auto* b) {
                const int at = a ? a->tag() : -1;
                const int bt = b ? b->tag() : -1;
                return at < bt;
            });
            out.surface.processed.reserve(faces.size());
            for (const auto* face : faces) {
                SurfaceFaceRow row;
                row.facePersistentId = face && face->tag() >= 0
                    ? persistentId(face->tag()) : 0;
                row.status = unsupported ? SurfaceFaceStatus::Refused
                                         : SurfaceFaceStatus::Rejected;
                row.failure = out.surface.setFailure;
                out.surface.processed.push_back(std::move(row));
            }
            out.surface.reportComplete = true;
            return out;
        }
        KernelGeometryProvider geometry(built.contexts);
        ScalableCandidateProvider candidates(
            built.contexts,geometry,constructionMode);
        if(cacheSession) {
            std::map<SharedCurveNodeKey,SharedCurveNodeBinding> ledger;
#if defined(PAACS_TOPO_TESTING)
            {
                const t2c_detail::SurfaceCacheAllocationZoneScope allocationZone(
                    t2c_detail::SurfaceCacheAllocationZone::AuthorityKey);
#endif
                for(const auto& binding:built.request.sharedBindings)
                    ledger.emplace(binding.key,binding);
#if defined(PAACS_TOPO_TESTING)
            }
#endif
            std::map<FacePersistentId,std::vector<std::uint8_t>> authorities;
            std::vector<FacePersistentId> canonicalFaces;
            std::vector<FacePersistentId> expectedDirtyFaces;
            std::uint64_t dependencyVisits=0;
            {
#if defined(PAACS_TOPO_TESTING)
            const t2c_detail::SurfaceCacheAllocationZoneScope dirtyAllocationZone(
                t2c_detail::SurfaceCacheAllocationZone::DirtySet);
#endif
            canonicalFaces.reserve(built.request.faces.size());
            expectedDirtyFaces.reserve(built.request.faces.size());
            const bool initiallyEmpty=cacheSession->empty();
            for(const auto& face:built.request.faces) {
                const auto context=built.contexts.find(face.chart);
                if(context==built.contexts.end())
                    throw std::logic_error("missing cache face context");
                const auto certificate=geometry.chartCertificate(face.chart);
#if defined(PAACS_TOPO_TESTING)
                {
                    const t2c_detail::SurfaceCacheAllocationZoneScope
                        authorityAllocationZone(
                            t2c_detail::SurfaceCacheAllocationZone::AuthorityKey);
#endif
                    authorities.emplace(face.faceStamp.persistentId,
                        t2c_detail::encodeSurfaceFaceCacheAuthority(
                            setStamp,face,certificate,ledger,destination,
                            context->second.lo,context->second.hi,
                            context->second.metricScale));
#if defined(PAACS_TOPO_TESTING)
                }
#endif
                canonicalFaces.push_back(face.faceStamp.persistentId);
                bool dirty=initiallyEmpty;
                if(dependencyVisits==
                   (std::numeric_limits<std::uint64_t>::max)())
                    throw std::length_error("surface dependency work overflow");
                ++dependencyVisits; // canonical face visit
                for(const auto& point:proposedFacePoints) {
                    if(dependencyVisits==
                       (std::numeric_limits<std::uint64_t>::max)())
                        throw std::length_error(
                            "surface dependency work overflow");
                    ++dependencyVisits;
                    if(point.facePersistentId==face.faceStamp.persistentId)
                        dirty=true;
                }
                const auto inspectChain=[&](const SurfaceChain& chain) {
                    if(dependencyVisits==
                       (std::numeric_limits<std::uint64_t>::max)())
                        throw std::length_error(
                            "surface dependency work overflow");
                    ++dependencyVisits;
                    if(!chain.owningCurve)return false;
                    bool affected=false;
                    for(const auto& split:proposedEdgeSplits) {
                        if(dependencyVisits==
                           (std::numeric_limits<std::uint64_t>::max)())
                            throw std::length_error(
                                "surface dependency work overflow");
                        ++dependencyVisits;
                        affected=affected||split.identity.parentArtifact==
                            *chain.owningCurve;
                    }
                    return affected;
                };
                {
                    for(const auto& chain:face.outer.pieces)
                        dirty=inspectChain(chain)||dirty;
                    for(const auto& hole:face.holes)
                        for(const auto& chain:hole.pieces)
                            dirty=inspectChain(chain)||dirty;
                    for(const auto& chain:face.embedded)
                        dirty=inspectChain(chain)||dirty;
                }
                if(dirty)expectedDirtyFaces.push_back(
                    face.faceStamp.persistentId);
            }
            std::sort(canonicalFaces.begin(),canonicalFaces.end());
            std::sort(expectedDirtyFaces.begin(),expectedDirtyFaces.end());
            if(!cacheSession->configureExpectedDirty(
                    canonicalFaces,expectedDirtyFaces,dependencyVisits)) {
                // The cached adapter transports the exact cache diagnostic
                // through its fixed report envelope.  It is still invoked so
                // every face receives one processed/unattempted outcome.
            }
            }
            const CheckedSizingCacheEvaluator evaluator(checkedSpec);
            out.surface=t2c_detail::buildAndPublishSurfaceSetCached(
                built.request,destination,geometry,candidates,{},authorities,
                *cacheSession,evaluator,candidates);
            const auto rows=cacheSession->attemptMetrics();
            out.metrics.assign(rows.begin(),rows.end());
        } else {
            out.surface = buildAndPublishSurfaceSet(
                built.request,destination,geometry,candidates);
            out.metrics = candidates.takeRows();
        }
        if (out.surface.setStatus != SurfaceSetStatus::Published &&
            candidates.recoveryFailure()) {
            // The candidate provider throws to stop construction. Retain the
            // structured recovery receipt alongside that provider failure.
            const auto copyRecovery = [&](SurfaceFailure& failure) noexcept {
                failure.bndReject = candidates.recoveryFailure()->reject;
                failure.bndRefuse = candidates.recoveryFailure()->refuse;
                failure.bndFail = candidates.recoveryFailure()->fail;
            };
            copyRecovery(out.surface.setFailure);
            for (auto& row : out.surface.processed)
                if (row.status != SurfaceFaceStatus::Ready)
                    copyRecovery(row.failure);
        }
        out.phase = phaseFor(out.surface.setFailure);
        if (sizingEvaluationFailed(checkedSpec)) {
            out.reason = ProductionSurfaceReason::SizingEvaluationFailed;
        } else if (out.surface.setStatus == SurfaceSetStatus::Published) {
            out.reason = ProductionSurfaceReason::None;
            out.phase = ProductionSurfacePhase::Publication;
        } else if (out.surface.setFailure.reason == SurfaceReason::AllocationFailed) {
            out.reason = ProductionSurfaceReason::AllocationFailed;
        } else if (out.surface.setFailure.reason ==
                   SurfaceReason::CandidateProviderError) {
            out.reason = ProductionSurfaceReason::CandidateConstructionFailure;
        } else if (out.surface.setFailure.reason ==
                       SurfaceReason::GeometryProviderError ||
                   out.surface.setFailure.reason ==
                       SurfaceReason::ProviderDomainRejected) {
            out.reason = ProductionSurfaceReason::GeometryProviderFailure;
        } else {
            out.reason = ProductionSurfaceReason::SurfaceAdapterStopped;
        }
        return out;
    } catch (const SizingEvaluationError&) {
        finalizeCacheException(
            t2c_detail::SurfaceFaceCacheReason::AuthorityChanged);
        out.reason = ProductionSurfaceReason::SizingEvaluationFailed;
        out.phase = ProductionSurfacePhase::Spacing;
        out.surface.setStatus = SurfaceSetStatus::Failed;
        out.surface.setFailure.reason = SurfaceReason::InternalError;
    } catch (const std::bad_alloc&) {
        finalizeCacheException(
            t2c_detail::SurfaceFaceCacheReason::AllocationFailed);
        out.reason = ProductionSurfaceReason::AllocationFailed;
        out.phase = ProductionSurfacePhase::Preflight;
        out.surface.setStatus = SurfaceSetStatus::Failed;
        out.surface.setFailure.reason = SurfaceReason::AllocationFailed;
    } catch (const std::exception&) {
        finalizeCacheException(t2c_detail::SurfaceFaceCacheReason::Invalid);
        out.reason = ProductionSurfaceReason::GeometryProviderFailure;
        out.phase = ProductionSurfacePhase::BoundaryMapping;
        out.surface.setStatus = SurfaceSetStatus::Failed;
        out.surface.setFailure.reason = SurfaceReason::GeometryProviderError;
    } catch (...) {
        finalizeCacheException(t2c_detail::SurfaceFaceCacheReason::Invalid);
        out.reason = ProductionSurfaceReason::GeometryProviderFailure;
        out.phase = ProductionSurfacePhase::BoundaryMapping;
        out.surface.setStatus = SurfaceSetStatus::Failed;
        out.surface.setFailure.reason = SurfaceReason::GeometryProviderError;
    }
    return out;
}

ProductionSurfaceOutcome meshFacesConstrained(
    const geo::Model& model, const sizing::MeshSpec& spec,
    mesh::Mesh& destination) noexcept {
    return meshFacesConstrained(model, spec, destination, {});
}

ProductionSurfaceOutcome meshFacesConstrained(
    const geo::Model& model, const sizing::MeshSpec& spec,
    mesh::Mesh& destination, const ProductionSurfaceOptions& options) noexcept {
    ProductionSurfaceOutcome out;
    const auto before = shapeOf(destination);
    try {
        const sizing::MeshSpec checkedSpec = detail::makeCheckedSpec(spec);
        const auto mode = options.recoverPlanarConstraints
            ? detail::SurfaceCandidateConstructionMode::ScalablePlanarRecovery
            : detail::SurfaceCandidateConstructionMode::ScalableDelaunay;
        auto built = buildRequest(model,checkedSpec,destination,{},{},nullptr,{},mode);
        if (built.reason != ProductionSurfaceReason::None) {
            out.status = built.reason == ProductionSurfaceReason::EmptyFaceSet
                             ? ProductionSurfaceStatus::Unsupported
                             : ProductionSurfaceStatus::Rejected;
            out.reason = built.reason;
            out.phase = ProductionSurfacePhase::Preflight;
            out.destinationUnchanged = shapeOf(destination).nodes == before.nodes &&
                                       shapeOf(destination).cells == before.cells;
            return out;
        }
        KernelGeometryProvider geometry(built.contexts);
        ScalableCandidateProvider candidates(built.contexts, geometry, mode,
                                              options.maximumConstraintWork);
        auto surface = buildAndPublishSurfaceSet(
            built.request, destination, geometry, candidates);
        // Candidate rows were fully allocated before publication. Move their
        // vector after the accepted append tail; copying here would re-open a
        // recoverable allocation point after global mutation.
        out.faceMetrics = candidates.takeRows();
        out.planarRecovery = candidates.planarRecovery();
        out.reportComplete = surface.reportComplete;
        out.processedFaces = surface.processed.size();
        out.unattemptedFaces = surface.unattempted.size();
        out.status = statusFor(surface);
        out.phase = phaseFor(surface.setFailure);
        if (detail::sizingEvaluationFailed(checkedSpec)) {
            out.status = ProductionSurfaceStatus::Failed;
            out.reason = ProductionSurfaceReason::SizingEvaluationFailed;
            out.phase = ProductionSurfacePhase::Spacing;
        } else if (surface.setStatus == SurfaceSetStatus::Published) {
            out.reason = ProductionSurfaceReason::None;
            out.phase = ProductionSurfacePhase::Publication;
            out.destinationUnchanged = false;
        } else if (candidates.recoveryFailure()) {
            out.reason =
                ProductionSurfaceReason::CandidateBoundaryRecoveryStopped;
            out.phase = ProductionSurfacePhase::CandidateBoundaryRecovery;
            if (candidates.recoveryFailure()->reject !=
                bnd::SteinerReject::None) {
                out.recoveryKind = ProductionSurfaceRecoveryKind::Rejected;
                out.recoveryDetail = static_cast<std::uint32_t>(
                    candidates.recoveryFailure()->reject);
            } else if (candidates.recoveryFailure()->refuse !=
                       bnd::SteinerRefuse::None) {
                out.recoveryKind = ProductionSurfaceRecoveryKind::Refused;
                out.recoveryDetail = static_cast<std::uint32_t>(
                    candidates.recoveryFailure()->refuse);
            } else if (candidates.recoveryFailure()->fail !=
                       bnd::SteinerFail::None) {
                out.recoveryKind = ProductionSurfaceRecoveryKind::Failed;
                out.recoveryDetail = static_cast<std::uint32_t>(
                    candidates.recoveryFailure()->fail);
            }
        } else if (surface.setFailure.reason ==
                   SurfaceReason::CandidateProviderError) {
            out.reason = ProductionSurfaceReason::CandidateConstructionFailure;
            out.phase = ProductionSurfacePhase::ScalableConstruction;
        } else if (surface.setFailure.reason ==
                   SurfaceReason::GeometryProviderError ||
                   surface.setFailure.reason ==
                       SurfaceReason::ProviderDomainRejected) {
            out.reason = ProductionSurfaceReason::GeometryProviderFailure;
            out.phase = ProductionSurfacePhase::BoundaryMapping;
        } else if (surface.setFailure.reason ==
                   SurfaceReason::AllocationFailed) {
            out.reason = ProductionSurfaceReason::AllocationFailed;
        } else {
            out.reason = ProductionSurfaceReason::SurfaceAdapterStopped;
        }
        out.retainedArtifact = std::move(surface.retainedArtifact);
        if (!out.published())
            out.destinationUnchanged =
                shapeOf(destination).nodes == before.nodes &&
                shapeOf(destination).cells == before.cells;
        return out;
    } catch (const detail::SizingEvaluationError&) {
        out.status = ProductionSurfaceStatus::Failed;
        out.reason = ProductionSurfaceReason::SizingEvaluationFailed;
        out.phase = ProductionSurfacePhase::Spacing;
    } catch (const std::bad_alloc&) {
        out.status = ProductionSurfaceStatus::Failed;
        out.reason = ProductionSurfaceReason::AllocationFailed;
        out.phase = ProductionSurfacePhase::Preflight;
    } catch (...) {
        out.status = ProductionSurfaceStatus::Failed;
        out.reason = ProductionSurfaceReason::GeometryProviderFailure;
        out.phase = ProductionSurfacePhase::BoundaryMapping;
    }
    out.destinationUnchanged = shapeOf(destination).nodes == before.nodes &&
                               shapeOf(destination).cells == before.cells;
    return out;
}

}  // namespace paacs::algo

#if defined(_MSC_VER)
#pragma warning(pop)
#endif
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
