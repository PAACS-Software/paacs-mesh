// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/Meshers.hpp"
#include "paacs/algo/PreOcctMeshingCore.hpp"
#include "paacs/algo/PreparedSweepMesher.hpp"
#include "paacs/algo/ProductionSurfaceMesher.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <set>

namespace paacs::algo {
namespace {
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point a) {
    return std::chrono::duration<double>(Clock::now() - a).count();
}
double length(Vec3 p) {
    return std::hypot(p.x, p.y, p.z);
}
bool finite(Vec3 p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
struct Stop {
    GeometryReason reason;
};
std::uint64_t add(std::uint64_t a, std::uint64_t b) {
    if (b > UINT64_MAX - a)
        throw Stop{GeometryReason::OutputBudget};
    return a + b;
}
std::uint64_t mul(std::uint64_t a, std::uint64_t b) {
    if (a && b > UINT64_MAX / a)
        throw Stop{GeometryReason::OutputBudget};
    return a * b;
}
std::uint64_t bytes(std::uint64_t n, std::uint64_t c, std::uint64_t conn) {
    return add(add(mul(n, sizeof(Vec3) + sizeof(mesh::Classification)), mul(c, sizeof(mesh::Cell))),
               mul(conn, sizeof(Index)));
}
std::uint64_t bytes(const mesh::Mesh &m) {
    return bytes(static_cast<std::uint64_t>(m.nodeCount()),
                 static_cast<std::uint64_t>(m.cellCount()),
                 static_cast<std::uint64_t>(m.connectivityCount()));
}
void limits(std::uint64_t n, std::uint64_t c, std::uint64_t b, const GeometryMeshingOptions &o) {
    if (n > o.maximumNodes || c > o.maximumCells || b > o.maximumOutputBytes)
        throw Stop{GeometryReason::OutputBudget};
}
void constructionWorkFloor(std::uint64_t nodes, std::uint64_t cells,
                           const GeometryMeshingOptions &o) {
    // At least one construction operation is necessary per output record.
    // This is a per-body floor, not a count of every arithmetic operation.
    if (nodes > o.maximumMeshingWork || cells > o.maximumMeshingWork - nodes)
        throw Stop{GeometryReason::MeshingBudget};
}
void cancelled(const GeometryMeshingOptions &o) {
    if (o.cancellation && o.cancellation->load(std::memory_order_relaxed))
        throw Stop{GeometryReason::Cancelled};
}
void validate(const GeometryMeshingOptions &o) {
    const auto ratio = [](double x) { return std::isfinite(x) && (x == 0 || x >= 1); };
    if (!(o.targetSize > 0) || !std::isfinite(o.targetSize) ||
        !(o.relativeRecognitionTolerance > 0 && o.relativeRecognitionTolerance < 1e-3) ||
        !(o.minimumElementQuality >= 0 && o.minimumElementQuality <= 1) || o.threads < 0 ||
        !(o.warningElementQuality >= 0 && o.warningElementQuality <= 1) ||
        o.strategy == GeometryStrategy::Unsupported ||
        static_cast<unsigned>(o.strategy) > static_cast<unsigned>(GeometryStrategy::Unsupported) ||
        (o.spatialField && !o.spatialField->valid()) || !ratio(o.warningElementSizeRatio) ||
        !ratio(o.maximumElementSizeRatio) || !o.edgeGeneralRecovery.valid())
        throw Stop{GeometryReason::InvalidRequest};
    cancelled(o);
}
struct Work {
    std::uint64_t count = 0;
    const GeometryMeshingOptions &o;
    void tick(std::uint64_t n = 1) {
        if (n > o.maximumAnalysisWork - count)
            throw Stop{GeometryReason::AnalysisBudget};
        count += n;
        if ((count & 1023U) == 0)
            cancelled(o);
    }
};
void hash(std::uint64_t &h, std::uint64_t v) {
    for (unsigned i = 0; i < 8; ++i) {
        h ^= (v >> (i * 8U)) & 255U;
        h *= 1099511628211ULL;
    }
}
using Pair = std::pair<int, int>;
Pair key(int a, int b) {
    return std::minmax(a, b);
}

class View final : public geo::Model {
  public:
    const geo::Model *original = nullptr;
    std::vector<const geo::GVertex *> vs;
    std::vector<const geo::GEdge *> es;
    std::vector<const geo::GFace *> fs;
    std::vector<const geo::GSolid *> ss;
    std::vector<const geo::GVertex *> vertices() const override { return vs; }
    std::vector<const geo::GEdge *> edges() const override { return es; }
    std::vector<const geo::GFace *> faces() const override { return fs; }
    std::vector<const geo::GSolid *> solids() const override { return ss; }
    const geo::GVertex *vertex(int t) const override { return original->vertex(t); }
    const geo::GEdge *edge(int t) const override { return original->edge(t); }
    const geo::GFace *face(int t) const override { return original->face(t); }
    const geo::GSolid *solid(int t) const override { return original->solid(t); }
};
struct FaceData {
    const geo::GFace *face = nullptr;
    std::vector<std::vector<int>> rings;
    std::set<int> vertices, edges;
};
struct Receipt {
    View view;
    GeometryBodyReport report;
    std::vector<FaceData> faces;
    std::map<Pair, int> edgeTags;
    std::map<int, int> topVertex, topEdge, sideFace, verticalEdge;
    std::size_t source = 0, target = 0;
    MappedHexBlockRequest mapped;
    Vec3 lo, hi;
    double tolerance = 0;
    GeometryBodyControl control;
    GeometryStrategy requestedStrategy = GeometryStrategy::Auto;
    std::vector<double> profileFractions;
};
std::vector<int> ringVertices(const geo::Model &model, const std::vector<int> &edges, Work &work) {
    if (edges.size() < 2)
        throw Stop{GeometryReason::UnsupportedChart};
    const auto *first = model.edge(edges.front());
    if (!first)
        throw Stop{GeometryReason::InvalidTopology};
    for (bool reverse : {false, true}) {
        std::vector<int> result;
        result.reserve(edges.size());
        const int start = reverse ? first->endVertexTag() : first->startVertexTag();
        int current = start;
        bool valid = true;
        for (int tag : edges) {
            work.tick();
            const auto *e = model.edge(tag);
            if (!e || e->tag() != tag) {
                valid = false;
                break;
            }
            result.push_back(current);
            if (e->startVertexTag() == current)
                current = e->endVertexTag();
            else if (e->endVertexTag() == current)
                current = e->startVertexTag();
            else {
                valid = false;
                break;
            }
        }
        if (valid && current == start) {
            std::set<int> unique(result.begin(), result.end());
            if (unique.size() != result.size())
                throw Stop{GeometryReason::UnsupportedChart};
            return result;
        }
    }
    throw Stop{GeometryReason::InvalidTopology};
}
Receipt inspect(const geo::Model &model, const geo::GSolid &solid, Work &work) {
    Receipt r;
    r.view.original = &model;
    r.view.ss = {&solid};
    r.report.solidTag = solid.tag();
    solid.bbox(r.lo, r.hi);
    if (!finite(r.lo) || !finite(r.hi) || !(r.hi.x > r.lo.x) || !(r.hi.y > r.lo.y) ||
        !(r.hi.z > r.lo.z))
        throw Stop{GeometryReason::NonfiniteGeometry};
    const double scale = length(r.hi - r.lo),
                 coordinate = std::max({std::abs(r.lo.x), std::abs(r.lo.y), std::abs(r.lo.z),
                                        std::abs(r.hi.x), std::abs(r.hi.y), std::abs(r.hi.z)});
    if (!std::isfinite(scale) || scale <= 0)
        throw Stop{GeometryReason::NonfiniteGeometry};
    r.tolerance = scale * work.o.relativeRecognitionTolerance +
                  coordinate * 32 * std::numeric_limits<double>::epsilon();
    auto tags = solid.boundaryFaceTags();
    std::sort(tags.begin(), tags.end());
    if (tags.empty() || std::adjacent_find(tags.begin(), tags.end()) != tags.end())
        throw Stop{GeometryReason::InvalidTopology};
    std::set<int> edges, vertices;
    std::map<int, unsigned> incidence;
    r.report.geometryHash = 1469598103934665603ULL;
    hash(r.report.geometryHash, static_cast<std::uint64_t>(solid.tag()));
    for (int tag : tags) {
        work.tick();
        const auto *f = model.face(tag);
        if (!f || tag < 0 || f->tag() != tag)
            throw Stop{GeometryReason::InvalidTopology};
        if (f->periodicU() || f->periodicV())
            throw Stop{GeometryReason::UnsupportedChart};
        auto loops = f->boundaryLoops();
        if (loops.empty())
            throw Stop{GeometryReason::InvalidTopology};
        FaceData fd;
        fd.face = f;
        for (const auto &loop : loops) {
            auto ring = ringVertices(model, loop, work);
            fd.vertices.insert(ring.begin(), ring.end());
            fd.rings.push_back(std::move(ring));
            for (int edge : loop) {
                if (!fd.edges.insert(edge).second)
                    throw Stop{GeometryReason::UnsupportedChart};
                edges.insert(edge);
                ++incidence[edge];
            }
        }
        vertices.insert(fd.vertices.begin(), fd.vertices.end());
        r.report.holes += loops.size() - 1;
        if (f->kind() != geo::SurfaceKind::Plane)
            ++r.report.curvedFaces;
        hash(r.report.geometryHash, static_cast<std::uint64_t>(tag));
        hash(r.report.geometryHash, static_cast<unsigned>(f->kind()));
        for (int e : fd.edges)
            hash(r.report.geometryHash, static_cast<std::uint64_t>(e));
        r.faces.push_back(std::move(fd));
        r.view.fs.push_back(f);
    }
    for (int tag : edges) {
        work.tick();
        const auto *e = model.edge(tag);
        if (!e || incidence[tag] != 2)
            throw Stop{GeometryReason::InvalidTopology};
        if (!r.edgeTags.emplace(key(e->startVertexTag(), e->endVertexTag()), tag).second)
            throw Stop{GeometryReason::UnsupportedChart};
        r.view.es.push_back(e);
    }
    for (int tag : vertices) {
        work.tick();
        const auto *v = model.vertex(tag);
        if (!v || tag < 0 || v->tag() != tag || !finite(v->position()))
            throw Stop{GeometryReason::InvalidTopology};
        const auto p = v->position();
        hash(r.report.geometryHash, static_cast<std::uint64_t>(tag));
        hash(r.report.geometryHash, std::bit_cast<std::uint64_t>(p.x));
        hash(r.report.geometryHash, std::bit_cast<std::uint64_t>(p.y));
        hash(r.report.geometryHash, std::bit_cast<std::uint64_t>(p.z));
        r.view.vs.push_back(v);
    }
    r.report.vertices = vertices.size();
    r.report.edges = edges.size();
    r.report.faces = tags.size();
    return r;
}
bool recognizeSweep(Receipt &r, Work &work) {
    if (r.report.curvedFaces)
        return false;
    for (const auto *e : r.view.es)
        if (e->kind() != geo::CurveKind::Line)
            return false;
    for (std::size_t a = 0; a < r.faces.size(); ++a)
        for (std::size_t b = r.control.sourceFaceTag < 0 ? a + 1 : 0; b < r.faces.size(); ++b) {
            work.tick();
            if (a == b ||
                (r.control.sourceFaceTag >= 0 && r.faces[a].face->tag() != r.control.sourceFaceTag))
                continue;
            const auto &source = r.faces[a];
            const auto &target = r.faces[b];
            if (source.rings.size() != target.rings.size() ||
                source.vertices.size() != target.vertices.size() ||
                source.edges.size() != target.edges.size() ||
                source.vertices.size() * 2 != r.view.vs.size() ||
                source.edges.size() + 2 != r.faces.size())
                continue;
            bool disjoint = true;
            for (int id : source.vertices)
                if (target.vertices.count(id))
                    disjoint = false;
            if (!disjoint)
                continue;
            const auto p0 = r.view.vertex(*source.vertices.begin())->position();
            for (int candidate : target.vertices) {
                const auto shift = r.view.vertex(candidate)->position() - p0;
                Vec2 l, h;
                source.face->paramBounds(l, h);
                const auto normal = source.face->normalAt({(l.u + h.u) * .5, (l.v + h.v) * .5});
                if (!finite(normal) || length(normal) == 0 ||
                    std::abs(shift.dot(normal) / length(normal)) <= r.tolerance)
                    continue;
                std::map<int, int> matched;
                std::set<int> used;
                bool valid = true;
                for (int i : source.vertices) {
                    int found = -1;
                    for (int j : target.vertices) {
                        work.tick();
                        if (length(r.view.vertex(i)->position() + shift -
                                   r.view.vertex(j)->position()) <= r.tolerance) {
                            if (found >= 0) {
                                valid = false;
                                break;
                            }
                            found = j;
                        }
                    }
                    if (!valid || found < 0 || !used.insert(found).second) {
                        valid = false;
                        break;
                    }
                    matched.emplace(i, found);
                }
                if (!valid)
                    continue;
                std::map<int, int> topEdges, sides, vertical;
                for (const auto [i, j] : matched) {
                    const auto it = r.edgeTags.find(key(i, j));
                    if (it == r.edgeTags.end()) {
                        valid = false;
                        break;
                    }
                    vertical.emplace(i, it->second);
                }
                if (!valid)
                    continue;
                std::set<int> sideTags;
                for (int edge : source.edges) {
                    const auto *e = r.view.edge(edge);
                    const int i = e->startVertexTag(), j = e->endVertexTag();
                    const auto ti = r.edgeTags.find(key(matched.at(i), matched.at(j)));
                    if (ti == r.edgeTags.end() || !target.edges.count(ti->second)) {
                        valid = false;
                        break;
                    }
                    topEdges.emplace(edge, ti->second);
                    const std::set<int> expected{i, j, matched.at(i), matched.at(j)};
                    int side = -1;
                    for (std::size_t f = 0; f < r.faces.size(); ++f) {
                        work.tick();
                        if (f == a || f == b)
                            continue;
                        const auto &face = r.faces[f];
                        if (face.rings.size() == 1 && face.rings[0].size() == 4 &&
                            face.vertices == expected) {
                            if (side >= 0) {
                                valid = false;
                                break;
                            }
                            side = face.face->tag();
                        }
                    }
                    if (!valid || side < 0 || !sideTags.insert(side).second) {
                        valid = false;
                        break;
                    }
                    sides.emplace(edge, side);
                }
                if (valid && sideTags.size() + 2 == r.faces.size()) {
                    r.source = a;
                    r.target = b;
                    r.topVertex = std::move(matched);
                    r.topEdge = std::move(topEdges);
                    r.sideFace = std::move(sides);
                    r.verticalEdge = std::move(vertical);
                    r.report.sourceFaceTag = source.face->tag();
                    r.report.targetFaceTag = target.face->tag();
                    r.report.sweep = shift;
                    return true;
                }
            }
        }
    return false;
}
std::uint32_t steps(double extent, double h) {
    const double n = std::ceil(extent / h);
    if (!std::isfinite(n) || n > UINT32_MAX)
        throw Stop{GeometryReason::OutputBudget};
    return static_cast<std::uint32_t>(std::max(1.0, n));
}
bool recognizeBlock(Receipt &r, const GeometryMeshingOptions &o) {
    if (r.faces.size() != 6 || r.view.vs.size() != 8 || r.view.es.size() != 12 ||
        r.faces[r.source].rings.size() != 1 || r.faces[r.source].rings[0].size() != 4)
        return false;
    auto ring = r.faces[r.source].rings[0];
    const auto point = [&](int id) { return r.view.vertex(id)->position(); };
    Vec3 origin = point(ring[0]), u = point(ring[1]) - origin, v = point(ring[3]) - origin,
         w = r.report.sweep;
    if (length(origin + u + v - point(ring[2])) > r.tolerance)
        return false;
    if (u.cross(v).dot(w) < 0) {
        std::swap(ring[1], ring[3]);
        std::swap(u, v);
    }
    if (!(u.cross(v).dot(w) > 0))
        return false;
    auto &g = r.mapped.geometry;
    g.origin = origin;
    g.axisU = u;
    g.axisV = v;
    g.axisW = w;
    g.cellsU = steps(length(u), r.report.targetSize);
    g.cellsV = steps(length(v), r.report.targetSize);
    g.cellsW = r.report.layers;
    r.mapped.growth.w = {StructuredAxisGrowthMode::EndToStartCellWidthRatio,
                         r.control.endToStartLayerWidthRatio, r.control.maximumAdjacentLayerRatio};
    r.mapped.growth.w.boundaryFractions = r.profileFractions;
    g.volumeClassification = {3, r.report.solidTag};
    g.maxNodes = o.maximumNodes;
    g.maxCells = o.maximumCells;
    std::array<int, 8> ids{};
    for (std::size_t i = 0; i < 4; ++i) {
        ids[i] = ring[i];
        ids[i + 4] = r.topVertex.at(ring[i]);
    }
    for (std::size_t i = 0; i < 8; ++i)
        g.vertexClassifications[i] = {0, ids[i]};
    const std::array<std::array<unsigned, 2>, 12> pairs{{{0, 1},
                                                         {1, 2},
                                                         {2, 3},
                                                         {3, 0},
                                                         {4, 5},
                                                         {5, 6},
                                                         {6, 7},
                                                         {7, 4},
                                                         {0, 4},
                                                         {1, 5},
                                                         {2, 6},
                                                         {3, 7}}};
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const auto p = pairs[i];
        const auto e = r.edgeTags.find(key(ids[p[0]], ids[p[1]]));
        if (e == r.edgeTags.end())
            return false;
        g.edgeClassifications[i] = {1, e->second};
    }
    const std::array<std::array<unsigned, 4>, 6> faces{
        {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}}};
    for (std::size_t i = 0; i < faces.size(); ++i) {
        std::set<int> wanted;
        for (auto k : faces[i])
            wanted.insert(ids[k]);
        int tag = -1;
        for (const auto &f : r.faces)
            if (f.rings.size() == 1 && f.vertices == wanted)
                tag = f.face->tag();
        if (tag < 0)
            return false;
        g.boundaryClassifications[i] = {2, tag};
    }
    r.mapped.executionThreads = o.threads;
    r.mapped.validation = MappedHexValidationMode::CertifiedAffine;
    return true;
}
void estimateBlock(Receipt &r, const GeometryMeshingOptions &o) {
    const auto &g = r.mapped.geometry;
    const std::uint64_t u = g.cellsU, v = g.cellsV, w = g.cellsW;
    const auto nodes = mul(mul(u + 1, v + 1), w + 1), volume = mul(mul(u, v), w),
               boundary = mul(2, add(add(mul(u, v), mul(v, w)), mul(u, w)));
    r.report.estimatedNodes = nodes;
    const bool tet = r.control.elements == GeometryElements::Tet4;
    r.report.estimatedCells = add(mul(volume, tet ? 6 : 1), mul(boundary, tet ? 2 : 1));
    r.report.estimatedOutputBytes = bytes(
        nodes, r.report.estimatedCells, add(mul(volume, tet ? 24 : 8), mul(boundary, tet ? 6 : 4)));
    r.report.exactCountEstimate = true;
    limits(nodes, r.report.estimatedCells, r.report.estimatedOutputBytes, o);
}
std::uint64_t sourceBound(const Receipt &r, const FaceData &f, const GeometryMeshingOptions &o) {
    Vec2 lo, hi;
    f.face->paramBounds(lo, hi);
    const auto metric = f.face->metricScale();
    const double target =
        r.report.usedSpatialField ? r.report.spatialMinimumBound : r.report.targetSize;
    const auto nx = steps((hi.u - lo.u) * metric.u, target),
               ny = steps((hi.v - lo.v) * metric.v, target * .8660254037844386);
    std::uint64_t count = add(mul(add(nx, 3), add(ny, 3)), 4);
    for (int e : f.edges)
        count = add(count, add(steps(r.view.edge(e)->length(), target), 1));
    if (count > o.maximumSourcePoints)
        throw Stop{GeometryReason::OutputBudget};
    return count;
}
std::vector<Receipt> analyze(const geo::Model &model, const GeometryMeshingOptions &o,
                             GeometryAnalysis &out) {
    validate(o);
    Work work{0, o};
    std::vector<Receipt> receipts;
    try {
        auto solids = model.solids();
        std::sort(solids.begin(), solids.end(), [](const auto *a, const auto *b) {
            return a && b ? a->tag() < b->tag() : bool(b);
        });
        if (solids.empty())
            throw Stop{GeometryReason::InvalidTopology};
        std::map<int, const GeometryBodyControl *> controls;
        for (const auto &c : o.bodyControls) {
            work.tick();
            if (c.solidTag < 0 || c.sourceFaceTag < -1 || !(c.targetSize >= 0) ||
                !std::isfinite(c.targetSize) || !(c.endToStartLayerWidthRatio > 0) ||
                !std::isfinite(c.endToStartLayerWidthRatio) ||
                !(c.maximumAdjacentLayerRatio >= 1) ||
                !std::isfinite(c.maximumAdjacentLayerRatio) ||
                static_cast<unsigned>(c.elements) > static_cast<unsigned>(GeometryElements::Tet4) ||
                (c.spatialField && !c.spatialField->valid()) ||
                (c.tetRefinement.enabled &&
                 (!(c.tetRefinement.maximumEdgeToTarget >= 1) ||
                  !std::isfinite(c.tetRefinement.maximumEdgeToTarget) ||
                  !(c.tetRefinement.maximumRadiusEdge >= 1) ||
                  !std::isfinite(c.tetRefinement.maximumRadiusEdge) ||
                  !c.tetRefinement.maximumInsertions || !c.tetRefinement.maximumWork)) ||
                (c.sizeProfile &&
                 (!c.sizeProfile->valid() || c.layers || c.endToStartLayerWidthRatio != 1)) ||
                (c.strategy && static_cast<unsigned>(*c.strategy) >=
                                   static_cast<unsigned>(GeometryStrategy::Unsupported)) ||
                !controls.emplace(c.solidTag, &c).second)
                throw Stop{GeometryReason::InvalidRequest};
        }
        std::size_t controlsUsed = 0;
        std::uint64_t retainedProfileBytes = 0;
        int previous = -1;
        for (const auto *s : solids) {
            if (!s || s->tag() < 0 || s->tag() == previous)
                throw Stop{GeometryReason::InvalidTopology};
            previous = s->tag();
            auto r = inspect(model, *s, work);
            const auto control = controls.find(s->tag());
            if (control != controls.end()) {
                r.control = *control->second;
                ++controlsUsed;
            }
            r.report.targetSize = r.control.targetSize > 0 ? r.control.targetSize : o.targetSize;
            r.report.elements = r.control.elements;
            r.requestedStrategy = r.control.strategy.value_or(o.strategy);
            if (!r.control.spatialField)
                r.control.spatialField = o.spatialField;
            if (r.control.spatialField) {
                if (r.control.sizeProfile || r.control.layers || r.control.sourceFaceTag >= 0 ||
                    r.control.endToStartLayerWidthRatio != 1 ||
                    (r.requestedStrategy != GeometryStrategy::Auto &&
                     r.requestedStrategy != GeometryStrategy::BoundaryConformingTet4))
                    throw Stop{GeometryReason::NoMatchingStrategy};
                const double h = r.control.spatialField->minimumOnBox(r.lo, r.hi);
                if (!(h > 0) || !std::isfinite(h))
                    throw Stop{GeometryReason::SpatialFieldRejected};
                r.report.usedSpatialField = true;
                r.report.spatialMinimumBound = std::min(h, r.report.targetSize);
                r.requestedStrategy = GeometryStrategy::BoundaryConformingTet4;
            }
            if (r.control.tetRefinement.enabled) {
                if (r.requestedStrategy == GeometryStrategy::Auto)
                    r.requestedStrategy = GeometryStrategy::BoundaryConformingTet4;
                if (r.requestedStrategy != GeometryStrategy::BoundaryConformingTet4)
                    throw Stop{GeometryReason::NoMatchingStrategy};
            }
            cancelled(o);
            for (const auto &other : receipts) {
                work.tick();
                const double t = std::max(r.tolerance, other.tolerance);
                const bool separated = r.hi.x < other.lo.x - t || other.hi.x < r.lo.x - t ||
                                       r.hi.y < other.lo.y - t || other.hi.y < r.lo.y - t ||
                                       r.hi.z < other.lo.z - t || other.hi.z < r.lo.z - t;
                if (!separated)
                    throw Stop{GeometryReason::UnresolvedBodyInterface};
            }
            const bool sweep = recognizeSweep(r, work);
            if (r.control.sizeProfile && !sweep)
                throw Stop{GeometryReason::NoMatchingStrategy};
            if (sweep) {
                if (r.control.sizeProfile) {
                    const auto &profile = *r.control.sizeProfile;
                    const double extent = length(r.report.sweep);
                    if (std::abs(profile.length() - extent) > extent * 1e-10)
                        throw Stop{GeometryReason::InvalidRequest};
                    auto bound = r.control.profileLimits;
                    if (!(bound.maximumAdjacentGrowth == 0 ||
                          (std::isfinite(bound.maximumAdjacentGrowth) &&
                           bound.maximumAdjacentGrowth >= 1)))
                        throw Stop{GeometryReason::InvalidRequest};
                    bound.maximumAdjacentGrowth =
                        std::min(r.control.maximumAdjacentLayerRatio,
                                 bound.maximumAdjacentGrowth == 0 ? profile.maximumGrowth()
                                                                  : bound.maximumAdjacentGrowth);
                    bound.maximumCells =
                        std::min({bound.maximumCells, o.maximumCells, o.maximumNodes,
                                  static_cast<std::uint64_t>(UINT32_MAX)});
                    bound.maximumBytes = std::min(bound.maximumBytes, o.maximumProfileWorkingBytes -
                                                                          retainedProfileBytes);
                    if (o.cancellation)
                        bound.cancellation = o.cancellation;
                    if (r.control.profileLimits.cancellation &&
                        r.control.profileLimits.cancellation->load(std::memory_order_relaxed))
                        throw Stop{GeometryReason::Cancelled};
                    auto planned = sizing::planProfileDistribution(profile, bound);
                    r.report.profileReason = planned.reason;
                    r.report.profileWorkingBytes = planned.workingBytes;
                    if (!planned.planned) {
                        out.bodies.push_back(r.report);
                        throw Stop{planned.reason == sizing::ProfileReason::Cancelled
                                       ? GeometryReason::Cancelled
                                   : planned.reason == sizing::ProfileReason::AllocationFailed
                                       ? GeometryReason::AllocationFailed
                                       : GeometryReason::ProfileRejected};
                    }
                    r.report.usedSizeProfile = true;
                    r.report.layers = static_cast<std::uint32_t>(planned.positions.size() - 1);
                    r.report.minimumLayerWidth = planned.minimumWidth;
                    r.report.maximumLayerWidth = planned.maximumWidth;
                    r.report.maximumAdjacentLayerGrowth = planned.maximumAdjacentGrowth;
                    r.report.maximumLayerWidthToTarget = planned.maximumWidthToTarget;
                    r.report.plateauCellsRemoved = planned.plateauCellsRemoved;
                    // Normalized arrays are shared with the mapped engine and
                    // moved into the sweep request: no extra retained copy.
                    r.profileFractions = std::move(planned.positions);
                    retainedProfileBytes = add(retainedProfileBytes,
                                               mul(r.profileFractions.capacity(), sizeof(double)));
                    if (retainedProfileBytes > o.maximumProfileWorkingBytes)
                        throw Stop{GeometryReason::ProfileRejected};
                    for (auto &x : r.profileFractions)
                        x /= profile.length();
                    r.profileFractions.front() = 0;
                    r.profileFractions.back() = 1;
                } else {
                    r.report.layers = r.control.layers
                                          ? r.control.layers
                                          : steps(length(r.report.sweep), r.report.targetSize);
                    r.report.endToStartLayerWidthRatio = r.control.endToStartLayerWidthRatio;
                    if (r.report.layers == 1 && r.control.endToStartLayerWidthRatio != 1)
                        throw Stop{GeometryReason::InvalidRequest};
                }
                // Validate grading without allocating a potentially huge schedule.
                const double logRatio = std::abs(std::log(r.control.endToStartLayerWidthRatio));
                if (r.report.layers > 1 && logRatio / (r.report.layers - 1) >
                                               std::log(r.control.maximumAdjacentLayerRatio) +
                                                   32 * std::numeric_limits<double>::epsilon())
                    throw Stop{GeometryReason::InvalidRequest};
            }
            const bool block = sweep && recognizeBlock(r, o);
            r.report.strategy = r.requestedStrategy == GeometryStrategy::Auto
                                    ? (block   ? GeometryStrategy::MappedHex8
                                       : sweep ? GeometryStrategy::LinearPrismSweep
                                               : GeometryStrategy::BoundaryConformingTet4)
                                    : r.requestedStrategy;
            if ((r.report.strategy == GeometryStrategy::MappedHex8 && !block) ||
                (r.report.strategy == GeometryStrategy::LinearPrismSweep && !sweep))
                throw Stop{GeometryReason::NoMatchingStrategy};
            if (r.report.strategy == GeometryStrategy::BoundaryConformingTet4 &&
                (r.control.sourceFaceTag >= 0 || r.control.layers ||
                 r.control.endToStartLayerWidthRatio != 1 || r.control.sizeProfile))
                throw Stop{GeometryReason::NoMatchingStrategy};
            if (r.report.strategy == GeometryStrategy::BoundaryConformingTet4)
                r.report.layers = 0;
            if (r.report.strategy == GeometryStrategy::MappedHex8)
                estimateBlock(r, o);
            else if (r.report.strategy == GeometryStrategy::LinearPrismSweep) {
                const auto n = sourceBound(r, r.faces[r.source], o);
                const auto layers = r.report.layers;
                r.report.estimatedNodes = mul(n, add(layers, 1));
                r.report.estimatedCells =
                    add(mul(mul(n, 2), layers), add(mul(n, 4), mul(n, layers)));
                if (r.control.elements == GeometryElements::Tet4)
                    r.report.estimatedCells = mul(r.report.estimatedCells, 3);
                r.report.estimatedOutputBytes =
                    bytes(r.report.estimatedNodes, r.report.estimatedCells,
                          mul(r.report.estimatedCells, 6));
                // Conservative source estimates cap preparation, not the final
                // output. Exact output budgets are checked after source build.
            } else {
                std::uint64_t candidateBound = 0;
                for (const auto &f : r.faces)
                    candidateBound = add(candidateBound, sourceBound(r, f, o));
                if (candidateBound > o.maximumSourcePoints)
                    throw Stop{GeometryReason::OutputBudget};
                const auto d = r.hi - r.lo;
                const double h =
                    r.report.usedSpatialField ? r.report.spatialMinimumBound : r.report.targetSize;
                const auto n =
                    mul(mul(add(steps(d.x, h), 2), add(steps(d.y, h), 2)), add(steps(d.z, h), 2));
                if (n > o.maximumNodes)
                    throw Stop{GeometryReason::OutputBudget};
            }
            out.bodies.push_back(r.report);
            receipts.push_back(std::move(r));
        }
        if (controlsUsed != controls.size())
            throw Stop{GeometryReason::InvalidRequest};
        out.work = work.count;
        out.complete = true;
        out.reason = GeometryReason::None;
        return receipts;
    } catch (...) {
        out.work = work.count;
        throw;
    }
}

mesh::Mesh sweepBody(Receipt &r, const GeometryMeshingOptions &o) {
    const auto &face = r.faces[r.source];
    View source;
    source.original = r.view.original;
    source.fs = {face.face};
    for (int v : face.vertices)
        source.vs.push_back(source.vertex(v));
    for (int e : face.edges)
        source.es.push_back(source.edge(e));
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(r.report.targetSize);
    spec.threads = o.threads;
    spec.elementsPerCircle = 0;
    mesh::Mesh surface;
    surface.setDedupTolerance(r.tolerance * .01);
    meshEdges(source, spec, surface);
    const auto built = meshFacesConstrained(source, spec, surface, {true, o.maximumMeshingWork});
    r.report.engineStatus = static_cast<std::uint32_t>(built.status);
    r.report.engineReason = static_cast<std::uint32_t>(built.reason);
    r.report.enginePhase = static_cast<std::uint32_t>(built.phase);
    r.report.recoveryKind = static_cast<std::uint32_t>(built.recoveryKind);
    r.report.recoveryDetail = built.recoveryDetail;
    if (!built.published() && built.planarRecovery.reason == PlanarRecoveryReason::WorkBudget)
        throw Stop{GeometryReason::MeshingBudget};
    if (!built.published())
        throw Stop{GeometryReason::SourceMeshingFailed};
    if (static_cast<std::uint64_t>(surface.nodeCount()) > o.maximumSourcePoints)
        throw Stop{GeometryReason::OutputBudget};
    cancelled(o);
    std::vector<Index> map(static_cast<std::size_t>(surface.nodeCount()), kInvalidIndex);
    mesh::Mesh cap;
    for (Index i = 0; i < surface.cellCount(); ++i)
        if (surface.cell(i).type == mesh::CellType::Tri3) {
            const auto *c = surface.cellNodes(i);
            for (int j = 0; j < 3; ++j)
                map[static_cast<std::size_t>(c[j])] = 0;
        }
    for (Index i = 0; i < surface.nodeCount(); ++i)
        if (map[static_cast<std::size_t>(i)] == 0)
            map[static_cast<std::size_t>(i)] = cap.addNode(surface.node(i), surface.nodeClass(i));
    for (Index i = 0; i < surface.cellCount(); ++i)
        if (surface.cell(i).type == mesh::CellType::Tri3) {
            const auto *c = surface.cellNodes(i);
            std::array<Index, 3> ids{map[static_cast<std::size_t>(c[0])],
                                     map[static_cast<std::size_t>(c[1])],
                                     map[static_cast<std::size_t>(c[2])]};
            if ((cap.node(ids[1]) - cap.node(ids[0]))
                    .cross(cap.node(ids[2]) - cap.node(ids[0]))
                    .dot(r.report.sweep) < 0)
                std::swap(ids[1], ids[2]);
            (void)cap.addCell(mesh::CellType::Tri3, ids, {2, face.face->tag()});
        }
    PreparedMixedSweepRequest request;
    request.layers = r.report.layers;
    request.cellToCellGrowth =
        request.layers == 1
            ? 1
            : std::exp(std::log(r.control.endToStartLayerWidthRatio) / (request.layers - 1));
    request.maximumAllowedGrowth = r.control.maximumAdjacentLayerRatio;
    request.layerFractions = std::move(r.profileFractions);
    request.direction = r.report.sweep;
    request.totalThickness = length(r.report.sweep);
    request.executionThreads = o.threads;
    request.fullAudit = false;
    request.maximumNodes = o.maximumNodes;
    request.maximumCells = o.maximumCells;
    request.volumeClassification = {3, r.report.solidTag};
    request.targetClassification = {2, r.report.targetFaceTag};
    request.sideClassification = {2, r.sideFace.begin()->second};
    for (Index i = 0; i < cap.nodeCount(); ++i) {
        const auto cls = cap.nodeClass(i);
        mesh::Classification top = request.targetClassification,
                             middle = request.volumeClassification;
        if (cls.dim == 0) {
            top = {0, r.topVertex.at(cls.tag)};
            middle = {1, r.verticalEdge.at(cls.tag)};
        } else if (cls.dim == 1) {
            top = {1, r.topEdge.at(cls.tag)};
            middle = {2, r.sideFace.at(cls.tag)};
        }
        request.targetNodeClassifications.push_back(top);
        request.intermediateNodeClassifications.push_back(middle);
    }
    for (Index i = 0; i < surface.cellCount(); ++i)
        if (surface.cell(i).type == mesh::CellType::Line2) {
            const auto *c = surface.cellNodes(i);
            const auto a = map[static_cast<std::size_t>(c[0])],
                       b = map[static_cast<std::size_t>(c[1])];
            if (a < 0 || b < 0)
                throw Stop{GeometryReason::SourceMeshingFailed};
            request.sideEdgeClassifications.push_back(
                {a, b, {2, r.sideFace.at(surface.cell(i).onEntity.tag)}});
        }
    const auto n = static_cast<std::uint64_t>(cap.nodeCount()),
               tri = static_cast<std::uint64_t>(cap.cellCount()),
               edge = static_cast<std::uint64_t>(request.sideEdgeClassifications.size()),
               layers = static_cast<std::uint64_t>(request.layers);
    const auto nodes = mul(n, add(layers, 1)), wedges = mul(tri, layers), quads = mul(edge, layers),
               triangles = mul(tri, 2), cells = add(wedges, add(quads, triangles));
    limits(nodes, cells,
           bytes(nodes, cells, add(mul(wedges, 6), add(mul(quads, 4), mul(triangles, 3)))), o);
    constructionWorkFloor(nodes, cells, o);
    mesh::Mesh result;
    const auto swept = generatePreparedMixedSweep(cap, request, result);
    r.report.engineStatus = static_cast<std::uint32_t>(swept.status);
    r.report.engineReason = static_cast<std::uint32_t>(swept.reason);
    r.report.enginePhase = static_cast<std::uint32_t>(swept.phase);
    if (!swept.published())
        throw Stop{GeometryReason::VolumeMeshingFailed};
    return result;
}
class CappedBodyField final : public sizing::BoundedSizeField {
  public:
    CappedBodyField(const sizing::BoundedSizeField *f, double cap) : field_(f), cap_(cap) {}
    bool valid() const noexcept override { return !field_ || field_->valid(); }
    std::uint64_t boxQueryWorkBound() const noexcept override {
        return field_ ? field_->boxQueryWorkBound() : 1;
    }
    double minimumOnBox(const Vec3 &lo, const Vec3 &hi) const noexcept override {
        const double h = field_ ? field_->minimumOnBox(lo, hi) : cap_;
        return (h > 0 && std::isfinite(h)) ? std::min(h, cap_)
                                           : std::numeric_limits<double>::quiet_NaN();
    }
    double eval(const Vec3 &p) const override {
        const double h = field_ ? field_->eval(p) : cap_;
        return (h > 0 && std::isfinite(h)) ? std::min(h, cap_)
                                           : std::numeric_limits<double>::quiet_NaN();
    }

  private:
    const sizing::BoundedSizeField *field_;
    double cap_;
};
mesh::Mesh generate(Receipt &r, const GeometryMeshingOptions &o) {
    mesh::Mesh result;
    if (r.report.strategy == GeometryStrategy::MappedHex8) {
        r.mapped.geometry.maxNodes = o.maximumNodes;
        r.mapped.geometry.maxCells = o.maximumCells;
        limits(r.report.estimatedNodes, r.report.estimatedCells, r.report.estimatedOutputBytes, o);
        constructionWorkFloor(r.report.estimatedNodes, r.report.estimatedCells, o);
        const auto made = generateMappedHexBlock(r.mapped, result);
        r.report.engineStatus = static_cast<std::uint32_t>(made.status);
        r.report.engineReason = static_cast<std::uint32_t>(made.reason);
        r.report.enginePhase = static_cast<std::uint32_t>(made.phase);
        if (!made.published())
            throw Stop{GeometryReason::VolumeMeshingFailed};
    } else if (r.report.strategy == GeometryStrategy::LinearPrismSweep)
        result = sweepBody(r, o);
    else {
        auto ownership = buildModelBoundaryOwnershipProvider(r.view);
        if (!ownership.built)
            throw Stop{GeometryReason::InvalidTopology};
        PreOcctResourceLimits bound;
        bound.maximumNodes = o.maximumNodes;
        bound.maximumTets = o.maximumCells;
        bound.maximumWorkUnits = o.maximumMeshingWork;
        bound.maximumReportRows = std::min<std::uint64_t>(250000, o.maximumMeshingWork);
        auto request = makePreOcctTetRequest(r.control.tetRefinement.enabled
                                                 ? PreOcctTetOutput::RefinedTet4
                                                 : PreOcctTetOutput::BoundaryConformingTet4,
                                             ownership.provider.stamp(), {}, bound);
        auto& conforming=request.pipeline.source.source.conformingRequest;
        conforming.usePlanarSurfaceRecovery=o.usePlanarSurfaceRecovery;
        conforming.edgeGeneralRecovery=o.edgeGeneralRecovery;
        if(conforming.edgeGeneralRecovery.enabled) {
            conforming.edgeGeneralRecovery.maximumTotalWorkUnits=std::min(
                conforming.edgeGeneralRecovery.maximumTotalWorkUnits,o.maximumMeshingWork);
            conforming.edgeGeneralRecovery.limits.maxWorkUnits=std::min(
                conforming.edgeGeneralRecovery.limits.maxWorkUnits,o.maximumMeshingWork);
        }
        if (r.control.tetRefinement.enabled) {
            auto &policy = request.pipeline.source.source.refinement;
            policy.maxEdgeToTarget = r.control.tetRefinement.maximumEdgeToTarget;
            policy.maxRadiusEdge = r.control.tetRefinement.maximumRadiusEdge;
            policy.maxInsertions =
                std::min(r.control.tetRefinement.maximumInsertions, o.maximumNodes);
            policy.maxWorkUnits =
                std::min(r.control.tetRefinement.maximumWork, o.maximumMeshingWork);
            policy.maxReportTopologyHandles = add(policy.maxWorkUnits, policy.maxInsertions);
        }
        sizing::MeshSpec spec;
        if (r.control.spatialField)
            spec.field =
                std::make_shared<CappedBodyField>(r.control.spatialField, r.report.targetSize);
        else
            spec.field = std::make_shared<sizing::ConstantField>(r.report.targetSize);
        spec.threads = o.threads;
        const auto made = meshPreparedTetModel(r.view, spec, request, result, ownership.provider);
        r.report.engineStatus = static_cast<std::uint32_t>(made.status);
        r.report.engineReason = static_cast<std::uint32_t>(made.reason);
        r.report.enginePhase = static_cast<std::uint32_t>(made.phase);
        const auto copySurfaceDiagnostic = [&](const SurfaceSetOutcome *receipt) {
            if (!receipt)
                return;
            const auto& surface = *receipt;
            r.report.hasSurfaceDiagnostic = true;
            r.report.surfaceStatus = static_cast<std::uint32_t>(surface.setStatus);
            r.report.surfaceReason = static_cast<std::uint32_t>(surface.setFailure.reason);
            r.report.surfaceProcessedFaces = surface.processed.size();
            r.report.surfaceUnattemptedFaces = surface.unattempted.size();
            for (const auto& row : surface.processed)
                if (row.status != SurfaceFaceStatus::Ready) {
                    r.report.surfaceFailureFacePersistentId = row.facePersistentId;
                    break;
                }
            const auto& failure = surface.setFailure;
            if (failure.bndReject != bnd::SteinerReject::None) {
                r.report.recoveryKind =
                    static_cast<std::uint32_t>(ProductionSurfaceRecoveryKind::Rejected);
                r.report.recoveryDetail = static_cast<std::uint32_t>(failure.bndReject);
            } else if (failure.bndRefuse != bnd::SteinerRefuse::None) {
                r.report.recoveryKind =
                    static_cast<std::uint32_t>(ProductionSurfaceRecoveryKind::Refused);
                r.report.recoveryDetail = static_cast<std::uint32_t>(failure.bndRefuse);
            } else if (failure.bndFail != bnd::SteinerFail::None) {
                r.report.recoveryKind =
                    static_cast<std::uint32_t>(ProductionSurfaceRecoveryKind::Failed);
                r.report.recoveryDetail = static_cast<std::uint32_t>(failure.bndFail);
            }
        };
        if (made.conforming) {
            r.report.hasBoundaryDiagnostic = true;
            r.report.boundaryStatus = static_cast<std::uint32_t>(made.conforming->status);
            r.report.boundaryReason = static_cast<std::uint32_t>(made.conforming->failure.reason);
            r.report.boundaryPhase = static_cast<std::uint32_t>(made.conforming->failure.phase);
            copySurfaceDiagnostic(made.conforming->surfaceOutcome.get());
        }
        if (made.refined) {
            r.report.hasRefinementDiagnostic = true;
            r.report.refinementStatus = static_cast<std::uint32_t>(made.refined->status);
            r.report.refinementReason = static_cast<std::uint32_t>(made.refined->failure.reason);
            r.report.refinementPhase = static_cast<std::uint32_t>(made.refined->failure.phase);
            if (made.refined->retainedResult)
                r.report.refinementInsertions =
                    made.refined->retainedResult->metrics.committedInsertions;
            if (made.refined->sourceReceipt) {
                r.report.hasBoundaryDiagnostic = true;
                r.report.boundaryStatus =
                    static_cast<std::uint32_t>(made.refined->sourceReceipt->status);
                r.report.boundaryReason =
                    static_cast<std::uint32_t>(made.refined->sourceReceipt->failure.reason);
                r.report.boundaryPhase =
                    static_cast<std::uint32_t>(made.refined->sourceReceipt->failure.phase);
                copySurfaceDiagnostic(made.refined->sourceReceipt->surfaceOutcome.get());
            }
        }
        if (!made.published() && made.conforming &&
            made.conforming->failure.edgeOutcome.reason==bnd::TetEdgeRecoveryReason::ResourceBudgetReached)
            throw Stop{GeometryReason::MeshingBudget};
        if (!made.published())
            throw Stop{GeometryReason::VolumeMeshingFailed};
    }
    if (r.control.elements == GeometryElements::Tet4 &&
        r.report.strategy != GeometryStrategy::BoundaryConformingTet4) {
        ConvexToTetOptions conversion;
        conversion.maximumNodes = o.maximumNodes;
        conversion.maximumCells = o.maximumCells;
        conversion.maximumOutputBytes = o.maximumOutputBytes;
        conversion.maximumWorkingBytes = o.maximumConversionWorkingBytes;
        conversion.auditLimits = {o.maximumAuditFaceUses, o.maximumAuditTemporaryBytes,
                                  o.cancellation};
        conversion.cancellation = o.cancellation;
        conversion.minimumQuality = o.minimumElementQuality;
        conversion.warningQuality = o.warningElementQuality;
        mesh::Mesh converted;
        r.report.conversion = convertConvexCellsToTet4(result, conversion, converted);
        if (!r.report.conversion.published) {
            const auto reason = r.report.conversion.reason;
            throw Stop{reason == ConvexToTetReason::Cancelled ? GeometryReason::Cancelled
                       : reason == ConvexToTetReason::AllocationFailed
                           ? GeometryReason::AllocationFailed
                       : reason == ConvexToTetReason::OutputBudget ? GeometryReason::OutputBudget
                       : reason == ConvexToTetReason::QualityRequirement
                           ? GeometryReason::QualityRequirementNotMet
                           : GeometryReason::ConversionFailed};
        }
        result.swap(converted);
    }
    limits(static_cast<std::uint64_t>(result.nodeCount()),
           static_cast<std::uint64_t>(result.cellCount()), bytes(result), o);
    if (r.control.spatialField || r.control.tetRefinement.enabled ||
        o.warningElementSizeRatio > 0 || o.maximumElementSizeRatio > 0) {
        MeshSizeAuditOptions audit;
        audit.sizeCap = r.report.targetSize;
        audit.warningEdgeToTargetRatio = o.warningElementSizeRatio;
        audit.maximumEdgeToTargetRatio = o.maximumElementSizeRatio;
        audit.maximumCells = o.maximumCells;
        audit.maximumSizingWork = o.maximumSizingAuditWork;
        audit.cancellation = o.cancellation;
        const CappedBodyField field(r.control.spatialField, r.report.targetSize);
        r.report.sizeAudit = auditMeshSize(result, field, audit);
        if (!r.report.sizeAudit.complete || !r.report.sizeAudit.requirementMet)
            throw Stop{r.report.sizeAudit.reason == MeshSizeAuditReason::Cancelled
                           ? GeometryReason::Cancelled
                       : r.report.sizeAudit.reason == MeshSizeAuditReason::SizeRequirement
                           ? GeometryReason::SizeRequirementNotMet
                       : r.report.sizeAudit.reason == MeshSizeAuditReason::BudgetExceeded
                           ? GeometryReason::MeshingBudget
                           : GeometryReason::SpatialFieldRejected};
    }
    return result;
}
void append(mesh::Mesh &dst, const mesh::Mesh &src) {
    const Index offset = dst.nodeCount();
    dst.reserveAppend(static_cast<std::size_t>(src.nodeCount()),
                      static_cast<std::size_t>(src.cellCount()),
                      static_cast<std::size_t>(src.connectivityCount()));
    for (Index i = 0; i < src.nodeCount(); ++i)
        (void)dst.addNode(src.node(i), src.nodeClass(i));
    for (Index i = 0; i < src.cellCount(); ++i) {
        const auto &c = src.cell(i);
        const int n = mesh::nodesPerCell(c.type);
        std::array<Index, 10> ids{};
        for (int j = 0; j < n; ++j)
            ids[static_cast<std::size_t>(j)] = src.cellNodes(i)[j] + offset;
        (void)dst.addCell(c.type, std::span<const Index>(ids.data(), static_cast<std::size_t>(n)),
                          c.onEntity);
    }
}
} // namespace

GeometryAnalysis analyzeGeometryStrategies(const geo::Model &model,
                                           const GeometryMeshingOptions &options) noexcept {
    GeometryAnalysis out;
    try {
        (void)analyze(model, options, out);
    } catch (const Stop &s) {
        out.reason = s.reason;
    } catch (const std::bad_alloc &) {
        out.reason = GeometryReason::AllocationFailed;
    } catch (...) {
        out.reason = GeometryReason::InternalError;
    }
    return out;
}
GeometryMeshingOutcome meshGeometryModel(const geo::Model &model,
                                         const GeometryMeshingOptions &options,
                                         mesh::Mesh &destination) noexcept {
    GeometryMeshingOutcome out;
    std::size_t active = 0;
    try {
        const auto start = Clock::now();
        auto receipts = analyze(model, options, out.analysis);
        out.analysisSeconds = seconds(start);
        mesh::Mesh staged;
        std::uint64_t sizingAuditWork = 0;
        for (active = 0; active < receipts.size(); ++active) {
            cancelled(options);
            auto remaining = options;
            remaining.maximumNodes -= static_cast<std::uint64_t>(staged.nodeCount());
            remaining.maximumCells -= static_cast<std::uint64_t>(staged.cellCount());
            remaining.maximumOutputBytes -= bytes(staged);
            remaining.maximumSizingAuditWork -= sizingAuditWork;
            const auto begin = Clock::now();
            mesh::Mesh body;
            try {
                body = generate(receipts[active], remaining);
            } catch (...) {
                receipts[active].report.generationSeconds = seconds(begin);
                out.generationSeconds += receipts[active].report.generationSeconds;
                out.analysis.bodies[active] = receipts[active].report;
                throw;
            }
            auto &report = receipts[active].report;
            sizingAuditWork += report.sizeAudit.sizingWork;
            report.generationSeconds = seconds(begin);
            report.outputNodes = static_cast<std::uint64_t>(body.nodeCount());
            report.outputCells = static_cast<std::uint64_t>(body.cellCount());
            report.outputBytes = bytes(body);
            out.generationSeconds += report.generationSeconds;
            out.analysis.bodies[active] = report;
            if (staged.nodeCount() == 0)
                staged.swap(body);
            else
                append(staged, body);
        }
        cancelled(options);
        out.outputBytes = bytes(staged);
        if (options.fullBoundaryAudit) {
            const auto begin = Clock::now();
            out.audit = mesh::auditLinearVolumeBoundary(staged, {options.maximumAuditFaceUses,
                                                                 options.maximumAuditTemporaryBytes,
                                                                 options.cancellation});
            out.auditSeconds = seconds(begin);
            if (!out.audit.valid)
                throw Stop{out.audit.cancelled          ? GeometryReason::Cancelled
                           : out.audit.allocationFailed ? GeometryReason::AllocationFailed
                           : out.audit.budgetExceeded   ? GeometryReason::OutputBudget
                                                        : GeometryReason::BoundaryAuditFailed};
        }
        cancelled(options);
        if (options.minimumElementQuality > 0 || options.warningElementQuality > 0) {
            const auto begin = Clock::now();
            out.quality = analyzeMixedMeshQuality(staged);
            out.qualitySeconds = seconds(begin);
            if (out.quality.allocationFailed)
                throw Stop{GeometryReason::AllocationFailed};
            out.qualityAssessed = true;
            out.qualityWarning = out.quality.minimumQuality < options.warningElementQuality;
            if (!out.quality.valid || out.quality.minimumQuality < options.minimumElementQuality)
                throw Stop{GeometryReason::QualityRequirementNotMet};
        }
        cancelled(options);
        destination.swap(staged);
        out.published = true;
        out.destinationUnchanged = false;
        out.reason = GeometryReason::None;
    } catch (const Stop &s) {
        out.reason = s.reason;
    } catch (const std::bad_alloc &) {
        out.reason = GeometryReason::AllocationFailed;
    } catch (...) {
        out.reason = GeometryReason::InternalError;
    }
    if (!out.published && active < out.analysis.bodies.size())
        out.analysis.bodies[active].reason = out.reason;
    return out;
}
const char *geometryStrategyName(GeometryStrategy s) noexcept {
    switch (s) {
    case GeometryStrategy::Auto:
        return "auto";
    case GeometryStrategy::MappedHex8:
        return "mapped_hex8";
    case GeometryStrategy::LinearPrismSweep:
        return "linear_prism_sweep";
    case GeometryStrategy::BoundaryConformingTet4:
        return "constrained_tet4";
    case GeometryStrategy::Unsupported:
        return "unsupported";
    }
    return "unknown";
}
const char *geometryReasonName(GeometryReason r) noexcept {
    switch (r) {
    case GeometryReason::None:
        return "none";
    case GeometryReason::InvalidRequest:
        return "invalid_request";
    case GeometryReason::InvalidTopology:
        return "invalid_topology";
    case GeometryReason::NonfiniteGeometry:
        return "nonfinite_geometry";
    case GeometryReason::UnsupportedChart:
        return "unsupported_chart";
    case GeometryReason::NoMatchingStrategy:
        return "no_matching_strategy";
    case GeometryReason::UnresolvedBodyInterface:
        return "unresolved_body_interface";
    case GeometryReason::AnalysisBudget:
        return "analysis_budget";
    case GeometryReason::MeshingBudget:
        return "meshing_budget";
    case GeometryReason::OutputBudget:
        return "output_budget";
    case GeometryReason::SourceMeshingFailed:
        return "source_meshing_failed";
    case GeometryReason::VolumeMeshingFailed:
        return "volume_meshing_failed";
    case GeometryReason::BoundaryAuditFailed:
        return "boundary_audit_failed";
    case GeometryReason::QualityRequirementNotMet:
        return "quality_requirement_not_met";
    case GeometryReason::Cancelled:
        return "cancelled";
    case GeometryReason::AllocationFailed:
        return "allocation_failed";
    case GeometryReason::InternalError:
        return "internal_error";
    case GeometryReason::ProfileRejected:
        return "profile_rejected";
    case GeometryReason::ConversionFailed:
        return "conversion_failed";
    case GeometryReason::SpatialFieldRejected:
        return "spatial_field_rejected";
    case GeometryReason::SizeRequirementNotMet:
        return "size_requirement_not_met";
    }
    return "unknown";
}
} // namespace paacs::algo
