// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/CurvedPrismLayerMesher.hpp"

#include "paacs/core/Parallel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <new>
#include <numeric>
#include <vector>

namespace paacs::algo {
namespace {
struct Triangle {
    Vec3 a, b, c;
    std::array<Index, 3> nodes{};
    Vec3 centroid{};
};
struct Bounds {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    void add(const Vec3 &p) {
        lo.x = std::min(lo.x, p.x);
        lo.y = std::min(lo.y, p.y);
        lo.z = std::min(lo.z, p.z);
        hi.x = std::max(hi.x, p.x);
        hi.y = std::max(hi.y, p.y);
        hi.z = std::max(hi.z, p.z);
    }
};
struct BvhNode {
    Bounds bounds{};
    std::uint32_t begin = 0, end = 0, left = 0, right = 0;
    bool leaf = true;
};
bool finite(const Vec3 &p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
bool validClass(mesh::Classification c, int dim) noexcept {
    return c.dim == dim && c.tag >= 0;
}
double totalThickness(const CurvedPrismLayerRequest &r) noexcept {
    long double w = r.firstLayerThickness, total = 0;
    for (std::uint32_t i = 0; i < r.layers; ++i) {
        total += w;
        w *= r.cellToCellGrowth;
    }
    return static_cast<double>(total);
}

bool extractTriangles(const mesh::Mesh &mesh, std::vector<Triangle> &triangles,
                      std::vector<std::vector<Index>> *nodeTriangles,
                      std::vector<std::vector<Index>> *adjacency) {
    if (nodeTriangles)
        nodeTriangles->assign(static_cast<std::size_t>(mesh.nodeCount()), {});
    if (adjacency)
        adjacency->assign(static_cast<std::size_t>(mesh.nodeCount()), {});
    triangles.reserve(static_cast<std::size_t>(mesh.cellCount()));
    for (Index ci = 0; ci < mesh.cellCount(); ++ci) {
        const auto &cell = mesh.cell(ci);
        if (cell.type != mesh::CellType::Tri3)
            return false;
        const Index *n = mesh.cellNodes(ci);
        if (n[0] < 0 || n[1] < 0 || n[2] < 0 || n[0] >= mesh.nodeCount() ||
            n[1] >= mesh.nodeCount() || n[2] >= mesh.nodeCount() || n[0] == n[1] || n[1] == n[2] ||
            n[2] == n[0])
            return false;
        const Vec3 a = mesh.node(n[0]), b = mesh.node(n[1]), c = mesh.node(n[2]);
        if (!finite(a) || !finite(b) || !finite(c))
            return false;
        if (!((b - a).cross(c - a).norm2() > 0.0))
            return false;
        triangles.push_back({a, b, c, {n[0], n[1], n[2]}, (a + b + c) / 3.0});
        const Index ti = static_cast<Index>(triangles.size() - 1);
        if (nodeTriangles)
            for (Index node : {n[0], n[1], n[2]})
                (*nodeTriangles)[static_cast<std::size_t>(node)].push_back(ti);
        if (adjacency) {
            const std::array<std::array<Index, 2>, 3> edges{
                {{n[0], n[1]}, {n[1], n[2]}, {n[2], n[0]}}};
            for (const auto &e : edges) {
                (*adjacency)[static_cast<std::size_t>(e[0])].push_back(e[1]);
                (*adjacency)[static_cast<std::size_t>(e[1])].push_back(e[0]);
            }
        }
    }
    if (adjacency)
        for (auto &row : *adjacency) {
            std::sort(row.begin(), row.end());
            row.erase(std::unique(row.begin(), row.end()), row.end());
        }
    return !triangles.empty();
}
Bounds triangleBounds(const Triangle &t) {
    Bounds b;
    b.add(t.a);
    b.add(t.b);
    b.add(t.c);
    return b;
}
Bounds merge(Bounds a, const Bounds &b) {
    a.add(b.lo);
    a.add(b.hi);
    return a;
}
std::uint32_t buildBvh(std::vector<Triangle> &tris, std::vector<BvhNode> &nodes,
                       std::uint32_t begin, std::uint32_t end) {
    const std::uint32_t id = static_cast<std::uint32_t>(nodes.size());
    nodes.push_back({});
    Bounds bounds, centers;
    for (std::uint32_t i = begin; i < end; ++i) {
        bounds = merge(bounds, triangleBounds(tris[i]));
        centers.add(tris[i].centroid);
    }
    nodes[id].bounds = bounds;
    nodes[id].begin = begin;
    nodes[id].end = end;
    if (end - begin <= 8)
        return id;
    const Vec3 span = centers.hi - centers.lo;
    int axis = span.y > span.x ? 1 : 0;
    if ((axis == 0 ? span.x : span.y) < span.z)
        axis = 2;
    const std::uint32_t middle = begin + (end - begin) / 2U;
    std::nth_element(tris.begin() + begin, tris.begin() + middle, tris.begin() + end,
                     [&](const Triangle &a, const Triangle &b) {
                         return axis == 0   ? a.centroid.x < b.centroid.x
                                : axis == 1 ? a.centroid.y < b.centroid.y
                                            : a.centroid.z < b.centroid.z;
                     });
    nodes[id].leaf = false;
    nodes[id].left = buildBvh(tris, nodes, begin, middle);
    nodes[id].right = buildBvh(tris, nodes, middle, end);
    return id;
}
bool rayBox(const Vec3 &o, const Vec3 &d, const Bounds &b, double maxT) noexcept {
    double lo = 0.0, hi = maxT;
    const auto slab = [&](double origin, double direction, double mn, double mx) {
        if (std::abs(direction) < 1e-30)
            return origin >= mn && origin <= mx;
        double a = (mn - origin) / direction, bv = (mx - origin) / direction;
        if (a > bv)
            std::swap(a, bv);
        lo = std::max(lo, a);
        hi = std::min(hi, bv);
        return lo <= hi;
    };
    return slab(o.x, d.x, b.lo.x, b.hi.x) && slab(o.y, d.y, b.lo.y, b.hi.y) &&
           slab(o.z, d.z, b.lo.z, b.hi.z);
}
double rayTriangle(const Vec3 &o, const Vec3 &d, const Triangle &t) noexcept {
    const Vec3 e1 = t.b - t.a, e2 = t.c - t.a, p = d.cross(e2);
    const double det = e1.dot(p);
    if (std::abs(det) < 1e-14)
        return (std::numeric_limits<double>::infinity)();
    const double inv = 1.0 / det;
    const Vec3 s = o - t.a;
    const double u = s.dot(p) * inv;
    if (u < 0 || u > 1)
        return (std::numeric_limits<double>::infinity)();
    const Vec3 q = s.cross(e1);
    const double v = d.dot(q) * inv;
    if (v < 0 || u + v > 1)
        return (std::numeric_limits<double>::infinity)();
    const double distance = e2.dot(q) * inv;
    return distance > 1e-12 ? distance : (std::numeric_limits<double>::infinity)();
}
bool incident(Index node, const Triangle &t) noexcept {
    return t.nodes[0] == node || t.nodes[1] == node || t.nodes[2] == node;
}
double nearestHit(Index node, const Vec3 &o, const Vec3 &d, const std::vector<Triangle> &tris,
                  const std::vector<BvhNode> &bvh, bool skipIncident, double limit) {
    if (bvh.empty())
        return limit;
    double best = limit;
    std::array<std::uint32_t, 128> stack{};
    std::size_t count = 0;
    stack[count++] = 0;
    while (count) {
        const auto id = stack[--count];
        const auto &n = bvh[id];
        if (!rayBox(o, d, n.bounds, best))
            continue;
        if (n.leaf) {
            for (std::uint32_t i = n.begin; i < n.end; ++i) {
                if (skipIncident && incident(node, tris[i]))
                    continue;
                best = std::min(best, rayTriangle(o, d, tris[i]));
            }
        } else {
            if (count + 2 > stack.size())
                return best;
            stack[count++] = n.left;
            stack[count++] = n.right;
        }
    }
    return best;
}
} // namespace

CurvedPrismOutcome generateCurvedPrismLayers(const mesh::Mesh &source,
                                             const CurvedPrismLayerRequest &request,
                                             mesh::Mesh &destination) noexcept {
    CurvedPrismOutcome out;
    if (source.nodeCount() == 0 || source.cellCount() == 0) {
        out.reason = CurvedPrismReason::EmptySource;
        return out;
    }
    if (request.layers == 0 || !std::isfinite(request.firstLayerThickness) ||
        request.firstLayerThickness <= 0 || !std::isfinite(request.cellToCellGrowth) ||
        request.cellToCellGrowth < 1 || request.cellToCellGrowth > request.maximumAllowedGrowth ||
        request.maximumAllowedGrowth < 1 || request.maximumTotalThickness <= 0 ||
        request.minimumTotalThickness <= 0 ||
        request.minimumTotalThickness > request.maximumTotalThickness ||
        request.collisionSafetyFactor <= 0 || request.collisionSafetyFactor >= 1 ||
        request.maximumAdjacentThicknessRatio < 1 || request.globalBackoffFactor <= 0 ||
        request.globalBackoffFactor >= 1) {
        out.reason = CurvedPrismReason::InvalidLayerControl;
        return out;
    }
    if (!validClass(request.volumeClassification, 3) ||
        !validClass(request.targetClassification, 2) ||
        !validClass(request.sideClassification, 2)) {
        out.reason = CurvedPrismReason::InvalidClassification;
        return out;
    }
    try {
        std::vector<Triangle> sourceTriangles;
        std::vector<std::vector<Index>> nodeTriangles, adjacency;
        if (!extractTriangles(source, sourceTriangles, &nodeTriangles, &adjacency)) {
            out.reason = CurvedPrismReason::UnsupportedSource;
            return out;
        }
        std::vector<Triangle> obstacles;
        const mesh::Mesh *collision = request.collisionSurface ? request.collisionSurface : &source;
        if (!extractTriangles(*collision, obstacles, nullptr, nullptr)) {
            out.reason = CurvedPrismReason::InvalidCollisionSurface;
            return out;
        }
        if (obstacles.size() > request.maximumCollisionTriangles) {
            out.reason = CurvedPrismReason::CollisionBudgetExceeded;
            return out;
        }
        std::vector<BvhNode> bvh;
        bvh.reserve(obstacles.size() * 2U);
        (void)buildBvh(obstacles, bvh, 0, static_cast<std::uint32_t>(obstacles.size()));
        const std::size_t nodeCount = static_cast<std::size_t>(source.nodeCount());
        std::vector<Vec3> normals(nodeCount, Vec3{});
        std::vector<Vec3> nextNormals(nodeCount);
        std::vector<std::uint8_t> feature(nodeCount, 0);
        const double cosFeature =
            std::cos(request.featureAngleDegrees * 3.14159265358979323846 / 180.0);
        for (std::size_t node = 0; node < nodeCount; ++node) {
            Vec3 sum{};
            std::vector<Vec3> local;
            for (Index ti : nodeTriangles[node]) {
                const auto &t = sourceTriangles[static_cast<std::size_t>(ti)];
                const Vec3 n = (t.b - t.a).cross(t.c - t.a);
                sum = sum + n;
                local.push_back(n.normalized());
            }
            if (!(sum.norm2() > 0)) {
                out.reason = CurvedPrismReason::InvalidConnectivity;
                return out;
            }
            normals[node] = sum.normalized();
            for (std::size_t i = 0; i < local.size(); ++i)
                for (std::size_t j = i + 1; j < local.size(); ++j)
                    if (local[i].dot(local[j]) < cosFeature)
                        feature[node] = 1;
        }
        for (std::uint32_t iteration = 0; iteration < request.normalSmoothingIterations;
             ++iteration) {
            nextNormals = normals;
            for (std::size_t node = 0; node < nodeCount; ++node) {
                if (feature[node])
                    continue;
                Vec3 sum = normals[node];
                for (Index other : adjacency[node])
                    sum = sum + normals[static_cast<std::size_t>(other)];
                if (sum.norm2() > 0)
                    nextNormals[node] = sum.normalized();
            }
            normals.swap(nextNormals);
        }
        const double requested = std::min(totalThickness(request), request.maximumTotalThickness);
        if (!std::isfinite(requested) || requested < request.minimumTotalThickness) {
            out.reason = CurvedPrismReason::LayerCollapsed;
            return out;
        }
        std::vector<double> thickness(nodeCount, requested), nextThickness(nodeCount);
        const bool skipIncident = collision == &source;
        core::parallelForWithThreads(0, nodeCount, request.executionThreads, [&](std::size_t node) {
            const double hit = nearestHit(
                static_cast<Index>(node), source.node(static_cast<Index>(node)), normals[node],
                obstacles, bvh, skipIncident, requested / request.collisionSafetyFactor);
            if (std::isfinite(hit))
                thickness[node] = std::min(requested, hit * request.collisionSafetyFactor);
        });
        for (std::uint32_t it = 0; it < request.thicknessSmoothingIterations; ++it) {
            nextThickness = thickness;
            for (std::size_t node = 0; node < nodeCount; ++node)
                for (Index other : adjacency[node])
                    nextThickness[node] =
                        std::min(nextThickness[node], thickness[static_cast<std::size_t>(other)] *
                                                          request.maximumAdjacentThicknessRatio);
            thickness.swap(nextThickness);
        }
        out.metrics.requestedThickness = requested;
        out.metrics.minimumThickness = *std::min_element(thickness.begin(), thickness.end());
        out.metrics.maximumThickness = *std::max_element(thickness.begin(), thickness.end());
        out.metrics.clippedNodes = std::count_if(thickness.begin(), thickness.end(), [&](double t) {
            return t < requested * (1 - 1e-12);
        });
        out.metrics.featureNodes = std::count(feature.begin(), feature.end(), 1);
        out.metrics.rayQueries = nodeCount;
        if (out.metrics.minimumThickness < request.minimumTotalThickness) {
            out.reason = CurvedPrismReason::LayerCollapsed;
            return out;
        }
        std::vector<Vec3> targets(nodeCount);
        for (std::size_t n = 0; n < nodeCount; ++n)
            targets[n] = source.node(static_cast<Index>(n)) + normals[n] * thickness[n];
        PreparedMixedSweepRequest sweep;
        sweep.layers = request.layers;
        sweep.cellToCellGrowth = request.cellToCellGrowth;
        sweep.maximumAllowedGrowth = request.maximumAllowedGrowth;
        sweep.targetNodes = targets;
        sweep.volumeClassification = request.volumeClassification;
        sweep.targetClassification = request.targetClassification;
        sweep.sideClassification = request.sideClassification;
        sweep.emitBoundaryFaces = request.emitBoundaryFaces;
        sweep.executionThreads = request.executionThreads;
        for (std::uint32_t backoff = 0; backoff <= request.maximumGlobalBackoffSteps; ++backoff) {
            mesh::Mesh candidate;
            const auto result = generatePreparedMixedSweep(source, sweep, candidate);
            if (result.published()) {
                candidate.swap(destination);
                out.published = true;
                out.reason = CurvedPrismReason::None;
                out.destinationUnchanged = false;
                out.metrics.backoffSteps = backoff;
                out.metrics.sourceNodes = source.nodeCount();
                out.metrics.sourceTriangles = sourceTriangles.size();
                out.metrics.collisionTriangles = obstacles.size();
                out.metrics.sweep = result.metrics;
                return out;
            }
            for (std::size_t n = 0; n < targets.size(); ++n) {
                thickness[n] *= request.globalBackoffFactor;
                targets[n] = source.node(static_cast<Index>(n)) + normals[n] * thickness[n];
            }
            if (*std::min_element(thickness.begin(), thickness.end()) <
                request.minimumTotalThickness)
                break;
            sweep.targetNodes = targets;
        }
        out.reason = CurvedPrismReason::SweepFailed;
        return out;
    } catch (const std::bad_alloc &) {
        out.reason = CurvedPrismReason::AllocationFailed;
    } catch (...) {
        out.reason = CurvedPrismReason::InternalError;
    }
    return out;
}

} // namespace paacs::algo
