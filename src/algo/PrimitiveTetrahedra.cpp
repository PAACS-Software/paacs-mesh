// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "PrimitiveMesherInternal.hpp"
#include "paacs/algo/Tet10FastCertificate.hpp"
#include "paacs/algo/Tet10JacobianCertificate.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/CertifiedMeshStorage.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/topo/ElementTraits.hpp"
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace paacs::algo::primitive_detail {
namespace {
using Clock = std::chrono::steady_clock;
using Storage = mesh::detail::CertifiedMeshStorage;
using Type = mesh::CellType;
constexpr std::uint64_t chunk = 1024;
struct Stop {
    PrimitiveReason reason;
};
void need(bool ok, PrimitiveReason reason) {
    if (!ok)
        throw Stop{reason};
}
void cancel(const PrimitiveMeshingOptions &o) {
    need(!(o.cancellation && o.cancellation->load(std::memory_order_relaxed)) &&
             !(o.auditLimits.cancellation &&
               o.auditLimits.cancellation->load(std::memory_order_relaxed)),
         PrimitiveReason::Cancelled);
}
std::uint64_t add(std::uint64_t a, std::uint64_t b) {
    need(b <= UINT64_MAX - a, PrimitiveReason::OutputBudget);
    return a + b;
}
std::uint64_t mul(std::uint64_t a, std::uint64_t b) {
    need(a == 0 || b <= UINT64_MAX / a, PrimitiveReason::OutputBudget);
    return a * b;
}
double elapsed(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now() - t).count();
}
struct Counts {
    std::uint64_t n = 0, t = 0, b = 0;
};
std::uint64_t payload(Counts c, bool quadratic = false) {
    return add(add(mul(c.n, sizeof(Vec3) + sizeof(mesh::Classification)),
                   mul(add(c.t, c.b), sizeof(mesh::Cell))),
               mul(add(mul(c.t, quadratic ? 10 : 4), mul(c.b, quadratic ? 6 : 3)), sizeof(Index)));
}
void outputBudget(Counts c, bool quadratic, const PrimitiveMeshingOptions &o) {
    need(c.n <= o.maximumNodes && c.n <= INT32_MAX && add(c.t, c.b) <= o.maximumCells &&
             payload(c, quadratic) <= o.maximumOutputBytes && payload(c, quadratic) <= SIZE_MAX,
         PrimitiveReason::OutputBudget);
}
void workBudget(std::uint64_t bytes, const PrimitiveMeshingOptions &o, PrimitiveMeshingResult &r) {
    r.plannedWorkingBytes = std::max(r.plannedWorkingBytes, bytes);
    need(bytes <= o.maximumWorkingBytes && bytes <= SIZE_MAX, PrimitiveReason::WorkingBudget);
}
std::uint64_t auditBudget(Counts c, const PrimitiveMeshingOptions &o) {
    const auto faces = add(mul(c.t, 4), c.b),
               bytes = add(mul(faces, 16), mul(c.n, sizeof(std::size_t)));
    need(faces <= o.auditLimits.maximumFaceUses && bytes <= o.auditLimits.maximumTemporaryBytes,
         PrimitiveReason::WorkingBudget);
    return bytes;
}
Vec3 midpoint(Vec3 a, Vec3 b) {
    return {std::midpoint(a.x, b.x), std::midpoint(a.y, b.y), std::midpoint(a.z, b.z)};
}
std::array<Index, 4> quad(const Index *ids, const topo::LocalFace &f) {
    return {ids[f.v[0]], ids[f.v[1]], ids[f.v[2]], ids[f.v[3]]};
}
std::array<std::array<Index, 3>, 2> triangles(std::array<Index, 4> q) {
    const auto first = static_cast<unsigned>(std::min_element(q.begin(), q.end()) - q.begin());
    return {{{q[first], q[(first + 1) % 4], q[(first + 2) % 4]},
             {q[first], q[(first + 2) % 4], q[(first + 3) % 4]}}};
}
// Pulling triangulation uses ONE global vertex ordering on every shared face.
// Six tetrahedra per hex; two per pyramid. No new nodes and no face hash table.
mesh::Mesh convert(const mesh::Mesh &src, const PrimitiveMeshingResult &base,
                   const PrimitiveMeshingOptions &o, PrimitiveMeshingResult &r, Counts &c) {
    c = {static_cast<std::uint64_t>(src.nodeCount()),
         add(mul(base.hexes, 6), mul(base.pyramids, 2)),
         add(mul(base.boundaryQuads, 2), base.boundaryTriangles)};
    outputBudget(c, false, o);
    workBudget(add(base.outputBytes, auditBudget(c, o)), o, r);
    Storage storage(c.n, c.t + c.b, 4 * c.t + 3 * c.b);
    for (Index i = 0; i < src.nodeCount(); ++i)
        storage.setNode(i, src.node(i), src.nodeClass(i));
    const auto volumes = base.hexes + base.pyramids;
    core::parallelChunksWithThreads(
        0, volumes, chunk, o.threads, [&](std::size_t, std::size_t lo, std::size_t hi) {
            cancel(o);
            for (std::size_t i = lo; i < hi; ++i) {
                const auto *ids = src.cellNodes(i);
                const auto type =
                    i < base.hexes ? topo::ElementType::Hex8 : topo::ElementType::Pyramid5;
                const auto &tr = topo::traits(type);
                const Index anchor = *std::min_element(ids, ids + tr.cornerCount);
                std::uint64_t ordinal =
                                  i < base.hexes ? 6 * i : 6 * base.hexes + 2 * (i - base.hexes),
                              emitted = 0;
                for (unsigned f = 0; f < tr.faceCount; ++f) {
                    const auto &face = tr.faces[f];
                    bool incident = false;
                    for (unsigned j = 0; j < face.n; ++j)
                        incident |= ids[face.v[j]] == anchor;
                    if (incident)
                        continue;
                    if (face.n == 3) {
                        storage.setCell(ordinal, 4 * ordinal, Type::Tet4,
                                        std::array<Index, 4>{anchor, ids[face.v[0]], ids[face.v[1]],
                                                             ids[face.v[2]]},
                                        {3, 1});
                        ++ordinal;
                        ++emitted;
                    } else
                        for (const auto &t : triangles(quad(ids, face))) {
                            storage.setCell(ordinal, 4 * ordinal, Type::Tet4,
                                            std::array<Index, 4>{anchor, t[0], t[1], t[2]}, {3, 1});
                            ++ordinal;
                            ++emitted;
                        }
                }
                need(emitted == (i < base.hexes ? 6u : 2u), PrimitiveReason::InternalError);
            }
        });
    core::parallelChunksWithThreads(
        0, base.boundaryQuads + base.boundaryTriangles, chunk, o.threads,
        [&](std::size_t, std::size_t lo, std::size_t hi) {
            cancel(o);
            for (std::size_t j = lo; j < hi; ++j) {
                const Index i = static_cast<Index>(volumes + j);
                const auto *ids = src.cellNodes(i);
                if (j < base.boundaryQuads) {
                    const auto ts = triangles({ids[0], ids[1], ids[2], ids[3]});
                    for (unsigned k = 0; k < 2; ++k)
                        storage.setCell(c.t + 2 * j + k, 4 * c.t + 3 * (2 * j + k), Type::Tri3,
                                        ts[k], src.cell(i).onEntity);
                } else {
                    const auto at = 2 * base.boundaryQuads + j - base.boundaryQuads;
                    storage.setCell(c.t + at, 4 * c.t + 3 * at, Type::Tri3,
                                    std::array<Index, 3>{ids[0], ids[1], ids[2]},
                                    src.cell(i).onEntity);
                }
            }
        });
    mesh::Mesh result;
    storage.publishTo(result);
    return result;
}
struct EdgeUse {
    std::uint64_t key, slot;
};
struct EdgeTable {
    std::vector<EdgeUse> uses;
    std::vector<Index> binding;
    std::uint64_t unique = 0;
};
std::uint64_t edgeKey(Index a, Index b) {
    if (a > b)
        std::swap(a, b);
    return (std::uint64_t(a) << 32) | std::uint64_t(b);
}
std::uint64_t edgeWork(Counts c) {
    return mul(add(mul(c.t, 6), mul(c.b, 3)), sizeof(EdgeUse) + sizeof(Index));
}
EdgeTable edgeTable(const mesh::Mesh &m, Counts c, const PrimitiveMeshingOptions &o) {
    EdgeTable e;
    const auto uses = 6 * c.t + 3 * c.b;
    e.uses.resize(uses);
    e.binding.resize(uses);
    const auto &te = topo::traits(topo::ElementType::Tet4);
    const auto &tri = topo::traits(topo::ElementType::Tri3);
    core::parallelChunksWithThreads(
        0, c.t + c.b, chunk, o.threads, [&](std::size_t, std::size_t lo, std::size_t hi) {
            cancel(o);
            for (std::size_t i = lo; i < hi; ++i) {
                const auto &tr = i < c.t ? te : tri;
                const auto *ids = m.cellNodes(i);
                const auto start = i < c.t ? 6 * i : 6 * c.t + 3 * (i - c.t);
                for (unsigned j = 0; j < tr.edgeCount; ++j)
                    e.uses[start + j] = {edgeKey(ids[tr.edges[j].a], ids[tr.edges[j].b]),
                                         start + j};
            }
        });
    // In-place introsort: bounded auxiliary stack, no per-edge allocation.
    std::sort(e.uses.begin(), e.uses.end(),
              [](auto a, auto b) { return a.key < b.key || (a.key == b.key && a.slot < b.slot); });
    cancel(o);
    std::uint64_t previous = UINT64_MAX;
    for (std::size_t i = 0; i < e.uses.size(); ++i) {
        if (i % chunk == 0)
            cancel(o);
        const auto &u = e.uses[i];
        if (u.key != previous) {
            ++e.unique;
            previous = u.key;
        }
        e.binding[u.slot] = static_cast<Index>(c.n + e.unique - 1);
    }
    return e;
}
mesh::Classification edgeClass(const mesh::Mesh &m, Counts c, const EdgeTable &e, std::size_t begin,
                               std::size_t end) {
    int first = 0, second = 0;
    for (auto i = begin; i < end; ++i)
        if (e.uses[i].slot >= 6 * c.t) {
            const auto cell = c.t + (e.uses[i].slot - 6 * c.t) / 3;
            const auto tag = m.cell(static_cast<Index>(cell)).onEntity.tag;
            if (first == 0)
                first = tag;
            else if (tag != first) {
                need(second == 0 || second == tag, PrimitiveReason::BoundaryAuditFailed);
                second = tag;
            }
        }
    if (!first)
        return {3, 1};
    if (!second)
        return {2, first};
    if (first > second)
        std::swap(first, second);
    if (first == 1 && (second == 2 || second == 3))
        return {1, 99 + second};
    if ((first == 2 || first == 3) && second == 4)
        return {1, 101 + first};
    throw Stop{PrimitiveReason::BoundaryAuditFailed};
}
struct Frame {
    Vec3 x, y, z;
    explicit Frame(const PrimitiveShape &s) {
        z = s.axis / std::hypot(s.axis.x, s.axis.y, s.axis.z);
        auto hint = s.xDirection / std::hypot(s.xDirection.x, s.xDirection.y, s.xDirection.z);
        x = hint - z * hint.dot(z);
        x = x / x.norm();
        y = z.cross(x);
    }
};
// Six triangular sectors. Ring k has 6*k equally spaced angular vertices at
// radius k/m. There is no separate square core/collar. Each swept triangle is
// split into three conforming tetrahedra using the same section vertex order.
mesh::Mesh uniformSweep(const PrimitiveShape &s, const PrimitiveMeshingOptions &o,
                        PrimitiveMeshingResult &r, Counts &c) {
    constexpr double pi = std::numbers::pi;
    const bool circular = s.kind == PrimitiveKind::Torus || s.kind == PrimitiveKind::Elbow;
    const bool closed = s.kind == PrimitiveKind::Torus;
    const double sweep = closed ? 2 * pi : s.sweepAngle;
    const double tolerance = o.maximumChordError > 0 ? o.maximumChordError : .0025 * s.radius;
    const double scale =
        std::max({s.radius, circular ? s.majorRadius + s.radius : s.height, std::abs(s.origin.x),
                  std::abs(s.origin.y), std::abs(s.origin.z)});
    const double roundoff = 128 * std::numeric_limits<double>::epsilon() * scale;
    need(tolerance > roundoff, PrimitiveReason::NumericallyUnresolved);
    auto count = [](double value) -> std::uint64_t {
        need(std::isfinite(value) && value >= 0 && value <= UINT32_MAX - 2.,
             PrimitiveReason::OutputBudget);
        return static_cast<std::uint64_t>(std::max(1., std::ceil(value)));
    };
    // 1-sinc(theta) <= theta^2/6 bounds the polygonal volume deficit.
    const auto sizeRings = count(s.radius / o.targetSize);
    std::uint64_t geometryRings;
    if (circular) {
        // L >= sweep*R*m/r gives major step <= beta/m, beta=r/R.
        // sin(x)^2 <= x^2 bounds the SUM of both sagittas by
        // r/m^2 * (pi^2/72 + beta*(1+beta)/8). This avoids imposing
        // an arbitrary half-tolerance on each circle and overrefining both.
        const double beta = s.radius / s.majorRadius;
        const auto chordRings = count(std::sqrt((s.radius / (tolerance - roundoff)) *
                                                (pi * pi / 72 + beta * (1 + beta) / 8)));
        // 1-sinc(a)*sinc(b) <= (a*a+b*b)/6 for these positive small steps.
        const auto volumeRings =
            count(std::sqrt((pi * pi / 9 + beta * beta) / (6 * o.maximumRelativeVolumeError)));
        geometryRings = std::max(chordRings, volumeRings);
    } else {
        // Preserve the rc.6 cylinder counts and hashes exactly.
        const double angle = std::min(
            4 * std::asin(std::sqrt(std::min(.5, (tolerance - roundoff) / (2 * s.radius)))),
            std::sqrt(6 * o.maximumRelativeVolumeError));
        geometryRings = count(pi / (3 * angle));
    }
    auto m = std::max(sizeRings, geometryRings);
    m = add(m, m % 2); // multiples of 12 boundary vertices; exact quarter count
    const double spacing = s.radius / static_cast<double>(m);
    const auto sizeLayers =
        !o.axialPositions.empty()
            ? o.axialPositions.size() - 1
            : count(circular ? sweep * (s.majorRadius / spacing) : s.height / spacing);
    const auto layers = std::max(sizeLayers, std::uint64_t(closed ? 8 : 1));
    const auto planes = closed ? layers : add(layers, 1);
    const auto sectionNodes = add(1, mul(3, mul(m, add(m, 1))));
    const auto sectionTriangles = mul(6, mul(m, m));
    const auto boundaryEdges = mul(6, m);
    const auto capTriangles = closed ? 0 : sectionTriangles;
    c = {mul(sectionNodes, planes), mul(3, mul(sectionTriangles, layers)),
         add(mul(2, capTriangles), mul(2, mul(boundaryEdges, layers)))};
    outputBudget(c, false, o);
    const auto sectionBytes = add(
        add(mul(sectionNodes, sizeof(Vec2)), mul(sectionTriangles, sizeof(std::array<Index, 3>))),
        circular ? mul(planes, sizeof(Vec2)) : 0);
    workBudget(sectionBytes, o, r);
    workBudget(auditBudget(c, o), o, r);
    if (o.elements == PrimitiveElements::Tet10) {
        // E=N+T+B/2-chi; chi is 1 for a ball, 0 for a solid torus.
        const Counts quadratic{add(mul(c.n, 2), add(c.t, c.b / 2)) - (closed ? 0 : 1), c.t, c.b};
        outputBudget(quadratic, true, o);
        workBudget(add(payload(c), edgeWork(c)), o, r);
    }
    if (circular) {
        r.uniformCircular = true;
        r.circularGeometryLimited = geometryRings > sizeRings;
        r.circularRadialSpacing = spacing;
        const auto chordFactor = 2 * std::sin(sweep / (2 * static_cast<double>(layers)));
        r.circularInnerSpacing = (s.majorRadius - s.radius) * chordFactor;
        r.circularOuterSpacing = (s.majorRadius + s.radius) * chordFactor;
    } else {
        r.uniformCylinder = o.axialPositions.empty();
        r.cylinderGeometryLimited = geometryRings > sizeRings;
        r.cylinderRadialSpacing = spacing;
        r.cylinderAxialSpacing = s.height / static_cast<double>(layers);
    }
    r.radialLayers = static_cast<std::uint32_t>(m);
    r.sweepLayers = static_cast<std::uint32_t>(layers);
    r.segmentsPerQuarter = static_cast<std::uint32_t>(boundaryEdges / 4);
    r.chordErrorBound =
        2 * s.radius * std::pow(std::sin(pi / (2 * static_cast<double>(boundaryEdges))), 2) +
        roundoff;
    if (circular)
        r.chordErrorBound += 2 * (s.majorRadius + s.radius) *
                             std::pow(std::sin(sweep / (4 * static_cast<double>(layers))), 2);
    need(r.chordErrorBound <= tolerance, PrimitiveReason::GeometryTolerance);
    r.analyticVolume = circular ? pi * s.radius * s.radius * s.majorRadius * sweep
                                : pi * s.radius * s.radius * s.height;
    need(r.analyticVolume > 0 && std::isfinite(r.analyticVolume),
         PrimitiveReason::NumericallyUnresolved);
    const auto ringNode = [](std::uint64_t k, std::uint64_t j) -> Index {
        return k == 0 ? 0 : static_cast<Index>(1 + 3 * k * (k - 1) + j % (6 * k));
    };
    std::vector<Vec2> section(sectionNodes);
    std::vector<std::array<Index, 3>> faces(sectionTriangles);
    std::uint64_t at = 0;
    for (std::uint64_t k = 1; k <= m; ++k) {
        cancel(o);
        const double radius = s.radius * (static_cast<double>(k) / static_cast<double>(m));
        for (std::uint64_t j = 0; j < 6 * k; ++j) {
            if (j % chunk == 0)
                cancel(o);
            const double a = (2 * pi) * (static_cast<double>(j) / static_cast<double>(6 * k));
            section[ringNode(k, j)] = {radius * std::cos(a), radius * std::sin(a)};
        }
        for (std::uint64_t sector = 0; sector < 6; ++sector)
            for (std::uint64_t j = 0; j < k; ++j) {
                if (j % chunk == 0)
                    cancel(o);
                faces[at++] = {ringNode(k, sector * k + j), ringNode(k, sector * k + j + 1),
                               ringNode(k - 1, sector * (k - 1) + j)};
                if (j + 1 < k)
                    faces[at++] = {ringNode(k - 1, sector * (k - 1) + j),
                                   ringNode(k, sector * k + j + 1),
                                   ringNode(k - 1, sector * (k - 1) + j + 1)};
            }
    }
    need(at == sectionTriangles, PrimitiveReason::InternalError);
    std::vector<Vec2> trig(circular ? planes : 0);
    for (std::uint64_t k = 0; k < trig.size(); ++k) {
        if (k % chunk == 0)
            cancel(o);
        const double a = sweep * (static_cast<double>(k) / static_cast<double>(layers));
        trig[k] = {std::cos(a), std::sin(a)};
    }
    const Frame frame(s);
    const auto outerStart = static_cast<std::uint64_t>(ringNode(m, 0));
    Storage out(c.n, c.t + c.b, 4 * c.t + 3 * c.b);
    core::parallelChunksWithThreads(
        0, c.n, chunk, o.threads, [&](std::size_t, std::size_t lo, std::size_t hi) {
            cancel(o);
            for (std::size_t i = lo; i < hi; ++i) {
                const auto level = i / sectionNodes, local = i % sectionNodes;
                const auto p = section[local];
                const bool side = local >= outerStart;
                const int cap = closed ? 0 : level == 0 ? 2 : level == layers ? 3 : 0;
                const mesh::Classification tag = side && cap ? mesh::Classification{1, 99 + cap}
                                                 : side      ? mesh::Classification{2, 1}
                                                 : cap       ? mesh::Classification{2, cap}
                                                             : mesh::Classification{3, 1};
                const Vec3 point =
                    circular ? s.origin +
                                   (frame.x * trig[level].u + frame.y * trig[level].v) *
                                       (s.majorRadius + p.u) -
                                   frame.z * p.v
                             : s.origin + frame.x * p.u + frame.y * p.v +
                                   frame.z * (o.axialPositions.empty()
                                                  ? s.height * (static_cast<double>(level) /
                                                                static_cast<double>(layers))
                                                  : o.axialPositions[level]);
                out.setNode(i, point, tag);
            }
        });
    core::parallelChunksWithThreads(
        0, sectionTriangles * layers, chunk, o.threads,
        [&](std::size_t, std::size_t lo, std::size_t hi) {
            cancel(o);
            for (std::size_t i = lo; i < hi; ++i) {
                const auto level = i / sectionTriangles;
                auto v = faces[i % sectionTriangles];
                bool odd = false;
                for (auto [a, b] : {std::pair{0, 1}, std::pair{1, 2}, std::pair{0, 1}})
                    if (v[a] > v[b]) {
                        std::swap(v[a], v[b]);
                        odd = !odd;
                    }
                for (auto &id : v)
                    id += static_cast<Index>(level * sectionNodes);
                // The section ordering (not layer IDs) selects the diagonal.
                // This also agrees across the periodic last-to-first seam.
                const auto n = static_cast<Index>(((level + 1) % planes) * sectionNodes) -
                               static_cast<Index>(level * sectionNodes);
                std::array<std::array<Index, 4>, 3> t{{{v[0], v[1], v[2], v[2] + n},
                                                       {v[0], v[1] + n, v[1], v[2] + n},
                                                       {v[0], v[0] + n, v[1] + n, v[2] + n}}};
                for (unsigned j = 0; j < 3; ++j) {
                    // Combinatorial orientation from the original CCW triangle,
                    // not a coordinate-based inversion repair.
                    if (odd)
                        std::swap(t[j][1], t[j][2]);
                    out.setCell(3 * i + j, 4 * (3 * i + j), Type::Tet4, t[j], {3, 1});
                }
            }
        });
    core::parallelChunksWithThreads(
        0, capTriangles + boundaryEdges * layers, chunk, o.threads,
        [&](std::size_t, std::size_t lo, std::size_t hi) {
            cancel(o);
            for (std::size_t i = lo; i < hi; ++i) {
                std::array<std::array<Index, 3>, 2> t{};
                if (i < capTriangles) {
                    const auto v = faces[i];
                    const auto top = static_cast<Index>(layers * sectionNodes);
                    t = {{{v[0], v[2], v[1]}, {v[0] + top, v[1] + top, v[2] + top}}};
                } else {
                    const auto side = i - capTriangles;
                    const auto level = side / boundaryEdges, j = side % boundaryEdges;
                    const auto a = ringNode(m, j) + static_cast<Index>(level * sectionNodes);
                    const auto b = ringNode(m, j + 1) + static_cast<Index>(level * sectionNodes);
                    const auto n = static_cast<Index>(((level + 1) % planes) * sectionNodes) -
                                   static_cast<Index>(level * sectionNodes);
                    // Preserve oriented boundary and the section-selected
                    // diagonal even when the upper IDs wrap back to zero.
                    t = a < b ? std::array<std::array<Index, 3>, 2>{{{a, b, b + n},
                                                                     {a, b + n, a + n}}}
                              : std::array<std::array<Index, 3>, 2>{
                                    {{b, b + n, a + n}, {b, a + n, a}}};
                }
                for (unsigned j = 0; j < 2; ++j)
                    out.setCell(c.t + 2 * i + j, 4 * c.t + 3 * (2 * i + j), Type::Tri3, t[j],
                                {2, i < capTriangles ? static_cast<int>(2 + j) : 1});
            }
        });
    mesh::Mesh result;
    out.publishTo(result);
    return result;
}
Vec3 projectCone(Vec3 p, mesh::Classification tag, const PrimitiveShape &s, const Frame &f) {
    if (tag.dim == 3)
        return p;
    const auto local = p - s.origin;
    double height = local.dot(f.z);
    if (tag.dim == 1 || tag.tag == 2)
        height = 0;
    double u = local.dot(f.x), v = local.dot(f.y);
    if (tag.dim == 1 || tag.tag == 1) {
        const double radius = s.radius * (1 - height / s.height), d = std::hypot(u, v);
        need(d > 0 && radius > 0, PrimitiveReason::NumericallyUnresolved);
        u *= radius / d;
        v *= radius / d;
    }
    return s.origin + f.x * u + f.y * v + f.z * height;
}
void nodesAndMidpoints(Storage &out, const mesh::Mesh &src, Counts c, const EdgeTable &e,
                       const PrimitiveMeshingOptions &o, const PrimitiveShape *cone) {
    for (std::uint64_t i = 0; i < c.n; ++i) {
        if (i % chunk == 0)
            cancel(o);
        out.setNode(i, src.node(static_cast<Index>(i)), src.nodeClass(static_cast<Index>(i)));
    }
    const Frame f(cone ? *cone : PrimitiveShape{});
    std::size_t nextPoll = 0;
    for (std::size_t i = 0; i < e.uses.size();) {
        if (i >= nextPoll) {
            cancel(o);
            nextPoll = i + chunk;
        }
        auto end = i + 1;
        while (end < e.uses.size() && e.uses[end].key == e.uses[i].key)
            ++end;
        const auto key = e.uses[i].key;
        const auto tag = edgeClass(src, c, e, i, end);
        auto point = midpoint(src.node(static_cast<Index>(key >> 32)),
                              src.node(static_cast<Index>(key & UINT32_MAX)));
        if (cone)
            point = projectCone(point, tag, *cone, f);
        out.setNode(e.binding[e.uses[i].slot], point, tag);
        i = end;
    }
}
struct Measure {
    double volume, quality;
};
Measure measure(Vec3 a, Vec3 b, Vec3 c, Vec3 d) {
    const auto ab = b - a, ac = c - a, ad = d - a;
    const double v = ab.cross(ac).dot(ad) / 6;
    const double sum =
        ab.norm2() + ac.norm2() + ad.norm2() + (c - b).norm2() + (d - b).norm2() + (d - c).norm2();
    return {v, v > 0 ? 12 * std::pow(3 * v, 2. / 3) / sum : 0};
}
// Four corner children and four cells in the central octahedron. Pick the
// best of its three diagonals from actual projected coordinates. Boundary
// triangles are identical for all three choices, so neighbours cannot crack.
std::array<std::array<unsigned, 4>, 8> children(const Storage &out,
                                                const std::array<Index, 10> &v) {
    std::array<std::array<unsigned, 4>, 8> best{
        {{0, 4, 6, 7}, {4, 1, 5, 8}, {6, 5, 2, 9}, {7, 8, 9, 3}}};
    constexpr unsigned poles[3][2]{{4, 9}, {5, 7}, {6, 8}};
    constexpr unsigned cycles[3][4]{{5, 6, 7, 8}, {4, 8, 9, 6}, {4, 5, 9, 7}};
    double score = -1;
    for (unsigned d = 0; d < 3; ++d) {
        auto candidate = best;
        double q = 1;
        for (unsigned j = 0; j < 4; ++j) {
            auto &t = candidate[4 + j];
            t = {poles[d][0], cycles[d][j], cycles[d][(j + 1) % 4], poles[d][1]};
            q = std::min(q, measure(out.node(v[t[0]]), out.node(v[t[1]]), out.node(v[t[2]]),
                                    out.node(v[t[3]]))
                                .quality);
        }
        if (q > score) {
            score = q;
            best = candidate;
        }
    }
    need(score > 0, PrimitiveReason::InvalidCell);
    return best;
}
mesh::Mesh refineCone(const mesh::Mesh &src, Counts &c, const PrimitiveShape &s,
                      const PrimitiveMeshingOptions &o, PrimitiveMeshingResult &r) {
    const Counts planned{add(mul(c.n, 2), add(c.t, c.b / 2)) - 1, mul(c.t, 8), mul(c.b, 4)};
    outputBudget(planned, false, o);
    workBudget(add(payload(c), edgeWork(c)), o, r);
    auto e = edgeTable(src, c, o);
    need(c.n + e.unique == planned.n, PrimitiveReason::BoundaryAuditFailed);
    Storage out(planned.n, planned.t + planned.b, 4 * planned.t + 3 * planned.b);
    nodesAndMidpoints(out, src, c, e, o, &s);
    core::parallelChunksWithThreads(
        0, c.t, chunk, o.threads, [&](std::size_t, std::size_t lo, std::size_t hi) {
            cancel(o);
            for (std::size_t i = lo; i < hi; ++i) {
                std::array<Index, 10> v{};
                std::copy_n(src.cellNodes(i), 4, v.begin());
                std::copy_n(e.binding.data() + 6 * i, 6, v.begin() + 4);
                const auto pattern = children(out, v);
                for (unsigned j = 0; j < 8; ++j) {
                    const auto &t = pattern[j];
                    const auto ordinal = 8 * i + j;
                    out.setCell(ordinal, 4 * ordinal, Type::Tet4,
                                std::array<Index, 4>{v[t[0]], v[t[1]], v[t[2]], v[t[3]]}, {3, 1});
                }
            }
        });
    core::parallelChunksWithThreads(
        0, c.b, chunk, o.threads, [&](std::size_t, std::size_t lo, std::size_t hi) {
            cancel(o);
            for (std::size_t i = lo; i < hi; ++i) {
                const auto *v = src.cellNodes(c.t + i);
                const auto *b = e.binding.data() + 6 * c.t + 3 * i;
                const std::array<std::array<Index, 3>, 4> t{{{v[0], b[0], b[2]},
                                                             {b[0], v[1], b[1]},
                                                             {b[2], b[1], v[2]},
                                                             {b[0], b[1], b[2]}}};
                for (unsigned j = 0; j < 4; ++j)
                    out.setCell(planned.t + 4 * i + j, 4 * planned.t + 3 * (4 * i + j), Type::Tri3,
                                t[j], src.cell(c.t + i).onEntity);
            }
        });
    mesh::Mesh result;
    out.publishTo(result);
    c = planned;
    return result;
}
mesh::Mesh gradedCone(const PrimitiveShape &s, const PrimitiveMeshingOptions &o,
                      PrimitiveMeshingResult &r, Counts &c) {
    const double tolerance = o.maximumChordError > 0 ? o.maximumChordError : .0025 * s.radius;
    // Choose a small macro fan and power-of-two refinement together. The
    // factor two bounds triangles spanning two adjacent radial levels.
    unsigned sectors = 0, levels = 0;
    double bestCost = std::numeric_limits<double>::infinity();
    for (unsigned candidate : {8u, 12u, 16u, 24u, 32u}) {
        const double chord =
            2 * s.radius * std::pow(std::sin(std::numbers::pi / (2 * candidate)), 2);
        const double refinement =
            std::max(std::hypot(s.height, s.radius) / o.targetSize,
                     2 * chord / std::min(tolerance, s.radius * o.maximumRelativeVolumeError / 4));
        if (!(refinement > 0) || !std::isfinite(refinement))
            continue;
        const double depth = std::max(0., std::ceil(std::log2(refinement)));
        if (depth > 20)
            continue;
        const double cost = candidate * std::exp2(3 * depth);
        if (cost < bestCost) {
            bestCost = cost;
            sectors = candidate;
            levels = static_cast<unsigned>(depth);
        }
    }
    need(sectors != 0, PrimitiveReason::OutputBudget);
    c = {sectors + 2, sectors, 2 * sectors};
    auto final = c;
    for (unsigned i = 0; i < levels; ++i) {
        const Counts next{add(mul(final.n, 2), add(final.t, final.b / 2)) - 1, mul(final.t, 8),
                          mul(final.b, 4)};
        outputBudget(next, false, o);
        workBudget(add(payload(final), edgeWork(final)), o, r);
        final = next;
    }
    // Depth zero still has a real output, and must obey every output budget.
    outputBudget(final, false, o);
    if (o.elements == PrimitiveElements::Tet10) {
        const Counts quadratic{add(mul(final.n, 2), add(final.t, final.b / 2)) - 1, final.t,
                               final.b};
        outputBudget(quadratic, true, o);
        workBudget(add(payload(final), edgeWork(final)), o, r);
    }
    workBudget(auditBudget(final, o), o, r);
    r.gradedCone = true;
    r.coneRefinementLevels = levels;
    r.segmentsPerQuarter = (sectors / 4) * (1u << levels);
    r.sweepLayers = 1u << levels;
    const Frame f(s);
    Storage initial(c.n, c.t + c.b, 4 * c.t + 3 * c.b);
    initial.setNode(0, s.origin, {2, 2});
    initial.setNode(1, s.origin + f.z * s.height, {0, 201});
    for (unsigned i = 0; i < sectors; ++i) {
        const double angle = 2 * std::numbers::pi * i / sectors;
        initial.setNode(
            i + 2, s.origin + (f.x * std::cos(angle) + f.y * std::sin(angle)) * s.radius, {1, 101});
        const Index a = i + 2, b = (i + 1) % sectors + 2;
        initial.setCell(i, 4 * i, Type::Tet4, std::array<Index, 4>{0, a, b, 1}, {3, 1});
        initial.setCell(c.t + 2 * i, 4 * c.t + 6 * i, Type::Tri3, std::array<Index, 3>{a, b, 1},
                        {2, 1});
        initial.setCell(c.t + 2 * i + 1, 4 * c.t + 6 * i + 3, Type::Tri3,
                        std::array<Index, 3>{0, b, a}, {2, 2});
    }
    mesh::Mesh mesh;
    initial.publishTo(mesh);
    for (unsigned i = 0; i < levels; ++i) {
        cancel(o);
        auto next = refineCone(mesh, c, s, o, r);
        mesh.swap(next);
    }
    const double scale = std::max(
        {s.radius, s.height, std::abs(s.origin.x), std::abs(s.origin.y), std::abs(s.origin.z)});
    const double roundoff = 128 * std::numeric_limits<double>::epsilon() * scale;
    r.chordErrorBound = 0;
    for (std::uint64_t i = c.t; i < c.t + c.b; ++i) {
        if (i % chunk == 0)
            cancel(o);
        if (mesh.cell(i).onEntity.tag != 1)
            continue;
        const auto *v = mesh.cellNodes(i);
        double maxRadius = 0, maxAngle = 0;
        std::array<Vec2, 3> uv{};
        std::array<bool, 3> apex{};
        for (unsigned j = 0; j < 3; ++j) {
            const auto p = mesh.node(v[j]) - s.origin;
            uv[j] = {p.dot(f.x), p.dot(f.y)};
            apex[j] = mesh.nodeClass(v[j]).dim == 0;
            maxRadius = std::max(maxRadius, std::hypot(uv[j].u, uv[j].v));
        }
        for (unsigned a = 0; a < 3; ++a)
            for (unsigned b = a + 1; b < 3; ++b)
                if (!apex[a] && !apex[b]) {
                    const double cross = uv[a].u * uv[b].v - uv[a].v * uv[b].u;
                    maxAngle = std::max(
                        maxAngle,
                        std::abs(std::atan2(cross, uv[a].u * uv[b].u + uv[a].v * uv[b].v)));
                }
        r.chordErrorBound =
            std::max(r.chordErrorBound, 2 * maxRadius * std::pow(std::sin(maxAngle / 4), 2));
    }
    r.chordErrorBound += roundoff;
    need(r.chordErrorBound <= tolerance * (1 + 1e-12), PrimitiveReason::GeometryTolerance);
    r.analyticVolume = std::numbers::pi * s.radius * s.radius * s.height / 3;
    return mesh;
}
struct Part {
    long double volume = 0;
    double minimumVolume = std::numeric_limits<double>::infinity();
    PrimitiveQualityReport quality{};
};
void validate(const mesh::Mesh &m, Counts c, const PrimitiveMeshingOptions &o,
              PrimitiveMeshingResult &r) {
    const auto partsCount = (c.t + chunk - 1) / chunk;
    workBudget(mul(partsCount, sizeof(Part)), o, r);
    std::vector<Part> parts(partsCount);
    core::parallelChunksWithThreads(
        0, c.t, chunk, o.threads, [&](std::size_t part, std::size_t lo, std::size_t hi) {
            cancel(o);
            auto &a = parts[part];
            for (std::size_t i = lo; i < hi; ++i) {
                const auto *v = m.cellNodes(i);
                const auto q = measure(m.node(v[0]), m.node(v[1]), m.node(v[2]), m.node(v[3]));
                need(q.volume > 0 && std::isfinite(q.volume) &&
                         predicates::orient3d(m.node(v[0]), m.node(v[1]), m.node(v[2]),
                                              m.node(v[3])) == predicates::Sign::Positive,
                     PrimitiveReason::InvalidCell);
                need(q.quality > 0 && std::isfinite(q.quality),
                     PrimitiveReason::NumericallyUnresolved);
                bool sharp = false;
                for (unsigned j = 0; j < 4; ++j)
                    sharp |= m.nodeClass(v[j]).dim <= 1;
                recordQuality(a.quality, static_cast<Index>(i), q.quality, sharp, o.warningQuality);
                a.volume += q.volume;
                a.minimumVolume = std::min(a.minimumVolume, q.volume);
            }
        });
    r.quality = {};
    long double volume = 0;
    r.minimumCellVolume = std::numeric_limits<double>::infinity();
    for (const auto &a : parts) {
        volume += a.volume;
        r.minimumCellVolume = std::min(r.minimumCellVolume, a.minimumVolume);
        mergeQuality(r.quality, a.quality);
    }
    r.meshVolume = static_cast<double>(volume);
    r.minimumQuality = r.quality.worstQuality;
    r.qualityWarning = r.quality.lowQualityCells != 0;
    r.relativeVolumeError = std::abs(r.meshVolume - r.analyticVolume) / r.analyticVolume;
    need(r.minimumQuality >= o.minimumQuality &&
             !(r.qualityWarning && o.qualityPolicy == PrimitiveQualityPolicy::Reject),
         PrimitiveReason::QualityRequirement);
    need(r.relativeVolumeError <= o.maximumRelativeVolumeError, PrimitiveReason::GeometryTolerance);
}
long double quadraticVolume(const std::array<Vec3, 10> &p) {
    const auto &traits = topo::traits(topo::ElementType::Tet10);
    long double point[10][3]{};
    for (unsigned i = 1; i < 10; ++i) {
        point[i][0] = static_cast<long double>(p[i].x) - p[0].x;
        point[i][1] = static_cast<long double>(p[i].y) - p[0].y;
        point[i][2] = static_cast<long double>(p[i].z) - p[0].z;
    }
    long double delta[6][3]{};
    for (unsigned e = 0; e < 6; ++e)
        for (unsigned a = 0; a < 3; ++a)
            delta[e][a] =
                point[4 + e][a] - (point[traits.edges[e].a][a] + point[traits.edges[e].b][a]) / 2;
    constexpr int derivative[4][3]{{-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    long double volume = 0;
    // Five-point degree-three tetrahedral quadrature, normalized weights sum 1.
    for (unsigned q = 0; q < 5; ++q) {
        long double l[4]{.25L, .25L, .25L, .25L};
        if (q) {
            for (auto &v : l)
                v = 1.L / 6;
            l[q - 1] = .5L;
        }
        long double j[3][3]{};
        for (unsigned a = 0; a < 3; ++a)
            for (unsigned d = 0; d < 3; ++d) {
                j[a][d] = point[d + 1][a];
                for (unsigned e = 0; e < 6; ++e) {
                    const auto edge = traits.edges[e];
                    j[a][d] +=
                        4 * delta[e][a] *
                        (l[edge.a] * derivative[edge.b][d] + l[edge.b] * derivative[edge.a][d]);
                }
            }
        const long double det = j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1]) -
                                j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0]) +
                                j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
        volume += (q ? .45L : -.8L) * det / 6;
    }
    return volume;
}
mesh::Mesh elevate(const mesh::Mesh &src, Counts &c, const PrimitiveMeshingOptions &o,
                   PrimitiveMeshingResult &r) {
    const auto sourceBytes = payload(c), edgesBytes = edgeWork(c);
    workBudget(add(sourceBytes, edgesBytes), o, r);
    auto e = edgeTable(src, c, o);
    const Counts next{add(c.n, e.unique), c.t, c.b};
    outputBudget(next, true, o);
    Storage storage(next.n, c.t + c.b, 10 * c.t + 6 * c.b);
    nodesAndMidpoints(storage, src, c, e, o, nullptr);
    core::parallelChunksWithThreads(
        0, c.t + c.b, chunk, o.threads, [&](std::size_t, std::size_t lo, std::size_t hi) {
            cancel(o);
            for (std::size_t i = lo; i < hi; ++i) {
                if (i < c.t) {
                    std::array<Index, 10> v{};
                    std::copy_n(src.cellNodes(i), 4, v.begin());
                    std::copy_n(e.binding.data() + 6 * i, 6, v.begin() + 4);
                    storage.setCell(i, 10 * i, Type::Tet10, v, src.cell(i).onEntity);
                } else {
                    std::array<Index, 6> v{};
                    std::copy_n(src.cellNodes(i), 3, v.begin());
                    std::copy_n(e.binding.data() + 6 * c.t + 3 * (i - c.t), 3, v.begin() + 3);
                    storage.setCell(i, 10 * c.t + 6 * (i - c.t), Type::Tri6, v,
                                    src.cell(i).onEntity);
                }
            }
        });
    mesh::Mesh result;
    storage.publishTo(result);
    // Audit actual elevated edge uses against the sorted complete source-edge
    // groups. Also require unique midside ids, distinct from all corners.
    for (std::size_t i = 0; i < e.uses.size(); ++i) {
        if (i % chunk == 0)
            cancel(o);
        const auto slot = e.uses[i].slot;
        const auto cell = slot < 6 * c.t ? slot / 6 : c.t + (slot - 6 * c.t) / 3;
        const unsigned corner = cell < c.t ? 4 : 3;
        const unsigned local =
            static_cast<unsigned>(slot < 6 * c.t ? slot % 6 : (slot - 6 * c.t) % 3);
        const auto id = result.cellNodes(cell)[corner + local];
        need(id == e.binding[slot] && id >= static_cast<Index>(c.n) && id < result.nodeCount(),
             PrimitiveReason::BoundaryAuditFailed);
    }
    r.quadraticConnectivityVerified = true;
    struct Integral {
        long double volume = 0;
        double minimum = std::numeric_limits<double>::infinity();
    };
    const auto integralBytes = mul((c.t + chunk - 1) / chunk, sizeof(Integral));
    workBudget(add(add(sourceBytes, edgesBytes), integralBytes), o, r);
    std::vector<Integral> integrals((c.t + chunk - 1) / chunk);
    core::parallelChunksWithThreads(
        0, c.t, chunk, o.threads, [&](std::size_t part, std::size_t lo, std::size_t hi) {
            cancel(o);
            for (std::size_t i = lo; i < hi; ++i) {
                std::array<Vec3, 10> p{};
                const auto *v = result.cellNodes(i);
                for (unsigned j = 0; j < 10; ++j)
                    p[j] = result.node(v[j]);
                const auto volume = quadraticVolume(p);
                need(volume > 0 && std::isfinite(volume), PrimitiveReason::InvalidCell);
                integrals[part].volume += volume;
                integrals[part].minimum =
                    std::min(integrals[part].minimum, static_cast<double>(volume));
                if (!positiveTet10JacobianFast(p)) {
                    const auto proof = certifyTet10Jacobian(p, {});
                    need(proof.status == Tet10CertificateStatus::Certified,
                         PrimitiveReason::NumericallyUnresolved);
                }
            }
        });
    long double volume = 0;
    r.minimumCellVolume = std::numeric_limits<double>::infinity();
    for (const auto &part : integrals) {
        volume += part.volume;
        r.minimumCellVolume = std::min(r.minimumCellVolume, part.minimum);
    }
    r.meshVolume = static_cast<double>(volume);
    r.relativeVolumeError = std::abs(r.meshVolume - r.analyticVolume) / r.analyticVolume;
    need(r.relativeVolumeError <= o.maximumRelativeVolumeError, PrimitiveReason::GeometryTolerance);
    r.certifiedTet10s = c.t;
    c = next;
    return result;
}
} // namespace
PrimitiveMeshingResult tetrahedra(const PrimitiveShape &s, const PrimitiveMeshingOptions &o,
                                  mesh::Mesh &destination) noexcept {
    PrimitiveMeshingResult r;
    try {
        if (!o.axialPositions.empty()) {
            r.profiledSizing = true;
            r.transverseTarget = o.targetSize;
            r.minimumAxialSpacing = std::numeric_limits<double>::infinity();
            double previous = 0;
            for (std::size_t i = 1; i < o.axialPositions.size(); ++i) {
                const double width = o.axialPositions[i] - o.axialPositions[i - 1];
                r.minimumAxialSpacing = std::min(r.minimumAxialSpacing, width);
                r.maximumAxialSpacing = std::max(r.maximumAxialSpacing, width);
                if (previous > 0)
                    r.maximumAxialGrowth = std::max(r.maximumAxialGrowth,
                                                    std::max(width / previous, previous / width));
                previous = width;
            }
        }
        Counts c;
        mesh::Mesh linear;
        const auto begin = Clock::now();
        if (o.profileStrategy == PrimitiveProfileStrategy::Volume) {
            const auto reason = profiledCylinderSource(s, o, r, linear);
            need(reason == PrimitiveReason::None, reason);
            c = {r.nodes, r.tet4s, r.boundaryTriangles};
        } else if (s.kind == PrimitiveKind::Cone && o.coneStrategy != PrimitiveConeStrategy::Swept)
            linear = gradedCone(s, o, r, c);
        else if (s.kind == PrimitiveKind::Cylinder &&
                 o.cylinderStrategy != PrimitiveCylinderStrategy::Swept)
            linear = uniformSweep(s, o, r, c);
        else if ((s.kind == PrimitiveKind::Torus || s.kind == PrimitiveKind::Elbow) &&
                 s.innerRadius == 0 && o.circularStrategy != PrimitiveCircularStrategy::Swept)
            linear = uniformSweep(s, o, r, c);
        else {
            auto baseOptions = o;
            baseOptions.elements = PrimitiveElements::HexDominant;
            baseOptions.coneStrategy = PrimitiveConeStrategy::Swept;
            baseOptions.cylinderStrategy = PrimitiveCylinderStrategy::Swept;
            baseOptions.circularStrategy = PrimitiveCircularStrategy::Swept;
            baseOptions.minimumQuality = 0;
            baseOptions.warningQuality = 0;
            baseOptions.qualityPolicy = PrimitiveQualityPolicy::Warn;
            mesh::Mesh baseMesh;
            const auto base = meshPrimitive(s, baseOptions, baseMesh);
            need(base.published, base.reason);
            r = base;
            r.published = false;
            r.destinationUnchanged = true;
            r.plannedWorkingBytes = base.plannedWorkingBytes;
            linear = convert(baseMesh, base, o, r, c);
        }
        r.nodes = c.n;
        r.hexes = 0;
        r.pyramids = 0;
        r.boundaryQuads = 0;
        r.boundaryTriangles = c.b;
        r.tet4s = c.t;
        r.tet10s = 0;
        r.outputBytes = payload(c);
        r.generationSeconds = elapsed(begin);
        const auto validation = Clock::now();
        validate(linear, c, o, r);
        r.validationSeconds = elapsed(validation);
        const auto audit = Clock::now();
        workBudget(auditBudget(c, o), o, r);
        auto limits = o.auditLimits;
        if (o.cancellation)
            limits.cancellation = o.cancellation;
        r.audit = mesh::auditLinearVolumeBoundary(linear, limits);
        r.auditSeconds = elapsed(audit);
        need(!r.audit.cancelled, PrimitiveReason::Cancelled);
        need(!r.audit.allocationFailed, PrimitiveReason::AllocationFailed);
        need(!r.audit.budgetExceeded, PrimitiveReason::WorkingBudget);
        need(r.audit.valid, PrimitiveReason::BoundaryAuditFailed);
        if (o.elements == PrimitiveElements::Tet10) {
            const auto start = Clock::now();
            auto quadratic = elevate(linear, c, o, r);
            r.validationSeconds += elapsed(start);
            linear.swap(quadratic);
            r.tet10s = c.t;
            r.tet4s = 0;
            r.nodes = c.n;
            r.outputBytes = payload(c, true);
        }
        cancel(o);
        destination.swap(linear);
        r.published = true;
        r.destinationUnchanged = false;
        r.reason = PrimitiveReason::None;
    } catch (const Stop &e) {
        r.reason = e.reason;
    } catch (const std::bad_alloc &) {
        r.reason = PrimitiveReason::AllocationFailed;
    } catch (...) {
        r.reason = PrimitiveReason::InternalError;
    }
    return r;
}
} // namespace paacs::algo::primitive_detail
