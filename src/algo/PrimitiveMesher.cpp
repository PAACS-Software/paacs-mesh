// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/PrimitiveMesher.hpp"
#include "PrimitiveMesherInternal.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/CertifiedMeshStorage.hpp"
#include "paacs/mesh/Hex8Jacobian.hpp"
#include "paacs/predicates/Predicates.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <new>
#include <vector>

namespace paacs::algo {
namespace {
using Clock = std::chrono::steady_clock;
constexpr double pi = std::numbers::pi;
constexpr std::uint64_t chunkSize = 1024;
struct Stop {
    PrimitiveReason reason;
};
void require(bool ok, PrimitiveReason reason) {
    if (!ok)
        throw Stop{reason};
}
bool finite(Vec3 p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
void cancelled(const PrimitiveMeshingOptions &o) {
    require(!(o.cancellation && o.cancellation->load(std::memory_order_relaxed)) &&
                !(o.auditLimits.cancellation &&
                  o.auditLimits.cancellation->load(std::memory_order_relaxed)),
            PrimitiveReason::Cancelled);
}
std::uint64_t add(std::uint64_t a, std::uint64_t b) {
    require(b <= (std::numeric_limits<std::uint64_t>::max)() - a, PrimitiveReason::OutputBudget);
    return a + b;
}
std::uint64_t mul(std::uint64_t a, std::uint64_t b) {
    require(a == 0 || b <= (std::numeric_limits<std::uint64_t>::max)() / a,
            PrimitiveReason::OutputBudget);
    return a * b;
}
std::uint32_t count(double n, std::uint32_t minimum = 1) {
    require(std::isfinite(n) && n <= 1'000'000 && n >= 0, PrimitiveReason::OutputBudget);
    return std::max(minimum, static_cast<std::uint32_t>(std::ceil(n)));
}
double seconds(Clock::time_point begin) {
    return std::chrono::duration<double>(Clock::now() - begin).count();
}
struct Section {
    std::vector<Vec2> points;
    std::vector<std::uint8_t> boundary; // 0 interior, 1 outer, 2 inner
    std::vector<std::array<Index, 4>> quads;
    std::vector<std::array<Index, 2>> edges; // outer CCW, inner CW
};
Section section(std::uint32_t n, std::uint32_t m, double inner) {
    Section s;
    const std::uint64_t perimeter = 4ULL * n;
    const bool hollow = inner > 0;
    const auto coreNodes = hollow ? perimeter : (n + 1ULL) * (n + 1ULL);
    s.points.resize(coreNodes + perimeter * m);
    s.boundary.resize(s.points.size());
    s.quads.reserve((hollow ? 0 : n * static_cast<std::uint64_t>(n)) + perimeter * m);
    s.edges.reserve(perimeter * (hollow ? 2 : 1));
    auto core = [n](std::uint32_t i, std::uint32_t j) { return Index(j * (n + 1ULL) + i); };
    auto ring0 = [&](std::uint64_t k) -> Index {
        if (hollow)
            return static_cast<Index>(k);
        const auto side = k / n;
        const auto u = static_cast<std::uint32_t>(k % n);
        if (side == 0)
            return core(u, 0);
        if (side == 1)
            return core(n, u);
        if (side == 2)
            return core(n - u, n);
        return core(0, n - u);
    };
    auto ring = [&](std::uint32_t j, std::uint64_t k) -> Index {
        return j == 0 ? ring0(k) : static_cast<Index>(coreNodes + (j - 1) * perimeter + k);
    };
    if (!hollow) {
        for (std::uint32_t j = 0; j <= n; ++j)
            for (std::uint32_t i = 0; i <= n; ++i)
                s.points[core(i, j)] = {double(i) / n - .5, double(j) / n - .5};
        for (std::uint32_t j = 0; j < n; ++j)
            for (std::uint32_t i = 0; i < n; ++i)
                s.quads.push_back({core(i, j), core(i + 1, j), core(i + 1, j + 1), core(i, j + 1)});
    }
    for (std::uint64_t k = 0; k < perimeter; ++k) {
        const double a = -3 * pi / 4 + 2 * pi * double(k) / double(perimeter);
        const Vec2 outer{std::cos(a), std::sin(a)};
        if (hollow) {
            s.points[k] = outer * inner;
            s.boundary[k] = 2;
        }
        const Vec2 first = s.points[ring0(k)];
        for (std::uint32_t j = 1; j <= m; ++j) {
            const double t = double(j) / m;
            const auto id = ring(j, k);
            s.points[id] = first * (1 - t) + outer * t;
            if (j == m)
                s.boundary[id] = 1;
            s.quads.push_back({ring(j - 1, k), id, ring(j, (k + 1) % perimeter),
                               ring(j - 1, (k + 1) % perimeter)});
        }
        s.edges.push_back({ring(m, k), ring(m, (k + 1) % perimeter)});
    }
    if (hollow)
        for (std::uint64_t k = 0; k < perimeter; ++k)
            s.edges.push_back({ring(0, (k + 1) % perimeter), ring(0, k)});
    return s;
}

struct Metrics {
    PrimitiveQualityReport quality{};
    long double volume = 0;
    double minimumVolume = std::numeric_limits<double>::infinity(), minimumQuality = 1;
};
void tet(Metrics &m, Vec3 a, Vec3 b, Vec3 c, Vec3 d) {
    const auto ab = b - a, ac = c - a, ad = d - a;
    const double volume = ab.cross(ac).dot(ad) / 6;
    require(volume > 0 && std::isfinite(volume), PrimitiveReason::InvalidCell);
    const double sum =
        ab.norm2() + ac.norm2() + ad.norm2() + (c - b).norm2() + (d - b).norm2() + (d - c).norm2();
    const double quality = 12 * std::pow(3 * volume, 2. / 3) / sum;
    require(std::isfinite(quality) && quality > 0, PrimitiveReason::NumericallyUnresolved);
    m.volume += volume;
    m.minimumQuality = std::min(m.minimumQuality, quality);
}
Metrics validate(const mesh::Mesh &mesh, std::uint64_t volumes,
                 const PrimitiveMeshingOptions &options) {
    std::vector<Metrics> parts((volumes + chunkSize - 1) / chunkSize);
    core::parallelChunksWithThreads(
        0, volumes, chunkSize, options.threads,
        [&](std::size_t ci, std::size_t lo, std::size_t hi) {
            cancelled(options);
            auto &part = parts[ci];
            for (std::size_t i = lo; i < hi; ++i) {
                const auto *ids = mesh.cellNodes(static_cast<Index>(i));
                std::array<Vec3, 8> p{};
                const bool hex = mesh.cell(static_cast<Index>(i)).type == mesh::CellType::Hex8;
                for (int j = 0; j < (hex ? 8 : 5); ++j)
                    p[j] = mesh.node(ids[j]);
                Metrics cell;
                if (hex) {
                    require(mesh::positiveHex8Jacobian(p), PrimitiveReason::NumericallyUnresolved);
                    constexpr std::array<std::array<unsigned, 4>, 6> t{{{0, 1, 2, 6},
                                                                        {0, 2, 3, 6},
                                                                        {0, 3, 7, 6},
                                                                        {0, 7, 4, 6},
                                                                        {0, 4, 5, 6},
                                                                        {0, 5, 1, 6}}};
                    for (const auto &v : t)
                        tet(cell, p[v[0]], p[v[1]], p[v[2]], p[v[3]]);
                } else {
                    for (unsigned j = 0; j < 4; ++j)
                        require(predicates::orient3d(p[j], p[(j + 1) % 4], p[(j + 2) % 4], p[4]) ==
                                    predicates::Sign::Positive,
                                PrimitiveReason::InvalidCell);
                    // Convex planar quad base + apex off its plane. Check
                    // both diagonals, then count just one volume partition.
                    Metrics alternate;
                    tet(alternate, p[1], p[2], p[3], p[4]);
                    tet(alternate, p[1], p[3], p[0], p[4]);
                    tet(cell, p[0], p[1], p[2], p[4]);
                    tet(cell, p[0], p[2], p[3], p[4]);
                    cell.minimumQuality = std::min(cell.minimumQuality, alternate.minimumQuality);
                }
                bool sharp = false;
                for (int j = 0; j < (hex ? 8 : 5); ++j)
                    sharp = sharp || mesh.nodeClass(ids[j]).dim <= 1;
                primitive_detail::recordQuality(part.quality, static_cast<Index>(i),
                                                cell.minimumQuality, sharp, options.warningQuality);
                part.volume += cell.volume;
                part.minimumVolume = std::min(part.minimumVolume, static_cast<double>(cell.volume));
                part.minimumQuality = std::min(part.minimumQuality, cell.minimumQuality);
            }
        });
    Metrics result;
    for (const auto &p : parts) {
        primitive_detail::mergeQuality(result.quality, p.quality);
        result.volume += p.volume;
        result.minimumVolume = std::min(result.minimumVolume, p.minimumVolume);
        result.minimumQuality = std::min(result.minimumQuality, p.minimumQuality);
    }
    return result;
}
} // namespace

PrimitiveMeshingResult meshPrimitive(const PrimitiveShape &s, const PrimitiveMeshingOptions &o,
                                     mesh::Mesh &destination) noexcept {
    PrimitiveMeshingResult out;
    try {
        cancelled(o);
        const bool circular = s.kind == PrimitiveKind::Torus || s.kind == PrimitiveKind::Elbow;
        const bool closed = s.kind == PrimitiveKind::Torus, cone = s.kind == PrimitiveKind::Cone;
        const bool hollow = s.kind == PrimitiveKind::Tube || (circular && s.innerRadius > 0);
        require(static_cast<unsigned>(s.kind) <= static_cast<unsigned>(PrimitiveKind::Elbow) &&
                    finite(s.origin) && finite(s.axis) && finite(s.xDirection) &&
                    std::isfinite(s.radius) && s.radius > 0 && std::isfinite(o.targetSize) &&
                    o.targetSize > 0 && std::isfinite(o.maximumChordError) &&
                    o.maximumChordError >= 0 && std::isfinite(o.maximumRelativeVolumeError) &&
                    o.maximumRelativeVolumeError > 0 && o.maximumRelativeVolumeError < 1 &&
                    std::isfinite(o.minimumQuality) && o.minimumQuality >= 0 &&
                    o.minimumQuality <= 1 && o.threads >= 1 && o.threads <= 1024,
                PrimitiveReason::InvalidRequest);
        const auto unit = [](Vec3 v) {
            const double n = std::hypot(v.x, v.y, v.z);
            require(n > 0 && std::isfinite(n), PrimitiveReason::InvalidRequest);
            return v / n;
        };
        const Vec3 z = unit(s.axis), hinted = unit(s.xDirection);
        const Vec3 projected = hinted - z * hinted.dot(z);
        require(projected.norm() >= 1e-8, PrimitiveReason::InvalidRequest);
        const Vec3 x = unit(projected), y = z.cross(x);
        const double r1 = cone ? 0 : s.kind == PrimitiveKind::Frustum ? s.endRadius : s.radius;
        const double rmax = std::max(s.radius, r1);
        require(std::isfinite(r1) && (cone || r1 > 0), PrimitiveReason::InvalidRequest);
        if (circular || s.kind == PrimitiveKind::Tube)
            require(std::isfinite(s.innerRadius) && s.innerRadius >= 0 &&
                        s.innerRadius < s.radius && (!hollow || s.innerRadius > 0),
                    PrimitiveReason::InvalidRequest);
        if (s.kind == PrimitiveKind::Tube)
            require(s.innerRadius > 0, PrimitiveReason::InvalidRequest);
        if (hollow)
            require(s.innerRadius / s.radius > 0 && s.innerRadius / s.radius < 1,
                    PrimitiveReason::NumericallyUnresolved);
        if (circular)
            require(std::isfinite(s.majorRadius) && s.majorRadius > s.radius,
                    PrimitiveReason::InvalidRequest);
        else
            require(std::isfinite(s.height) && s.height > 0, PrimitiveReason::InvalidRequest);
        const double angle = closed ? 2 * pi : s.sweepAngle;
        if (s.kind == PrimitiveKind::Elbow)
            require(std::isfinite(angle) && angle > 0 && angle < 2 * pi,
                    PrimitiveReason::InvalidRequest);
        require(
            static_cast<unsigned>(o.elements) <= static_cast<unsigned>(PrimitiveElements::Tet10) &&
                static_cast<unsigned>(o.coneStrategy) <=
                    static_cast<unsigned>(PrimitiveConeStrategy::Graded) &&
                static_cast<unsigned>(o.cylinderStrategy) <=
                    static_cast<unsigned>(PrimitiveCylinderStrategy::Uniform) &&
                static_cast<unsigned>(o.circularStrategy) <=
                    static_cast<unsigned>(PrimitiveCircularStrategy::Uniform) &&
                static_cast<unsigned>(o.profileStrategy) <=
                    static_cast<unsigned>(PrimitiveProfileStrategy::Volume) &&
                (o.profileStrategy != PrimitiveProfileStrategy::Volume ||
                 (o.sizeProfile && s.kind == PrimitiveKind::Cylinder)) &&
                static_cast<unsigned>(o.qualityPolicy) <=
                    static_cast<unsigned>(PrimitiveQualityPolicy::Reject) &&
                std::isfinite(o.warningQuality) && o.warningQuality >= 0 && o.warningQuality <= 1 &&
                (o.coneStrategy != PrimitiveConeStrategy::Graded ||
                 (cone && o.elements != PrimitiveElements::HexDominant)) &&
                (o.cylinderStrategy != PrimitiveCylinderStrategy::Uniform ||
                 (s.kind == PrimitiveKind::Cylinder &&
                  o.elements != PrimitiveElements::HexDominant)) &&
                (o.circularStrategy != PrimitiveCircularStrategy::Uniform ||
                 (circular && !hollow && o.elements != PrimitiveElements::HexDominant)),
            PrimitiveReason::InvalidRequest);
        if (o.sizeProfile || !o.axialPositions.empty()) {
            require((s.kind == PrimitiveKind::Cylinder || s.kind == PrimitiveKind::Tube) &&
                        (!o.sizeProfile || o.axialPositions.empty()),
                    PrimitiveReason::InvalidRequest);
            if (o.sizeProfile)
                return primitive_detail::profiled(s, o, destination);
            require(o.axialPositions.size() >= 2 && o.axialPositions.size() <= 1'000'001 &&
                        o.axialPositions.front() == 0 && o.axialPositions.back() == s.height &&
                        std::isfinite(o.maximumAxialGrowth) && o.maximumAxialGrowth >= 1,
                    PrimitiveReason::InvalidRequest);
            double previous = 0;
            out.minimumAxialSpacing = std::numeric_limits<double>::infinity();
            for (std::size_t k = 1; k < o.axialPositions.size(); ++k) {
                if (k % 1024 == 0)
                    cancelled(o);
                const double h = o.axialPositions[k] - o.axialPositions[k - 1];
                require(std::isfinite(h) && h > 0, PrimitiveReason::InvalidRequest);
                const double g = previous > 0 ? std::max(h / previous, previous / h) : 1;
                require(std::isfinite(g) && g <= o.maximumAxialGrowth * (1 + 1e-10), PrimitiveReason::InvalidRequest);
                out.minimumAxialSpacing = std::min(out.minimumAxialSpacing, h);
                out.maximumAxialSpacing = std::max(out.maximumAxialSpacing, h);
                out.maximumAxialGrowth = std::max(out.maximumAxialGrowth, g);
                previous = h;
            }
            out.profiledSizing = true;
            out.transverseTarget = o.targetSize;
        }
        if (o.elements != PrimitiveElements::HexDominant)
            return primitive_detail::tetrahedra(s, o, destination);
        out.radialLayers = count((hollow ? s.radius - s.innerRadius : .5 * rmax) / o.targetSize);
        const double requestedTolerance =
            o.maximumChordError > 0 ? o.maximumChordError : .0025 * rmax;
        // A globally small sagitta can still swallow a thin wall. Keep the
        // combined angular error below a quarter of one radial layer.
        const double tolerance =
            hollow
                ? std::min(requestedTolerance, (s.radius - s.innerRadius) / (4 * out.radialLayers))
                : requestedTolerance;
        const double minorTolerance = circular ? tolerance * .5 : tolerance;
        // sagitta = 2*r*sin(angle/4)^2; asin avoids cancellation for tiny errors.
        auto maxAngle = [](double radius, double error) {
            return std::min(pi / 4, 4 * std::asin(std::sqrt(std::min(.5, error / (2 * radius)))));
        };
        // 1-sinc(a) <= a*a/6. Split the volume-error allowance between
        // the two circles of a torus; satisfy tighter volume requests by
        // refining before allocation, rather than only rejecting afterward.
        const double volumeAngle = std::sqrt(6 * o.maximumRelativeVolumeError / (circular ? 2 : 1));
        const double minorAngle = std::min(maxAngle(rmax, minorTolerance), volumeAngle);
        out.segmentsPerQuarter =
            count(std::max(pi * rmax / (2 * o.targetSize), pi / (2 * minorAngle)), 2);
        out.sweepLayers =
            circular
                ? count(
                      std::max(angle * (s.majorRadius + s.radius) / o.targetSize,
                               angle / std::min(maxAngle(s.majorRadius + s.radius, tolerance * .5),
                                                volumeAngle)),
                      closed ? 8 : 1)
            : !o.axialPositions.empty() ? static_cast<std::uint32_t>(o.axialPositions.size() - 1)
                                        : count(std::hypot(s.height, r1 - s.radius) / o.targetSize);
        const std::uint64_t n = out.segmentsPerQuarter, m = out.radialLayers, l = out.sweepLayers,
                            p = 4 * n;
        const auto sectionNodes = add(hollow ? p : mul(n + 1, n + 1), mul(p, m));
        const auto sectionCells = add(hollow ? 0 : mul(n, n), mul(p, m));
        const auto boundaryEdges = mul(p, hollow ? 2 : 1);
        const auto planes = closed || cone ? l : l + 1;
        out.nodes = add(mul(sectionNodes, planes), cone ? 1 : 0);
        out.hexes = mul(sectionCells, cone ? l - 1 : l);
        out.pyramids = cone ? sectionCells : 0;
        out.boundaryTriangles = cone ? boundaryEdges : 0;
        out.boundaryQuads =
            add(mul(boundaryEdges, cone ? l - 1 : l), mul(sectionCells, closed ? 0
                                                                        : cone ? 1
                                                                               : 2));
        const auto volumes = add(out.hexes, out.pyramids);
        const auto boundaries = add(out.boundaryQuads, out.boundaryTriangles),
                   cells = add(volumes, boundaries);
        const auto conn = add(add(mul(out.hexes, 8), mul(out.pyramids, 5)),
                              add(mul(out.boundaryQuads, 4), mul(out.boundaryTriangles, 3)));
        out.outputBytes = add(add(mul(out.nodes, sizeof(Vec3) + sizeof(mesh::Classification)),
                                  mul(cells, sizeof(mesh::Cell))),
                              mul(conn, sizeof(Index)));
        require(out.nodes <= o.maximumNodes && cells <= o.maximumCells &&
                    out.outputBytes <= o.maximumOutputBytes &&
                    out.nodes <=
                        static_cast<std::uint64_t>((std::numeric_limits<std::int32_t>::max)()) &&
                    out.outputBytes <=
                        static_cast<std::uint64_t>((std::numeric_limits<std::size_t>::max)()) &&
                    conn <= static_cast<std::uint64_t>((std::numeric_limits<Index>::max)()),
                PrimitiveReason::OutputBudget);
        const auto faceUses = add(add(mul(out.hexes, 6), mul(out.pyramids, 5)), boundaries);
        const auto auditBytes = add(mul(faceUses, 16), mul(out.nodes, sizeof(std::size_t)));
        out.plannedWorkingBytes = add(
            auditBytes, add(add(mul(sectionNodes, sizeof(Vec2) + sizeof(std::uint8_t)),
                                mul(sectionCells, 4 * sizeof(Index))),
                            add(mul(boundaryEdges, 2 * sizeof(Index)),
                                add(mul(planes, sizeof(Vec2)),
                                    mul((volumes + chunkSize - 1) / chunkSize, sizeof(Metrics))))));
        require(faceUses <= o.auditLimits.maximumFaceUses &&
                    auditBytes <= o.auditLimits.maximumTemporaryBytes &&
                    out.plannedWorkingBytes <= o.maximumWorkingBytes,
                PrimitiveReason::WorkingBudget);
        const double minorStep = pi / (2 * n);
        out.chordErrorBound = 2 * rmax * std::pow(std::sin(minorStep / 4), 2);
        if (circular)
            out.chordErrorBound +=
                2 * (s.majorRadius + s.radius) * std::pow(std::sin(angle / (4 * l)), 2);
        const double coordinateScale =
            std::max({std::abs(s.origin.x), std::abs(s.origin.y), std::abs(s.origin.z), rmax,
                      circular ? s.majorRadius : s.height});
        const double roundoffBound = 128 * std::numeric_limits<double>::epsilon() * coordinateScale;
        out.chordErrorBound += roundoffBound;
        require(std::isfinite(out.chordErrorBound) &&
                    out.chordErrorBound <= tolerance * (1 + 1e-12),
                PrimitiveReason::GeometryTolerance);
        const double area =
            pi * (s.radius * s.radius - (hollow ? s.innerRadius * s.innerRadius : 0));
        out.analyticVolume =
            circular ? area * s.majorRadius * angle
            : hollow ? area * s.height
                     : pi * s.height * (s.radius * s.radius + s.radius * r1 + r1 * r1) / 3;
        require(std::isfinite(out.analyticVolume) && out.analyticVolume > 0,
                PrimitiveReason::NumericallyUnresolved);
        mesh::Mesh local;
        const auto generationBegin = Clock::now();
        {
            const Section a = section(out.segmentsPerQuarter, out.radialLayers,
                                      hollow ? s.innerRadius / s.radius : 0);
            std::vector<Vec2> trig(planes);
            if (circular)
                for (std::uint64_t k = 0; k < planes; ++k)
                    trig[k] = {std::cos(angle * double(k) / double(l)),
                               std::sin(angle * double(k) / double(l))};
            mesh::detail::CertifiedMeshStorage storage(out.nodes, cells, conn);
            core::parallelChunksWithThreads(
                0, out.nodes, chunkSize, o.threads,
                [&](std::size_t, std::size_t lo, std::size_t hi) {
                    cancelled(o);
                    for (std::size_t ordinal = lo; ordinal < hi; ++ordinal) {
                        if (cone && ordinal + 1 == out.nodes) {
                            storage.setNode(ordinal, s.origin + z * s.height, {0, 201});
                            continue;
                        }
                        const auto k = ordinal / sectionNodes, j = ordinal % sectionNodes;
                        const double t = o.axialPositions.empty() ? double(k) / double(l)
                                                                  : o.axialPositions[k] / s.height,
                                     radius = s.radius * (1 - t) + r1 * t;
                        const auto uv = a.points[j];
                        Vec3 point;
                        if (circular)
                            point = s.origin +
                                    (x * trig[k].u + y * trig[k].v) *
                                        (s.majorRadius + s.radius * uv.u) -
                                    z * (s.radius * uv.v);
                        else
                            point = s.origin + x * (radius * uv.u) + y * (radius * uv.v) +
                                    z * (s.height * t);
                        require(finite(point), PrimitiveReason::NumericallyUnresolved);
                        mesh::Classification tag{3, 1};
                        const auto b = a.boundary[j];
                        if (b)
                            tag = {2, b == 1 ? 1 : 4};
                        if (!closed && (k == 0 || (!cone && k == l)))
                            tag =
                                b ? mesh::Classification{1, (b == 1 ? 101 : 103) + (k == l ? 1 : 0)}
                                  : mesh::Classification{2, k == 0 ? 2 : 3};
                        storage.setNode(ordinal, point, tag);
                    }
                });
            core::parallelChunksWithThreads(
                0, volumes, chunkSize, o.threads, [&](std::size_t, std::size_t lo, std::size_t hi) {
                    cancelled(o);
                    for (std::size_t ordinal = lo; ordinal < hi; ++ordinal) {
                        const auto layer = ordinal / sectionCells, quad = ordinal % sectionCells;
                        const auto &q = a.quads[quad];
                        const Index base = static_cast<Index>(layer * sectionNodes);
                        if (ordinal < out.hexes) {
                            const Index top =
                                static_cast<Index>(((layer + 1) % planes) * sectionNodes);
                            storage.setCell(ordinal, ordinal * 8, mesh::CellType::Hex8,
                                            std::array<Index, 8>{
                                                base + q[0], base + q[1], base + q[2], base + q[3],
                                                top + q[0], top + q[1], top + q[2], top + q[3]},
                                            {3, 1});
                        } else
                            storage.setCell(ordinal, out.hexes * 8 + (ordinal - out.hexes) * 5,
                                            mesh::CellType::Pyramid5,
                                            std::array<Index, 5>{base + q[0], base + q[1],
                                                                 base + q[2], base + q[3],
                                                                 static_cast<Index>(out.nodes - 1)},
                                            {3, 1});
                    }
                });
            const auto volumeConn = out.hexes * 8 + out.pyramids * 5;
            const auto sideQuads = boundaryEdges * (cone ? l - 1 : l);
            core::parallelChunksWithThreads(
                0, out.boundaryQuads, chunkSize, o.threads,
                [&](std::size_t, std::size_t lo, std::size_t hi) {
                    cancelled(o);
                    for (std::size_t i = lo; i < hi; ++i) {
                        std::array<Index, 4> q{};
                        int tag = 1;
                        if (i < sideQuads) {
                            const auto layer = i / boundaryEdges, edge = i % boundaryEdges;
                            const auto &e = a.edges[edge];
                            const Index base = static_cast<Index>(layer * sectionNodes),
                                        top = static_cast<Index>(((layer + 1) % planes) *
                                                                 sectionNodes);
                            q = {base + e[0], base + e[1], top + e[1], top + e[0]};
                            tag = edge < p ? 1 : 4;
                        } else {
                            const auto cap = (i - sideQuads) / sectionCells,
                                       j = (i - sideQuads) % sectionCells;
                            const auto &v = a.quads[j];
                            if (cap == 0) {
                                q = {v[3], v[2], v[1], v[0]};
                                tag = 2;
                            } else {
                                const Index b = static_cast<Index>(l * sectionNodes);
                                q = {b + v[0], b + v[1], b + v[2], b + v[3]};
                                tag = 3;
                            }
                        }
                        storage.setCell(volumes + i, volumeConn + i * 4, mesh::CellType::Quad4, q,
                                        {2, tag});
                    }
                });
            for (std::uint64_t i = 0; i < out.boundaryTriangles; ++i) {
                if (i % chunkSize == 0)
                    cancelled(o);
                const auto &e = a.edges[i];
                const Index b = static_cast<Index>((l - 1) * sectionNodes);
                storage.setCell(
                    volumes + out.boundaryQuads + i, volumeConn + out.boundaryQuads * 4 + i * 3,
                    mesh::CellType::Tri3,
                    std::array<Index, 3>{b + e[0], b + e[1], static_cast<Index>(out.nodes - 1)},
                    {2, 1});
            }
            storage.publishTo(local);
        }
        out.generationSeconds = seconds(generationBegin);
        const auto validationBegin = Clock::now();
        const auto quality = validate(local, volumes, o);
        out.meshVolume = static_cast<double>(quality.volume);
        out.minimumCellVolume = quality.minimumVolume;
        out.minimumQuality = quality.minimumQuality;
        out.quality = quality.quality;
        out.qualityWarning = out.quality.lowQualityCells != 0;
        out.relativeVolumeError =
            std::abs(out.meshVolume - out.analyticVolume) / out.analyticVolume;
        require(out.minimumQuality >= o.minimumQuality &&
                    !(out.qualityWarning && o.qualityPolicy == PrimitiveQualityPolicy::Reject),
                PrimitiveReason::QualityRequirement);
        require(out.relativeVolumeError <= o.maximumRelativeVolumeError,
                PrimitiveReason::GeometryTolerance);
        out.validationSeconds = seconds(validationBegin);
        cancelled(o);
        const auto auditBegin = Clock::now();
        auto limits = o.auditLimits;
        if (o.cancellation)
            limits.cancellation = o.cancellation;
        out.audit = mesh::auditLinearVolumeBoundary(local, limits);
        out.auditSeconds = seconds(auditBegin);
        require(!out.audit.cancelled, PrimitiveReason::Cancelled);
        require(!out.audit.allocationFailed, PrimitiveReason::AllocationFailed);
        require(!out.audit.budgetExceeded, PrimitiveReason::WorkingBudget);
        require(out.audit.valid, PrimitiveReason::BoundaryAuditFailed);
        cancelled(o);
        destination.swap(local);
        out.published = true;
        out.destinationUnchanged = false;
        out.reason = PrimitiveReason::None;
    } catch (const Stop &s) {
        out.reason = s.reason;
    } catch (const std::bad_alloc &) {
        out.reason = PrimitiveReason::AllocationFailed;
    } catch (...) {
        out.reason = PrimitiveReason::InternalError;
    }
    return out;
}
const char *primitiveKindName(PrimitiveKind k) noexcept {
    switch (k) {
    case PrimitiveKind::Cylinder:
        return "cylinder";
    case PrimitiveKind::Cone:
        return "cone";
    case PrimitiveKind::Frustum:
        return "frustum";
    case PrimitiveKind::Tube:
        return "tube";
    case PrimitiveKind::Torus:
        return "torus";
    case PrimitiveKind::Elbow:
        return "elbow";
    }
    return "unknown";
}
const char *primitiveReasonName(PrimitiveReason r) noexcept {
    switch (r) {
    case PrimitiveReason::None:
        return "none";
    case PrimitiveReason::InvalidRequest:
        return "invalid-request";
    case PrimitiveReason::OutputBudget:
        return "output-budget";
    case PrimitiveReason::WorkingBudget:
        return "working-budget";
    case PrimitiveReason::Cancelled:
        return "cancelled";
    case PrimitiveReason::NumericallyUnresolved:
        return "numerically-unresolved";
    case PrimitiveReason::InvalidCell:
        return "invalid-cell";
    case PrimitiveReason::BoundaryAuditFailed:
        return "boundary-audit-failed";
    case PrimitiveReason::GeometryTolerance:
        return "geometry-tolerance";
    case PrimitiveReason::QualityRequirement:
        return "quality-requirement";
    case PrimitiveReason::AllocationFailed:
        return "allocation-failed";
    case PrimitiveReason::InternalError:
        return "internal-error";
    }
    return "unknown";
}
const char *primitiveElementsName(PrimitiveElements e) noexcept {
    switch (e) {
    case PrimitiveElements::HexDominant:
        return "hex-dominant";
    case PrimitiveElements::Tet4:
        return "tet4";
    case PrimitiveElements::Tet10:
        return "tet10";
    }
    return "unknown";
}
} // namespace paacs::algo
