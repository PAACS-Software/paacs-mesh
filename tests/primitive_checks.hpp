// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "geometry_checks.hpp"
#include "primitive_fixtures.hpp"
#include <bit>
#include <map>
#include <set>
namespace primitive_test {
using geometry_test::require;
inline std::uint64_t hash(const mesh::Mesh &m) {
    std::uint64_t h = 1469598103934665603ULL;
    auto add = [&](std::uint64_t v) {
        for (int j = 0; j < 8; ++j) {
            h ^= v & 255;
            h *= 1099511628211ULL;
            v >>= 8;
        }
    };
    add(m.nodeCount());
    add(m.cellCount());
    for (Index i = 0; i < m.nodeCount(); ++i) {
        const auto p = m.node(i);
        add(std::bit_cast<std::uint64_t>(p.x));
        add(std::bit_cast<std::uint64_t>(p.y));
        add(std::bit_cast<std::uint64_t>(p.z));
        add(m.nodeClass(i).dim);
        add(m.nodeClass(i).tag);
    }
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto &c = m.cell(i);
        add(static_cast<unsigned>(c.type));
        add(c.onEntity.dim);
        add(c.onEntity.tag);
        for (int j = 0; j < mesh::nodesPerCell(c.type); ++j)
            add(m.cellNodes(i)[j]);
    }
    return h;
}
struct Check {
    double volume = 0, maximumSurfaceError = 0;
    long long euler = 0;
};
inline Check check(const algo::PrimitiveShape &s, const mesh::Mesh &m,
                   const algo::PrimitiveMeshingResult &r) {
    using K = algo::PrimitiveKind;
    Check out;
    out.volume = geometry_test::positiveVolume(m);
    const bool circular = s.kind == K::Torus || s.kind == K::Elbow;
    const bool hollow = s.kind == K::Tube || (circular && s.innerRadius > 0);
    const double end = s.kind == K::Cone ? 0 : s.kind == K::Frustum ? s.endRadius : s.radius;
    const Vec3 z = s.axis.normalized(), x = (s.xDirection - z * s.xDirection.dot(z)).normalized(),
               y = z.cross(x);
    const auto local = [&](Vec3 p) {
        p = p - s.origin;
        return Vec3{p.dot(x), p.dot(y), p.dot(z)};
    };
    const double scale =
        std::max({s.radius, circular ? s.majorRadius : s.height, std::abs(s.origin.x),
                  std::abs(s.origin.y), std::abs(s.origin.z)});
    const double eps = 2e-8 * std::max(s.radius, circular ? s.majorRadius : s.height) +
                       64 * std::numeric_limits<double>::epsilon() * scale;
    auto radial = [&](Vec3 p) {
        return circular ? std::hypot(std::hypot(p.x, p.y) - s.majorRadius, p.z)
                        : std::hypot(p.x, p.y);
    };
    auto radiusAt = [&](Vec3 p) {
        return circular ? s.radius : s.radius + (end - s.radius) * p.z / s.height;
    };
    const double theta = 2 * std::numbers::pi / (4 * r.segmentsPerQuarter);
    // Independent analytic oracle: Simpson integration is exact for the
    // quadratic area of an axial frustum; Pappus gives a swept annulus.
    const long double outer = s.radius, inner = hollow ? s.innerRadius : 0;
    const long double mid = (outer + end) / 2;
    const long double area = std::numbers::pi_v<long double> * (outer - inner) * (outer + inner);
    const long double analytic =
        circular ? area * s.majorRadius *
                       (s.kind == K::Torus ? 2 * std::numbers::pi_v<long double> : s.sweepAngle)
        : hollow ? area * s.height
                 : std::numbers::pi_v<long double> * s.height *
                       (outer * outer + 4 * mid * mid + static_cast<long double>(end) * end) / 6;
    require(std::abs(r.analyticVolume - analytic) <= 2e-12L * analytic,
            "independent analytic volume");
    double expected = static_cast<double>(analytic) * std::sin(theta) / theta;
    if (circular) {
        const double a = (s.kind == K::Torus ? 2 * std::numbers::pi : s.sweepAngle) / r.sweepLayers;
        expected *= std::sin(a) / a;
    }
    if (!r.gradedCone)
        require(std::abs(out.volume - expected) <= 2e-8 * r.analyticVolume,
                "independent polygonal volume formula");
    require(std::abs(out.volume - r.meshVolume) <= 2e-8 * r.analyticVolume,
            "independent cell volume");
    const double nodeTolerance = eps + (r.tet10s ? r.chordErrorBound : 0);
    for (Index i = 0; i < m.nodeCount(); ++i) {
        const auto p = local(m.node(i));
        const double radius = radial(p);
        require(radius <= radiusAt(p) + nodeTolerance, "node outside primitive");
        if (hollow)
            require(radius >= s.innerRadius - nodeTolerance, "node inside bore");
        if (!circular)
            require(p.z >= -eps && p.z <= s.height + eps, "axial bounds");
        const auto c = m.nodeClass(i);
        if (c.dim == 2 && (c.tag == 1 || c.tag == 4))
            require(std::abs(radius - (c.tag == 4 ? s.innerRadius : radiusAt(p))) <= nodeTolerance,
                    "boundary node distance");
        if (c.dim == 1)
            require(std::abs(radius - (c.tag >= 103 ? s.innerRadius : radiusAt(p))) <=
                        nodeTolerance,
                    "rim node distance");
        if (c.dim == 0)
            require(s.kind == K::Cone && c.tag == 201 && (p - Vec3{0, 0, s.height}).norm() <= eps,
                    "apex classification");
    }
    std::map<std::pair<Index, Index>, int> edges;
    if (hollow)
        for (Index i = 0; i < m.cellCount(); ++i) {
            if (m.cell(i).onEntity.dim != 3)
                continue;
            Vec3 center{};
            const auto n = mesh::nodesPerCell(m.cell(i).type);
            for (int j = 0; j < n; ++j)
                center = center + m.node(m.cellNodes(i)[j]) / double(n);
            const auto p = local(center);
            const double radius = radial(p);
            require(radius >= s.innerRadius - eps && radius <= s.radius + eps,
                    "thin-wall cell center outside true analytic solid");
        }
    std::set<Index> nodes;
    long long faces = 0;
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto &c = m.cell(i);
        if (c.onEntity.dim != 2)
            continue;
        ++faces;
        const auto *ids = m.cellNodes(i);
        const auto n = c.type == mesh::CellType::Tri6 ? 3 : mesh::nodesPerCell(c.type);
        Vec3 center{};
        for (int j = 0; j < n; ++j) {
            nodes.insert(ids[j]);
            ++edges[std::minmax(ids[j], ids[(j + 1) % n])];
            center = center + m.node(ids[j]) / double(n);
        }
        const auto cp = local(center);
        const Vec3 normal =
            (m.node(ids[1]) - m.node(ids[0])).cross(m.node(ids[2]) - m.node(ids[0]));
        Vec3 outward;
        if (c.onEntity.tag == 1 || c.onEntity.tag == 4) {
            for (int j = 0; j < n; ++j) {
                const auto p = local(m.node(ids[j]));
                require(std::abs(radial(p) - (c.onEntity.tag == 4 ? s.innerRadius : radiusAt(p))) <=
                            eps,
                        "face corner on correct surface");
            }
            const double expectedRadius = c.onEntity.tag == 4 ? s.innerRadius : radiusAt(cp);
            const double error = std::abs(radial(cp) - expectedRadius);
            out.maximumSurfaceError = std::max(out.maximumSurfaceError, error);
            require(error <= r.chordErrorBound + eps, "face center chord-error bound");
            if (circular) {
                const double rho = std::hypot(cp.x, cp.y);
                outward = x * (cp.x * (rho - s.majorRadius) / rho) +
                          y * (cp.y * (rho - s.majorRadius) / rho) + z * cp.z;
            } else
                outward = x * cp.x + y * cp.y -
                          z * ((end - s.radius) * std::hypot(cp.x, cp.y) / s.height);
            if (c.onEntity.tag == 4)
                outward = outward * (-1);
        } else if (circular) {
            const double a = c.onEntity.tag == 2 ? 0 : s.sweepAngle;
            outward = (x * (-std::sin(a)) + y * std::cos(a)) * (c.onEntity.tag == 2 ? -1 : 1);
        } else
            outward = z * (c.onEntity.tag == 2 ? -1 : 1);
        require(normal.dot(outward) > 0, "outward boundary orientation");
    }
    for (const auto &[edge, uses] : edges) {
        (void)edge;
        require(uses == 2, "boundary seam/crack/nonmanifold edge");
    }
    out.euler = static_cast<long long>(nodes.size()) - static_cast<long long>(edges.size()) + faces;
    require(out.euler == ((s.kind == K::Torus || hollow) ? 0 : 2), "boundary Euler characteristic");
    return out;
}
} // namespace primitive_test
